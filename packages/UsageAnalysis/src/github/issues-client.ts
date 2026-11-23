import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

import {GITHUB_ORG, GitHubIssue, GitHubIssueRef, ISSUE_PRIORITIES, ISSUE_REPOS, IssueRef, IssueResolution,
  isIssueRepo} from './issue';

const GRAPHQL_URL = 'https://api.github.com/graphql';
const TOKEN_CREDENTIAL = 'githubToken';
const BOARD_PROJECT = 17;
const CACHE_SIZE = 1000;
const ISSUES_PER_LOOKUP = 50;
const KEYS_PER_SEARCH = 5;
const SEARCH_RESULTS_LIMIT = 1000;
const MOVED_FROM = /^>\s*Moved from \[(GROK-\d+)\]/m;
const REPORT_MARKER = /^datagrok-report: \S+\/([\w-]+)\s*$/m;
const REPORT_ISSUES_TTL_MS = 10 * 60 * 1000;

const ISSUE_FIELDS = `number title url state stateReason body @include(if: $withBody)
  repository { name }
  milestone { title }
  issueType { name }
  labels(first: 20) { nodes { name } }
  assignees(first: 10) { nodes { login } }
  issueFieldValues(first: 20) { nodes {
    ... on IssueFieldSingleSelectValue { name field { ... on IssueFieldSingleSelect { name } } }
    ... on IssueFieldMultiSelectValue { options { name } field { ... on IssueFieldMultiSelect { name } } }
  } }
  projectItems(first: 10) { nodes {
    project { number }
    status: fieldValueByName(name: "Status") { ... on ProjectV2ItemFieldSingleSelectValue { name } }
  } }`;

const SEARCH_QUERY = `query($q: String!, $cursor: String, $withBody: Boolean!) {
  search(type: ISSUE, query: $q, first: 100, after: $cursor) {
    issueCount
    pageInfo { hasNextPage endCursor }
    nodes { ... on Issue { ${ISSUE_FIELDS} } }
  }
}`;

const REPORT_ISSUES_QUERY = `query($q: String!, $cursor: String) {
  search(type: ISSUE, query: $q, first: 100, after: $cursor) {
    pageInfo { hasNextPage endCursor }
    nodes { ... on Issue { number url state body repository { name } } }
  }
}`;

export interface ReportIssue {
  key: string;
  url: string;
  open: boolean;
}

interface ReportIssueNode {
  number: number;
  url: string;
  state: string;
  body: string;
  repository: Named;
}

interface ReportIssuesData {
  search: {pageInfo: {hasNextPage: boolean; endCursor: string | null}; nodes: (ReportIssueNode | Record<string, never>)[]};
}

interface GraphQLError {
  message: string;
  type?: string;
}

interface GraphQLResponse<T> {
  data?: T | null;
  errors?: GraphQLError[];
}

interface Named {
  name: string;
}

interface IssueFieldValueNode {
  field?: Named;
  name?: string;
  options?: Named[];
}

interface IssueNode {
  number: number;
  title: string;
  url: string;
  state: 'OPEN' | 'CLOSED';
  stateReason: string | null;
  body?: string;
  repository: Named;
  milestone: {title: string} | null;
  issueType: Named | null;
  labels: {nodes: Named[]};
  assignees: {nodes: {login: string}[]};
  issueFieldValues: {nodes: IssueFieldValueNode[]};
  projectItems: {nodes: {project: {number: number}, status: Partial<Named> | null}[]};
}

interface SearchData {
  search: {
    issueCount: number;
    pageInfo: {hasNextPage: boolean, endCursor: string | null};
    nodes: (IssueNode | Record<string, never>)[];
  };
}

type LookupData = Record<string, {issue: IssueNode | null} | null>;

export class GitHubIssuesClient {
  pkg: DG.Package;
  token: Promise<string> | null = null;
  cache = new DG.LruCache<string, GitHubIssue | null>(CACHE_SIZE);
  reportIssuesCache: {at: number; issues: Promise<Map<string, ReportIssue>>} | null = null;

  constructor(pkg: DG.Package) {
    this.pkg = pkg;
  }

  async milestoneIssues(milestone: string): Promise<GitHubIssue[]> {
    const scope = `${repoScope()} is:issue milestone:"${milestone.replace(/"/g, '')}"`;
    const [open, closed] = await Promise.all([
      this.search(`${scope} is:open`, false), this.search(`${scope} is:closed`, false)]);
    return [...open, ...closed].map(toIssue);
  }

  /** Issues ticket-relay filed for reports, by report id. Newest first up to GitHub's 1000 search results;
   *  an open issue wins over a closed one. */
  reportIssues(): Promise<Map<string, ReportIssue>> {
    if (!this.reportIssuesCache || Date.now() - this.reportIssuesCache.at > REPORT_ISSUES_TTL_MS) {
      const issues = this.searchReportIssues();
      issues.catch(() => this.reportIssuesCache = null);
      this.reportIssuesCache = {at: Date.now(), issues};
    }
    return this.reportIssuesCache.issues;
  }

  async searchReportIssues(): Promise<Map<string, ReportIssue>> {
    const q = `repo:${GITHUB_ORG}/core is:issue in:body "datagrok-report" sort:created-desc`;
    const byReport = new Map<string, ReportIssue>();
    let cursor: string | null = null;
    do {
      const data: ReportIssuesData = await this.graphql<ReportIssuesData>(REPORT_ISSUES_QUERY, {q, cursor});
      const search = data.search;
      for (const node of search.nodes) {
        const id = 'body' in node ? REPORT_MARKER.exec(node.body)?.[1] : undefined;
        if (!id || byReport.get(id)?.open)
          continue;
        byReport.set(id, {key: `${node.repository.name}#${node.number}`, url: node.url, open: node.state === 'OPEN'});
      }
      cursor = search.pageInfo.hasNextPage ? search.pageInfo.endCursor : null;
    } while (cursor);
    return byReport;
  }

  async findIssues(refs: IssueRef[]): Promise<Map<string, GitHubIssue>> {
    const missing = new Map(refs.filter((r) => !this.cache.has(r.key)).map((r) => [r.key, r]));
    const found = await this.lookUp([...missing.values()]);
    for (const key of missing.keys())
      this.cache.set(key, found.get(key) ?? null);
    const result = new Map<string, GitHubIssue>();
    for (const ref of refs) {
      const issue = found.get(ref.key) ?? this.cache.get(ref.key);
      if (issue)
        result.set(ref.key, issue);
    }
    return result;
  }

  async lookUp(refs: IssueRef[]): Promise<Map<string, GitHubIssue>> {
    const numbered = refs.flatMap((r) => r.kind === 'github' ? [r] : []);
    const jiraKeys = refs.flatMap((r) => r.kind === 'jira' ? [r.key] : []);
    const [byNumber, byJiraKey] = await Promise.all([this.issuesByNumber(numbered), this.issuesByJiraKey(jiraKeys)]);
    return new Map([...byNumber, ...byJiraKey]);
  }

  async issuesByNumber(refs: GitHubIssueRef[]): Promise<Map<string, GitHubIssue>> {
    const found = new Map<string, GitHubIssue>();
    for (const chunk of chunks(refs, ISSUES_PER_LOOKUP)) {
      const aliases = chunk.map((r, i) => `i${i}: repository(owner: "${GITHUB_ORG}", name: "${r.repo}") ` +
        `{ issue(number: ${r.number}) { ${ISSUE_FIELDS} } }`);
      const data = await this.graphql<LookupData>(`query($withBody: Boolean!) { ${aliases.join('\n')} }`,
        {withBody: false}, (e) => e.type === 'NOT_FOUND');
      chunk.forEach((ref, i) => {
        const node = data[`i${i}`]?.issue;
        if (node)
          found.set(ref.key, toIssue(node));
      });
    }
    return found;
  }

  async issuesByJiraKey(keys: string[]): Promise<Map<string, GitHubIssue>> {
    const found = new Map<string, GitHubIssue>();
    for (const chunk of chunks(keys, KEYS_PER_SEARCH)) {
      const phrases = chunk.map((k) => `"Moved from ${k}"`).join(' OR ');
      for (const node of await this.search(`${repoScope()} is:issue in:body ${phrases}`, true)) {
        const key = MOVED_FROM.exec(node.body ?? '')?.[1];
        if (key && chunk.includes(key))
          found.set(key, toIssue(node));
      }
    }
    return found;
  }

  async search(q: string, withBody: boolean): Promise<IssueNode[]> {
    const nodes: IssueNode[] = [];
    let cursor: string | null = null;
    do {
      const {search}: SearchData = await this.graphql<SearchData>(SEARCH_QUERY, {q, cursor, withBody});
      if (search.issueCount > SEARCH_RESULTS_LIMIT)
        throw new Error(`GitHub search matches ${search.issueCount} issues, more than it returns: ${q}`);
      nodes.push(...search.nodes.filter(isIssueNode));
      cursor = search.pageInfo.hasNextPage ? search.pageInfo.endCursor : null;
    } while (cursor);
    return nodes;
  }

  async graphql<T>(query: string, variables: Record<string, unknown>,
    tolerated: (e: GraphQLError) => boolean = () => false): Promise<T> {
    const resp = await grok.dapi.fetchProxy(GRAPHQL_URL, {
      method: 'POST',
      headers: {'Authorization': `Bearer ${await this.readToken()}`, 'Content-Type': 'application/json'},
      body: JSON.stringify({query, variables}),
    });
    const text = await resp.text();
    if (resp.status === 401)
      this.token = null;
    if (!resp.ok)
      throw new Error(`GitHub API: HTTP ${resp.status}, ${errorMessage(text)}`);
    const json: GraphQLResponse<T> = JSON.parse(text);
    const errors = (json.errors ?? []).filter((e) => !tolerated(e));
    if (errors.length > 0 || !json.data)
      throw new Error(`GitHub API: ${errors.map((e) => e.message).join('; ') || 'empty response'}`);
    return json.data;
  }

  readToken(): Promise<string> {
    if (!this.token) {
      this.token = this.pkg.getCredentials().then((credentials: DG.Credentials | null) => {
        const token: unknown = credentials?.parameters[TOKEN_CREDENTIAL];
        if (typeof token !== 'string' || token === '') {
          throw new Error(
            `GitHub token is not set: add "${TOKEN_CREDENTIAL}" to the ${this.pkg.name} package credentials`);
        }
        return token;
      });
      this.token.catch(() => this.token = null);
    }
    return this.token;
  }
}

function errorMessage(body: string): string {
  try {
    const message: unknown = JSON.parse(body).message;
    if (typeof message === 'string')
      return message;
  } catch {}
  return body.slice(0, 200);
}

function repoScope(): string {
  return ISSUE_REPOS.map((r) => `repo:${GITHUB_ORG}/${r}`).join(' ');
}

function chunks<T>(items: T[], size: number): T[][] {
  return Array.from({length: Math.ceil(items.length / size)}, (_, i) => items.slice(i * size, (i + 1) * size));
}

function isIssueNode(node: IssueNode | Record<string, never>): node is IssueNode {
  return 'number' in node;
}

function toIssue(node: IssueNode): GitHubIssue {
  const repo = node.repository.name;
  if (!isIssueRepo(repo))
    throw new Error(`GitHub returned an issue of an unexpected repository: ${repo}`);
  const fields = new Map(node.issueFieldValues.nodes
    .flatMap((v): [string, IssueFieldValueNode][] => v.field ? [[v.field.name, v]] : []));
  const priority = fields.get('Priority')?.name;
  return {
    key: `${repo}#${node.number}`,
    repo,
    number: node.number,
    url: node.url,
    title: node.title.trim(),
    resolution: resolution(node),
    boardStatus: node.projectItems.nodes.find((p) => p.project.number === BOARD_PROJECT)?.status?.name ?? null,
    priority: ISSUE_PRIORITIES.find((p) => p === priority) ?? null,
    clients: (fields.get('Client')?.options ?? []).map((o) => o.name),
    type: node.issueType?.name ?? null,
    labels: node.labels.nodes.map((l) => l.name),
    assignees: node.assignees.nodes.map((a) => a.login),
    milestone: node.milestone?.title ?? null,
  };
}

function resolution(node: IssueNode): IssueResolution {
  if (node.state === 'OPEN')
    return 'Unresolved';
  return node.stateReason === 'COMPLETED' ? 'Done' : 'Won\'t fix';
}
