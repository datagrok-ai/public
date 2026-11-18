export const GITHUB_ORG = 'datagrok-ai';
export const ISSUE_REPOS = ['core', 'public'] as const;
export type IssueRepo = typeof ISSUE_REPOS[number];

export const ISSUE_PRIORITIES = ['Blocker', 'Urgent', 'High', 'Medium', 'Low'] as const;
export type IssuePriority = typeof ISSUE_PRIORITIES[number];

export type IssueResolution = 'Unresolved' | 'Done' | 'Won\'t fix';

export interface GitHubIssue {
  key: string;
  repo: IssueRepo;
  number: number;
  url: string;
  title: string;
  resolution: IssueResolution;
  boardStatus: string | null;
  priority: IssuePriority | null;
  clients: string[];
  type: string | null;
  labels: string[];
  assignees: string[];
  milestone: string | null;
}

export interface GitHubIssueRef {
  kind: 'github';
  key: string;
  repo: IssueRepo;
  number: number;
}

export interface JiraIssueRef {
  kind: 'jira';
  key: string;
}

export type IssueRef = GitHubIssueRef | JiraIssueRef;

export function isOpen(issue: GitHubIssue): boolean {
  return issue.resolution === 'Unresolved';
}

export function isFixed(issue: GitHubIssue): boolean {
  return issue.resolution === 'Done';
}

function priorityRank(priority: IssuePriority | null): number {
  return priority === null ? ISSUE_PRIORITIES.length : ISSUE_PRIORITIES.indexOf(priority);
}

export function compareIssues(a: GitHubIssue, b: GitHubIssue): number {
  return Number(isOpen(b)) - Number(isOpen(a)) ||
    priorityRank(a.priority) - priorityRank(b.priority) ||
    a.repo.localeCompare(b.repo) || a.number - b.number;
}

export function ticketVerdict(refs: IssueRef[], issues: Map<string, GitHubIssue>): string {
  if (refs.length === 0)
    return '';
  const found = refs.flatMap((ref) => issues.get(ref.key) ?? []);
  if (found.length === 0)
    return 'Not on GitHub';
  const unfixed = found.filter((i) => !isFixed(i));
  if (unfixed.length === 0)
    return 'Fixed';
  const priority = ISSUE_PRIORITIES.find((p) => unfixed.some((i) => i.priority === p)) ?? 'No priority';
  return `${unfixed.length < found.length ? 'Partially Fixed' : 'Wasn\'t Fixed'} (${priority})`;
}

export function isIssueRepo(name: string): name is IssueRepo {
  return ISSUE_REPOS.some((r) => r === name);
}

const REPO_PATTERN = ISSUE_REPOS.join('|');
const REF_PATTERN = new RegExp(
  `\\b(GROK-\\d+)\\b|(?:https://github\\.com/)?${GITHUB_ORG}/(${REPO_PATTERN})(?:#|/issues/)(\\d+)` +
  `|\\b(${REPO_PATTERN})#(\\d+)\\b`, 'g');

export function findIssueRefs(text: string | null): IssueRef[] {
  const refs = new Map<string, IssueRef>();
  for (const [, jiraKey, urlRepo, urlNumber, shortRepo, shortNumber] of (text ?? '').matchAll(REF_PATTERN)) {
    const repo = urlRepo ?? shortRepo;
    if (jiraKey)
      refs.set(jiraKey, {kind: 'jira', key: jiraKey});
    else if (isIssueRepo(repo)) {
      const number = Number(urlNumber ?? shortNumber);
      refs.set(`${repo}#${number}`, {kind: 'github', key: `${repo}#${number}`, repo, number});
    }
  }
  return [...refs.values()];
}
