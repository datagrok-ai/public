/// The maintainer reports (build-plan.md WO-8, conventions.md §10, §11.2): what has no owner, what a
/// home still cites that is gone, which features have no test or document, which folders want a feature
/// of their own, and what a branch touches. They read the JSONL, the manifest and the build's own
/// reports — never the index — so they answer wherever a build ran, with or without `kuzu`. The one
/// exception is `replay`, which scores the tiers of `tests-for` against the history and so needs the index.
import * as fs from 'fs';
import * as path from 'path';
import {spawnSync} from 'child_process';
import {TypeSystem, Issue, isSubtype} from './types';
import {Row, compare} from './normalize';
import {helpPage, kebab} from './ids';
import {Graph} from './build/emitter';
import {readJsonl, dataFile, readManifest} from './generation';
import {Section, coverageNote} from './answer';
import {KuzuConnection, run, quote} from './kuzu';
import {resolveTargets, linkedTests} from './ops';

export type ReportName = 'orphans' | 'stale' | 'coverage' | 'proposed' | 'media' | 'diff' | 'replay';
export type ReportFormat = 'table' | 'json' | 'md';
export type Repo = 'core' | 'public' | 'both';

export const REPORT_NAMES: ReportName[] = ['orphans', 'stale', 'coverage', 'proposed', 'media', 'diff', 'replay'];
/** The five `build` writes; `diff` needs a base revision and `replay` the index, so only the verb can produce them. */
export const BUILD_REPORTS: ReportName[] = ['orphans', 'stale', 'coverage', 'proposed', 'media'];
const MEDIA_ROWS = 50;
export const REPOS: Repo[] = ['core', 'public', 'both'];

/** Where a folder earns a feature of its own, and the id shape that would be proposed for it. */
const PROPOSED_PARENTS = ['public/packages', 'public/libraries', 'core/client/d4/lib/src/viewers', 'core/server/datlas/lib/src/services'];
/** Check codes that mean a home still points at something that is gone (homes.ts `checkCitations`). */
const STALE_CODES = ['missing-cited-path', 'missing-doc-link', 'citation-escape', 'bad-anchor', 'missing-path'];
/** A `help-url` that names a Datagrok help page (samples.ts `HELP_URL`); anything else links outside the repo. */
const HELP_LINK = /^(?:https?:\/\/(?:[\w-]+\.)*datagrok\.ai)?\/help\//;
const ORPHAN_GROUPS = 50;
const ORPHAN_FILES = 5;
/** A test file by its path, for the files the index holds no test of (change-tests/plan.md § Validation); a fixture
 * under a test folder is neither a test nor a source. */
const TEST_PATH = /(?:^|\/)(?:tests?|__tests__)\/|\.test\.ts$|_test\.dart$/;
const FIXTURE_PATH = /\/fixtures?\//;
const SOURCE_PATH = /\.(?:ts|dart)$/;
const REPLAY_SUBJECT = 60;
const REPLAY_FILES = 3;

export interface Report {
  name: ReportName;
  title: string;
  /** One line: what the report counted. */
  summary: string;
  /** What it could not see, the Dart pass first. */
  notes: string[];
  sections: Section[];
}

export interface Ownership {
  ambiguous: {file: string, features: string[], rung: number}[];
  orphans: {file: string, loc: number}[];
  resolved_by_chain: {file: string, owner: string, over: string[]}[];
  /** What the build actually saw, the denominator every ownership number is a fraction of (membership.ts). */
  inventory?: Record<string, number>;
}

/** The slice of a build the reports read: rows by id and by type, edge rows by file name, and the build's own records. */
export interface GraphData {
  nodes: Map<string, Row>;
  byType: Map<string, Row[]>;
  edges: Map<string, Row[]>;
  ownership?: Ownership;
  /** `reports/problems.json`: free text per problem kind. */
  problems: Record<string, string[]>;
  /** `reports/home-issues.json`: the check issues behind a partial homes source, structured. */
  homeIssues: Issue[];
  sources: Record<string, string>;
  /** The revisions the graph was built from; `diff` compares them with the working tree. */
  revisions: Record<string, string>;
  repoRoot: string;
}

export interface ReportOptions {
  system: TypeSystem;
  repoRoot: string;
  /** `report diff` only: the revision to compare HEAD against. */
  base?: string;
}

/** Only what [name] reads: the declaration rows alone are 30 MB of the monorepo's build. */
function needs(name: ReportName, system: TypeSystem): {nodes: string[], edges: string[]} {
  const concrete = (test: (type: string) => boolean) => [...system.nodes.values()].filter((t) => !t.abstract && test(t.name)).map((t) => t.name);
  const features = concrete((t) => isSubtype(system, t, 'feature'));
  switch (name) {
    case 'orphans':
      return {nodes: ['source-file', 'package', 'library'], edges: ['is-implemented-in', 'participates-in']};
    case 'stale':
      return {nodes: ['ticket', 'declaration', 'doc-page', ...concrete((t) => !!system.nodes.get(t)!.members.help_url)], edges: ['tracked-in', 'defines-concept']};
    case 'coverage':
      return {nodes: [...features, 'test'], edges: ['tests', 'covers', 'automates', 'documents', 'owner', 'part-of']};
    case 'proposed':
      return {nodes: ['source-file'], edges: ['is-implemented-in']};
    case 'media':
      return {nodes: ['media'], edges: ['embeds', 'thumbnail']};
    default:
      return {nodes: [...features, 'test'], edges: ['is-implemented-in', 'participates-in', 'documents', 'tests', 'covers', 'automates', 'owner']};
  }
}

/** The build output under [kgDir], loaded for one report. */
export function readGraph(kgDir: string, repoRoot: string, system: TypeSystem, name: ReportName): GraphData {
  const want = needs(name, system);
  const data = empty(repoRoot);
  const manifest = readManifest(kgDir);
  data.sources = manifest?.sources ?? {};
  data.revisions = manifest?.revisions ?? {};
  data.ownership = readJson<Ownership>(path.join(kgDir, 'reports', 'ownership.json'));
  data.problems = readJson<Record<string, string[]>>(path.join(kgDir, 'reports', 'problems.json')) ?? {};
  data.homeIssues = readJson<Issue[]>(path.join(kgDir, 'reports', 'home-issues.json')) ?? [];
  for (const type of want.nodes)
    for (const row of readJsonl(dataFile(kgDir, 'nodes', type))) {
      data.nodes.set(String(row.id), row);
      push(data.byType, type, row);
    }
  for (const edge of want.edges) data.edges.set(edge, [...readJsonl(dataFile(kgDir, 'edges', edge))] as Row[]);
  return data;
}

/** The same slice straight from the build that produced it, so `build` writes the reports without re-reading 100 MB. */
export function fromGraph(graph: Graph, repoRoot: string, sources: Record<string, string>, revisions: Record<string, string> = {}): GraphData {
  const data = empty(repoRoot);
  data.sources = sources;
  data.revisions = revisions;
  data.ownership = graph.reports.ownership as Ownership | undefined;
  data.problems = graph.details;
  data.homeIssues = graph.reports['home-issues'] as Issue[] ?? [];
  for (const row of graph.nodes) {
    data.nodes.set(String(row.id), row);
    push(data.byType, String(row.type), row);
  }
  for (const row of graph.edges) push(data.edges, row.type === 'ref' ? String(row.name) : String(row.type), row);
  return data;
}

export function makeReport(name: ReportName, data: GraphData, options: ReportOptions): Report {
  const report = name === 'orphans' ? orphans(data) : name === 'stale' ? stale(data, options)
    : name === 'coverage' ? coverage(data, options) : name === 'proposed' ? proposed(data) : name === 'media' ? media(data) : diff(data, options);
  const note = coverageNote(data.sources);
  if (note) report.notes.unshift(note);
  return report;
}

/** `reports/<name>.json` and `reports/<name>.md` for everything but `diff`; returns what it wrote. */
export function writeReports(kgDir: string, data: GraphData, options: ReportOptions): string[] {
  const dir = path.join(kgDir, 'reports');
  fs.mkdirSync(dir, {recursive: true});
  const written: string[] = [];
  for (const name of BUILD_REPORTS) {
    const report = makeReport(name, data, options);
    fs.writeFileSync(path.join(dir, `${name}.json`), `${JSON.stringify(report, null, 2)}\n`);
    fs.writeFileSync(path.join(dir, `${name}.md`), markdown(report));
    written.push(`${name}.json`, `${name}.md`);
  }
  return written;
}


/**
 * Files no feature owns, grouped by the package or the core sub-project they sit in. A group's orphan count means
 * nothing on its own, so each row carries the files that are owned and the files that only take part, and the
 * summary carries the inventory this build observed — the denominator the numbers are a fraction of (review 3 #10).
 * The `owner` column is the group's own owner, from the `package.json` author: no feature is not no owner.
 */
function orphans(data: GraphData): Report {
  const owned = new Set((data.edges.get('is-implemented-in') ?? []).map((e) => String(e.to)));
  const participating = new Set((data.edges.get('participates-in') ?? []).map((e) => String(e.from)));
  const groups = new Map<string, {orphans: {file: string, loc: number}[], owned: number, participating: number, files: number}>();
  const group = (file: string) => {
    const key = groupOf(file);
    let entry = groups.get(key);
    if (!entry) groups.set(key, entry = {orphans: [], owned: 0, participating: 0, files: 0});
    return entry;
  };
  for (const row of data.byType.get('source-file') ?? []) {
    const entry = group(String(row.path ?? ''));
    entry.files++;
    if (owned.has(String(row.id))) entry.owned++;
    else if (participating.has(String(row.id))) entry.participating++;
  }
  for (const orphan of data.ownership?.orphans ?? []) {
    const entry = group(orphan.file);
    entry.orphans.push(orphan);
  }
  const owners = new Map([...data.byType.get('package') ?? [], ...data.byType.get('library') ?? []]
    .map((row) => [String(row.path), String(row.owner ?? '')] as [string, string]));
  const rows = [...groups].map(([name, e]) => ({
    group: name, owner: owners.get(name) ?? '', files: e.files, owned: e.owned, participating: e.participating, orphans: e.orphans.length,
    loc: e.orphans.reduce((sum, f) => sum + f.loc, 0),
    largest: [...e.orphans].sort((a, b) => b.loc - a.loc || compare(a.file, b.file)).slice(0, ORPHAN_FILES)
      .map((f) => `${f.file.slice(name.length + 1)} (${f.loc})`).join(', '),
  })).filter((r) => r.orphans).sort((a, b) => b.loc - a.loc || compare(a.group, b.group));
  const loc = rows.reduce((sum, r) => sum + r.loc, 0);
  const inventory = data.ownership?.inventory;
  const total = data.ownership?.orphans.length ?? 0;
  const of = inventory ? ` of ${inventory.observed_files} observed (${inventory.owned_files} owned, ${inventory.participating_files} participating)` : '';
  const lines = inventory ? ` of ${inventory.observed_loc}` : '';
  return {
    name: 'orphans', title: 'Orphan files',
    summary: `${total} files${of} have no owner, ${loc}${lines} lines, in ${rows.length} groups; the ${Math.min(ORPHAN_GROUPS, rows.length)} largest groups below.`,
    notes: [
      ...(data.ownership ? [] : ['no reports/ownership.json: run grok kg build']),
      ...(inventory ? [] : ['no inventory in reports/ownership.json: the counts have no denominator']),
      'the inventory covers the files the extractors observed, not every file in the repositories',
    ],
    sections: [{title: 'groups', rows: rows.slice(0, ORPHAN_GROUPS)}],
  };
}

/**
 * The media backlog (conventions.md §5.7): what pages show that nothing describes, descriptions the file has outgrown,
 * proposals a person has not reviewed, files no page shows, embeds of files that are gone, embeds with no alt text,
 * one file committed under several paths, the largest files, and what is rated unfit and still shown.
 */
function media(data: GraphData): Report {
  const items = (data.byType.get('media') ?? []).filter((m) => m.status !== 'proposed');
  const shown = new Map<string, Row[]>();
  for (const e of data.edges.get('embeds') ?? []) {
    const list = shown.get(String(e.to));
    if (list) list.push(e);
    else shown.set(String(e.to), [e]);
  }
  const thumbnails = new Set((data.edges.get('thumbnail') ?? []).map((e) => String(e.to)));
  const pages = (m: Row) => new Set((shown.get(String(m.id)) ?? []).map((e) => String(e.from))).size;
  const bytes = (m: Row) => Number(m.bytes ?? 0);
  const file = (m: Row) => ({id: m.id, format: m.format, bytes: bytes(m), pages: pages(m)});
  const byPages = (a: Row, b: Row) => pages(b) - pages(a) || bytes(b) - bytes(a) || compare(String(a.id), String(b.id));
  const byBytes = (a: Row, b: Row) => bytes(b) - bytes(a) || compare(String(a.id), String(b.id));
  const undescribed = items.filter((m) => m.description === undefined && shown.has(String(m.id))).sort(byPages);
  const stale = items.filter((m) => typeof m.described_blob === 'string' && m.blob !== undefined && m.described_blob !== m.blob).sort(byPages);
  const proposals = items.filter((m) => m.description !== undefined && m.reviewed !== true).sort(byPages);
  const unreferenced = items.filter((m) => m.path !== undefined && !shown.has(String(m.id)) && !thumbnails.has(String(m.id))).sort(byBytes);
  const noAlt = [...shown.values()].flat().filter((e) => e.form === 'image' && !e.alt).map((e) => ({page: e.from, line: e.line, media: e.to}))
    .sort((a, b) => compare(String(a.page), String(b.page)) || Number(a.line) - Number(b.line));
  const blobs = new Map<string, Row[]>();
  for (const m of items)
    if (typeof m.blob === 'string') {
      const list = blobs.get(m.blob);
      if (list) list.push(m);
      else blobs.set(m.blob, [m]);
    }
  const duplicates = [...blobs.values()].filter((l) => l.length > 1).map((l) => ({blob: l[0].blob, bytes: bytes(l[0]), copies: l.length,
    paths: l.map((m) => String(m.path)).sort(compare).join(', ')})).sort((a, b) => b.bytes * b.copies - a.bytes * a.copies || compare(String(a.blob), String(b.blob)));
  const unfit = items.filter((m) => m.quality === 'unfit' && shown.has(String(m.id))).sort(byPages);
  const total = items.reduce((sum, m) => sum + bytes(m), 0);
  const section = (title: string, rows: Record<string, unknown>[]) => ({title, rows: rows.slice(0, MEDIA_ROWS), total: rows.length});
  return {
    name: 'media', title: 'Media backlog',
    summary: `${items.length} media (${Math.round(total / 1048576)} MB): ${undescribed.length} shown but undescribed, ${stale.length} described before the file changed, ` +
      `${proposals.length} awaiting review, ${unreferenced.length} shown nowhere, ${(data.problems.broken_embeds ?? []).length} broken embeds, ${noAlt.length} images without alt text, ` +
      `${duplicates.length} files committed more than once, ${unfit.length} rated unfit and still shown.`,
    notes: [
      ...((data.problems.untracked_media ?? []).length ? [`${data.problems.untracked_media.length} media files under the roots are not tracked by git; they are listed, not indexed`] : []),
      'a record beside the file (media.yaml) carries the description; grok kg enrich media proposes one for what is undescribed or stale',
    ],
    sections: [
      section('undescribed', undescribed.map(file)),
      section('stale', stale.map((m) => ({...file(m), described_blob: String(m.described_blob).slice(0, 12), blob: String(m.blob).slice(0, 12)}))),
      section('awaiting review', proposals.map((m) => ({...file(m), described_by: m.described_by, quality: m.quality}))),
      section('unreferenced', unreferenced.map((m) => ({id: m.id, format: m.format, bytes: bytes(m), path: m.path}))),
      section('broken embeds', (data.problems.broken_embeds ?? []).map((line) => ({embed: line}))),
      section('images without alt text', noAlt),
      section('duplicates', duplicates),
      section('largest', [...items].sort(byBytes).slice(0, MEDIA_ROWS).map(file)),
      section('unfit but shown', unfit.map((m) => ({...file(m), quality_notes: m.quality_notes}))),
      section('untracked', (data.problems.untracked_media ?? []).map((path) => ({path}))),
    ],
  };
}

/** `public/packages/Chem`, `public/libraries/utils`, `core/client/d4`, `public/js-api`: the first two segments of a repo path. */
function groupOf(file: string): string {
  const segments = file.split('/');
  const deep = segments[0] === 'core' || ['public/packages', 'public/libraries'].includes(segments.slice(0, 2).join('/'));
  return segments.slice(0, deep ? 3 : 2).join('/');
}

/** What the graph still points at and cannot reach: a gone path, an unknown ticket, a help page, a spec, a declaration. */
function stale(data: GraphData, options: ReportOptions): Report {
  const rows: Record<string, unknown>[] = [];
  for (const issue of data.homeIssues) {
    if (!STALE_CODES.includes(issue.code)) continue;
    rows.push({kind: 'citation', source: issue.line ? `${issue.file}:${issue.line}` : issue.file, target: issue.target ?? '',
      reason: `${issue.code}: ${issue.message}`});
  }
  // without the snapshot a ticket with no external provenance proves nothing: the backlog was never consulted
  const backlog = /^ok\b/.test(data.sources.backlog ?? '');
  for (const edge of data.edges.get('tracked-in') ?? []) {
    const ticket = data.nodes.get(String(edge.to));
    if (ticket && ticket.provenance === 'external') continue;
    rows.push({kind: 'ticket', source: String(edge.from), target: String(edge.to),
      reason: backlog ? 'not in the backlog snapshot' : `unknown: the backlog snapshot is ${data.sources.backlog ?? 'absent'}`});
  }
  for (const row of data.nodes.values()) {
    if (typeof row.help_url !== 'string' || !HELP_LINK.test(row.help_url)) continue;
    const page = helpPage(options.repoRoot, row.help_url);
    if (page && data.nodes.has(`doc:${page}`)) continue;
    rows.push({kind: 'help-url', source: String(row.id), target: row.help_url, reason: page ? `${page} is not a doc-page` : 'names no page under public/help'});
  }
  for (const problem of data.problems.unresolved_ids ?? []) {
    const realized = /^(.+?): realized_as (\S+) (.+)$/.exec(problem);
    if (realized) rows.push({kind: 'realized-as', source: realized[1], target: realized[2], reason: realized[3]});
    const help = /^(.+?): help-url (\S+) (.+)$/.exec(problem);
    if (help && HELP_LINK.test(help[2])) rows.push({kind: 'help-url', source: help[1], target: help[2], reason: help[3]});
  }
  if (data.sources.dart === 'ok')
    for (const edge of data.edges.get('defines-concept') ?? []) {
      const decl = data.nodes.get(String(edge.from));
      if (!decl || decl.provenance !== 'annotation' || !String(decl.path ?? '').endsWith('.dart')) continue;
      rows.push({kind: 'declaration', source: String(edge.to), target: String(edge.from), reason: 'defined_by names a declaration the Dart pass does not have'});
    }
  const kinds = new Map<string, number>();
  for (const row of rows) kinds.set(String(row.kind), (kinds.get(String(row.kind)) ?? 0) + 1);
  return {
    name: 'stale', title: 'Stale references',
    summary: `${rows.length} stale references${kinds.size ? `: ${[...kinds].sort(([a], [b]) => compare(a, b)).map(([k, n]) => `${n} ${k}`).join(', ')}` : ''}.`,
    notes: [
      ...(backlog ? [] : ['the backlog snapshot was not read, so a tracked-in ticket is reported as unknown, not as absent']),
      ...(data.sources.dart === 'ok' ? [] : ['the Dart pass is not ok, so defined_by declarations are not checked']),
    ],
    sections: [{title: 'references', rows: rows.sort((a, b) => compare(String(a.kind), String(b.kind)) || compare(String(a.source), String(b.source)) || compare(String(a.target), String(b.target)))}],
  };
}

/**
 * One row per feature: who owns it, what really runs against it, what documents it, and whether it says what it is.
 * A test count that mixes a skipped test, a test whose name is built at run time and a test that runs is not a
 * coverage number (review 3 #10), so each of them is its own column, and so is the coverage a feature only inherits
 * from the features under it. A stub — a feature with no home of its own — is marked rather than counted as a gap.
 */
function coverage(data: GraphData, options: ReportOptions): Report {
  const features = [...data.nodes.values()].filter((r) => isSubtype(options.system, String(r.type), 'feature'));
  const documents = count(data.edges.get('documents'), 'to');
  const owners = new Map((data.edges.get('owner') ?? []).map((e) => [String(e.from), String(e.to)]));
  const scenarios = new Map<string, string[]>();
  for (const edge of data.edges.get('covers') ?? []) push(scenarios, String(edge.to), String(edge.from));
  const automated = count(data.edges.get('automates'), 'to');
  const children = new Map<string, string[]>();
  for (const edge of data.edges.get('part-of') ?? []) push(children, String(edge.to), String(edge.from));
  const direct = new Map<string, string[]>();
  for (const edge of data.edges.get('tests') ?? []) push(direct, String(edge.to), String(edge.from));
  const split = (tests: string[]) => {
    const out = {runnable: 0, skipped: 0, dynamic: 0};
    for (const id of tests) {
      const test = data.nodes.get(id);
      if (test?.skipped === true) out.skipped++;
      else if (test?.dynamic === true) out.dynamic++;
      else out.runnable++;
    }
    return out;
  };
  const rows = features.map((f) => {
    const id = String(f.id);
    const own = direct.get(id) ?? [];
    const counts = split(own);
    const covered = scenarios.get(id) ?? [];
    return {
      feature: id, owner: (f.owner as string | undefined) ?? owners.get(id) ?? '', status: (f.status as string | undefined) ?? '',
      stub: !f.home, tests_runnable: counts.runnable, tests_skipped: counts.skipped, tests_dynamic: counts.dynamic,
      inherited: descendants(id, children).reduce((sum, d) => sum + (direct.get(d)?.length ?? 0), 0),
      scenarios: covered.length, scenario_automated: covered.reduce((sum, s) => sum + (automated.get(s) ?? 0), 0),
      user_help: (f.user_help as string | undefined) ?? '', developer_help: (f.developer_help as string | undefined) ?? '',
      no_tests: !own.length && !covered.length, no_docs: !documents.get(id) && !f.user_help && !f.developer_help, no_description: !f.description,
    };
  }).sort((a, b) => compare(a.feature, b.feature));
  const gaps = (key: 'no_tests' | 'no_docs' | 'no_description') => rows.filter((r) => !r.stub && r[key]).length;
  const authored = rows.filter((r) => !r.stub).length;
  return {
    name: 'coverage', title: 'Feature coverage',
    summary: `${rows.length} features, ${rows.length - authored} of them stubs; of the ${authored} with a home, ${gaps('no_tests')} have no test or scenario, ` +
      `${gaps('no_docs')} no document beside the home, ${gaps('no_description')} no first paragraph.`,
    notes: ['tests_runnable excludes a skipped test and one whose name is built at run time; inherited counts the tests of the features under this one'],
    sections: [{title: 'features', rows}],
  };
}

/** Everything under a feature in the part-of tree, the feature itself excluded. */
function descendants(id: string, children: Map<string, string[]>): string[] {
  const out: string[] = [];
  const queue = [...children.get(id) ?? []];
  while (queue.length) {
    const next = queue.shift()!;
    if (next === id || out.includes(next)) continue;
    out.push(next);
    queue.push(...children.get(next) ?? []);
  }
  return out;
}

/**
 * Folders that carry code no feature claims: candidates for a home of their own (conventions.md §10). One owned file
 * used to hide a whole folder, so a partly claimed area never appeared; the ranking is the unowned remainder instead
 * (review 3 #10), and a folder already fully owned is the only one that drops out.
 */
function proposed(data: GraphData): Report {
  const owned = new Set((data.edges.get('is-implemented-in') ?? []).map((e) => String(e.to)));
  const folders = new Map<string, {files: number, loc: number, unowned: number, unowned_loc: number}>();
  for (const row of data.byType.get('source-file') ?? []) {
    const folder = folderOf(String(row.path ?? ''));
    if (!folder) continue;
    const entry = folders.get(folder) ?? {files: 0, loc: 0, unowned: 0, unowned_loc: 0};
    const loc = Number(row.loc ?? 0);
    entry.files++;
    entry.loc += loc;
    if (!owned.has(String(row.id))) {
      entry.unowned++;
      entry.unowned_loc += loc;
    }
    folders.set(folder, entry);
  }
  const rows = [...folders].filter(([, e]) => e.unowned)
    .map(([folder, e]) => ({path: folder, files: e.files, loc: e.loc, unowned_files: e.unowned, unowned_loc: e.unowned_loc, suggested: suggestedId(folder)}))
    .sort((a, b) => b.unowned_loc - a.unowned_loc || compare(a.path, b.path));
  const partial = rows.filter((r) => r.unowned_files < r.files).length;
  return {
    name: 'proposed', title: 'Proposed features',
    summary: `${rows.length} folders carry code no feature owns, of ${folders.size} under ${PROPOSED_PARENTS.join(', ')}; ` +
      `${partial} of them ${partial === 1 ? 'is' : 'are'} partly owned already; ${rows.filter((r) => r.suggested).length} with a suggested id.`,
    notes: ['a folder is only seen through the source files the build extracted; the ranking is the unowned remainder'],
    sections: [{title: 'folders', rows}],
  };
}

/** The folder directly under one of the proposed parents that [file] sits in; a file lying in the parent itself has none. */
function folderOf(file: string): string | undefined {
  for (const parent of PROPOSED_PARENTS) {
    if (!file.startsWith(`${parent}/`)) continue;
    const rest = file.slice(parent.length + 1).split('/');
    return rest.length > 1 ? `${parent}/${rest[0]}` : undefined;
  }
  return undefined;
}

function suggestedId(folder: string): string {
  const name = folder.split('/').pop()!;
  if (folder.startsWith('public/packages/')) return `domains/${kebab(name)}`;
  if (folder.startsWith('core/client/d4/lib/src/viewers/')) return `visualize/viewers/${kebab(name)}`;
  return '';
}

/** How a changed file reaches a feature, strongest first. */
const RELATIONS = ['owns', 'home', 'documents', 'participates'];

interface Hit {
  relation: string;
  confidence: number;
}

function stronger(a: Hit, b: Hit): boolean {
  const rank = (h: Hit) => RELATIONS.indexOf(h.relation);
  return rank(a) !== rank(b) ? rank(a) < rank(b) : a.confidence > b.confidence;
}

/** What a branch touches: the features behind the changed files, and everything that tests them. */
function diff(data: GraphData, options: ReportOptions): Report {
  const {files, problems} = changedFiles(options.repoRoot, options.base!);
  const owns = new Map<string, Row[]>();
  for (const edge of data.edges.get('is-implemented-in') ?? []) push(owns, String(edge.to), edge);
  const parts = new Map<string, Row[]>();
  for (const edge of data.edges.get('participates-in') ?? []) push(parts, String(edge.from), edge);
  const documents = new Map<string, Row[]>();
  for (const edge of data.edges.get('documents') ?? []) push(documents, String(edge.from), edge);
  const homes = new Map<string, string>();
  for (const row of data.nodes.values())
    if (typeof row.home === 'string' && isSubtype(options.system, String(row.type), 'feature')) homes.set(row.home, String(row.id));

  const fileRows: Record<string, unknown>[] = [];
  const features = new Map<string, Hit>();
  let unowned = 0;
  let orphanedDeletions = 0;
  for (const {file, change} of files) {
    const hits = new Map<string, Hit>();
    const touch = (feature: string, relation: string, confidence: number) => {
      const best = hits.get(feature);
      if (!best || stronger({relation, confidence}, best)) hits.set(feature, {relation, confidence});
    };
    for (const edge of owns.get(`file:${file}`) ?? []) touch(String(edge.from), 'owns', Number(edge.confidence ?? 1));
    for (const edge of parts.get(`file:${file}`) ?? []) touch(String(edge.to), 'participates', Number(edge.confidence ?? 1));
    // a home document has no file: node of its own, and a help page reaches its feature through documents
    const home = homes.get(file);
    if (home) touch(home, 'home', 1);
    for (const edge of documents.get(`doc:${file}`) ?? []) touch(String(edge.to), 'documents', Number(edge.confidence ?? 1));
    if (!hits.size) {
      // a deleted file is gone from the working tree: the graph is the only place its owner can still be read
      if (change === 'deleted') orphanedDeletions++;
      else unowned++;
      fileRows.push({file, feature: '', relation: change === 'deleted' ? 'deleted' : '', confidence: '', change});
      continue;
    }
    for (const [feature, hit] of [...hits].sort(([a], [b]) => compare(a, b))) {
      fileRows.push({file, feature, ...hit, change});
      const best = features.get(feature);
      if (!best || stronger(hit, best)) features.set(feature, hit);
    }
  }
  const ids = new Set(features.keys());
  const owners = new Map((data.edges.get('owner') ?? []).map((e) => [String(e.from), String(e.to)]));
  const featureRows = [...features].sort(([a], [b]) => compare(a, b)).map(([feature, hit]) => ({
    feature, name: String(data.nodes.get(feature)?.name ?? ''), owner: (data.nodes.get(feature)?.owner as string | undefined) ?? owners.get(feature) ?? '',
    relation: hit.relation, confidence: hit.confidence,
  }));
  const testRows = (data.edges.get('tests') ?? []).filter((e) => ids.has(String(e.to)))
    .map((e) => ({test: String(e.from), feature: String(e.to), level: String(data.nodes.get(String(e.from))?.level ?? '')})).sort((a, b) => compare(a.test, b.test));
  const scenarioRows = (data.edges.get('covers') ?? []).filter((e) => ids.has(String(e.to)))
    .map((e) => ({scenario: String(e.from), feature: String(e.to)})).sort((a, b) => compare(a.scenario, b.scenario));
  const covered = new Set(scenarioRows.map((r) => r.scenario));
  const automationRows = (data.edges.get('automates') ?? []).filter((e) => covered.has(String(e.to)))
    .map((e) => ({test: String(e.from), scenario: String(e.to)})).sort((a, b) => compare(a.test, b.test));
  const deleted = files.filter((f) => f.change === 'deleted').length;
  return {
    name: 'diff', title: `Features touched by ${options.base}...HEAD`,
    summary: `${files.length} files changed (${deleted} deleted), ${unowned} of them in no feature; ${featureRows.length} features touched; ` +
      `${testRows.length} tests, ${scenarioRows.length} scenarios and ${automationRows.length} automations cover them.`,
    notes: [
      ...problems,
      ...(orphanedDeletions ? [`${orphanedDeletions} deleted file${orphanedDeletions === 1 ? '' : 's'} had no owner in this graph either; they are listed as deleted`] : []),
      ...staleGraph(data, options.repoRoot),
    ],
    sections: [
      {title: 'features', rows: featureRows},
      {title: 'files', rows: fileRows},
      {title: 'tests', rows: testRows},
      {title: 'scenarios', rows: scenarioRows},
      {title: 'automations', rows: automationRows},
    ],
  };
}

interface Change {
  file: string;
  change: 'added' | 'modified' | 'deleted' | 'renamed' | 'changed';
}

/** What `git diff --name-status` calls a change, in the words the report uses. */
const CHANGES: Record<string, Change['change']> = {A: 'added', M: 'modified', D: 'deleted', R: 'renamed', C: 'added', T: 'changed'};

/**
 * The changed files of both repositories, the public ones under the `public/` prefix the graph gives them. A monorepo
 * revision means nothing inside the submodule, so the public baseline is the gitlink that revision recorded
 * (`git rev-parse <base>:public`) and only falls back to the same string when the gitlink cannot be read (review 3
 * #10). Output is NUL-delimited: a path may contain anything but a NUL, quoting included.
 */
function changedFiles(repoRoot: string, base: string): {files: Change[], problems: string[]} {
  const found = new Map<string, Change>();
  const problems: string[] = [];
  const roots = new Set<string>();
  for (const [cwd, prefix] of [[repoRoot, ''], [path.join(repoRoot, 'public'), 'public/']]) {
    // a checkout where public/ is a plain folder rather than the submodule answers the same diff twice
    const top = spawnSync('git', ['-C', cwd, 'rev-parse', '--show-toplevel'], {encoding: 'utf8'});
    if (top.status !== 0 || roots.has(top.stdout.trim())) continue;
    roots.add(top.stdout.trim());
    let from = base;
    if (prefix) {
      const link = gitlink(repoRoot, base);
      if (link) from = link;
      else problems.push(`public/: ${base}:public names no gitlink, so ${base} is used in the submodule as well`);
    }
    const r = spawnSync('git', ['-C', cwd, 'diff', '--name-status', '-z', `${from}...HEAD`], {encoding: 'utf8', maxBuffer: 32 * 1024 * 1024});
    if (r.status !== 0) {
      problems.push(`${prefix || 'the monorepo'}: git diff ${from}...HEAD failed: ${(r.stderr ?? '').trim().split('\n')[0] || `exit ${r.status}`}`);
      continue;
    }
    // NUL-separated fields, a status then its path — and for a rename or a copy, the old path then the new one
    const fields = r.stdout.split('\0').filter((f) => f !== '');
    for (let i = 0; i < fields.length;) {
      const status = fields[i++];
      const paths = /^[RC]/.test(status) ? [fields[i++], fields[i++]] : [fields[i++]];
      const file = paths[paths.length - 1];
      // the submodule's own gitlink is not a file of either repository
      if (file && file !== 'public') found.set(prefix + file, {file: prefix + file, change: CHANGES[status[0]] ?? 'modified'});
    }
  }
  return {files: [...found.values()].sort((a, b) => compare(a.file, b.file)), problems};
}

/** The public commit [base] recorded: mode 160000 is the gitlink, anything else is a folder that only looks like one. */
function gitlink(repoRoot: string, base: string): string | undefined {
  const r = spawnSync('git', ['-C', repoRoot, 'ls-tree', base, '--', 'public'], {encoding: 'utf8'});
  const m = r.status === 0 ? /^160000 commit ([0-9a-f]{7,40})/.exec(r.stdout.trim()) : null;
  return m ? m[1] : undefined;
}

/** A graph built from other commits than the ones checked out answers about files that are no longer there. */
function staleGraph(data: GraphData, repoRoot: string): string[] {
  const notes: string[] = [];
  for (const [name, dir] of [['reddata', repoRoot], ['public', path.join(repoRoot, 'public')]]) {
    const built = data.revisions[name];
    if (!built) continue;
    const head = spawnSync('git', ['-C', dir, 'rev-parse', 'HEAD'], {encoding: 'utf8'});
    if (head.status !== 0 || !head.stdout.trim() || head.stdout.trim() === built) continue;
    notes.push(`${name}: the graph was built at ${built.slice(0, 10)} and the working tree is at ${head.stdout.trim().slice(0, 10)}; ` +
      'ownership may be out of date — run grok kg build');
  }
  return notes;
}

/** One commit of the history `replay` scores: its files under the prefix the graph gives them. */
export interface Commit {
  sha: string;
  subject: string;
  files: string[];
}

/** The last [n] non-merge commits of the monorepo, of `public/`, or of both, newest first, the public paths prefixed. */
export function readLog(repoRoot: string, n: number, repo: Repo): Commit[] {
  const out: Commit[] = [];
  for (const [name, cwd, prefix] of [['core', repoRoot, ''], ['public', path.join(repoRoot, 'public'), 'public/']] as const) {
    if (repo !== 'both' && repo !== name) continue;
    const r = spawnSync('git', ['-C', cwd, '-c', 'core.quotepath=false', 'log', '--no-merges', '--name-only', '--format=%x01%H%x00%s', '-n', String(n)],
      {encoding: 'utf8', maxBuffer: 64 * 1024 * 1024});
    if (r.status !== 0) throw new Error(`${name}: git log failed: ${(r.stderr ?? '').trim().split('\n')[0] || `exit ${r.status}`}`);
    for (const block of r.stdout.split('\x01').slice(1)) {
      const [head, ...lines] = block.split('\n');
      const [sha, subject] = head.split('\0');
      out.push({sha, subject: subject ?? '', files: lines.map((l) => l.trim()).filter((f) => f && f !== 'public').map((f) => prefix + f)});
    }
  }
  return out;
}

/**
 * The yardstick of the notation (change-tests/plan.md § Validation): for each commit that changed a source file and
 * a test file, the immediate and reachable tiers of its source files against the current index, and whether the
 * test files it changed are in them. A test file is one the index holds a test of, or one named like one; a commit
 * none of whose source files the index knows cannot be scored and is counted in the notes instead.
 */
export async function replay(conn: KuzuConnection, commits: Commit[], sources?: Record<string, string>): Promise<Report> {
  const started = Date.now();
  const paths = [...new Set(commits.flatMap((c) => c.files))];
  const frameworks = new Map<string, string[]>();
  const tests = await run(conn, `MATCH (t:Artifact) WHERE t.${quote('type')} = 'test' AND t.${quote('path')} IN $paths ` +
    `RETURN DISTINCT t.${quote('path')} AS path, t.${quote('framework')} AS framework`, {paths});
  for (const r of tests.rows) push(frameworks, String(r.path), String(r.framework ?? ''));
  const isTest = (f: string) => frameworks.has(f) || TEST_PATH.test(f);
  const candidates = commits.map((c) => {
    const files = c.files.filter((f) => SOURCE_PATH.test(f) && !FIXTURE_PATH.test(f));
    return {commit: c, sources: files.filter((f) => !isTest(f)), tests: files.filter(isTest)};
  }).filter((c) => c.sources.length && c.tests.length);
  const resolved = await resolveTargets(conn, [...new Set(candidates.flatMap((c) => c.sources))]);
  const known = (f: string) => String(resolved.get(f)?.id ?? '').startsWith('file:');
  const perFramework = new Map<string, {test_files: number, immediate: number, linked: number}>();
  const misses: Record<string, unknown>[] = [];
  const cache = new Map<string, {immediate: Set<string>, linked: Set<string>}>();
  let considered = 0, hitImmediate = 0, hitLinked = 0, unscored = 0, partial = 0;
  for (const {commit, sources: changed, tests: testFiles} of candidates) {
    const ids = changed.filter(known).map((f) => String(resolved.get(f)!.id)).sort();
    if (!ids.length) {
      unscored++;
      continue;
    }
    if (ids.length < changed.length) partial++;
    considered++;
    let tiers = cache.get(ids.join(' '));
    if (!tiers) {
      const linked = await linkedTests(conn, ids.map((id) => ({id})));
      const immediate = new Set(linked.immediate.map((r) => String(r.path)));
      cache.set(ids.join(' '), tiers = {immediate, linked: new Set([...immediate, ...linked.reachable.map((r) => String(r.path))])});
    }
    const missed = testFiles.filter((f) => !tiers!.immediate.has(f));
    if (!missed.length) hitImmediate++;
    if (testFiles.every((f) => tiers!.linked.has(f))) hitLinked++;
    for (const f of testFiles)
      for (const framework of frameworks.get(f) ?? ['']) {
        const row = perFramework.get(framework) ?? {test_files: 0, immediate: 0, linked: 0};
        row.test_files++;
        if (tiers.immediate.has(f)) row.immediate++;
        if (tiers.linked.has(f)) row.linked++;
        perFramework.set(framework, row);
      }
    if (missed.length)
      misses.push({commit: commit.sha.slice(0, 10), subject: commit.subject.slice(0, REPLAY_SUBJECT), sources: few(changed), missed: few(missed)});
  }
  const rate = (n: number) => considered ? `${Math.round(100 * n / considered)}% (${n}/${considered})` : 'n/a';
  const frameworkRows = [...perFramework].sort(([a], [b]) => compare(a, b)).map(([framework, r]) => ({framework: framework || '(no test in the index)', ...r,
    immediate_rate: `${Math.round(100 * r.immediate / r.test_files)}%`, linked_rate: `${Math.round(100 * r.linked / r.test_files)}%`}));
  const note = coverageNote(sources);
  return {
    name: 'replay', title: `History replay over ${commits.length} commits`,
    summary: `${considered} of ${commits.length} commits changed a source file the index knows and a test file; every changed test file is in the ` +
      `immediate tier for ${rate(hitImmediate)}, in immediate or reachable for ${rate(hitLinked)}; ${((Date.now() - started) / 1000).toFixed(1)} s.`,
    notes: [
      ...(note ? [note] : []),
      `${unscored} commit${unscored === 1 ? '' : 's'} with a source and a test file skipped: none of the source files is in the index (new or deleted ` +
        `since the build); ${partial} of the scored commits have some source files the index does not know`,
    ],
    sections: [
      {title: 'frameworks', rows: frameworkRows},
      {title: 'misses', rows: misses},
    ],
  };
}

/** The first few paths and how many more there are: a commit can touch hundreds. */
function few(files: string[]): string {
  return `${files.slice(0, REPLAY_FILES).join(', ')}${files.length > REPLAY_FILES ? ` (+${files.length - REPLAY_FILES} more)` : ''}`;
}

function empty(repoRoot: string): GraphData {
  return {nodes: new Map(), byType: new Map(), edges: new Map(), problems: {}, homeIssues: [], sources: {}, revisions: {}, repoRoot};
}

function count(rows: Row[] | undefined, key: 'from' | 'to'): Map<string, number> {
  const out = new Map<string, number>();
  for (const row of rows ?? []) out.set(String(row[key]), (out.get(String(row[key])) ?? 0) + 1);
  return out;
}

function push<T>(map: Map<string, T[]>, key: string, value: T): void {
  const list = map.get(key);
  if (list) list.push(value);
  else map.set(key, [value]);
}

/** The PR-comment rendering: a title, the summary, the notes as a quote, and one table per section. */
export function markdown(report: Report): string {
  const out = [`# ${report.title}`, '', report.summary, ''];
  for (const note of report.notes) out.push(`> ${note}`, '');
  for (const section of report.sections) {
    out.push(`## ${section.title} (${section.rows.length})`, '');
    out.push(...table(section.rows), '');
  }
  return `${out.join('\n').replace(/\n+$/, '')}\n`;
}

/** A GitHub table with padded cells, or `_none_` when the section is empty. */
function table(rows: Record<string, unknown>[]): string[] {
  if (!rows.length) return ['_none_'];
  const keys = [...new Set(rows.flatMap((r) => Object.keys(r)))];
  const cells = rows.map((row) => keys.map((k) => cell(row[k])));
  const widths = keys.map((k, i) => Math.max(k.length, ...cells.map((c) => c[i].length)));
  const line = (values: string[]) => `| ${values.map((v, i) => v.padEnd(widths[i])).join(' | ')} |`;
  return [line(keys), `|${widths.map((w) => '-'.repeat(w + 2)).join('|')}|`, ...cells.map(line)];
}

function cell(value: unknown): string {
  if (value === null || value === undefined) return '';
  return String(Array.isArray(value) ? value.join(', ') : value).replace(/\|/g, '\\|').replace(/\s*\n\s*/g, ' ');
}

function readJson<T>(file: string): T | undefined {
  return fs.existsSync(file) ? JSON.parse(fs.readFileSync(file, 'utf8')) as T : undefined;
}

