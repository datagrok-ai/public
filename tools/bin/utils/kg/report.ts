/// The maintainer reports (build-plan.md WO-8, conventions.md §10, §11.2): what has no owner, what a
/// home still cites that is gone, which features have no test or document, which folders want a feature
/// of their own, and what a branch touches. They read the JSONL, the manifest and the build's own
/// reports — never the index — so they answer wherever a build ran, with or without `kuzu`.
import * as fs from 'fs';
import * as path from 'path';
import {spawnSync} from 'child_process';
import {TypeSystem, isSubtype} from './types';
import {Row} from './build/normalize';
import {Graph} from './build/emitter';
import {Manifest} from './build/write';
import {coverageNote} from './ops';
import {helpPage} from './build/extract/ts/samples';
import {kebab} from './build/extract/markers';
import {printOutput} from '../server-output';
import type {Section} from './ops';

export type ReportName = 'orphans' | 'stale' | 'coverage' | 'proposed' | 'diff';
export type ReportFormat = 'table' | 'json' | 'md';

export const REPORT_NAMES: ReportName[] = ['orphans', 'stale', 'coverage', 'proposed', 'diff'];
/** The four `build` writes; `diff` needs a base revision, so only the verb can produce it. */
export const BUILD_REPORTS: ReportName[] = ['orphans', 'stale', 'coverage', 'proposed'];

/** Where a folder earns a feature of its own, and the id shape that would be proposed for it. */
const PROPOSED_PARENTS = ['public/packages', 'public/libraries', 'core/client/d4/lib/src/viewers', 'core/server/datlas/lib/src/services'];
/** Check codes that mean a home still points at something that is gone (homes.ts `checkCitations`). */
const STALE_CODES = ['missing-cited-path', 'missing-doc-link', 'citation-escape', 'bad-anchor', 'missing-path'];
/** A `help-url` that names a Datagrok help page (samples.ts `HELP_URL`); anything else links outside the repo. */
const HELP_LINK = /^(?:https?:\/\/(?:[\w-]+\.)*datagrok\.ai)?\/help\//;
const ORPHAN_GROUPS = 50;
const ORPHAN_FILES = 5;
/** `file:line: code: message`, the shape the homes extractor records a check issue in. */
const HOME_ISSUE = /^([^:\s]+?)(?::(\d+))?: ([a-z][a-z-]*): (.*)$/;

export interface Report {
  name: ReportName;
  title: string;
  /** One line: what the report counted. */
  summary: string;
  /** What it could not see, the Dart batch first. */
  notes: string[];
  sections: Section[];
}

export interface Ownership {
  ambiguous: {file: string, features: string[], rung: number}[];
  orphans: {file: string, loc: number}[];
  resolved_by_chain: {file: string, owner: string, over: string[]}[];
}

/** The slice of a build the reports read: rows by id and by type, edge rows by file name, and the build's own records. */
export interface GraphData {
  nodes: Map<string, Row>;
  byType: Map<string, Row[]>;
  edges: Map<string, Row[]>;
  ownership?: Ownership;
  /** `reports/problems.json`: free text per problem kind. */
  problems: Record<string, string[]>;
  sources: Record<string, string>;
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
      return {nodes: [], edges: []};
    case 'stale':
      return {nodes: ['ticket', 'declaration', 'doc-page', ...concrete((t) => !!system.nodes.get(t)!.members.help_url)], edges: ['tracked-in', 'defines-concept']};
    case 'coverage':
      return {nodes: features, edges: ['tests', 'covers', 'automates', 'documents', 'owner']};
    case 'proposed':
      return {nodes: ['source-file'], edges: ['is-implemented-in']};
    default:
      return {nodes: features, edges: ['is-implemented-in', 'participates-in', 'documents', 'tests', 'covers', 'automates', 'owner']};
  }
}

/** The build output under [kgDir], loaded for one report. */
export function readGraph(kgDir: string, repoRoot: string, system: TypeSystem, name: ReportName): GraphData {
  const want = needs(name, system);
  const data = empty(repoRoot);
  const manifest = readJson<Manifest>(path.join(kgDir, 'manifest.json'));
  data.sources = manifest?.sources ?? {};
  data.ownership = readJson<Ownership>(path.join(kgDir, 'reports', 'ownership.json'));
  data.problems = readJson<Record<string, string[]>>(path.join(kgDir, 'reports', 'problems.json')) ?? {};
  for (const type of want.nodes)
    for (const row of readJsonl(path.join(kgDir, 'data', 'nodes', `${type}.jsonl`))) {
      data.nodes.set(String(row.id), row);
      push(data.byType, type, row);
    }
  for (const edge of want.edges) data.edges.set(edge, readJsonl(path.join(kgDir, 'data', 'edges', `${edge}.jsonl`)));
  return data;
}

/** The same slice straight from the build that produced it, so `build` writes the reports without re-reading 100 MB. */
export function fromGraph(graph: Graph, repoRoot: string, sources: Record<string, string>): GraphData {
  const data = empty(repoRoot);
  data.sources = sources;
  data.ownership = graph.reports.ownership as Ownership | undefined;
  data.problems = graph.details;
  for (const row of graph.nodes) {
    data.nodes.set(String(row.id), row);
    push(data.byType, String(row.type), row);
  }
  for (const row of graph.edges) push(data.edges, row.type === 'ref' ? String(row.name) : String(row.type), row);
  return data;
}

export function makeReport(name: ReportName, data: GraphData, options: ReportOptions): Report {
  const report = name === 'orphans' ? orphans(data) : name === 'stale' ? stale(data, options)
    : name === 'coverage' ? coverage(data, options) : name === 'proposed' ? proposed(data) : diff(data, options);
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

export function printReport(report: Report, format: ReportFormat): void {
  if (format === 'json') {
    printOutput(report, 'json');
    return;
  }
  if (format === 'md') {
    process.stdout.write(markdown(report));
    return;
  }
  console.log(report.summary);
  for (const note of report.notes) console.log(note);
  for (const section of report.sections) {
    console.log(`\n${section.title} (${section.rows.length})`);
    printOutput(section.rows, 'table');
  }
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

/** Files no feature owns, grouped by the package or the core sub-project they sit in. */
function orphans(data: GraphData): Report {
  const groups = new Map<string, {file: string, loc: number}[]>();
  for (const orphan of data.ownership?.orphans ?? []) push(groups, groupOf(orphan.file), orphan);
  const rows = [...groups].map(([group, files]) => ({
    group, files: files.length, loc: files.reduce((sum, f) => sum + f.loc, 0),
    largest: [...files].sort((a, b) => b.loc - a.loc || compare(a.file, b.file)).slice(0, ORPHAN_FILES)
      .map((f) => `${f.file.slice(group.length + 1)} (${f.loc})`).join(', '),
  })).sort((a, b) => b.loc - a.loc || compare(a.group, b.group));
  const loc = rows.reduce((sum, r) => sum + r.loc, 0);
  return {
    name: 'orphans', title: 'Orphan files',
    summary: `${data.ownership?.orphans.length ?? 0} files with no owner in ${rows.length} groups, ${loc} lines; the ${Math.min(ORPHAN_GROUPS, rows.length)} largest groups below.`,
    notes: data.ownership ? [] : ['no reports/ownership.json: run grok kg build'],
    sections: [{title: 'groups', rows: rows.slice(0, ORPHAN_GROUPS)}],
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
  for (const issue of data.problems.home_issues ?? []) {
    const m = HOME_ISSUE.exec(issue);
    if (!m || !STALE_CODES.includes(m[3])) continue;
    rows.push({kind: 'citation', source: m[2] ? `${m[1]}:${m[2]}` : m[1], target: /'([^']+)'/.exec(m[4])?.[1] ?? '', reason: `${m[3]}: ${m[4]}`});
  }
  for (const edge of data.edges.get('tracked-in') ?? []) {
    const ticket = data.nodes.get(String(edge.to));
    if (ticket && ticket.provenance === 'external') continue;
    rows.push({kind: 'ticket', source: String(edge.from), target: String(edge.to), reason: 'not in the backlog snapshot'});
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
      rows.push({kind: 'declaration', source: String(edge.to), target: String(edge.from), reason: 'defined_by names a declaration the Dart batch does not have'});
    }
  const kinds = new Map<string, number>();
  for (const row of rows) kinds.set(String(row.kind), (kinds.get(String(row.kind)) ?? 0) + 1);
  return {
    name: 'stale', title: 'Stale references',
    summary: `${rows.length} stale references${kinds.size ? `: ${[...kinds].sort(([a], [b]) => compare(a, b)).map(([k, n]) => `${n} ${k}`).join(', ')}` : ''}.`,
    notes: [
      ...(data.sources.backlog === 'missing' ? ['the backlog snapshot was not read, so tracked-in tickets are not checked'] : []),
      ...(data.sources.dart === 'ok' ? [] : ['the Dart batch is not ok, so defined_by declarations are not checked']),
    ],
    sections: [{title: 'references', rows: rows.sort((a, b) => compare(String(a.kind), String(b.kind)) || compare(String(a.source), String(b.source)) || compare(String(a.target), String(b.target)))}],
  };
}

/** One row per feature: who owns it, what tests it, what documents it, and whether it says what it is. */
function coverage(data: GraphData, options: ReportOptions): Report {
  const features = [...data.nodes.values()].filter((r) => isSubtype(options.system, String(r.type), 'feature'));
  const tests = count(data.edges.get('tests'), 'to');
  const documents = count(data.edges.get('documents'), 'to');
  const owners = new Map((data.edges.get('owner') ?? []).map((e) => [String(e.from), String(e.to)]));
  const scenarios = new Map<string, string[]>();
  for (const edge of data.edges.get('covers') ?? []) push(scenarios, String(edge.to), String(edge.from));
  const automated = count(data.edges.get('automates'), 'to');
  const rows = features.map((f) => {
    const id = String(f.id);
    const covered = scenarios.get(id) ?? [];
    const scenarioTests = covered.reduce((sum, s) => sum + (automated.get(s) ?? 0), 0);
    return {
      feature: id, owner: (f.owner as string | undefined) ?? owners.get(id) ?? '', status: (f.status as string | undefined) ?? '',
      tests: tests.get(id) ?? 0, scenarios: covered.length, automated: scenarioTests,
      no_tests: !(tests.get(id) ?? 0) && !covered.length, no_docs: !documents.get(id), no_description: !f.description,
    };
  }).sort((a, b) => compare(a.feature, b.feature));
  const gaps = (key: 'no_tests' | 'no_docs' | 'no_description') => rows.filter((r) => r[key]).length;
  return {
    name: 'coverage', title: 'Feature coverage',
    summary: `${rows.length} features: ${gaps('no_tests')} with no test or scenario, ${gaps('no_docs')} with no document beside the home, ${gaps('no_description')} with no first paragraph.`,
    notes: [], sections: [{title: 'features', rows}],
  };
}

/** Folders that carry code no feature claims: candidates for a home of their own (conventions.md §10). */
function proposed(data: GraphData): Report {
  const owned = new Set((data.edges.get('is-implemented-in') ?? []).map((e) => String(e.to)));
  const folders = new Map<string, {files: number, loc: number, owned: boolean}>();
  for (const row of data.byType.get('source-file') ?? []) {
    const folder = folderOf(String(row.path ?? ''));
    if (!folder) continue;
    const entry = folders.get(folder) ?? {files: 0, loc: 0, owned: false};
    entry.files++;
    entry.loc += Number(row.loc ?? 0);
    entry.owned ||= owned.has(String(row.id));
    folders.set(folder, entry);
  }
  const rows = [...folders].filter(([, e]) => !e.owned)
    .map(([folder, e]) => ({path: folder, files: e.files, loc: e.loc, suggested: suggestedId(folder)}))
    .sort((a, b) => b.loc - a.loc || compare(a.path, b.path));
  return {
    name: 'proposed', title: 'Proposed features',
    summary: `${rows.length} folders with code and no owner, of ${folders.size} under ${PROPOSED_PARENTS.join(', ')}; ${rows.filter((r) => r.suggested).length} with a suggested id.`,
    notes: ['a folder is only seen through the source files the build extracted'],
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
  for (const file of files) {
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
      unowned++;
      continue;
    }
    for (const [feature, hit] of [...hits].sort(([a], [b]) => compare(a, b))) {
      fileRows.push({file, feature, ...hit});
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
    .map((e) => ({test: String(e.from), feature: String(e.to), kind: String(e.kind ?? '')})).sort((a, b) => compare(a.test, b.test));
  const scenarioRows = (data.edges.get('covers') ?? []).filter((e) => ids.has(String(e.to)))
    .map((e) => ({scenario: String(e.from), feature: String(e.to)})).sort((a, b) => compare(a.scenario, b.scenario));
  const covered = new Set(scenarioRows.map((r) => r.scenario));
  const automationRows = (data.edges.get('automates') ?? []).filter((e) => covered.has(String(e.to)))
    .map((e) => ({test: String(e.from), scenario: String(e.to)})).sort((a, b) => compare(a.test, b.test));
  return {
    name: 'diff', title: `Features touched by ${options.base}...HEAD`,
    summary: `${files.length} files changed, ${unowned} of them in no feature; ${featureRows.length} features touched; ` +
      `${testRows.length} tests, ${scenarioRows.length} scenarios and ${automationRows.length} automations cover them.`,
    notes: problems,
    sections: [
      {title: 'features', rows: featureRows},
      {title: 'files', rows: fileRows},
      {title: 'tests', rows: testRows},
      {title: 'scenarios', rows: scenarioRows},
      {title: 'automations', rows: automationRows},
    ],
  };
}

/** `git diff --name-only <base>...HEAD` in the monorepo and in the public submodule, whose paths the graph prefixes. */
function changedFiles(repoRoot: string, base: string): {files: string[], problems: string[]} {
  const files = new Set<string>();
  const problems: string[] = [];
  const roots = new Set<string>();
  for (const [cwd, prefix] of [[repoRoot, ''], [path.join(repoRoot, 'public'), 'public/']]) {
    // a checkout where public/ is a plain folder rather than the submodule answers the same diff twice
    const top = spawnSync('git', ['-C', cwd, 'rev-parse', '--show-toplevel'], {encoding: 'utf8'});
    if (top.status !== 0 || roots.has(top.stdout.trim())) continue;
    roots.add(top.stdout.trim());
    const r = spawnSync('git', ['-C', cwd, 'diff', '--name-only', `${base}...HEAD`], {encoding: 'utf8', maxBuffer: 32 * 1024 * 1024});
    if (r.status !== 0) {
      problems.push(`${prefix || 'the monorepo'}: git diff ${base}...HEAD failed: ${(r.stderr ?? '').trim().split('\n')[0] || `exit ${r.status}`}`);
      continue;
    }
    for (const line of r.stdout.split('\n')) {
      const file = line.trim();
      // the submodule's own gitlink is not a file of either repository
      if (file && file !== 'public') files.add(prefix + file);
    }
  }
  return {files: [...files].sort(compare), problems};
}

function empty(repoRoot: string): GraphData {
  return {nodes: new Map(), byType: new Map(), edges: new Map(), problems: {}, sources: {}, repoRoot};
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

function readJson<T>(file: string): T | undefined {
  return fs.existsSync(file) ? JSON.parse(fs.readFileSync(file, 'utf8')) as T : undefined;
}

function readJsonl(file: string): Row[] {
  return fs.existsSync(file) ? fs.readFileSync(file, 'utf8').split('\n').filter(Boolean).map((l) => JSON.parse(l) as Row) : [];
}

/** Code-point order, the same on every platform and locale (write.ts). */
function compare(a: string, b: string): number {
  return a < b ? -1 : a > b ? 1 : 0;
}
