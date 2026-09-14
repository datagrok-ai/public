/// The bounded operations agents get before raw Cypher (conventions.md §11.1, build-plan.md WO-7):
/// impact, tests-for, explain and find, each a read-only template over the index. They answer with
/// sections of rows, and say up front when the manifest reports a source they could not read.
import {OutputFormat, printOutput} from '../server-output';
import {KuzuConnection, run, quote} from './kuzu';

export interface Section {
  title: string;
  rows: Record<string, unknown>[];
  /** Rows before paging, when the section was paged: it is counted whole first, so a header can say `50 of 394`. */
  total?: number;
  /** The whole header line of a section that came back empty, saying why it did. */
  empty?: string;
}

export interface OpsResult {
  op: string;
  target: Record<string, unknown> | null;
  /** The Dart-coverage clause, when the manifest says the batch is not `ok`. */
  note?: string;
  sections: Section[];
}

export interface OpsOptions {
  limit: number;
}

export const DEFAULT_LIMIT = 50;
const HOP_TARGETS = 20;
const SEARCH_SCAN = 500;
/** An op prints wider cells than `grok s`: its lists are ids, and half an id is worse than none. */
const CELL_BUDGET = 80;
const NO_AFFECTS = 'features: none derivable (no affects edges from tickets; the Jira Feature field is not populated yet)';
const UNOWNED = 'features (0): no home document owns this file; see grok kg report proposed';
/** `find` is the vocabulary search (conventions.md §11.1): the authored types first, everything else after. */
const VOCABULARY = [['feature', 'concept'], ['scenario', 'initiative', 'package', 'library']];

/** A `file:` or `doc:` node carries the path it was made from; a home document is a `doc:` node with no file of its own. */
const NODE_PATH = /^(?:file|doc):(.+)$/;

/** Ids may be written with or without the sigil; a path names the source file it belongs to. */
function bareId(arg: string): string {
  return arg.trim().replace(/^~/, '');
}

/** Every section is counted whole and paged after, so the header can tell a page from the total. */
function section(title: string, rows: Record<string, unknown>[], limit: number, empty?: string): Section {
  const paged: Section = {title, rows: rows.slice(0, limit), total: rows.length};
  if (!rows.length && empty) paged.empty = empty;
  return paged;
}

export async function resolveTarget(conn: KuzuConnection, arg: string): Promise<Record<string, unknown> | null> {
  const id = bareId(arg);
  const p = id.replace(/\\/g, '/');
  for (const candidate of [id, `file:${p}`, `doc:${p}`]) {
    const {rows} = await run(conn, `MATCH (n) WHERE n.${quote('id')} = $id ` +
      `RETURN n.${quote('id')} AS id, label(n) AS root, n.${quote('type')} AS type, n.${quote('name')} AS name, n.${quote('status')} AS status`, {id: candidate});
    if (rows.length) return rows[0];
  }
  return null;
}

/** The index describes itself: `find` needs the tables, and the columns each of them really has. */
async function nodeTables(conn: KuzuConnection): Promise<string[]> {
  const {rows} = await run(conn, 'CALL show_tables() RETURN name, type');
  return rows.filter((r) => r.type === 'NODE').map((r) => String(r.name));
}

async function columnsOf(conn: KuzuConnection, table: string): Promise<string[]> {
  const {rows} = await run(conn, `CALL table_info('${table}') RETURN name`);
  return rows.map((r) => String(r.name));
}

/**
 * The features a node belongs to, strongest relation first, the way report.ts's diff reads a changed file: the features
 * that own it, the one whose home it is, the ones it documents, and the ones it takes part in.
 */
async function featuresOf(conn: KuzuConnection, target: Record<string, unknown>): Promise<Record<string, unknown>[]> {
  const id = String(target.id);
  if (target.root === 'Feature') return [{feature: id, relation: 'self', name: target.name, status: target.status ?? null}];
  const select = (relation: string) => `RETURN f.${quote('id')} AS feature, '${relation}' AS relation, f.${quote('name')} AS name, f.${quote('status')} AS status`;
  const owns = await run(conn, `MATCH (f:Feature)-[:${quote('IS_IMPLEMENTED_IN')}]->(n) WHERE n.${quote('id')} = $id ${select('owns')}`, {id});
  const path = NODE_PATH.exec(id)?.[1];
  const home = path ? await run(conn, `MATCH (f:Feature) WHERE f.${quote('home')} = $path ${select('home')}`, {path}) : {rows: []};
  const documents = await run(conn, `MATCH (n)-[:${quote('DOCUMENTS')}]->(f:Feature) WHERE n.${quote('id')} = $id ${select('documents')}`, {id});
  const part = await run(conn, `MATCH (n)-[:${quote('PARTICIPATES_IN')}]->(f:Feature) WHERE n.${quote('id')} = $id ${select('participates')}`, {id});
  const features: Record<string, unknown>[] = [];
  for (const row of [...owns.rows, ...home.rows, ...documents.rows, ...part.rows])
    if (!features.some((f) => f.feature === row.feature)) features.push(row);
  return features;
}

/** What a change to this file, declaration or feature reaches: the features that own it, their evidence and their work. */
export async function impact(conn: KuzuConnection, target: Record<string, unknown>, options: OpsOptions): Promise<OpsResult> {
  const id = String(target.id);
  const sections: Section[] = [];
  const features = await featuresOf(conn, target);
  sections.push(section('features', features, options.limit, NODE_PATH.test(id) ? UNOWNED : undefined));
  const ids = features.map((f) => String(f.feature));
  if (ids.length) {
    const owners = await run(conn, `MATCH (f:Feature)-[:${quote('owner')}]->(p) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, p.${quote('id')} AS owner, p.${quote('name')} AS name`, {ids});
    sections.push(section('owners', owners.rows, options.limit));
    sections.push(section('evidence', await evidence(conn, ids), options.limit));
    const docs = await run(conn, `MATCH (d)-[:${quote('DOCUMENTS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, d.${quote('id')} AS document, d.${quote('type')} AS type`, {ids});
    sections.push(section('documents', docs.rows, options.limit));
    const work = await run(conn, `MATCH (t)-[:${quote('AFFECTS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, t.${quote('id')} AS ticket, 'affects' AS relation, t.${quote('state')} AS state`, {ids});
    const tracked = await run(conn, `MATCH (f:Feature)-[:${quote('TRACKED_IN')}]->(t) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, t.${quote('id')} AS ticket, 'tracked-in' AS relation, t.${quote('state')} AS state`, {ids});
    sections.push(section('work', [...work.rows, ...tracked.rows], options.limit));
  }
  const callers = await callersOf(conn, id);
  if (callers) sections.push(section('callers', callers, options.limit));
  return {op: 'impact', target, sections};
}

/** Who reaches this: a declaration through uses and calls, a source file through its imports and its declarations. */
async function callersOf(conn: KuzuConnection, id: string): Promise<Record<string, unknown>[] | null> {
  if (/^(decl|func|ep):/.test(id)) {
    const {rows} = await run(conn, `MATCH (c)-[e:${quote('USES')}|${quote('CALLS')}]->(n) WHERE n.${quote('id')} = $id ` +
      `RETURN c.${quote('id')} AS caller, label(e) AS via, n.${quote('id')} AS target, e.${quote('count')} AS count ORDER BY count DESC`, {id});
    return rows;
  }
  if (!id.startsWith('file:')) return null;
  const importers = await run(conn, `MATCH (c)-[:${quote('IMPORTS')}]->(n) WHERE n.${quote('id')} = $id ` +
    `RETURN c.${quote('id')} AS caller, 'IMPORTS' AS via, n.${quote('id')} AS target, null AS count`, {id});
  const users = await run(conn, `MATCH (c)-[e:${quote('USES')}]->(d)<-[:${quote('DECLARES')}]-(n) WHERE n.${quote('id')} = $id ` +
    `RETURN c.${quote('id')} AS caller, 'USES' AS via, d.${quote('id')} AS target, e.${quote('count')} AS count ORDER BY count DESC`, {id});
  return [...importers.rows, ...users.rows];
}

/** The tests, scenarios and automations of a feature and everything under it. */
export async function testsFor(conn: KuzuConnection, target: Record<string, unknown>, options: OpsOptions): Promise<OpsResult> {
  const roots = (await featuresOf(conn, target)).map((r) => String(r.feature));
  const sections: Section[] = [];
  if (!roots.length)
    return {op: 'tests-for', target, sections: [section('features', [], options.limit, NODE_PATH.test(String(target.id)) ? UNOWNED : undefined)]};
  const descendants = await run(conn, `MATCH (d:Feature)-[:${quote('PART_OF')}*0..5]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
    `RETURN DISTINCT d.${quote('id')} AS feature, d.${quote('name')} AS name, d.${quote('status')} AS status`, {ids: roots});
  const ids = descendants.rows.map((r) => String(r.feature));
  sections.push(section('features', descendants.rows, options.limit));
  const tests = await run(conn, `MATCH (t)-[:${quote('TESTS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
    `RETURN t.${quote('framework')} AS framework, t.${quote('level')} AS level, t.${quote('id')} AS test, f.${quote('id')} AS feature, t.${quote('skipped')} AS skipped ` +
    `ORDER BY framework, test`, {ids});
  sections.push(section('tests', tests.rows, options.limit,
    `tests (0): no test carries ~${roots.join(', ~')} and no owned file contains tests`));
  const scenarios = await run(conn, `MATCH (s)-[:${quote('COVERS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
    `RETURN s.${quote('id')} AS scenario, f.${quote('id')} AS feature, s.${quote('manual_only')} AS manual_only, s.${quote('priority')} AS priority ` +
    `ORDER BY scenario`, {ids});
  sections.push(section('scenarios', scenarios.rows, options.limit));
  const covered = scenarios.rows.map((r) => String(r.scenario));
  const automations = covered.length ? await run(conn, `MATCH (t)-[:${quote('AUTOMATES')}]->(s) WHERE s.${quote('id')} IN $ids ` +
    `RETURN t.${quote('framework')} AS framework, t.${quote('id')} AS test, s.${quote('id')} AS scenario ORDER BY framework, test`, {ids: covered}) : {rows: []};
  sections.push(section('automations', automations.rows, options.limit));
  return {op: 'tests-for', target, sections};
}

/** Everything one hop away from a node, grouped by edge type and direction. */
export async function explain(conn: KuzuConnection, target: Record<string, unknown>, options: OpsOptions): Promise<OpsResult> {
  const id = String(target.id);
  const node = await run(conn, `MATCH (n) WHERE n.${quote('id')} = $id RETURN n`, {id});
  const properties = Object.entries((node.rows[0]?.n ?? {}) as Record<string, unknown>)
    .filter(([k, v]) => k !== '_label' && v !== null && v !== undefined && !(Array.isArray(v) && !v.length))
    .map(([property, value]) => ({property, value}));
  const sections: Section[] = [section('properties', properties, options.limit)];
  if (target.type === 'release') sections.push(await shipped(conn, id, options.limit));
  const out = await run(conn, `MATCH (n)-[e]->(m) WHERE n.${quote('id')} = $id RETURN label(e) AS edge, m.${quote('id')} AS other, m.${quote('type')} AS type`, {id});
  const into = await run(conn, `MATCH (n)<-[e]-(m) WHERE n.${quote('id')} = $id RETURN label(e) AS edge, m.${quote('id')} AS other, m.${quote('type')} AS type`, {id});
  sections.push(section('edges', [...group(out.rows, 'out'), ...group(into.rows, 'in')], options.limit));
  return {op: 'explain', target, sections};
}

/** What a release shipped: the features its tickets affect, whether the ticket targets it or a commit in it resolves the ticket. */
async function shipped(conn: KuzuConnection, id: string, limit: number): Promise<Section> {
  const select = `RETURN DISTINCT f.${quote('id')} AS feature, f.${quote('name')} AS name, f.${quote('status')} AS status`;
  const targeted = await run(conn, `MATCH (t)-[:${quote('TARGETS_RELEASE')}]->(r) WHERE r.${quote('id')} = $id ` +
    `MATCH (t)-[:${quote('AFFECTS')}]->(f:Feature) ${select}`, {id});
  const included = await run(conn, `MATCH (r)-[:${quote('INCLUDES')}]->(c)-[:${quote('RESOLVES')}]->(t)-[:${quote('AFFECTS')}]->(f:Feature) ` +
    `WHERE r.${quote('id')} = $id ${select}`, {id});
  const rows = [...targeted.rows];
  for (const row of included.rows)
    if (!rows.some((r) => r.feature === row.feature)) rows.push(row);
  return section('features', rows, limit, NO_AFFECTS);
}

/** One row per edge type and direction: how many, and the first few targets. */
function group(rows: Record<string, unknown>[], direction: 'in' | 'out'): Record<string, unknown>[] {
  const byEdge = new Map<string, Record<string, unknown>[]>();
  for (const row of rows) {
    const key = String(row.edge);
    let list = byEdge.get(key);
    if (!list) byEdge.set(key, list = []);
    list.push(row);
  }
  return [...byEdge].sort(([a], [b]) => a < b ? -1 : 1).map(([edge, list]) => ({
    edge, direction, count: list.length,
    targets: list.slice(0, HOP_TARGETS).map((r) => String(r.other)).join(', ') + (list.length > HOP_TARGETS ? `, +${list.length - HOP_TARGETS} more` : ''),
  }));
}

/** The vocabulary search: ids, names, aliases, descriptions and keywords, table by table, exact first. */
export async function find(conn: KuzuConnection, text: string, options: OpsOptions): Promise<OpsResult> {
  const q = text.trim().toLowerCase();
  const rows: Record<string, unknown>[] = [];
  for (const table of await nodeTables(conn)) {
    const columns = await columnsOf(conn, table);
    const where = [
      ...['id', 'name', 'description'].filter((c) => columns.includes(c)).map((c) => `lower(n.${quote(c)}) CONTAINS q`),
      // kuzu 0.11.3 crashes on `any(x IN list ...)` over a column that is null on some rows, and on a parameter
      // inside a lambda, and on a column the table does not have: list_filter, a WITH binding and the real columns.
      ...['aliases', 'keywords'].filter((c) => columns.includes(c)).map((c) => `size(list_filter(n.${quote(c)}, x -> lower(x) CONTAINS q)) > 0`),
    ].join(' OR ');
    const found = await run(conn, `WITH $q AS q MATCH (n:${quote(table)}) WHERE ${where} ` +
      `RETURN n.${quote('id')} AS id, n.${quote('type')} AS type, n.${quote('name')} AS name, n.${quote('status')} AS status LIMIT ${SEARCH_SCAN}`, {q});
    rows.push(...found.rows);
  }
  const rank = (row: Record<string, unknown>) => {
    const id = String(row.id).toLowerCase(), name = String(row.name ?? '').toLowerCase();
    if (id === q || name === q) return 0;
    if (id.startsWith(q) || name.startsWith(q) || id.split('/').pop() === q) return 1;
    return 2;
  };
  const vocabulary = (row: Record<string, unknown>) => {
    const tier = VOCABULARY.findIndex((types) => types.includes(String(row.type)));
    return tier < 0 ? VOCABULARY.length : tier;
  };
  const sorted = rows.map((row) => ({...row, match: rank(row)}))
    .sort((a, b) => vocabulary(a) - vocabulary(b) || a.match - b.match || (String(a.id) < String(b.id) ? -1 : 1));
  return {op: 'find', target: {text}, sections: [section('matches', sorted, options.limit)]};
}

async function evidence(conn: KuzuConnection, ids: string[]): Promise<Record<string, unknown>[]> {
  const rows: Record<string, unknown>[] = [];
  for (const [edge, kind] of [['TESTS', 'test'], ['COVERS', 'scenario']]) {
    const found = await run(conn, `MATCH (a)-[:${quote(edge)}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, a.${quote('id')} AS artifact, '${kind}' AS kind, a.${quote('framework')} AS framework`, {ids});
    rows.push(...found.rows);
  }
  const automated = await run(conn, `MATCH (t)-[:${quote('AUTOMATES')}]->(s)-[:${quote('COVERS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
    `RETURN f.${quote('id')} AS feature, t.${quote('id')} AS artifact, 'automation' AS kind, t.${quote('framework')} AS framework`, {ids});
  return [...rows, ...automated.rows];
}

/** What the manifest says about the Dart batch; every op leads with it unless the batch is `ok`. */
export function coverageNote(sources: Record<string, string> | undefined): string | undefined {
  const status = sources?.dart;
  if (status === 'ok') return undefined;
  return status === undefined || status === 'missing' ? 'Dart coverage unknown (no kg-dart batch)' : `Dart coverage ${status} (the kg-dart batch is ${status})`;
}

export function printOps(result: OpsResult, output: OutputFormat): void {
  if (output === 'json') {
    printOutput(result, 'json');
    return;
  }
  if (result.note) console.log(result.note);
  if (result.target && result.op !== 'find') console.log(`${result.op} ${result.target.id}${result.target.name ? ` (${result.target.name})` : ''}`);
  for (const section of result.sections) {
    if (!section.rows.length) {
      console.log(`\n${section.empty ?? `${section.title} (0)`}`);
      continue;
    }
    const total = section.total ?? section.rows.length;
    const count = section.rows.length < total ? `${section.rows.length} of ${total}; --limit to see more` : `${total}`;
    console.log(`\n${section.title} (${count})`);
    printOutput(section.rows, output, CELL_BUDGET);
  }
}
