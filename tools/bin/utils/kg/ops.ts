/// The bounded operations agents get before raw Cypher (conventions.md §11.1, build-plan.md WO-7):
/// impact, tests-for, explain and find, each a read-only template over the index. They answer with
/// sections of rows, and say up front when the manifest reports a source they could not read.
import {OutputFormat, printOutput} from '../server-output';
import {KuzuConnection, run, quote} from './kuzu';

export interface Section {
  title: string;
  rows: Record<string, unknown>[];
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
/** `find` is the vocabulary search (conventions.md §11.1): the authored types first, everything else after. */
const VOCABULARY = [['feature', 'concept'], ['scenario', 'initiative', 'package', 'library']];

/** A `file:` or `doc:` node carries the path it was made from; a home document is a `doc:` node with no file of its own. */
const NODE_PATH = /^(?:file|doc):(.+)$/;

/** Ids may be written with or without the sigil; a path names the source file it belongs to. */
function bareId(arg: string): string {
  return arg.trim().replace(/^~/, '');
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
  sections.push({title: 'features', rows: features});
  const ids = features.map((f) => String(f.feature));
  if (ids.length) {
    const owners = await run(conn, `MATCH (f:Feature)-[:${quote('owner')}]->(p) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, p.${quote('id')} AS owner, p.${quote('name')} AS name`, {ids});
    sections.push({title: 'owners', rows: owners.rows});
    sections.push({title: 'evidence', rows: await evidence(conn, ids, options.limit)});
    const docs = await run(conn, `MATCH (d)-[:${quote('DOCUMENTS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, d.${quote('id')} AS document, d.${quote('type')} AS type LIMIT ${options.limit}`, {ids});
    sections.push({title: 'documents', rows: docs.rows});
    const work = await run(conn, `MATCH (t)-[:${quote('AFFECTS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, t.${quote('id')} AS ticket, 'affects' AS relation, t.${quote('state')} AS state LIMIT ${options.limit}`, {ids});
    const tracked = await run(conn, `MATCH (f:Feature)-[:${quote('TRACKED_IN')}]->(t) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, t.${quote('id')} AS ticket, 'tracked-in' AS relation, t.${quote('state')} AS state LIMIT ${options.limit}`, {ids});
    sections.push({title: 'work', rows: [...work.rows, ...tracked.rows]});
  }
  const callers = await callersOf(conn, id, options.limit);
  if (callers) sections.push({title: 'callers', rows: callers});
  return {op: 'impact', target, sections};
}

/** Who reaches this: a declaration through uses and calls, a source file through its imports and its declarations. */
async function callersOf(conn: KuzuConnection, id: string, limit: number): Promise<Record<string, unknown>[] | null> {
  if (/^(decl|func|ep):/.test(id)) {
    const {rows} = await run(conn, `MATCH (c)-[e:${quote('USES')}|${quote('CALLS')}]->(n) WHERE n.${quote('id')} = $id ` +
      `RETURN c.${quote('id')} AS caller, label(e) AS via, n.${quote('id')} AS target, e.${quote('count')} AS count ORDER BY count DESC LIMIT ${limit}`, {id});
    return rows;
  }
  if (!id.startsWith('file:')) return null;
  const importers = await run(conn, `MATCH (c)-[:${quote('IMPORTS')}]->(n) WHERE n.${quote('id')} = $id ` +
    `RETURN c.${quote('id')} AS caller, 'IMPORTS' AS via, n.${quote('id')} AS target, null AS count LIMIT ${limit}`, {id});
  const users = await run(conn, `MATCH (c)-[e:${quote('USES')}]->(d)<-[:${quote('DECLARES')}]-(n) WHERE n.${quote('id')} = $id ` +
    `RETURN c.${quote('id')} AS caller, 'USES' AS via, d.${quote('id')} AS target, e.${quote('count')} AS count ORDER BY count DESC LIMIT ${limit}`, {id});
  return [...importers.rows, ...users.rows];
}

/** The tests, scenarios and automations of a feature and everything under it. */
export async function testsFor(conn: KuzuConnection, target: Record<string, unknown>, options: OpsOptions): Promise<OpsResult> {
  const roots = (await featuresOf(conn, target)).map((r) => String(r.feature));
  const sections: Section[] = [];
  if (!roots.length) return {op: 'tests-for', target, sections: [{title: 'features', rows: []}]};
  const descendants = await run(conn, `MATCH (d:Feature)-[:${quote('PART_OF')}*0..5]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
    `RETURN DISTINCT d.${quote('id')} AS feature, d.${quote('name')} AS name, d.${quote('status')} AS status`, {ids: roots});
  const ids = descendants.rows.map((r) => String(r.feature));
  sections.push({title: 'features', rows: descendants.rows});
  const tests = await run(conn, `MATCH (t)-[:${quote('TESTS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
    `RETURN t.${quote('framework')} AS framework, t.${quote('level')} AS level, t.${quote('id')} AS test, f.${quote('id')} AS feature, t.${quote('skipped')} AS skipped ` +
    `ORDER BY framework, test LIMIT ${options.limit}`, {ids});
  sections.push({title: 'tests', rows: tests.rows});
  const scenarios = await run(conn, `MATCH (s)-[:${quote('COVERS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
    `RETURN s.${quote('id')} AS scenario, f.${quote('id')} AS feature, s.${quote('manual_only')} AS manual_only, s.${quote('priority')} AS priority ` +
    `ORDER BY scenario LIMIT ${options.limit}`, {ids});
  sections.push({title: 'scenarios', rows: scenarios.rows});
  const covered = scenarios.rows.map((r) => String(r.scenario));
  const automations = covered.length ? await run(conn, `MATCH (t)-[:${quote('AUTOMATES')}]->(s) WHERE s.${quote('id')} IN $ids ` +
    `RETURN t.${quote('framework')} AS framework, t.${quote('id')} AS test, s.${quote('id')} AS scenario ORDER BY framework, test LIMIT ${options.limit}`, {ids: covered}) : {rows: []};
  sections.push({title: 'automations', rows: automations.rows});
  return {op: 'tests-for', target, sections};
}

/** Everything one hop away from a node, grouped by edge type and direction. */
export async function explain(conn: KuzuConnection, target: Record<string, unknown>, options: OpsOptions): Promise<OpsResult> {
  const id = String(target.id);
  const node = await run(conn, `MATCH (n) WHERE n.${quote('id')} = $id RETURN n`, {id});
  const properties = Object.entries((node.rows[0]?.n ?? {}) as Record<string, unknown>)
    .filter(([k, v]) => k !== '_label' && v !== null && v !== undefined && !(Array.isArray(v) && !v.length))
    .map(([property, value]) => ({property, value}));
  const sections: Section[] = [{title: 'properties', rows: properties}];
  const out = await run(conn, `MATCH (n)-[e]->(m) WHERE n.${quote('id')} = $id RETURN label(e) AS edge, m.${quote('id')} AS other, m.${quote('type')} AS type`, {id});
  const into = await run(conn, `MATCH (n)<-[e]-(m) WHERE n.${quote('id')} = $id RETURN label(e) AS edge, m.${quote('id')} AS other, m.${quote('type')} AS type`, {id});
  sections.push({title: 'edges', rows: [...group(out.rows, 'out'), ...group(into.rows, 'in')].slice(0, options.limit)});
  return {op: 'explain', target, sections};
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
  return {op: 'find', target: {text}, sections: [{title: 'matches', rows: sorted.slice(0, options.limit)}]};
}

async function evidence(conn: KuzuConnection, ids: string[], limit: number): Promise<Record<string, unknown>[]> {
  const rows: Record<string, unknown>[] = [];
  for (const [edge, kind] of [['TESTS', 'test'], ['COVERS', 'scenario']]) {
    const found = await run(conn, `MATCH (a)-[:${quote(edge)}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, a.${quote('id')} AS artifact, '${kind}' AS kind, a.${quote('framework')} AS framework LIMIT ${limit}`, {ids});
    rows.push(...found.rows);
  }
  const automated = await run(conn, `MATCH (t)-[:${quote('AUTOMATES')}]->(s)-[:${quote('COVERS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
    `RETURN f.${quote('id')} AS feature, t.${quote('id')} AS artifact, 'automation' AS kind, t.${quote('framework')} AS framework LIMIT ${limit}`, {ids});
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
    console.log(`\n${section.title} (${section.rows.length})`);
    printOutput(section.rows, output);
  }
}
