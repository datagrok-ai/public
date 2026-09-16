/// The bounded operations agents get before raw Cypher (conventions.md §11.1, build-plan.md WO-7):
/// impact, tests-for, explain and find, each a read-only template over the index, answering with the
/// sections of `answer.ts`.
import {KuzuConnection, run, quote} from './kuzu';
import {Section, Answer, EdgeGroup, section} from './answer';

export interface OpsOptions {
  limit: number;
  /** The manifest's `edge_groups`, for `explain` to head its edges by folder. */
  groups?: Record<string, string[]>;
  /** `find`: how many substring matches one table may contribute. The exact query is not bounded by it. */
  scan?: number;
}

export const DEFAULT_LIMIT = 50;
const HOP_TARGETS = 20;
const SEARCH_SCAN = 500;
const NO_AFFECTS = 'shipped (0): none derivable (no affects edges from tickets; the Jira Feature field is not populated yet)';
const UNOWNED = 'features (0): no home document owns this file; see grok kg report proposed';
/** How many nodes the containment expansion may reach before it answers from what it has. */
const CONTAINED_CAP = 5000;
/** `find` is the vocabulary search (conventions.md §11.1): the authored types first, everything else after. */
const VOCABULARY = [['feature', 'concept'], ['scenario', 'initiative', 'package', 'library']];
/** How a node reaches a feature, strongest first; `report diff` ranks the same way. */
const RELATIONS = ['self', 'owns', 'home', 'documents', 'participates'];
/** A node that holds others: a package or a library declares its contents, a file its declarations. */
const CONTAINER = /^(?:pkg|lib|file):/;

/** A `file:` or `doc:` node carries the path it was made from; a home document is a `doc:` node with no file of its own. */
const NODE_PATH = /^(?:file|doc):(.+)$/;
/** Reference properties materialize as rel tables named after the property (§7.6); they are in no folder and share one heading. */
const REFERENCE_GROUP = 'reference';

/** The group headings in order: the manifest's `edge_groups` folders, then the reference predicates. */
export function groupOrder(groups: Record<string, string[]> = {}): string[] {
  return [...Object.keys(groups), REFERENCE_GROUP];
}

/** The group heading a rel table belongs under: its folder per the manifest's `edge_groups`, `reference` for the lowercase
 * tables a reference property makes (§7.5), nothing over a generation that predates `edge_groups`. */
export function groupOf(label: string, groups: Record<string, string[]> = {}): string {
  for (const [group, labels] of Object.entries(groups))
    if (labels.includes(label)) return group;
  return label === label.toLowerCase() ? REFERENCE_GROUP : '';
}

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

/** One step of a reasoning path, rendered the way it is read: `→ declares →` out of a node, `→ owner ←` into it. */
function hop(edge: string, direction: 'in' | 'out'): string {
  return `→ ${edge} ${direction === 'out' ? '→' : '←'}`;
}

/**
 * The target, what it contains and what contains it — one bounded step each way (review 3 #9). A package or a file
 * declares its functions, files and declarations; anything declared belongs to the file or package that declares it
 * and to the functions it implements. Features are then looked for over the whole set, so a change to a declaration
 * reaches the feature that claims it and a package answers for the code it holds.
 */
async function contained(conn: KuzuConnection, target: Record<string, unknown>): Promise<Map<string, string[]>> {
  const id = String(target.id);
  const reached = new Map<string, string[]>([[id, [id]]]);
  const add = (from: string, edge: string, direction: 'in' | 'out', to: string) => {
    if (reached.has(to) || reached.size >= CONTAINED_CAP) return;
    reached.set(to, [...reached.get(from)!, hop(edge, direction), to]);
  };
  const step = async (ids: string[], pattern: string, edge: string, direction: 'in' | 'out') => {
    if (!ids.length) return;
    const {rows} = await run(conn, `MATCH ${pattern} WHERE n.${quote('id')} IN $ids ` +
      `RETURN n.${quote('id')} AS anchor, m.${quote('id')} AS other ORDER BY anchor, other`, {ids});
    for (const row of rows) add(String(row.anchor), edge, direction, String(row.other));
  };
  if (CONTAINER.test(id)) await step([id], `(n)-[:${quote('DECLARES')}]->(m)`, 'declares', 'out');
  const declared = [...reached.keys()].filter((k) => !CONTAINER.test(k));
  await step(declared, `(m)-[:${quote('DECLARES')}]->(n)`, 'declares', 'in');
  await step(declared, `(n)-[:${quote('IMPLEMENTS')}]->(m)`, 'implements', 'out');
  return reached;
}

/**
 * The features the target belongs to through containment, strongest relation first: the features that own it or
 * anything it contains, the one whose home it is, the ones it documents, and the ones it takes part in. Each row
 * carries the chain that produced it, so the answer can be checked rather than trusted.
 */
async function featuresOf(conn: KuzuConnection, target: Record<string, unknown>): Promise<Record<string, unknown>[]> {
  const id = String(target.id);
  if (target.root === 'Feature')
    return [{feature: id, relation: 'self', name: target.name, status: target.status ?? null, via: id, path: [id]}];
  const reached = await contained(conn, target);
  const ids = [...reached.keys()];
  const select = (anchor: string) => `RETURN ${anchor} AS anchor, f.${quote('id')} AS feature, f.${quote('name')} AS name, f.${quote('status')} AS status`;
  const owns = await run(conn, `MATCH (f:Feature)-[:${quote('IS_IMPLEMENTED_IN')}]->(n) WHERE n.${quote('id')} IN $ids ${select(`n.${quote('id')}`)}`, {ids});
  const paths = ids.map((i) => NODE_PATH.exec(i)?.[1]).filter((p): p is string => !!p);
  const home = paths.length
    ? await run(conn, `MATCH (f:Feature) WHERE f.${quote('home')} IN $paths ${select(`f.${quote('home')}`)}`, {paths}) : {rows: []};
  const documents = await run(conn, `MATCH (n)-[:${quote('DOCUMENTS')}]->(f:Feature) WHERE n.${quote('id')} IN $ids ${select(`n.${quote('id')}`)}`, {ids});
  const part = await run(conn, `MATCH (n)-[:${quote('PARTICIPATES_IN')}]->(f:Feature) WHERE n.${quote('id')} IN $ids ${select(`n.${quote('id')}`)}`, {ids});
  const byFeature = new Map<string, Record<string, unknown>>();
  const consider = (rows: Record<string, unknown>[], relation: string, edge: string, direction: 'in' | 'out', byPath: boolean) => {
    for (const row of rows) {
      const anchor = byPath ? `doc:${row.anchor}` : String(row.anchor);
      const chain = reached.get(anchor) ?? reached.get(`file:${row.anchor}`);
      if (!chain) continue;
      const path = [...chain, hop(edge, direction), String(row.feature)];
      const candidate = {feature: row.feature, relation, name: row.name, status: row.status ?? null, via: path.join(' '), path};
      const best = byFeature.get(String(row.feature));
      if (!best || stronger(candidate, best)) byFeature.set(String(row.feature), candidate);
    }
  };
  consider(owns.rows, 'owns', 'is-implemented-in', 'in', false);
  consider(home.rows, 'home', 'home', 'in', true);
  consider(documents.rows, 'documents', 'documents', 'out', false);
  consider(part.rows, 'participates', 'participates-in', 'out', false);
  return [...byFeature.values()].sort((a, b) =>
    RELATIONS.indexOf(String(a.relation)) - RELATIONS.indexOf(String(b.relation)) ||
    (a.path as string[]).length - (b.path as string[]).length ||
    (String(a.feature) < String(b.feature) ? -1 : 1));
}

/** A shorter chain beats a longer one; at equal length the stronger relation wins. */
function stronger(a: Record<string, unknown>, b: Record<string, unknown>): boolean {
  const length = (r: Record<string, unknown>) => (r.path as string[]).length;
  if (length(a) !== length(b)) return length(a) < length(b);
  const rank = (r: Record<string, unknown>) => RELATIONS.indexOf(String(r.relation));
  return rank(a) !== rank(b) ? rank(a) < rank(b) : String(a.via) < String(b.via);
}

/** What a change to this file, declaration or feature reaches: the features that own it, their evidence and their work. */
export async function impact(conn: KuzuConnection, target: Record<string, unknown>, options: OpsOptions): Promise<Answer> {
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
  else sections.push(await packageOwners(conn, id, options.limit));
  const callers = await callersOf(conn, id);
  if (callers) sections.push(section('callers', callers, options.limit));
  return {op: 'impact', target, sections};
}

/**
 * The owner a file has when no feature owns it: the package or library that declares it answers for it through the
 * author of its `package.json` (conventions.md §5.1). One `declares` hop for a file, two for what a file declares.
 */
async function packageOwners(conn: KuzuConnection, id: string, limit: number): Promise<Section> {
  const {rows} = await run(conn, `MATCH (c)-[:${quote('DECLARES')}*1..2]->(n), (c)-[:${quote('owner')}]->(p) WHERE n.${quote('id')} = $id ` +
    `RETURN DISTINCT c.${quote('id')} AS package, p.${quote('id')} AS owner, p.${quote('name')} AS name`, {id});
  return section('owners', rows.map((r) => ({...r,
    via: [id, hop('declares', 'in'), String(r.package), hop('owner', 'out'), String(r.owner)].join(' ')})), limit);
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
export async function testsFor(conn: KuzuConnection, target: Record<string, unknown>, options: OpsOptions): Promise<Answer> {
  const found = await featuresOf(conn, target);
  const roots = found.map((r) => String(r.feature));
  const sections: Section[] = [];
  if (!roots.length)
    return {op: 'tests-for', target, sections: [section('features', [], options.limit, NODE_PATH.test(String(target.id)) ? UNOWNED : undefined),
      await packageOwners(conn, String(target.id), options.limit)]};
  const via = new Map(found.map((r) => [String(r.feature), r]));
  const descendants = await run(conn, `MATCH (d:Feature)-[:${quote('PART_OF')}*0..5]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
    `RETURN DISTINCT d.${quote('id')} AS feature, d.${quote('name')} AS name, d.${quote('status')} AS status, f.${quote('id')} AS root`, {ids: roots});
  const reached = new Map<string, Record<string, unknown>>();
  for (const r of descendants.rows) {
    if (reached.has(String(r.feature))) continue;
    const root = via.get(String(r.root))!;
    const path = r.feature === r.root ? root.path as string[] : [...root.path as string[], hop('part-of', 'in'), String(r.feature)];
    reached.set(String(r.feature), {feature: r.feature, name: r.name, status: r.status, via: path.join(' '), path});
  }
  const ids = [...reached.keys()];
  sections.push(section('features', [...reached.values()], options.limit));
  const tests = await run(conn, `MATCH (t)-[:${quote('TESTS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
    `RETURN t.${quote('framework')} AS framework, t.${quote('level')} AS level, t.${quote('id')} AS test, f.${quote('id')} AS feature, ` +
    `t.${quote('skipped')} AS skipped, t.${quote('dynamic')} AS dynamic ORDER BY framework, test`, {ids});
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

/** Everything one hop away from a node, grouped by edge type and direction, with the evidence behind each group. */
export async function explain(conn: KuzuConnection, target: Record<string, unknown>, options: OpsOptions): Promise<Answer> {
  const id = String(target.id);
  const node = await run(conn, `MATCH (n) WHERE n.${quote('id')} = $id RETURN n`, {id});
  const properties = Object.entries((node.rows[0]?.n ?? {}) as Record<string, unknown>)
    .filter(([k, v]) => k !== '_label' && v !== null && v !== undefined && !(Array.isArray(v) && !v.length))
    .map(([property, value]) => ({property, value}));
  const sections: Section[] = [section('properties', properties, options.limit)];
  if (target.type === 'release') sections.push(...await releaseSections(conn, node.rows[0]?.n as Record<string, unknown> ?? {}, id, options.limit));
  const select = `RETURN label(e) AS edge, m.${quote('id')} AS other, m.${quote('name')} AS name, m.${quote('type')} AS type, ` +
    `e.${quote('derived_by')} AS derived_by, e.${quote('confidence')} AS confidence, e.${quote('evidence')} AS evidence`;
  const out = await run(conn, `MATCH (n)-[e]->(m) WHERE n.${quote('id')} = $id ${select} ORDER BY edge, other`, {id});
  const into = await run(conn, `MATCH (n)<-[e]-(m) WHERE n.${quote('id')} = $id ${select} ORDER BY edge, other`, {id});
  const order = groupOrder(options.groups);
  const rank = (row: Record<string, unknown>) => {
    const at = order.indexOf(String(row.group));
    return at < 0 ? order.length : at;
  };
  const edges = [...group(out.rows, 'out', options.groups), ...group(into.rows, 'in', options.groups)].sort((a, b) => rank(a) - rank(b));
  sections.push(section('edges', edges, options.limit));
  return {op: 'explain', target, sections};
}

/**
 * A release answers three different questions and used to answer them as one (review 3 #5): what was planned for it
 * (`fix-version`), what it actually carries (its commits and the tickets whose fixes were picked), and what it shipped
 * in features — which only a released, non-dry-run record can support.
 */
async function releaseSections(conn: KuzuConnection, node: Record<string, unknown>, id: string, limit: number): Promise<Section[]> {
  const kind = (k: string) => `MATCH (t)-[e:${quote('TARGETS_RELEASE')}]->(r) WHERE r.${quote('id')} = $id AND e.${quote('kind')} = '${k}' ` +
    `RETURN t.${quote('id')} AS ticket, t.${quote('name')} AS name, t.${quote('state')} AS state, e.${quote('confidence')} AS confidence ORDER BY ticket`;
  const targeted = await run(conn, kind('fix-version'), {id});
  const picked = await run(conn, kind('picked'), {id});
  const commits = await run(conn, `MATCH (r)-[:${quote('INCLUDES')}]->(c) WHERE r.${quote('id')} = $id ` +
    `RETURN c.${quote('id')} AS item, 'includes' AS relation, c.${quote('name')} AS name ORDER BY item`, {id});
  const included = [...commits.rows,
    ...picked.rows.map((r) => ({item: r.ticket, relation: 'picked', name: r.name}))];
  const sections = [
    section('targeted', targeted.rows, limit, 'targeted (0): no ticket carries this release as its fix version'),
    section('included', included, limit, 'included (0): the record names no commit and no pick'),
  ];
  const blocked: string[] = [];
  if (node.dry_run === true) blocked.push('the record is a dry run');
  if (node.state !== 'released') blocked.push(`release state ${node.state ?? 'unknown'}`);
  if (blocked.length) {
    sections.push({title: 'shipped', rows: [], total: 0, empty: `shipped (0): not derivable (${blocked.join('; ')})`});
    return sections;
  }
  const select = `RETURN DISTINCT f.${quote('id')} AS feature, f.${quote('name')} AS name, f.${quote('status')} AS status ORDER BY feature`;
  const affected = await run(conn, `MATCH (t)-[:${quote('TARGETS_RELEASE')}]->(r) WHERE r.${quote('id')} = $id ` +
    `MATCH (t)-[:${quote('AFFECTS')}]->(f:Feature) ${select}`, {id});
  const resolved = await run(conn, `MATCH (r)-[:${quote('INCLUDES')}]->(c)-[:${quote('RESOLVES')}]->(t)-[:${quote('AFFECTS')}]->(f:Feature) ` +
    `WHERE r.${quote('id')} = $id ${select}`, {id});
  const rows = [...affected.rows];
  for (const row of resolved.rows)
    if (!rows.some((r) => r.feature === row.feature)) rows.push(row);
  sections.push(section('shipped', rows, limit, NO_AFFECTS));
  return sections;
}

/** One row per edge type and direction: how many, how they were derived, the evidence of the first, and the first few targets. */
function group(rows: Record<string, unknown>[], direction: 'in' | 'out', groups?: Record<string, string[]>): EdgeGroup[] {
  const byEdge = new Map<string, Record<string, unknown>[]>();
  for (const row of rows) {
    const key = String(row.edge);
    let list = byEdge.get(key);
    if (!list) byEdge.set(key, list = []);
    list.push(row);
  }
  return [...byEdge].sort(([a], [b]) => a < b ? -1 : 1).map(([edge, list]) => {
    const confidences = list.map((r) => Number(r.confidence)).filter((c) => !Number.isNaN(c)).sort((a, b) => a - b);
    const text = (v: unknown) => v === null || v === undefined ? undefined : String(v);
    return {
      group: groupOf(edge, groups), edge, direction, count: list.length,
      derived_by: [...new Set(list.map((r) => String(r.derived_by ?? '')).filter(Boolean))].sort(),
      confidence: confidences.length ? [confidences[0], confidences[confidences.length - 1]] as [number, number] : null,
      evidence: (list.find((r) => Array.isArray(r.evidence) && r.evidence.length)?.evidence as string[] | undefined) ?? [],
      targets: list.slice(0, HOP_TARGETS).map((r) => ({id: String(r.other), name: text(r.name), type: text(r.type)})),
    };
  });
}

/**
 * The vocabulary search. An exact id, name or alias is asked for in its own unbounded query, so that the scan cap on
 * the substring query can no longer drop it (review 3 #9); both sets are then ranked together — authored types first,
 * exact before prefix before substring — and only the ranked whole is paged.
 */
export async function find(conn: KuzuConnection, text: string, options: OpsOptions): Promise<Answer> {
  const q = text.trim().toLowerCase();
  const select = `RETURN n.${quote('id')} AS id, n.${quote('type')} AS type, n.${quote('name')} AS name, n.${quote('status')} AS status`;
  const exact: Record<string, unknown>[] = [];
  const scanned: Record<string, unknown>[] = [];
  for (const table of await nodeTables(conn)) {
    const columns = await columnsOf(conn, table);
    // kuzu 0.11.3 crashes on `any(x IN list ...)` over a column that is null on some rows, and on a parameter
    // inside a lambda, and on a column the table does not have: list_filter, a WITH binding and the real columns.
    const alias = (op: string) => ['aliases', 'keywords'].filter((c) => columns.includes(c))
      .map((c) => `size(list_filter(n.${quote(c)}, x -> lower(x) ${op} q)) > 0`);
    const hit = (op: string) => [...['id', 'name'].filter((c) => columns.includes(c)).map((c) => `lower(n.${quote(c)}) ${op} q`), ...alias(op)].join(' OR ');
    const both = await Promise.all([
      run(conn, `WITH $q AS q MATCH (n:${quote(table)}) WHERE ${hit('=')} ${select}`, {q}),
      run(conn, `WITH $q AS q MATCH (n:${quote(table)}) WHERE ${hit('CONTAINS')}` +
        `${columns.includes('description') ? ` OR lower(n.${quote('description')}) CONTAINS q` : ''} ${select} ORDER BY id LIMIT ${options.scan ?? SEARCH_SCAN}`, {q}),
    ]);
    exact.push(...both[0].rows);
    scanned.push(...both[1].rows);
  }
  const matched = new Set(exact.map((r) => String(r.id)));
  const rank = (row: Record<string, unknown>) => {
    if (matched.has(String(row.id))) return 0;
    const id = String(row.id).toLowerCase(), name = String(row.name ?? '').toLowerCase();
    if (id.startsWith(q) || name.startsWith(q) || id.split('/').pop() === q) return 1;
    return 2;
  };
  const vocabulary = (row: Record<string, unknown>) => {
    const tier = VOCABULARY.findIndex((types) => types.includes(String(row.type)));
    return tier < 0 ? VOCABULARY.length : tier;
  };
  const seen = new Set<string>();
  const rows = [...exact, ...scanned].filter((r) => !seen.has(String(r.id)) && seen.add(String(r.id)))
    .map((row): Record<string, unknown> => ({...row, match: rank(row)}))
    .sort((a, b) => vocabulary(a) - vocabulary(b) || (a.match as number) - (b.match as number) || (String(a.id) < String(b.id) ? -1 : 1));
  return {op: 'find', target: {text}, sections: [section('matches', rows, options.limit)]};
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
