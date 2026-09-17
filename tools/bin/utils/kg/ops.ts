/// The bounded operations agents get before raw Cypher (conventions.md §11.1, build-plan.md WO-7):
/// impact, tests-for, explain and find, each a read-only template over the index, answering with the
/// sections of `answer.ts`.
import * as fs from 'fs';
import * as path from 'path';
import {KuzuConnection, run, quote} from './kuzu';
import {Section, Answer, EdgeGroup, section} from './answer';
import {unitOf, testUnitsOf} from './ids';

export type Tier = 'immediate' | 'reachable' | 'feature';
export const TIERS: Tier[] = ['immediate', 'reachable', 'feature'];
/** `linked` is the quick set: the two tiers the change reaches through structure, without the feature's blast radius. */
export const LINKED: Tier[] = ['immediate', 'reachable'];
export const TIER_CHOICES = `${TIERS.join(', ')}, linked (${LINKED.join(',')}), all or a comma list of the first three`;

export interface OpsOptions {
  limit: number;
  /** The manifest's `edge_groups`, for `explain` to head its edges by folder. */
  groups?: Record<string, string[]>;
  /** `find`: how many substring matches one table may contribute. The exact query is not bounded by it. */
  scan?: number;
  /** `tests-for`: the tiers to answer with; every tier when absent. */
  tiers?: Tier[];
  /** The checkout the `run` rows may look at for a runner's config file; without it the package folder answers. */
  repoRoot?: string;
}

/** One row of the `changes` section: a path the caller asked about, known to the index or not. */
export interface Change {
  path: string;
  known: boolean;
  repo: string;
}

type Chain = string[];
type Target = Record<string, unknown>;

export const DEFAULT_LIMIT = 50;
/** How far the importer walk goes beyond the direct importers (plan.md § Tiers). */
const REACH_DEPTH = 4;
/** An entry file until the packages extractor marks them (`entry`): every test imports it, so the walk stops there. */
const ENTRY = /\/(?:package|package-test|package-api)\.ts$/;
const DEVTOOLS = 'public/packages/DevTools';
const RUNNERS = ['dg', 'xamgle', 'dart', 'playwright', 'node'];
/** The runners whose whole suite is one command in the unit, taken when at least half of its test files are selected. */
const WHOLE_SUITE = ['dg', 'dart', 'playwright', 'node'];
const WHOLE_NAMES = 3;
/** DevTools lists a client test under `Core: <segments of its category>`. */
const CORE_CATEGORY = 'Core';
const DECLARED = /^(?:decl|func|ep):/;
/** How many levels of barrels (a file that only re-exports) a target is expanded through. */
const BARREL_DEPTH = 3;
/** A file of a test folder that declares no test is a helper; a test importing it reaches what the helper uses (§8.1). */
const TEST_FOLDER = /(?:^|\/)(?:test|tests|__tests__)\//;
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

/** `--tier` as written: a tier, `linked`, `all` or a comma list of tiers, in tier order; nothing when a word is not a tier. */
export function parseTiers(text: string): Tier[] | undefined {
  const words = text.split(',').map((s) => s.trim().toLowerCase()).filter(Boolean);
  if (!words.length) return undefined;
  if (words.length === 1 && words[0] === 'all') return [...TIERS];
  if (words.length === 1 && words[0] === 'linked') return [...LINKED];
  if (words.some((w) => !TIERS.includes(w as Tier))) return undefined;
  return TIERS.filter((t) => words.includes(t));
}

/** Ids may be written with or without the sigil; a path names the source file it belongs to. */
function bareId(arg: string): string {
  return arg.trim().replace(/^~/, '');
}

export async function resolveTarget(conn: KuzuConnection, arg: string): Promise<Target | null> {
  return (await resolveTargets(conn, [arg])).get(arg) ?? null;
}

/** Every argument at once: the id as written, then the file and the doc a path names, the first that exists. */
export async function resolveTargets(conn: KuzuConnection, args: string[]): Promise<Map<string, Target | null>> {
  const candidates = (arg: string) => {
    const id = bareId(arg);
    const p = id.replace(/\\/g, '/');
    return [id, `file:${p}`, `doc:${p}`];
  };
  const {rows} = await run(conn, `MATCH (n) WHERE n.${quote('id')} IN $ids ` +
    `RETURN n.${quote('id')} AS id, label(n) AS root, n.${quote('type')} AS type, n.${quote('name')} AS name, n.${quote('status')} AS status`,
  {ids: [...new Set(args.flatMap(candidates))]});
  const found = new Map(rows.map((r) => [String(r.id), r]));
  return new Map(args.map((arg) => [arg, candidates(arg).map((c) => found.get(c)).find((r) => r) ?? null]));
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

/**
 * What a change to this file, declaration or feature reaches: the features that own it, their evidence and their
 * work, the documents citing it (immediate) and documenting its features, and who imports or calls it — the direct
 * importers of a file and everything that reaches it through them (plan.md § Tiers).
 */
export async function impact(conn: KuzuConnection, target: Record<string, unknown>, options: OpsOptions): Promise<Answer> {
  const id = String(target.id);
  const sections: Section[] = [];
  const features = await featuresOf(conn, target);
  sections.push(section('features', features, options.limit, NODE_PATH.test(id) ? UNOWNED : undefined));
  const ids = features.map((f) => String(f.feature));
  const contents = [...(await contained(conn, target)).keys()];
  const cites = await run(conn, `MATCH (d)-[:${quote('MENTIONS')}]->(n) WHERE n.${quote('id')} IN $ids ` +
    `RETURN d.${quote('id')} AS document, d.${quote('type')} AS type, n.${quote('id')} AS target ORDER BY document, target`, {ids: contents});
  if (ids.length) {
    const owners = await run(conn, `MATCH (f:Feature)-[:${quote('owner')}]->(p) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, p.${quote('id')} AS owner, p.${quote('name')} AS name`, {ids});
    sections.push(section('owners', owners.rows, options.limit));
    sections.push(section('evidence', await evidence(conn, ids), options.limit));
    sections.push(section('cites', cites.rows, options.limit));
    const docs = await run(conn, `MATCH (d)-[:${quote('DOCUMENTS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, d.${quote('id')} AS document, d.${quote('type')} AS type`, {ids});
    sections.push(section('documents', docs.rows, options.limit));
    const work = await run(conn, `MATCH (t)-[:${quote('AFFECTS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, t.${quote('id')} AS ticket, 'affects' AS relation, t.${quote('state')} AS state`, {ids});
    const tracked = await run(conn, `MATCH (f:Feature)-[:${quote('TRACKED_IN')}]->(t) WHERE f.${quote('id')} IN $ids ` +
      `RETURN f.${quote('id')} AS feature, t.${quote('id')} AS ticket, 'tracked-in' AS relation, t.${quote('state')} AS state`, {ids});
    sections.push(section('work', [...work.rows, ...tracked.rows], options.limit));
  }
  else {
    sections.push(await packageOwners(conn, id, options.limit));
    sections.push(section('cites', cites.rows, options.limit));
  }
  const answer: Answer = {op: 'impact', target, sections};
  if (id.startsWith('file:')) {
    const files = new Map<string, Chain>([[id, [id]]]);
    const notes = await barrels(conn, files);
    if (notes.length) answer.notes = notes;
    const {direct, reachable} = await importers(conn, files);
    const users = await run(conn, `MATCH (c)-[e:${quote('USES')}]->(d)<-[:${quote('DECLARES')}]-(n) WHERE n.${quote('id')} = $id ` +
      `RETURN c.${quote('id')} AS caller, 'USES' AS via, d.${quote('id')} AS target, e.${quote('count')} AS count ORDER BY count DESC`, {id});
    sections.push(section('importers', [...[...direct.keys()].map((caller) => ({caller, via: 'IMPORTS', target: id, count: null})), ...users.rows], options.limit));
    sections.push(section('reachable', [...reachable].map(([importer, chain]) => ({importer, via: chain.join(' ')})), options.limit));
  }
  else if (DECLARED.test(id)) {
    const callers = await run(conn, `MATCH (c)-[e:${quote('USES')}|${quote('CALLS')}]->(n) WHERE n.${quote('id')} = $id ` +
      `RETURN c.${quote('id')} AS caller, label(e) AS via, n.${quote('id')} AS target, e.${quote('count')} AS count ORDER BY count DESC`, {id});
    sections.push(section('callers', callers.rows, options.limit));
  }
  return answer;
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

/** The source files the targets stand for, each with the chain from its target: a file is itself, a declaration
 * is the file that declares it. */
async function sourceFiles(conn: KuzuConnection, targets: Target[]): Promise<Map<string, Chain>> {
  const files = new Map<string, Chain>();
  const declared: string[] = [];
  for (const t of targets) {
    const id = String(t.id);
    if (id.startsWith('file:')) files.set(id, [id]);
    else if (DECLARED.test(id)) declared.push(id);
  }
  if (declared.length) {
    const {rows} = await run(conn, `MATCH (f)-[:${quote('DECLARES')}]->(d) WHERE d.${quote('id')} IN $ids AND f.${quote('id')} STARTS WITH 'file:' ` +
      `RETURN d.${quote('id')} AS decl, f.${quote('id')} AS file ORDER BY decl, file`, {ids: declared});
    for (const r of rows)
      if (!files.has(String(r.file))) files.set(String(r.file), [String(r.decl), hop('declares', 'in'), String(r.file)]);
  }
  return files;
}

/** A barrel — a file that declares nothing and only re-exports — stands for the files it re-exports, and a barrel of
 * barrels (`js-api/src/dataframe.ts` → `dataframe/index.ts` → the classes) for theirs, BARREL_DEPTH levels down; they
 * join [files] with the barrel's chain, and one note per barrel says how many. */
async function barrels(conn: KuzuConnection, files: Map<string, Chain>): Promise<string[]> {
  if (!files.size || !(await columnsOf(conn, 'IMPORTS')).includes('reexport')) return [];
  const rootOf = new Map<string, string>();
  const added = new Map<string, number>();
  let frontier = [...files.keys()];
  for (let depth = 0; depth < BARREL_DEPTH && frontier.length; depth++) {
    const declaring = new Set((await run(conn, `MATCH (b)-[:${quote('DECLARES')}]->(d) WHERE b.${quote('id')} IN $ids AND d.${quote('id')} STARTS WITH 'decl:' ` +
      `RETURN DISTINCT b.${quote('id')} AS barrel`, {ids: frontier})).rows.map((r) => String(r.barrel)));
    const {rows} = await run(conn, `MATCH (b)-[e:${quote('IMPORTS')}]->(m) WHERE b.${quote('id')} IN $ids ` +
      `RETURN b.${quote('id')} AS barrel, m.${quote('id')} AS target, e.${quote('reexport')} AS reexport ORDER BY barrel, target`, {ids: frontier});
    const exports = new Map<string, {all: boolean, targets: string[]}>();
    for (const r of rows) {
      const barrel = String(r.barrel);
      if (declaring.has(barrel)) continue;
      const entry = exports.get(barrel) ?? {all: true, targets: []};
      entry.all = entry.all && r.reexport === true;
      entry.targets.push(String(r.target));
      exports.set(barrel, entry);
    }
    frontier = [];
    for (const [barrel, {all, targets}] of exports) {
      if (!all) continue;
      const root = rootOf.get(barrel) ?? barrel;
      for (const t of targets.filter((t) => t.startsWith('file:') && !files.has(t))) {
        files.set(t, [...files.get(barrel)!, hop('imports', 'out'), t]);
        rootOf.set(t, root);
        added.set(root, (added.get(root) ?? 0) + 1);
        frontier.push(t);
      }
    }
  }
  return [...added].map(([barrel, n]) => `barrel: ${barrel.slice('file:'.length)} expanded to ${n} re-exported file${n === 1 ? '' : 's'}`);
}

/**
 * Who imports these files, hop by hop: the direct importers, then everything that reaches them within REACH_DEPTH
 * hops without passing through an entry file (plan.md § Tiers). Breadth first, so the chain kept is the shortest.
 */
async function importers(conn: KuzuConnection, files: Map<string, Chain>): Promise<{direct: Map<string, Chain>, reachable: Map<string, Chain>}> {
  const flagged = (await columnsOf(conn, 'Component')).includes('entry');
  const direct = new Map<string, Chain>();
  const reachable = new Map<string, Chain>();
  const seen = new Set(files.keys());
  let frontier = files;
  for (let depth = 1; depth <= REACH_DEPTH && frontier.size; depth++) {
    const {rows} = await run(conn, `MATCH (c)-[:${quote('IMPORTS')}]->(n) WHERE n.${quote('id')} IN $ids ` +
      `RETURN c.${quote('id')} AS importer, n.${quote('id')} AS imported${flagged ? `, c.${quote('entry')} AS entry` : ''} ORDER BY imported, importer`,
    {ids: [...frontier.keys()]});
    const next = new Map<string, Chain>();
    for (const r of rows) {
      const id = String(r.importer);
      if (seen.has(id)) continue;
      seen.add(id);
      const chain = [...frontier.get(String(r.imported))!, hop('imports', 'in'), id];
      (depth === 1 ? direct : reachable).set(id, chain);
      if (!(flagged ? r.entry === true : ENTRY.test(id))) next.set(id, chain);
    }
    frontier = next;
  }
  return {direct, reachable};
}

/** The test files named after the targets' files: the `mirrors` edge the build draws (conventions.md §8.1). */
async function mirrors(conn: KuzuConnection, files: Map<string, Chain>): Promise<Map<string, Chain>> {
  const found = new Map<string, Chain>();
  if (!files.size) return found;
  const {rows} = await run(conn, `MATCH (t)-[:${quote('MIRRORS')}]->(f) WHERE f.${quote('id')} IN $ids ` +
    `RETURN t.${quote('id')} AS test, f.${quote('id')} AS file ORDER BY file, test`, {ids: [...files.keys()]});
  for (const r of rows)
    if (!found.has(String(r.test))) found.set(String(r.test), [...files.get(String(r.file))!, hop('mirrors', 'in'), String(r.test)]);
  return found;
}

/** The files using what the targets declare: a declaration target directly, a file target through its declarations. A
 * test reaches a use through a helper it imports; a helper most of the unit's test files import (the setup of a suite)
 * is a hub, and what it reaches is returned apart, for the reachable tier. */
async function users(conn: KuzuConnection, targets: Target[], files: Map<string, Chain>): Promise<{found: Map<string, Chain>, hub: Map<string, Chain>}> {
  const found = new Map<string, Chain>();
  const hub = new Map<string, Chain>();
  const declared = targets.map((t) => String(t.id)).filter((id) => DECLARED.test(id));
  if (declared.length) {
    const {rows} = await run(conn, `MATCH (c)-[:${quote('USES')}]->(d) WHERE d.${quote('id')} IN $ids AND c.${quote('id')} STARTS WITH 'file:' ` +
      `RETURN c.${quote('id')} AS user, d.${quote('id')} AS decl ORDER BY decl, user`, {ids: declared});
    for (const r of rows) found.set(String(r.user), [String(r.decl), hop('uses', 'in'), String(r.user)]);
  }
  const fileIds = [...files.keys()].filter((id) => !DECLARED.test(files.get(id)![0]));
  if (fileIds.length) {
    const {rows} = await run(conn, `MATCH (c)-[:${quote('USES')}]->(d)<-[:${quote('DECLARES')}]-(n) WHERE n.${quote('id')} IN $ids AND c.${quote('id')} STARTS WITH 'file:' ` +
      `RETURN DISTINCT c.${quote('id')} AS user, d.${quote('id')} AS decl, n.${quote('id')} AS file ORDER BY file, decl, user`, {ids: fileIds});
    for (const r of rows)
      if (!found.has(String(r.user))) found.set(String(r.user), [...files.get(String(r.file))!, hop('declares', 'out'), String(r.decl), hop('uses', 'in'), String(r.user)]);
  }
  const helpers = [...found.keys()].filter((id) => TEST_FOLDER.test(id));
  if (helpers.length) {
    const {rows} = await run(conn, `MATCH (t)-[:${quote('IMPORTS')}]->(h) WHERE h.${quote('id')} IN $ids ` +
      `OPTIONAL MATCH (h)-[:${quote('DECLARES')}]->(x:Artifact) RETURN t.${quote('id')} AS test, h.${quote('id')} AS helper, count(x) AS declared ORDER BY helper, test`, {ids: helpers});
    const byHelper = new Map<string, string[]>();
    for (const r of rows) {
      const test = String(r.test), helper = String(r.helper);
      if (Number(r.declared) || found.has(test) || unitOf(test.slice('file:'.length)) !== unitOf(helper.slice('file:'.length))) continue;
      byHelper.set(helper, [...byHelper.get(helper) ?? [], test]);
    }
    const sizes = await suiteSizes(conn);
    for (const [helper, tests] of byHelper) {
      const unit = unitOf(helper.slice('file:'.length));
      const total = [...sizes].filter(([k]) => k.endsWith(` ${unit}`)).reduce((n, [, v]) => n + v, 0);
      const into = tests.length >= HUB_MIN && tests.length * 2 >= total ? hub : found;
      for (const test of tests) into.set(test, [...found.get(helper)!, hop('imports', 'in'), test]);
    }
  }
  return {found, hub};
}

const HUB_MIN = 3;

const TEST_COLUMNS = `t.${quote('framework')} AS framework, t.${quote('level')} AS level, t.${quote('id')} AS test, t.${quote('name')} AS name, ` +
  `t.${quote('path')} AS path, t.${quote('category')} AS category`;

/** The unit of the changed file a chain starts from: its first `file:` node (a declaration target leads with its file). */
function unitOfChain(chain: Chain): string | undefined {
  const file = chain.find((s) => s.startsWith('file:'));
  return file === undefined ? undefined : unitOf(file.slice('file:'.length));
}

/** The tests declared in these files, one row each with the file's chain, the feature it tests when it tests one. Without
 * a [tier] the link decides: a test file standing for the unit of the changed file (`testUnitsOf`) is immediate, any other
 * reachable. */
async function testsIn(conn: KuzuConnection, files: Map<string, Chain>, tier?: Tier): Promise<Record<string, unknown>[]> {
  if (!files.size) return [];
  const paths = [...files.keys()].map((id) => id.slice('file:'.length));
  const {rows} = await run(conn, `MATCH (t:Artifact) WHERE t.${quote('type')} = 'test' AND t.${quote('path')} IN $paths ` +
    `OPTIONAL MATCH (t)-[:${quote('TESTS')}]->(f:Feature) RETURN ${TEST_COLUMNS}, f.${quote('id')} AS feature, ` +
    `t.${quote('skipped')} AS skipped, t.${quote('dynamic')} AS dynamic ORDER BY path, test`, {paths});
  return rows.map((r) => {
    const chain = files.get(`file:${r.path}`)!;
    const unit = unitOfChain(chain);
    const path = String(r.path);
    const same = unit === undefined ? unitOf(path) === undefined : testUnitsOf(path, r.category === null ? undefined : String(r.category)).includes(unit);
    return {...r, feature: r.feature ?? null, via: chain.join(' '), tier: tier ?? (same ? 'immediate' : 'reachable')};
  });
}

/** The tests the changed files declare through the `declares` edge, once the build emits it (plan.md § Graph changes);
 * until then the table has no file → test pair and the path match in `testsIn` is the whole answer. */
async function declaredTests(conn: KuzuConnection, files: Map<string, Chain>): Promise<Record<string, unknown>[]> {
  const pairs = await run(conn, `CALL show_connection('DECLARES') RETURN *`);
  if (!files.size || !pairs.rows.some((r) => r['source table name'] === 'Component' && r['destination table name'] === 'Artifact')) return [];
  const {rows} = await run(conn, `MATCH (n)-[:${quote('DECLARES')}]->(t:Artifact) WHERE n.${quote('id')} IN $ids AND t.${quote('type')} = 'test' ` +
    `OPTIONAL MATCH (t)-[:${quote('TESTS')}]->(f:Feature) RETURN n.${quote('id')} AS file, ${TEST_COLUMNS}, f.${quote('id')} AS feature, ` +
    `t.${quote('skipped')} AS skipped, t.${quote('dynamic')} AS dynamic ORDER BY file, test`, {ids: [...files.keys()]});
  return rows.map(({file, ...r}) => ({...r, feature: r.feature ?? null, via: [...files.get(String(file))!, hop('declares', 'out'), String(r.test)].join(' '), tier: 'immediate'}));
}

/** The immediate and reachable tiers of these targets (plan.md § Tiers): the tests in the source files, in the files
 * importing them, in their mirror test files and in the files using what they declare — immediate in the unit of the
 * changed file, reachable from another unit — then the tests reached through the import walk; a test met twice keeps
 * its first chain unless a later link lifts it to immediate. `seen` holds every test id met, so the feature tier can
 * leave them out. */
export async function linkedTests(conn: KuzuConnection, targets: Target[], seen = new Set<string>()):
  Promise<{files: Map<string, Chain>, immediate: Record<string, unknown>[], reachable: Record<string, unknown>[], notes: string[]}> {
  const files = await sourceFiles(conn, targets);
  const notes = await barrels(conn, files);
  const {direct, reachable} = await importers(conn, files);
  const linked = new Map<string, Record<string, unknown>>();
  const used = await users(conn, targets, files);
  for (const row of [...await declaredTests(conn, files), ...await testsIn(conn, files), ...await testsIn(conn, direct),
    ...await testsIn(conn, await mirrors(conn, files)), ...await testsIn(conn, used.found), ...await testsIn(conn, used.hub, 'reachable')]) {
    const best = linked.get(String(row.test));
    if (!best || (best.tier === 'reachable' && row.tier === 'immediate')) linked.set(String(row.test), row);
  }
  for (const test of linked.keys()) seen.add(test);
  const walked = (await testsIn(conn, reachable, 'reachable')).filter((r) => !seen.has(String(r.test)));
  for (const r of walked) seen.add(String(r.test));
  const rows = [...linked.values()];
  return {files, immediate: rows.filter((r) => r.tier === 'immediate'), reachable: [...rows.filter((r) => r.tier === 'reachable'), ...walked], notes};
}

/** How many test files each unit holds per framework, so the run rows know when a selection is most of a suite. */
async function suiteSizes(conn: KuzuConnection): Promise<Map<string, number>> {
  const {rows} = await run(conn, `MATCH (t:Artifact) WHERE t.${quote('type')} = 'test' RETURN DISTINCT t.${quote('framework')} AS framework, t.${quote('path')} AS path`);
  const out = new Map<string, number>();
  for (const r of rows) {
    const unit = unitOf(String(r.path));
    if (unit) out.set(`${r.framework} ${unit}`, (out.get(`${r.framework} ${unit}`) ?? 0) + 1);
  }
  return out;
}

/**
 * The tests for one target or several (plan.md § Tiers): the immediate tier — tests in the changed files, in the
 * files importing them, in their mirror test files and in the files using what they declare — then the reachable
 * tier through the import walk, then everything the owning features and their subtrees carry, with the scenarios
 * and automations of those features, and one `run` row per runner invocation over the tests the `tiers` select.
 */
export async function testsFor(conn: KuzuConnection, target: Target | Target[], options: OpsOptions, changes?: Change[]): Promise<Answer> {
  const targets = Array.isArray(target) ? target : [target];
  const tiers = options.tiers ?? TIERS;
  const wanted = (t: Tier) => tiers.includes(t);
  const sections: Section[] = [];
  if (changes) sections.push(section('changes', changes, options.limit));
  const byFeature = new Map<string, Record<string, unknown>>();
  for (const t of targets)
    for (const row of await featuresOf(conn, t)) {
      const best = byFeature.get(String(row.feature));
      if (!best || stronger(row, best)) byFeature.set(String(row.feature), row);
    }
  const found = [...byFeature.values()];
  const roots = found.map((r) => String(r.feature));
  const unowned = targets.length === 1 && NODE_PATH.test(String(targets[0].id)) ? UNOWNED : undefined;
  const seen = new Set<string>();
  const {files, immediate, reachable: reached, notes} = await linkedTests(conn, targets, seen);
  const none = files.size ? undefined : ' (0): none of the targets is a source file or a declaration';
  const featureTier: Record<string, unknown>[] = [];
  const scenarios: Record<string, unknown>[] = [];
  const automations: Record<string, unknown>[] = [];
  const features = new Map<string, Record<string, unknown>>();
  if (roots.length) {
    const via = new Map(found.map((r) => [String(r.feature), r]));
    const descendants = await run(conn, `MATCH (d:Feature)-[:${quote('PART_OF')}*0..5]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
      `RETURN DISTINCT d.${quote('id')} AS feature, d.${quote('name')} AS name, d.${quote('status')} AS status, f.${quote('id')} AS root`, {ids: roots});
    for (const r of descendants.rows) {
      if (features.has(String(r.feature))) continue;
      const root = via.get(String(r.root))!;
      const path = r.feature === r.root ? root.path as string[] : [...root.path as string[], hop('part-of', 'in'), String(r.feature)];
      features.set(String(r.feature), {feature: r.feature, name: r.name, status: r.status, via: path.join(' '), path});
    }
    const ids = [...features.keys()];
    const tests = await run(conn, `MATCH (t)-[:${quote('TESTS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
      `RETURN ${TEST_COLUMNS}, f.${quote('id')} AS feature, t.${quote('skipped')} AS skipped, t.${quote('dynamic')} AS dynamic ORDER BY framework, test`, {ids});
    for (const r of tests.rows) {
      if (seen.has(String(r.test))) continue;
      seen.add(String(r.test));
      featureTier.push({...r, via: [...features.get(String(r.feature))!.path as string[], hop('tests', 'in'), String(r.test)].join(' '), tier: 'feature'});
    }
    const covers = await run(conn, `MATCH (s)-[:${quote('COVERS')}]->(f:Feature) WHERE f.${quote('id')} IN $ids ` +
      `RETURN s.${quote('id')} AS scenario, f.${quote('id')} AS feature, s.${quote('manual_only')} AS manual_only, s.${quote('priority')} AS priority ` +
      `ORDER BY scenario`, {ids});
    scenarios.push(...covers.rows);
    const covered = scenarios.map((r) => String(r.scenario));
    if (covered.length) automations.push(...(await run(conn, `MATCH (t)-[:${quote('AUTOMATES')}]->(s) WHERE s.${quote('id')} IN $ids ` +
      `RETURN t.${quote('framework')} AS framework, t.${quote('id')} AS test, s.${quote('id')} AS scenario ORDER BY framework, test`, {ids: covered})).rows);
  }
  sections.push(section('features', [...features.values()], options.limit, unowned));
  if (!roots.length && targets.length === 1) sections.push(await packageOwners(conn, String(targets[0].id), options.limit));
  if (wanted('immediate')) sections.push(section('immediate', immediate, options.limit, none && `immediate${none}`));
  if (wanted('reachable')) sections.push(section('reachable', reached, options.limit, none && `reachable${none}`));
  if (wanted('feature')) {
    sections.push(section('feature', featureTier, options.limit, !roots.length ? `feature (0): no feature owns the target${targets.length === 1 ? '' : 's'}`
      : `feature (0): ${seen.size ? 'every test of' : 'no test carries'} ~${roots.join(', ~')}${seen.size ? ' is already in the immediate or reachable tier' : ' and no owned file contains tests'}`));
    sections.push(section('scenarios', scenarios, options.limit));
    sections.push(section('automations', automations, options.limit));
  }
  const selected = [...(wanted('immediate') ? immediate : []), ...(wanted('reachable') ? reached : []), ...(wanted('feature') ? featureTier : [])];
  const runs = runRows(selected, options.repoRoot, selected.length ? await suiteSizes(conn) : undefined);
  sections.push(section('run', runs.rows, options.limit));
  const answer: Answer = Array.isArray(target) ? {op: 'tests-for', targets, sections} : {op: 'tests-for', target, sections};
  if (notes.length || runs.notes.length) answer.notes = [...notes, ...runs.notes];
  return answer;
}

/** A shell argument: double-quoted, with the quotes inside escaped. */
function arg(text: string): string {
  return `"${text.replace(/"/g, '\\"')}"`;
}

/** The folder holding [file] up from the test's folder, inside [repoRoot]; none when the checkout is not at hand. */
function nearest(repoRoot: string | undefined, testPath: string, file: string): string | undefined {
  if (!repoRoot) return undefined;
  for (let dir = path.posix.dirname(testPath); dir !== '.' && dir !== '/'; dir = path.posix.dirname(dir))
    if (fs.existsSync(path.join(repoRoot, dir, file))) return dir;
  return undefined;
}

/**
 * One row per runner invocation over the selected tests (plan.md § `grok test --recent`): the DG runner takes one
 * category with one `--test` name or the whole category, the Dart runner (`pub run test`: core is Dart 1.x, which
 * has no `dart test`) takes files and one `-n`, Playwright its spec files and one `--grep`, vitest its files. A client
 * test runs in DevTools under `Core: <its category segments>`. A dynamic test has no runnable name: a DG category runs
 * whole for it, a client one skips it and says so in the row. When [suites] says a unit holds at most twice as many test files of a
 * framework as the selection, the rows of that unit collapse into its whole suite and the row says so.
 */
export function runRows(tests: Record<string, unknown>[], repoRoot?: string, suites?: Map<string, number>): {rows: Record<string, unknown>[], notes: string[]} {
  const groups = new Map<string, {framework: string, cwd: string, category: string, files: string[], names: string[], dynamic: boolean, skipped: number, whole?: string}>();
  const notes = new Set<string>();
  const selected = new Map<string, Set<string>>();
  for (const t of tests) {
    const unit = unitOf(String(t.path));
    if (unit && WHOLE_SUITE.includes(String(t.framework))) selected.set(`${t.framework} ${unit}`, (selected.get(`${t.framework} ${unit}`) ?? new Set()).add(String(t.path)));
  }
  const whole = new Map<string, string>();
  for (const [key, files] of selected) {
    const total = suites?.get(key);
    if (total && files.size * 2 >= total) whole.set(key, `whole suite: ${files.size} of ${total} test files selected`);
  }
  for (const t of tests) {
    const framework = String(t.framework);
    const p = String(t.path);
    let cwd: string | undefined;
    let category = '';
    if (framework === 'dg') {
      cwd = unitOf(p);
      category = String(t.category ?? '');
    }
    else if (framework === 'xamgle') {
      cwd = DEVTOOLS;
      category = [CORE_CATEGORY, ...String(t.category ?? '').split('|').map((s) => s.trim()).filter(Boolean)].join(': ');
    }
    else if (framework === 'dart') cwd = unitOf(p);
    else if (framework === 'playwright') cwd = nearest(repoRoot, p, 'playwright.config.ts') ?? unitOf(p);
    else if (framework === 'node') cwd = nearest(repoRoot, p, 'package.json') ?? unitOf(p);
    if (cwd === undefined) notes.add(RUNNERS.includes(framework) ? `run: no package folder for ${p}` : `run: no runner known for framework '${framework}'`);
    const suite = whole.get(`${framework} ${unitOf(p)}`);
    const key = [framework, cwd ?? '', suite ? '' : category].join(' ');
    let group = groups.get(key);
    if (!group) groups.set(key, group = {framework, cwd: cwd ?? '', category, files: [], names: [], dynamic: false, skipped: 0, whole: suite});
    if (framework === 'xamgle' && t.dynamic === true) {
      group.skipped++;
      continue;
    }
    const file = cwd === undefined ? p : path.posix.relative(cwd, p);
    if (!group.files.includes(file)) group.files.push(file);
    group.names.push(String(t.name));
    if (t.dynamic === true) group.dynamic = true;
  }
  const rows: Record<string, unknown>[] = [];
  for (const g of groups.values()) {
    const note = g.whole ?? (g.skipped ? `${g.skipped} dynamic test${g.skipped === 1 ? '' : 's'} skipped: registered in a loop, no runnable name` : undefined);
    if (!g.names.length) {
      notes.add(`run: ${g.category}: ${note}`);
      continue;
    }
    const one = g.names.length === 1 && !g.dynamic && !g.whole ? g.names[0] : undefined;
    const files = g.whole ? '' : ` ${g.files.join(' ')}`;
    const command = !g.cwd ? '?'
      : g.framework === 'dart' ? `pub run test${files}${one ? ` -n ${arg(one)}` : ''}`
        : g.framework === 'playwright' ? `npx playwright test${files}${one ? ` --grep ${arg(one)}` : ''}`
          : g.framework === 'node' ? `npx vitest run${files}`
            : g.whole ? 'grok test' : `grok test --category ${arg(g.category)}${one ? ` --test ${arg(one)}` : ''}`;
    const names = g.whole ? g.names.slice(0, WHOLE_NAMES) : g.names;
    rows.push(note ? {framework: g.framework, cwd: g.cwd, command, tests: g.names.length, names, note}
      : {framework: g.framework, cwd: g.cwd, command, tests: g.names.length, names});
  }
  return {rows, notes: [...notes]};
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
