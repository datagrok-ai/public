/// The one seam between the canonical JSONL and the graph index (conventions.md §11.1, build-plan.md
/// WO-7): the optional `kuzu` binding, the DDL derived from the type system, and the bulk load through
/// temporary CSVs. Nothing else in `grok kg` may import `kuzu`, and nothing here may be needed to build:
/// when the binding is absent `loadKuzu()` returns null and the caller says so.
///
/// One rule the 0.11.3 binding imposes on every caller, and it takes the whole process down: **close every
/// `QueryResult` before the connection and the database**. One left open when the database closes segfaults
/// at native teardown — on win32-x64 `kg build` wrote all its output and then exited 139 — so `run` and
/// `exec` close theirs, and nothing else may call `conn.query` directly. It is also what killed a vitest
/// worker that read the type system while a database was open; with the results closed, that order is safe.
import * as fs from 'fs';
import * as path from 'path';
import {createRequire} from 'module';
import {TypeSystem, Member, NodeType, pascal, graphLabel} from './types';

export interface KuzuQueryResult {
  getAll(): Promise<Record<string, unknown>[]>;
  getColumnNames(): Promise<string[]>;
  close(): void;
}

export interface KuzuPreparedStatement {
  isSuccess(): boolean;
  getErrorMessage(): string;
}

export interface KuzuConnection {
  query(statement: string): Promise<KuzuQueryResult | KuzuQueryResult[]>;
  prepare(statement: string): Promise<KuzuPreparedStatement>;
  execute(prepared: KuzuPreparedStatement, params: Record<string, unknown>): Promise<KuzuQueryResult>;
  close(): Promise<void>;
}

export interface KuzuDatabase {
  close(): Promise<void>;
}

/** What `grok kg` uses of the binding; the local shape keeps `kuzu` out of the type graph. */
export interface KuzuModule {
  VERSION: string;
  Database: new (databasePath?: string, bufferManagerSize?: number, enableCompression?: boolean, readOnly?: boolean) => KuzuDatabase;
  Connection: new (database: KuzuDatabase, numThreads?: number) => KuzuConnection;
}

export const KUZU_VERSION = '0.11.3';
export const MISSING_KUZU = `kuzu ${KUZU_VERSION} is not installed; the JSONL under .kg/data is complete. ` +
  `Install with: npm install kuzu@${KUZU_VERSION}`;

/** Kuzu identifiers are case-insensitive, so one lowercase name is one table. */
export interface Column {
  name: string;
  type: string;
}

export interface NodeTable {
  name: string;
  root: string;
  columns: Column[];
  /** Concrete node types stored here, in `nodes/<type>.jsonl` order. */
  types: string[];
}

export interface RelTable {
  name: string;
  pairs: {from: string, to: string}[];
  columns: Column[];
  /** `data/edges/<name>.jsonl` files that feed this table: an edge type name, a reference property name, or both. */
  groups: string[];
}

export interface Ddl {
  nodes: NodeTable[];
  rels: RelTable[];
  statements: string[];
}

export interface TableRows {
  table: string;
  rows: number;
}

export interface LoadResult {
  db: string;
  ms: number;
  bytes: number;
  nodes: TableRows[];
  rels: TableRows[];
  /** Rows whose list values CSV cannot carry (see `LIST_UNSAFE`), inserted one by one instead. */
  parameterized: number;
  /** Peak resident memory of this process during the load, in MB; recorded in the manifest as telemetry. */
  memoryMb: number;
  /** `<platform>-<arch>` the index was written on; Kuzu's format is not portable between them. */
  platform: string;
}

const BUILD_ORDER = ['id', 'type', 'types'];
const EDGE_HEAD = ['from', 'to'];
/** Joins a root pair into one map key; no table name contains it. */
const SEPARATOR = '\u0000';
/** A list item CSV can carry: Kuzu's list reader splits on `,`, honours `[](){}` and quotes, and unescapes nothing. */
const LIST_UNSAFE = /^$|^["'\s]|[\s]$|[,[\]{}]/;

/** Kuzu takes 80% of free memory by default; a CLI over a 100 MB graph needs a fraction of that. The load
 * needs room for the tables it builds, a reader only for what it touches — `--memory <MB>` and
 * `KG_KUZU_MEMORY` override both. */
export const BUILD_MEMORY_MB = 2048;
export const READ_MEMORY_MB = 512;
/** Graph identifiers reach Cypher inside backticks; nothing else may (conventions.md §7.1). */
const IDENTIFIER = /^[A-Za-z_][A-Za-z0-9_]*$/;
/** How often the load samples its own resident memory for the high-water mark. */
const SAMPLE_MS = 250;

const nodeRequire = typeof require === 'function' ? require : createRequire(path.join(process.cwd(), 'grok.js'));

/** The binding, or null when it is not installed: the index is optional, the JSONL is not. */
export function loadKuzu(): KuzuModule | null {
  try {
    return nodeRequire('kuzu') as KuzuModule;
  }
  catch {
    return null;
  }
}

/** The index schema: one node table per root, one rel table per concrete edge type and per reference property. */
export function ddl(system: TypeSystem): Ddl {
  const nodes = nodeTables(system);
  const rels = relTables(system);
  const seen = new Map<string, string>();
  for (const name of [...nodes.map((t) => t.name), ...rels.map((t) => t.name)]) {
    const first = seen.get(name.toLowerCase());
    if (first) throw new Error(`kuzu tables '${first}' and '${name}' differ only in case; rename one of them (conventions.md §7.5)`);
    seen.set(name.toLowerCase(), name);
  }
  const statements = [
    ...nodes.map((t) => `CREATE NODE TABLE ${quote(t.name)}(${t.columns.map((c) => `${column(c)}${c.name === 'id' ? ' PRIMARY KEY' : ''}`).join(', ')})`),
    ...rels.map((t) => `CREATE REL TABLE ${quote(t.name)}(${t.pairs.map((p) => `FROM ${quote(p.from)} TO ${quote(p.to)}`).join(', ')}` +
      `${t.columns.length ? `, ${t.columns.map(column).join(', ')}` : ''})`),
  ];
  return {nodes, rels, statements};
}

function nodeTables(system: TypeSystem): NodeTable[] {
  const tables: NodeTable[] = [];
  for (const root of system.roots) {
    const subtree = [...system.nodes.values()].filter((n) => n.root === root);
    const columns = new Map<string, Column>(BUILD_ORDER.map((n) => [n, {name: n, type: n === 'types' ? 'STRING[]' : 'STRING'}]));
    const members = new Map<string, {member: Member, type: NodeType}>();
    for (const type of subtree)
      for (const member of Object.values(type.members)) {
        const type2 = columnType(member);
        const held = members.get(member.name);
        if (!held) {
          members.set(member.name, {member, type});
          columns.set(member.name, {name: member.name, type: type2});
          continue;
        }
        const merged = widen(columns.get(member.name)!.type, type2);
        if (!merged)
          throw new Error(`kuzu table ${pascal(root)}: ${held.type.name} declares ${member.name} as ${held.member.spec}, ` +
            `${type.name} as ${member.spec}; they cannot share a column`);
        columns.get(member.name)!.type = merged;
      }
    for (const m of system.buildFields.node) columns.set(m.name, {name: m.name, type: columnType(m)});
    tables.push({name: pascal(root), root, columns: order(columns, system), types: subtree.filter((t) => !t.abstract).map((t) => t.name)});
  }
  return tables;
}

/** id, type, types, the base fields in their reserved order, then the subtree's own members and the build fields. */
function order(columns: Map<string, Column>, system: TypeSystem): Column[] {
  const head = [...BUILD_ORDER, ...system.reservedNodeFields, ...system.buildFields.node.map((m) => m.name)];
  const rest = [...columns.keys()].filter((n) => !head.includes(n)).sort();
  return [...head.filter((n, i) => head.indexOf(n) === i), ...rest].map((n) => columns.get(n)).filter((c): c is Column => c !== undefined);
}

function relTables(system: TypeSystem): RelTable[] {
  const tables = new Map<string, RelTable>();
  const add = (name: string, group: string, pairs: {from: string, to: string}[], columns: Column[]) => {
    const key = name.toLowerCase();
    const held = tables.get(key);
    if (!held) {
      tables.set(key, {name, pairs, columns, groups: [group]});
      return;
    }
    if (!held.groups.includes(group)) held.groups.push(group);
    for (const pair of pairs)
      if (!held.pairs.some((p) => p.from === pair.from && p.to === pair.to)) held.pairs.push(pair);
    for (const column of columns) {
      const there = held.columns.find((c) => c.name === column.name);
      if (!there) held.columns.push(column);
      else {
        const merged = widen(there.type, column.type);
        if (!merged) throw new Error(`kuzu table ${held.name}: ${column.name} is ${there.type} and ${column.type} in ${held.groups.join(' and ')}`);
        there.type = merged;
      }
    }
  };
  const build = () => system.buildFields.edge.map((m) => ({name: m.name, type: columnType(m)}));
  for (const edge of system.edges.values()) {
    if (edge.abstract) continue;
    const pairs = cross(rootsOf(system, edge.from), rootsOf(system, edge.to));
    add(graphLabel(edge.name), edge.name, pairs, [...build(), ...Object.values(edge.properties).map((m) => ({name: m.name, type: columnType(m)}))]);
  }
  for (const type of system.nodes.values())
    for (const member of Object.values(type.own)) {
      if (member.kind !== 'ref') continue;
      const from = type.abstract ? [...system.nodes.values()].filter((n) => n.chain.includes(type.name)) : [type];
      add(member.name, member.name, cross(rootsOf(system, from.map((t) => t.name)), rootsOf(system, member.refs!)), build());
    }
  return [...tables.values()].sort((a, b) => a.name < b.name ? -1 : a.name > b.name ? 1 : 0);
}

/** The node tables a `from:`/`to:` union admits; `node` is every root. */
function rootsOf(system: TypeSystem, types: string[]): string[] {
  const roots = new Set<string>();
  for (const t of types) {
    if (t === 'node') system.roots.forEach((r) => roots.add(r));
    else {
      const root = system.nodes.get(t)?.root;
      if (root) roots.add(root);
    }
  }
  return system.roots.filter((r) => roots.has(r)).map(pascal);
}

function cross(from: string[], to: string[]): {from: string, to: string}[] {
  return from.flatMap((f) => to.map((t) => ({from: f, to: t})));
}

/** Two declarations of one column: equal types stand, a scalar pair widens to STRING, anything else is a conflict. */
function widen(a: string, b: string): string | null {
  if (a === b) return a;
  return a.endsWith('[]') || b.endsWith('[]') ? null : 'STRING';
}

export function columnType(member: Member): string {
  const base = member.kind !== 'scalar' ? 'STRING' : member.scalar === 'number' ? 'DOUBLE' : member.scalar === 'boolean' ? 'BOOLEAN' : 'STRING';
  return member.list ? `${base}[]` : base;
}

function column(c: Column): string {
  return `${quote(c.name)} ${c.type}`;
}

/** Every identifier is quoted: `order`, `key`, `from` and `type` are Cypher keywords. A name that is not a
 * schema identifier is refused here rather than escaped: quoting is not a place to accept arbitrary text. */
export function quote(identifier: string): string {
  if (!IDENTIFIER.test(identifier))
    throw new Error(`'${identifier}' is not a graph identifier (letters, digits and _, not starting with a digit); ` +
      'table and column names come from the type system (conventions.md §7.1)');
  return `\`${identifier}\``;
}

/** A single-quoted Cypher literal: a checkout path can hold an apostrophe, and a Windows one backslashes.
 * Kuzu 0.11.3 reads `\'`, not the SQL `''` — a doubled quote is a parser error there. */
export function literal(value: string): string {
  return `'${value.replace(/\\/g, '\\\\').replace(/'/g, "\\'")}'`;
}

/** `--memory <MB>` first, then `KG_KUZU_MEMORY`, then [fallback]; anything unreadable is ignored. */
export function memoryMb(option: unknown, fallback: number): number {
  for (const value of [option, process.env.KG_KUZU_MEMORY]) {
    const mb = Number(value);
    if (value !== undefined && value !== '' && Number.isFinite(mb) && mb > 0) return Math.round(mb);
  }
  return fallback;
}

export async function open(kgDir: string, readonly = true, mb?: number): Promise<{db: KuzuDatabase, conn: KuzuConnection} | null> {
  const kuzu = loadKuzu();
  if (!kuzu) return null;
  const pool = memoryMb(mb, readonly ? READ_MEMORY_MB : BUILD_MEMORY_MB) * 1024 * 1024;
  const db = new kuzu.Database(path.join(kgDir, 'kg.kuzu'), pool, true, readonly);
  return {db, conn: new kuzu.Connection(db)};
}

export interface QueryRows {
  columns: string[];
  rows: Record<string, unknown>[];
}

/** One statement, or several: the last result is the answer. Values come back as plain JSON. */
export async function run(conn: KuzuConnection, cypher: string, params?: Record<string, unknown>): Promise<QueryRows> {
  const results = many(params ? await conn.execute(await prepare(conn, cypher), params) : await conn.query(cypher));
  try {
    const last = results[results.length - 1];
    if (!last) return {columns: [], rows: []};
    const columns = await last.getColumnNames();
    const rows = (await last.getAll()).map(plain) as Record<string, unknown>[];
    return {columns, rows};
  }
  finally {
    for (const result of results) result.close();
  }
}

/** A statement whose rows nobody reads; its result is closed all the same (see the rules above). */
async function exec(conn: KuzuConnection, statement: string): Promise<void> {
  for (const result of many(await conn.query(statement))) result.close();
}

function many(result: KuzuQueryResult | KuzuQueryResult[]): KuzuQueryResult[] {
  return Array.isArray(result) ? result : [result];
}

/** A statement that will not bind reports it here, not by throwing on execute. */
async function prepare(conn: KuzuConnection, cypher: string): Promise<KuzuPreparedStatement> {
  const prepared = await conn.prepare(cypher);
  if (!prepared.isSuccess()) throw new Error(prepared.getErrorMessage());
  return prepared;
}

/** Kuzu hands back BigInt for INT64 and Date for timestamps; neither survives JSON. */
function plain(value: unknown): unknown {
  if (typeof value === 'bigint') return Number.isSafeInteger(Number(value)) ? Number(value) : value.toString();
  if (value instanceof Date) return value.toISOString();
  if (Array.isArray(value)) return value.map(plain);
  if (value && typeof value === 'object') return Object.fromEntries(Object.entries(value).filter(([k]) => k !== '_id').map(([k, v]) => [k, plain(v)]));
  return value;
}

/** Writes `<kgDir>/kg.kuzu` from the graph in `<kgDir>/data`, through CSVs under `<kgDir>/tmp`. The caller
 * builds a generation of its own, so nothing here is ever loaded over a database a reader may hold. */
export async function load(kgDir: string, system: TypeSystem, mb?: number): Promise<LoadResult> {
  const kuzu = loadKuzu();
  if (!kuzu) throw new Error(MISSING_KUZU);
  const started = Date.now();
  const schema = ddl(system);
  const dbPath = path.join(kgDir, 'kg.kuzu');
  const tmp = path.join(kgDir, 'tmp');
  for (const p of [dbPath, `${dbPath}.wal`, `${dbPath}.tmp`, tmp]) fs.rmSync(p, {recursive: true, force: true});
  fs.mkdirSync(tmp, {recursive: true});
  const db = new kuzu.Database(dbPath, memoryMb(mb, BUILD_MEMORY_MB) * 1024 * 1024);
  const conn = new kuzu.Connection(db);
  const chains = new Map<string, string[]>([...system.nodes].map(([name, type]) => [name, type.chain]));
  let peak = process.memoryUsage().rss;
  const watch = setInterval(() => peak = Math.max(peak, process.memoryUsage().rss), SAMPLE_MS);
  watch.unref();
  try {
    for (const statement of schema.statements) await exec(conn, statement);
    const table = new Map<string, string>();
    const nodes: TableRows[] = [];
    let parameterized = 0;
    for (const t of schema.nodes) {
      const csv = new CsvTable(path.join(tmp, `${t.name}.csv`), presentColumns(t, system, kgDir));
      const late: Record<string, unknown>[] = [];
      for (const type of t.types)
        for (const row of readJsonl(path.join(kgDir, 'data', 'nodes', `${type}.jsonl`))) {
          row.types = chains.get(String(row.type)) ?? [String(row.type)];
          table.set(String(row.id), t.name);
          if (!csv.write(row)) late.push(row);
        }
      csv.end();
      if (csv.rows) await exec(conn, `COPY ${quote(t.name)} (${csv.columns.map((c) => quote(c.name)).join(', ')}) FROM ${literal(posix(csv.file))} (HEADER=true, PARALLEL=false)`);
      for (const row of late) await insertNode(conn, t, row);
      parameterized += late.length;
      nodes.push({table: t.name, rows: csv.rows + late.length});
    }
    const rels: TableRows[] = [];
    for (const t of schema.rels) {
      const columns = [...EDGE_HEAD.map((name) => ({name, type: 'STRING'})), ...t.columns];
      const files = new Map<string, CsvTable>();
      const late: Record<string, unknown>[] = [];
      let rows = 0;
      for (const group of t.groups)
        for (const row of readJsonl(path.join(kgDir, 'data', 'edges', `${group}.jsonl`))) {
          const pair = {from: table.get(String(row.from)), to: table.get(String(row.to))};
          if (!pair.from || !pair.to) continue;
          const key = [pair.from, pair.to].join(SEPARATOR);
          let csv = files.get(key);
          if (!csv) files.set(key, csv = new CsvTable(path.join(tmp, `${t.name}-${pair.from}-${pair.to}.csv`), columns));
          if (csv.write(row)) rows++;
          else late.push({...row, from_table: pair.from, to_table: pair.to});
        }
      for (const [key, csv] of files) {
        csv.end();
        const [from, to] = key.split(SEPARATOR);
        if (csv.rows) await exec(conn, `COPY ${quote(t.name)} FROM ${literal(posix(csv.file))} (HEADER=true, PARALLEL=false, from=${literal(from)}, to=${literal(to)})`);
      }
      for (const row of late) await insertRel(conn, t, row);
      parameterized += late.length;
      rels.push({table: t.name, rows: rows + late.length});
    }
    await exec(conn, 'CHECKPOINT');
    peak = Math.max(peak, process.memoryUsage().rss);
    return {db: dbPath, ms: Date.now() - started, bytes: sizeOf(dbPath), nodes, rels, parameterized,
      memoryMb: Math.ceil(peak / 1048576), platform: `${process.platform}-${process.arch}`};
  }
  finally {
    clearInterval(watch);
    await conn.close();
    await db.close();
    fs.rmSync(tmp, {recursive: true, force: true});
  }
}

/** The columns of [table] the types that have a file can actually fill. */
function presentColumns(table: NodeTable, system: TypeSystem, kgDir: string): Column[] {
  const present = new Set<string>([...BUILD_ORDER, ...system.buildFields.node.map((m) => m.name)]);
  for (const type of table.types) {
    if (!fs.existsSync(path.join(kgDir, 'data', 'nodes', `${type}.jsonl`))) continue;
    for (const name of Object.keys(system.nodes.get(type)!.members)) present.add(name);
  }
  return table.columns.filter((c) => present.has(c.name));
}

export function* readJsonl(file: string): Generator<Record<string, unknown>> {
  if (!fs.existsSync(file)) return;
  for (const line of fs.readFileSync(file, 'utf8').split('\n'))
    if (line) yield JSON.parse(line);
}

/** Rows a CSV cannot carry go in one by one, so a list value is never mangled by the reader. */
async function insertNode(conn: KuzuConnection, table: NodeTable, row: Record<string, unknown>): Promise<void> {
  const columns = table.columns.filter((c) => row[c.name] !== undefined);
  const cypher = `CREATE (n:${quote(table.name)} {${columns.map((c) => `${quote(c.name)}: $${c.name}`).join(', ')}})`;
  await run(conn, cypher, Object.fromEntries(columns.map((c) => [c.name, parameter(row[c.name], c.type)])));
}

async function insertRel(conn: KuzuConnection, table: RelTable, row: Record<string, unknown>): Promise<void> {
  const columns = table.columns.filter((c) => row[c.name] !== undefined);
  const properties = columns.length ? ` {${columns.map((c) => `${quote(c.name)}: $${c.name}`).join(', ')}}` : '';
  const cypher = `MATCH (a:${quote(String(row.from_table))}), (b:${quote(String(row.to_table))}) WHERE a.${quote('id')} = $from_id AND b.${quote('id')} = $to_id ` +
    `CREATE (a)-[:${quote(table.name)}${properties}]->(b)`;
  const params: Record<string, unknown> = {from_id: row.from, to_id: row.to};
  for (const c of columns) params[c.name] = parameter(row[c.name], c.type);
  await run(conn, cypher, params);
}

/** A Record member is stored as its JSON, the way the CSV writes it. */
function parameter(value: unknown, type: string): unknown {
  return type.startsWith('STRING') && value !== null && typeof value === 'object' && !Array.isArray(value) ? JSON.stringify(value) : value;
}

/** One CSV per table (per root pair for a rel table), written as the rows stream past. */
class CsvTable {
  rows = 0;
  private buffer: string[] = [];

  constructor(readonly file: string, readonly columns: Column[]) {
    this.buffer.push(`${columns.map((c) => c.name).join(',')}\n`);
  }

  /** False when a value needs a parameter: the caller inserts that row itself. */
  write(row: Record<string, unknown>): boolean {
    const cells: string[] = [];
    for (const c of this.columns) {
      const cell = csvCell(row[c.name]);
      if (cell === null) return false;
      cells.push(cell);
    }
    this.buffer.push(`${cells.join(',')}\n`);
    this.rows++;
    if (this.buffer.length > 4096) this.flush();
    return true;
  }

  end(): void {
    this.flush();
  }

  private flush(): void {
    fs.appendFileSync(this.file, this.buffer.join(''));
    this.buffer = [];
  }
}

/** The CSV text of one value, or null when the CSV reader would not give it back unchanged. */
function csvCell(value: unknown): string | null {
  if (value === undefined || value === null) return '';
  if (typeof value === 'number' || typeof value === 'boolean') return String(value);
  if (Array.isArray(value)) {
    const items = value.map(String);
    return items.some((i) => LIST_UNSAFE.test(i)) ? null : text(`[${items.join(',')}]`);
  }
  return text(typeof value === 'object' ? JSON.stringify(value) : String(value));
}

/** RFC 4180 quoting; the reader runs with PARALLEL=false, so a newline inside a value is fine. */
function text(value: string): string {
  return `"${value.replace(/"/g, '""')}"`;
}

function posix(p: string): string {
  return p.split(path.sep).join('/');
}

function sizeOf(p: string): number {
  if (!fs.existsSync(p)) return 0;
  const stat = fs.statSync(p);
  if (!stat.isDirectory()) return stat.size;
  return fs.readdirSync(p).reduce((total, f) => total + sizeOf(path.join(p, f)), 0);
}
