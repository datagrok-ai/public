/// Entity schema snapshots, diffs and migration scaffolds for plugin manifests
/// (core/docs/features/core-tables-ems/MIGRATIONS.md §3-5). Mirrors the server's
/// DomainSnapshot (snapshot.dart): a snapshot built here from a manifest hashes
/// identically to the one the server records at deploy. Pure — no I/O.
import {createHash} from 'crypto';

export const SNAPSHOT_FORMAT = 1;

export interface SnapshotColumn {
  name: string;
  type: string;
  dbName?: string;
  ref?: string;
  onDelete?: string;
  required?: true;
  unique?: true;
  isName?: true;
  semType?: string;
  min?: number;
  max?: number;
  choices?: string[];
  default?: string;
  friendlyName?: string;
  description?: string;
  format?: string;
}

export interface SnapshotRelation {
  name: string;
  via: string;
  target: string;
  viaSelf: string;
  viaTarget: string;
  allowCreate?: false;
}

export interface SnapshotTable {
  pgTable?: string;
  securityMode?: string;
  promotion?: string;
  defaultRowVisibility?: string;
  delegate?: string;
  softDelete?: false;
  audit?: false;
  extensible?: true;
  idempotency?: true;
  ginIndex?: true;
  businessKey?: string[];
  nameColumn?: string;
  friendlyName?: string;
  description?: string;
  singularName?: string;
  pluralName?: string;
  filters?: any[];
  relations?: SnapshotRelation[];
  schemas?: string[];
  columns: SnapshotColumn[];
}

export interface Snapshot {
  format: number;
  schema: string;
  version?: string;
  hash?: string;
  extensibleTables?: true;
  propertySchemas?: {[name: string]: SnapshotColumn[]};
  tables: {[name: string]: SnapshotTable};
}

export interface Change {
  kind: string;
  table?: string;
  column?: string;
  detail: string;
  physical: boolean;
  auto: boolean;
  before?: any;
  after?: any;
}

const sqlTypes: {[type: string]: string} = {
  string: 'text', int: 'int', float: 'float8', bool: 'bool', datetime: 'timestamp without time zone',
  string_list: 'text[]', ref: 'uuid', user: 'uuid', group: 'uuid', file: 'text',
};
const coreRefTables: {[type: string]: string} = {user: 'users', group: 'groups'};

/** Every table knob the diff compares, in report order (snapshot.dart `_knobs`). */
const knobs = ['pgTable', 'securityMode', 'promotion', 'defaultRowVisibility', 'delegate', 'softDelete',
  'audit', 'extensible', 'idempotency', 'ginIndex', 'readOnly', 'entityBacked', 'systemColumns',
  'businessKey', 'nameColumn', 'friendlyName', 'description', 'singularName', 'pluralName', 'schemas'];
/** Registry-only knobs the plugin deploy engine refuses to change. */
const refusedKnobs = ['securityMode', 'promotion', 'defaultRowVisibility', 'delegate', 'softDelete'];
/** Knobs that add or drop an index (and a column, for idempotency). */
const indexKnobs = ['businessKey', 'idempotency', 'ginIndex'];
const columnMeta = ['semType', 'min', 'max', 'choices', 'default', 'friendlyName', 'description', 'format', 'isName'];

/** Compact JSON with recursively sorted keys and integral numbers printed as integers —
 * the one form both implementations hash. Keys sort by UTF-16 code units (default `<`). */
export function canonical(v: any): string {
  if (v === null || v === undefined)
    return 'null';
  if (typeof v === 'boolean' || typeof v === 'number')
    return String(v);
  if (typeof v === 'string')
    return JSON.stringify(v);
  if (Array.isArray(v))
    return `[${v.map(canonical).join(',')}]`;
  if (typeof v === 'object') {
    const keys = Object.keys(v).filter((k) => v[k] !== undefined).sort();
    return `{${keys.map((k) => `${JSON.stringify(k)}:${canonical(v[k])}`).join(',')}}`;
  }
  return JSON.stringify(String(v));
}

/** Everything the hash covers — and nothing a rename or a version bump touches. */
export function snapshotBody(s: Snapshot): any {
  const body: any = {};
  if (s.extensibleTables === true)
    body.extensibleTables = true;
  if (s.propertySchemas != null && Object.keys(s.propertySchemas).length > 0)
    body.propertySchemas = s.propertySchemas;
  body.tables = s.tables;
  return body;
}

export function hashOf(s: Snapshot): string {
  return createHash('sha256').update(canonical(snapshotBody(s)), 'utf8').digest('hex');
}

/** The canonical column map: every field at its engine default is omitted. */
function columnOf(name: string, c: any): SnapshotColumn {
  const m: SnapshotColumn = {name: name, type: c.type};
  if (c.dbName != null && c.dbName !== name)
    m.dbName = c.dbName;
  if (c.ref != null)
    m.ref = c.ref;
  if (c.onDelete != null)
    m.onDelete = c.onDelete;
  if (c.required === true)
    m.required = true;
  if (c.unique === true)
    m.unique = true;
  if (c.isName === true)
    m.isName = true;
  if (c.semType != null)
    m.semType = c.semType;
  if (typeof c.min === 'number')
    m.min = c.min;
  if (typeof c.max === 'number')
    m.max = c.max;
  if (Array.isArray(c.choices) && c.choices.length > 0)
    m.choices = c.choices.map(String);
  if (c.default != null)
    m.default = String(c.default);
  if (c.friendlyName != null)
    m.friendlyName = c.friendlyName;
  if (c.description != null)
    m.description = c.description;
  if (c.format != null)
    m.format = c.format;
  return m;
}

function columnsOf(columns: any): SnapshotColumn[] {
  return Object.keys(columns ?? {}).map((name) => columnOf(name, columns[name]));
}

/** The server's display-name normalization: trimmed, empty collapses to absent. */
function displayName(v: any): string | undefined {
  const s = typeof v === 'string' ? v.trim() : '';
  return s.length > 0 ? s : undefined;
}

/** The server's primary-name resolution: the one `isName` column (relational or from an
 * attached property schema), else a string column literally named `name`. */
function nameColumnOf(table: any, manifest: any): string | undefined {
  const groups: any[] = [table.columns ?? {}];
  for (const s of table.schemas ?? [])
    groups.push(manifest.propertySchemas?.[s] ?? {});
  for (const g of groups) {
    for (const name of Object.keys(g)) {
      if (g[name].isName === true)
        return name;
    }
  }
  for (const g of groups) {
    if (g.name?.type === 'string')
      return 'name';
  }
  return undefined;
}

/** One side of a junction: the declared column, otherwise the single ref column of
 * [via] pointing at [to]. */
function relationSide(manifest: any, tableName: string, relName: string, via: string,
  explicit: string | undefined, to: string): string {
  if (explicit != null)
    return explicit;
  const columns = manifest.tables[via]?.columns ?? {};
  const candidates = Object.keys(columns).filter((c) => columns[c].type === 'ref' && columns[c].ref === to);
  if (candidates.length !== 1) {
    throw new Error(`relation '${tableName}.${relName}': junction table '${via}' has ` +
      `${candidates.length === 0 ? 'no' : 'more than one'} ref column targeting '${to}'`);
  }
  return candidates[0];
}

function relationsOf(manifest: any, tableName: string, table: any): SnapshotRelation[] {
  const declared = table.relations ?? {};
  return Object.keys(declared).map((name) => {
    const r = declared[name];
    const m: SnapshotRelation = {
      name: name, via: r.via, target: r.target,
      viaSelf: relationSide(manifest, tableName, name, r.via, r.viaSelf, tableName),
      viaTarget: relationSide(manifest, tableName, name, r.via, r.viaTarget, r.target),
    };
    if (r.allowCreate === false)
      m.allowCreate = false;
    return m;
  });
}

/** The server's filter-descriptor normalization: `{column, type?, bins?, label?}`. */
function filtersOf(filters: any): any[] | undefined {
  if (!Array.isArray(filters) || filters.length === 0)
    return undefined;
  return filters.map((f) => {
    const m: any = {column: f.column};
    if (f.type != null)
      m.type = String(f.type);
    if (Number.isInteger(f.bins))
      m.bins = f.bins;
    const label = displayName(f.label);
    if (label != null)
      m.label = label;
    return m;
  });
}

/** The canonical table map: every knob at its engine default is omitted, so two
 * declarations that mean the same thing serialize identically. */
function tableOf(manifest: any, name: string, t: any): SnapshotTable {
  // Key order mirrors snapshot.dart's tableOf (columns last), so sealed files match too.
  const m: Partial<SnapshotTable> = {};
  const securityMode = t.securityMode ?? 'table';
  if (securityMode !== 'table')
    m.securityMode = securityMode;
  if (securityMode === 'row' && (t.promotion ?? 'lazy') !== 'lazy')
    m.promotion = t.promotion;
  if ((t.defaultRowVisibility ?? 'table') !== 'table')
    m.defaultRowVisibility = t.defaultRowVisibility;
  if (t.delegate != null)
    m.delegate = t.delegate;
  if (t.softDelete === false)
    m.softDelete = false;
  if (t.audit === false)
    m.audit = false;
  if (t.extensible === true)
    m.extensible = true;
  if (t.idempotency === true)
    m.idempotency = true;
  if (t.ginIndex === true)
    m.ginIndex = true;
  if (Array.isArray(t.businessKey) && t.businessKey.length > 0)
    m.businessKey = t.businessKey.map(String);
  const nameColumn = nameColumnOf(t, manifest);
  if (nameColumn != null)
    m.nameColumn = nameColumn;
  if (t.friendlyName != null)
    m.friendlyName = t.friendlyName;
  if (t.description != null)
    m.description = t.description;
  const singular = displayName(t.singularName);
  if (singular != null)
    m.singularName = singular;
  const plural = displayName(t.pluralName);
  if (plural != null)
    m.pluralName = plural;
  const filters = filtersOf(t.filters);
  if (filters != null)
    m.filters = filters;
  const relations = relationsOf(manifest, name, t);
  if (relations.length > 0)
    m.relations = relations;
  if (Array.isArray(t.schemas) && t.schemas.length > 0)
    m.schemas = t.schemas.map(String).sort();
  m.columns = columnsOf(t.columns);
  return m as SnapshotTable;
}

/** Builds the snapshot of a parsed (and validated) `schema.json` manifest, applying the
 * same normalization the server's manifest parser does before it records one. */
export function buildSnapshot(manifest: any): Snapshot {
  const s: Snapshot = {format: SNAPSHOT_FORMAT, schema: manifest.name, version: manifest.version, tables: {}};
  if (manifest.extensible?.tables === true)
    s.extensibleTables = true;
  const schemas = manifest.propertySchemas ?? {};
  const schemaNames = Object.keys(schemas).sort();
  if (schemaNames.length > 0) {
    s.propertySchemas = {};
    for (const n of schemaNames)
      s.propertySchemas[n] = columnsOf(schemas[n]);
  }
  for (const n of Object.keys(manifest.tables ?? {}).sort())
    s.tables[n] = tableOf(manifest, n, manifest.tables[n]);
  s.hash = hashOf(s);
  return s;
}

/** The sealed-file form: stable key order, sorted tables, two-space indent, trailing newline. */
export function serialize(s: Snapshot): string {
  const out: any = {format: s.format, schema: s.schema, version: s.version, hash: s.hash ?? hashOf(s)};
  const body = snapshotBody(s);
  if (body.extensibleTables)
    out.extensibleTables = true;
  if (body.propertySchemas) {
    out.propertySchemas = {};
    for (const n of Object.keys(body.propertySchemas).sort())
      out.propertySchemas[n] = body.propertySchemas[n];
  }
  out.tables = {};
  for (const n of Object.keys(s.tables).sort())
    out.tables[n] = s.tables[n];
  return JSON.stringify(out, null, 2) + '\n';
}

/** A snapshot with nothing declared — the baseline when no sealed or recorded one exists. */
export function emptySnapshot(schema: string): Snapshot {
  return {format: SNAPSHOT_FORMAT, schema: schema, tables: {}};
}

function show(v: any): string {
  return v == null ? '(none)' : typeof v === 'string' ? `"${v}"` : canonical(v);
}

export function describeChange(c: Change): string {
  const path = c.table == null ? '' : c.column == null ? c.table : `${c.table}.${c.column}`;
  return `${c.kind} ${path === '' ? '' : `${path}: `}${c.detail}` +
    `${c.physical ? ' [physical]' : ''}${c.auto ? '' : ' [manual]'}`;
}

function change(kind: string, detail: string,
  o: {table?: string, column?: string, physical?: boolean, auto?: boolean, before?: any, after?: any}): Change {
  return {kind, detail, table: o.table, column: o.column, physical: o.physical ?? false, auto: o.auto ?? true,
    before: o.before, after: o.after};
}

function diffColumn(changes: Change[], t: string, n: string, a: SnapshotColumn, b: SnapshotColumn): void {
  if (a.type !== b.type) {
    changes.push(change('change-type', `type ${a.type} -> ${b.type}`,
      {table: t, column: n, physical: true, auto: false, before: a, after: b}));
  } else if (a.ref !== b.ref) {
    changes.push(change('change-ref', `ref ${show(a.ref)} -> ${show(b.ref)}`,
      {table: t, column: n, physical: true, auto: false, before: a, after: b}));
  }
  if ((a.required === true) !== (b.required === true)) {
    changes.push(change('change-required', `required ${a.required === true} -> ${b.required === true}`,
      {table: t, column: n, physical: true, auto: false, before: a, after: b}));
  }
  if ((a.unique === true) !== (b.unique === true)) {
    changes.push(change('change-unique', `unique ${a.unique === true} -> ${b.unique === true}`,
      {table: t, column: n, physical: true, auto: b.unique === true, before: a, after: b}));
  }
  if (a.onDelete !== b.onDelete) {
    changes.push(change('change-on-delete', `onDelete ${show(a.onDelete)} -> ${show(b.onDelete)}`,
      {table: t, column: n, auto: false, before: a, after: b}));
  }
  const meta: string[] = [];
  for (const k of columnMeta) {
    if (canonical((a as any)[k]) !== canonical((b as any)[k]))
      meta.push(`${k} ${show((a as any)[k])} -> ${show((b as any)[k])}`);
  }
  if (meta.length > 0)
    changes.push(change('column-meta', meta.join(', '), {table: t, column: n, before: a, after: b}));
}

function diffTable(changes: Change[], t: string, a: SnapshotTable, b: SnapshotTable): void {
  for (const k of knobs) {
    if (canonical((a as any)[k]) === canonical((b as any)[k]))
      continue;
    const index = indexKnobs.includes(k);
    changes.push(change('table-knob', `${k} ${show((a as any)[k])} -> ${show((b as any)[k])}`,
      {table: t, physical: index, auto: !index && !refusedKnobs.includes(k), before: a, after: b}));
  }
  if (canonical(a.filters) !== canonical(b.filters))
    changes.push(change('filters', 'declared filters changed', {table: t}));
  if (canonical(a.relations) !== canonical(b.relations))
    changes.push(change('relations', 'declared relations changed', {table: t}));

  const before: {[name: string]: SnapshotColumn} = {};
  for (const c of a.columns ?? [])
    before[c.name] = c;
  const after: {[name: string]: SnapshotColumn} = {};
  for (const c of b.columns ?? [])
    after[c.name] = c;

  const dropped: Change[] = [];
  const added: Change[] = [];
  for (const n of Object.keys(before)) {
    if (after[n] != null)
      continue;
    const ch = change('drop-column', 'dropped', {table: t, column: n, physical: true, auto: false, before: before[n]});
    dropped.push(ch);
    changes.push(ch);
  }
  for (const n of Object.keys(after)) {
    const nc = after[n];
    if (before[n] == null) {
      const required = nc.required === true && nc.default == null;
      const ch = change('add-column', `added (${nc.type}${required ? ', required without a default' : ''})`,
        {table: t, column: n, physical: true, auto: !required, after: nc});
      added.push(ch);
      changes.push(ch);
    } else
      diffColumn(changes, t, n, before[n], nc);
  }
  if (dropped.length === 1 && added.length === 1 && dropped[0].before.type === added[0].after.type) {
    const hint = ` (rename? ${dropped[0].column} -> ${added[0].column})`;
    dropped[0].detail += hint;
    added[0].detail += hint;
  }

  const kept = Object.keys(before).filter((n) => after[n] != null);
  const keptAfter = Object.keys(after).filter((n) => before[n] != null);
  if (kept.join(',') !== keptAfter.join(','))
    changes.push(change('reorder-columns', 'column order changed', {table: t}));
}

/** The ordered change list between two snapshots, classified with the plugin rules
 * (MIGRATIONS.md §4): retiring a table or column is a physical drop the deployer refuses. */
export function diff(from: Snapshot, to: Snapshot): Change[] {
  const changes: Change[] = [];
  const fromExt = from.extensibleTables === true, toExt = to.extensibleTables === true;
  if (fromExt !== toExt)
    changes.push(change('schema-extensible', `extensible tables ${fromExt} -> ${toExt}`, {}));

  const fromTables = Object.keys(from.tables).sort();
  const toTables = Object.keys(to.tables).sort();
  for (const t of fromTables) {
    if (to.tables[t] == null)
      changes.push(change('drop-table', 'dropped', {table: t, physical: true, auto: false, before: from.tables[t]}));
  }
  for (const t of toTables) {
    if (from.tables[t] == null)
      changes.push(change('add-table', 'created', {table: t, physical: true, after: to.tables[t]}));
    else
      diffTable(changes, t, from.tables[t], to.tables[t]);
  }

  const fromSchemas = from.propertySchemas ?? {}, toSchemas = to.propertySchemas ?? {};
  const schemas = new Set([...Object.keys(fromSchemas), ...Object.keys(toSchemas)]);
  for (const s of [...schemas].sort()) {
    const a = fromSchemas[s], b = toSchemas[s];
    if (a == null || b == null || canonical(a) !== canonical(b))
      changes.push(change('property-schema', a == null ? 'added' : b == null ? 'removed' : 'changed', {table: s}));
  }
  return changes;
}

// --- SQL scaffold (DomainDdlGenerator vocabulary) ---------------------------------------

function physicalName(c: SnapshotColumn): string {
  return c.dbName ?? c.name;
}

function columnDefinition(c: SnapshotColumn): string {
  return `${physicalName(c)} ${sqlTypes[c.type]}${c.required === true ? ' NOT NULL' : ''}`;
}

function sqlLiteral(c: SnapshotColumn): string {
  const v = c.default;
  if (v == null)
    return 'NULL';
  if (c.type === 'int' || c.type === 'float')
    return v;
  if (c.type === 'bool')
    return v.toLowerCase() === 'true' ? 'true' : 'false';
  return `'${v.replace(/'/g, "''")}'`;
}

/** The physical table a reference column's FK points at; a cross-schema `<Schema>.<table>`
 * ref resolves through the registry at deploy time, so it gets no FK here. */
function refTarget(c: SnapshotColumn, schema: string): string | null {
  if (coreRefTables[c.type] != null)
    return `public.${coreRefTables[c.type]}`;
  if (c.type !== 'ref' || c.ref == null || c.ref.includes('.'))
    return null;
  return `${schema}.${c.ref}`;
}

function addConstraint(qt: string, name: string, definition: string): string {
  return ['DO $$', 'BEGIN',
    `  IF NOT EXISTS (SELECT 1 FROM pg_constraint WHERE conname = '${name}') THEN`,
    `    ALTER TABLE ${qt} ADD CONSTRAINT ${name} ${definition};`,
    '  END IF;', 'END$$'].join('\n');
}

function uniqueIndex(schema: string, table: string, c: SnapshotColumn): string {
  return `CREATE UNIQUE INDEX IF NOT EXISTS ux_${table}_${physicalName(c)} ` +
    `ON ${schema}.${table} (${physicalName(c)}) WHERE NOT is_deleted`;
}

function businessKeyIndex(schema: string, table: string, key: string[] | undefined): string[] {
  return key == null || key.length === 0 ? [] :
    [`CREATE UNIQUE INDEX IF NOT EXISTS ux_${table}_business_key ON ${schema}.${table} (${key.join(', ')}) ` +
      'WHERE NOT is_deleted'];
}

/** Indexes and FKs of [columns] (plus the table-level ones when [t] is given). */
function tableConstraints(schema: string, table: string, columns: SnapshotColumn[], t?: SnapshotTable): string[] {
  const qt = `${schema}.${table}`;
  const res: string[] = [];
  for (const c of columns) {
    if (c.unique === true)
      res.push(uniqueIndex(schema, table, c));
  }
  if (t != null) {
    res.push(...businessKeyIndex(schema, table, t.businessKey));
    if (t.idempotency === true) {
      res.push(`CREATE UNIQUE INDEX IF NOT EXISTS ux_${table}_idempotency_key ON ${qt} (idempotency_key) ` +
        'WHERE idempotency_key IS NOT NULL');
    }
  }
  for (const c of columns) {
    const target = refTarget(c, schema);
    if (target != null) {
      res.push(addConstraint(qt, `fk_${table}_${physicalName(c)}`,
        `FOREIGN KEY (${physicalName(c)}) REFERENCES ${target}(id)`));
    }
    if (c.type === 'ref' || coreRefTables[c.type] != null)
      res.push(`CREATE INDEX IF NOT EXISTS ix_${table}_${physicalName(c)} ON ${qt} (${physicalName(c)})`);
  }
  if (t?.ginIndex === true)
    res.push(`CREATE INDEX IF NOT EXISTS ix_${table}_data ON ${qt} USING gin (data jsonb_path_ops)`);
  return res;
}

function createTable(schema: string, table: string, t: SnapshotTable): string[] {
  const cols = [
    'id uuid PRIMARY KEY DEFAULT uuid_generate_v4()',
    'is_deleted bool NOT NULL DEFAULT false',
    'version int NOT NULL DEFAULT 1',
    'created_on timestamp without time zone NOT NULL DEFAULT now()',
    'updated_on timestamp without time zone NOT NULL DEFAULT now()',
    'author_id uuid',
  ];
  if (t.idempotency === true)
    cols.push('idempotency_key uuid');
  cols.push("data jsonb NOT NULL DEFAULT '{}'");
  for (const c of t.columns ?? [])
    cols.push(columnDefinition(c));
  return [`CREATE TABLE IF NOT EXISTS ${schema}.${table} (\n  ${cols.join(',\n  ')}\n)`,
    ...tableConstraints(schema, table, t.columns ?? [], t)];
}

function addColumn(schema: string, table: string, c: SnapshotColumn): string[] {
  let def = columnDefinition(c);
  if (c.required === true && c.default != null)
    def += ` DEFAULT ${sqlLiteral(c)}`;
  return [`ALTER TABLE ${schema}.${table} ADD COLUMN IF NOT EXISTS ${def}`,
    ...tableConstraints(schema, table, [c])];
}

function swapFk(schema: string, table: string, a: SnapshotColumn, b: SnapshotColumn): string[] {
  const res = [`ALTER TABLE ${schema}.${table} DROP CONSTRAINT IF EXISTS fk_${table}_${physicalName(a)}`,
    ...tableConstraints(schema, table, [b]).filter((s) => s.includes('FOREIGN KEY'))];
  if (b.type === 'ref' && b.ref != null && b.ref.includes('.'))
    res.push(`-- TODO: the cross-schema target of "${b.ref}" resolves through the registry`);
  return res;
}

function knobStatements(schema: string, c: Change): {up: string[], down: string[]} {
  const t = c.table!;
  const qt = `${schema}.${t}`;
  const before: SnapshotTable = c.before, after: SnapshotTable = c.after;
  if (c.detail.startsWith('businessKey')) {
    return {
      up: [`DROP INDEX IF EXISTS ${schema}.ux_${t}_business_key`, ...businessKeyIndex(schema, t, after.businessKey)],
      down: [`DROP INDEX IF EXISTS ${schema}.ux_${t}_business_key`, ...businessKeyIndex(schema, t, before.businessKey)],
    };
  }
  if (c.detail.startsWith('idempotency')) {
    const add = [`ALTER TABLE ${qt} ADD COLUMN IF NOT EXISTS idempotency_key uuid`,
      `CREATE UNIQUE INDEX IF NOT EXISTS ux_${t}_idempotency_key ON ${qt} (idempotency_key) ` +
        'WHERE idempotency_key IS NOT NULL'];
    const drop = [`DROP INDEX IF EXISTS ${schema}.ux_${t}_idempotency_key`,
      `ALTER TABLE ${qt} DROP COLUMN IF EXISTS idempotency_key`];
    return after.idempotency === true ? {up: add, down: drop} : {up: drop, down: add};
  }
  if (c.detail.startsWith('ginIndex')) {
    const add = [`CREATE INDEX IF NOT EXISTS ix_${t}_data ON ${qt} USING gin (data jsonb_path_ops)`];
    const drop = [`DROP INDEX IF EXISTS ${schema}.ix_${t}_data`];
    return after.ginIndex === true ? {up: add, down: drop} : {up: drop, down: add};
  }
  return {up: [], down: []};
}

/** The up and down statements of one physical change; empty = nothing could be inferred. */
function render(schema: string, c: Change): {up: string[], down: string[]} {
  const t = c.table!;
  const qt = `${schema}.${t}`;
  switch (c.kind) {
    case 'add-table':
      return {up: createTable(schema, t, c.after), down: ['-- data loss', `DROP TABLE IF EXISTS ${qt}`]};
    case 'drop-table':
      return {up: ['-- data loss', `DROP TABLE IF EXISTS ${qt}`],
        down: ['-- data is not restored', ...createTable(schema, t, c.before)]};
    case 'add-column': {
      const mc: SnapshotColumn = c.after;
      const up = addColumn(schema, t, mc);
      if (mc.required === true && mc.default == null)
        up.unshift('-- TODO: a required column with no default cannot be added to a populated table');
      return {up, down: [`ALTER TABLE ${qt} DROP COLUMN IF EXISTS ${physicalName(mc)}`]};
    }
    case 'drop-column': {
      const mc: SnapshotColumn = c.before;
      return {up: ['-- data loss', `ALTER TABLE ${qt} DROP COLUMN IF EXISTS ${physicalName(mc)}`],
        down: ['-- data is not restored', ...addColumn(schema, t, mc)]};
    }
    case 'change-type': {
      const a: SnapshotColumn = c.before, b: SnapshotColumn = c.after;
      const ta = sqlTypes[a.type], tb = sqlTypes[b.type];
      const review = '-- review: values that do not survive the cast fail the statement';
      return {
        up: [review, `ALTER TABLE ${qt} ALTER COLUMN ${physicalName(b)} TYPE ${tb} USING ${physicalName(b)}::${tb}`],
        down: [review, `ALTER TABLE ${qt} ALTER COLUMN ${physicalName(a)} TYPE ${ta} USING ${physicalName(a)}::${ta}`],
      };
    }
    case 'change-ref':
      return {up: swapFk(schema, t, c.before, c.after), down: swapFk(schema, t, c.after, c.before)};
    case 'change-required': {
      const b: SnapshotColumn = c.after;
      const col = physicalName(b);
      const setNotNull = [b.default != null
        ? `UPDATE ${qt} SET ${col} = ${sqlLiteral(b)} WHERE ${col} IS NULL`
        : `-- TODO: backfill NULLs in ${qt}.${col} before SET NOT NULL`,
      `ALTER TABLE ${qt} ALTER COLUMN ${col} SET NOT NULL`];
      const dropNotNull = [`ALTER TABLE ${qt} ALTER COLUMN ${col} DROP NOT NULL`];
      return b.required === true ? {up: setNotNull, down: dropNotNull} : {up: dropNotNull, down: setNotNull};
    }
    case 'change-unique': {
      const b: SnapshotColumn = c.after;
      const create = uniqueIndex(schema, t, b);
      const drop = `DROP INDEX IF EXISTS ${schema}.ux_${t}_${physicalName(b)}`;
      return b.unique === true ? {up: [create], down: [drop]} : {up: [drop], down: [create]};
    }
    case 'table-knob':
      return knobStatements(schema, c);
  }
  return {up: [], down: []};
}

export interface ScaffoldOptions {
  from?: string;
  to?: string;
  title?: string;
  generatedBy?: string;
}

function shortHash(h: string | undefined): string {
  return h == null ? '(none)' : h.substring(0, 12);
}

function script(schema: string, changes: Change[], rendered: Change[], isDown: boolean, o: ScaffoldOptions): string {
  const lines: string[] = [];
  const from = isDown ? o.to : o.from, to = isDown ? o.from : o.to;
  lines.push(`-- ${o.title ?? 'Schema migration'}${isDown ? ' (down)' : ''}`);
  // The deployer runs a script only when this transition matches the recorded snapshot.
  lines.push(`-- ems-migration: from=${from ?? 'none'} to=${to ?? 'none'}`);
  lines.push(`-- Generated by ${o.generatedBy ?? 'grok schema migrate'} from snapshot ` +
    `${shortHash(from)} to ${shortHash(to)}.`);
  lines.push('-- Review every statement before committing.');
  const info = changes.filter((c) => !c.physical);
  if (info.length > 0) {
    lines.push('-- Registry-only changes (no SQL; converged by the next deploy):');
    for (const c of info)
      lines.push(`--   * ${describeChange(c)}`);
  }
  const applied = changes.filter((c) => c.physical && c.auto);
  if (!isDown && applied.length > 0) {
    lines.push('-- Applied by the deployer itself (no SQL here; reversed in the down script):');
    for (const c of applied)
      lines.push(`--   * ${describeChange(c)}`);
  }
  if (rendered.length === 0)
    lines.push(isDown ? '-- No physical change.' : '-- No manual statement.');
  for (const c of rendered) {
    lines.push('');
    lines.push(`-- ${describeChange(c)}`);
    const statements = render(schema, c)[isDown ? 'down' : 'up'];
    if (statements.length === 0)
      lines.push('-- TODO: no statement could be inferred for this change.');
    for (const s of statements)
      lines.push(s.endsWith(';') || s.startsWith('--') ? s : `${s};`);
  }
  return lines.join('\n') + '\n';
}

/** Renders a diff as an up and a down script (MIGRATIONS.md §5). Only the physical changes
 * the deployer refuses get an up statement — it applies the additive ones itself, and a
 * script duplicating them would race the deploy; the down script reverses every physical
 * change, in reverse order. Registry-only changes are listed in the header. */
export function scaffold(changes: Change[], schema: string, options: ScaffoldOptions = {}): {up: string, down: string} {
  const physical = changes.filter((c) => c.physical);
  return {
    up: script(schema, changes, physical.filter((c) => !c.auto), false, options),
    down: script(schema, changes, [...physical].reverse(), true, options),
  };
}
