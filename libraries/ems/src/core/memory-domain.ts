import {uuid4} from 'datagrok-api/u2core';
import {Filters, Access, MemoryFrame, Rows, isEmpty, Component} from '@datagrok-libraries/u2';
import type {FilterCondition, FilterGroup, IProperty, AccessData, FieldAccess,
  MemoryColumn} from '@datagrok-libraries/u2';
import {DomainBackendError, isUnscoped, toPlural, toSingular, readableInfo} from './domain-backend.js';
import type {AuditEntryLike, DomainBackend, DomainBatchOptionsLike, DomainBatchReportLike,
  DomainBatchValidationLike, DomainBatchValidationRowLike, DomainFrameLike, DomainProbeLike,
  DomainQueryLike, DomainReadScope, DomainSupportLike, DomainTableInfoLike, DomainTableLike,
  DomainTransactionOpLike, DomainTransactionResultLike} from './domain-backend.js';
import type {EditState} from './edit-state.js';
import {MemoryEditState, isJsonCell, jsonCell} from './edit-state.js';

/** A history line as kept: its `seq` is where it stands in the table's history. */
type AuditLine = Omit<AuditEntryLike, 'seq'>;

export interface MemoryColumnJson {
  type: string;
  required?: boolean;
  ref?: string;
  choices?: string[];
  min?: number;
  max?: number;
  friendlyName?: string;
  description?: string;
  editor?: string;
  isName?: boolean;
  searchable?: boolean;
  /** A ref column's grammar filter over its target (`$<column>` = this table's sibling). */
  filter?: string;
  /** A ref column's referential action, as the server's `onDelete`. */
  onDelete?: 'cascade' | 'setnull' | 'restrict';
  semType?: string;
  default?: unknown;
}

export interface MemoryTableJson {
  columns: Record<string, MemoryColumnJson>;
  businessKey?: string[];
  friendlyName?: string;
  singularName?: string;
  pluralName?: string;
  /** `{name: {check: sql}}` or `{name: {expr: grammar, message?}}` — only `expr` entries reach `info`. */
  constraints?: Record<string, {check?: string, expr?: string, message?: string}>;
  permissions?: string[] | Record<string, {description?: string}>;
  /** A self-referencing table: exactly one `ref` column must target it, and that column becomes
   * `info.parentColumn` — what `ancestors` walks. */
  hierarchy?: boolean;
}

/** The `schema.json` subset the backend reads — a real one satisfies it. */
export interface MemorySchemaJson {
  name: string;
  tables: Record<string, MemoryTableJson>;
}

export interface MemoryDomainOptions {
  /** Initial rows per table name; a row without an id gets one. */
  rows?: Record<string, Record<string, unknown>[]>;
  /** The access every table answers (default: everything allowed, every declared column
   * editable, the system columns readonly, every declared permission granted). */
  access?: AccessData;
  /** The user id stamped into `author_id` on insert (default `'me'`). */
  author?: string;
}

type Row = Record<string, unknown>;

/** The system columns every table carries (js-api `DOMAIN_SYSTEM_COLUMNS`). */
const SYSTEM: [string, string][] = [['id', 'string'], ['version', 'int'], ['created_on', 'datetime'],
  ['updated_on', 'datetime'], ['author_id', 'string']];

/** Column type → property type, as the registry's `rowProperties` answers them. */
const PROPERTY_TYPES: Record<string, string> = {
  string: 'string', int: 'int', float: 'double', bool: 'bool', datetime: 'datetime',
  string_list: 'string_list', ref: 'string', user: 'string', group: 'string', file: 'string', json: 'map',
};

const CORE_REF_SEM_TYPES: Record<string, string> = {user: 'User', group: 'Group'};

/** The store's soft-delete flag — the server's own column, projected to queries as
 * `Rows.DELETED`; not a declared property, so no form or filter sees it. */
const IS_DELETED = 'is_deleted';

const TEXT_OPERATORS = new Set(['like', '!like', 'starts', 'ends', 'matches', '!matches']);
const TEXT = '::text';

/** `DomainBackend` over a schema.json in memory, meaning what the server means: filters through the filter
 * feature's own `toMask`, a transaction ordered as the server orders it and landing on every table or none. */
export class MemoryDomainBackend implements DomainBackend {
  private readonly _tables = new Map<string, MemoryTable>();
  private readonly _schema: string;

  constructor(schema: MemorySchemaJson, options: MemoryDomainOptions = {}) {
    this._schema = schema.name;
    for (const [name, table] of Object.entries(schema.tables)) {
      this._tables.set(`${schema.name}.${name}`,
        new MemoryTable(this, schema.name, name, table, options.rows?.[name] ?? [], options.access, options.author));
    }
    for (const child of this._tables.values()) {
      for (const [column, target] of Object.entries(child.refs)) {
        this._tables.get(target)?.info.childTables.push({schema: child.schema, table: child.name, fkColumn: column,
          label: child.properties.find((p) => p.name === column)?.friendlyName ?? column,
          onDelete: schema.tables[child.name].columns[column].onDelete ?? 'restrict'});
      }
    }
  }

  table(address: string): Promise<DomainTableLike> {
    const table = this._tables.get(address);
    return table === undefined ?
      Promise.reject(new DomainBackendError('not-found', `Unknown table "${address}"`)) : Promise.resolve(table);
  }

  /** The table handle without the round trip — for tests that seed and inspect. */
  tableSync(address: string): MemoryTable | undefined {
    return this._tables.get(address);
  }

  /** Every writer's batch concatenated into one transaction, the results sliced back; every
   * writer learns every draft id the batch resolved, so a child's reference to another
   * writer's draft becomes the real id in its frame too. */
  async saveAll(edits: EditState[]): Promise<void> {
    if (edits.some((edit) => edit.isSaving.peek()))
      throw new DomainBackendError('validation', 'The batch is already being saved.');
    const writers = edits.map((edit) => {
      if (!(edit instanceof MemoryEditState))
        throw new Error('MemoryDomainBackend.saveAll: the writer is not a memory one');
      return edit;
    });
    const drafts = new Set(writers.flatMap((edit) => [...edit.draftIds()]));
    const parts = writers.map((edit) => ({edit, pending: edit.buildOps(drafts)}));
    // every participant closed for the whole transaction, as the platform session closes its
    // editors (`domains-session.ts`): an edit made meanwhile is not in the batch being sent
    for (const part of parts)
      part.edit.setSaving(true);
    try {
      const pending = parts.flatMap((p) => p.pending);
      const results = await this.transaction(pending.map((x) => x.op));
      const resolved = MemoryEditState.assignedOf(pending, results);
      let offset = 0;
      for (const part of parts) {
        part.edit.applyResults(part.pending, results.slice(offset, offset + part.pending.length), resolved);
        offset += part.pending.length;
      }
    } finally {
      for (const part of parts)
        part.edit.setSaving(false);
    }
  }

  dropAccess(): void {}

  /** The server's `/transaction`: ops may target any table (`<schema>.<table>`, or a bare name in
   * this schema); they run in a stable topological order — a `$ref` use after the insert that
   * declares it, a child table's delete before its parent's — and land all together or not at
   * all; `results[i]` answers the op at request index `i`. */
  async transaction(ops: DomainTransactionOpLike[]): Promise<DomainTransactionResultLike[]> {
    const tables = ops.map((op, i) => this._target(op.table, i));
    const declared = new Map<string, number>();
    for (const [i, op] of ops.entries()) {
      if (op.ref === undefined)
        continue;
      if (declared.has(op.ref))
        throw new DomainBackendError('bad-ref', `Operation ${i}: duplicate reference "${op.ref}"`);
      declared.set(op.ref, i);
    }
    const before = ops.map((op, i) => {
      const uses = new Set<string>();
      MemoryDomainBackend._uses(op.values, uses);
      if (op.op !== 'insert')
        MemoryDomainBackend._uses(op.id, uses);
      const deps = new Set<number>();
      for (const name of uses) {
        const at = declared.get(name);
        if (at === undefined)
          throw new DomainBackendError('bad-ref', `Operation ${i}: unknown reference "$${name}"`);
        if (at !== i)
          deps.add(at);
      }
      return deps;
    });
    const refers = (from: number, to: number): boolean => {
      const [b, a] = [tables[from], tables[to]];
      if (b !== a)
        return Object.values(b.refs).includes(a.address) && !Object.values(a.refs).includes(b.address);
      const row = a.rows.find((r) => r.id === ops[from].id);
      return row !== undefined && Object.entries(a.refs)
        .some(([column, target]) => target === a.address && row[column] === ops[to].id);
    };
    // one pass, read in both directions: a child's DELETE runs before its parent's, while a
    // restore is the mirror — the restore arm below vetoes a child coming back under a still
    // deleted parent, so the PARENT's restore runs first
    for (const op of ['delete', 'restore'] as const) {
      for (let i = 0; i < ops.length; i++) {
        for (let j = 0; j < ops.length; j++) {
          if (i === j || ops[i].op !== op || ops[j].op !== op)
            continue;
          if (op === 'delete' ? refers(i, j) : refers(j, i))
            before[j].add(i);
        }
      }
    }
    const order: number[] = [];
    const done = new Set<number>();
    while (order.length < ops.length) {
      const next = ops.findIndex((_, i) => !done.has(i) && [...before[i]].every((j) => done.has(j)));
      if (next < 0) {
        const stuck = ops.map((_, i) => i).filter((i) => !done.has(i));
        const cycle = [...declared].filter(([, i]) => !done.has(i));
        const among = cycle.length > 0 ? cycle.map(([name]) => `"${name}"`).join(', ') :
          stuck.map((i) => `${tables[i].address} "${ops[i].id}"`).join(', ');
        throw new DomainBackendError('bad-ref',
          `Operation ${cycle[0]?.[1] ?? stuck[0]}: circular reference among refs ${among}`);
      }
      order.push(next);
      done.add(next);
    }

    const copies = new Map<MemoryTable, Row[]>();
    // `rowsOf` is also how an op READS another table — a delete's FK scan, a restore's parent
    // check — so the tables it copied are not the tables the transaction wrote
    const written = new Set<MemoryTable>();
    const rowsOf = (table: MemoryTable): Row[] => {
      let rows = copies.get(table);
      if (rows === undefined)
        copies.set(table, rows = table.rows.map((row) => ({...row})));
      return rows;
    };
    const refs = new Map<string, string>();
    const resolve = (v: unknown, index: number): unknown => {
      if (Array.isArray(v))
        return v.map((x) => resolve(x, index));
      if (typeof v !== 'string' || !v.startsWith('$'))
        return v;
      if (v.startsWith('$$'))
        return v.slice(1);
      const name = v.slice(1);
      if (!refs.has(name))
        throw new DomainBackendError('bad-ref', `Operation ${index}: unknown reference "$${name}"`);
      return refs.get(name);
    };
    const results: DomainTransactionResultLike[] = new Array(ops.length);
    const audit: [MemoryTable, AuditLine][] = [];
    const tx = uuid4();
    const ts = new Date().toISOString();
    const parentsLive = (table: MemoryTable, row: Row, columns: string[]): void => {
      for (const column of columns) {
        const parentTable = this._tables.get(table.refs[column]);
        const parent = parentTable === undefined || row[column] == null ? undefined :
          rowsOf(parentTable).find((r) => r.id === row[column]);
        if (parent?.[IS_DELETED] === true) {
          throw new DomainBackendError('restrict',
            `Column "${column}" references a deleted row in "${parentTable!.name}" ("${String(row[column])}")`);
        }
      }
    };
    for (const index of order) {
      const op = ops[index];
      const table = tables[index];
      const rows = rowsOf(table);
      const values = Object.fromEntries(Object.entries(op.values ?? {})
        .map(([k, v]) => [k, table.decode(k, resolve(v, index), index)]));
      const id = op.op === 'insert' ? undefined : String(resolve(op.id, index));
      const at = rows.findIndex((row) => row.id === id);
      const entry = (name: string, rowId: string, prior: Row | null, next: Row | null): void => {
        audit.push([table, {id: rowId, tx_id: tx, op: name, actor_id: null, ts, before: prior, after: next}]);
      };
      if (op.op === 'insert') {
        const row = table.stamp(table.withDefaults(values), 1);
        table.check(row, index);
        table.checkKey(rows, row, index);
        parentsLive(table, row, Object.keys(table.refs));
        rows.push(row);
        if (op.ref !== undefined)
          refs.set(op.ref, row.id as string);
        results[index] = MemoryDomainBackend._stamped(row, 1);
        entry(op.op, row.id as string, null, {...row});
      } else if (at < 0)
        throw new DomainBackendError('not-found', `Operation ${index}: no row "${id}"`);
      else if (op.op === 'restore' && rows[at][IS_DELETED] !== true)
        throw new DomainBackendError('not-found', `Operation ${index}: no deleted row "${id}"`);
      else if (op.expectedVersion !== undefined && rows[at].version !== op.expectedVersion) {
        throw new DomainBackendError('version-conflict',
          `Operation ${index}: row "${id}" is at version ${rows[at].version}, expected ${op.expectedVersion}`);
      } else if (op.op === 'delete') {
        // the server's FK veto the child-first ordering exists to beat; a child already in the
        // trash holds nothing back
        for (const child of this._tables.values()) {
          for (const [column, target] of Object.entries(child.refs)) {
            if (target === table.address &&
                rowsOf(child).some((row) => row[column] === id && row[IS_DELETED] !== true))
              throw new DomainBackendError('validation', `Operation ${index}: row "${id}" is referenced by ${child.name}.${column}`);
          }
        }
        // soft delete: the row stays in the store and leaves every query that excludes deleted
        entry(op.op, id!, {...rows[at]}, null);
        rows[at] = table.stamp({...rows[at], [IS_DELETED]: true}, (rows[at].version as number) + 1);
        results[index] = {id};
      } else if (op.op === 'restore') {
        parentsLive(table, rows[at], Object.keys(table.refs));
        const row = table.stamp({...rows[at], [IS_DELETED]: false}, (rows[at].version as number) + 1);
        entry('undelete', id!, {...rows[at]}, {...row});
        rows[at] = row;
        results[index] = MemoryDomainBackend._stamped(row, row.version as number);
      } else {
        const moved = Object.keys(op.expected ?? {})
          .filter((c) => rows[at][c] !== resolve(op.expected![c], index));
        if (moved.length > 0) {
          throw new DomainBackendError('version-conflict',
            `Operation ${index}: row "${id}" changed since it was read: ${moved.join(', ')}`);
        }
        const row = table.stamp({...rows[at], ...values}, (rows[at].version as number) + 1);
        table.check(row, index, Object.keys(values));
        table.checkKey(rows, row, index, Object.keys(values));
        parentsLive(table, row, Object.keys(values).filter((column) => column in table.refs));
        entry(op.op, id!, {...rows[at]}, {...row});
        rows[at] = row;
        results[index] = MemoryDomainBackend._stamped(row, row.version as number);
      }
      written.add(table);
    }
    for (const [table, rows] of copies) {
      table.rows.splice(0, table.rows.length, ...rows);
      if (written.has(table))
        table.seq++;
    }
    for (const [table, line] of audit)
      table.history.push(line);
    return results;
  }

  // more than the server's {id, version}: the platform editor re-reads the row for the stamped columns
  private static _stamped(row: Row, version: number): DomainTransactionResultLike {
    return {id: row.id as string, version, created_on: row.created_on as string,
      updated_on: row.updated_on as string, author_id: row.author_id as string};
  }

  private _target(name: string, index: number): MemoryTable {
    const table = this._tables.get(name.includes('.') ? name : `${this._schema}.${name}`);
    if (table === undefined)
      throw new DomainBackendError('not-found', `Operation ${index}: unknown table "${name}"`);
    return table;
  }

  private static _uses(v: unknown, into: Set<string>): void {
    if (Array.isArray(v)) {
      for (const x of v)
        MemoryDomainBackend._uses(x, into);
    } else if (typeof v === 'object' && v !== null) {
      for (const x of Object.values(v))
        MemoryDomainBackend._uses(x, into);
    } else if (typeof v === 'string' && v.startsWith('$') && !v.startsWith('$$'))
      into.add(v.slice(1));
  }
}

export class MemoryTable implements DomainTableLike {
  /** The server's `maxUpdateWhereRows` (= `maxDeleteWhereRows`). */
  static readonly maxUpdateWhereRows = 1000;
  /** The server's `maxBatchReportRows`: a batch report itemizes at most this many rows. */
  static readonly maxBatchReportRows = 1000;
  /** How deep `ancestors` walks before it gives up — the server's recursion cap. */
  static readonly maxPathDepth = 64;

  readonly address: string;
  readonly properties: IProperty[];
  readonly info: DomainTableInfoLike;
  private readonly _columns: IProperty[];
  private readonly _info: DomainTableInfoLike;
  private readonly _hidden: string[];
  readonly support: DomainSupportLike;
  /** Installed only for a hierarchy table, as the seam's rule says — see {@link _ancestors}. */
  readonly ancestors?: DomainTableLike['ancestors'];
  /** The store — what a query copies from and a transaction writes to. */
  readonly rows: Row[];
  /** Every ref column → its target address. */
  readonly refs: Record<string, string>;
  /** What every transaction that touched this table recorded, oldest first. */
  readonly history: AuditLine[] = [];
  /** One bump per transaction that wrote this table's rows — the memory twin of
   * `domain_tables.change_seq` (consolidation 1-5), so `probe` answers an unscoped read with a
   * token instead of a scan, exactly as `DgDomainTable` does. */
  seq = 0;

  private readonly _access: AccessData;
  private readonly _author: string;

  constructor(private readonly _backend: MemoryDomainBackend, readonly schema: string, readonly name: string,
    json: MemoryTableJson, rows: Row[], access?: AccessData, author = 'me') {
    this.address = `${schema}.${name}`;
    this._author = author;
    const columns = Object.entries(json.columns);
    this._columns = [
      ...SYSTEM.map(([column, type]) => MemoryTable._property(column, type,
        {get: (row) => row[column], semType: column === 'author_id' ? 'User' : undefined})),
      ...columns.map(([column, c]) => MemoryTable._column(schema, column, c)),
    ];
    this.refs = Object.fromEntries(columns.filter(([, c]) => c.type === 'ref' && c.ref)
      .map(([column, c]) => [column, MemoryTable._address(schema, c.ref!)]));
    const named = columns.find(([, c]) => c.isName)?.[0] ??
      columns.find(([n, c]) => n === 'name' && c.type === 'string')?.[0];
    const singular = json.singularName ?? toSingular(name);
    const searchable = columns.filter(([, c]) => c.searchable === true).map(([n]) => n);
    const permissions = Array.isArray(json.permissions) ? json.permissions : Object.keys(json.permissions ?? {});
    this._info = {
      nameColumn: named ?? null, businessKey: json.businessKey ?? [], rowAddress: 'businessKey',
      singularName: singular, pluralName: json.pluralName ?? toPlural(singular),
      ...(json.friendlyName === undefined ? {} : {friendlyName: json.friendlyName}),
      searchableColumns: searchable.length > 0 ? searchable : named === undefined ? [] : [named],
      constraints: Object.entries(json.constraints ?? {}).filter(([, c]) => typeof c.expr === 'string')
        .map(([n, c]) => ({name: n, expr: c.expr!, ...(c.message === undefined ? {} : {message: c.message})})),
      refFilters: Object.fromEntries(columns.filter(([, c]) => c.filter).map(([n, c]) => [n, c.filter!])),
      permissions,
      childTables: [],
      hierarchy: json.hierarchy === true,
      parentColumn: json.hierarchy === true ? MemoryTable._parentColumn(this.address, this.refs) : null,
    };
    // `watch` is false and means it: there are no subscriptions here, and declaring that is the
    // point — a control gates on what the backend says, not on which backend it is
    this.support = {systemColumns: SYSTEM.map(([column]) => column), writes: true, deleted: true,
      restore: true, audit: true, ancestors: this._info.hierarchy === true, probe: true, version: true,
      watch: false, updateWhere: true, captions: true, transaction: true, concurrency: 'version', filters: 'full',
      // what `_plan` honours: the upsert merge needs a business key to match on
      batch: {upsert: this._info.businessKey.length > 0, partial: true, validate: true, skipDuplicates: true}};
    if (this.support.ancestors)
      this.ancestors = (id) => this._ancestors(id);
    this.rows = rows.map((row) => this.stamp({...row}, 1));
    const granted = Object.fromEntries(permissions.map((p) => [p, true]));
    this._access = access ? {can: {...granted, ...access.can}, fields: access.fields} : {
      can: {view: true, insert: true, edit: true, delete: true, share: true, ...granted},
      fields: Object.fromEntries(this._columns.map((p): [string, FieldAccess] =>
        [p.name!, SYSTEM.some(([column]) => column === p.name) ? 'readonly' : 'editable'])),
    };
    const readable = (name: string) => Access.from(this._access).field(name) !== 'hidden';
    this.properties = this._columns.filter((p) => readable(p.name!));
    this._hidden = this._columns.map((p) => p.name!)
      .filter((name) => !readable(name) && !SYSTEM.some(([column]) => column === name));
    this.info = readableInfo(this._info, this._columns, readable);
  }

  access(): Promise<AccessData> {
    return Promise.resolve({can: {...this._access.can}, fields: {...this._access.fields}});
  }

  async query(spec: DomainQueryLike = {}): Promise<Row[]> {
    this._checkColumns(spec);
    const captions = this._captionTargets(spec.captions);
    let rows = await this._where(spec);
    if (spec.sort)
      rows = MemoryTable._sorted(rows, spec.sort);
    const offset = spec.offset ?? 0;
    rows = rows.slice(offset, spec.limit === undefined ? undefined : offset + spec.limit);
    return rows.map((row) => {
      const out: Row = spec.columns ? Object.fromEntries([...SYSTEM.map(([c]) => c), ...spec.columns]
        .filter((c) => c in row).map((c) => [c, row[c]])) : {...row};
      // as the server off row mode: edit and delete are the table's answer, share is not
      // carried (null); a seed row carrying its own boolean is a row-mode table
      if (spec.withAccess) {
        for (const [column, capability] of Access.ROW_COLUMNS)
          out[column] = row[column] ?? (capability === 'share' ? null : this._access.can[capability] === true);
      }
      if (spec.deleted !== undefined && spec.deleted !== 'exclude')
        out[Rows.DELETED] = row[IS_DELETED] === true;
      for (const name of this._hidden)
        delete out[name];
      // the server's LEFT JOIN under the target's own View predicate: a null FK and a target the
      // caller may not see are the same null, and nothing says which
      for (const [name, target] of captions) {
        const hit = typeof row[name] !== 'string' ? undefined :
          target.rows.find((r) => r.id === row[name]);
        out[Rows.caption(name)] = hit === undefined ? null : target._nameOf(hit);
      }
      return out;
    });
  }

  private _checkColumns(spec: DomainQueryLike): void {
    for (const key of (spec.sort ?? '').split(',')) {
      if (this._hidden.includes(key.trim().replace(/^!/, '')))
        throw new DomainBackendError('filter', `Unknown or inaccessible sort column "${key.trim()}"`);
    }
    for (const c of spec.columns ?? []) {
      const declared = !SYSTEM.some(([s]) => s === c) && this.properties.some((p) => p.name === c) &&
        this._access.fields[c] !== 'hidden';
      if (!declared && !(spec.withAccess === true && Access.ROW_COLUMNS.some(([name]) => name === c)) &&
          !(spec.deleted !== undefined && spec.deleted !== 'exclude' && c === Rows.DELETED))
        throw new DomainBackendError('filter', `Unknown or inaccessible column "${c}"`);
    }
  }

  // the server's `_compileCaptionExpand` refusals word for word: unknown, hidden and non-ref share one (no oracle)
  private _captionTargets(names: string[] | undefined): [string, MemoryTable][] {
    if (names === undefined)
      return [];
    const seen = new Set<string>();
    return names.map((name): [string, MemoryTable] => {
      if (name.includes('.'))
        throw new DomainBackendError('filter', `Nested caption "${name}" is not supported`);
      if (seen.has(name))
        throw new DomainBackendError('filter', `Duplicate caption "${name}"`);
      seen.add(name);
      const target = this._backend.tableSync(this.refs[name] ?? '');
      if (target === undefined || this._access.fields[name] === 'hidden')
        throw new DomainBackendError('filter', `Unknown or inaccessible caption column "${name}"`);
      return [name, target];
    });
  }

  /** The rows as a `MemoryFrame` with `MemoryEditState` as its writer — the same shape the dg
   * backend answers with the js-api editor; an appended page lands in the same frame. */
  async frame(spec: DomainQueryLike): Promise<DomainFrameLike> {
    const [rows, access] = await Promise.all([this.query(spec), this.access()]);
    const columns: MemoryColumn[] = [
      ...this._columns.map((p) => ({name: p.name!,
        type: isJsonCell(p) ? 'string' : p.propertyType ?? p.type ?? 'string', semType: p.semType})),
      {name: Rows.STATE, type: 'string'},
      ...(spec.withAccess ? Access.ROW_COLUMNS.map(([name]) => ({name, type: 'bool'})) : []),
      ...(spec.deleted === undefined || spec.deleted === 'exclude' ? [] : [{name: Rows.DELETED, type: 'bool'}]),
      ...(spec.captions ?? []).map((name) => ({name: Rows.caption(name), type: 'string'})),
    ];
    const df = new MemoryFrame(columns, rows.map((row) => this._frameRow(row)));
    // a frame over deleted rows is read-only until they are restored — the same upper bound
    // `DomainSource.access` publishes, so the writer and the controls agree
    const bound = Access.from(access);
    const edit = new MemoryEditState(this, df,
      spec.deleted === undefined || spec.deleted === 'exclude' ? bound : bound.narrow({edit: false, insert: false}));
    return {
      df,
      edit,
      append: async (page) => {
        const more = await this.query(page);
        const held = new Set(df.rows.map((row) => row.id));
        df.rows.push(...more.filter((row) => row.id === undefined || !held.has(row.id))
          .map((row) => this._frameRow(row)));
        df.onRowsAdded.fire(undefined);
        return more.length;
      },
      dispose: () => edit.dispose(),
    };
  }

  private _frameRow(row: Row): Row {
    for (const prop of this.properties) {
      if (prop.name! in row)
        row[prop.name!] = jsonCell(prop, row[prop.name!]);
    }
    return row;
  }

  async count(scope: DomainReadScope = {}): Promise<number> {
    return (await this._where(scope)).length;
  }

  /** The server's live probe mirrored: the matching rows counted and the newest `updated_on`
   * among them, in one call — `aggregate({measures: [{fn: 'count'}, {fn: 'max', column:
   * 'updated_on'}]})` under the same scope. An unscoped read is answered by {@link seq} alone
   * with `count: -1`, the branch `DgDomainTable._probe` takes over the change token. */
  async probe(scope: DomainReadScope = {}): Promise<DomainProbeLike> {
    if (isUnscoped(scope))
      return {count: -1, last: String(this.seq)};
    const rows = await this._where(scope);
    let last: string | null = null;
    for (const row of rows) {
      const updated = row.updated_on;
      if (typeof updated === 'string' && (last === null || updated > last))
        last = updated;
    }
    return {count: rows.length, last};
  }

  transaction(ops: DomainTransactionOpLike[]): Promise<DomainTransactionResultLike[]> {
    return this._backend.transaction(ops);
  }

  audit(id: string, page: {limit?: number, before?: number} = {}): Promise<AuditEntryLike[]> {
    const lines = this.history.map((line, seq) => ({...line, seq}))
      .filter((line) => line.id === id && (page.before === undefined || line.seq < page.before));
    return Promise.resolve(lines.slice(-Math.min(Math.max(page.limit ?? 100, 1), 1000)));
  }

  /** The soft delete undone, as the server's `POST …/{id}/restore`: a row whose ref column points
   * at a row still in the trash is refused `restrict` naming that column — restoring it would
   * leave a live row referring to a deleted one — and a row that is not in the trash is not found.
   * The route IS the transaction op, so the veto, the stamp and the `'undelete'` audit line are
   * one rule. */
  async restore(id: string): Promise<void> {
    await this.transaction([{op: 'restore', table: this.address, id}]);
  }

  /** The filtered bulk edit, as the server's `POST …/{table}/update`: a non-empty filter selects
   * the LIVE rows the caller may edit that `values` would change (so a rerun moves on), `values`
   * is checked against the writable columns before anything is written, and at most `limit` rows
   * (the server's cap of {@link maxUpdateWhereRows} at most) are patched inside ONE transaction —
   * `hasMore` says the filter matched past the cap.
   * Any refusal rolls the whole batch back. */
  async updateWhere(filter: DomainQueryLike['filter'], values: Record<string, unknown>,
    options: {limit?: number} = {}): Promise<{updated: number, hasMore: boolean}> {
    if ((filter ?? '') === '' || (Array.isArray(filter) && filter.length === 0))
      throw new DomainBackendError('validation', `${this.address}: a filter is required`);
    if (Object.keys(values).length === 0)
      throw new DomainBackendError('validation', `${this.address}: updateWhere requires non-empty values`);
    const writable = new Set(this._columns.map((p) => p.name!)
      .filter((name) => this._access.fields[name] === 'editable'));
    for (const column of Object.keys(values)) {
      if (!writable.has(column))
        throw new DomainBackendError('validation', `${this.address}: column "${column}" cannot be updated`);
    }
    // as the server clamps it: a limit of 0 would report `hasMore` forever and write nothing
    const cap = Math.max(1, Math.min(options.limit ?? MemoryTable.maxUpdateWhereRows,
      MemoryTable.maxUpdateWhereRows));
    const plain = (v: unknown) => v instanceof Date ? v.toISOString() : v ?? null;
    // the Edit predicate and the values narrow the selection silently, as the server's do
    const matched = (await this._where({filter})).filter((row) => this._canEdit(row) &&
      Object.entries(values).some(([column, value]) => !Component.sameValue(plain(row[column]), plain(value))));
    const ids = matched.slice(0, cap).map((row) => String(row.id));
    const sent = MemoryTable._literal(values);
    await this.transaction(ids.map((id) => ({op: 'update' as const, table: this.address, id, values: sent})));
    return {updated: ids.length, hasMore: matched.length > cap};
  }

  /** The bulk upload, as the server's `POST …/{table}/batch` runs it (`batch_loader.dart`): the
   * payload's columns are checked against column security before anything else, every row is
   * validated, and then business-key duplicates — inside the batch and against the live rows —
   * are skipped and reported, or reported as errors under `errorOnDuplicate`; `'upsert'` merges
   * the matches instead of skipping them. `allOrNothing` (the default) answers the report with
   * `error` set and writes nothing, the way the platform client hands the server's envelope back;
   * otherwise the good rows land in ONE transaction and the bad ones are reported per row. */
  async batch(rows: Record<string, unknown>[],
    options: DomainBatchOptionsLike = {}): Promise<DomainBatchReportLike> {
    const {errors, duplicate, ops, posted, failed} = this._plan(rows, options);
    const cap = MemoryTable.maxBatchReportRows;
    if (options.allOrNothing !== false && errors.size > 0) {
      return {error: 'validation', inserted: 0, updated: 0, skipped: 0, errorCount: errors.size,
        rows: failed.slice(0, cap)};
    }
    const results = await this.transaction(ops.map((op) => ({...op, values: MemoryTable._literal(op.values!)})));
    const landed = posted.map((index, i) => ({index, id: results[i].id ?? null,
      status: ops[i].op === 'insert' ? 'inserted' : 'updated'}));
    return {
      inserted: ops.filter((op) => op.op === 'insert').length,
      updated: ops.filter((op) => op.op === 'update').length,
      skipped: duplicate.size,
      errorCount: errors.size,
      rows: [...failed,
        ...[...duplicate.entries()].sort(([a], [b]) => a - b).map(([index, existingId]) =>
          ({index, id: existingId, status: 'duplicate', ...(existingId === null ? {} : {existingId})})),
        ...landed.filter((r) => r.status === 'updated'), ...landed.filter((r) => r.status === 'inserted')]
        .slice(0, cap),
    };
  }

  /** The dry run, as the server's `POST …/{table}/batch?validateOnly=true`: the plan {@link batch}
   * would apply, reported and thrown away. The same code decides both, so a preview and the commit
   * cannot disagree; the verdict is about the table as it is now. */
  async validate(rows: Record<string, unknown>[],
    options: DomainBatchOptionsLike = {}): Promise<DomainBatchValidationLike> {
    const {errors, duplicate, ops, posted} = this._plan(rows, options);
    const op = new Map(posted.map((index, i) => [index, ops[i]]));
    const lines: DomainBatchValidationRowLike[] = rows.map((_, index) => {
      const problems = errors.get(index);
      if (problems !== undefined)
        return {index, predicted: 'error' as const, errors: problems};
      if (duplicate.has(index)) {
        const existingId = duplicate.get(index);
        return {index, predicted: 'skip' as const, ...(existingId === null ? {} : {existingId})};
      }
      const planned = op.get(index)!;
      return planned.op === 'insert' ? {index, predicted: 'insert' as const} :
        {index, predicted: 'update' as const, existingId: planned.id!};
    });
    return {validateOnly: true, rowCount: rows.length,
      willInsert: ops.filter((x) => x.op === 'insert').length,
      willUpdate: ops.filter((x) => x.op === 'update').length,
      willSkip: duplicate.size, errorCount: errors.size, rows: lines};
  }

  // shared by the dry run and the commit: resolves the ops, writes nothing
  private _plan(rows: Record<string, unknown>[], options: DomainBatchOptionsLike) {
    const upsert = options.mode === 'upsert';
    const key = this._info.businessKey;
    const writable = new Set(this._columns.map((p) => p.name!)
      .filter((name) => this._access.fields[name] === 'editable'));
    for (const column of new Set(rows.flatMap((row) => Object.keys(row)))) {
      if (!this._columns.some((p) => p.name === column))
        throw new DomainBackendError('validation', `Unknown column "${column}"`);
      if (!writable.has(column))
        throw new DomainBackendError('validation', `Column "${column}" is not writable`);
    }
    if (upsert) {
      if (key.length === 0)
        throw new DomainBackendError('validation', `${this.address} declares no business key to upsert by`);
      for (const column of key) {
        if (!rows.every((row) => column in row)) {
          throw new DomainBackendError('validation',
            `Upsert requires business key column "${column}" in the payload`);
        }
      }
    }
    const errors = new Map<number, {column?: string, code?: string, message: string}[]>();
    const addError = (index: number, column: string, code: string, message: string): void => {
      const list = errors.get(index) ?? [];
      list.push({column, code, message});
      errors.set(index, list);
    };
    for (const [index, row] of rows.entries()) {
      for (const prop of this._columns) {
        const problem = MemoryEditState.problemOf(prop, row[prop.name!]);
        if (problem !== null)
          addError(index, prop.name!, 'invalid-value', problem);
      }
    }
    const spell = (values: Row) => key.map((column) => String(values[column] ?? '')).join('\u0000');
    const live = new Map<string, string>();
    for (const row of this.rows) {
      if (row[IS_DELETED] !== true)
        live.set(spell(row), String(row.id));
    }
    const duplicate = new Map<number, string | null>();
    const seen = new Set<string>();
    for (const [index, row] of key.length === 0 ? [] : [...rows.entries()]) {
      // an invalid row is reported as an error and lands nowhere: it is not a duplicate, and it
      // does not take the key away from a later good row
      if (errors.has(index))
        continue;
      const spelled = spell(row);
      // the first occurrence wins: the merge may not touch one target row twice
      const clash = seen.has(spelled) ? 'Duplicate business key in batch' :
        !upsert && live.has(spelled) ? 'Duplicate business key' : null;
      seen.add(spelled);
      if (clash === null)
        continue;
      if (options.errorOnDuplicate === true)
        addError(index, key.join(','), 'unique', clash);
      else
        duplicate.set(index, live.get(spelled) ?? null);
    }
    const failed = [...errors.keys()].sort((a, b) => a - b).map((index) =>
      ({index, id: null, status: 'error', errors: errors.get(index)!}));
    const ops: DomainTransactionOpLike[] = [];
    const posted: number[] = [];
    for (const [index, values] of rows.entries()) {
      if (errors.has(index) || duplicate.has(index))
        continue;
      const id = upsert ? live.get(spell(values)) : undefined;
      posted.push(index);
      ops.push(id === undefined ? {op: 'insert', table: this.address, values} :
        {op: 'update', table: this.address, id, values});
    }
    return {errors, duplicate, ops, posted, failed};
  }

  // the seed alone may be deleted, so a row opened from a trash list still has its breadcrumb
  private async _ancestors(id: string): Promise<{id: string, name: string}[]> {
    const parent = this._info.parentColumn!;
    const visible = (key: unknown) => typeof key !== 'string' ? undefined :
      this.rows.find((row) => row.id === key && row[IS_DELETED] !== true);
    // no-oracle: a row the caller cannot see answers no path, never that it exists
    let row = this.rows.find((r) => r.id === id);
    if (row === undefined)
      return [];
    const chain: {id: string, name: string}[] = [];
    const seen = new Set<string>([id]);
    while (chain.length < MemoryTable.maxPathDepth) {
      const next = visible(row[parent]);
      if (next === undefined || seen.has(next.id as string))
        break;
      seen.add(next.id as string);
      chain.push({id: next.id as string, name: this._nameOf(next)});
      row = next;
    }
    return chain.reverse();
  }

  /** The schema's rules as the server refuses on them; an update is checked over the columns it sends alone. */
  check(row: Row, index: number, columns?: string[]): void {
    for (const prop of this._columns) {
      if (columns !== undefined && !columns.includes(prop.name!))
        continue;
      const problem = MemoryEditState.problemOf(prop, row[prop.name!]);
      if (problem !== null)
        throw new DomainBackendError('validation', `Operation ${index}: column "${prop.name}": ${problem}`);
    }
  }

  /** The business key's unique index over the live rows; an update is checked when it sends a key column. */
  checkKey(rows: Row[], row: Row, index: number, columns?: string[]): void {
    const key = this._info.businessKey;
    if (key.length > 0 && (columns === undefined || key.some((c) => columns.includes(c))) &&
        rows.some((r) => r.id !== row.id && r[IS_DELETED] !== true &&
          key.every((c) => row[c] != null && r[c] === row[c])))
      throw new DomainBackendError('validation', `Operation ${index}: duplicate business key (${key.join(', ')})`);
  }

  // a bulk edit and an upload carry values, never `$` references
  private static _literal(values: Row): Row {
    const escape = (v: unknown): unknown => Array.isArray(v) ? v.map(escape) :
      typeof v === 'string' && v.startsWith('$') ? `$${v}` : v;
    return Object.fromEntries(Object.entries(values).map(([k, v]) => [k, escape(v)]));
  }

  /** A written value as the server coerces it, for the cells the frame holds as text. */
  decode(column: string, value: unknown, index: number): unknown {
    const prop = this.properties.find((p) => p.name === column);
    if (prop === undefined || !isJsonCell(prop) || value === null || value === undefined)
      return value;
    if (value === '')
      return null;
    const type = prop.propertyType ?? prop.type;
    let decoded: unknown;
    try {
      decoded = JSON.parse(typeof value === 'string' ? value : JSON.stringify(value));
    } catch {
      decoded = undefined;
    }
    if (type === 'string_list' && Array.isArray(decoded) && !decoded.some((x) => typeof x === 'object' && x !== null))
      return decoded.map((x) => `${x}`);
    if (type === 'map' && typeof decoded === 'object' && decoded !== null && !Array.isArray(decoded))
      return decoded;
    throw new DomainBackendError('validation',
      `Operation ${index}: Cannot convert "${typeof value === 'string' ? value : JSON.stringify(value)}" to ${type}`);
  }

  withDefaults(values: Row): Row {
    for (const prop of this._columns) {
      if (isEmpty(values[prop.name!]) && prop.defaultValue != null)
        values[prop.name!] = prop.defaultValue;
    }
    return values;
  }

  stamp(row: Row, version: number): Row {
    const now = new Date().toISOString();
    row.id ??= uuid4();
    row.version = version;
    row.created_on ??= now;
    row.updated_on = now;
    row.author_id ??= this._author;
    return row;
  }

  private _canEdit(row: Row): boolean {
    const own = row[`${Access.ROW_PREFIX}edit`];
    return typeof own === 'boolean' ? own : this._access.can.edit === true;
  }

  private _nameOf(row: Row): string {
    const column = this.info.nameColumn;
    return String((column === null ? undefined : row[column]) ?? row.id);
  }

  private async _where(scope: DomainReadScope): Promise<Row[]> {
    const {filter, search, deleted = 'exclude'} = scope;
    let rows = this.rows;
    if (filter !== undefined && filter !== '') {
      const root = this._resolved(Filters.bind(MemoryTable._tree(filter), {me: this._author}));
      Filters.walk(root, (n) => {
        if (!Filters.isGroup(n) && this._hidden.includes(n.property.replace(TEXT, '')))
          throw new DomainBackendError('filter', `Unknown or inaccessible column "${n.property.replace(TEXT, '')}"`);
      });
      const frame = Filters.recordsFrame(this.rows, this._columns);
      const texts = Filters.recordsFrame(this.rows,
        this._columns.map((p) => ({...p, propertyType: 'string', type: 'string'})));
      const mask = await Filters.toMask({...frame, column: (name) => name.endsWith(TEXT) ?
        texts.column(name.slice(0, -TEXT.length)) : frame.column(name)}, root);
      rows = rows.filter((_, i) => mask.get(i));
    }
    if (deleted !== 'include')
      rows = rows.filter((row) => (row[IS_DELETED] === true) === (deleted === 'only'));
    if (search) {
      const columns = this.info.searchableColumns;
      if (columns.length === 0)
        throw new DomainBackendError('validation', `Table "${this.name}" has no searchable column`);
      const q = search.toLowerCase();
      rows = rows.filter((row) => columns.some((c) => String(row[c] ?? '').toLowerCase().includes(q)));
    }
    return rows;
  }

  // `under` becomes its ids, `fuzzy` a substring match; a text operator reads any column as text (`::text`)
  private _resolved(root: FilterGroup): FilterGroup {
    const cast = (name: string) => this._columns.some((p) => p.name === name && p.propertyType !== 'string_list');
    return {...root, nodes: root.nodes.map((node) => 'nodes' in node ? this._resolved(node) :
      node.operator === 'under' ? this._subtree(node) :
        node.operator === 'fuzzy' ? {...node, operator: 'like'} :
          TEXT_OPERATORS.has(node.operator) && cast(node.property) ? {...node, property: node.property + TEXT} : node)};
  }

  // walks the target's live rows, as the server does: a deleted branch truncates the subtree instead of leaking it
  private _subtree(condition: FilterCondition): FilterCondition {
    const property = condition.property;
    const target = property === 'id' ? this : this._backend.tableSync(this.refs[property] ?? '');
    const parent = target === undefined ? null : target.info.parentColumn;
    if (parent === null || parent === undefined)
      throw new DomainBackendError('filter', `Unknown or inaccessible column "${property}"`);
    const ids = new Set<string>([String(condition.value)]);
    for (let grew = true; grew;) {
      grew = false;
      for (const row of target!.rows) {
        if (row[IS_DELETED] !== true && typeof row.id === 'string' && !ids.has(row.id) &&
            ids.has(String(row[parent]))) {
          ids.add(row.id);
          grew = true;
        }
      }
    }
    return {...condition, operator: 'in', value: [...ids]};
  }

  private static _tree(filter: NonNullable<DomainQueryLike['filter']>): FilterGroup {
    if (typeof filter !== 'string')
      return Filters.fromDomainTree(filter);
    const {root, problems} = Filters.parse(filter);
    if (problems.length > 0)
      throw new DomainBackendError('validation', problems[0].message);
    return root;
  }

  private static _sorted(rows: Row[], sort: string): Row[] {
    const keys = sort.split(',').map((k) => k.trim()).filter((k) => k !== '')
      .map((k) => k.startsWith('!') ? {column: k.slice(1), dir: -1} : {column: k, dir: 1});
    return [...rows].sort((a, b) => {
      for (const {column, dir} of keys) {
        const x = a[column] as string | number | null;
        const y = b[column] as string | number | null;
        if (x === y)
          continue;
        // as Postgres orders: nulls last ascending, first descending
        if (x === null || x === undefined)
          return dir;
        if (y === null || y === undefined)
          return -dir;
        return (x < y ? -1 : 1) * dir;
      }
      return 0;
    });
  }

  private static _parentColumn(address: string, refs: Record<string, string>): string {
    const self = Object.entries(refs).filter(([, target]) => target === address).map(([column]) => column);
    if (self.length !== 1) {
      throw new DomainBackendError('validation', `${address}.hierarchy: invalid-hierarchy — exactly one ` +
        `ref column must target the table itself, found ${self.length}`);
    }
    return self[0];
  }

  private static _address(schema: string, ref: string): string {
    return ref.includes('.') ? ref : `${schema}.${ref}`;
  }

  private static _column(schema: string, name: string, c: MemoryColumnJson): IProperty {
    const ref = c.type === 'ref' && c.ref ? MemoryTable._address(schema, c.ref) : undefined;
    const semType = c.semType ?? CORE_REF_SEM_TYPES[c.type] ?? ref ?? (c.type === 'file' ? 'File' : undefined);
    return MemoryTable._property(name, PROPERTY_TYPES[c.type] ?? 'string', {
      get: (row) => row[name],
      set: (row, value) => row[name] = value,
      nullable: c.required !== true,
      choices: c.choices, min: c.min, max: c.max, friendlyName: c.friendlyName, description: c.description,
      editor: c.editor, semType, defaultValue: c.default,
    });
  }

  private static _property(name: string, type: string, extra: Partial<IProperty>): IProperty {
    const prop: Record<string, unknown> = {name, propertyType: type, type};
    for (const [key, value] of Object.entries(extra)) {
      if (value !== undefined)
        prop[key] = value;
    }
    return prop as IProperty;
  }
}
