/* `DomainBackend` over a schema.json held in memory (GOAL step-back "A local in-memory domain
   backend"): the same seam the platform fills, answered without a server — the gallery, the
   headless tests, an agent iterating on an app. Filters are evaluated with the filter feature's
   own `toMask`, so a query means here what it means on the server; a transaction is ordered the
   way the server orders it (forward `$ref`s, child-first deletes) and lands on every table it
   touches or on none. */
import {BitArray} from 'datagrok-api/u2core';
import {Filters} from '../core/filter/index.js';
import {uuid4} from '../core/uuid.js';
import {notify} from '../components/display/notify.js';
import type {FilterGroup} from '../core/filter/model.js';
import {INT_NULL, FLOAT_NULL} from '../core/filter/evaluate.js';
import type {MaskColumnLike, MaskFrameLike} from '../core/filter/evaluate.js';
import type {IProperty} from '../core/property-like.js';
import type {AccessData, FieldAccess} from '../core/access.js';
import {Access} from '../core/access.js';
import {DomainBackendError} from './domain-backend.js';
import type {AuditEntryLike, DomainBackend, DomainDeletedMode, DomainFrameLike, DomainQueryLike,
  DomainTableInfoLike, DomainTableLike, DomainTransactionOpLike,
  DomainTransactionResultLike} from './domain-backend.js';
import type {EditState} from './edit-state.js';
import {MemoryEditState} from './edit-state.js';
import {MemoryFrame} from './memory-frame.js';
import type {MemoryColumn} from './memory-frame.js';
import {Rows} from './rows-like.js';

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
          label: child.properties.find((p) => p.name === column)?.friendlyName ?? column});
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
  async saveAll(edits: EditState[]): Promise<boolean> {
    if (edits.some((edit) => edit.isSaving.peek())) {
      notify.warning('The batch is already being saved.');
      return false;
    }
    const parts = edits.map((edit) => ({edit: edit as MemoryEditState, pending: (edit as MemoryEditState).buildOps()}));
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
      return true;
    } finally {
      for (const part of parts)
        part.edit.setSaving(false);
    }
  }

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
    for (const [i, child] of ops.entries()) {
      for (const [j, parent] of ops.entries()) {
        if (i === j || child.op !== 'delete' || parent.op !== 'delete')
          continue;
        const [b, a] = [tables[i], tables[j]];
        if (b === a) {
          const row = a.rows.find((r) => r.id === child.id);
          if (row !== undefined && Object.entries(a.refs)
            .some(([column, target]) => target === a.address && row[column] === parent.id))
            before[j].add(i);
        }
        else if (Object.values(b.refs).includes(a.address) && !Object.values(a.refs).includes(b.address))
          before[j].add(i);
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
    const audit: [MemoryTable, AuditEntryLike][] = [];
    const tx = uuid4();
    const ts = new Date().toISOString();
    for (const index of order) {
      const op = ops[index];
      const table = tables[index];
      const rows = rowsOf(table);
      const values = Object.fromEntries(Object.entries(op.values ?? {}).map(([k, v]) => [k, resolve(v, index)]));
      const id = op.op === 'insert' ? undefined : String(resolve(op.id, index));
      const at = rows.findIndex((row) => row.id === id);
      const entry = (rowId: string, prior: Row | null, next: Row | null): void => {
        audit.push([table, {id: rowId, tx_id: tx, op: op.op, actor_id: null, ts, before: prior, after: next}]);
      };
      if (op.op === 'insert') {
        const row = table.stamp(values, 1);
        table.check(row, index);
        rows.push(row);
        if (op.ref !== undefined)
          refs.set(op.ref, row.id as string);
        results[index] = MemoryDomainBackend._stamped(row, 1);
        entry(row.id as string, null, {...row});
      } else if (at < 0)
        throw new DomainBackendError('not-found', `Operation ${index}: no row "${id}"`);
      else if (op.op === 'delete') {
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
        entry(id!, {...rows[at]}, null);
        rows[at] = table.stamp({...rows[at], [IS_DELETED]: true}, (rows[at].version as number) + 1);
        results[index] = {id};
      } else {
        if (op.expectedVersion !== undefined && rows[at].version !== op.expectedVersion) {
          throw new DomainBackendError('version-conflict',
            `Operation ${index}: row "${id}" is at version ${rows[at].version}, expected ${op.expectedVersion}`);
        }
        const row = table.stamp({...rows[at], ...values}, (rows[at].version as number) + 1);
        table.check(row, index);
        entry(id!, {...rows[at]}, {...row});
        rows[at] = row;
        results[index] = MemoryDomainBackend._stamped(row, row.version as number);
      }
    }
    for (const [table, rows] of copies)
      table.rows.splice(0, table.rows.length, ...rows);
    for (const [table, line] of audit)
      table.history.push(line);
    return results;
  }

  /** What an op answers with: the id and version the server answers, plus the system columns the
   * stamp wrote — the platform editor re-reads the row for those, and a writer here has it. */
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

  /** Every `$name` in a value, lists included; `$$` is a literal. */
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
  readonly address: string;
  readonly properties: IProperty[];
  readonly info: DomainTableInfoLike;
  /** The store — what a query copies from and a transaction writes to. */
  readonly rows: Row[];
  /** Every ref column → its target address. */
  readonly refs: Record<string, string>;
  /** What every transaction that touched this table recorded, oldest first. */
  readonly history: AuditEntryLike[] = [];

  private readonly _access: AccessData;
  private readonly _author: string;

  constructor(private readonly _backend: MemoryDomainBackend, readonly schema: string, readonly name: string,
    json: MemoryTableJson, rows: Row[], access?: AccessData, author = 'me') {
    this.address = `${schema}.${name}`;
    this._author = author;
    const columns = Object.entries(json.columns);
    this.properties = [
      ...SYSTEM.map(([column, type]) => MemoryTable._property(column, type,
        {get: (row) => row[column], semType: column === 'author_id' ? 'User' : undefined})),
      ...columns.map(([column, c]) => MemoryTable._column(schema, column, c)),
    ];
    this.refs = Object.fromEntries(columns.filter(([, c]) => c.type === 'ref' && c.ref)
      .map(([column, c]) => [column, MemoryTable._address(schema, c.ref!)]));
    const named = columns.find(([, c]) => c.isName)?.[0] ??
      columns.find(([n, c]) => n === 'name' && c.type === 'string')?.[0];
    const singular = json.singularName ?? json.friendlyName?.replace(/s$/, '') ?? name.replace(/_/g, ' ');
    const searchable = columns.filter(([, c]) => c.searchable === true).map(([n]) => n);
    const permissions = Array.isArray(json.permissions) ? json.permissions : Object.keys(json.permissions ?? {});
    this.info = {
      nameColumn: named ?? null, businessKey: json.businessKey ?? [],
      singularName: singular, pluralName: json.pluralName ?? json.friendlyName ?? `${singular}s`,
      searchableColumns: searchable.length > 0 ? searchable : named === undefined ? [] : [named],
      constraints: Object.entries(json.constraints ?? {}).filter(([, c]) => typeof c.expr === 'string')
        .map(([n, c]) => ({name: n, expr: c.expr!, ...(c.message === undefined ? {} : {message: c.message})})),
      refFilters: Object.fromEntries(columns.filter(([, c]) => c.filter).map(([n, c]) => [n, c.filter!])),
      permissions,
      childTables: [],
    };
    this.rows = rows.map((row) => this.stamp({...row}, 1));
    const granted = Object.fromEntries(permissions.map((p) => [p, true]));
    this._access = access ? {can: {...granted, ...access.can}, fields: access.fields} : {
      can: {view: true, insert: true, edit: true, delete: true, share: true, ...granted},
      fields: Object.fromEntries(this.properties.map((p): [string, FieldAccess] =>
        [p.name!, SYSTEM.some(([column]) => column === p.name) ? 'readonly' : 'editable'])),
    };
  }

  access(): Promise<AccessData> {
    return Promise.resolve({can: {...this._access.can}, fields: {...this._access.fields}});
  }

  async query(spec: DomainQueryLike = {}): Promise<Row[]> {
    let rows = await this._where(spec.filter, spec.search, spec.deleted);
    if (spec.sort)
      rows = MemoryTable._sorted(rows, spec.sort);
    const offset = spec.offset ?? 0;
    rows = rows.slice(offset, spec.limit === undefined ? undefined : offset + spec.limit);
    return rows.map((row) => {
      const out: Row = spec.columns ?
        Object.fromEntries(spec.columns.filter((c) => c in row).map((c) => [c, row[c]])) : {...row};
      // as the server off row mode: edit and delete are the table's answer, share is not
      // carried (null); a seed row carrying its own boolean is a row-mode table
      if (spec.withAccess) {
        for (const [column, capability] of Access.ROW_COLUMNS)
          out[column] = row[column] ?? (capability === 'share' ? null : this._access.can[capability] === true);
      }
      if (spec.deleted !== undefined && spec.deleted !== 'exclude')
        out[Rows.DELETED] = row[IS_DELETED] === true;
      return out;
    });
  }

  /** The rows as a `MemoryFrame` with `MemoryEditState` as its writer — the same shape the dg
   * backend answers with the js-api editor; an appended page lands in the same frame. */
  async frame(spec: DomainQueryLike): Promise<DomainFrameLike> {
    const [rows, access] = await Promise.all([this.query(spec), this.access()]);
    const columns: MemoryColumn[] = [
      ...this.properties.map((p) => ({name: p.name!, type: p.propertyType ?? p.type ?? 'string', semType: p.semType})),
      {name: Rows.STATE, type: 'string'},
      ...(spec.withAccess ? Access.ROW_COLUMNS.map(([name]) => ({name, type: 'bool'})) : []),
      ...(spec.deleted === undefined || spec.deleted === 'exclude' ? [] : [{name: Rows.DELETED, type: 'bool'}]),
    ];
    const df = new MemoryFrame(columns, rows);
    const edit = new MemoryEditState(this, df, Access.from(access));
    return {
      df,
      edit,
      append: async (page) => {
        const more = await this.query(page);
        df.rows.push(...more);
        df.onRowsAdded.fire(undefined);
        return more.length;
      },
      dispose: () => edit.dispose(),
    };
  }

  async count(filter?: DomainQueryLike['filter'], search?: string,
    deleted?: DomainDeletedMode): Promise<number> {
    return (await this._where(filter, search, deleted)).length;
  }

  transaction(ops: DomainTransactionOpLike[]): Promise<DomainTransactionResultLike[]> {
    return this._backend.transaction(ops);
  }

  audit(id: string): Promise<AuditEntryLike[]> {
    return Promise.resolve(this.history.filter((line) => line.id === id));
  }

  /** The soft delete undone, as the server's `POST …/{id}/restore`: a row whose ref column points
   * at a row still in the trash is refused naming that column — restoring it would leave a live
   * row referring to a deleted one. Restoring a live row does nothing. */
  async restore(id: string): Promise<void> {
    const at = this.rows.findIndex((row) => row.id === id);
    if (at < 0)
      throw new DomainBackendError('not-found', `${this.address}: no row "${id}"`);
    const row = this.rows[at];
    if (row[IS_DELETED] !== true)
      return;
    for (const [column, target] of Object.entries(this.refs)) {
      const parent = this._backend.tableSync(target)?.rows.find((r) => r.id === row[column]);
      if (parent?.[IS_DELETED] === true) {
        throw new DomainBackendError('validation',
          `Cannot restore row "${id}": ${column} refers to the deleted ${target} "${parent.id}"`);
      }
    }
    const restored = this.stamp({...row, [IS_DELETED]: false}, (row.version as number) + 1);
    this.rows[at] = restored;
    this.history.push({id, tx_id: uuid4(), op: 'undelete', actor_id: null,
      ts: restored.updated_on as string, before: {...row}, after: {...restored}});
  }

  /** The schema's rules — required, choices, min, max — as the server's `_validateRow` refuses on. */
  check(row: Row, index: number): void {
    for (const prop of this.properties) {
      const problem = MemoryEditState.problemOf(prop, row[prop.name!]);
      if (problem !== null)
        throw new DomainBackendError('validation', `Operation ${index}: column "${prop.name}": ${problem}`);
    }
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

  private async _where(filter: DomainQueryLike['filter'], search?: string,
    deleted: DomainDeletedMode = 'exclude'): Promise<Row[]> {
    let rows = this.rows;
    if (filter !== undefined && filter !== '') {
      const mask = await Filters.toMask(this._frameLike(), MemoryTable._tree(filter));
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

  private static _tree(filter: NonNullable<DomainQueryLike['filter']>): FilterGroup {
    if (typeof filter !== 'string')
      return Filters.fromDomainTree(filter);
    const {root, problems} = Filters.parse(filter);
    if (problems.length > 0)
      throw new DomainBackendError('validation', problems[0].message);
    return root;
  }

  /** The store as `toMask` reads a frame: raw arrays per column, built on demand. */
  private _frameLike(): MaskFrameLike {
    const rows = this.rows;
    const columns = new Map<string, MaskColumnLike>();
    return {
      rowCount: rows.length,
      column: (name) => {
        const prop = this.properties.find((p) => p.name === name);
        if (prop === undefined)
          return null;
        let column = columns.get(name);
        if (column === undefined)
          columns.set(name, column = MemoryTable._maskColumn(name, prop.propertyType ?? prop.type ?? 'string', rows));
        return column;
      },
    };
  }

  private static _maskColumn(name: string, type: string, rows: Row[]): MaskColumnLike {
    const n = rows.length;
    const cell = (i: number) => rows[i][name];
    const isNull = (v: unknown) => v === null || v === undefined;
    const number = (v: unknown, nil: number) => isNull(v) ? nil : Number(v);
    const column = (raw: () => ArrayLike<number>, extra: Partial<MaskColumnLike> = {}): MaskColumnLike =>
      ({name, type, length: n, getRawData: raw, ...extra});
    switch (Filters.kindOf({name, type})) {
      case Filters.KIND.INT:
        return column(() => Int32Array.from(rows, (r) => number(r[name], INT_NULL)));
      case Filters.KIND.FLOAT:
        return column(() => Float64Array.from(rows, (r) => number(r[name], FLOAT_NULL)));
      case Filters.KIND.DATE_TIME:
        return column(() => Float64Array.from(rows, (r) => {
          const v = r[name];
          return isNull(v) ? FLOAT_NULL : (v instanceof Date ? v.getTime() : Date.parse(String(v))) * 1000;
        }));
      case Filters.KIND.BOOL:
        return column(() => BitArray.create(n, (i) => cell(i) === true).getBuffer());
      case Filters.KIND.BIG_INT: case Filters.KIND.STRING_LIST:
        return column(() => new Int32Array(0), {get: cell});
      default: {
        const text = (v: unknown) => isNull(v) ? '' : String(v);
        const categories = [...new Set(rows.map((r) => text(r[name])))];
        const index = new Map(categories.map((c, i) => [c, i]));
        return column(() => Int32Array.from(rows, (r) => index.get(text(r[name]))!), {categories});
      }
    }
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
