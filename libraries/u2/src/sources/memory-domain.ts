/* `DomainBackend` over a schema.json held in memory (GOAL step-back "A local in-memory domain
   backend"): the same seam the platform fills, answered without a server — the gallery, the
   headless tests, an agent iterating on an app. Filters are evaluated with the filter feature's
   own `toMask`, so a query means here what it means on the server. */
import {BitArray} from 'datagrok-api/u2core';
import {Filters} from '../core/filter/index.js';
import type {FilterGroup} from '../core/filter/model.js';
import {INT_NULL, FLOAT_NULL} from '../core/filter/evaluate.js';
import type {MaskColumnLike, MaskFrameLike} from '../core/filter/evaluate.js';
import type {IProperty} from '../core/property-like.js';
import type {AccessData, FieldAccess} from '../core/access.js';
import {Access} from '../core/access.js';
import {DomainBackendError} from './domain-backend.js';
import type {DomainBackend, DomainFrameLike, DomainQueryLike, DomainTableInfoLike, DomainTableLike,
  DomainTransactionOpLike, DomainTransactionResultLike} from './domain-backend.js';
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
  semType?: string;
  default?: unknown;
}

export interface MemoryTableJson {
  columns: Record<string, MemoryColumnJson>;
  businessKey?: string[];
  friendlyName?: string;
  singularName?: string;
  pluralName?: string;
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
   * editable, the system columns readonly). */
  access?: AccessData;
  /** The user id stamped into `author_id` on insert (default `'me'`). */
  author?: string;
}

/** The system columns every table carries (js-api `DOMAIN_SYSTEM_COLUMNS`). */
const SYSTEM: [string, string][] = [['id', 'string'], ['version', 'int'], ['created_on', 'datetime'],
  ['updated_on', 'datetime'], ['author_id', 'string']];

/** Column type → property type, as the registry's `rowProperties` answers them. */
const PROPERTY_TYPES: Record<string, string> = {
  string: 'string', int: 'int', float: 'double', bool: 'bool', datetime: 'datetime',
  string_list: 'string_list', ref: 'string', user: 'string', group: 'string', file: 'string', json: 'map',
};

const CORE_REF_SEM_TYPES: Record<string, string> = {user: 'User', group: 'Group'};

export class MemoryDomainBackend implements DomainBackend {
  private readonly _tables = new Map<string, MemoryTable>();

  constructor(schema: MemorySchemaJson, options: MemoryDomainOptions = {}) {
    for (const [name, table] of Object.entries(schema.tables)) {
      this._tables.set(`${schema.name}.${name}`,
        new MemoryTable(schema.name, name, table, options.rows?.[name] ?? [], options.access, options.author));
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
}

export class MemoryTable implements DomainTableLike {
  readonly address: string;
  readonly properties: IProperty[];
  readonly info: DomainTableInfoLike;
  /** The store — what a query copies from and a transaction writes to. */
  readonly rows: Record<string, unknown>[];

  private readonly _access: AccessData;
  private readonly _author: string;

  constructor(schema: string, name: string, json: MemoryTableJson, rows: Record<string, unknown>[],
    access?: AccessData, author = 'me') {
    this.address = `${schema}.${name}`;
    this._author = author;
    const columns = Object.entries(json.columns);
    this.properties = [
      ...SYSTEM.map(([column, type]) => MemoryTable._property(column, type,
        {get: (row) => row[column], semType: column === 'author_id' ? 'User' : undefined})),
      ...columns.map(([column, c]) => MemoryTable._column(schema, column, c)),
    ];
    const named = columns.find(([, c]) => c.isName)?.[0] ??
      columns.find(([n, c]) => n === 'name' && c.type === 'string')?.[0];
    const singular = json.singularName ?? json.friendlyName?.replace(/s$/, '') ?? name.replace(/_/g, ' ');
    this.info = {nameColumn: named ?? null, businessKey: json.businessKey ?? [],
      singularName: singular, pluralName: json.pluralName ?? json.friendlyName ?? `${singular}s`};
    this.rows = rows.map((row) => this._stamp({...row}, 1));
    this._access = access ?? {
      can: {view: true, insert: true, edit: true, delete: true, share: true},
      fields: Object.fromEntries(this.properties.map((p): [string, FieldAccess] =>
        [p.name!, SYSTEM.some(([column]) => column === p.name) ? 'readonly' : 'editable'])),
    };
  }

  access(): Promise<AccessData> {
    return Promise.resolve({can: {...this._access.can}, fields: {...this._access.fields}});
  }

  async query(spec: DomainQueryLike = {}): Promise<Record<string, unknown>[]> {
    let rows = await this._where(spec.filter);
    if (spec.sort)
      rows = MemoryTable._sorted(rows, spec.sort);
    const offset = spec.offset ?? 0;
    rows = rows.slice(offset, spec.limit === undefined ? undefined : offset + spec.limit);
    return rows.map((row) => {
      const out: Record<string, unknown> = spec.columns ?
        Object.fromEntries(spec.columns.filter((c) => c in row).map((c) => [c, row[c]])) : {...row};
      // as the server off row mode: edit and delete are the table's answer, share is not
      // carried (null); a seed row carrying its own boolean is a row-mode table
      if (spec.withAccess) {
        for (const [column, capability] of Access.ROW_COLUMNS)
          out[column] = row[column] ?? (capability === 'share' ? null : this._access.can[capability] === true);
      }
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

  async count(filter?: DomainQueryLike['filter']): Promise<number> {
    return (await this._where(filter)).length;
  }

  /** All or nothing: every op is applied to a copy of the store, which replaces it at the end. */
  async transaction(ops: DomainTransactionOpLike[]): Promise<DomainTransactionResultLike[]> {
    const rows = this.rows.map((row) => ({...row}));
    const refs = new Map<string, string>();
    const declared = new Set(ops.map((op) => op.ref).filter((ref): ref is string => ref !== undefined));
    const resolve = (v: unknown, index: number): unknown => {
      if (Array.isArray(v))
        return v.map((x) => resolve(x, index));
      if (typeof v !== 'string' || !v.startsWith('$'))
        return v;
      if (v.startsWith('$$'))
        return v.slice(1);
      const name = v.slice(1);
      if (refs.has(name))
        return refs.get(name);
      throw new DomainBackendError('bad-ref', declared.has(name) ?
        `Operation ${index}: forward reference "$${name}" — refs may only point to earlier operations` :
        `Operation ${index}: unknown reference "$${name}"`);
    };
    const results: DomainTransactionResultLike[] = [];
    for (const [index, op] of ops.entries()) {
      const values = Object.fromEntries(Object.entries(op.values ?? {}).map(([k, v]) => [k, resolve(v, index)]));
      const at = rows.findIndex((row) => row.id === op.id);
      if (op.op === 'insert') {
        const row = this._stamp(values, 1);
        this._check(row, index);
        rows.push(row);
        if (op.ref !== undefined)
          refs.set(op.ref, row.id as string);
        results.push({id: row.id as string, version: 1});
      } else if (at < 0)
        throw new DomainBackendError('not-found', `Operation ${index}: no row "${op.id}"`);
      else if (op.op === 'delete') {
        rows.splice(at, 1);
        results.push({id: op.id});
      } else {
        if (op.expectedVersion !== undefined && rows[at].version !== op.expectedVersion) {
          throw new DomainBackendError('version-conflict',
            `Operation ${index}: row "${op.id}" is at version ${rows[at].version}, expected ${op.expectedVersion}`);
        }
        const row = this._stamp({...rows[at], ...values}, (rows[at].version as number) + 1);
        this._check(row, index);
        rows[at] = row;
        results.push({id: op.id, version: row.version as number});
      }
    }
    this.rows.splice(0, this.rows.length, ...rows);
    return results;
  }

  private async _where(filter: DomainQueryLike['filter']): Promise<Record<string, unknown>[]> {
    if (filter === undefined || filter === '')
      return this.rows;
    const root = MemoryTable._tree(filter);
    const mask = await Filters.toMask(this._frameLike(), root);
    return this.rows.filter((_, i) => mask.get(i));
  }

  private static _tree(filter: NonNullable<DomainQueryLike['filter']>): FilterGroup {
    if (typeof filter !== 'string')
      return Filters.fromDomainTree(filter);
    const {root, problems} = Filters.parse(filter);
    if (problems.length > 0)
      throw new DomainBackendError('validation', problems[0].message);
    return root;
  }

  /** The schema's rules — required, choices, min, max — as the server's `_validateRow` refuses on. */
  private _check(row: Record<string, unknown>, index: number): void {
    for (const prop of this.properties) {
      const problem = MemoryEditState.problemOf(prop, row[prop.name!]);
      if (problem !== null)
        throw new DomainBackendError('validation', `Operation ${index}: column "${prop.name}": ${problem}`);
    }
  }

  private _stamp(row: Record<string, unknown>, version: number): Record<string, unknown> {
    const now = new Date().toISOString();
    row.id ??= crypto.randomUUID();
    row.version = version;
    row.created_on ??= now;
    row.updated_on = now;
    row.author_id ??= this._author;
    return row;
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

  private static _maskColumn(name: string, type: string, rows: Record<string, unknown>[]): MaskColumnLike {
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

  private static _sorted(rows: Record<string, unknown>[], sort: string): Record<string, unknown>[] {
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

  private static _column(schema: string, name: string, c: MemoryColumnJson): IProperty {
    const ref = c.type === 'ref' && c.ref ? (c.ref.includes('.') ? c.ref : `${schema}.${c.ref}`) : undefined;
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
