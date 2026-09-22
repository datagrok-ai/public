/* The manifest editor's model: the draft envelope the server answers (`POST /domains/schemas/draft`
   → `{manifest, inventory, diagnostics}`) plus the user's edits, projected back into the manifest
   to submit by `toJSON()`. Every table and column is identified by its REMOTE name — the logical
   name is what the user edits — and every declaration key the editor does not expose is carried
   through untouched (Astra 4: deep preservation). Refs follow one RULE, never an edit: a
   relation the warehouse reported as `ref` becomes `{type: 'ref', ref}` while its target table
   is included, and the column keeps its scalar type otherwise. The rules the server checks are
   mirrored in `ManifestRules`; the dry run remains the authority. */
import {computed, signal, Signal, ReadonlySignal} from '../../../core/signals.js';
import {ManifestRules} from './manifest-rules.js';

export interface ManifestColumnJson {
  type: string;
  /** The remote column, where it differs from the logical name. */
  column?: string;
  ref?: string;
  required?: boolean;
  isName?: boolean;
  searchable?: boolean;
  friendlyName?: string;
  [key: string]: unknown;
}

export interface ManifestTableJson {
  /** The remote table, where it differs from the logical name. */
  table?: string;
  businessKey?: string[];
  friendlyName?: string;
  writable?: boolean;
  columns: Record<string, ManifestColumnJson>;
  [key: string]: unknown;
}

export interface ManifestStorageJson {
  kind: string;
  connection?: string;
  schema?: string;
  catalog?: string;
  writable?: boolean;
}

export interface ManifestJson {
  name: string;
  version?: string;
  storage?: ManifestStorageJson;
  tables: Record<string, ManifestTableJson>;
  [key: string]: unknown;
}

/** A remote table the draft looked at: bindable ones carry their key, the rest the reason. */
export interface InventoryTable {
  remote: string;
  logical?: string;
  bindable: boolean;
  key?: string[];
  code?: string;
  message?: string;
}

/** A remote column the draft reports on — one it could not bind (`code`), or one it only
 * annotates with its warehouse type. */
export interface InventoryColumn {
  table: string;
  remote: string;
  dbType?: string;
  code?: string;
  message?: string;
}

/** A foreign key the warehouse reported, and whether the draft made a ref of it. */
export interface InventoryRelation {
  table: string;
  column: string;
  targetTable: string;
  targetColumn: string;
  status: 'ref' | 'plain';
  code?: string;
  message?: string;
}

export interface DraftInventory {
  tables?: InventoryTable[];
  columns?: InventoryColumn[];
  relations?: InventoryRelation[];
}

export interface ManifestDiagnostic {
  code: string;
  message: string;
  /** A manifest path (`tables.order.columns.ship_via`); absent for a schema-wide finding. */
  path?: string;
}

export interface DraftEnvelope {
  manifest: ManifestJson;
  inventory?: DraftInventory;
  diagnostics?: ManifestDiagnostic[];
}

export type ManifestSelection =
  | {kind: 'schema'}
  | {kind: 'table', table: string}
  | {kind: 'column', table: string, column: string};

export interface TableView {
  remote: string;
  logical: string;
  friendlyName: string;
  included: boolean;
  /** False for a view, a keyless table, a table whose key the platform cannot hold — or one the
   * draft holds no declaration for. */
  bindable: boolean;
  reason?: string;
  code?: string;
  /** The key's remote column names. */
  key: string[];
  readOnly: boolean;
}

export interface RelationView {
  table: string;
  column: string;
  targetTable: string;
  targetColumn: string;
  /** What the warehouse and the draft said. */
  status: 'ref' | 'plain';
  /** Whether the column is a ref right now — the rule's answer. */
  ref: boolean;
  targetIncluded: boolean;
  targetLogical: string | null;
  reason: string;
  /** Including the target would make it a ref. */
  canFix: boolean;
}

export interface ColumnView {
  table: string;
  remote: string;
  logical: string;
  /** The type the manifest will carry: `ref` under the rule, the scalar type otherwise. */
  type: string;
  scalarType: string;
  dbType?: string;
  included: boolean;
  isKey: boolean;
  required: boolean;
  isName: boolean;
  searchable: boolean;
  supported: boolean;
  reason?: string;
  code?: string;
  relation?: RelationView;
}

interface ColumnState {
  remote: string;
  logical: string;
  decl: ManifestColumnJson;
  included: boolean;
}

interface TableState {
  remote: string;
  logical: string;
  decl: ManifestTableJson | null;
  included: boolean;
  key: string[];
  columns: Map<string, ColumnState>;
  inventory?: InventoryTable;
}

export class ManifestModel {
  /** The manifest's `name` — what the schema is registered as. */
  readonly name: Signal<string>;
  /** Schema metadata the create envelope carries beside the manifest (Astra 9). */
  readonly friendlyName: Signal<string>;
  readonly tables: ReadonlySignal<TableView[]>;
  readonly relations: ReadonlySignal<RelationView[]>;
  readonly writable: ReadonlySignal<boolean>;
  readonly manifest: ReadonlySignal<ManifestJson>;
  /** Bumped by every edit — what a view re-renders on without rebuilding the manifest. */
  readonly revision: ReadonlySignal<number>;

  private readonly _rev = signal(0);
  private readonly _draft: ManifestJson;
  private readonly _tables = new Map<string, TableState>();
  private readonly _order: string[] = [];
  private readonly _unsupported = new Map<string, InventoryColumn>();
  private readonly _dbTypes = new Map<string, string>();
  private readonly _relations: InventoryRelation[];
  private readonly _columnViews = new Map<string, ReadonlySignal<ColumnView[]>>();
  private _writable: boolean;

  constructor(draft: DraftEnvelope, options: {friendlyName?: string} = {}) {
    this._draft = draft.manifest;
    this.name = signal(draft.manifest.name);
    this.friendlyName = signal(options.friendlyName ?? '');
    this.revision = this._rev;
    this._writable = draft.manifest.storage?.writable === true;
    this._relations = draft.inventory?.relations ?? [];
    for (const [logical, decl] of Object.entries(draft.manifest.tables))
      this._addTable(decl.table ?? logical, logical, decl);
    for (const t of draft.inventory?.tables ?? []) {
      const state = this._tables.get(t.remote);
      if (state === undefined)
        this._addTable(t.remote, t.logical ?? t.remote, null, t);
      else
        state.inventory = t;
    }
    for (const c of draft.inventory?.columns ?? []) {
      const key = ManifestModel._key(c.table, c.remote);
      if (c.dbType !== undefined)
        this._dbTypes.set(key, c.dbType);
      if (c.code !== undefined)
        this._unsupported.set(key, c);
    }
    this.tables = computed(() => {
      this._rev.value;
      return this._order.map((remote) => this._tableView(this._tables.get(remote)!));
    });
    this.relations = computed(() => {
      this._rev.value;
      return this._relations.map((r) => this._relationView(r));
    });
    this.writable = computed(() => {
      this._rev.value;
      return this._writable;
    });
    this.manifest = computed(() => {
      this._rev.value;
      return this._toJSON(this.name.value);
    });
  }

  /** The storage the draft was made over — read-only here, the vocabulary of the editor. */
  get storage(): ManifestStorageJson | undefined {
    return this._draft.storage;
  }

  table(remote: string): TableView | undefined {
    const state = this._tables.get(remote);
    return state === undefined ? undefined : this._tableView(state);
  }

  /** The table's columns in declaration order, the unsupported ones after. */
  columns(table: string): ReadonlySignal<ColumnView[]> {
    let view = this._columnViews.get(table);
    if (view === undefined) {
      view = computed(() => {
        this._rev.value;
        const state = this._tables.get(table);
        if (state === undefined)
          return [];
        const views = [...state.columns.values()].map((c) => this._columnView(state, c));
        for (const u of this._unsupported.values()) {
          if (u.table === table && !state.columns.has(u.remote))
            views.push(this._unsupportedView(state, u));
        }
        return views;
      });
      this._columnViews.set(table, view);
    }
    return view;
  }

  column(table: string, remote: string): ColumnView | undefined {
    return this.columns(table).peek().find((c) => c.remote === remote);
  }

  /** The relation a column carries, if the warehouse reported one. */
  relationOf(table: string, column: string): RelationView | undefined {
    return this.relations.peek().find((r) => r.table === table && r.column === column);
  }

  includeTable(remote: string, on: boolean): void {
    const state = this._tables.get(remote);
    if (state === undefined || state.decl === null || state.included === on)
      return;
    this._mutate(() => state.included = on);
  }

  /** Every bindable table in, or every table out. */
  includeTables(on: boolean): void {
    this._mutate(() => {
      for (const state of this._tables.values()) {
        if (state.decl !== null)
          state.included = on;
      }
    });
  }

  /** A key column stays: the row id encodes it. */
  includeColumn(table: string, remote: string, on: boolean): void {
    const state = this._tables.get(table);
    const column = state?.columns.get(remote);
    if (state === undefined || column === undefined || state.key.includes(remote) || column.included === on)
      return;
    this._mutate(() => column.included = on);
  }

  checkSchemaName(name: string): string | null {
    return ManifestRules.checkSchemaName(name);
  }

  checkTableName(remote: string, logical: string): string | null {
    const problem = ManifestRules.checkTableName(logical);
    if (problem !== null)
      return problem;
    for (const other of this._tables.values()) {
      if (other.remote !== remote && other.decl !== null && other.logical === logical)
        return `Another table is already named "${logical}"`;
    }
    return null;
  }

  checkColumnName(table: string, remote: string, logical: string): string | null {
    const problem = ManifestRules.checkColumnName(logical);
    if (problem !== null)
      return problem;
    const state = this._tables.get(table);
    for (const other of state?.columns.values() ?? []) {
      if (other.remote !== remote && other.logical === logical)
        return `Another column of ${state!.logical} is already named "${logical}"`;
    }
    return null;
  }

  /** Applies the name when it passes, else answers the problem and keeps the name it had. */
  renameTable(remote: string, logical: string): string | null {
    const state = this._tables.get(remote);
    if (state === undefined)
      return null;
    const problem = this.checkTableName(remote, logical);
    if (problem === null && state.logical !== logical)
      this._mutate(() => state.logical = logical);
    return problem;
  }

  renameColumn(table: string, remote: string, logical: string): string | null {
    const column = this._tables.get(table)?.columns.get(remote);
    if (column === undefined)
      return null;
    const problem = this.checkColumnName(table, remote, logical);
    if (problem === null && column.logical !== logical)
      this._mutate(() => column.logical = logical);
    return problem;
  }

  setFriendlyName(table: string, name: string): void {
    const decl = this._tables.get(table)?.decl;
    if (decl === undefined || decl === null)
      return;
    this._mutate(() => {
      if (name === '')
        delete decl.friendlyName;
      else
        decl.friendlyName = name;
    });
  }

  /** At most one name column per table; null clears it. Only a string column names a row. */
  setNameColumn(table: string, remote: string | null): void {
    this._setSingle(table, remote, 'isName');
  }

  /** At most one searchable column per table; null clears it. */
  setSearchable(table: string, remote: string | null): void {
    this._setSingle(table, remote, 'searchable');
  }

  /** A key column is always required. */
  setRequired(table: string, remote: string, on: boolean): void {
    const state = this._tables.get(table);
    const column = state?.columns.get(remote);
    if (state === undefined || column === undefined || state.key.includes(remote))
      return;
    this._mutate(() => ManifestModel._flag(column.decl, 'required', on));
  }

  /** Opts one table out of writes under a writable storage. */
  setReadOnly(table: string, on: boolean): void {
    const decl = this._tables.get(table)?.decl;
    if (decl === undefined || decl === null)
      return;
    this._mutate(() => {
      if (on)
        decl.writable = false;
      else
        delete decl.writable;
    });
  }

  setWritable(on: boolean): void {
    if (this._writable !== on)
      this._mutate(() => this._writable = on);
  }

  /** The manifest to submit: the draft's declarations with the edits applied and the rule's refs. */
  toJSON(): ManifestJson {
    return this.manifest.peek();
  }

  /** The node a manifest path names (`tables.<t>`, `tables.<t>.columns.<c>`, `tables.<t>.<key>`),
   * by the CURRENT logical names; the schema for a root path or one no node owns. */
  resolvePath(path: string | undefined): ManifestSelection {
    const parts = (path ?? '').split('.');
    if (parts[0] !== 'tables' || parts.length < 2)
      return {kind: 'schema'};
    const table = [...this._tables.values()].find((t) => t.decl !== null && t.logical === parts[1]);
    if (table === undefined)
      return {kind: 'schema'};
    if (parts[2] === 'columns' && parts.length >= 4) {
      const column = [...table.columns.values()].find((c) => c.logical === parts[3]);
      if (column !== undefined)
        return {kind: 'column', table: table.remote, column: column.remote};
    }
    return {kind: 'table', table: table.remote};
  }

  static sameSelection(a: ManifestSelection, b: ManifestSelection): boolean {
    return a.kind === b.kind && (a.kind === 'schema' ||
      (a.table === (b as {table: string}).table &&
        (a.kind === 'table' || a.column === (b as {column: string}).column)));
  }

  private _addTable(remote: string, logical: string, decl: ManifestTableJson | null,
    inventory?: InventoryTable): void {
    const columns = new Map<string, ColumnState>();
    const keyByLogical = new Map<string, string>();
    for (const [name, c] of Object.entries(decl?.columns ?? {})) {
      const columnRemote = c.column ?? name;
      keyByLogical.set(name, columnRemote);
      columns.set(columnRemote, {remote: columnRemote, logical: name, decl: {...c}, included: true});
    }
    const key = decl?.businessKey?.map((k) => keyByLogical.get(k) ?? k) ?? inventory?.key ?? [];
    this._tables.set(remote, {remote, logical, decl: decl === null ? null : {...decl}, included: decl !== null,
      key, columns, inventory});
    this._order.push(remote);
  }

  private _mutate(edit: () => void): void {
    edit();
    this._rev.value = this._rev.peek() + 1;
  }

  private _setSingle(table: string, remote: string | null, flag: 'isName' | 'searchable'): void {
    const state = this._tables.get(table);
    if (state === undefined)
      return;
    const target = remote === null ? undefined : state.columns.get(remote);
    if (target !== undefined && this._columnView(state, target).type !== 'string')
      return;
    this._mutate(() => {
      for (const column of state.columns.values())
        ManifestModel._flag(column.decl, flag, column === target);
    });
  }

  private static _flag(decl: ManifestColumnJson, flag: 'required' | 'isName' | 'searchable', on: boolean): void {
    if (on)
      decl[flag] = true;
    else
      delete decl[flag];
  }

  private _tableView(state: TableState): TableView {
    const bindable = state.decl !== null;
    const inventory = state.inventory;
    return {
      remote: state.remote, logical: state.logical,
      friendlyName: state.decl?.friendlyName ?? '',
      included: bindable && state.included, bindable,
      reason: bindable ? undefined : inventory?.message ?? 'the draft holds no declaration for this table',
      code: bindable ? undefined : inventory?.code,
      key: state.key, readOnly: state.decl?.writable === false,
    };
  }

  private _relationView(r: InventoryRelation): RelationView {
    const target = this._tables.get(r.targetTable);
    const targetIncluded = target !== undefined && target.decl !== null && target.included;
    const ref = r.status === 'ref' && targetIncluded;
    const canFix = r.status === 'ref' && !targetIncluded && target !== undefined && target.decl !== null;
    const reason = ref ? `ref: ${target!.logical}` :
      r.status === 'ref' ? 'target not included — stays a plain value' :
        r.message ?? 'stays a plain value';
    return {table: r.table, column: r.column, targetTable: r.targetTable, targetColumn: r.targetColumn,
      status: r.status, ref, targetIncluded, targetLogical: target?.logical ?? null, reason, canFix};
  }

  private _columnView(state: TableState, column: ColumnState): ColumnView {
    const relation = this._relations.find((r) => r.table === state.remote && r.column === column.remote);
    const view = relation === undefined ? undefined : this._relationView(relation);
    const scalarType = this._scalarType(state, column);
    const isKey = state.key.includes(column.remote);
    return {
      table: state.remote, remote: column.remote, logical: column.logical,
      type: view?.ref === true ? 'ref' : scalarType, scalarType,
      dbType: this._dbTypes.get(ManifestModel._key(state.remote, column.remote)),
      included: column.included, isKey,
      required: isKey || column.decl.required === true,
      isName: column.decl.isName === true, searchable: column.decl.searchable === true,
      supported: true, relation: view,
    };
  }

  private _unsupportedView(state: TableState, u: InventoryColumn): ColumnView {
    return {
      table: state.remote, remote: u.remote, logical: u.remote, type: u.dbType ?? '', scalarType: u.dbType ?? '',
      dbType: u.dbType, included: false, isKey: false, required: false, isName: false, searchable: false,
      supported: false, reason: u.message, code: u.code,
    };
  }

  /** The scalar type of a column — for a drafted ref, the type of the key it points at. */
  private _scalarType(state: TableState, column: ColumnState): string {
    const decl = column.decl;
    if (decl.type !== 'ref')
      return decl.type;
    const relation = this._relations.find((r) => r.table === state.remote && r.column === column.remote);
    const target = relation === undefined ? undefined : this._tables.get(relation.targetTable);
    const key = target === undefined || relation === undefined ? undefined : target.columns.get(relation.targetColumn);
    return key === undefined || key.decl.type === 'ref' ? 'string' : key.decl.type;
  }

  private _toJSON(name: string): ManifestJson {
    const tables: Record<string, ManifestTableJson> = {};
    for (const remote of this._order) {
      const state = this._tables.get(remote)!;
      if (state.decl === null || !state.included)
        continue;
      const columns: Record<string, ManifestColumnJson> = {};
      for (const column of state.columns.values()) {
        if (column.included)
          columns[column.logical] = this._columnJSON(state, column);
      }
      const {table: _table, businessKey: _key, writable, filters, delegate, ...rest} = state.decl;
      const decl: ManifestTableJson = {...rest, columns};
      if (state.logical !== state.remote)
        decl.table = state.remote;
      if (state.key.length > 0)
        decl.businessKey = state.key.map((k) => ManifestModel._current(state, k)!);
      const named = (draftName: string): string | null =>
        ManifestModel._current(state, state.decl!.columns[draftName]?.column ?? draftName);
      if (Array.isArray(filters)) {
        const kept: unknown[] = [];
        for (const f of filters as {column?: unknown}[]) {
          const [head, ...tail] = String(f.column ?? '').split('.');
          const current = named(head);
          if (current !== null)
            kept.push({...f, column: [current, ...tail].join('.')});
        }
        if (kept.length > 0)
          decl.filters = kept;
      }
      const delegateNow = typeof delegate === 'string' ? named(delegate) : null;
      if (delegateNow !== null)
        decl.delegate = delegateNow;
      if (this._writable && writable === false)
        decl.writable = false;
      tables[state.logical] = decl;
    }
    const {storage, tables: _tables, name: _name, ...rest} = this._draft;
    const json: ManifestJson = {...rest, name, tables};
    if (storage !== undefined) {
      const {writable: _writable, ...storageRest} = storage;
      json.storage = this._writable ? {...storageRest, writable: true} : storageRest;
    }
    return json;
  }

  private _columnJSON(state: TableState, column: ColumnState): ManifestColumnJson {
    const {type: _type, ref: _ref, column: _column, ...rest} = column.decl;
    const relation = this._relations.find((r) => r.table === state.remote && r.column === column.remote);
    const view = relation === undefined ? undefined : this._relationView(relation);
    const json: ManifestColumnJson = view?.ref === true ? {type: 'ref', ref: view.targetLogical!, ...rest} :
      {type: this._scalarType(state, column), ...rest};
    if (column.logical !== column.remote)
      json.column = column.remote;
    return json;
  }

  /** The logical name a declaration key naming a column carries out: the current one, the name
   * itself where no column is declared under it, null once the column is excluded. */
  private static _current(state: TableState, remote: string): string | null {
    const column = state.columns.get(remote);
    return column === undefined ? remote : column.included ? column.logical : null;
  }

  private static _key(table: string, column: string): string {
    return `${table}\u0000${column}`;
  }
}

/** One row of the access grid: who, and which of View / Edit / Delete. */
export interface AccessGrant {
  /** `'schema'` — every table — or a table's remote name. */
  scope: string;
  group: string;
  view: boolean;
  edit: boolean;
  delete: boolean;
}

/** Who may see a column: everyone who may see the row (null), or only these groups. */
export interface ColumnVisibility {
  table: string;
  column: string;
  groups: string[] | null;
}

export interface AccessJson {
  grants: AccessGrant[];
  visibility: ColumnVisibility[];
}

export type AccessCapability = 'view' | 'edit' | 'delete';

/** Access is not in the manifest — declared grants are refused for a user schema — so the dialog
 * applies what this holds after Create through the grants API and the column restrictions.
 * Tables and columns are keyed by REMOTE name, as the manifest model keys them; the editor's
 * `plan()` resolves them to the logical names the API takes. */
export class AccessModel {
  readonly grants: Signal<AccessGrant[]>;
  readonly visibility: Signal<ColumnVisibility[]>;

  constructor(json: Partial<AccessJson> = {}) {
    this.grants = signal(json.grants?.map((g) => ({...g})) ?? []);
    this.visibility = signal(json.visibility?.map((v) =>
      ({...v, groups: v.groups === null ? null : [...v.groups]})) ?? []);
  }

  grantsOf(scope: string): ReadonlySignal<AccessGrant[]> {
    return computed(() => this.grants.value.filter((g) => g.scope === scope));
  }

  /** Replaces the scope's rows wholesale — what the access grid hands back. */
  setGrants(scope: string, rows: Omit<AccessGrant, 'scope'>[]): void {
    this.grants.value = [...this.grants.peek().filter((g) => g.scope !== scope),
      ...rows.map((r) => ({...r, scope}))];
  }

  addGroup(scope: string, group: string): void {
    if (this.grants.peek().some((g) => g.scope === scope && g.group === group))
      return;
    this.grants.value = [...this.grants.peek(), {scope, group, view: true, edit: false, delete: false}];
  }

  removeGroup(scope: string, group: string): void {
    this.grants.value = this.grants.peek().filter((g) => !(g.scope === scope && g.group === group));
  }

  setGrant(scope: string, group: string, capability: AccessCapability, on: boolean): void {
    this.grants.value = this.grants.peek().map((g) =>
      g.scope === scope && g.group === group ? {...g, [capability]: on} : g);
  }

  visibilityOf(table: string, column: string): string[] | null {
    return this.visibility.peek().find((v) => v.table === table && v.column === column)?.groups ?? null;
  }

  setVisibility(table: string, column: string, groups: string[] | null): void {
    const rest = this.visibility.peek().filter((v) => !(v.table === table && v.column === column));
    this.visibility.value = groups === null ? rest : [...rest, {table, column, groups: [...groups]}];
  }

  toJSON(): AccessJson {
    return {grants: this.grants.peek().map((g) => ({...g})),
      visibility: this.visibility.peek().map((v) => ({...v, groups: v.groups === null ? null : [...v.groups]}))};
  }
}
