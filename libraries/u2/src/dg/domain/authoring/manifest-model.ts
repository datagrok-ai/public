/* The manifest editor's model: the draft envelope the server answers (`POST /domains/schemas/draft`
   → `{manifest, inventory, diagnostics}`) plus the user's edits, projected back into the manifest
   to submit by `toJSON()`. Every table and column is identified by its REMOTE name — the logical
   name is what the user edits — and every declaration key the editor does not expose is carried
   through untouched (Astra 4: deep preservation). A declared type — `ref` included — is the
   authority; the inventory only explains, and is consulted when the USER changes inclusion:
   excluding a ref's target demotes the ref to the plain type the inventory knows (no such type:
   the column drops out with a reason), including a target promotes back the refs the inventory
   stamped `ref`. Remote names are vocabulary of the external storage alone. The rules the server
   checks are mirrored in `ManifestRules`; the dry run remains the authority. */
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
  /** The platform type the draft mapped the warehouse type to — what a demoted ref becomes. */
  type?: string;
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
  /** Whether the column is a ref onto this target right now. */
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
  /** The type the manifest will carry — as declared, until an inclusion change demoted or
   * promoted it; '' for a column dropped by a demotion. */
  type: string;
  dbType?: string;
  included: boolean;
  isKey: boolean;
  required: boolean;
  isName: boolean;
  searchable: boolean;
  supported: boolean;
  /** Why the column cannot be included: unsupported, or a demoted ref no plain type serves. */
  reason?: string;
  code?: string;
  relation?: RelationView;
}

interface ColumnState {
  remote: string;
  logical: string;
  /** The declaration as drafted or registered: the authority on the type. */
  decl: ManifestColumnJson;
  included: boolean;
  /** The type carried now; '' once a demotion found no plain type for it. */
  type: string;
  /** The ref's target table by REMOTE name — from the declaration, else the inventory's relation. */
  target?: string;
  demoted?: string;
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
  private readonly _plainTypes = new Map<string, string>();
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
      if (c.type !== undefined)
        this._plainTypes.set(key, c.type);
      if (c.code !== undefined)
        this._unsupported.set(key, c);
    }
    for (const state of this._tables.values()) {
      for (const column of state.columns.values()) {
        const relation = this._relations.find((r) => r.table === state.remote && r.column === column.remote);
        column.target = column.decl.type === 'ref' ? this._byLogical(column.decl.ref)?.remote :
          relation?.status === 'ref' ? relation.targetTable : undefined;
      }
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

  /** Excluding a table demotes the refs onto it; including one promotes them back. */
  includeTable(remote: string, on: boolean): void {
    const state = this._tables.get(remote);
    if (state === undefined || state.decl === null || state.included === on)
      return;
    this._mutate(() => {
      state.included = on;
      this._follow(state, on);
    });
  }

  /** Every bindable table in, or every table out. */
  includeTables(on: boolean): void {
    this._mutate(() => {
      const bindable = [...this._tables.values()].filter((s) => s.decl !== null);
      for (const state of bindable)
        state.included = on;
      for (const state of bindable)
        this._follow(state, on);
    });
  }

  /** A key column stays: the row id encodes it; a demoted column has no type to come back with. */
  includeColumn(table: string, remote: string, on: boolean): void {
    const state = this._tables.get(table);
    const column = state?.columns.get(remote);
    if (state === undefined || column === undefined || state.key.includes(remote) || column.included === on ||
        column.demoted !== undefined)
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
      columns.set(columnRemote, {remote: columnRemote, logical: name, decl: {...c}, included: true, type: c.type});
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

  private _byLogical(logical: string | undefined): TableState | undefined {
    return [...this._tables.values()].find((t) => t.decl !== null && t.logical === logical);
  }

  /** The refs onto [target] follow its inclusion: out, they become the plain type the inventory
   * knows or drop out with a reason; in, they are refs again. */
  private _follow(target: TableState, on: boolean): void {
    for (const state of this._tables.values()) {
      for (const column of state.columns.values()) {
        if (column.target !== target.remote)
          continue;
        if (on) {
          column.type = 'ref';
          if (column.demoted !== undefined) {
            column.demoted = undefined;
            column.included = true;
          }
        } else if (column.type === 'ref') {
          const plain = this._plainType(state, column);
          column.type = plain ?? '';
          if (plain === null) {
            column.demoted = `${target.remote} is not included and the draft names no plain type for this column`;
            column.included = false;
          }
        }
      }
    }
  }

  /** The type a ref carries as a plain value: what the inventory mapped the column to, else the
   * type of the key the inventory's relation points at; null where the inventory is silent. */
  private _plainType(state: TableState, column: ColumnState): string | null {
    const mapped = this._plainTypes.get(ManifestModel._key(state.remote, column.remote));
    if (mapped !== undefined)
      return mapped;
    const relation = this._relations.find((r) => r.table === state.remote && r.column === column.remote);
    const key = relation === undefined ? undefined :
      this._tables.get(relation.targetTable)?.columns.get(relation.targetColumn);
    return key === undefined || key.decl.type === 'ref' ? null : key.decl.type;
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
    const column = this._tables.get(r.table)?.columns.get(r.column);
    const targetIncluded = target !== undefined && target.decl !== null && target.included;
    const ref = column !== undefined && column.type === 'ref' && column.target === r.targetTable;
    const canFix = r.status === 'ref' && !ref && !targetIncluded && target !== undefined && target.decl !== null;
    const reason = ref ? `ref: ${target!.logical}` :
      column?.demoted !== undefined ? 'target not included — left out' :
        r.status === 'ref' ? 'target not included — stays a plain value' :
          r.message ?? 'stays a plain value';
    return {table: r.table, column: r.column, targetTable: r.targetTable, targetColumn: r.targetColumn,
      status: r.status, ref, targetIncluded, targetLogical: target?.logical ?? null, reason, canFix};
  }

  private _columnView(state: TableState, column: ColumnState): ColumnView {
    const relation = this._relations.find((r) => r.table === state.remote && r.column === column.remote);
    const view = relation === undefined ? undefined : this._relationView(relation);
    const isKey = state.key.includes(column.remote);
    return {
      table: state.remote, remote: column.remote, logical: column.logical, type: column.type,
      dbType: this._dbTypes.get(ManifestModel._key(state.remote, column.remote)),
      included: column.included, isKey,
      required: isKey || column.decl.required === true,
      isName: column.decl.isName === true, searchable: column.decl.searchable === true,
      supported: true, reason: column.demoted, relation: view,
    };
  }

  private _unsupportedView(state: TableState, u: InventoryColumn): ColumnView {
    return {
      table: state.remote, remote: u.remote, logical: u.remote, type: u.dbType ?? '',
      dbType: u.dbType, included: false, isKey: false, required: false, isName: false, searchable: false,
      supported: false, reason: u.message, code: u.code,
    };
  }

  private _toJSON(name: string): ManifestJson {
    const external = this._draft.storage?.kind === 'external';
    const tables: Record<string, ManifestTableJson> = {};
    for (const remote of this._order) {
      const state = this._tables.get(remote)!;
      if (state.decl === null || !state.included)
        continue;
      const columns: Record<string, ManifestColumnJson> = {};
      for (const column of state.columns.values()) {
        if (column.included)
          columns[column.logical] = this._columnJSON(column, external);
      }
      const {table: _table, businessKey: _key, writable, filters, delegate, ...rest} = state.decl;
      const decl: ManifestTableJson = {...rest, columns};
      if (external && state.logical !== state.remote)
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

  /** The declaration with the current type: a ref names its target's CURRENT logical name (the
   * declared one where the target is not in the manifest). */
  private _columnJSON(column: ColumnState, external: boolean): ManifestColumnJson {
    const {type: _type, ref, column: _column, ...rest} = column.decl;
    const target = column.target === undefined ? undefined : this._tables.get(column.target);
    const json: ManifestColumnJson = column.type === 'ref' ? {type: 'ref', ref: target?.logical ?? ref, ...rest} :
      {type: column.type, ...rest};
    if (external && column.logical !== column.remote)
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

/** The schema — every included table — or one table by its remote name. */
export type AccessScope = {kind: 'schema'} | {kind: 'table', table: string};

/** A group by id, with the label it is shown under (a duplicate name carries a disambiguator). */
export interface AccessPrincipal {
  id: string;
  label: string;
}

/** One row of the access grid: who, and which of View / Edit / Delete. */
export interface AccessGrant {
  scope: AccessScope;
  group: AccessPrincipal;
  view: boolean;
  edit: boolean;
  delete: boolean;
}

/** Who may see a column: everyone who may see the row (null), or only these groups. */
export interface ColumnVisibility {
  table: string;
  column: string;
  groups: AccessPrincipal[] | null;
}

export interface AccessJson {
  grants: AccessGrant[];
  visibility: ColumnVisibility[];
}

export type AccessCapability = 'view' | 'edit' | 'delete';

/** Access is not in the manifest — declared grants are refused for a user schema — so the dialog
 * applies what this holds after Create through the table grants and the column restrictions: a
 * schema-scope row is the same grant on every included table, there is no schema-wide row
 * access. Tables and columns are keyed by REMOTE name, as the manifest model keys them; the
 * editor's `plan()` resolves them to the logical names the API takes. */
export class AccessModel {
  readonly grants: Signal<AccessGrant[]>;
  readonly visibility: Signal<ColumnVisibility[]>;

  constructor(json: Partial<AccessJson> = {}) {
    const copy = AccessModel._copy(json);
    this.grants = signal(copy.grants);
    this.visibility = signal(copy.visibility);
  }

  static sameScope(a: AccessScope, b: AccessScope): boolean {
    return a.kind === b.kind && (a.kind === 'schema' || a.table === (b as {table: string}).table);
  }

  /** A bare string names a group that is its own label. */
  static principal(group: string | AccessPrincipal): AccessPrincipal {
    return typeof group === 'string' ? {id: group, label: group} : group;
  }

  grantsOf(scope: AccessScope): ReadonlySignal<AccessGrant[]> {
    return computed(() => this.grants.value.filter((g) => AccessModel.sameScope(g.scope, scope)));
  }

  /** Replaces the scope's rows wholesale — what the access grid hands back. */
  setGrants(scope: AccessScope, rows: Omit<AccessGrant, 'scope'>[]): void {
    this.grants.value = [...this.grants.peek().filter((g) => !AccessModel.sameScope(g.scope, scope)),
      ...rows.map((r) => ({...r, scope}))];
  }

  addGroup(scope: AccessScope, group: string | AccessPrincipal): void {
    const principal = AccessModel.principal(group);
    if (this._grant(scope, principal.id) !== undefined)
      return;
    this.grants.value = [...this.grants.peek(), {scope, group: principal, view: true, edit: false, delete: false}];
  }

  removeGroup(scope: AccessScope, groupId: string): void {
    this.grants.value = this.grants.peek()
      .filter((g) => !(AccessModel.sameScope(g.scope, scope) && g.group.id === groupId));
  }

  setGrant(scope: AccessScope, groupId: string, capability: AccessCapability, on: boolean): void {
    this.grants.value = this.grants.peek().map((g) =>
      AccessModel.sameScope(g.scope, scope) && g.group.id === groupId ? {...g, [capability]: on} : g);
  }

  visibilityOf(table: string, column: string): AccessPrincipal[] | null {
    return this.visibility.peek().find((v) => v.table === table && v.column === column)?.groups ?? null;
  }

  setVisibility(table: string, column: string, groups: (string | AccessPrincipal)[] | null): void {
    const rest = this.visibility.peek().filter((v) => !(v.table === table && v.column === column));
    this.visibility.value = groups === null ? rest :
      [...rest, {table, column, groups: groups.map((g) => AccessModel.principal(g))}];
  }

  toJSON(): AccessJson {
    return AccessModel._copy({grants: this.grants.peek(), visibility: this.visibility.peek()});
  }

  private static _copy(json: Partial<AccessJson>): AccessJson {
    return {
      grants: json.grants?.map((g) => ({...g, scope: {...g.scope}, group: {...g.group}})) ?? [],
      visibility: json.visibility?.map((v) =>
        ({...v, groups: v.groups === null ? null : v.groups.map((g) => ({...g}))})) ?? [],
    };
  }

  private _grant(scope: AccessScope, groupId: string): AccessGrant | undefined {
    return this.grants.peek().find((g) => AccessModel.sameScope(g.scope, scope) && g.group.id === groupId);
  }
}
