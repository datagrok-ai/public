/* The manifest editor's model: the draft envelope the server answers (`POST /domains/schemas/draft`
   → `{manifest, inventory, diagnostics}`) plus the user's edits, projected back into the manifest
   to submit by `toJSON()`. Every table and column is identified by its REMOTE name — the logical
   name is what the user edits — and every declaration key the editor does not expose is carried
   through untouched (Astra 4: deep preservation). A declared type — `ref` included — is the
   authority; the inventory only explains, and is consulted when the USER changes inclusion:
   excluding a ref's target demotes the ref to the plain type the inventory knows (no such type:
   the column drops out with a reason), including a target promotes back the refs the inventory
   stamped `ref`. Remote names are vocabulary of the external storage alone. The rules the server
   checks are mirrored in `ManifestRules`; the dry run remains the authority.

   Edit mode (a `baseline`): the REGISTERED manifest is the immutable comparison baseline and the
   working copy starts equal to it; a fresh draft over the same schema supplies only candidates —
   tables and columns the baseline does not carry, unchecked — and the catalog facts that mark
   drift on registered items (missing remotely, a changed type or key, "unknown" when the catalog
   could not be read). Registered names are locked by baseline membership, declared refs stay
   authoritative (a new warehouse foreign key is a suggestion), and every edit is a field on one
   flat surface, which is what the change list, the apply patch and the three-way rebase read. */
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
  kind: 'domain' | 'external';
  connection?: string;
  schema?: string;
  catalog?: string;
  writable?: boolean;
}

export interface ManifestJson {
  name: string;
  version?: string;
  /** The schema row's creation time, beside `version` on a registered manifest: an apply sends
   * both back, so a schema deleted and re-created cannot take a stale edit. */
  incarnation?: string;
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

/** How a registered item stands against the catalog read for this edit. */
export interface DriftView {
  kind: 'missing' | 'unbindable' | 'type' | 'key' | 'unknown';
  /** Validate is blocked while the item is kept. */
  blocks: boolean;
  reason: string;
}

export interface TableView {
  remote: string;
  logical: string;
  friendlyName: string;
  included: boolean;
  /** False for a view, a keyless table, a table whose key the platform cannot hold. */
  bindable: boolean;
  /** False for a bindable table the draft holds no declaration for — it was read over other
   * tables; nothing here can include it. */
  drafted: boolean;
  reason?: string;
  code?: string;
  /** The key's remote column names. */
  key: string[];
  readOnly: boolean;
  /** In the registered manifest (edit mode): its logical name and remote mapping are locked. */
  registered: boolean;
  drift?: DriftView;
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
  /** A foreign key the registered manifest does not declare — a ref on request ({@link ManifestModel.setRef}). */
  suggested: boolean;
  /** A ref this edit made of a suggestion; it can be a plain value again. */
  promoted: boolean;
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
  /** Every foreign key the warehouse reported on the column. */
  relations: RelationView[];
  registered: boolean;
  drift?: DriftView;
}

/** One line of the Review change list; `id` is stable across re-plans, by LOGICAL names (a
 * candidate's id follows its logical name, so renaming one renames its lines), so a server
 * plan's `lost` (keyed by dropped table and column) maps onto it. */
export interface ManifestChange {
  id: string;
  text: string;
  /** The table's logical name, where the change is about a table or a column. */
  table?: string;
  column?: string;
  /** A registered table or column leaves the manifest: the server purge names what goes with it. */
  removes: boolean;
}

/** One user edit replayed over a newer baseline: the change, what it went from and to, and the
 * server's value where the same field changed there too. */
export interface RebaseOp extends ManifestChange {
  from: unknown;
  to: unknown;
  server?: unknown;
}

export interface RebaseReport {
  applied: RebaseOp[];
  conflicts: RebaseOp[];
  dropped: RebaseOp[];
}

/** The apply body's manifest half: whole descriptors of the changed tables, the dropped ones,
 * and the schema metadata that changed. */
export interface ManifestPatch {
  friendlyName?: string;
  description?: string;
  storage?: {writable: boolean};
  tables?: Record<string, ManifestTableJson>;
  dropTables?: string[];
}

export interface ManifestModelOptions {
  friendlyName?: string;
  description?: string;
  takenNames?: string[];
  /** Edit mode: the registered manifest, the immutable baseline the draft annotates and extends. */
  baseline?: ManifestJson;
}

interface ColumnState {
  remote: string;
  logical: string;
  /** The declaration as drafted or registered: the authority on the type. */
  decl: ManifestColumnJson;
  included: boolean;
  /** The type carried now; '' once a demotion found no plain type for it. */
  type: string;
  /** The ref's target table by REMOTE name — from the declaration, else (create) the inventory's relation. */
  target?: string;
  demoted?: string;
  registered: boolean;
}

interface TableState {
  remote: string;
  logical: string;
  decl: ManifestTableJson | null;
  included: boolean;
  key: string[];
  columns: Map<string, ColumnState>;
  inventory?: InventoryTable;
  registered: boolean;
  drift?: DriftView;
}

interface FieldDiff {
  key: string;
  from: unknown;
  to: unknown;
}

/** A draft whose catalog read did not answer, or answered without listing the schema. */
const CATALOG_UNKNOWN_CODES = ['external-unreachable', 'external-schema-unlisted', 'external-schema-missing'];
const UNKNOWN_DRIFT: DriftView = {kind: 'unknown', blocks: false,
  reason: 'the warehouse catalog could not be read — shown as registered'};
const SEP = '\u0000';

/** The flat field surface both models expose for the change list and the rebase: JSON scalars
 * under string keys, compared by content. */
class Fields {
  /** JSON with the keys of every object sorted: equal content reads equal whatever the key order. */
  static stable(x: unknown): string {
    return JSON.stringify(x, (_k, v: unknown) =>
      v !== null && typeof v === 'object' && !Array.isArray(v) ?
        Object.fromEntries(Object.entries(v as Record<string, unknown>).sort(([a], [b]) => a < b ? -1 : a > b ? 1 : 0)) : v);
  }

  static same(a: unknown, b: unknown): boolean {
    return Fields.stable(a) === Fields.stable(b);
  }

  /** The fields whose value differs, in [b]'s order; a field only one side holds is no change. */
  static diff(a: Map<string, unknown>, b: Map<string, unknown>): FieldDiff[] {
    const out: FieldDiff[] = [];
    for (const [key, to] of b) {
      if (a.has(key) && !Fields.same(a.get(key), to))
        out.push({key, from: a.get(key), to});
    }
    return out;
  }

  /** Replays [ops], made against [before], over a state reloaded to [after]: an op on a field
   * that no longer `exists` is dropped; one whose field the server changed to another value is
   * a conflict; the rest go through `apply` and count as applied once the surface agrees.
   * Inclusion goes first — a demotion or a promotion follows its target's inclusion. */
  static replay(ops: FieldDiff[], before: Map<string, unknown>, after: Map<string, unknown>, io: {
    exists: (key: string) => boolean,
    /** The value a field reads as where the reloaded state holds no row for it. */
    fallback: (key: string) => unknown,
    apply: (d: FieldDiff) => void,
    surface: () => Map<string, unknown>,
    describe: (d: FieldDiff) => ManifestChange,
  }): RebaseReport {
    const report: RebaseReport = {applied: [], conflicts: [], dropped: []};
    const inclusion = (d: FieldDiff): boolean => d.key.endsWith('].included');
    const now = (key: string): unknown => {
      const s = io.surface();
      return s.has(key) ? s.get(key) : io.fallback(key);
    };
    for (const d of [...ops.filter(inclusion), ...ops.filter((d) => !inclusion(d))]) {
      const op = {...io.describe(d), from: d.from, to: d.to};
      if (!io.exists(d.key)) {
        report.dropped.push(op);
        continue;
      }
      const server = after.has(d.key) ? after.get(d.key) : io.fallback(d.key);
      if (!Fields.same(server, before.get(d.key)) && !Fields.same(server, d.to)) {
        report.conflicts.push({...op, server});
        continue;
      }
      if (!Fields.same(now(d.key), d.to))
        io.apply(d);
      (Fields.same(now(d.key), d.to) ? report.applied : report.dropped).push(op);
    }
    return report;
  }
}

export class ManifestModel {
  /** The manifest's `name` — what the schema is registered as. */
  readonly name: Signal<string>;
  /** Schema metadata the create envelope carries beside the manifest (Astra 9); in edit mode the
   * DomainSchema entity's, changed through the apply. */
  readonly friendlyName: Signal<string>;
  readonly description: Signal<string>;
  readonly tables: ReadonlySignal<TableView[]>;
  readonly relations: ReadonlySignal<RelationView[]>;
  readonly writable: ReadonlySignal<boolean>;
  readonly manifest: ReadonlySignal<ManifestJson>;
  /** Registered items kept while the catalog says they cannot bind — what gates Validate in edit
   * mode, addressed by manifest path like the dry run's findings. */
  readonly blockers: ReadonlySignal<ManifestDiagnostic[]>;
  /** Bumped by every edit — what a view re-renders on without rebuilding the manifest. */
  readonly revision: ReadonlySignal<number>;

  private readonly _rev = signal(0);
  private _baseline: ManifestJson | null;
  /** The manifest whose storage and unexposed keys `toJSON` carries: the draft's, or the baseline. */
  private _source!: ManifestJson;
  private readonly _tables = new Map<string, TableState>();
  private readonly _order: string[] = [];
  private readonly _unsupported = new Map<string, InventoryColumn>();
  private readonly _dbTypes = new Map<string, string>();
  private readonly _plainTypes = new Map<string, string>();
  private readonly _listed = new Map<string, number>();
  private _relations: InventoryRelation[] = [];
  private readonly _columnViews = new Map<string, ReadonlySignal<ColumnView[]>>();
  private _writable = false;
  private _catalog: 'known' | 'unknown' = 'known';
  /** The surface as loaded — the baseline every change and rebase is read against. */
  private _loaded = new Map<string, unknown>();
  /** `name` was set by hand and no longer follows the friendly name. */
  private _detached = false;
  private readonly _taken: Set<string>;

  constructor(draft: DraftEnvelope | null, options: ManifestModelOptions = {}) {
    this._baseline = options.baseline ?? null;
    this.name = signal('');
    this.friendlyName = signal(options.friendlyName ?? '');
    this.description = signal(options.description ?? '');
    this._taken = new Set([...ManifestRules.RESERVED_SCHEMA_NAMES, ...(options.takenNames ?? [])]);
    this.revision = this._rev;
    this._load(draft);
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
    this.blockers = computed(() => {
      this._rev.value;
      return this._blockers();
    });
  }

  /** The storage the draft was made over — read-only here, the vocabulary of the editor. */
  get storage(): ManifestStorageJson | undefined {
    return this._source.storage;
  }

  /** Over a registered baseline: names locked by membership, candidates from the draft, drift. */
  get editing(): boolean {
    return this._baseline !== null;
  }

  /** Whether the catalog behind this edit was read: `unknown` marks nothing missing. */
  get catalog(): 'known' | 'unknown' {
    return this._catalog;
  }

  /** The edit tokens an apply sends back, from the registered manifest. */
  get version(): string | undefined {
    return this._baseline?.version;
  }

  get incarnation(): string | undefined {
    return this._baseline?.incarnation;
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

  /** A suggested foreign key becomes a ref (its target included), or a ref this edit made becomes
   * the plain value again; a declared ref is authoritative and stays. */
  setRef(table: string, remote: string, on: boolean): void {
    const state = this._tables.get(table);
    const column = state?.columns.get(remote);
    if (state === undefined || column === undefined)
      return;
    if (on) {
      const relation = this._relations.find((r) => r.table === table && r.column === remote && r.status === 'ref');
      const target = relation === undefined ? undefined : this._tables.get(relation.targetTable);
      if (column.type === 'ref' || target === undefined || target.decl === null || !target.included)
        return;
      this._mutate(() => {
        column.type = 'ref';
        column.target = relation!.targetTable;
      });
    } else {
      const target = column.target === undefined ? undefined : this._tables.get(column.target);
      if (column.type !== 'ref' || column.decl.type === 'ref' || target === undefined)
        return;
      const plain = this._plainType(state, column, target);
      if (plain === null)
        return;
      this._mutate(() => {
        column.type = plain;
        column.target = undefined;
      });
    }
  }

  checkSchemaName(name: string): string | null {
    return ManifestRules.checkSchemaName(name);
  }

  /** The friendly name; the identifier follows it ({@link proposeName}) until set by hand — and
   * never over a registered schema, whose identifier is locked. */
  setSchemaFriendlyName(friendlyName: string): void {
    this.friendlyName.value = friendlyName;
    if (!this._detached && !this.editing)
      this.name.value = this.proposeName(friendlyName);
  }

  setDescription(description: string): void {
    this.description.value = description;
  }

  /** The identifier by hand; an empty one lets it follow the friendly name again (a view
   * re-renders even when that is the name it already had). Locked over a registered schema. */
  setSchemaName(name: string): void {
    if (this.editing)
      return;
    this._detached = name !== '';
    const next = name === '' ? this.proposeName(this.friendlyName.peek()) : name;
    if (next === this.name.peek())
      this._rev.value = this._rev.peek() + 1;
    else
      this.name.value = next;
  }

  /** A registered name {@link proposeName} steps past from now on — found by probing the registry,
   * which is not listed. */
  markTaken(name: string): void {
    this._taken.add(name);
  }

  /** The harmonized identifier, free of the registered and the reserved names: `_2`, `_3`… on a collision. */
  proposeName(friendlyName: string): string {
    const base = ManifestRules.identifier(friendlyName);
    let candidate = base;
    for (let n = 2; base !== '' && this._taken.has(candidate); n++)
      candidate = `${base.slice(0, ManifestRules.MAX_SCHEMA_NAME_LENGTH - `_${n}`.length)}_${n}`;
    return candidate;
  }

  /** A suggested friendly name, numbered like the identifier {@link proposeName} steps it to. */
  proposeFriendlyName(label: string): string {
    const name = this.proposeName(label);
    return name === ManifestRules.identifier(label) ? label : `${label} ${name.slice(name.lastIndexOf('_') + 1)}`;
  }

  checkTableName(remote: string, logical: string): string | null {
    const state = this._tables.get(remote);
    if (state?.registered === true && state.logical !== logical)
      return `${state.logical} is registered under that name — a logical name is for life`;
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
    const state = this._tables.get(table);
    const column = state?.columns.get(remote);
    if (column?.registered === true && column.logical !== logical)
      return `${state!.logical}.${column.logical} is registered under that name — a logical name is for life`;
    const problem = ManifestRules.checkColumnName(logical);
    if (problem !== null)
      return problem;
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

  /** Whether rows of the table may be written: a writable storage, the table not opted out. */
  writes(table: TableView): boolean {
    return this.writable.peek() && !table.readOnly;
  }

  /** Opts one table out of writes; the opt-out stays, dormant, while the storage is read-only. */
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

  // ─────────────────────── edit mode: the change list, the patch, the rebase ───────────────────────

  /** The edits since the load, in manifest order, worded for Review. Column lines of a table
   * that is out now say nothing beyond its removal. */
  changes(): ManifestChange[] {
    const now = this._surface();
    const out: ManifestChange[] = [];
    for (const d of Fields.diff(this._loaded, now)) {
      const field = ManifestModel._field(d.key);
      if (field.table !== undefined && field.column !== undefined && now.get(`table[${field.table}].included`) !== true)
        continue;
      if (field.table !== undefined && field.column === undefined && field.name !== 'included' &&
          now.get(`table[${field.table}].included`) !== true)
        continue;
      out.push(this._describe(d));
    }
    return out;
  }

  /** The apply body's manifest half (edit mode); empty when nothing changed. A changed table
   * travels as its whole descriptor; the comparison is by content, not by key order or a
   * redundant remote name. */
  patch(): ManifestPatch {
    const baseline = this._baseline!;
    const patch: ManifestPatch = {};
    if (this.friendlyName.peek() !== (this._loaded.get('schema.friendlyName') as string))
      patch.friendlyName = this.friendlyName.peek();
    if (this.description.peek() !== (this._loaded.get('schema.description') as string))
      patch.description = this.description.peek();
    if (this._writable !== (this._loaded.get('schema.writable') as boolean))
      patch.storage = {writable: this._writable};
    const json = this.toJSON();
    const tables: Record<string, ManifestTableJson> = {};
    for (const [logical, decl] of Object.entries(json.tables)) {
      const before = baseline.tables[logical];
      if (before === undefined || !Fields.same(ManifestModel.canonicalTable(before, logical), ManifestModel.canonicalTable(decl, logical)))
        tables[logical] = decl;
    }
    if (Object.keys(tables).length > 0)
      patch.tables = tables;
    const dropped = Object.keys(baseline.tables).filter((logical) => json.tables[logical] === undefined);
    if (dropped.length > 0)
      patch.dropTables = dropped;
    return patch;
  }

  /** The three-way merge after a version conflict: the edits since the load, as field-level ops,
   * replayed over the new baseline. An op whose field the server changed as well is a conflict
   * (the server's value stands, the user's is reported); one on an item the new baseline and
   * draft no longer hold is dropped and reported. */
  rebase(baseline: ManifestJson, draft: DraftEnvelope | null,
    options: {friendlyName: string, description: string}): RebaseReport {
    const before = this._loaded;
    const ops = Fields.diff(before, this._surface());
    this._baseline = baseline;
    this.friendlyName.value = options.friendlyName;
    this.description.value = options.description;
    this._load(draft);
    const after = this._loaded;
    // an edit on a table that is a candidate again (unregistered on the server, still in the warehouse)
    // has nowhere to go unless the replayed inclusion took it in
    const exists = (key: string): boolean => {
      const f = ManifestModel._field(key);
      const state = f.table === undefined ? undefined : this._tables.get(f.table);
      return after.has(key) && (f.table === undefined || (f.column === undefined && f.name === 'included') ||
        state !== undefined && (state.registered || state.included));
    };
    const report = Fields.replay(ops, before, after, {exists, fallback: () => undefined,
      apply: (d) => this._apply(d), surface: () => this._surface(), describe: (d) => this._describe(d)});
    this._rev.value = this._rev.peek() + 1;
    return report;
  }

  private _apply(d: FieldDiff): void {
    const f = ManifestModel._field(d.key);
    const to = d.to;
    if (f.table === undefined) {
      if (f.name === 'friendlyName')
        this.setSchemaFriendlyName(to as string);
      else if (f.name === 'description')
        this.setDescription(to as string);
      else if (f.name === 'writable')
        this.setWritable(to as boolean);
      return;
    }
    if (f.column === undefined) {
      if (f.name === 'included')
        this.includeTable(f.table, to as boolean);
      else if (f.name === 'logical')
        this.renameTable(f.table, to as string);
      else if (f.name === 'friendlyName')
        this.setFriendlyName(f.table, to as string);
      else if (f.name === 'readOnly')
        this.setReadOnly(f.table, to as boolean);
      return;
    }
    const column = f.column;
    if (f.name === 'included')
      this.includeColumn(f.table, column, to as boolean);
    else if (f.name === 'logical')
      this.renameColumn(f.table, column, to as string);
    else if (f.name === 'required')
      this.setRequired(f.table, column, to as boolean);
    else if (f.name === 'isName' || f.name === 'searchable') {
      const flagged = this._tables.get(f.table)?.columns.get(column)?.decl[f.name] === true;
      if (to === true)
        this._setSingle(f.table, column, f.name);
      else if (flagged)
        this._setSingle(f.table, null, f.name);
    } else if (f.name === 'ref')
      this.setRef(f.table, column, to !== null);
  }

  /** The Review wording of one field change, by the current logical names. */
  private _describe(d: FieldDiff): ManifestChange {
    const f = ManifestModel._field(d.key);
    const onOff = (on: unknown): string => on === true ? 'on' : 'off';
    if (f.table === undefined) {
      const text = f.name === 'writable' ? `Writable: ${onOff(d.to)}` :
        f.name === 'friendlyName' ? `Friendly name: "${String(d.to)}"` : `Description: "${String(d.to)}"`;
      return {id: `schema:${f.name}`, text, removes: false};
    }
    const state = this._tables.get(f.table);
    const table = state?.logical ?? f.table;
    if (f.column === undefined) {
      const id = f.name === 'included' ? `table:${table}` : `table:${table}:${f.name}`;
      const text = f.name === 'included' ? (d.to === true ? `Table ${table} added` : `Table ${table} removed`) :
        f.name === 'logical' ? `Table ${String(d.from)} renamed to ${String(d.to)}` :
          f.name === 'friendlyName' ? `Table ${table}: friendly name "${String(d.to)}"` :
            `Table ${table}: read-only ${onOff(d.to)}`;
      return {id, text, table, removes: f.name === 'included' && d.to === false && state?.registered === true};
    }
    const cs = state?.columns.get(f.column);
    const column = cs?.logical ?? f.column;
    const at = `${table}.${column}`;
    const id = f.name === 'included' ? `column:${at}` : `column:${at}:${f.name}`;
    let text: string;
    if (f.name === 'included')
      text = d.to === true ? `Column ${at} added` : `Column ${at} removed`;
    else if (f.name === 'logical')
      text = `Column ${table}.${String(d.from)} renamed to ${String(d.to)}`;
    else if (f.name === 'ref') {
      const target = (remote: unknown): string => this._tables.get(String(remote))?.logical ?? String(remote);
      text = d.to === null ? `Column ${at}: ref to ${target(d.from)} demoted to a plain ${cs?.type ?? 'value'}` :
        `Column ${at}: now a ref to ${target(d.to)}`;
    } else {
      const flag = f.name === 'isName' ? 'name column' : f.name;
      text = `Column ${at}: ${flag} ${onOff(d.to)}`;
    }
    return {id, text, table, column, removes: f.name === 'included' && d.to === false && cs?.registered === true};
  }

  /** Every editable field as one flat map — what the change list, the dirty check and the rebase
   * compare. Remote names key the tables and columns; values are JSON scalars. */
  private _surface(): Map<string, unknown> {
    const s = new Map<string, unknown>();
    s.set('schema.friendlyName', this.friendlyName.peek());
    s.set('schema.description', this.description.peek());
    s.set('schema.writable', this._writable);
    for (const remote of this._order) {
      const state = this._tables.get(remote)!;
      if (state.decl === null)
        continue;
      const t = `table[${remote}]`;
      s.set(`${t}.included`, state.included);
      s.set(`${t}.logical`, state.logical);
      s.set(`${t}.friendlyName`, state.decl.friendlyName ?? '');
      s.set(`${t}.readOnly`, state.decl.writable === false);
      for (const column of state.columns.values()) {
        const c = `column[${remote}${SEP}${column.remote}]`;
        s.set(`${c}.included`, column.included);
        s.set(`${c}.logical`, column.logical);
        s.set(`${c}.required`, column.decl.required === true);
        s.set(`${c}.isName`, column.decl.isName === true);
        s.set(`${c}.searchable`, column.decl.searchable === true);
        s.set(`${c}.ref`, column.type === 'ref' ? column.target ?? null : null);
      }
    }
    return s;
  }

  private static _field(key: string): {table?: string, column?: string, name: string} {
    const dot = key.lastIndexOf('.');
    const name = key.slice(dot + 1);
    const item = key.slice(0, dot);
    if (item.startsWith('table['))
      return {table: item.slice(6, -1), name};
    if (item.startsWith('column[')) {
      const [table, column] = item.slice(7, -1).split(SEP);
      return {table, column, name};
    }
    return {name};
  }

  /** A descriptor without the remote names that repeat the logical ones — the registry emits
   * them either way; the editor only where they differ. */
  /** A table descriptor by content: a remote name equal to its logical one is dropped, so a
   * registered descriptor and a sent one compare alike. */
  static canonicalTable(decl: ManifestTableJson, logical: string): ManifestTableJson {
    const {table, columns, ...rest} = decl;
    const out: ManifestTableJson = {...rest, columns: {}};
    if (table !== undefined && table !== logical)
      out.table = table;
    for (const [name, c] of Object.entries(columns)) {
      const {column, ...cRest} = c;
      out.columns[name] = column !== undefined && column !== name ? {...cRest, column} : cRest;
    }
    return out;
  }

  // ─────────────────────── loading ───────────────────────

  private _load(draft: DraftEnvelope | null): void {
    const baseline = this._baseline;
    const source = baseline ?? draft!.manifest;
    this._source = source;
    this.name.value = source.name;
    this._writable = source.storage?.writable === true;
    this._tables.clear();
    this._order.length = 0;
    this._unsupported.clear();
    this._dbTypes.clear();
    this._plainTypes.clear();
    this._listed.clear();
    this._relations = draft?.inventory?.relations ?? [];
    this._catalog = draft !== null && !(draft.diagnostics ?? []).some((d) => CATALOG_UNKNOWN_CODES.includes(d.code)) ?
      'known' : 'unknown';
    for (const [logical, decl] of Object.entries(source.tables))
      this._addTable(decl.table ?? logical, logical, decl, source, {registered: baseline !== null, included: true});
    if (baseline !== null && draft !== null)
      this._addCandidates(draft.manifest);
    for (const t of draft?.inventory?.tables ?? []) {
      const state = this._tables.get(t.remote);
      if (state === undefined)
        this._addTable(t.remote, t.logical ?? t.remote, null, source, {registered: false, included: false, inventory: t});
      else
        state.inventory = t;
    }
    for (const c of draft?.inventory?.columns ?? []) {
      const key = ManifestModel._key(c.table, c.remote);
      this._listed.set(c.table, (this._listed.get(c.table) ?? 0) + 1);
      if (c.dbType !== undefined)
        this._dbTypes.set(key, c.dbType);
      if (c.type !== undefined)
        this._plainTypes.set(key, c.type);
      if (c.code !== undefined)
        this._unsupported.set(key, c);
    }
    // create: a foreign key the draft qualified is the target a plain column follows in and out
    if (baseline === null) {
      for (const state of this._tables.values()) {
        for (const column of state.columns.values()) {
          if (column.target !== undefined)
            continue;
          const relation = this._relations.find((r) => r.table === state.remote && r.column === column.remote);
          if (relation?.status === 'ref')
            column.target = relation.targetTable;
        }
      }
    } else {
      // a read that lists none of the registered tables is a grants or catalog problem, not a warehouse
      // that lost every one of them
      const registered = [...this._tables.values()].filter((t) => t.registered);
      if (registered.length > 0 && registered.every((t) => t.inventory === undefined))
        this._catalog = 'unknown';
      this._markTableDrift();
    }
    this._loaded = this._surface();
  }

  /** The draft's tables and columns the baseline does not carry, as unchecked candidates under
   * names free of the registered ones. */
  private _addCandidates(draft: ManifestJson): void {
    for (const [logical, decl] of Object.entries(draft.tables)) {
      const remote = decl.table ?? logical;
      const state = this._tables.get(remote);
      if (state === undefined) {
        const free = ManifestModel._free(logical, [...this._tables.values()].filter((t) => t.decl !== null).map((t) => t.logical));
        this._addTable(remote, free, decl, draft, {registered: false, included: false});
        continue;
      }
      for (const [name, c] of Object.entries(decl.columns)) {
        const columnRemote = c.column ?? name;
        if (state.columns.has(columnRemote))
          continue;
        const free = ManifestModel._free(name, [...state.columns.values()].map((x) => x.logical));
        // the registered table already has its name and searchable columns
        const {isName: _isName, searchable: _searchable, ...rest} = c;
        state.columns.set(columnRemote, {remote: columnRemote, logical: free, decl: rest, included: false,
          type: c.type, target: ManifestModel._refTarget(draft, c), registered: false});
      }
    }
  }

  private static _free(name: string, taken: string[]): string {
    let candidate = name;
    for (let n = 2; taken.includes(candidate); n++)
      candidate = `${name}_${n}`;
    return candidate;
  }

  /** The remote table a declared ref points at, by the naming of the manifest it came from. */
  private static _refTarget(source: ManifestJson, c: ManifestColumnJson): string | undefined {
    if (c.type !== 'ref' || c.ref === undefined)
      return undefined;
    return source.tables[c.ref]?.table ?? c.ref;
  }

  private _addTable(remote: string, logical: string, decl: ManifestTableJson | null, source: ManifestJson,
    options: {registered: boolean, included: boolean, inventory?: InventoryTable}): void {
    const columns = new Map<string, ColumnState>();
    const keyByLogical = new Map<string, string>();
    for (const [name, c] of Object.entries(decl?.columns ?? {})) {
      const columnRemote = c.column ?? name;
      keyByLogical.set(name, columnRemote);
      columns.set(columnRemote, {remote: columnRemote, logical: name, decl: {...c}, included: true, type: c.type,
        target: ManifestModel._refTarget(source, c), registered: options.registered});
    }
    const key = decl?.businessKey?.map((k) => keyByLogical.get(k) ?? k) ?? options.inventory?.key ?? [];
    this._tables.set(remote, {remote, logical, decl: decl === null ? null : {...decl},
      included: decl !== null && options.included, key, columns, inventory: options.inventory,
      registered: options.registered});
    this._order.push(remote);
  }

  /** What the catalog says about every registered table: nothing when it could not be read. */
  private _markTableDrift(): void {
    for (const state of this._tables.values()) {
      if (!state.registered)
        continue;
      if (this._catalog === 'unknown') {
        state.drift = UNKNOWN_DRIFT;
        continue;
      }
      const inv = state.inventory;
      if (inv === undefined || inv.code === 'external-table-missing') {
        state.drift = {kind: 'missing', blocks: true, reason: 'registered; missing remotely — the warehouse no longer has this table'};
        continue;
      }
      if (!inv.bindable) {
        state.drift = {kind: 'unbindable', blocks: true, reason: inv.message ?? 'registered; the warehouse table can no longer be bound'};
        continue;
      }
      if (inv.key !== undefined && [...inv.key].sort().join(',') !== [...state.key].sort().join(',')) {
        state.drift = {kind: 'key', blocks: false,
          reason: `the warehouse key is ${inv.key.join(', ')}; registered with ${state.key.join(', ')} — the row ids depend ` +
            'on it, and a structural change makes the dry run report external-key-mismatch'};
      }
    }
  }

  /** What the catalog says about a registered column, read live: float over int is compatible
   * under a read-only binding alone, so the verdict follows the Writable switch. */
  private _columnDrift(state: TableState, column: ColumnState): DriftView | undefined {
    if (!column.registered)
      return undefined;
    if (this._catalog === 'unknown')
      return UNKNOWN_DRIFT;
    const inv = state.inventory;
    if (inv === undefined || !inv.bindable || (this._listed.get(state.remote) ?? 0) === 0)
      return undefined;
    const key = ManifestModel._key(state.remote, column.remote);
    const dbType = this._dbTypes.get(key);
    const unsupported = this._unsupported.get(key);
    const mapped = this._plainTypes.get(key);
    if (dbType === undefined && mapped === undefined && unsupported === undefined)
      return {kind: 'missing', blocks: true, reason: 'registered; missing remotely — the warehouse table no longer has this column'};
    if (unsupported?.code === 'external-column-type')
      return {kind: 'type', blocks: true, reason: unsupported.message ?? `${dbType ?? ''} is not supported`};
    if (mapped === undefined)
      return undefined;
    const target = column.target === undefined ? undefined : this._tables.get(column.target);
    const declared = column.decl.type !== 'ref' ? column.decl.type :
      target === undefined ? undefined : this._keyType(state, column, target);
    if (declared === undefined || declared === mapped)
      return undefined;
    const compatible = declared === 'float' && mapped === 'int' && !this._writable;
    return {kind: 'type', blocks: !compatible,
      reason: `declared ${declared}; the warehouse column is ${dbType ?? mapped} (${mapped})` +
        (compatible ? ' — read as float' : ' — the types do not agree')};
  }

  private _blockers(): ManifestDiagnostic[] {
    const out: ManifestDiagnostic[] = [];
    for (const remote of this._order) {
      const state = this._tables.get(remote)!;
      if (!state.included || !state.registered)
        continue;
      if (state.drift?.blocks === true) {
        out.push({path: `tables.${state.logical}`, message: state.drift.reason,
          code: state.drift.kind === 'missing' ? 'external-table-missing' : state.inventory?.code ?? 'external-table-drift'});
        continue;
      }
      for (const column of state.columns.values()) {
        const drift = this._columnDrift(state, column);
        if (column.included && drift?.blocks === true) {
          out.push({path: `tables.${state.logical}.columns.${column.logical}`, message: drift.reason,
            code: drift.kind === 'missing' ? 'external-column-missing' : 'external-column-type'});
        }
      }
    }
    return out;
  }

  private _mutate(edit: () => void): void {
    edit();
    this._rev.value = this._rev.peek() + 1;
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
          const plain = this._plainType(state, column, target);
          if (plain !== null)
            column.type = plain;
          // a key column has no type to drop out with: the ref stands, and the dry run names it
          else if (!state.key.includes(column.remote)) {
            column.type = '';
            column.demoted = `${target.remote} is not included and the draft names no plain type for this column`;
            column.included = false;
          }
        }
      }
    }
  }

  /** The type a ref carries as a plain value: what the inventory mapped the column to, else the
   * type of the key column it points at; null where neither says. */
  private _plainType(state: TableState, column: ColumnState, target: TableState): string | null {
    const mapped = this._plainTypes.get(ManifestModel._key(state.remote, column.remote));
    return mapped ?? this._keyType(state, column, target) ?? null;
  }

  /** The declared type of the key column a ref points at — the relation's, or the target's
   * one-column key; undefined where neither says or the key is itself a ref. */
  private _keyType(state: TableState, column: ColumnState, target: TableState): string | undefined {
    const relation = this._relations.find((r) => r.table === state.remote && r.column === column.remote &&
      r.targetTable === target.remote);
    const keyColumn = relation?.targetColumn ?? (target.key.length === 1 ? target.key[0] : undefined);
    const key = keyColumn === undefined ? undefined : target.columns.get(keyColumn);
    return key === undefined || key.decl.type === 'ref' ? undefined : key.decl.type;
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
    const drafted = state.decl !== null;
    const inventory = state.inventory;
    const bindable = drafted || inventory?.bindable === true;
    return {
      remote: state.remote, logical: state.logical,
      friendlyName: state.decl?.friendlyName ?? '',
      included: drafted && state.included, bindable, drafted,
      reason: !bindable ? inventory?.message ?? 'the draft holds no declaration for this table' :
        drafted ? state.drift?.reason : 'not in this draft — it was read over other tables',
      code: !bindable ? inventory?.code : drafted ? undefined : 'external-table-not-drafted',
      key: state.key, readOnly: state.decl?.writable === false,
      registered: state.registered, drift: state.drift,
    };
  }

  private _relationView(r: InventoryRelation): RelationView {
    const target = this._tables.get(r.targetTable);
    const column = this._tables.get(r.table)?.columns.get(r.column);
    const targetIncluded = target !== undefined && target.decl !== null && target.included;
    const ref = column !== undefined && column.type === 'ref' && column.target === r.targetTable;
    const canFix = r.status === 'ref' && !ref && !targetIncluded && target !== undefined && target.decl !== null &&
      column?.target === r.targetTable;
    const suggested = this.editing && r.status === 'ref' && !ref && targetIncluded && column !== undefined &&
      column.type !== 'ref' && column.target === undefined;
    const promoted = ref && column.decl.type !== 'ref' && this.editing;
    const reason = promoted ? `ref: ${target!.logical} — new in this edit` :
      ref ? `ref: ${target!.logical}` :
        column?.demoted !== undefined ? 'target not included — left out' :
          suggested ? 'a foreign key the registered manifest does not declare' :
            r.status === 'ref' ? 'target not included — stays a plain value' :
              r.message ?? 'stays a plain value';
    return {table: r.table, column: r.column, targetTable: r.targetTable, targetColumn: r.targetColumn,
      status: r.status, ref, targetIncluded, targetLogical: target?.logical ?? null, reason, canFix, suggested, promoted};
  }

  private _columnView(state: TableState, column: ColumnState): ColumnView {
    const relations = this._relations.filter((r) => r.table === state.remote && r.column === column.remote)
      .map((r) => this._relationView(r));
    const isKey = state.key.includes(column.remote);
    const drift = this._columnDrift(state, column);
    return {
      table: state.remote, remote: column.remote, logical: column.logical, type: column.type,
      dbType: this._dbTypes.get(ManifestModel._key(state.remote, column.remote)),
      included: column.included, isKey,
      required: isKey || column.decl.required === true,
      isName: column.decl.isName === true, searchable: column.decl.searchable === true,
      supported: true, reason: column.demoted ?? drift?.reason, relations,
      registered: column.registered, drift,
    };
  }

  private _unsupportedView(state: TableState, u: InventoryColumn): ColumnView {
    return {
      table: state.remote, remote: u.remote, logical: u.remote, type: u.dbType ?? '',
      dbType: u.dbType, included: false, isKey: false, required: false, isName: false, searchable: false,
      supported: false, reason: u.message, code: u.code, relations: [], registered: false,
    };
  }

  private _toJSON(name: string): ManifestJson {
    const external = this._source.storage?.kind === 'external';
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
      // a dormant opt-out is kept: the parser accepts it under a read-only storage
      if (writable === false)
        decl.writable = false;
      tables[state.logical] = decl;
    }
    const {storage, tables: _tables, name: _name, version, incarnation: _incarnation, ...rest} = this._source;
    const json: ManifestJson = {...rest, name, tables};
    if (version !== undefined)
      json.version = version;
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
    const json: ManifestColumnJson = column.type === 'ref' ? {type: 'ref', ref: target?.logical ?? ref ?? column.target, ...rest} :
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
    return `${table}${SEP}${column}`;
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
  /** Permissions the grid does not edit — Share and custom ones — held as loaded, kept as they are. */
  other?: string[];
}

/** Who may see a column: everyone who may see the row (null), or only these groups. */
export interface ColumnVisibility {
  table: string;
  column: string;
  groups: AccessPrincipal[] | null;
  /** Edit mode: who may also edit it, and the permissions beyond View and Edit, as loaded. */
  edit?: AccessPrincipal[];
  other?: {group: AccessPrincipal, permission: string}[];
  /** Restricted, but its groups could not be read: shown as such, never edited. */
  unknown?: boolean;
}

export interface AccessJson {
  grants: AccessGrant[];
  visibility: ColumnVisibility[];
}

export type AccessCapability = 'view' | 'edit' | 'delete';

/** `GET /domains/schemas/{s}/access`: per table (by logical name) its entity, the direct grants
 * as complete permission sets where the caller holds Share (null elsewhere), and per column its
 * restriction state with the View / Edit / other grants where readable. */
export interface AccessSnapshotGroup {
  id: string;
  friendlyName: string;
  personal?: boolean;
}

export interface AccessSnapshotGrant {
  group: AccessSnapshotGroup;
  permissions: string[];
}

export interface AccessSnapshotTarget {
  entityId: string;
  canShare: boolean;
  grants: AccessSnapshotGrant[] | null;
}

export interface AccessSnapshotTable extends AccessSnapshotTarget {
  remote: string;
  coreSchema: AccessSnapshotTarget;
}

export interface AccessSnapshotColumn {
  state: 'unrestricted' | 'restricted';
  canShare: boolean;
  schemaId?: string;
  view?: AccessSnapshotGroup[] | null;
  edit?: AccessSnapshotGroup[] | null;
  other?: {group: AccessSnapshotGroup, permission: string}[] | null;
}

export interface AccessSnapshot {
  tables: Record<string, AccessSnapshotTable>;
  columns: Record<string, AccessSnapshotColumn>;
  version?: string;
  incarnation?: string;
}

/** One access op of the delta: its target by REMOTE names, and the Review line — by logical
 * names, its id stable — the editor and the rebase report reuse. */
export interface AccessTriple {
  table: string;
  group: AccessPrincipal;
  permission: string;
  change: ManifestChange;
}

export interface AccessRestriction {
  table: string;
  column: string;
  grant: {group: AccessPrincipal, permission: string}[];
  revoke: {group: AccessPrincipal, permission: string}[];
  change: ManifestChange;
}

/** The exact permission-triple deltas since the load. */
export interface AccessDelta {
  grant: AccessTriple[];
  revoke: AccessTriple[];
  restrict: AccessRestriction[];
  unrestrict: {table: string, column: string, change: ManifestChange}[];
}

export interface AccessEditOptions {
  /** The logical name of a table (or of its column) by remote names — the manifest model's, so a
   * candidate resolves too; without it the registered manifest's naming serves. */
  names?: (table: string, column?: string) => string | undefined;
  /** The author's own personal group: restricting a column keeps them on it, explicitly. */
  author?: AccessPrincipal;
}

const TABLE_PERMISSIONS: [AccessCapability, string][] = [['view', 'View'], ['edit', 'Edit'], ['delete', 'Delete']];

/** Access is not in the manifest — declared grants are refused for a user schema — so the dialog
 * applies what this holds after Create through the table grants and the column restrictions: a
 * schema-scope row is the same grant on every included table, there is no schema-wide row
 * access. Tables and columns are keyed by REMOTE name, as the manifest model keys them; the
 * editor's `plan()` resolves them to the logical names the API takes.
 *
 * In edit mode ({@link AccessModel.edit}) the rows are the snapshot's direct grants per table,
 * complete permission sets kept; there are no schema-scope rows ("every table" is a bulk edit
 * over explicit tables), a target the caller may not share is read-only, and {@link delta}
 * answers the exact triples the apply carries. */
export class AccessModel {
  readonly grants: Signal<AccessGrant[]>;
  readonly visibility: Signal<ColumnVisibility[]>;

  private _loaded: AccessJson | null = null;
  private _locked: Set<string> | 'all' = new Set();
  private readonly _entities = new Map<string, string>();
  private readonly _principals = new Map<string, AccessPrincipal>();
  private readonly _logical = new Map<string, string>();
  private _names: ((table: string, column?: string) => string | undefined) | undefined;
  private _author: AccessPrincipal | undefined;

  constructor(json: Partial<AccessJson> = {}) {
    const copy = AccessModel._copy(json);
    this.grants = signal(copy.grants);
    this.visibility = signal(copy.visibility);
  }

  /** Over a registered schema's snapshot (null: it could not be read — everything stays as it is);
   * [manifest] is the registered manifest, whose naming maps the snapshot's logical names to the
   * remote ones the model keys by. */
  static edit(snapshot: AccessSnapshot | null, manifest: ManifestJson, options: AccessEditOptions = {}): AccessModel {
    const model = new AccessModel();
    model._names = options.names;
    model._author = options.author;
    model._loadSnapshot(snapshot, manifest);
    return model;
  }

  static sameScope(a: AccessScope, b: AccessScope): boolean {
    return a.kind === b.kind && (a.kind === 'schema' || a.table === (b as {table: string}).table);
  }

  /** A bare string names a group that is its own label. */
  static principal(group: string | AccessPrincipal): AccessPrincipal {
    return typeof group === 'string' ? {id: group, label: group} : group;
  }

  /** Loaded from a snapshot: no schema-scope rows, read-only targets, exact deltas. */
  get editing(): boolean {
    return this._loaded !== null;
  }

  /** The registry entity behind a table, from the snapshot. */
  entityOf(table: string): string | undefined {
    return this._entities.get(table);
  }

  /** Whether the target's grants may be edited here: always, until a snapshot says the caller
   * cannot share it — or could not be read at all. */
  canEdit(scope: AccessScope): boolean {
    return scope.kind === 'schema' || this._locked !== 'all' && !this._locked.has(scope.table);
  }

  canEditColumn(table: string, column: string): boolean {
    return this._locked !== 'all' && !this._locked.has(`columns:${table}`) && !this._locked.has(`${table}${SEP}${column}`);
  }

  grantsOf(scope: AccessScope): ReadonlySignal<AccessGrant[]> {
    return computed(() => this.grants.value.filter((g) => AccessModel.sameScope(g.scope, scope)));
  }

  /** Replaces the scope's rows wholesale — what the access grid hands back. */
  setGrants(scope: AccessScope, rows: Omit<AccessGrant, 'scope'>[]): void {
    const kept = new Map(this.grants.peek().filter((g) => AccessModel.sameScope(g.scope, scope))
      .map((g) => [g.group.id, g.other]));
    this.grants.value = [...this.grants.peek().filter((g) => !AccessModel.sameScope(g.scope, scope)),
      ...rows.map((r) => AccessModel._row(scope, {...r, other: r.other ?? kept.get(r.group.id)}))];
    for (const r of rows)
      this._principals.set(r.group.id, r.group);
  }

  addGroup(scope: AccessScope, group: string | AccessPrincipal): void {
    const principal = AccessModel.principal(group);
    this._principals.set(principal.id, principal);
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

  /** The rows held alike on every editable one of [tables] — the bulk view over "every table":
   * a group on each of them, a capability where each grants it. */
  everyTable(tables: string[]): Omit<AccessGrant, 'scope'>[] {
    const editable = tables.filter((t) => this.canEdit({kind: 'table', table: t}));
    if (editable.length === 0)
      return [];
    const rows = editable.map((t) => this.grantsOf({kind: 'table', table: t}).peek());
    return rows[0].filter((g) => rows.every((r) => r.some((x) => x.group.id === g.group.id))).map((g) => {
      const all = rows.map((r) => r.find((x) => x.group.id === g.group.id)!);
      return {group: g.group, view: all.every((x) => x.view), edit: all.every((x) => x.edit),
        delete: all.every((x) => x.delete)};
    });
  }

  /** The bulk edit: what changed against the bulk view ({@link everyTable}) written onto every
   * editable one of [tables] — a capability toggled lands on each, a group added joins each, a
   * group gone loses on each the capabilities the view showed as held alike; what the view did
   * not show or change — a table's wider grant, Share — stays. */
  setEveryTable(tables: string[], rows: Omit<AccessGrant, 'scope'>[]): void {
    const editable = tables.filter((t) => this.canEdit({kind: 'table', table: t}));
    const before = this.everyTable(tables);
    const gone = before.filter((g) => !rows.some((r) => r.group.id === g.group.id));
    let grants = this.grants.peek();
    for (const t of editable) {
      const scope: AccessScope = {kind: 'table', table: t};
      const write = (was: Omit<AccessGrant, 'scope'>, r: Omit<AccessGrant, 'scope'>): void => {
        const at = grants.findIndex((g) => AccessModel.sameScope(g.scope, scope) && g.group.id === r.group.id);
        const next = {...grants[at]};
        for (const [capability] of TABLE_PERMISSIONS) {
          if (was[capability] !== r[capability])
            next[capability] = r[capability];
        }
        const empty = !next.view && !next.edit && !next.delete && (next.other ?? []).length === 0;
        grants = empty ? grants.filter((_g, i) => i !== at) : grants.map((g, i) => i === at ? next : g);
      };
      for (const was of gone)
        write(was, {group: was.group, view: false, edit: false, delete: false});
      for (const r of rows) {
        const was = before.find((g) => g.group.id === r.group.id);
        if (was !== undefined)
          write(was, r);
        else if (!grants.some((g) => AccessModel.sameScope(g.scope, scope) && g.group.id === r.group.id))
          grants = [...grants, AccessModel._row(scope, r)];
      }
    }
    for (const r of rows)
      this._principals.set(r.group.id, r.group);
    this.grants.value = grants;
  }

  visibilityOf(table: string, column: string): AccessPrincipal[] | null {
    return this.visibility.peek().find((v) => v.table === table && v.column === column)?.groups ?? null;
  }

  /** The column's visibility entry, unknown state included. */
  columnOf(table: string, column: string): ColumnVisibility | undefined {
    return this.visibility.peek().find((v) => v.table === table && v.column === column);
  }

  /** Who may see the column; in edit mode a group newly let in may also edit it where it holds
   * Edit on the table, one the snapshot or this edit already knew keeps the Edit it had, one
   * taken out edits it no more, and the rest is kept as loaded. */
  setVisibility(table: string, column: string, groups: (string | AccessPrincipal)[] | null): void {
    const previous = this.columnOf(table, column);
    const rest = this.visibility.peek().filter((v) => !(v.table === table && v.column === column));
    if (groups === null) {
      this.visibility.value = rest;
      return;
    }
    const principals = groups.map((g) => AccessModel.principal(g));
    for (const p of principals)
      this._principals.set(p.id, p);
    const entry: ColumnVisibility = {table, column, groups: principals};
    if (this.editing) {
      const loaded = this._loaded!.visibility.find((v) => v.table === table && v.column === column);
      const known = new Map<string, boolean>();
      for (const v of [loaded, previous]) {
        for (const g of v?.groups ?? [])
          known.set(g.id, (v!.edit ?? []).some((e) => e.id === g.id));
      }
      const writes = new Set(this.grantsOf({kind: 'table', table}).peek().filter((g) => g.edit).map((g) => g.group.id));
      entry.edit = principals.filter((p) => known.get(p.id) ?? writes.has(p.id));
      const other = previous?.other ?? loaded?.other;
      if (other !== undefined)
        entry.other = other;
    }
    this.visibility.value = [...rest, entry];
  }

  toJSON(): AccessJson {
    return AccessModel._copy({grants: this.grants.peek(), visibility: this.visibility.peek()});
  }

  /** The permission triples that turn the loaded rows into the current ones — nothing on a
   * read-only target, nothing a group holds beyond what the grid edits. A column restricted for
   * the first time keeps its author on it (View, and Edit where they may edit the table) as
   * explicit triples, unless they are among its groups already. */
  delta(): AccessDelta {
    const loaded = this._loaded ?? {grants: [], visibility: []};
    const out: AccessDelta = {grant: [], revoke: [], restrict: [], unrestrict: []};
    const now = this.grants.peek();
    const tables: string[] = [];
    for (const g of [...now, ...loaded.grants]) {
      if (g.scope.kind === 'table' && !tables.includes(g.scope.table))
        tables.push(g.scope.table);
    }
    for (const table of tables) {
      const scope: AccessScope = {kind: 'table', table};
      if (!this.canEdit(scope))
        continue;
      const before = loaded.grants.filter((g) => AccessModel.sameScope(g.scope, scope));
      const after = now.filter((g) => AccessModel.sameScope(g.scope, scope));
      const groups = [...after.map((g) => g.group), ...before.map((g) => g.group)]
        .filter((g, i, all) => all.findIndex((x) => x.id === g.id) === i);
      for (const group of groups) {
        const was = before.find((g) => g.group.id === group.id);
        const is = after.find((g) => g.group.id === group.id);
        for (const [capability, permission] of TABLE_PERMISSIONS) {
          const held = was?.[capability] === true;
          const wanted = is?.[capability] === true;
          if (wanted && !held)
            out.grant.push({table, group, permission, change: this._tableChange('grant', table, group, permission)});
          else if (held && !wanted)
            out.revoke.push({table, group, permission, change: this._tableChange('revoke', table, group, permission)});
        }
      }
    }
    const columns: ColumnVisibility[] = [...this.visibility.peek(), ...loaded.visibility];
    const seen = new Set<string>();
    for (const c of columns) {
      const key = `${c.table}${SEP}${c.column}`;
      if (seen.has(key) || !this.canEditColumn(c.table, c.column))
        continue;
      seen.add(key);
      const was = loaded.visibility.find((v) => v.table === c.table && v.column === c.column);
      const is = this.columnOf(c.table, c.column);
      if (was?.unknown === true || is?.unknown === true)
        continue;
      if (is === undefined) {
        out.unrestrict.push({table: c.table, column: c.column,
          change: this._columnChange('unrestrict', c.table, c.column, ' visible to everyone again')});
        continue;
      }
      const triples = (v: ColumnVisibility | undefined): {group: AccessPrincipal, permission: string}[] =>
        [...(v?.groups ?? []).map((group) => ({group, permission: 'View'})),
          ...(v?.edit ?? []).map((group) => ({group, permission: 'Edit'}))];
      const held = triples(was);
      const wanted = triples(is);
      const author = this._author;
      const kept = was === undefined && author !== undefined && !is.groups!.some((g) => g.id === author.id);
      if (kept) {
        wanted.push({group: author, permission: 'View'});
        if (this._grant({kind: 'table', table: c.table}, author.id)?.edit === true)
          wanted.push({group: author, permission: 'Edit'});
      }
      const same = (a: {group: AccessPrincipal, permission: string}, b: {group: AccessPrincipal, permission: string}): boolean =>
        a.group.id === b.group.id && a.permission === b.permission;
      const grant = wanted.filter((t) => !held.some((h) => same(h, t)));
      const revoke = held.filter((h) => !wanted.some((t) => same(h, t)));
      if (was !== undefined && grant.length === 0 && revoke.length === 0)
        continue;
      const who = (permission: string): string => {
        const labels = wanted.filter((t) => t.permission === permission && !(kept && t.group.id === author!.id))
          .map((t) => t.group.label);
        const you = kept && wanted.some((t) => t.permission === permission && t.group.id === author?.id);
        return labels.length === 0 ? (you ? 'you alone' : 'nobody else') : `${labels.join(', ')}${you ? ' and you' : ''}`;
      };
      const tail = was === undefined ?
        ` restricted — visible to ${who('View')}${wanted.some((t) => t.permission === 'Edit') ? `; Edit for ${who('Edit')}` : ''}` :
        `: ${[...grant.map((t) => `${t.permission} for ${t.group.label}`),
          ...revoke.map((t) => `${t.permission} revoked from ${t.group.label}`)].join(', ')}`;
      out.restrict.push({table: c.table, column: c.column, grant, revoke,
        change: this._columnChange('restrict', c.table, c.column, tail)});
    }
    return out;
  }

  /** The Review line of a table op, its id by the LOGICAL table name. */
  private _tableChange(kind: 'grant' | 'revoke', table: string, group: AccessPrincipal, permission: string): ManifestChange {
    const t = this._name(table);
    return {id: `access:${kind}:${t}:${group.id}:${permission}`, table: t, removes: false,
      text: kind === 'grant' ? `${t}: ${permission} for ${group.label}` : `${t}: ${permission} revoked from ${group.label}`};
  }

  private _columnChange(kind: 'restrict' | 'unrestrict', table: string, column: string, tail: string): ManifestChange {
    const t = this._name(table);
    const c = this._name(table, column);
    return {id: `access:${kind}:${t}.${c}`, table: t, column: c, removes: false, text: `${t}.${c}${tail}`};
  }

  private _name(table: string, column?: string): string {
    return this._names?.(table, column) ?? this._logical.get(column === undefined ? table : `${table}${SEP}${column}`) ??
      column ?? table;
  }

  /** The three-way merge of the access rows: the edits since the load replayed over a fresh
   * snapshot ({@link ManifestModel.rebase}); [known] answers whether a table (and a column of it)
   * is still in the manifest, so an edit on a table the snapshot never held — one this edit adds
   * — is replayed rather than dropped. */
  rebase(snapshot: AccessSnapshot | null, manifest: ManifestJson,
    known: (table: string, column?: string) => boolean): RebaseReport {
    const before = this._loadedSurface();
    const ops = Fields.diff(before, this._surface());
    this._loadSnapshot(snapshot, manifest);
    const after = this._loadedSurface();
    const stillThere = (key: string): boolean => {
      const f = AccessModel._field(key);
      return after.has(key) || known(f.table, f.column);
    };
    return Fields.replay(ops, before, after, {exists: stillThere, fallback: AccessModel._default,
      apply: (d) => this._apply(d), surface: () => this._surface(), describe: (d) => this._describe(d)});
  }

  /** What a field reads as with no row behind it: no grant, visible to everyone. */
  private static _default(key: string): unknown {
    return key.startsWith('grant[') ? false : null;
  }

  private _apply(d: FieldDiff): void {
    const f = AccessModel._field(d.key);
    if (f.column !== undefined) {
      const ids = d.to as string[] | null;
      if (this.canEditColumn(f.table, f.column))
        this.setVisibility(f.table, f.column, ids === null ? null : ids.map((id) => this._principal(id)));
      return;
    }
    const scope: AccessScope = {kind: 'table', table: f.table};
    if (!this.canEdit(scope))
      return;
    const group = f.group!;
    if (d.to === true && this._grant(scope, group) === undefined)
      this.addGroup(scope, this._principal(group));
    this.setGrant(scope, group, f.name as AccessCapability, d.to === true);
    const row = this._grant(scope, group);
    if (row !== undefined && !row.view && !row.edit && !row.delete && (row.other ?? []).length === 0)
      this.removeGroup(scope, group);
  }

  /** `grant[<table>][<group>].<capability>` and `visibility[<table>\0<column>].view` (sorted ids,
   * or null for everyone) — the fields the rebase compares. The current rows, plus every loaded
   * row that is gone at its default; {@link _loadedSurface} is the mirror image, so the two
   * share one key set and a removed row is a change like any other. */
  private _surface(): Map<string, unknown> {
    const loaded = this._loaded ?? {grants: [], visibility: []};
    return AccessModel._fill(AccessModel._rows(this.grants.peek(), this.visibility.peek()),
      AccessModel._rows(loaded.grants, loaded.visibility));
  }

  private _loadedSurface(): Map<string, unknown> {
    const loaded = this._loaded ?? {grants: [], visibility: []};
    return AccessModel._fill(AccessModel._rows(loaded.grants, loaded.visibility),
      AccessModel._rows(this.grants.peek(), this.visibility.peek()));
  }

  private static _rows(grants: AccessGrant[], visibility: ColumnVisibility[]): Map<string, unknown> {
    const s = new Map<string, unknown>();
    for (const g of grants) {
      if (g.scope.kind !== 'table')
        continue;
      for (const [capability] of TABLE_PERMISSIONS)
        s.set(`grant[${g.scope.table}][${g.group.id}].${capability}`, g[capability]);
    }
    for (const v of visibility)
      s.set(`visibility[${v.table}${SEP}${v.column}].view`, v.unknown === true ? 'unknown' : v.groups!.map((g) => g.id).sort());
    return s;
  }

  private static _fill(s: Map<string, unknown>, keys: Map<string, unknown>): Map<string, unknown> {
    for (const key of keys.keys()) {
      if (!s.has(key))
        s.set(key, AccessModel._default(key));
    }
    return s;
  }

  private static _field(key: string): {table: string, group?: string, column?: string, name: string} {
    const dot = key.lastIndexOf('.');
    const name = key.slice(dot + 1);
    const item = key.slice(0, dot);
    if (item.startsWith('grant[')) {
      const [table, group] = item.slice(6, -1).split('][');
      return {table, group, name};
    }
    const [table, column] = item.slice(11, -1).split(SEP);
    return {table, column, name};
  }

  private _describe(d: FieldDiff): ManifestChange {
    const f = AccessModel._field(d.key);
    if (f.column !== undefined) {
      const to = d.to as string[] | null;
      return to === null ? this._columnChange('unrestrict', f.table, f.column, ' visible to everyone again') :
        this._columnChange('restrict', f.table, f.column,
          ` restricted — visible to ${to.length === 0 ? 'nobody else' : to.map((id) => this._principal(id).label).join(', ')}`);
    }
    const permission = TABLE_PERMISSIONS.find(([c]) => c === f.name)![1];
    return this._tableChange(d.to === true ? 'grant' : 'revoke', f.table, this._principal(f.group!), permission);
  }

  private _principal(id: string): AccessPrincipal {
    return this._principals.get(id) ?? {id, label: id};
  }

  private _loadSnapshot(snapshot: AccessSnapshot | null, manifest: ManifestJson): void {
    const grants: AccessGrant[] = [];
    const visibility: ColumnVisibility[] = [];
    this._entities.clear();
    this._logical.clear();
    for (const [logical, t] of Object.entries(manifest.tables)) {
      const remote = t.table ?? logical;
      this._logical.set(remote, logical);
      for (const [name, c] of Object.entries(t.columns))
        this._logical.set(`${remote}${SEP}${c.column ?? name}`, name);
    }
    this._locked = snapshot === null ? 'all' : new Set();
    const principal = (g: AccessSnapshotGroup): AccessPrincipal => {
      const p = {id: g.id, label: g.friendlyName};
      this._principals.set(p.id, p);
      return p;
    };
    if (snapshot !== null) {
      const locked = this._locked as Set<string>;
      const remoteOf = (logical: string): string =>
        snapshot.tables[logical]?.remote ?? manifest.tables[logical]?.table ?? logical;
      for (const [logical, t] of Object.entries(snapshot.tables)) {
        const remote = remoteOf(logical);
        this._entities.set(remote, t.entityId);
        if (!t.canShare || t.grants === null)
          locked.add(remote);
        if (!t.coreSchema.canShare)
          locked.add(`columns:${remote}`);
        for (const g of t.grants ?? []) {
          const other = g.permissions.filter((p) => !TABLE_PERMISSIONS.some(([, name]) => name === p));
          grants.push(AccessModel._row({kind: 'table', table: remote}, {group: principal(g.group),
            view: g.permissions.includes('View'), edit: g.permissions.includes('Edit'),
            delete: g.permissions.includes('Delete'), other: other.length === 0 ? undefined : other}));
        }
      }
      for (const [key, c] of Object.entries(snapshot.columns)) {
        const dot = key.indexOf('.');
        const tl = key.slice(0, dot);
        const cl = key.slice(dot + 1);
        const table = remoteOf(tl);
        const column = manifest.tables[tl]?.columns[cl]?.column ?? cl;
        if (!c.canShare)
          locked.add(`${table}${SEP}${column}`);
        if (c.state !== 'restricted')
          continue;
        const entry: ColumnVisibility = {table, column, groups: (c.view ?? []).map(principal),
          edit: (c.edit ?? []).map(principal)};
        if (c.other !== undefined && c.other !== null && c.other.length > 0)
          entry.other = c.other.map((o) => ({group: principal(o.group), permission: o.permission}));
        if (c.view === null || c.view === undefined)
          entry.unknown = true;
        visibility.push(entry);
      }
    }
    this.grants.value = grants;
    this.visibility.value = visibility;
    this._loaded = AccessModel._copy({grants, visibility});
  }

  private static _row(scope: AccessScope, r: Omit<AccessGrant, 'scope'>): AccessGrant {
    const row: AccessGrant = {scope, group: r.group, view: r.view, edit: r.edit, delete: r.delete};
    if (r.other !== undefined)
      row.other = r.other;
    return row;
  }

  private static _copy(json: Partial<AccessJson>): AccessJson {
    return {
      grants: json.grants?.map((g) => ({...g, scope: {...g.scope}, group: {...g.group},
        ...(g.other === undefined ? {} : {other: [...g.other]})})) ?? [],
      visibility: json.visibility?.map((v) => ({...v, groups: v.groups === null ? null : v.groups.map((g) => ({...g})),
        ...(v.edit === undefined ? {} : {edit: v.edit.map((g) => ({...g}))}),
        ...(v.other === undefined ? {} : {other: v.other.map((o) => ({group: {...o.group}, permission: o.permission}))}),
      })) ?? [],
    };
  }

  private _grant(scope: AccessScope, groupId: string): AccessGrant | undefined {
    return this.grants.peek().find((g) => AccessModel.sameScope(g.scope, scope) && g.group.id === groupId);
  }
}
