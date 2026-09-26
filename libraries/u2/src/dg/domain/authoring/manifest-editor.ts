/* The reusable manifest editor: the tree and the context panel side by side over one
   `ManifestModel` and one `AccessModel`, the field offer decided by the editor context.
   Platform-free — the dialog lane (PowerPack) drafts, dry-runs, creates and applies the access;
   this control only edits and answers `plan()` (create) or `editPlan()` (edit: the apply body
   and the change list). Diagnostics the dry run returns are set on `diagnostics` and land on the
   row and the panel that own their manifest path. */
import {Control} from '../../../core/component.js';
import {signal, Signal, ReadonlySignal} from '../../../core/signals.js';
import {Splitter} from '../../../components/containers/splitter.js';
import type {EditorContext, FieldOffer} from './editor-context.js';
import {AccessModel, ManifestModel} from './manifest-model.js';
import type {AccessJson, AccessPrincipal, AccessSnapshot, DraftEnvelope, ManifestChange, ManifestDiagnostic,
  ManifestJson, ManifestPatch, ManifestSelection, RebaseReport} from './manifest-model.js';
import {ManifestTree} from './manifest-tree.js';
import {ManifestContextPanel} from './manifest-panel.js';
import type {PrincipalPicker} from './manifest-panel.js';

export interface ManifestEditorOptions {
  context: EditorContext;
  /** The groups the access pickers offer up front — by id with a label, or a bare name that is both. */
  groups?: (string | AccessPrincipal)[];
  /** The look-up that finds any other group or user (the platform's picker, in the dialog). */
  principalPicker?: PrincipalPicker;
  /** The schema's friendly name, carried beside the manifest; in edit mode the DomainSchema entity's. */
  friendlyName?: string;
  /** The DomainSchema entity's description (edit mode). */
  description?: string;
  /** Access rows to start from (an editor reopened on what it produced). */
  access?: Partial<AccessJson>;
  /** Edit mode: the registered manifest (`GET …/manifest`, with `version` and `incarnation`). */
  baseline?: ManifestJson;
  /** Edit mode: the access snapshot (`GET …/access`); null when it could not be read — access
   * is then shown as far as it goes and never changed. */
  snapshot?: AccessSnapshot | null;
  /** Edit mode: the author's own personal group — restricting a column keeps them on it explicitly. */
  author?: AccessPrincipal;
  /** The left pane's share of the width (0.42 by default). */
  treeShare?: number;
  /** Why the Writable switch cannot be turned on here (a DML right the author lacks); the switch
   * is then disabled with this hint. */
  writableDisabled?: string;
}

/** A grant to apply after Create, by the names the API takes: one per table and group, the
 * schema-scope rows fanned out over every included table. */
export interface PlannedGrant {
  /** The table's logical name. */
  table: string;
  group: AccessPrincipal;
  view: boolean;
  edit: boolean;
  delete: boolean;
}

/** A column restriction to apply after Create, by logical names. */
export interface PlannedRestriction {
  table: string;
  column: string;
  groups: AccessPrincipal[];
}

/** Everything the dialog submits: the create envelope and what to apply afterwards. Grants and
 * restrictions of tables and columns the manifest no longer carries are left out. */
export interface ManifestPlan {
  name: string;
  friendlyName: string;
  manifest: ManifestJson;
  grants: PlannedGrant[];
  restrictions: PlannedRestriction[];
}

/** The `access` section of the apply body: exact permission triples by logical names, group ids;
 * a column op names the restriction state it was made from, which the server checks under the
 * lock (`access-conflict`). */
export interface ApplyAccess {
  grant: {table: string, group: string, permission: string}[];
  revoke: {table: string, group: string, permission: string}[];
  restrict: {table: string, column: string, from: 'restricted' | 'unrestricted', grant: {group: string, permission: string}[],
    revoke: {group: string, permission: string}[]}[];
  unrestrict: {table: string, column: string, from: 'restricted'}[];
}

/** `POST /domains/schemas/{s}/apply` — the edit tokens, the manifest patch and the access deltas;
 * a no-op edit carries the tokens alone. */
export interface ApplyPayload extends ManifestPatch {
  ifVersion: string;
  ifIncarnation?: string;
  access?: ApplyAccess;
}

export interface EditPlan {
  payload: ApplyPayload;
  /** The Review list, in manifest order then access; empty for a no-op edit. */
  changes: ManifestChange[];
  /** Registered items kept while the catalog says they cannot bind: Validate waits on them. */
  blockers: ManifestDiagnostic[];
}

export class ManifestEditor extends Control {
  readonly model: ManifestModel;
  readonly access: AccessModel;
  readonly tree: ManifestTree;
  readonly panel: ManifestContextPanel;
  readonly selected: ReadonlySignal<ManifestSelection>;
  /** The dry run's findings, addressed by manifest path; set them and the rows light up. */
  readonly diagnostics: Signal<ManifestDiagnostic[]>;
  readonly context: EditorContext;
  /** Edit mode: the snapshot and the manifest were read at different schema versions — reload
   * before editing, or the access rows describe another state than the tables. */
  stale = false;

  private readonly _author: AccessPrincipal | undefined;

  constructor(draft: DraftEnvelope | null, options: ManifestEditorOptions) {
    super();
    this.context = options.context;
    const editing = options.context.mode === 'edit';
    if (editing && options.baseline === undefined)
      throw new Error('u2: the manifest editor\'s edit mode needs the registered manifest (options.baseline)');
    if (!editing && draft === null)
      throw new Error('u2: the manifest editor needs a draft');
    this._author = options.author;
    this.model = new ManifestModel(draft, {friendlyName: options.friendlyName, description: options.description,
      baseline: editing ? options.baseline : undefined});
    this.access = editing ? this._accessOf(options.snapshot ?? null, options.baseline!) : new AccessModel(options.access);
    this.diagnostics = signal<ManifestDiagnostic[]>(draft?.diagnostics ?? []);
    this.root.classList.add('u2-manifest-editor');
    this.root.dataset.u2 = 'manifest-editor';
    const editable = options.context.mode !== 'view';
    this.tree = this.runInScope(() => new ManifestTree(this.model, {editable, diagnostics: this.diagnostics}));
    this.selected = this.tree.selected;
    this.panel = this.runInScope(() => new ManifestContextPanel(this.model, {selected: this.selected,
      access: this.access, context: options.context, groups: options.groups?.map((g) => AccessModel.principal(g)),
      principalPicker: options.principalPicker, diagnostics: this.diagnostics,
      writableDisabled: options.writableDisabled}));
    const share = options.treeShare ?? 0.42;
    const splitter = this.runInScope(() => new Splitter([this.tree, this.panel],
      {direction: 'horizontal', sizes: [share, 1 - share], minSize: 200}));
    this.root.append(splitter.root);
  }

  get offer(): FieldOffer {
    return this.panel.offer;
  }

  /** Opens the path to the node and selects it — where a diagnostic's row is taken to. */
  select(selection: ManifestSelection): Promise<void> {
    return this.tree.select(selection);
  }

  /** The create envelope plus the access to apply, resolved to the manifest's logical names. A
   * schema-scope row and a table row for the same group merge into one grant per table; Edit and
   * Delete only where the binding and the table are writable — the lock the panel shows. */
  plan(): ManifestPlan {
    if (this.model.editing)
      throw new Error('u2: plan() is the create envelope — an edit answers editPlan()');
    const model = this.model;
    const manifest = model.toJSON();
    const included = model.tables.peek().filter((t) => t.included);
    const logicalOf = (remote: string): string | null => included.find((t) => t.remote === remote)?.logical ?? null;
    const grants = new Map<string, PlannedGrant>();
    for (const g of this.access.grants.peek()) {
      const scope = g.scope;
      const tables = scope.kind === 'schema' ? included : included.filter((t) => t.remote === scope.table);
      for (const t of tables) {
        const writes = model.writes(t);
        const key = `${t.logical}\u0000${g.group.id}`;
        const merged = grants.get(key) ?? {table: t.logical, group: g.group, view: false, edit: false, delete: false};
        grants.set(key, {...merged, view: merged.view || g.view, edit: merged.edit || writes && g.edit,
          delete: merged.delete || writes && g.delete});
      }
    }
    const restrictions: PlannedRestriction[] = [];
    for (const v of this.access.visibility.peek()) {
      const table = logicalOf(v.table);
      const column = model.column(v.table, v.column);
      if (v.groups !== null && table !== null && column !== undefined && column.included)
        restrictions.push({table, column: column.logical, groups: [...v.groups]});
    }
    return {name: model.name.peek(), friendlyName: model.friendlyName.peek(), manifest, grants: [...grants.values()],
      restrictions};
  }

  /** The apply body and the Review list of an edit: the tokens from the baseline, changed tables
   * as whole descriptors, the dropped ones, the metadata that changed, and the access deltas as
   * exact triples by logical names — none on a table the same body drops (the server purge owns
   * them), on a column that is out, or on a key column (it cannot be restricted). */
  editPlan(): EditPlan {
    const model = this.model;
    if (!model.editing)
      throw new Error('u2: editPlan() needs the editor\'s edit mode — a create answers plan()');
    if (model.version === undefined)
      throw new Error('u2: the registered manifest carries no version — an apply cannot name what it was edited against');
    const payload: ApplyPayload = {ifVersion: model.version, ...model.patch()};
    if (model.incarnation !== undefined)
      payload.ifIncarnation = model.incarnation;
    const changes = model.changes();
    const included = (remote: string): boolean => model.table(remote)?.included === true;
    const restrictable = (table: string, column: string): boolean => {
      const c = model.column(table, column);
      return c !== undefined && c.included && !c.isKey;
    };
    const delta = this.access.delta();
    const access: ApplyAccess = {grant: [], revoke: [], restrict: [], unrestrict: []};
    for (const kind of ['grant', 'revoke'] as const) {
      for (const op of delta[kind].filter((o) => included(o.table))) {
        access[kind].push({table: op.change.table!, group: op.group.id, permission: op.permission});
        changes.push(op.change);
      }
    }
    const triple = (t: {group: AccessPrincipal, permission: string}): {group: string, permission: string} =>
      ({group: t.group.id, permission: t.permission});
    for (const op of delta.restrict.filter((o) => included(o.table) && restrictable(o.table, o.column))) {
      access.restrict.push({table: op.change.table!, column: op.change.column!, from: op.from, grant: op.grant.map(triple),
        revoke: op.revoke.map(triple)});
      changes.push(op.change);
    }
    for (const op of delta.unrestrict.filter((o) => included(o.table) && restrictable(o.table, o.column))) {
      access.unrestrict.push({table: op.change.table!, column: op.change.column!, from: op.from});
      changes.push(op.change);
    }
    if (access.grant.length + access.revoke.length + access.restrict.length + access.unrestrict.length > 0)
      payload.access = access;
    return {payload, changes, blockers: model.blockers.peek()};
  }

  /** After a `version-conflict`: the user's edits replayed over the schema as it is now — the
   * registered manifest, its access snapshot and a fresh draft. Same-field conflicts and edits on
   * items that vanished are reported, never applied silently; the tree and the panel rebuild. */
  rebase(baseline: ManifestJson, snapshot: AccessSnapshot | null, draft: DraftEnvelope | null,
    options: {friendlyName: string, description: string}): RebaseReport {
    const manifest = this.model.rebase(baseline, draft, options);
    const access = this.access.rebase(snapshot, baseline, (table, column) => this._known(table, column));
    this.stale = ManifestEditor._stale(snapshot, baseline);
    this.diagnostics.value = draft?.diagnostics ?? [];
    return {applied: [...manifest.applied, ...access.applied], conflicts: [...manifest.conflicts, ...access.conflicts],
      dropped: [...manifest.dropped, ...access.dropped]};
  }

  private _accessOf(snapshot: AccessSnapshot | null, baseline: ManifestJson): AccessModel {
    this.stale = ManifestEditor._stale(snapshot, baseline);
    return AccessModel.edit(snapshot, baseline, {author: this._author, names: (table, column) => {
      const t = this.model.table(table);
      return t === undefined || !t.drafted ? undefined : column === undefined ? t.logical : this.model.column(table, column)?.logical;
    }});
  }

  /** An access edit still has a target: a table (and its column) the apply will carry — registered,
   * or a candidate this edit includes; a table unregistered on the server but still in the
   * warehouse is a candidate again, and an edit on it has nowhere to go. */
  private _known(table: string, column?: string): boolean {
    const t = this.model.table(table);
    if (t === undefined || !t.drafted || !(t.registered || t.included))
      return false;
    if (column === undefined)
      return true;
    const c = this.model.column(table, column);
    return c !== undefined && (c.registered || c.included);
  }

  private static _stale(snapshot: AccessSnapshot | null, baseline: ManifestJson): boolean {
    return snapshot !== null && (snapshot.version !== baseline.version || snapshot.incarnation !== baseline.incarnation);
  }
}
