/**
 * The `domains` facade — the ONE import a plugin needs to build UI over domain
 * tables:
 *
 * ```ts
 * import {domains} from '@datagrok-libraries/domain-ui';
 *
 * const issues = await domains.table('grit.issue');       // the async boundary
 * const saved = await issues.formDialog({values: {project_id: project.id}});
 * grok.shell.preview = domains.formView(issues.form());
 * ```
 *
 * {@link domains.table} is the ONLY await: it resolves the typed client, the
 * registry metadata and the caller's capabilities once, and every widget factory on
 * the resulting {@link DomainTable} handle is SYNCHRONOUS — data loads inside the
 * widget afterwards (see `DomainForm.ready`). Strings stay on the one-shot dialog
 * openers (`pick` / `create` / `edit`), which are async by nature and keep no state.
 *
 * @module facade
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {discardFunc, saveFunc} from './actions';
import {AppView, AppViewOptions, DomainAppView, DomainAppViewOptions, DomainEntityAppView,
  DomainEntityAppViewOptions} from './app-view';
import {DomainGrid, DomainGridOptions} from './domain-grid';
import {EntityListOptions, EntityListWidget} from './entity-list';
import {acquireDomainContext, editorsOf, IDomainTableContext} from './frame-editor';
import {DomainForm, DomainFormOptions} from './form';
import {applyDomainUiStyles} from './styles';

/** Options of {@link domains.view} / {@link domains.formView}. */
export interface DomainWidgetViewOptions extends AppViewOptions {}

/** Options of {@link domains.dialog} / `DomainTable.formDialog`. */
export interface DomainDialogOptions {
  /** Dialog caption; a sensible default comes from the widget. */
  title?: string;
  /** Caption of the affirmative button (default 'OK'). */
  okText?: string;
  width?: number;
  height?: number;
}


/**
 * A prefetched handle on one domain table: the typed client plus the registry
 * metadata and the caller's capabilities, resolved ONCE by {@link domains.table}.
 *
 * Everything built from it is synchronous, which is the whole point — a page that
 * opens three widgets over the same table pays for one round of metadata, not
 * three. Data-plane reads are passed through so a handle is also enough to read
 * back what a form just wrote.
 *
 * **Capabilities are a snapshot.** Affordances (a read-only form, a missing Save)
 * come from the capabilities as they were when the handle was acquired; a later
 * grant change — or a `grok.dapi.domains.invalidateUiCaches()` — reaches only
 * widgets built from a NEWLY acquired handle. Re-acquire to re-gate.
 */
export class DomainTable<TRow = any, TInsert = DG.DomainRowInsert<TRow>,
  TColumn extends string = string,
  TExpand extends {[key: string]: {}} = {[key: string]: {}}> implements IDomainTableContext {
  private constructor(
    /** The typed data client — every read/write member of `grok.dapi.domains.table()`. */
    public readonly client: DG.DomainTableClient<TRow, TInsert, TColumn, TExpand>,
    /** Registry {@link DG.Property} metadata of the declared columns. */
    public readonly properties: DG.Property[],
    /** Display identity, security mode, audit flag, FK-inverted child tables. */
    public readonly info: DG.DomainTableInfo,
    /** The caller's effective capabilities, snapshot at acquisition. */
    public readonly capabilities: DG.DomainTableCapabilities) {}

  /** See {@link domains.table}. */
  static async acquire<TRow = any, TInsert = DG.DomainRowInsert<TRow>,
    TColumn extends string = string,
    TExpand extends {[key: string]: {}} = {[key: string]: {}}>(
    nameOrClient: string | DG.DomainTableClient<TRow, TInsert, TColumn, TExpand>):
      Promise<DomainTable<TRow, TInsert, TColumn, TExpand>> {
    const client = typeof nameOrClient === 'string'
      ? grok.dapi.domains.table<TRow, TInsert, TColumn, TExpand>(nameOrClient) : nameOrClient;
    const context = await acquireDomainContext(client);
    return new DomainTable<TRow, TInsert, TColumn, TExpand>(
      client, context.properties, context.info, context.capabilities);
  }

  /** `'<schema>.<table>'` — the row entity type and semType. */
  get table(): string { return `${this.client.schema}.${this.client.table}`; }

  /** Domain schema name. */
  get schema(): string { return this.client.schema; }

  /** Table name within the schema. */
  get name(): string { return this.client.table; }

  /** The property form of one row — a new one by default, an existing one with
   * `{row}` or `{id}`. Synchronous; the row loads inside the widget. */
  form(options?: DomainFormOptions): DomainForm {
    return new DomainForm(this, options);
  }

  /** {@link form} in a dialog: OK saves (a failed save keeps the dialog open),
   * Cancel/Esc discards silently. Resolves to whether a row was saved. */
  formDialog(options?: DomainFormOptions & DomainDialogOptions): Promise<boolean> {
    const form = this.form(options);
    const title = options?.title ?? (form.isEditing
      ? `Edit ${this.info.singularName}` : `New ${this.info.singularName}`);
    return domains.dialog(form, Object.assign({}, options, {title: title}));
  }

  /** An editable grid over the table (or over `options.query`): batch editing
   * with markers, one-transaction save, permission gating down to per-column
   * writability. Synchronous; the rows load inside the widget
   * ({@link DomainGrid.ready}). */
  grid(options?: DomainGridOptions): DomainGrid {
    return new DomainGrid(this, options);
  }

  /** The list of the table's rows — cards / brief / grid, search, New, per-item
   * commands — as an embeddable widget. Synchronous, loads inside
   * ({@link EntityListWidget.ready}). */
  list(options?: EntityListOptions): EntityListWidget {
    return new EntityListWidget(this, options);
  }

  /** {@link list} as a PAGE: the ribbon, the `DomainQuery` ⇄ URL round trip, the
   * row pages it opens, and the unsaved-changes gate. */
  listView(options?: DomainAppViewOptions): DomainAppView {
    return new DomainAppView(this.client, options, this);
  }

  /** The page of ONE row: its form, a tab per child table, history, and the
   * actions this caller may perform. Takes the row or its id. */
  entityView(row: string | DG.DomainRow, options?: DomainEntityAppViewOptions): DomainEntityAppView {
    return typeof row === 'string'
      ? new DomainEntityAppView(this.client, row, options, this)
      : DomainEntityAppView.forRow(this.client, row, options, this);
  }

  /** THE app: the list page, the row page behind every item, deep links, the
   * gate — `grok.shell.addView((await domains.table('grit.issue')).app())`.
   * {@link listView} under a name that says what it is. */
  app(options?: DomainAppViewOptions): DomainAppView {
    return this.listView(options);
  }

  // ─────────────────────── data-plane passthroughs ─────────────────────────

  /** {@link DG.DomainTableClient.query} — the fluent builder, or a spec. */
  query(): DG.DomainQueryBuilder<TRow, TColumn, TExpand, TRow, DG.DataFrame>;
  query(spec: DG.DomainQuerySpec<TColumn, keyof TExpand & string>): Promise<TRow[]>;
  query(spec?: DG.DomainQuerySpec<TColumn, keyof TExpand & string>): any {
    return spec === undefined ? this.client.query() : this.client.query(spec);
  }

  /** {@link DG.DomainTableClient.queryDf} — the same query as a typed DataFrame. */
  queryDf(spec?: DG.DomainQuerySpec<TColumn, keyof TExpand & string>): Promise<DG.DataFrame> {
    return this.client.queryDf(spec ?? {});
  }

  /** One row by id. */
  get(id: string): Promise<TRow> { return this.client.get(id); }

  /** How many rows match [filter]. */
  count(filter?: DG.DomainFilter<TColumn>): Promise<number> { return this.client.count(filter); }
}


/** The actions a composed page performs as a WHOLE (they drive every editor on it),
 * so its ribbon shows them once instead of once per child — see
 * {@link DomainWidgetView.buildRibbon}. */
const PAGE_ACTIONS = ['Save', 'Discard'];

/**
 * A page hosting any `DG.Widget`s — the container of {@link domains.view}.
 *
 * It is an {@link AppView}, so it REUSES the shipped page machinery rather than
 * repeating it: the status bar with the pending-change count and Save/Discard, the
 * unsaved-changes gate ({@link AppView.confirmDiscard}) and the `beforeunload`
 * handler. What it adds is the composition: children's roots mounted, their editors
 * rolled up into ONE dirty state, and their `getFunctions()` merged onto the ribbon.
 */
export class DomainWidgetView extends AppView {
  private readonly _children: DG.Widget[];

  constructor(children: DG.Widget[], options?: DomainWidgetViewOptions) {
    super(options);
    applyDomainUiStyles();
    this._children = children ?? [];
    this.box = true;
    this.root.appendChild(ui.box(ui.divV(this.widgets.map((c) => c.root), 'domain-ui-page')));
    this.name = options?.name ?? this._defaultName();
    for (const child of this.widgets) {
      // A widget that loads (a form) gets its editor after construction: the status
      // bar has to pick it up when it appears, not only when the page was built.
      const ready = (child as any).ready;
      if (ready != null && typeof ready.then === 'function')
        ready.then(() => {
          this.rebindEditors();
          this.syncUrl();
        });
    }
    this.rebindEditors();
    this.buildRibbon();
    this.syncUrl();
  }

  get type(): string { return 'domain-page'; }

  /** The widgets this page hosts, in layout order. NB not `Widget.children`, which
   * stays the platform's own DOM-derived structural view. The dirty rollup, the
   * nested status and the merged functions all come from here — see
   * {@link AppView.editors} / {@link AppView.getWidgetStatus}. */
  get widgets(): DG.Widget[] { return this._children ?? []; }

  detach(): void {
    for (const child of this.widgets ?? [])
      child.detach();
    super.detach();
  }

  /**
   * The children's actions on the ribbon (§F.1.3 — the form's Save/Reset show up
   * on the page for free), deduplicated: ONE page-level Save and Discard, which
   * drive every editor of the page ({@link AppView.save} / {@link AppView.discard})
   * exactly as the status bar's buttons do, plus each child's EXTRA actions —
   * Reset, AddRow, Refresh, a custom action — bound to the child that offers them
   * and named by it when more than one child offers the same.
   *
   * A composed page would otherwise repeat the whole vocabulary per child
   * (`Save, Discard, Reset, Save, Discard, AddRow, ...`). This is presentation
   * only: {@link AppView.getFunctions} and every child's `getFunctions()` keep
   * reporting the full set, which is what the machine surface reads.
   */
  protected buildRibbon(): void {
    const offered = (name: string) =>
      this.widgets.some((c) => c.getFunctions().some((f) => f.name === name));
    const items: HTMLElement[] = [];
    if (offered('Save'))
      items.push(_actionButton(saveFunc(), this));
    if (offered('Discard'))
      items.push(_actionButton(discardFunc(), this));

    const extras: {func: DG.Func, child: DG.Widget}[] = [];
    for (const child of this.widgets)
      for (const func of child.getFunctions())
        if (!PAGE_ACTIONS.includes(func.name))
          extras.push({func: func, child: child});
    for (const extra of extras) {
      const ambiguous = extras.some((e) => e !== extra && e.func.name === extra.func.name);
      items.push(_actionButton(extra.func, extra.child,
        ambiguous ? `${extra.func.friendlyName ?? extra.func.name} (${_childLabel(extra.child)})` : undefined));
    }
    this.setRibbonPanels(items.length === 0 ? [] : [items]);
  }

  /** Page identity only (the v1 ruling): the row an edit form on the page shows.
   * Field values never enter the URL. */
  protected syncUrl(): void {
    const form = this.widgets.find((c) => c instanceof DomainForm) as DomainForm | undefined;
    const id = form?.rowId;
    this.setUrlParams(id == null ? {} : {entity: id});
  }

  private _defaultName(): string {
    const form = this.widgets.find((c) => c instanceof DomainForm) as DomainForm | undefined;
    return form == null ? 'Page' : form.table;
  }
}


/** The `domains` facade — see the module doc. */
export namespace domains {

  /**
   * THE async boundary: resolves the typed client, the registry properties, the
   * table info and the caller's capabilities for [nameOrClient] (a
   * `'<schema>.<table>'` address, or a client you already hold), and hands back the
   * handle every widget factory is synchronous on.
   */
  export function table<TRow = any, TInsert = DG.DomainRowInsert<TRow>,
    TColumn extends string = string,
    TExpand extends {[key: string]: {}} = {[key: string]: {}}>(
    nameOrClient: string | DG.DomainTableClient<TRow, TInsert, TColumn, TExpand>):
      Promise<DomainTable<TRow, TInsert, TColumn, TExpand>> {
    return DomainTable.acquire<TRow, TInsert, TColumn, TExpand>(nameOrClient);
  }

  /**
   * A page hosting [children] — one dirty state, one status bar, one unsaved-changes
   * gate, and the children's actions merged onto the ribbon:
   *
   * ```ts
   * grok.shell.addView(domains.view([projects.form({row}), issues.grid({...})]));
   * ```
   */
  export function view(children: DG.Widget[], options?: DomainWidgetViewOptions): DomainWidgetView {
    return new DomainWidgetView(children, options);
  }

  /** {@link view} of ONE widget — the routed form page of §F.1.3:
   * `grok.shell.preview = domains.formView(issues.form())`. */
  export function formView(widget: DG.Widget, options?: DomainWidgetViewOptions): DomainWidgetView {
    return view([widget], options);
  }

  /**
   * Hosts [widget] in a dialog. OK saves through the widget's editors and closes
   * only when the save LANDED (a rejected save keeps the dialog open with its
   * markers); Cancel, Esc and a programmatic close discard silently (the dialog
   * ruling) and resolve `false`.
   *
   * Exactly one answer reaches the caller on every path — the first one — so no
   * promise is ever left pending.
   */
  export function dialog(widget: DG.Widget, options?: DomainDialogOptions): Promise<boolean> {
    return new Promise<boolean>((resolve) => {
      let decided = false;
      let closed: {unsubscribe(): void} | null = null;
      const decide = (saved: boolean) => {
        closed?.unsubscribe();
        closed = null;
        if (decided)
          return;
        decided = true;
        widget.detach();
        resolve(saved);
      };
      const dlg = DG.Dialog.create({title: options?.title ?? _titleOf(widget)});
      dlg.add(widget.root);
      // A custom button, not `onOK`: the platform's OK closes the dialog before the
      // save is known, and a failed save must leave the user with their values.
      dlg.addButton(options?.okText ?? 'OK', async () => {
        if (!(await _saveWidget(widget)))
          return;
        decide(true);
        dlg.close();
      }, 0);
      dlg.onCancel(() => decide(false));
      // Covers the X, Esc and a programmatic close — and fires for the decided
      // paths too, where the guard keeps the first answer.
      closed = dlg.onClose.subscribe(() => decide(false));
      dlg.show({width: options?.width, height: options?.height});
    });
  }

  /** The platform's lookup picker for [table] — the target's Domain View in a
   * dialog; resolves to the picked row, or null. */
  export function pick(table: string): Promise<DG.DomainRow | null> {
    return DG.DomainObjectHandler.pickRow(table);
  }

  /** The platform's create dialog for [table]; resolves to whether a row was saved. */
  export function create(table: string): Promise<boolean> {
    return DG.DomainObjectHandler.createRow(table);
  }

  /** The platform's edit dialog for [row]; resolves to whether it was saved. */
  export function edit(row: DG.DomainRow): Promise<boolean> {
    return DG.DomainObjectHandler.editRow(row);
  }

  /** The standard 409 reload/overwrite dialog, naming [subject]. */
  export function conflictDialog(subject: string): Promise<'reload' | 'overwrite' | null> {
    return DG.DomainObjectHandler.showConflictDialog(subject);
  }

  /** The row's audit trail as an embeddable element. */
  export function auditPane(row: DG.DomainRow): HTMLElement {
    return DG.DomainObjectHandler.auditPane(row);
  }

  /** The read-only grants pane of a REGISTRY entity (a schema or table id). */
  export function grantsPane(entityId: string, name?: string,
    options?: {readOnly?: boolean}): HTMLElement {
    return DG.DomainObjectHandler.grantsPane(entityId, name, options);
  }

  /** Opens a `/domains` address through the platform's own navigation. */
  export function open(path: string): void {
    DG.DomainObjectHandler.openPath(path);
  }

  /** Permalink to [row]'s Entity View; null for an unsaved row. */
  export function deepLink(row: DG.DomainRow | null): string | null {
    return DG.DomainObjectHandler.deepLink(row);
  }
}

/** A ribbon button running [func] on [widget], captioned by the function unless the
 * page needs to say WHICH widget it acts on. */
function _actionButton(func: DG.Func, widget: DG.Widget, caption?: string): HTMLElement {
  return ui.button(caption ?? func.friendlyName ?? func.name,
    () => func.apply({widget: widget}), `${func.description ?? ''}`);
}

/** How a ribbon button names the child it acts on when several offer the same
 * action: the table the child edits, or its widget type. */
function _childLabel(child: DG.Widget): string {
  const table = (child as any).table ?? editorsOf(child)[0]?.table;
  return typeof table === 'string' ? table.split('.').pop()! : child.type;
}

/** Saves everything [widget] answers for: its own `save()` when it has one, its
 * editors otherwise. */
async function _saveWidget(widget: any): Promise<boolean> {
  if (typeof widget?.save === 'function')
    return await widget.save() === true;
  for (const editor of editorsOf(widget))
    if (editor.isDirty && !(await editor.save()))
      return false;
  return true;
}

function _titleOf(widget: any): string {
  return widget instanceof DomainForm
    ? (widget.isEditing ? 'Edit' : 'New') + ` ${widget.table}` : `${widget?.type ?? 'Widget'}`;
}
