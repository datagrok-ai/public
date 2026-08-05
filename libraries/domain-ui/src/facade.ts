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

import {AppView, AppViewOptions} from './app-view';
import {DomainFrameEditor, editorsOf, IEditorHost} from './frame-editor';
import {DomainForm, DomainFormOptions, IDomainTableContext} from './form';
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
export class DomainTable<TRow = any> implements IDomainTableContext {
  private constructor(
    /** The typed data client — every read/write member of `grok.dapi.domains.table()`. */
    public readonly client: DG.DomainTableClient<TRow>,
    /** Registry {@link DG.Property} metadata of the declared columns. */
    public readonly properties: DG.Property[],
    /** Display identity, security mode, audit flag, FK-inverted child tables. */
    public readonly info: DG.DomainTableInfo,
    /** The caller's effective capabilities, snapshot at acquisition. */
    public readonly capabilities: DG.DomainTableCapabilities) {}

  /** See {@link domains.table}. */
  static async acquire<TRow = any>(
    nameOrClient: string | DG.DomainTableClient<TRow>): Promise<DomainTable<TRow>> {
    const client = typeof nameOrClient === 'string'
      ? grok.dapi.domains.table<TRow>(nameOrClient) : nameOrClient;
    const address = `${client.schema}.${client.table}`;
    const [properties, info, capabilities] = await Promise.all([
      grok.dapi.domains.registry.rowProperties(address),
      grok.dapi.domains.registry.tableInfo(address),
      client.capabilities(),
    ]);
    return new DomainTable<TRow>(client, properties, info, capabilities);
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

  // ─────────────────────── data-plane passthroughs ─────────────────────────

  /** {@link DG.DomainTableClient.query} — the fluent builder, or a spec. */
  query(): DG.DomainQueryBuilder<TRow>;
  query(spec: DG.DomainQuerySpec): Promise<TRow[]>;
  query(spec?: DG.DomainQuerySpec): any {
    return spec === undefined ? this.client.query() : this.client.query(spec);
  }

  /** {@link DG.DomainTableClient.queryDf} — the same query as a typed DataFrame. */
  queryDf(spec?: DG.DomainQuerySpec): Promise<DG.DataFrame> {
    return this.client.queryDf(spec ?? {});
  }

  /** One row by id. */
  get(id: string): Promise<TRow> { return this.client.get(id); }

  /** How many rows match [filter]. */
  count(filter?: DG.DomainFilter): Promise<number> { return this.client.count(filter); }
}


/**
 * A page hosting any `DG.Widget`s — the container of {@link domains.view}.
 *
 * It is an {@link AppView}, so it REUSES the shipped page machinery rather than
 * repeating it: the status bar with the pending-change count and Save/Discard, the
 * unsaved-changes gate ({@link AppView.confirmDiscard}) and the `beforeunload`
 * handler. What it adds is the composition: children's roots mounted, their editors
 * rolled up into ONE dirty state, and their `getFunctions()` merged onto the ribbon.
 */
export class DomainWidgetView extends AppView implements IEditorHost {
  /** The widgets this page hosts, in layout order. NB not `Widget.children`, which
   * stays the platform's own DOM-derived structural view. */
  readonly widgets: DG.Widget[];

  constructor(children: DG.Widget[], options?: DomainWidgetViewOptions) {
    super(options);
    applyDomainUiStyles();
    this.widgets = children ?? [];
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

  /** Every editor on the page: the ones tracked directly plus the children's, read
   * live (a child's editor may appear after the page was built). */
  get editors(): DomainFrameEditor[] {
    const all = super.editors.slice();
    for (const child of this.widgets ?? [])
      for (const editor of editorsOf(child))
        if (!all.includes(editor))
          all.push(editor);
    return all;
  }

  /** The children's actions, merged and de-duplicated — what the ribbon offers and
   * what a caller (or the AI) invokes with `f.apply({widget: page})`. */
  getFunctions(): DG.Func[] {
    const merged: DG.Func[] = [];
    for (const child of this.widgets)
      for (const func of child.getFunctions())
        if (!merged.some((f) => f.name === func.name))
          merged.push(func);
    return merged;
  }

  /** The page's structure: every child's status, nested under its own name — so a
   * form's inputs are addressable as `DomainForm[0].title` from the page. */
  getWidgetStatus(): DG.IWidgetStatus {
    const parts: {[name: string]: Element} = {page: this.root};
    const inputs: DG.IInputStatus[] = [];
    const descriptions: string[] = [];
    let error: string | null = null;
    this.widgets.forEach((child, index) => {
      const prefix = `${child.type}[${index}]`;
      const status = child.getWidgetStatus();
      for (const name of Object.keys(status.parts ?? {}))
        parts[`${prefix}.${name}`] = status.parts[name];
      for (const input of status.inputs ?? [])
        inputs.push(Object.assign({}, input, {name: `${prefix}.${input.name}`}));
      if (status.description != null)
        descriptions.push(`${prefix}: ${status.description}`);
      if (error == null && status.error != null)
        error = status.error;
    });
    const changes = this.editors.reduce((n, e) => n + e.changeCount, 0);
    descriptions.unshift(`The "${this.name}" page hosts ${this.widgets.length} ` +
      `widget${this.widgets.length === 1 ? '' : 's'}, with ${changes} unsaved ` +
      `change${changes === 1 ? '' : 's'}.`);
    return {
      parts: parts,
      hitAreas: {},
      shortcuts: {},
      events: [],
      description: descriptions.join(' '),
      error: error,
      inputs: inputs.length > 0 ? inputs : undefined,
    };
  }

  /** Saves every pending batch on the page — what the shared `Save` Func calls. */
  save(): Promise<boolean> { return this.saveAll(); }

  /** Drops every pending batch on the page. */
  discard(): void { this.discardAll(); }

  /** {@link discard} — a page has no separate "back to the opened values". */
  reset(): void { this.discardAll(); }

  detach(): void {
    for (const child of this.widgets ?? [])
      child.detach();
    super.detach();
  }

  /** The children's actions on the ribbon, each bound to the child that offers it
   * (§F.1.3 — the form's Save/Reset show up on the page's ribbon for free). */
  protected buildRibbon(): void {
    const items: HTMLElement[] = [];
    for (const child of this.widgets)
      for (const func of child.getFunctions())
        items.push(ui.button(func.friendlyName ?? func.name,
          () => func.apply({widget: child}), `${func.description ?? ''}`));
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
  export function table<TRow = any>(
    nameOrClient: string | DG.DomainTableClient<TRow>): Promise<DomainTable<TRow>> {
    return DomainTable.acquire<TRow>(nameOrClient);
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
