/**
 * The app view hierarchy: {@link AppView} (any page of a Datagrok app),
 * {@link DomainAppView} (the list page of one domain table) and
 * {@link DomainEntityAppView} (the page of one row). Defaults all the way down —
 * a subclass that overrides NOTHING is a working browse/CRUD app:
 *
 * ```ts
 * //name: Issues
 * //meta.role: app
 * //input: string path {meta.url: true; optional: true}
 * //output: view result
 * export function issuesApp(): DG.ViewBase {
 *   return new DomainAppView(grok.dapi.domains.table('grit.issue'));
 * }
 * ```
 *
 * Routing rides the platform's own `#app` + `meta.url` convention: no new
 * registration mechanism, no router. The page publishes its state as URL
 * parameters (the `DomainQuery` ones, plus the reserved `view=` and `entity=`)
 * and restores itself from them, so a pasted link opens the exact view.
 *
 * @module app-view
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DomainFrameEditor} from './frame-editor';
import {DomainGrid} from './domain-grid';
import {EntityListWidget, EntityListMode, EntityListOptions} from './entity-list';
import {domainHandler} from './handler';
import {confirmDiscardChanges} from './unsaved';
import {applyDomainUiStyles} from './styles';

/** Reserved URL parameter carrying the list's render mode. */
export const VIEW_MODE_PARAM = 'view';
/** Reserved URL parameter carrying the id of the row an entity page shows. */
export const ENTITY_PARAM = 'entity';

/** Options common to every {@link AppView}. */
export interface AppViewOptions {
  /** View name (the tab caption); each subclass has a sensible default. */
  name?: string;
  helpUrl?: string;
}

/**
 * Base of every page of a domain app: a {@link DG.ViewBase} with a URL contract,
 * a ribbon, a status bar, and — the part apps get wrong on their own — ONE place
 * that decides what happens to unsaved changes.
 *
 * **Unsaved changes.** The editing controls ({@link DomainFrameEditor},
 * {@link DomainGrid}) are refresh-agnostic by design: rebuilding discards
 * pending edits and they never reconcile. This class is the counterpart of that
 * contract: {@link track} registers an editor, {@link confirmDiscard} prompts
 * save / discard / cancel before anything that would rebuild, and the status bar
 * shows the pending count with Save and Discard. Subclasses call
 * {@link confirmDiscard} — they never rebuild behind it.
 */
export class AppView extends DG.ViewBase {
  /** The app's FuncCall, captured at construction — what makes a page opened
   * later (an entity page) keep the app's URL prefix and breadcrumbs. */
  readonly appCall: DG.FuncCall | undefined;

  private readonly _statusPanel = ui.div([], 'domain-ui-status-bar');
  private readonly _statusText = ui.divText('');
  private readonly _statusButtons: HTMLElement;
  private readonly _tracked: DomainFrameEditor[] = [];
  private _editorSubs: {unsubscribe(): void}[] = [];

  constructor(options?: AppViewOptions) {
    super(null, '', true);
    applyDomainUiStyles();
    this.appCall = grok.functions.getCurrentCall();
    this.name = options?.name ?? 'App';
    if (options?.helpUrl != null)
      this.helpUrl = options.helpUrl;
    this._statusButtons = ui.divH([
      ui.button('Save', () => this.saveAll(), 'Save every pending change'),
      ui.button('Discard', () => this.discardAll(), 'Discard every pending change'),
    ]);
    this._statusPanel.appendChild(this._statusText);
    this._statusPanel.appendChild(this._statusButtons);
    this.statusBarPanels = [this._statusPanel as HTMLDivElement];
    this.syncStatusBar();
  }

  get type(): string { return 'domain-app'; }

  /** Every editor whose pending changes this page answers for. Subclasses widen
   * it with the ones their widgets own (and rebuild them freely — this is read
   * on demand, never cached). */
  get editors(): DomainFrameEditor[] { return this._tracked; }

  /** Whether anything on the page has unsaved changes. */
  get isDirty(): boolean { return this.editors.some((e) => e.isDirty); }

  /** Registers an editor the page owns directly (a detail grid, a form). */
  track(editor: DomainFrameEditor): void {
    if (editor != null && !this._tracked.includes(editor))
      this._tracked.push(editor);
    this.rebindEditors();
  }

  /**
   * THE gate in front of anything that would discard pending changes — a query
   * change, a mode switch, a refresh, leaving the page. Resolves to whether the
   * caller may proceed; see {@link confirmDiscardChanges} for the outcomes
   * (including the mid-save refusal: an editor is closed while its transaction
   * is in flight, so nothing is prompted and nothing proceeds).
   */
  confirmDiscard(action: string = 'continue'): Promise<boolean> {
    return confirmDiscardChanges(this.editors, {action: action, subject: this.name});
  }

  /** Saves every pending batch (one transaction each). */
  async saveAll(): Promise<boolean> {
    for (const editor of this.editors)
      if (editor.isDirty && !(await editor.save()))
        return false;
    return true;
  }

  /** Drops every pending batch. */
  discardAll(): void {
    for (const editor of this.editors)
      if (editor.isDirty)
        editor.discard();
  }

  /** Re-subscribes the status bar to the CURRENT set of editors — call it after
   * anything that rebuilds them (a list reloading its grid). */
  protected rebindEditors(): void {
    // The editors' observables never complete on detach, so every subscription
    // is released by hand, here and in detach().
    for (const sub of this._editorSubs)
      sub.unsubscribe();
    this._editorSubs = [];
    for (const editor of this.editors) {
      this._editorSubs.push(editor.onDirtyChanged.subscribe(() => this.syncStatusBar()));
      this._editorSubs.push(editor.onSavingChanged.subscribe(() => this.syncStatusBar()));
    }
    this.syncStatusBar();
  }

  /** 'N unsaved changes' + Save / Discard while anything is pending; the
   * buttons close while a save is in flight (the editors refuse both then). */
  protected syncStatusBar(): void {
    const editors = this.editors.filter((e) => e.isDirty);
    const changes = editors.reduce((n, e) => n + e.changeCount, 0);
    const saving = this.editors.some((e) => e.isSaving);
    ui.setDisplay(this._statusPanel, changes > 0 || saving);
    ui.setDisplay(this._statusButtons, !saving);
    this._statusText.textContent = saving ? 'Saving...'
      : `${changes} unsaved change${changes === 1 ? '' : 's'}`;
  }

  // ─────────────────────── routing ─────────────────────────

  /** The page's URL parameters — the ones a deep link carries. The platform
   * hands a view only the path on navigation, so the query string is read from
   * the address bar (which the router has already updated). */
  protected urlParams(): {[key: string]: string} {
    const params: {[key: string]: string} = {};
    if (typeof window === 'undefined')
      return params;
    new URLSearchParams(window.location.search).forEach((value, key) => params[key] = value);
    return params;
  }

  /** Publishes [params] as the page's URL, under whatever prefix the app has
   * (`/apps/<Package>/<App>`) — the platform composes it. */
  protected setUrlParams(params: {[key: string]: string}): void {
    const query = new URLSearchParams();
    for (const key of Object.keys(params))
      if (params[key] != null && params[key] !== '')
        query.set(key, params[key]);
    const s = query.toString();
    this.path = s === '' ? '' : `?${s}`;
  }

  /** The page's own path, prefix included and parameters stripped. */
  protected get pagePath(): string { return `${this.path ?? ''}`.split('?')[0]; }

  /** Restores the page from {@link urlParams}. Override with the page's own
   * state; the base implementation does nothing. */
  restoreFromUrl(): void { }

  /** The platform hands back / forward navigation to the view that accepts the
   * path; app URLs re-run the app function instead, and both end up in
   * {@link restoreFromUrl}. */
  handlePath(_urlPath: string): void { this.restoreFromUrl(); }

  acceptsPath(urlPath: string): boolean {
    const own = this.pagePath;
    return own !== '' && `${urlPath}`.toLowerCase().startsWith(own.toLowerCase());
  }

  detach(): void {
    for (const sub of this._editorSubs)
      sub.unsubscribe();
    this._editorSubs = [];
    super.detach();
  }
}


/** Options of {@link DomainAppView} — the list options plus the view's own. */
export interface DomainAppViewOptions extends AppViewOptions, EntityListOptions {
  /** What opening an item does (default: the table's {@link DomainEntityAppView},
   * kept inside the app). */
  onOpen?: (row: DG.DomainRow) => void;
}

/**
 * The main page of a domain app: an {@link EntityListWidget} over one table, a
 * ribbon, the `DomainQuery` ⇄ URL round trip, and the unsaved-changes policy of
 * {@link AppView}.
 *
 * ```ts
 * const view = new DomainAppView(grok.dapi.domains.table('grit.issue'));
 * grok.shell.addView(view);
 * ```
 *
 * The page loads asynchronously (capabilities, registry metadata, rows), so the
 * constructor returns immediately — which is what an `#app` function needs.
 *
 * **URL.** `DomainQuery.toUrlParams()` (`filters[0]`, `orderBy[0]`, `limit`...)
 * plus the reserved {@link VIEW_MODE_PARAM}; opening a row pushes a
 * {@link DomainEntityAppView} carrying {@link ENTITY_PARAM}, and a link with that
 * parameter opens straight into it. Nothing UI-only ever enters the DomainQuery.
 */
export class DomainAppView extends AppView {
  readonly client: DG.DomainTableClient;
  readonly handler: DG.DomainObjectHandler;
  /** The list — null until the page has loaded (and while it reloads). */
  list: EntityListWidget | null = null;

  protected readonly options: DomainAppViewOptions;
  private readonly _host = ui.box();
  private _query: DG.DomainQuery;
  private _listSub: {unsubscribe(): void} | null = null;

  constructor(client: DG.DomainTableClient, options?: DomainAppViewOptions) {
    super(options);
    this.client = client;
    this.options = options ?? {};
    this.handler = domainHandler(this.table);
    this.name = options?.name ?? this.client.table;
    this.box = true;
    this.root.appendChild(this._host);
    this._query = new DG.DomainQuery({schema: client.schema, table: client.table});
    this._host.appendChild(ui.divText('Loading...'));
    this.load();
  }

  /** `'<schema>.<table>'`. */
  get table(): string { return `${this.client.schema}.${this.client.table}`; }

  /** What the page is showing, serializable: the shape a deep link, a saved
   * filter or an "open in Table View" carries. */
  get query(): DG.DomainQuery { return this._query; }

  get editors(): DomainFrameEditor[] {
    // Read during the base constructor (before this class's fields exist), so
    // every access is optional.
    return super.editors.concat(this.list?.editors ?? []);
  }

  /** Builds (or rebuilds) the page from the URL. */
  async load(): Promise<void> {
    const params = this.urlParams();
    this._query = this._queryFromParams(params);
    const entity = params[ENTITY_PARAM];
    try {
      const list = await EntityListWidget.create(this.client, Object.assign({}, this.options, {
        query: this._spec(this._query),
        mode: this._modeFromParams(params),
        // ONE gate for the whole page: the list's own rebuilds (search, mode,
        // query) prompt for every editor the page owns, not just its grid's.
        confirmDiscard: () => this.confirmDiscard('reload the list'),
        onOpen: (row: DG.DomainRow) => this.open(row),
      }));
      this.list?.detach();
      this.list = list;
      this.name = this.options.name ?? list.info.pluralName;
      this._host.innerHTML = '';
      this._host.appendChild(list.root);
      // A reload rebuilds the grid (and its editor) — the status bar follows
      // whatever the list owns NOW.
      this._listSub?.unsubscribe();
      this._listSub = list.onRefreshed.subscribe(() => {
        this.rebindEditors();
        // A mode switch is reserved UI state: it belongs in the URL, so a deep
        // link reopens the page the way the user left it.
        this.syncUrl();
      });
      this.rebindEditors();
      this.buildRibbon();
      this.syncUrl();
      if (entity != null && entity !== '')
        this.open(entity);
    } catch (e: any) {
      this._host.innerHTML = '';
      this._host.appendChild(ui.divText(`Cannot open ${this.table}: ${e?.message ?? e}`));
      grok.shell.error(`${this.table}: ${e?.message ?? e}`);
    }
  }

  /** Ribbon of the list page. Override to add app commands — call `super` to
   * keep Refresh and Open in Table View. */
  protected buildRibbon(): void {
    this.setRibbonPanels([[
      ui.iconFA('sync-alt', () => this.refresh(), 'Reload'),
      ui.iconFA('external-link-alt', () => this.openInTableView(),
        'Open the current query in a Table View'),
    ]]);
  }

  /** Applies a new query: prompts for pending changes, reloads, and republishes
   * the URL. Resolves to whether it was applied. */
  async setQuery(query: DG.DomainQuery): Promise<boolean> {
    if (this.list == null)
      return false;
    if (!(await this.list.setQuery(this._spec(query))))
      return false;
    this._query = query;
    this.syncUrl();
    return true;
  }

  /** Reloads the list (prompting for pending changes first). */
  async refresh(): Promise<void> {
    await this.list?.refresh();
  }

  /** Opens the platform's Table View on the current query — the full result set,
   * past the list's row cap, as a recorded, refreshable DataFrame. */
  async openInTableView(): Promise<void> {
    try {
      await this.query.run();
    } catch (e: any) {
      grok.shell.error(`${e?.message ?? e}`);
    }
  }

  /** Opens a row's page inside the app (the default of {@link EntityListWidget}'s
   * open action). Pass the row or its id. */
  open(row: DG.DomainRow | string): DomainEntityAppView | null {
    if (this.options.onOpen != null && typeof row !== 'string') {
      this.options.onOpen(row);
      return null;
    }
    const view = typeof row === 'string'
      ? new DomainEntityAppView(this.client, row, {helpUrl: this.helpUrl ?? undefined})
      : DomainEntityAppView.forRow(this.client, row, {helpUrl: this.helpUrl ?? undefined});
    // Inheriting the prefix keeps the entity page inside the app's URL — a page
    // opened by a click is created outside the app function's call.
    view.parentView = this;
    if (this.appCall != null)
      view.parentCall = this.appCall;
    grok.shell.addView(view);
    return view;
  }

  /** The list's state as URL parameters — the query plus the reserved mode. */
  protected urlState(): {[key: string]: string} {
    const params = this._query.toUrlParams();
    const mode = this.list?.mode;
    if (mode != null && mode !== 'cards')
      params[VIEW_MODE_PARAM] = mode;
    return params;
  }

  /** Publishes {@link urlState}. */
  protected syncUrl(): void {
    this.setUrlParams(this.urlState());
  }

  restoreFromUrl(): void {
    this.load();
  }

  detach(): void {
    this._listSub?.unsubscribe();
    this._listSub = null;
    this.list?.detach();
    super.detach();
  }

  private _queryFromParams(params: {[key: string]: string}): DG.DomainQuery {
    try {
      return DG.DomainQuery.fromUrlParams(this.client.schema, this.client.table, params);
    } catch (e: any) {
      grok.shell.warning(`Ignoring the URL query: ${e?.message ?? e}`);
      return new DG.DomainQuery({schema: this.client.schema, table: this.client.table});
    }
  }

  private _modeFromParams(params: {[key: string]: string}): EntityListMode | undefined {
    const mode = params[VIEW_MODE_PARAM];
    return mode === 'cards' || mode === 'brief' || mode === 'grid' ? mode : this.options.mode;
  }

  /** The query as a spec for the list. A query that declares no limit leaves the
   * spec's out too, so the list applies its own row cap instead of the
   * `DomainQuery` ceiling. */
  private _spec(query: DG.DomainQuery): DG.DomainQuerySpec {
    if (query.isAggregate)
      throw new Error('an aggregate DomainQuery cannot back a list — use openInTableView()');
    const spec = query.toSpec();
    if (query.limit == null)
      delete spec.limit;
    return spec;
  }
}


/** Options of {@link DomainEntityAppView}. */
export interface DomainEntityAppViewOptions extends AppViewOptions {
  /** Make the detail (child-table) grids editable — master-detail editing, each
   * grid saving its own transaction (default true). */
  editableDetails?: boolean;
}

/**
 * The page of ONE domain row: its name, the platform's row actions, a property
 * form, a tab per child table (an editable {@link DomainGrid} scoped to the row)
 * and its history.
 *
 * ```ts
 * grok.shell.addView(new DomainEntityAppView(grok.dapi.domains.table('grit.issue'), id));
 * ```
 *
 * Everything on it is permission-gated: actions come from the row's server-truth
 * permissions, the form from the caller's writable columns, the detail grids from
 * the child tables' own capabilities.
 */
export class DomainEntityAppView extends AppView {
  readonly client: DG.DomainTableClient;
  readonly handler: DG.DomainObjectHandler;
  /** The row — null until it has loaded. */
  row: DG.DomainRow | null = null;

  protected readonly options: DomainEntityAppViewOptions;
  private readonly _host = ui.box();
  private readonly _id: string;
  private _capabilities: DG.DomainTableCapabilities | null = null;
  private _details: DomainGrid[] = [];
  /** Values as loaded — what {@link save} diffs the form's writes against. */
  private _loaded: {[key: string]: any} = {};

  constructor(client: DG.DomainTableClient, id: string, options?: DomainEntityAppViewOptions) {
    super(options);
    this.client = client;
    this.options = options ?? {};
    this._id = id;
    this.handler = domainHandler(`${client.schema}.${client.table}`);
    this.name = options?.name ?? id;
    this.box = true;
    this.root.appendChild(this._host);
    this._host.appendChild(ui.divText('Loading...'));
    this.load();
  }

  /** The page of a row already in hand: named from it right away, so the tab
   * caption is right before the page has re-read the row from the server (which
   * it does, for its current version and the columns a list did not fetch). */
  static forRow(client: DG.DomainTableClient, row: DG.DomainRow,
    options?: DomainEntityAppViewOptions): DomainEntityAppView {
    const view = new DomainEntityAppView(client, row.id, options);
    view.name = options?.name ?? row.displayName;
    return view;
  }

  get table(): string { return `${this.client.schema}.${this.client.table}`; }

  /** Effective capabilities of the current user on the row's table (null until
   * the page has loaded). */
  get capabilities(): DG.DomainTableCapabilities | null { return this._capabilities; }

  get editors(): DomainFrameEditor[] {
    // Read during the base constructor, before `_details` is initialized.
    return super.editors.concat((this._details ?? []).map((g) => g.editor));
  }

  /** Loads the row and builds the page. */
  async load(): Promise<void> {
    try {
      const [values, capabilities, info] = await Promise.all([
        this.client.get(this._id) as Promise<{[key: string]: any}>,
        this.client.capabilities(),
        grok.dapi.domains.registry.tableInfo(this.table),
      ]);
      this._capabilities = capabilities;
      this._loaded = Object.assign({}, values);
      this.row = this.handler.rowFrom(values);
      this.name = this.options.name ?? this.row.displayName;
      this._host.innerHTML = '';
      this._host.appendChild(ui.divV([
        ui.h1(this.handler.getCaption(this.row as any)),
        ui.box(await this._tabs(info)),
      ], 'domain-ui-entity-page'));
      this.setUrlParams({[ENTITY_PARAM]: this._id});
      await this.buildRibbon();
    } catch (e: any) {
      this._host.innerHTML = '';
      this._host.appendChild(ui.divText(`Cannot open ${this.table} ${this._id}: ${e?.message ?? e}`));
      grok.shell.error(`${this.table}: ${e?.message ?? e}`);
    }
  }

  /** The row's own actions (Edit, Clone, Delete, Share, Watch, History, Copy
   * link), exactly the ones the current user may perform. */
  protected async buildRibbon(): Promise<void> {
    const row = this.row;
    if (row == null)
      return;
    const icons: HTMLElement[] = [ui.iconFA('arrow-left', () => this.close(), 'Back')];
    for (const action of await this.handler.getRibbonActions(row as any))
      icons.push(ui.iconFA(action.icon ?? 'bolt', () => this._run(action), action.name));
    this.setRibbonPanels([icons]);
  }

  /** Re-reads the row from the server and rebuilds the page. Prompts first when
   * a detail grid has pending changes. */
  async reload(): Promise<void> {
    if (!(await this.confirmDiscard('reload this page')))
      return;
    this._detachDetails();
    await this.load();
  }

  /**
   * Writes what the form changed, as one update guarded by the row's version.
   * A conflict goes through the platform's standard reload/overwrite dialog —
   * the same one the grid and the platform's own editors use.
   */
  async save(): Promise<boolean> {
    const row = this.row;
    if (row == null || this._capabilities == null)
      return false;
    const values = row.values;
    const changes: {[key: string]: any} = {};
    for (const name of this._capabilities.writableColumns)
      if (`${values[name] ?? ''}` !== `${this._loaded[name] ?? ''}`)
        changes[name] = values[name] ?? null;
    if (Object.keys(changes).length === 0) {
      grok.shell.info('No changes to save');
      return false;
    }
    return await this._update(changes, row.version);
  }

  detach(): void {
    this._detachDetails();
    super.detach();
  }

  // ─────────────────────── internals ─────────────────────────

  private async _update(changes: {[key: string]: any}, version?: number): Promise<boolean> {
    try {
      await this.client.update(this._id, changes as any,
        version == null ? undefined : {version: version});
      grok.shell.info('Saved');
      // Through the gate: rebuilding the page rebuilds the detail grids, whose
      // pending edits are none of this save's business.
      await this.reload();
      return true;
    } catch (e: any) {
      if (!(e instanceof DG.DomainVersionConflictError)) {
        grok.shell.error(`${e?.message ?? e}`);
        return false;
      }
      const decision = await DG.DomainObjectHandler.showConflictDialog(this.name);
      if (decision === 'reload') {
        await this.reload();
        return false;
      }
      if (decision === 'overwrite')
        return await this._update(changes, e.currentVersion ?? undefined);
      return false;
    }
  }

  private async _tabs(info: DG.DomainTableInfo): Promise<HTMLElement> {
    const row = this.row!;
    const tabs = ui.tabControl();
    tabs.addPane('Details', () => ui.wait(async () => await this._detailsPane()));
    for (const tab of await this.handler.getDetailTabs(row as any))
      tabs.addPane(tab.name, () => ui.wait(async () => await this._detailPane(tab)));
    if (info.audit)
      tabs.addPane('History', () => DG.DomainObjectHandler.auditPane(row));
    return tabs.root;
  }

  /** The property form of the row, plus Save / Revert when the caller may edit
   * it. The form itself is the platform's reflective one: writable columns only,
   * inputs and validation from the registry. */
  private async _detailsPane(): Promise<HTMLElement> {
    const form = await this.handler.renderEditor(this.row as any);
    if (this._capabilities?.canEdit !== true)
      return form;
    return ui.divV([form, ui.buttonsInput([
      ui.button('Save', () => this.save(), 'Save the changes'),
      ui.button('Revert', () => this.reload(), 'Reload the row, discarding the changes'),
    ])]);
  }

  /** One child table, scoped to this row: an editable grid whose new rows come
   * prefilled with the foreign key, saving its own transaction. */
  private async _detailPane(tab: DG.DomainDetailTab): Promise<HTMLElement> {
    const client = grok.dapi.domains.table(tab.table);
    const grid = await DomainGrid.create(client, {
      query: {filter: tab.filter},
      editable: this.options.editableDetails ?? true,
      defaults: {[tab.fkColumn]: this.row?.id},
    });
    this._details.push(grid);
    this.rebindEditors();
    return grid.root;
  }

  private async _run(action: DG.DomainAction): Promise<void> {
    const result = await action.run();
    // Delete closes the page (the row it shows is gone); everything else may
    // have changed the row, so the page re-reads it.
    if (action.name === 'Delete') {
      if (result !== false)
        this.close();
      return;
    }
    if (action.name !== 'Copy link' && action.name !== 'Open' && action.name !== 'History')
      await this.reload();
  }

  private _detachDetails(): void {
    for (const grid of this._details)
      grid.detach();
    this._details = [];
    this.rebindEditors();
  }
}
