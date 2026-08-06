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
 * export function issuesApp(_path?: string): DG.ViewBase {
 *   return new DomainAppView(grok.dapi.domains.table('grit.issue'));
 * }
 * ```
 *
 * The declared `path` input is what makes the app URL-addressable; the view reads
 * the deep link itself (`restoreFromUrl`), so the parameter stays unused — but a
 * function whose signature does not match its annotation is a `grok check`
 * finding waiting to happen.
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

import {acquireDomainContext, DomainFrameEditor, editorsOf, IDomainTableContext,
  IEditorHost} from './frame-editor';
import {DomainGrid} from './domain-grid';
import {EntityListWidget, EntityListMode, EntityListOptions,
  ENTITY_LIST_ROW_CAP} from './entity-list';
import {DomainForm} from './form';
import {actionChangesRow, DomainActionsProvider, DomainRowAction, domainHandler,
  isDeleteAction, isOpenAction, rowActions} from './handler';
import {discardFunc, saveFunc} from './actions';
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
export class AppView extends DG.ViewBase implements IEditorHost {
  /** The app's FuncCall, captured at construction — what makes a page opened
   * later (an entity page) keep the app's URL prefix and breadcrumbs. */
  readonly appCall: DG.FuncCall | undefined;

  private readonly _statusPanel = ui.div([], 'domain-ui-status-bar');
  private readonly _statusText = ui.divText('');
  private readonly _statusButtons: HTMLElement;
  private readonly _tracked: DomainFrameEditor[] = [];
  private _editorSubs: {unsubscribe(): void}[] = [];
  /** The prompt currently up, if any — see {@link confirmDiscard}. */
  private _prompting: Promise<boolean> | null = null;
  private _beforeUnload: ((e: BeforeUnloadEvent) => any) | null = null;

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

  /** Distinct per page kind: the platform keys saved layouts and view lookups on
   * it, and a row page is not a list page. */
  get type(): string { return 'domain-app'; }

  /** The widgets this page hosts, in layout order — where its editors, its
   * status and its functions come from. Subclasses override; NB this is not
   * `Widget.children`, which stays the platform's own DOM-derived view. */
  get widgets(): DG.Widget[] { return []; }

  /** Every editor whose pending changes this page answers for: the ones tracked
   * directly plus every hosted widget's, read on demand (a child's editor may
   * appear after the page was built, and rebuilt children must not go stale). */
  get editors(): DomainFrameEditor[] {
    // Read during the base constructor, before the subclass fields exist.
    const all = (this._tracked ?? []).slice();
    for (const child of this.widgets ?? [])
      for (const editor of editorsOf(child))
        if (!all.includes(editor))
          all.push(editor);
    return all;
  }

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
   *
   * Re-entrant by design: two gestures racing each other (a double-clicked Back,
   * a Back landing on top of a search) must ask ONCE and both act on the same
   * answer — a second dialog would offer a save of a batch the first one has
   * already discarded.
   */
  confirmDiscard(action: string = 'continue'): Promise<boolean> {
    if (this._prompting != null)
      return this._prompting;
    const done = () => this._prompting = null;
    const prompting = confirmDiscardChanges(this.editors, {action: action, subject: this.name})
      .then((proceed) => {
        done();
        return proceed;
      }, (e) => {
        done();
        throw e;
      });
    this._prompting = prompting;
    return prompting;
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

  /** {@link saveAll} — what the shared `Save` Func calls on a page. */
  save(): Promise<boolean> { return this.saveAll(); }

  /** {@link discardAll} — what the shared `Discard` Func calls on a page. */
  discard(): void { this.discardAll(); }

  /** {@link discardAll} — a page has no separate "back to the opened values". */
  reset(): void { this.discardAll(); }

  // ─────────────────────── the machine surface ─────────────────────────

  /**
   * The page's structure: every hosted widget's status, nested under its own
   * name — so a form's inputs are addressable as `DomainForm[0].title` from the
   * page — plus a summary naming what the page hosts and what is pending.
   */
  getWidgetStatus(): DG.IWidgetStatus {
    const parts: {[name: string]: Element} = {page: this.root};
    const inputs: DG.IInputStatus[] = [];
    const descriptions: string[] = [];
    let error: string | null = null;
    const children = this.widgets ?? [];
    children.forEach((child, index) => {
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
    descriptions.unshift(`The "${this.name}" page hosts ${children.length} ` +
      `widget${children.length === 1 ? '' : 's'}, with ${changes} unsaved ` +
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

  /** The hosted widgets' actions, merged and de-duplicated by name — what the
   * ribbon offers and what a caller (or the AI) invokes with
   * `f.apply({widget: page})`. */
  getFunctions(): DG.Func[] {
    const merged: DG.Func[] = [];
    for (const child of this.widgets ?? [])
      for (const func of child.getFunctions())
        if (!merged.some((f) => f.name === func.name))
          merged.push(func);
    // A page that tracks editors of its own (no widget answers for them) can
    // still be saved and discarded as a whole.
    if (this.editors.length > 0 && !merged.some((f) => f.name === 'Save'))
      merged.push(saveFunc(), discardFunc());
    return merged;
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
    this._syncBeforeUnload(changes > 0);
  }

  /** Closing the tab, reloading the page and following a link OUT of the client
   * never reach {@link confirmDiscard} — the browser owns that gesture, and the
   * only thing it takes is a `beforeunload` handler while something is pending.
   * Registered exactly while the page is dirty (an always-on one makes every
   * reload of a clean page ask). */
  private _syncBeforeUnload(dirty: boolean): void {
    if (typeof window === 'undefined')
      return;
    if (dirty && this._beforeUnload == null) {
      this._beforeUnload = (e: BeforeUnloadEvent) => {
        // The wording is the browser's own — the text a page supplies is ignored.
        e.preventDefault();
        e.returnValue = '';
        return '';
      };
      window.addEventListener('beforeunload', this._beforeUnload);
    }
    else if (!dirty && this._beforeUnload != null) {
      window.removeEventListener('beforeunload', this._beforeUnload);
      this._beforeUnload = null;
    }
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

  /** Where a back / forward navigation onto this page's URL lands. Reachable
   * from the router only on the page that {@link acceptsPath} — every other page
   * is routed BY it (see {@link DomainAppView.handlePath}); callers (and tests)
   * may still invoke it directly. */
  handlePath(_urlPath: string): void { this.restoreFromUrl(); }

  /**
   * ONE router per app: pages do not claim URLs.
   *
   * The platform hands a path to the FIRST open view that accepts it
   * (`_openViewsHandler`, `core/client/xamgle/lib/src/features/routing.dart`), and
   * every page of one app shares the app's URL prefix — so if each page claimed
   * the prefix, which one handled a Back would depend on the order the views
   * happen to sit in. The app's ROOT page ({@link DomainAppView}) is the only one
   * that says yes, and it decides what the URL means.
   */
  acceptsPath(_urlPath: string): boolean { return false; }

  /** Whether [urlPath] is inside this page's own URL prefix — what the app's
   * root page accepts on behalf of every page of the app. */
  protected ownsPath(urlPath: string): boolean {
    const own = this.pagePath;
    return own !== '' && `${urlPath}`.toLowerCase().startsWith(own.toLowerCase());
  }

  detach(): void {
    for (const sub of this._editorSubs)
      sub.unsubscribe();
    this._editorSubs = [];
    this._syncBeforeUnload(false);
    super.detach();
  }
}


/** Options of {@link DomainAppView} — the list options plus the view's own.
 * `confirmDiscard` is NOT among them: the page IS the gate, and it passes its own
 * to the list so one prompt covers every editor on the page. */
export interface DomainAppViewOptions extends AppViewOptions, Omit<EntityListOptions, 'confirmDiscard'> {
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
  /** The row pages this page has opened, by row id — what makes a back / forward
   * onto `?entity=` REUSE a page instead of stacking another copy of it. */
  private readonly _pages = new Map<string, DomainEntityAppView>();
  private _query: DG.DomainQuery;
  private _listSub: {unsubscribe(): void} | null = null;
  private _context: IDomainTableContext | null;
  /** Two loads can overlap (a Back landing on a slow one): the newest wins and
   * the older one drops whatever it built. */
  private _loading = 0;

  /** [context] is the prefetched handle a `domains.table(...).listView()` passes
   * in; without one the page resolves the table's metadata itself while it
   * loads. */
  constructor(client: DG.DomainTableClient, options?: DomainAppViewOptions,
    context?: IDomainTableContext) {
    super(options);
    this.client = client;
    this.options = options ?? {};
    this._context = context ?? null;
    this.handler = domainHandler(this.table);
    this.name = options?.name ?? this.client.table;
    this.box = true;
    this.root.appendChild(this._host);
    this._query = this._seedQuery();
    this._host.appendChild(ui.divText('Loading...'));
    this.load();
  }

  /** `'<schema>.<table>'`. */
  get table(): string { return `${this.client.schema}.${this.client.table}`; }

  /** What the page is showing, serializable: the shape a deep link, a saved
   * filter or an "open in Table View" carries. */
  get query(): DG.DomainQuery { return this._query; }

  /** The list, once it has loaded — the one widget this page hosts. */
  get widgets(): DG.Widget[] { return this.list == null ? [] : [this.list]; }

  /**
   * Builds (or rebuilds) the page from the URL.
   *
   * Rebuilding REPLACES the list, and with it whatever its grid had pending, so
   * a rebuild of a live list goes through the gate FIRST: a CANCEL leaves the
   * page exactly as it was — the current (dirty) list stays mounted, its batch
   * untouched, and the URL is republished from what is actually on screen,
   * because a browser Back has already changed it.
   */
  async load(): Promise<void> {
    const generation = ++this._loading;
    const params = this.urlParams();
    // Applied only once the load has committed: a refused one must not leave the
    // page reporting a query it does not show.
    const query = this._queryFromParams(params);
    const entity = params[ENTITY_PARAM];
    try {
      // The first load has nothing to discard; every later one asks for the whole
      // page, not just the list's grid.
      if (this.list != null && !(await this.confirmDiscard('reload the list'))) {
        this.syncUrl();
        return;
      }
      const context = await this._resolveContext();
      if (generation !== this._loading)
        return;
      const list = new EntityListWidget(context, Object.assign({}, this.options, {
        query: this._spec(query),
        mode: this._modeFromParams(params),
        // ONE gate for the whole page: the list's own rebuilds (search, mode,
        // query) prompt for every editor the page owns, not just its grid's.
        confirmDiscard: () => this.confirmDiscard('reload the list'),
        onOpen: (row: DG.DomainRow) => this.open(row),
      }));
      await list.ready;
      if (generation !== this._loading) {
        list.detach();
        return;
      }
      this._query = query;
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
      if (generation !== this._loading)
        return;
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

  /**
   * Opens a row's page inside the app (the default of {@link EntityListWidget}'s
   * open action). Pass the row or its id.
   *
   * FIND-OR-ACTIVATE: a row that already has a page gets that page back, not a
   * second copy of it — clicking the same row twice, and every back / forward
   * onto its URL, land on the one page (with whatever it had pending).
   */
  open(row: DG.DomainRow | string): DomainEntityAppView | null {
    if (this.options.onOpen != null && typeof row !== 'string') {
      this.options.onOpen(row);
      return null;
    }
    const id = typeof row === 'string' ? row : row.id;
    const existing = this._pages.get(id);
    if (existing != null) {
      this.activate(existing);
      return existing;
    }
    const options: DomainEntityAppViewOptions = {
      helpUrl: this.helpUrl ?? undefined,
      actions: this.options.actions,
      onClosed: () => this._pages.delete(id),
    };
    const view = typeof row === 'string'
      ? new DomainEntityAppView(this.client, row, options, this._context ?? undefined)
      : DomainEntityAppView.forRow(this.client, row, options, this._context ?? undefined);
    // Inheriting the prefix keeps the entity page inside the app's URL — a page
    // opened by a click is created outside the app function's call.
    view.parentView = this;
    if (this.appCall != null)
      view.parentCall = this.appCall;
    this._pages.set(id, view);
    grok.shell.addView(view);
    this.activate(view);
    return view;
  }

  /**
   * THE router of the app: every page of it shares this prefix, and this is the
   * one view the platform hands the path to ({@link acceptsPath}).
   *
   * `entity=` names a row page: it is found or opened, never stacked. Without it
   * the URL is the list's own, so the list restores itself from it (through the
   * gate, as every rebuild does).
   */
  handlePath(_urlPath: string): void {
    const entity = this.urlParams()[ENTITY_PARAM];
    if (entity == null || entity === '')
      this.restoreFromUrl();
    else
      this.open(entity)?.restoreFromUrl();
  }

  /** Accepts the whole app prefix — on behalf of every page of the app, which is
   * what makes back / forward independent of the order the views sit in. */
  acceptsPath(urlPath: string): boolean { return this.ownsPath(urlPath); }

  /**
   * Makes [view] the current one AFTER the router is done: it sets the current
   * view to whichever view accepted the path — this one — the moment
   * {@link handlePath} returns, so an activation made inside it would be undone.
   * A microtask lands right after that, and before the browser paints.
   */
  protected activate(view: DG.ViewBase): void {
    Promise.resolve().then(() => {
      if (!view.closing)
        grok.shell.v = view as any;
    });
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
    this._pages.clear();
    super.detach();
  }

  /** The prefetched table metadata: the handle's when the page was built from
   * one, resolved once (and kept) otherwise. */
  private async _resolveContext(): Promise<IDomainTableContext> {
    if (this._context == null)
      this._context = await acquireDomainContext(this.client);
    return this._context;
  }

  /** The page's query for the URL it was opened with: the deep link when it
   * carries one, `options.query` otherwise — so `new DomainAppView(client,
   * {query})` (the "my campaigns" page) lists what it was asked for, while a
   * pasted link still wins over the default. */
  private _queryFromParams(params: {[key: string]: string}): DG.DomainQuery {
    let query: DG.DomainQuery;
    try {
      query = DG.DomainQuery.fromUrlParams(this.client.schema, this.client.table, params);
    } catch (e: any) {
      grok.shell.warning(`Ignoring the URL query: ${e?.message ?? e}`);
      return this._seedQuery();
    }
    return Object.keys(query.toUrlParams()).length === 0 ? this._seedQuery() : query;
  }

  /** `options.query` as a {@link DG.DomainQuery} — the inverse of
   * {@link DG.DomainQuery.toSpec}, so the seeded page publishes the URL a deep
   * link would restore it from. */
  private _seedQuery(): DG.DomainQuery {
    const query = new DG.DomainQuery({schema: this.client.schema, table: this.client.table});
    const spec = this.options?.query;
    if (spec == null)
      return query;
    if (spec.filter != null)
      // A grammar string travels as itself; a condition node or tree as the
      // per-element JSON escape hatch `toSpec()` decodes back.
      query.filters = [typeof spec.filter === 'string' ? spec.filter : JSON.stringify(spec.filter)];
    if (spec.sort != null && spec.sort !== '')
      query.orderBy = `${spec.sort}`.split(',');
    if (spec.columns != null)
      query.columns = spec.columns.slice();
    if (spec.expand != null)
      query.joins = spec.expand.slice();
    if (spec.limit != null)
      query.limit = spec.limit;
    if (spec.offset != null)
      query.offset = spec.offset;
    return query;
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
  /** Customizes the row's actions — see {@link EntityListOptions.actions}. */
  actions?: DomainActionsProvider;
  /** Called when the page is closed — how the app that opened it forgets it
   * (see {@link DomainAppView.open}). */
  onClosed?: () => void;
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
  private _context: IDomainTableContext | null;
  private _details: DomainGrid[] = [];
  /** The property form of the row — THE writer of everything the Details pane
   * edits (dirty state, validation, the 409 flow and the save come from its
   * one-row editor). Null until the page has loaded. */
  private _form: DomainForm | null = null;

  constructor(client: DG.DomainTableClient, id: string, options?: DomainEntityAppViewOptions,
    context?: IDomainTableContext) {
    super(options);
    this.client = client;
    this.options = options ?? {};
    this._context = context ?? null;
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
    options?: DomainEntityAppViewOptions, context?: IDomainTableContext): DomainEntityAppView {
    const view = new DomainEntityAppView(client, row.id, options, context);
    view.name = options?.name ?? row.displayName;
    return view;
  }

  get type(): string { return 'domain-entity-app'; }

  get table(): string { return `${this.client.schema}.${this.client.table}`; }

  /** Effective capabilities of the current user on the row's table (null until
   * the page has loaded). */
  get capabilities(): DG.DomainTableCapabilities | null {
    return this._context?.capabilities ?? null;
  }

  /** The property form of the row, once the page has loaded. */
  get form(): DomainForm | null { return this._form; }

  /** The form and the detail grids — everything on the page that edits. */
  get widgets(): DG.Widget[] {
    const children: DG.Widget[] = this._form == null ? [] : [this._form];
    return children.concat(this._details ?? []);
  }

  /** Loads the row and builds the page. */
  async load(): Promise<void> {
    try {
      const [values, context] = await Promise.all([
        this.client.get(this._id) as Promise<{[key: string]: any}>,
        this._resolveContext(),
      ]);
      this.row = this.handler.rowFrom(values);
      this.name = this.options.name ?? this.row.displayName;
      // Before the tabs (and the ribbon they precede): the form is the page's
      // writer, and what it offers has to be on the ribbon whether or not the
      // Details pane has been rendered yet.
      this._buildForm(context);
      this._host.innerHTML = '';
      this._host.appendChild(ui.divV([
        ui.h1(this.handler.getCaption(this.row as any)),
        ui.box(await this._tabs(context.info)),
      ], 'domain-ui-entity-page'));
      this.setUrlParams({[ENTITY_PARAM]: this._id});
      await this.buildRibbon();
    } catch (e: any) {
      this._host.innerHTML = '';
      this._host.appendChild(ui.divText(`Cannot open ${this.table} ${this._id}: ${e?.message ?? e}`));
      grok.shell.error(`${this.table}: ${e?.message ?? e}`);
    }
  }

  /** The prefetched table metadata: the app's when the page was opened from one,
   * resolved once (and kept) otherwise. */
  private async _resolveContext(): Promise<IDomainTableContext> {
    if (this._context == null)
      this._context = await acquireDomainContext(this.client);
    return this._context;
  }

  /**
   * The row's own actions (Edit, Clone, Delete, Share, Watch, History, Copy
   * link), exactly the ones the current user may perform — WITHOUT the Open one,
   * which opens the very row this page is showing — plus the form's own
   * functions (Save / Discard / Reset) when the caller may edit it.
   */
  protected async buildRibbon(): Promise<void> {
    const row = this.row;
    if (row == null)
      return;
    const icons: HTMLElement[] = [ui.iconFA('arrow-left', () => this.back(), 'Back')];
    for (const action of await rowActions(this.handler, row, this.options.actions))
      if (!isOpenAction(action))
        icons.push(ui.iconFA(action.icon ?? 'bolt', () => this.runAction(action), action.name));
    const panels: HTMLElement[][] = [icons];
    const form = this._form;
    if (form != null) {
      const buttons = form.getFunctions().map((func) =>
        ui.button(func.friendlyName ?? func.name, () => func.apply({widget: form}),
          `${func.description ?? ''}`));
      if (buttons.length > 0)
        panels.push(buttons);
    }
    this.setRibbonPanels(panels);
  }

  /** Leaves the page — through the gate, because closing it drops every detail
   * grid and `detach()` cannot ask anything. Resolves to whether it closed. */
  async back(): Promise<boolean> {
    if (!(await this.confirmDiscard('leave this page')))
      return false;
    this.close();
    return true;
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
   * Writes what the Details form changed — one transaction, guarded by the row's
   * version, with the standard 409 reload/overwrite dialog and the per-field
   * markers, all from the form's own one-row editor. Nothing here diffs values
   * by hand any more.
   */
  async save(): Promise<boolean> {
    if (this._form == null)
      return false;
    if (!this._form.isDirty) {
      grok.shell.info('No changes to save');
      return false;
    }
    return await this._form.save();
  }

  /** Puts the Details form back to the values the row was loaded with. */
  discard(): void { this._form?.discard(); }

  detach(): void {
    this._detachDetails();
    this._form?.detach();
    this._form = null;
    super.detach();
    this.options.onClosed?.();
  }

  /** Back / forward onto this page's URL: it is already showing that row, so the
   * restore is a re-read of it (through the gate, like every other rebuild). */
  restoreFromUrl(): void {
    this.reload();
  }

  // ─────────────────────── internals ─────────────────────────

  private async _tabs(info: DG.DomainTableInfo): Promise<HTMLElement> {
    const row = this.row!;
    const tabs = ui.tabControl();
    tabs.addPane('Details', () => this._form!.root);
    // waitBox, not wait: a canvas grid has no intrinsic size, so its host has to
    // BE a sized box — inside a plain div it renders at the viewer's default
    // 400x300 whatever the pane's width.
    for (const tab of await this.handler.getDetailTabs(row as any))
      tabs.addPane(tab.name, () => ui.waitBox(async () => await this._detailPane(tab)));
    if (info.audit)
      tabs.addPane('History', () => DG.DomainObjectHandler.auditPane(row));
    return tabs.root;
  }

  /**
   * The property form of the row — a {@link DomainForm} in edit mode, so the
   * page inherits ONE writer for it: dirty state and the pending count on the
   * status bar, per-field validation from the registry, reference pickers, the
   * standard 409 flow, and a save that is one transaction. Save / Discard ride
   * the ribbon (the form's own functions) and the page's gate.
   */
  private _buildForm(context: IDomainTableContext): void {
    this._form?.detach();
    const form = new DomainForm(context, {id: this._id});
    this._form = form;
    // The editor appears when the row has loaded; the status bar has to pick it
    // up then, not only when the page was built.
    form.ready.then(() => {
      if (this._form === form)
        this.rebindEditors();
    });
    this.rebindEditors();
  }

  /** One child table, scoped to this row: an editable grid whose new rows come
   * prefilled with the foreign key, saving its own transaction. */
  private async _detailPane(tab: DG.DomainDetailTab): Promise<HTMLElement> {
    const client = grok.dapi.domains.table(tab.table);
    // An explicit cap, the list's own: without one the server applies its default
    // page size (100) and the pane quietly shows a fraction of the children.
    const [grid, total] = await Promise.all([
      DomainGrid.create(client, {
        query: {filter: tab.filter, limit: ENTITY_LIST_ROW_CAP},
        editable: this.options.editableDetails ?? true,
        defaults: {[tab.fkColumn]: this.row?.id},
      }),
      client.count(tab.filter).catch(() => null),
    ]);
    this._details.push(grid);
    this.rebindEditors();
    const shown = grid.dataFrame.rowCount;
    if (total == null || total <= shown)
      return ui.box(grid.root);
    // Same words as the list's banner — the same situation. The wrapper is TAGGED:
    // a plain `ui.divV` around a `ui.box` hands the box ui.css's default 400x300
    // (`div.ui-div > div.ui-box`), which is the collapse the tab host itself was
    // fixed for — and it only shows up past the row cap, which is why it survived.
    return ui.divV([
      ui.divText(`Showing the first ${shown} of ${total} rows — refine the filter to ` +
        'narrow the results.', 'grok-domains-cap-banner'),
      ui.box(grid.root),
    ], 'domain-ui-detail-pane');
  }

  /**
   * Runs a row action the way the ribbon does: Delete closes the page (the row
   * it shows is gone), anything else that MAY have changed the row re-reads it —
   * through the gate, so a detail grid's pending edits are never dropped
   * silently. A cancelled action (`run()` returning false) changes nothing.
   */
  async runAction(action: DomainRowAction): Promise<void> {
    const result = await action.run();
    // Delete closes the page (the row it shows is gone); everything else may
    // have changed the row, so the page re-reads it.
    if (isDeleteAction(action)) {
      if (result !== false)
        this._closeDeleted();
      return;
    }
    // A CANCELLED action changed nothing: reloading behind it would pop the
    // unsaved-changes prompt for a dialog the user has just dismissed.
    if (result !== false && actionChangesRow(action))
      await this.reload();
  }

  /** The row is gone: its child rows went with it, so their pending edits can
   * never be saved. They are dropped EXPLICITLY, and said out loud — letting
   * `detach()` swallow them is exactly the silence the gate exists to prevent. */
  private _closeDeleted(): void {
    const changes = this.editors.reduce((n, e) => n + (e.isDirty ? e.changeCount : 0), 0);
    if (changes > 0) {
      this.discardAll();
      grok.shell.warning(
        `Discarded ${changes} unsaved change${changes === 1 ? '' : 's'} — the row was deleted.`);
    }
    this.close();
  }

  private _detachDetails(): void {
    for (const grid of this._details)
      grid.detach();
    this._details = [];
    this.rebindEditors();
  }
}
