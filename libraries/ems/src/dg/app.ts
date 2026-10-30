import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {Control, signal, computed, batch, ReadonlySignal, Signal, button, div, divH, divV, span,
  icon, Filters, allowedActions, actionsMenu, iconButton, ButtonGroup, Breadcrumbs, Rows, Menu, LabelColumn, Reason,
  text} from '@datagrok-libraries/u2';
import type {FilterGroup, IWidgetStatus, Action, RowView, Scope} from '@datagrok-libraries/u2';
import type {DomainSource} from '../core/domain-source.js';
import {capital, nounOf} from '../core/domain-backend.js';
import type {DomainTableInfoLike} from '../core/domain-backend.js';
import {SharedSession, confirmDiscard} from '../core/session.js';
import {domains, DomainTable} from './facade.js';
import {DomainList} from './list.js';
import type {DomainListMode} from './list.js';
import {SAVE_SHORTCUTS, isSaveKey, DomainForm} from './form.js';
import {DomainFilters} from './filters.js';
import {SavedSearches} from './saved-searches.js';
import type {DomainFiltersOptions} from './filters.js';
import type {DomainChildrenOptions} from './children.js';
import {DomainAddress} from './address.js';
import {DomainErrors} from './errors.js';
import {DgDomainTable} from './backend.js';
import {DomainDataTable} from './grid.js';
import {domainOutput, outputKey} from './output.js';
import {PropertyDialog} from './authoring/property-dialog.js';
import {ViewSync} from '@datagrok-libraries/u2/dg';
import type {ViewState} from '@datagrok-libraries/u2/dg';

export type DomainAppPage = 'list' | 'entity';

export type DomainAppMode = 'live' | 'trash';

/** What the status bar says while a live source is behind the server ({@link DomainApp.refresh}). */
const STALE = 'Data changed — Refresh';

/** Below this width the ribbon has no room for the presets. */
const COMPACT_WIDTH = 700;

/** The trash reads newest-deleted first: a soft delete is a write, so it stamps `updated_on`. */
const TRASH_SORT = '!updated_on';

/** What the app's status bar says, as VALUES rather than one composed sentence: a view spec renders
 * the slots it wants instead of parsing text. */
export interface DomainAppStatus {
  /** What the page holds: "50 of 2,077 issues", "3 issues", "1 deleted issue", or the row's caption on
   * the entity page. `''` while there is nothing to say yet. */
  count: string;
  /** How many rows a bulk action would act on; 0 on the entity page and with nothing picked. */
  selected: number;
  /** Something true that is not a failure: "Data changed — Refresh", "Filter not applied",
   * "3 unsaved changes", "2 restores pending". Null when there is nothing to say. */
  notice: string | null;
  /** Why the pending changes are stuck, or why the load failed. Null when nothing is wrong. */
  problem: string | null;
}

/** The ribbon as slots rather than positions: the `ribbon` option adds to the slot it means, and a
 * view spec composes them without counting groups. */
export interface DomainRibbon {
  /** New, Save, Discard, the ⋯ menu, Refresh — the app's own verbs. */
  main: (Control | HTMLElement)[];
  /** What acts on the collection: the search box. */
  tools: (Control | HTMLElement)[];
  /** A query switch ({@link DomainApp.presets}); empty by default. */
  presets: (Control | HTMLElement)[];
}

export interface DomainAppOptions {
  table: DomainTable;
  /** The view path the app's paths extend, e.g. `/apps/Grit/Issues`. */
  base: string;
  /** The list's initial query. */
  query?: string;
  /** The list's initial order, `'!number,title'` (`!` descending). */
  sort?: string;
  pageSize?: number;
  /** Whether the list follows the server (`DomainSourceOptions.live`), ON by default: it probes
   * the table every {@link liveMs} and reloads while the session is clean, marking the page
   * {@link DomainApp.stale} while it is not. The entity page is one row and is not polled. */
  live?: boolean;
  liveMs?: number;
  mode?: DomainListMode;
  itemHeight?: number;
  /** What every draft the app creates starts with — the preset a view that opened the app carries
   * (a location's id on a containers app). Threaded into the entity page's draft source. */
  defaults?: Record<string, unknown>;
  /** The form's columns, in this order; every column by default. */
  include?: string[];
  exclude?: string[];
  /** The entity page's panes under the form, on by default: the child collections (a tab per
   * referring table; options narrow them) and the row's history — where the table keeps one
   * (`support.audit`); a table that keeps none has no pane. */
  children?: boolean | DomainChildrenOptions;
  history?: boolean;
  /** The query box over the list (default), or the condition builder; `placement: 'left'` puts it in a pane
   * beside the list. */
  filters?: DomainFiltersOptions & {placement?: 'top' | 'left'};
  /** Key → the name of an action in `table.actions`, run over the current row while the app has
   * the focus: `{'Ctrl+Shift+C': 'Close'}` (modifiers in the order Ctrl, Alt, Shift; Cmd counts
   * as Ctrl). */
  shortcuts?: Record<string, string>;
  /** The ribbon from the one the app builds, once:
   * `(base, app) => ({...base, presets: [app.presets(['Mine', 'assignee = $me'])]})`. */
  ribbon?: (base: DomainRibbon, app: DomainApp) => DomainRibbon;
  /** The Actions menu from the one the app builds, on every open — extend `base` rather than replace it. */
  menu?: (base: Action[], app: DomainApp) => Action[];
  /** What the app wires beyond its options, run once it is built under its owner `self`: validators,
   * effects over its pages, a hosted filter on `listSource`. */
  setup?: (app: DomainApp, self: Scope) => void;
  /** The list page as a table (the ribbon's list/table switch): its row height, and a column's cell where the app
   * draws it better than its type. */
  grid?: {rowHeight?: number, cells?: Record<string, (row: RowView) => HTMLElement | string>};
}

/** The app shape over one table: a list page and an entity page under one `SharedSession`, the URL following
 * the page, and the unsaved gate in front of every move between them. `DomainTable.app()` hosts it in a view
 * (`ViewHost`). An app's own parts are options, not a subclass: `ribbon`, `menu`, `presets()`, `shortcuts`, `setup`. */
export class DomainApp extends Control {
  /** The `entity` that stands for a draft: `?entity=new` opens a create page. */
  static readonly NEW = 'new';
  /** Every app alive — what find-or-activate looks through. */
  static readonly live = new Set<DomainApp>();

  readonly table: DomainTable;
  /** The view path the app's paths extend — the shell's own app route once the view is docked
   * ({@link rebase}), the option it was built with until then. */
  get base(): string { return this._base.value; }
  /** The session every source the app makes joins — the ambient one it was built under, else its own. */
  readonly session: SharedSession;
  readonly page: ReadonlySignal<DomainAppPage>;
  /** The row the entity page shows, {@link NEW} for a draft; null on the list page. */
  readonly entity: ReadonlySignal<string | null>;
  readonly list: DomainList;
  readonly listSource: DomainSource;
  /** One row per entity page — a source of one sharing the session, replaced with the entity. */
  readonly entitySource: ReadonlySignal<DomainSource | null>;
  readonly form: ReadonlySignal<DomainForm | null>;
  /** The entity form's label column: content lined up with the form reads its `width`. */
  readonly labelColumn = new LabelColumn();
  readonly breadcrumbs: Breadcrumbs;
  /** The entity page below the form — where the children and history panes go. */
  readonly panes: HTMLElement;
  /** `base?q=…&search=…` on the list page; on the entity page `base/<key>` where the app
   * addresses a row by a path segment ({@link entityPath}), else `base?entity=<id>`. */
  readonly path: ReadonlySignal<string>;
  /** What the status bar says, slot by slot: the count, the selection, a notice, a problem. */
  readonly status: ReadonlySignal<DomainAppStatus>;
  /** The session's pending changes when there are any, else the page's source. */
  readonly summary: ReadonlySignal<string>;
  /** Whether the page's collection is behind the server — a `live` source saw it move while this
   * session had unsaved changes, so it did not reload. {@link refresh} is the way forward. */
  readonly stale: ReadonlySignal<boolean>;
  /** What is wrong with the filter in the box while the rows do not answer it — said in the
   * status bar, where a tooltip on a red box is not looked at. */
  readonly filterProblem: ReadonlySignal<string | null>;
  /** The query box over the list. It is a ROW OF THE LIST PAGE, not a ribbon item: the shell's
   * ribbon is one fixed 32px line with `overflow: hidden`, and a query is longer than that line
   * has to give. */
  readonly filters: DomainFilters;
  /** Which collection the list answers: the live rows, or the table's trash — where every path
   * the app emits carries `?trash=1`. Written by the ⋯ menu ({@link setMode}) and read back from
   * the address bar by {@link open}. */
  readonly mode: ReadonlySignal<DomainAppMode>;
  /** {@link mode} as the one question most of the app asks of it. */
  readonly trash: ReadonlySignal<boolean>;
  /** See {@link DomainAppOptions.shortcuts}. */
  shortcuts: Record<string, string>;
  /** The view is too narrow for the whole ribbon: the presets fold into the Actions menu. */
  readonly compact = signal(false);
  /** The list page's rows as the list or as a table — the ribbon's switch. */
  readonly view = signal<'list' | 'table'>('list');
  /** The table view's columns: the Output and Aggregations under the filters. */
  readonly output: ReturnType<typeof domainOutput>;

  private readonly _mode = signal<DomainAppMode>('live');
  private readonly _count = signal('');
  private readonly _filterProblem = signal<string | null>(null);
  private readonly _base: Signal<string>;
  private readonly _page = signal<DomainAppPage>('list');
  private readonly _entity = signal<string | null>(null);
  private readonly _entitySource = signal<DomainSource | null>(null);
  private readonly _form = signal<DomainForm | null>(null);
  private readonly _formHost: HTMLElement;
  private _panes: Control[] = [];
  private _ribbon: DomainRibbon | undefined;
  private _unloadGuard: ((e: BeforeUnloadEvent) => void) | undefined;
  private _deleting: DomainSource | null = null;
  private readonly _presets: Action[] = [];
  private readonly _sync: ViewSync;
  private _opening: {address: string, done: Promise<boolean>} | undefined;

  private static readonly _byView = new WeakMap<DG.ViewBase, DomainApp>();
  private static readonly _views = new WeakMap<DomainApp, DG.ViewBase>();

  constructor(private readonly _options: DomainAppOptions) {
    super();
    const table = _options.table;
    this.table = table;
    this._base = signal(_options.base);
    this.session = SharedSession.ambient ?? new SharedSession();
    this.page = this._page;
    this.entity = this._entity;
    this.entitySource = this._entitySource;
    this.form = this._form;
    this.mode = this._mode;
    this.trash = computed(() => this._mode.value === 'trash');
    this.filterProblem = this._filterProblem;
    this.shortcuts = _options.shortcuts ?? {};
    this.root.classList.add('u2-domain-app');
    this.root.dataset.u2 = 'domain-app';
    const onKeyDown = (e: KeyboardEvent) => this._shortcut(e);
    this.root.addEventListener('keydown', onKeyDown);
    this.own(() => this.root.removeEventListener('keydown', onKeyDown));
    // the chrome the view's host places sits outside the app root, so Save has to be caught above both
    const onSave = (e: KeyboardEvent) => this._save(e);
    document.addEventListener('keydown', onSave);
    this.own(() => document.removeEventListener('keydown', onSave));
    DomainApp.live.add(this);
    this.own(() => DomainApp.live.delete(this));

    this.listSource = this._source({query: _options.query ?? '', sort: _options.sort, pageSize: _options.pageSize,
      live: _options.live ?? true, liveMs: _options.liveMs});
    this.own(() => this.listSource.dispose());
    this.list = this.build(() => new DomainList(this.listSource, {mode: _options.mode,
      itemHeight: _options.itemHeight, create: () => void this.create()}));
    // Enter on a row is the way into the entity page
    let seen = this.listSource.activate.peek();
    this.effect(() => {
      const bumped = this.listSource.activate.value;
      if (bumped === seen)
        return;
      seen = bumped;
      const row = this.listSource.currentRow.peek();
      if (row !== null && !Rows.isDraft(row))
        void this.goTo('entity', row.id);
    });

    const title = DomainApp.titleOf(table.info);
    this.breadcrumbs = this.build(() => new Breadcrumbs({items: [title], onClick: (index) => {
      // in the trash the root crumb reads "Issues › Trash" and is the way back out of it
      if (index === 0)
        void (this.trash.peek() ? this.setTrash(false) : this.goTo('list'));
    }}));
    this._formHost = div([], 'u2-domain-app-form u2-form-fill');
    this.panes = div([], 'u2-domain-app-panes');
    this.panes.dataset.u2Part = 'panes';
    this.filters = this.build(() => domains.filters(this.listSource, _options.filters));
    // the rows answer the filter that was applied, not the one being written: while the box is
    // invalid they are shown as stale and the count stands down
    this.effect(() => {
      const problem = this.filters.problem.value;
      this._filterProblem.value = problem;
      this.list.stale.value = problem !== null;
    });
    const plural = nounOf(table.info);
    const banner = this.build(() => divH([icon('info-circle'),
      span(`Deleted ${plural}. Restore what you need, then Save.`)], 'u2-domain-app-trash-banner'));
    banner.dataset.u2Part = 'trash-banner';
    this.effect(() => banner.hidden = !this.trash.value);
    this.output = this.build(() => domainOutput(this.listSource, table, _options.exclude ?? [],
      computed(() => this.view.value === 'table')));
    const listPage = divV([banner, ...(_options.filters?.placement === 'left' ?
      [this.filters.pane(this.list.root, [this.output.root])] : [this.filters.root, this.output.root, this.list.root])],
    'u2-domain-app-list');
    this._tableView();
    listPage.dataset.u2Part = 'list-page';
    const entityPage = divV([this._formHost, this.panes], 'u2-domain-app-entity');
    entityPage.dataset.u2Part = 'entity-page';
    // above both pages: the entity page always says which row, the list page only in the trash —
    // on a plain list the table's name is the view title already
    this.root.append(this.breadcrumbs.root, listPage, entityPage);
    this.effect(() => {
      const page = this._page.value;
      listPage.hidden = page !== 'list';
      entityPage.hidden = page !== 'entity';
    });
    this.effect(() => {
      const trash = this.trash.value;
      const entity = this._page.value === 'entity';
      this.breadcrumbs.root.hidden = !entity && !trash;
      if (trash) {
        this.breadcrumbs.setItems([title, 'Trash']);
        return;
      }
      const source = this._entitySource.value;
      // the row is a keyed proxy: the cells the save wrote back (the autonumbered name among
      // them) move under it without `currentRow` changing, so the caption follows the rows tick
      source?.rows.items.value;
      const row = source?.currentRow.value ?? null;
      // a row the table does not have is named by the address that asked for it
      const caption = row !== null ? this._caption(source!, row) :
        source?.state.value === 'ready' && !source.isDraft ? this._entity.peek() ?? '…' : '…';
      this.breadcrumbs.setItems([title, caption]);
    });
    // the entity page over a draft: once saved, the page is the row's
    this.effect(() => {
      const source = this._entitySource.value;
      const row = source?.currentRow.value ?? null;
      if (row !== null && source!.isDraft && this._entity.peek() === DomainApp.NEW && !Rows.isDraft(row))
        this._entity.value = row.id;
    });
    // and once the draft is discarded — by the Discard button or by a gate's DISCARD — the page
    // has no row left, so the app goes back to the list. `had` keeps the load's own empty window
    // (the frame arrives before `newRow` fills it) from counting as a discard; the move is
    // deferred past the discard that is running, which reads dirty until it has settled. A saved
    // delete of the row the draft became empties it too, mid-batch: that move is the next effect's.
    let had = false;
    this.effect(() => {
      const source = this._entitySource.value;
      if (source === null || !source.isDraft) {
        had = false;
        return;
      }
      const rows = source.rows.items.value.length;
      if (rows === 0 && had && !this.session.isSaving.peek())
        queueMicrotask(() => void this.goTo('list'));
      had = rows > 0;
    });

    // a row deleted from its own page: once the save has taken it, back to the list
    this.effect(() => {
      const source = this._entitySource.value;
      const empty = source !== null && source.state.value === 'ready' && source.rows.items.value.length === 0;
      const settled = !this.session.isDirty.value && !this.session.isSaving.value;
      if (source !== this._deleting)
        this._deleting = null;
      if (this._deleting === null || !settled)
        return;
      this._deleting = null;
      if (empty)
        queueMicrotask(() => void this.goTo('list'));
    });

    this.path = computed(() => {
      const entity = this._entity.value;
      if (this._page.value === 'entity' && entity !== null) {
        const address = this._entityAddress(entity);
        const source = this._entitySource.value;
        // the draft sentinel is spelled as is; a row's key is escaped so `~new` is the row "new"
        const segment = DomainApp._isDraft(source, source?.currentRow.value ?? null) ? address :
          DomainAddress.segment(address);
        return this.entityPath ? `${this.base}/${segment}` : `${this.base}?entity=${segment}`;
      }
      const q = this.listSource.query.value;
      const search = this.listSource.search.value;
      const sort = this.listSource.sort.value;
      let path = this._withTrash(Filters.queryPath(this.base, typeof q === 'string' ? q : Filters.format(q)));
      if (search !== '')
        path = `${path}${path.includes('?') ? '&' : '?'}search=${encodeURIComponent(search)}`;
      // only an order the user picked; the trash's is checked too, as the mode flips before the sort
      if (!this.trash.value && sort !== TRASH_SORT && sort !== (this._options.sort ?? ''))
        path = `${path}${path.includes('?') ? '&' : '?'}sort=${encodeURIComponent(sort)}`;
      return path;
    });
    this._sync = new ViewSync({
      path: this.path,
      draftEntity: DomainApp.NEW,
      state: () => this._state(),
      open: (address) => this.open(address),
      owns: (path) => DomainAddress.restOf(path, this._bases) !== null,
      isCurrentView: () => this._isCurrentView(),
    }, {owner: this});
    // the source folds a refusal, the pending changes and the poll's news into its one summary line, which are
    // slots of their own here: the count is taken while the source says nothing else, and stands until then
    this.effect(() => {
      const entity = this._page.value === 'entity';
      const source = entity ? this._entitySource.value : this.listSource;
      if (source === null) {
        this._count.value = '';
        return;
      }
      // the rows tick is read for the same reason the breadcrumb reads it — a save writes the
      // server's values into the cells of a proxy `currentRow` never stops pointing at
      source.rows.items.value;
      // the entity page shows one row, not a collection: it says which row, never "1 substance"
      const row = entity ? source.currentRow.value : null;
      if (row !== null && !Rows.isDraft(row))
        this._count.value = this._caption(source, row);
      else if (source.state.value === 'error' || entity && !source.isDraft && source.state.value === 'ready')
        this._count.value = '';
      else if (source.error.value === undefined && source.changeCount.value === 0)
        this._count.value = source.summary.value;
    });
    this.status = computed(() => {
      const entity = this._page.value === 'entity';
      const source = entity ? this._entitySource.value : this.listSource;
      const dirty = this.session.isDirty.value;
      const pending = dirty ? this._pendingSummary() : null;
      const saved = entity ? null : SavedSearches.title(this.listSource);
      const count = saved === null || this._count.value === '' ? this._count.value : `${saved} · ${this._count.value}`;
      // what a bulk action would act on, beside the count — the platform grid's own convention
      const selected = entity ? 0 : this.listSource.selection.value.length;
      const problem = source !== null && source.error.value !== undefined ? source.summary.value : null;
      // the rows answer the filter that was applied, not the one in the box: while that one is
      // invalid they are shown stale, and a refusal of the old one is not the news
      if (!entity && this.list.stale.value && !dirty)
        return {count, selected, notice: null, problem: this._filterProblem.value ?? 'Filter not applied'};
      // a live source saw the collection move while this session had changes to protect: what is
      // pending still matters, and so does that the rows under it are no longer the server's
      if (source !== null && source.stale.value)
        return {count, selected, notice: pending === null ? STALE : `${pending} — ${STALE}`, problem};
      return {count, selected, notice: pending, problem};
    });
    this.summary = computed(() => DomainApp.summaryOf(this.status.value));
    this.stale = computed(() =>
      (this._page.value === 'entity' ? this._entitySource.value : this.listSource)?.stale.value ?? false);
    this._wireMode();
    this.own(() => this._close());
    const fit = new ResizeObserver(() => {
      const width = this.root.clientWidth;
      if (width > 0)
        this.compact.value = width < COMPACT_WIDTH;
    });
    fit.observe(this.root);
    this.own(() => fit.disconnect());
    this.build(() => _options.setup?.(this, this.scope));
  }

  /** Re-points the app at the route the shell actually mounted its view at (`/apps/Stockroom`):
   * every path the app produces moves with it, and find-or-activate looks the table up under the
   * new base. See `DomainTable.app()`, which resolves it from the docked view. */
  rebase(base: string): void {
    if (base === '' || base === this._base.peek())
      return;
    this._base.value = base;
  }

  /** A trailing segment (`…/<table>/<keyOrId>`) or `?entity=<id>` → the entity page, `?q=` and
   * `?search=` → the list under them, nothing → the bare list; through the gate, so a Back with
   * unsaved changes asks first. A path is authoritative — `open('')` is the list; with no argument
   * the address bar is read instead, which is how a cold deep link reaches the app (the router
   * hands an app func only the path under its root, `routing.dart` ~:351, never the query). This
   * is the restore side of the URL, so no move it makes costs a history entry. Asked again for the
   * address it is restoring — a Back reaches both the app func and the view's own Back handler —
   * it joins that move: one gate. */
  async open(path?: string): Promise<boolean> {
    const address = path ?? `${location.pathname}${location.search}`;
    if (this._opening?.address !== address) {
      const opening = {address, done: this._sync.restore(() => this._open(address))};
      const clear = () => {
        if (this._opening === opening)
          this._opening = undefined;
      };
      this._opening = opening;
      void opening.done.then(clear, clear);
    }
    return this._opening!.done;
  }

  /** Moves between the pages through the gate: resolves to whether the move was made — a cancel
   * keeps the page and its unsaved changes. */
  async goTo(page: 'list'): Promise<boolean>;
  async goTo(page: 'entity', entity: string): Promise<boolean>;
  async goTo(page: DomainAppPage, entity: string | null = null): Promise<boolean> {
    return this._goTo(page, entity);
  }

  /** The create page — a pristine draft; `goTo('entity', 'new')` is the row keyed "new". */
  async create(): Promise<boolean> {
    return this._goTo('entity', DomainApp.NEW, null, true);
  }

  /** {@link DomainAddress.entityPath} over the app's base. */
  get entityPath(): boolean {
    return DomainAddress.entityPath(this.base);
  }

  /** How the row is addressed in a `/domains` path: its id where the table addresses rows by id
   * (`info.rowAddress`), else {@link DomainAddress.keyOf} over the business key — the raw key, which
   * {@link DomainAddress.segment} spells in the path. */
  keyOf(row: RowView): string {
    const info = this.table.info;
    return info.rowAddress === 'id' ? row.id : DomainAddress.keyOf(row, info.businessKey);
  }

  /** The ribbon the view's host places, in named slots: New, Save, Discard, the ⋯ menu and Refresh in
   * `main`, the search box in `tools`, the `ribbon` option's query switch in `presets`. Built once,
   * owned by the app. */
  ribbon(): DomainRibbon {
    if (this._ribbon !== undefined)
      return this._ribbon;
    const singular = nounOf(this.table.info, true);
    const add = new Control(this.build(() =>
      iconButton('plus', () => void this.create(), {tooltip: `New ${singular}`})));
    add.root.dataset.u2 = 'new-button';
    // in the trash the slot stays, so SAVE and DISCARD do not move
    add.effect(() => {
      const trash = this.trash.value;
      add.root.hidden = !trash && !this.listSource.access.value.can('insert');
      add.root.classList.toggle('u2-domain-app-invisible', trash);
    });
    const save = domains.saveButton(this.session);
    const pristine = () => !this.session.isDirty.value && this._page.value === 'entity' &&
      DomainApp._isDraft(this._entitySource.value, this._entitySource.value?.currentRow.value ?? null);
    const discard = domains.discardButton(this.session, {leave: {when: pristine, run: () => void this.goTo('list')}});
    // nothing on a read-only page is saved — the trash until a restore is staged
    save.effect(() => save.root.hidden = discard.root.hidden = !this.session.isDirty.value &&
      ((this._page.value === 'entity' ? this._entitySource.value : this.listSource)?.readOnly.value ?? false));
    const more = this._actionsMenu();
    // the status bar is text: the way to act on "Data changed" is here, and only while it stands
    const reload = new Control(button('Refresh', () => void this.refresh()));
    reload.root.dataset.u2 = 'refresh-button';
    reload.effect(() => reload.root.hidden = !this.stale.value);
    const search = domains.search(this.listSource);
    const sort = this._sortMenu();
    const switcher = this._viewSwitch();
    const tools: Control[] = [search, sort, switcher];
    if (this.table.table instanceof DgDomainTable)
      tools.push(this._propertiesButton());
    for (const control of [add, save, discard, more, reload, ...tools])
      this.own(() => control.dispose());
    this._listOnly(search);
    this._listOnly(sort, () => this.trash.value);
    this._listOnly(switcher);
    const base: DomainRibbon = {main: [add, save, discard, more, reload], tools, presets: []};
    return this._ribbon = this._options.ribbon?.(base, this) ?? base;
  }

  // the table takes the list's place, rebuilt when its columns change; the aggregations' values repaint it in place
  private _tableView(): void {
    let shown: DomainDataTable | null = null;
    this.effect(() => {
      const table = this.view.value === 'table';
      const columns = this.output.columns.value;
      const extraColumns = this.output.aggregations.value.map((a) => ({name: outputKey(a),
        header: outputKey(a), align: 'right' as const, render: (row: RowView) => {
          const value = this.output.values.peek().get(String(row.id))?.[outputKey(a)];
          return typeof value === 'number' ? DomainDataTable.number(value) : text(value);
        }}));
      this.list.root.hidden = table;
      this.output.root.hidden = !table;
      // the pane's filters stop taking its height, so Output sits right under them
      this.filters.root.style.flex = table ? 'none' : '';
      if (!table)
        return;
      this.output.show();
      const view = shown = new DomainDataTable(this.listSource, {columns, extraColumns,
        rowHeight: this._options.grid?.rowHeight, cells: this._options.grid?.cells});
      this.list.root.after(view.root);
      return () => {
        view.root.remove();
        shown = null;
      };
    });
    this.effect(() => {
      this.output.values.value;
      shown?.table.peek()?.refresh();
    });
  }

  private _viewSwitch(): Control {
    const group = this.build(() => new ButtonGroup({toggle: 'single', iconOnly: true, density: 'ribbon', items: [
      {id: 'list', icon: 'list', tooltip: 'Show as list', onClick: () => this.view.value = 'list'},
      {id: 'table', icon: 'table', tooltip: 'Show as table', onClick: () => this.view.value = 'table'}]}));
    group.effect(() => group.selected.value = [this.view.value]);
    return group;
  }

  private _propertiesButton(): Control {
    const [schema, table] = this.table.address.split('.');
    const button = new Control(this.build(() => iconButton('pencil-ruler',
      () => void domains.editProperties(schema, table).catch((e) => DomainErrors.report(e)), {tooltip: 'Edit properties…'})));
    button.root.hidden = true;
    void PropertyDialog.extensible(schema).then((tables) => button.root.hidden = !tables.has(table));
    return button;
  }

  /** The ribbon as `ViewHost.create` takes it — `[main, tools, presets]`, empty slots dropped.
   * @deprecated build or override {@link ribbon}; this is the shape the platform wants. */
  ribbonGroups(): (Control | HTMLElement)[][] {
    const {main, tools, presets} = this.ribbon();
    return [main, tools, presets].filter((group) => group.length > 0);
  }

  /** The Actions menu: the presets on a {@link compact} ribbon, the entity page's row actions, import
   * and bulk edit, Show deleted — then as the `menu` option changes it. Nothing in the trash: the
   * breadcrumb is the way out. */
  menuActions(): Action[] {
    const trash = this.trash.peek();
    // what the caller may do (`requires`) and what the table can do at all are different
    // questions: a permission hides an action, an undeclared support leaves it out entirely
    const support = this.table.table.support;
    const table: Action[] = [];
    if (!trash && support.writes && !this.listSource.readOnly.value) {
      table.push({text: 'Import…', icon: 'upload', requires: 'insert', run: () =>
        void domains.import(this.table).then((report) => report === null ? null : this.listSource.refresh())});
    }
    if (!trash && support.updateWhere) {
      table.push({text: 'Bulk edit…', icon: 'edit', requires: 'edit',
        run: () => void domains.bulkEdit(this.listSource)});
    }
    const mode: Action[] = trash || !support.deleted || !support.restore ? [] : [{
      text: 'Show deleted', icon: 'history', requires: 'delete', run: () => void this.setTrash(true)}];
    const presets = this.compact.peek() && this._page.peek() === 'list' && !trash ? this._presets : [];
    const base = [presets, this._rowActions(), table, mode].filter((section) => section.length > 0)
      .flatMap((section, i) => i === 0 ? section : [{...section[0], separator: true}, ...section.slice(1)]);
    return this._options.menu?.(base, this) ?? base;
  }

  /** Moves between the live collection and the trash — ONE switch: the list's `deleted` mode, its
   * sort, the breadcrumb and the `?trash=1` in every path all read it. Through the gate, so pending
   * changes are never dropped silently; resolves to whether the move was made. */
  async setMode(mode: DomainAppMode): Promise<boolean> {
    if (mode === this._mode.peek())
      return true;
    const on = mode === 'trash';
    if (!await confirmDiscard(this.session, {action: on ? 'open the trash' : 'leave the trash'}))
      return false;
    // one move, one history entry: the page and the mode settle together
    batch(() => {
      // the trash is a list view: an entity page open over a live row has nothing to show in it
      if (on && this._page.peek() === 'entity')
        this._show('list', null);
      this._mode.value = mode;
    });
    return true;
  }

  /** {@link setMode} as the ⋯ menu and an app's own code ask for it. */
  async setTrash(on: boolean): Promise<boolean> {
    return this.setMode(on ? 'trash' : 'live');
  }

  /** Reads the page's collection again — what the status bar's "Refresh" stands for while the
   * page is {@link stale}. Through the gate: a reload drops whatever is pending. */
  async refresh(): Promise<boolean> {
    const source = this._page.peek() === 'entity' ? this._entitySource.peek() : this.listSource;
    if (source === null || !await confirmDiscard(this.session, {action: 'reload the rows'}))
      return false;
    await source.refresh();
    return true;
  }

  /** A query switch for the ribbon — one button per preset, the pressed one writing the list's
   * query (`$me` is the current user's id); the one matching the query in force is shown pressed.
   * A press goes through the gate, like the filters. */
  presets(...entries: [label: string, query: string][]): ButtonGroup {
    const source = this.listSource;
    const bound = (at: number) => this._bound(entries[at][1]);
    const apply = (at: number) => {
      const q = source.query.peek();
      // `ButtonGroup`'s `single` is a radio group and never deselects: a press on the preset in
      // force clears the query here instead, and the effect below drops the pressed look
      const tree = Filters.equals(bound(at), DomainFilters.treeOf(q, source.schema)) ?
        Filters.group('and') : bound(at);
      source.query.value = typeof q === 'string' ? Filters.format(tree) : tree;
    };
    const press = (at: number) =>
      void confirmDiscard(this.session, {action: 'change the filter'}).then((ok) => ok && apply(at));
    const group = this.build(() => new ButtonGroup({toggle: 'single',
      items: entries.map(([label], at) => ({id: String(at), text: label, onClick: () => press(at)}))}));
    group.root.dataset.u2 = 'domain-presets';
    group.root.setAttribute('aria-label', 'Show');
    group.root.prepend(span('Show:', 'u2-domain-presets-label'));
    this._presets.push(...entries.map(([label], at): Action => ({text: label, run: () => press(at),
      get check() {
        return group.selected.peek().includes(String(at));
      }})));
    group.effect(() => {
      // the schema arrives with the load: read the state so the presets are bound again once it has
      source.state.value;
      const current = DomainFilters.treeOf(source.query.value, source.schema);
      // the box shows a preset as it is written ("assignee = $me"), not the id the query carries
      for (let i = 0; i < entries.length; i++)
        DomainFilters.showAs(source, Filters.format(bound(i)), entries[i][1]);
      const at = entries.findIndex((_, i) => Filters.equals(bound(i), current));
      group.selected.value = at < 0 ? [] : [String(at)];
    });
    this._listOnly(group, () => this.compact.value);
    return group;
  }

  getWidgetStatus(): IWidgetStatus {
    return {...super.getWidgetStatus(), shortcuts: {...SAVE_SHORTCUTS, ...this.shortcuts}};
  }

  /** The status as one line — what a plain status bar and the tests read. */
  static summaryOf(status: DomainAppStatus): string {
    const head = status.selected === 0 ? status.count : `${status.count} · ${status.selected} selected`;
    return [head, status.problem ?? status.notice].filter((part) => part !== null && part !== '').join(' — ');
  }

  /** One status-bar panel per slot (the view's host places them side by side): the count, the selection,
   * the notice, the problem — so a refusal never takes the count off the line. */
  static statusPanels(app: DomainApp): HTMLElement[] {
    const reads: ((status: DomainAppStatus) => string)[] = [
      (status) => status.count,
      (status) => status.selected === 0 ? '' : `${status.selected} selected`,
      (status) => status.notice ?? '',
      (status) => app._shortProblem(status),
    ];
    return ['count', 'selected', 'notice', 'problem'].map((name, at) => {
      const el = span('', `u2-domain-status-${name}`);
      el.dataset.u2Part = name;
      app.effect(() => {
        const status = app.status.value;
        const text = reads[at](status);
        el.textContent = text;
        el.hidden = text === '';
        el.classList.toggle('u2-domain-status-after', reads.slice(0, at).some((read) => read(status) !== ''));
        el.title = name === 'problem' && text !== (status.problem ?? '') ? status.problem ?? '' : '';
      });
      return el;
    });
  }

  // the list page shows the problem in full: here it is a few words (the reason is the tooltip), or none
  private _shortProblem(status: DomainAppStatus): string {
    const problem = status.problem ?? '';
    if (problem === '' || this._page.value !== 'list')
      return problem;
    if (problem === this._filterProblem.value)
      return '';
    if (this.listSource.state.value !== 'error')
      return problem;
    const q = this.listSource.query.value;
    return (typeof q === 'string' ? q.trim() !== '' : q.nodes.length > 0) ? 'Filter refused' : 'Rows not loaded';
  }

  private _pendingSummary(): string {
    const dirty = this.session.sources.value.filter((s) => s.isDirty.value);
    // a source's summary leads with its refusal, which is the `problem` slot already
    return dirty.length === 1 && dirty[0].error.value === undefined ?
      dirty[0].summary.value : this.session.summary.value;
  }

  private _caption(source: DomainSource, row: RowView): string {
    const shown = this.table.renderer.caption(source.current(row));
    // an emptied name leaves a saved row nameless (a blank, or the bare id): it keeps the name the row holds
    return shown === '' || shown === row.id ? this.table.renderer.caption(row) : shown;
  }

  /** How a table is named in a view title, a breadcrumb root or a child tab: the display name
   * the schema declares, else its plural name made presentable ('order_line' → 'Order lines'). */
  static titleOf(info: Pick<DomainTableInfoLike, 'pluralName' | 'friendlyName'>): string {
    return info.friendlyName ?? capital((info.pluralName || 'rows').replace(/_/g, ' '));
  }

  /** Plain-key shortcuts stay out of text entry; chords with Ctrl/Alt still fire there. */
  static isEditable(target: EventTarget | null): boolean {
    return target instanceof HTMLElement &&
      (target.isContentEditable || target instanceof HTMLInputElement || target instanceof HTMLTextAreaElement ||
        target instanceof HTMLSelectElement);
  }

  /** The key a keyboard event stands for, in the form {@link shortcuts} is keyed by. */
  static keyOf(e: KeyboardEvent): string {
    const key = e.key.length === 1 ? e.key.toUpperCase() : e.key;
    return [e.ctrlKey || e.metaKey ? 'Ctrl' : '', e.altKey ? 'Alt' : '', e.shiftKey ? 'Shift' : '', key]
      .filter((part) => part !== '').join('+');
  }

  // a `DomainForm` inside the app sees the event first and marks it handled, so the two never save twice
  private _save(e: KeyboardEvent): void {
    if (e.defaultPrevented || !isSaveKey(e))
      return;
    if (!this.root.contains(e.target as Node) && !this._isCurrentView())
      return;
    e.preventDefault();
    if (this.session.isDirty.peek())
      void this.session.save();
  }

  // the chrome the view's host places outside the app root is the app's own only while the shell shows its view
  private _isCurrentView(): boolean {
    const view = DomainApp._views.get(this);
    return view !== undefined && (grok.shell.v as DG.ViewBase | undefined)?.dart === view.dart;
  }

  private _shortcut(e: KeyboardEvent): void {
    if (!e.ctrlKey && !e.metaKey && !e.altKey && DomainApp.isEditable(e.target))
      return;
    const pressed = DomainApp.keyOf(e).toLowerCase();
    const name = Object.entries(this.shortcuts).find(([key]) => key.toLowerCase() === pressed)?.[1];
    const source = this._page.peek() === 'entity' ? this._entitySource.peek() : this.listSource;
    const row = source?.currentRow.peek() ?? null;
    if (name === undefined || source === null || row === null)
      return;
    const action = allowedActions(this.table.actions.for(row), {access: source.access.peek(), row})
      .find((a) => a.text === name);
    if (action === undefined || Reason.now(action.disabled) !== null)
      return;
    e.preventDefault();
    action.run();
  }

  // the menu is built on every open, so its items read the page and the mode in force
  private _actionsMenu(): Control {
    const control = this.build(() => {
      const el = iconButton('ellipsis-v', () => {}, {tooltip: 'More actions',
        menu: () => actionsMenu(this.menuActions(), {access: this.listSource.access.peek()})});
      el.dataset.u2 = 'actions-menu';
      return new Control(el);
    });
    control.effect(() => {
      this._page.value;
      this.compact.value;
      this.trash.value;
      this._entitySource.value?.currentRow.value;
      control.root.hidden = allowedActions(this.menuActions(), {access: this.listSource.access.value}).length === 0;
    });
    return control;
  }

  // through the gate: the reload a new order makes drops what is pending
  private _sortMenu(): Control {
    const source = this.listSource;
    const write = (sort: string) => void confirmDiscard(this.session, {action: 'change the order'})
      .then((ok) => ok && (source.sort.value = sort));
    const control = this.build(() => {
      const el = iconButton('sort-amount-down', () => {}, {tooltip: 'Sort by', menu: () => {
        const current = source.sort.peek().split(',')[0] ?? '';
        const desc = current.startsWith('!');
        const column = desc ? current.slice(1) : current;
        const menu = new Menu();
        for (const prop of source.schema.properties.filter((p) => p.ref === undefined && p.name !== 'id')) {
          menu.item(DomainTable.captionOf(prop), () => write(`${desc ? '!' : ''}${prop.name}`),
            {check: prop.name === column});
        }
        menu.separator();
        const unsorted = column === '' ? 'Pick a column to sort by' : null;
        menu.item('Ascending', () => write(column), {check: column !== '' && !desc, disabled: unsorted});
        menu.item('Descending', () => write(`!${column}`), {check: desc, disabled: unsorted});
        return menu;
      }});
      el.dataset.u2 = 'sort-menu';
      return new Control(el);
    });
    return control;
  }

  private _rowActions(): Action[] {
    const source = this._page.peek() === 'entity' ? this._entitySource.peek() : null;
    const row = source?.currentRow.peek() ?? null;
    if (source === null || row === null || Rows.isDraft(row))
      return [];
    return allowedActions([...this.table.actions.for(row), DomainList.deleteAction(() => {
      this._deleting = source;
      source.edit.peek()?.markDeleted(row.id);
    })], {access: source.access.peek(), row});
  }

  private _wireMode(): void {
    let saved: string | null = null;
    this.effect(() => {
      const on = this.mode.value === 'trash';
      this.listSource.deleted.value = on ? 'only' : 'exclude';
      // what a trash list is read for is what went in last; leaving gives the user's order back
      if (on) {
        saved ??= this.listSource.sort.peek();
        this.listSource.sort.value = TRASH_SORT;
      }
      else if (saved !== null) {
        this.listSource.sort.value = saved;
        saved = null;
      }
      this.root.classList.toggle('u2-domain-app-trash', on);
    });
  }

  private _withTrash(path: string): string {
    return this.trash.value ? `${path}${path.includes('?') ? '&' : '?'}trash=1` : path;
  }

  private _listOnly(control: Control, folded: () => boolean = () => false): void {
    this.effect(() => control.root.hidden = this._page.value !== 'list' || folded());
  }

  // parsed against the loaded schema, so a column-typed value round-trips the way the filters' own do
  private _bound(query: string) {
    const schema = this.listSource.schema;
    return Filters.bind(Filters.parse(query, schema).root, {me: grok.shell.user.id}, schema);
  }

  /** The browser's own gate: `beforeunload` is armed while anything is pending — a batch in
   * flight included, whose rows read clean until its write-back lands. */
  guardUnload(): void {
    if (this._unloadGuard !== undefined)
      return;
    const guard = (e: BeforeUnloadEvent) => {
      e.preventDefault();
      e.returnValue = '';
    };
    this._unloadGuard = guard;
    this.effect(() => {
      if (this.session.isDirty.value || this.session.isSaving.value)
        window.addEventListener('beforeunload', guard);
      else
        window.removeEventListener('beforeunload', guard);
    });
    this.own(() => window.removeEventListener('beforeunload', guard));
  }

  /** The app behind a view — filled by {@link register}. */
  static of(view: DG.ViewBase): DomainApp | undefined {
    return DomainApp._byView.get(view);
  }

  static viewOf(app: DomainApp): DG.ViewBase | undefined {
    return DomainApp._views.get(app);
  }

  /** The base path of the app open over a table address: the one the shell shows, else the last. */
  static baseOf(address: string): string | undefined {
    return DomainApp._pick((app) => app.table.address === address)?.base;
  }

  /** Ties an app to the view showing it (`DomainTable.app()`): {@link of} answers for it and
   * {@link activate} brings the view to the front, until the unregister runs. */
  static register(app: DomainApp, view: DG.ViewBase): () => void {
    DomainApp._byView.set(view, app);
    DomainApp._views.set(app, view);
    return () => {
      DomainApp._byView.delete(view);
      DomainApp._views.delete(app);
    };
  }

  /** Find-or-activate: the live app under `base` shows the entity and its view comes to the
   * front; false when no app is open there. */
  static activate(base: string | undefined, entity: string): boolean {
    const app = DomainApp._at(base);
    if (app === undefined)
      return false;
    void app.goTo('entity', entity);
    DomainApp._front(app);
    return true;
  }

  /** What `domains.route` reaches an app that is already open with: the whole address under its
   * base — the row segment, `?q=`, `?search=` — and its view comes to the front. Answers THAT
   * view, which the resolver hands back to the platform: an address the platform is told nothing
   * was found for falls through to the Dart domain view. */
  static openAt(base: string | undefined, address: string): DG.ViewBase | undefined {
    const app = DomainApp._at(base);
    if (app === undefined)
      return undefined;
    void app.open(address);
    DomainApp._front(app);
    return DomainApp._views.get(app);
  }

  private static _at(base: string | undefined): DomainApp | undefined {
    return base === undefined ? undefined : DomainApp._pick((app) => app.base === base);
  }

  private static _pick(test: (app: DomainApp) => boolean): DomainApp | undefined {
    const apps = [...DomainApp.live].filter(test);
    return apps.find((app) => app._isCurrentView()) ?? apps[apps.length - 1];
  }

  private static _front(app: DomainApp): void {
    const view = DomainApp._views.get(app);
    if (view !== undefined)
      grok.shell.v = view;
  }

  private async _open(path?: string): Promise<boolean> {
    const address = path ?? `${location.pathname}${location.search}`;
    const at = address.indexOf('?');
    const params = new URLSearchParams(at < 0 ? '' : address.slice(at + 1));
    // trash mode is the list's own view of the table, and a path is authoritative about it
    this._mode.value = params.get('trash') === '1' ? 'trash' : 'live';
    const segment = DomainAddress.segmentOf(at < 0 ? address : address.slice(0, at), this._bases);
    // the draft sentinel is read off the raw spelling, before the escape comes off: `~new` is the
    // row keyed "new" in a segment and in `?entity=` alike
    const raw = segment ?? params.get('entity');
    const draft = raw === DomainApp.NEW;
    const entity = raw === null ? null : DomainAddress.unescape(raw);
    if (entity !== null && entity !== '') {
      // a table addressing rows by id never has its segment split: the segment IS the id
      const info = this.table.info;
      return this._goTo('entity', entity, segment === null || draft || info.rowAddress === 'id' ? null :
        DomainAddress.keyQuery(entity, info.businessKey, this.listSource.schema), draft);
    }
    if (!await this._goTo('list', null))
      return false;
    const q = params.get('q') ?? '';
    const search = params.get('search') ?? '';
    const trash = this._mode.peek() === 'trash';
    const sort = trash ? this.listSource.sort.peek() : params.get('sort') ?? this._options.sort ?? '';
    const current = this.listSource.query.peek();
    // already on the list, the gate `goTo` skipped: the query is changing under pending changes
    if ((typeof current === 'string' ? current : Filters.format(current)) === q &&
        this.listSource.search.peek() === search && this.listSource.sort.peek() === sort)
      return true;
    if (!await confirmDiscard(this.session, {action: 'change the filter'}))
      return false;
    batch(() => {
      this.listSource.query.value = q;
      this.listSource.search.value = search;
      this.listSource.sort.value = sort;
    });
    return true;
  }

  // `draft` is the create page: an entity of "new" without it is the row so keyed
  private async _goTo(page: DomainAppPage, entity: string | null,
    query: FilterGroup | null = null, draft = false): Promise<boolean> {
    if (page === this._page.peek() && (page === 'list' || entity === this._entity.peek()))
      return true;
    if (!await confirmDiscard(this.session, {action: 'leave this page'}) || this.scope.isDisposed)
      return false;
    this._show(page, entity, query, draft);
    return true;
  }

  private _entityAddress(entity: string): string {
    const row = this._entitySource.value?.currentRow.value ?? null;
    if (row === null || Rows.isDraft(row))
      return entity;
    return this.entityPath ? this.keyOf(row) : row.id;
  }

  // both routes: a `/domains/…` link still reaches an app that has rebased onto `/apps/…`
  private get _bases(): string[] {
    return [this.base, this._options.base];
  }

  private _state(): ViewState {
    const q = this.listSource.query.peek();
    const entity = this._entity.peek();
    const source = this._entitySource.peek();
    // the row keyed "new" is not the draft: spelled as its segment, it stays a move of its own
    const draft = DomainApp._isDraft(source, source?.currentRow.peek() ?? null);
    return {page: this._page.peek(), entity: entity === null || draft ? entity : DomainAddress.escape(entity),
      search: this.listSource.search.peek(), query: typeof q === 'string' ? q : Filters.format(q),
      mode: this._mode.peek()};
  }

  private static _isDraft(source: DomainSource | null, row: RowView | null): boolean {
    return source !== null && source.isDraft && (row === null || Rows.isDraft(row));
  }

  private _show(page: DomainAppPage, entity: string | null, query: FilterGroup | null = null,
    draft = false): void {
    this._close();
    if (page === 'entity' && entity !== null) {
      const source = this.build(() => draft ? this._source({draft: true, defaults: this._options.defaults}) :
        this._source({query: query ?? Filters.group('and', [Filters.cond('id', '=', entity)]), pageSize: 1}));
      // a source of one has no list to make its row current
      source.effect(() => {
        const rows = source.rows.items.value;
        if (source.state.value === 'ready' && rows.length > 0 && source.currentRow.peek() === null)
          source.currentRow.value = rows[0];
      });
      const {children, history, include, exclude} = this._options;
      const form = this.build(() => new DomainForm(source, {system: 'footer', include, exclude,
        labelColumn: this.labelColumn, empty: `${capital(nounOf(this.table.info, true))} "${entity}" was not found.`}));
      this._panes = this.build(() => [
        ...(children === false ? [] : [domains.children(source, children === true ? undefined : children)]),
        ...(history === false || !this.table.table.support.audit ? [] :
          [domains.history(source, undefined, {include, exclude})]),
      ]);
      this._formHost.replaceChildren(form.root);
      this.panes.replaceChildren(...this._panes.map((pane) => pane.root));
      batch(() => {
        this._entitySource.value = source;
        this._form.value = form;
      });
    }
    batch(() => {
      this._entity.value = page === 'entity' ? entity : null;
      this._page.value = page;
    });
  }

  private _close(): void {
    const form = this._form.peek();
    const source = this._entitySource.peek();
    batch(() => {
      this._form.value = null;
      this._entitySource.value = null;
    });
    for (const pane of this._panes)
      pane.dispose();
    this._panes = [];
    form?.dispose();
    source?.dispose();
    this._formHost.replaceChildren();
    this.panes.replaceChildren();
  }

  private _source(options: Parameters<DomainTable['source']>[0]): DomainSource {
    return SharedSession.runWith(this.session, () => this.table.source(options));
  }
}
