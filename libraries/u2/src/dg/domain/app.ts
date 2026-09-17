/* The app shape over one table (GOAL "What it looks like"): a list page and an entity page under
   one `SharedSession`, the URL following the page (`?q=` the list's query, `?entity=` the row),
   `open(path)` the way back from the address bar, and the unsaved gate in front of every move
   between them. `DomainTable.app()` (WO 2-9) wraps it in an `appView`; the pieces are public so a
   hand-built page composes the same ribbon and status. Made to be subclassed (`app({app: Sub})`):
   `ribbon()` is overridden, `presets()` builds a query switch for it, `shortcuts` maps keys onto
   the table's actions. */
import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {Control} from '../../core/component.js';
import {signal, computed, batch, ReadonlySignal, Signal} from '../../core/signals.js';
import {button, div, divV, span} from '../../core/elements.js';
import {Filters} from '../../core/filter/index.js';
import type {FilterGroup} from '../../core/filter/index.js';
import type {IWidgetStatus} from '../../core/widget-like.js';
import {allowedActions, actionsMenu} from '../../components/actions/actions.js';
import type {Action} from '../../components/actions/actions.js';
import {iconButton} from '../../components/actions/buttons.js';
import {ButtonGroup} from '../../components/actions/button-group.js';
import {Breadcrumbs} from '../../components/navigation/breadcrumbs.js';
import type {DomainSource} from '../../sources/domain-source.js';
import type {DomainTableInfoLike} from '../../sources/domain-backend.js';
import {SharedSession, confirmDiscard} from '../../sources/session.js';
import {Rows} from '../../sources/rows-like.js';
import type {RowView} from '../../sources/rows-like.js';
import {domains, DomainTable} from './index.js';
import {DomainList} from './list.js';
import type {DomainListMode} from './list.js';
import {SAVE_SHORTCUTS, isSaveKey, DomainForm} from './form.js';
import {DomainFilters} from './filters.js';
import type {DomainChildrenOptions} from './children.js';
import {DomainAddress} from './address.js';
import {ViewSync} from './view-sync.js';
import type {ViewState} from './view-sync.js';

export type DomainAppPage = 'list' | 'entity';

export type DomainAppMode = 'live' | 'trash';

/** What the status bar says while a live source is behind the server ({@link DomainApp.refresh}). */
const STALE = 'Data changed — Refresh';

/** The trash reads newest-deleted first: a soft delete is a write, so it stamps `updated_on`. */
const TRASH_SORT = '!updated_on';

/** What the app's status bar says, as VALUES rather than one composed sentence: a phase-4 view
 * spec renders the slots it wants instead of parsing text, and a refusal no longer takes the row
 * count off the line (U37). */
export interface DomainAppStatus {
  /** What the page holds: "50 of 2,077", "3 issues", "1 deleted issue", or the row's caption on
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

/** The ribbon as slots rather than positions: a subclass adds to the slot it means, and a phase-4
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
  pageSize?: number;
  /** Whether the list follows the server (`DomainSourceOptions.live`), ON by default: it probes
   * the table every {@link liveMs} and reloads while the session is clean, marking the page
   * {@link DomainApp.stale} while it is not. The entity page is one row and is not polled. */
  live?: boolean;
  liveMs?: number;
  mode?: DomainListMode;
  /** What every draft the app creates starts with — the preset a view that opened the app carries
   * (a location's id on a containers app). Threaded into the entity page's draft source. */
  defaults?: Record<string, unknown>;
  /** The form's columns, in this order; every column by default. */
  include?: string[];
  /** The entity page's panes under the form, on by default: the child collections (a tab per
   * referring table; options narrow them) and the row's history. */
  children?: boolean | DomainChildrenOptions;
  history?: boolean;
  /** Key → the name of an action in `table.actions`, run over the current row while the app has
   * the focus: `{'Ctrl+Shift+C': 'Close'}` (modifiers in the order Ctrl, Alt, Shift; Cmd counts
   * as Ctrl). */
  shortcuts?: Record<string, string>;
}

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
  /** See {@link DomainAppOptions.shortcuts}; a subclass may declare its own. */
  shortcuts: Record<string, string>;

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
  private readonly _sync: ViewSync;

  private static readonly _byView = new WeakMap<DG.ViewBase, DomainApp>();
  private static readonly _views = new WeakMap<DomainApp, DG.ViewBase>();
  private static readonly _bases = new Map<string, string>();

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
    // the chrome `appView` places sits outside the app root, so Save has to be caught above both
    const onSave = (e: KeyboardEvent) => this._save(e);
    document.addEventListener('keydown', onSave);
    this.own(() => document.removeEventListener('keydown', onSave));
    DomainApp.live.add(this);
    DomainApp._bases.set(table.address, this.base);
    this.own(() => {
      DomainApp.live.delete(this);
      if (DomainApp._bases.get(table.address) === this.base)
        DomainApp._bases.delete(table.address);
    });

    this.listSource = this._source({query: _options.query ?? '', pageSize: _options.pageSize,
      live: _options.live ?? true, liveMs: _options.liveMs});
    this.own(() => this.listSource.dispose());
    this.list = this.runInScope(() => new DomainList(this.listSource, {mode: _options.mode}));
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
    this.breadcrumbs = this.runInScope(() => new Breadcrumbs({items: [title], onClick: (index) => {
      // in the trash the root crumb reads "Issues › Trash" and is the way back out of it; the ⋯
      // menu's "Exit trash" must not be the only one
      if (index === 0)
        void (this.trash.peek() ? this.setTrash(false) : this.goTo('list'));
    }}));
    this._formHost = div([], 'u2-domain-app-form');
    this.panes = div([], 'u2-domain-app-panes');
    this.panes.dataset.u2Part = 'panes';
    this.filters = this.runInScope(() => domains.filters(this.listSource));
    // the rows answer the filter that was applied, not the one being written: while the box is
    // invalid they are shown as stale and the count stands down
    this.effect(() => {
      const problem = this.filters.problem.value;
      this._filterProblem.value = problem;
      this.list.stale.value = problem !== null;
    });
    const listPage = divV([this.filters.root, this.list.root], 'u2-domain-app-list');
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
      const caption = row === null ? '…' : table.renderer.caption(row);
      this.breadcrumbs.setItems([title, caption]);
    });
    // the entity page over a draft: once saved, the page is the row's
    this.effect(() => {
      const row = this._entitySource.value?.currentRow.value ?? null;
      if (row !== null && this._entity.peek() === DomainApp.NEW && !Rows.isDraft(row))
        this._entity.value = row.id;
    });
    // and once the draft is discarded — by the Discard button or by a gate's DISCARD — the page
    // has no row left, so the app goes back to the list. `had` keeps the load's own empty window
    // (the frame arrives before `newRow` fills it) from counting as a discard; the move is
    // deferred past the discard that is running, which reads dirty until it has settled.
    let had = false;
    this.effect(() => {
      const source = this._entitySource.value;
      if (source === null || !source.isDraft) {
        had = false;
        return;
      }
      const rows = source.rows.items.value.length;
      if (rows === 0 && had)
        queueMicrotask(() => void this.goTo('list'));
      had = rows > 0;
    });

    this.path = computed(() => {
      const entity = this._entity.value;
      if (this._page.value === 'entity' && entity !== null) {
        const address = this._entityAddress(entity);
        return this.entityPath ? `${this.base}/${encodeURIComponent(address)}` :
          `${this.base}?entity=${encodeURIComponent(address)}`;
      }
      const q = this.listSource.query.value;
      const search = this.listSource.search.value;
      const path = this._withTrash(Filters.queryPath(this.base, typeof q === 'string' ? q : Filters.format(q)));
      return search === '' ? path : `${path}${path.includes('?') ? '&' : '?'}search=${encodeURIComponent(search)}`;
    });
    this._sync = new ViewSync({
      path: this.path,
      draftEntity: DomainApp.NEW,
      state: () => this._state(),
      open: (address) => this.open(address),
      owns: (path) => DomainAddress.restOf(path, this._bases) !== null,
      isCurrentView: () => this._isCurrentView(),
    }, this);
    // What the page HOLDS, and nothing else: a refusal, the pending changes and the poll's news
    // are slots of their own, while the source folds all of them into its one summary line — so
    // the count is taken while the source is saying nothing else, and stands until it can be
    // taken again. A failed load holds nothing, and says nothing.
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
        this._count.value = this.table.renderer.caption(row);
      else if (source.state.value === 'error')
        this._count.value = '';
      else if (source.error.value === undefined && source.changeCount.value === 0)
        this._count.value = source.summary.value;
    });
    this.status = computed(() => {
      const entity = this._page.value === 'entity';
      const source = entity ? this._entitySource.value : this.listSource;
      const dirty = this.session.isDirty.value;
      const pending = dirty ? this._pendingSummary() : null;
      const count = this._count.value;
      // what a bulk action would act on, beside the count — the platform grid's own convention
      const selected = entity ? 0 : this.listSource.selection.value.length;
      const problem = source !== null && source.error.value !== undefined ? source.summary.value : null;
      // the rows answer the filter that was applied, not the one in the box: while that one is
      // invalid they are shown stale, and a refusal of the old one is not the news
      if (!entity && this.list.stale.value && !dirty)
        return {count, selected, notice: this._filterProblem.value ?? 'Filter not applied', problem};
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
  }

  /** Re-points the app at the route the shell actually mounted its view at (`/apps/Stockroom`):
   * every path the app produces moves with it, and find-or-activate looks the table up under the
   * new base. See `DomainTable.app()`, which resolves it from the docked view. */
  rebase(base: string): void {
    if (base === '' || base === this._base.peek())
      return;
    if (DomainApp._bases.get(this.table.address) === this._base.peek())
      DomainApp._bases.set(this.table.address, base);
    this._base.value = base;
  }

  /** A trailing segment (`…/<table>/<keyOrId>`) or `?entity=<id>` → the entity page, `?q=` and
   * `?search=` → the list under them, nothing → the bare list; through the gate, so a Back with
   * unsaved changes asks first. A path is authoritative — `open('')` is the list; with no argument
   * the address bar is read instead, which is how a cold deep link reaches the app (the router
   * hands an app func only the path under its root, `routing.dart` ~:351, never the query). This
   * is the restore side of the URL, so no move it makes costs a history entry. */
  async open(path?: string): Promise<boolean> {
    return this._sync.restore(() => this._open(path));
  }

  /** Moves between the pages through the gate: resolves to whether the move was made — a cancel
   * keeps the page and its unsaved changes. */
  async goTo(page: 'list'): Promise<boolean>;
  async goTo(page: 'entity', entity: string): Promise<boolean>;
  async goTo(page: DomainAppPage, entity: string | null = null): Promise<boolean> {
    return this._goTo(page, entity);
  }

  /** {@link DomainAddress.entityPath} over the app's base. */
  get entityPath(): boolean {
    return DomainAddress.entityPath(this.base);
  }

  /** {@link DomainAddress.keyOf} over the table's business key. */
  keyOf(row: RowView): string {
    return DomainAddress.keyOf(row, this.table.info.businessKey);
  }

  /** The ribbon `appView` places, in named slots: New, Save, Discard, the ⋯ menu and Refresh in
   * `main`, the search box in `tools`, a subclass's query switch in `presets`. Built once, owned
   * by the app. */
  ribbon(): DomainRibbon {
    if (this._ribbon !== undefined)
      return this._ribbon;
    const add = new Control(button('New', () => void this.goTo('entity', DomainApp.NEW)));
    add.root.dataset.u2 = 'new-button';
    add.effect(() => add.root.hidden = !this.listSource.access.value.can('insert'));
    const save = domains.saveButton(this.session);
    const discard = domains.discardButton(this.session);
    const more = this._actionsMenu();
    // the status bar is text: the way to act on "Data changed" is here, and only while it stands
    const reload = new Control(button('Refresh', () => void this.refresh()));
    reload.root.dataset.u2 = 'refresh-button';
    reload.effect(() => reload.root.hidden = !this.stale.value);
    const search = domains.search(this.listSource);
    for (const control of [add, save, discard, more, reload, search])
      this.own(() => control.dispose());
    this._listOnly(search);
    return this._ribbon = {main: [add, save, discard, more, reload], tools: [search], presets: []};
  }

  /** The ribbon as `appView` takes it — `[main, tools, presets]`, empty slots dropped.
   * @deprecated build or override {@link ribbon}; this is the shape the platform wants. */
  ribbonGroups(): (Control | HTMLElement)[][] {
    const {main, tools, presets} = this.ribbon();
    return [main, tools, presets].filter((group) => group.length > 0);
  }

  /** The table-wide actions of the ⋯ menu: the import wizard, the bulk edit over what the list
   * holds, and the trash. Nothing writes in the trash, so only the way out of it is offered there. */
  menuActions(): Action[] {
    const trash = this.trash.peek();
    // what the caller may do (`requires`) and what the table can do at all are different
    // questions: a permission hides an action, an undeclared support leaves it out entirely
    const support = this.table.table.support;
    const out: Action[] = [];
    if (!trash && support.writes) {
      out.push({name: 'Import…', icon: 'upload', requires: 'insert', run: () =>
        void domains.import(this.table).then((report) => report === null ? null : this.listSource.refresh())});
      out.push({name: 'Bulk edit…', icon: 'edit', requires: 'edit',
        run: () => void domains.bulkEdit(this.listSource)});
    }
    if (support.deleted && support.restore) {
      out.push({name: trash ? 'Exit trash' : 'Trash', icon: trash ? 'arrow-left' : 'trash-alt',
        requires: 'delete', run: () => void this.setTrash(!trash)});
    }
    return out;
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

  /** {@link setMode} as the ⋯ menu and a subclass ask for it. */
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
    const group = this.runInScope(() => new ButtonGroup({toggle: 'single', density: 'ribbon',
      items: entries.map(([label], at) => ({id: String(at), label, onClick: () =>
        void confirmDiscard(this.session, {action: 'change the filter'}).then((ok) => ok && apply(at))}))}));
    group.root.dataset.u2 = 'domain-presets';
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
    this._listOnly(group);
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

  /** One status-bar panel per slot (`appView` places them side by side): the count, the selection,
   * the notice, the problem — so a refusal never takes the count off the line. */
  static statusPanels(app: DomainApp): HTMLElement[] {
    const slot = (name: string, read: (status: DomainAppStatus) => string): HTMLElement => {
      const el = span('', `u2-domain-status-${name}`);
      el.dataset.u2Part = name;
      app.effect(() => {
        const text = read(app.status.value);
        el.textContent = text;
        el.hidden = text === '';
      });
      return el;
    };
    return [
      slot('count', (status) => status.count),
      slot('selected', (status) => status.selected === 0 ? '' : `${status.selected} selected`),
      slot('notice', (status) => status.notice ?? ''),
      slot('problem', (status) => status.problem ?? ''),
    ];
  }

  /** What the pending changes are called: one dirty source words them in its own terms ("1
   * restore pending", "2 deletions pending"), several can only be counted. */
  private _pendingSummary(): string {
    const dirty = this.session.sources.value.filter((s) => s.isDirty.value);
    // a source's summary leads with its refusal, which is the `problem` slot already
    return dirty.length === 1 && dirty[0].error.value === undefined ?
      dirty[0].summary.value : this.session.summary.value;
  }

  /** How a table is named in a view title, a breadcrumb root or a child tab: the display name
   * the schema declares, else its plural name made presentable ('order_line' → 'Order lines'). */
  static titleOf(info: DomainTableInfoLike): string {
    const plural = (info.pluralName || 'rows').replace(/_/g, ' ');
    return info.friendlyName ?? `${plural.charAt(0).toUpperCase()}${plural.slice(1)}`;
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

  /** Ctrl+S / Ctrl+Enter → the session's Save, from everywhere the app reaches: its pages and the
   * ribbon, the toolbox and the status bar `appView` placed outside them. A `DomainForm` inside
   * the app sees the event first and marks it handled, so the two never save twice. */
  private _save(e: KeyboardEvent): void {
    if (e.defaultPrevented || !isSaveKey(e))
      return;
    if (!this.root.contains(e.target as Node) && !this._isCurrentView())
      return;
    e.preventDefault();
    if (this.session.isDirty.peek())
      void this.session.save();
  }

  /** Whether the shell shows this app's view — what makes its chrome, outside the app root, the
   * app's own. A hand-built page that never registered a view answers false. */
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
      .find((a) => a.name === name);
    if (action === undefined || action.enabled === false)
      return;
    e.preventDefault();
    action.run();
  }

  /** The ⋯ button: the menu is built on every open, so its items read the mode in force.
   * Permission ⇒ hidden — the button is gone while the access allows nothing in it. */
  private _actionsMenu(): Control {
    const control = this.runInScope(() => {
      const el = iconButton('ellipsis-v', () => {}, {tooltip: 'More actions'});
      el.dataset.u2 = 'actions-menu';
      el.title = 'More actions';
      el.addEventListener('click', () =>
        actionsMenu(this.menuActions(), {access: this.listSource.access.peek()}).show({anchor: el}));
      return new Control(el);
    });
    control.effect(() => control.root.hidden =
      allowedActions(this.menuActions(), {access: this.listSource.access.value}).length === 0);
    return control;
  }

  /** Trash mode is one flip of the list source's `deleted` mode: the search box, the filters and
   * the list stay bound to the one source, and its access narrows itself to read-only. */
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

  private _listOnly(...controls: Control[]): void {
    this.effect(() => {
      const list = this._page.value === 'list';
      for (const control of controls)
        control.root.hidden = !list;
    });
  }

  /** A preset's query with `$me` bound, as a tree — parsed against the table's schema once it is
   * loaded, so a column-typed value round-trips the way the filters' own do. */
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

  /** The base path of the app last opened over a table address. */
  static baseOf(address: string): string | undefined {
    return DomainApp._bases.get(address);
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
   * was found for falls through to the Dart domain view, which is how Back off a `/domains`
   * address used to close the app it should have restored. */
  static openAt(base: string | undefined, address: string): DG.ViewBase | undefined {
    const app = DomainApp._at(base);
    if (app === undefined)
      return undefined;
    void app.open(address);
    DomainApp._front(app);
    return DomainApp._views.get(app);
  }

  private static _at(base: string | undefined): DomainApp | undefined {
    return base === undefined ? undefined : [...DomainApp.live].find((app) => app.base === base);
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
    const entity = segment ?? params.get('entity');
    if (entity !== null && entity !== '') {
      return this._goTo('entity', entity, segment === null ? null :
        DomainAddress.keyQuery(segment, this.table.info.businessKey, this.listSource.schema));
    }
    if (!await this._goTo('list', null))
      return false;
    const q = params.get('q') ?? '';
    const search = params.get('search') ?? '';
    const current = this.listSource.query.peek();
    // already on the list, the gate `goTo` skipped: the query is changing under pending changes
    if ((typeof current === 'string' ? current : Filters.format(current)) === q &&
        this.listSource.search.peek() === search)
      return true;
    if (!await confirmDiscard(this.session, {action: 'change the filter'}))
      return false;
    batch(() => {
      this.listSource.query.value = q;
      this.listSource.search.value = search;
    });
    return true;
  }

  /** `query` is how the row was addressed when it was not by its id — the business key an address
   * segment spelled. */
  private async _goTo(page: DomainAppPage, entity: string | null,
    query: FilterGroup | null = null): Promise<boolean> {
    if (page === this._page.peek() && (page === 'list' || entity === this._entity.peek()))
      return true;
    if (!await confirmDiscard(this.session, {action: 'leave this page'}))
      return false;
    this._show(page, entity, query);
    return true;
  }

  /** How the entity page addresses its row: the loaded row spells it — the business key where the
   * app addresses rows by a segment, the id otherwise — and until the row is in, the address the
   * page was opened by. */
  private _entityAddress(entity: string): string {
    const row = this._entitySource.value?.currentRow.value ?? null;
    if (row === null || Rows.isDraft(row))
      return entity;
    return this.entityPath ? this.keyOf(row) : row.id;
  }

  /** Both of the app's routes: the one the shell mounted the view at, and the one it was built
   * with — a `/domains/…` link still reaches an app that has rebased onto `/apps/…`. */
  private get _bases(): string[] {
    return [this.base, this._options.base];
  }

  private _state(): ViewState {
    const q = this.listSource.query.peek();
    return {page: this._page.peek(), entity: this._entity.peek(), search: this.listSource.search.peek(),
      query: typeof q === 'string' ? q : Filters.format(q), mode: this._mode.peek()};
  }

  private _show(page: DomainAppPage, entity: string | null, query: FilterGroup | null = null): void {
    this._close();
    if (page === 'entity' && entity !== null) {
      const draft = entity === DomainApp.NEW;
      const source = draft ? this._source({draft: true, defaults: this._options.defaults}) :
        this._source({query: query ?? Filters.group('and', [Filters.cond('id', '=', entity)]), pageSize: 1});
      // a source of one has no list to make its row current
      source.effect(() => {
        const rows = source.rows.items.value;
        if (source.state.value === 'ready' && rows.length > 0 && source.currentRow.peek() === null)
          source.currentRow.value = rows[0];
      });
      const {children, history, include} = this._options;
      const form = this.runInScope(() => new DomainForm(source, {system: 'footer', include,
        empty: `${this.table.info.singularName || 'Row'} "${entity}" was not found.`}));
      this._panes = this.runInScope(() => [
        ...(children === false ? [] : [domains.children(source, children === true ? undefined : children)]),
        ...(history === false ? [] : [domains.history(source)]),
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
