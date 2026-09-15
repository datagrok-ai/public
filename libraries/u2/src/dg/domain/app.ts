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
import {button, div, divV} from '../../core/elements.js';
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
import {DomainForm, SAVE_SHORTCUTS, isSaveKey} from './form.js';
import {DomainFilters} from './filters.js';
import type {DomainChildrenOptions} from './children.js';

export type DomainAppPage = 'list' | 'entity';

/** A row id as the platform spells one — what tells an id from a business key in an address. */
const ID = /^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i;

/** What a user move changes: a path that changes without it is a rewrite of the same page. */
interface AppState {page: DomainAppPage; entity: string | null; query: string; search: string}

export interface DomainAppOptions {
  table: DomainTable;
  /** The view path the app's paths extend, e.g. `/apps/Grit/Issues`. */
  base: string;
  /** The list's initial query. */
  query?: string;
  pageSize?: number;
  mode?: DomainListMode;
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
  /** The session's pending changes when there are any, else the page's source. */
  readonly summary: ReadonlySignal<string>;
  /** Trash mode: the list answers the table's deleted rows — read-only, Restore on each — and
   * every path it emits carries `?trash=1`. Written by the ⋯ menu ({@link setTrash}) and read
   * back from the address bar by {@link open}. */
  readonly trash: ReadonlySignal<boolean>;
  /** See {@link DomainAppOptions.shortcuts}; a subclass may declare its own. */
  shortcuts: Record<string, string>;

  private readonly _trash = signal(false);
  private readonly _base: Signal<string>;
  private readonly _page = signal<DomainAppPage>('list');
  private readonly _entity = signal<string | null>(null);
  private readonly _entitySource = signal<DomainSource | null>(null);
  private readonly _form = signal<DomainForm | null>(null);
  private readonly _formHost: HTMLElement;
  private _panes: Control[] = [];
  private _ribbon: (Control | HTMLElement)[][] | undefined;
  private _unloadGuard: ((e: BeforeUnloadEvent) => void) | undefined;
  /** Nonzero while a move is being restored from the address bar ({@link open}), which is where
   * `popstate` and `handlePath` arrive: the entry exists already, so nothing is pushed. */
  private _restoring = 0;

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
    this.trash = this._trash;
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
    // Back restores the app itself: the platform router re-parses the URL on `popstate` but never
    // asks a JS-hosted view (`handlePath` is not called on the stand — probed 2026-09-15), so the
    // address the entry carries would only move the URL. A router that does call `handlePath`
    // lands on the same `open`, which is idempotent.
    const onPop = () => this._restoreFromUrl();
    window.addEventListener('popstate', onPop);
    this.own(() => window.removeEventListener('popstate', onPop));
    DomainApp.live.add(this);
    DomainApp._bases.set(table.address, this.base);
    this.own(() => {
      DomainApp.live.delete(this);
      if (DomainApp._bases.get(table.address) === this.base)
        DomainApp._bases.delete(table.address);
    });

    this.listSource = this._source({query: _options.query ?? '', pageSize: _options.pageSize});
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
      if (index === 0)
        void this.goTo('list');
    }}));
    this._formHost = div([], 'u2-domain-app-form');
    this.panes = div([], 'u2-domain-app-panes');
    this.panes.dataset.u2Part = 'panes';
    const listPage = div([this.list.root], 'u2-domain-app-list');
    listPage.dataset.u2Part = 'list-page';
    const entityPage = divV([this.breadcrumbs.root, this._formHost, this.panes], 'u2-domain-app-entity');
    entityPage.dataset.u2Part = 'entity-page';
    this.root.append(listPage, entityPage);
    this.effect(() => {
      const page = this._page.value;
      listPage.hidden = page !== 'list';
      entityPage.hidden = page !== 'entity';
    });
    this.effect(() => {
      if (this._trash.value) {
        this.breadcrumbs.setItems([title, 'Trash']);
        return;
      }
      const row = this._entitySource.value?.currentRow.value ?? null;
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
    // One history entry per user move. The shell mirrors `view.path` onto the address bar with
    // `replaceState` (`routing.dart` onViewUrlChanged, a synchronous stream), so a move costs no
    // entry of its own — the app pushes it, and the shell's replace then rewrites that same entry.
    // The effect reads `path` and nothing else: a dependency on the state signals themselves would
    // put it behind `appView`'s mirror in the flush order, too late to push.
    let state = this._state();
    this.effect(() => {
      const path = this.path.value;
      const previous = state;
      state = this._state();
      if (this._restoring === 0 && DomainApp._moved(previous, state) && typeof history !== 'undefined')
        history.pushState(null, '', path);
    });
    this.summary = computed(() => {
      const entity = this._page.value === 'entity';
      const source = entity ? this._entitySource.value : this.listSource;
      // the rows and the count answer the filter that was applied, not the one in the box: while
      // that one is invalid the count stands down, and a refusal of the old one is not the news
      if (!entity && this.list.stale.value && !this.session.isDirty.value)
        return '—';
      // a refusal outranks the change count: the summary says why the pending changes are stuck
      if (source !== null && source.error.value !== undefined)
        return source.summary.value;
      if (this.session.isDirty.value)
        return this.session.summary.value;
      // the entity page shows one row, not a collection: it says which row, never "1 substance"
      const row = entity ? source?.currentRow.value ?? null : null;
      if (row !== null && !Rows.isDraft(row))
        return this.table.renderer.caption(row);
      return source?.summary.value ?? '';
    });
    this._wireTrash();
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
    this._restoring++;
    try {
      return await this._open(path);
    } finally {
      this._restoring--;
    }
  }

  /** Moves between the pages through the gate: resolves to whether the move was made — a cancel
   * keeps the page and its unsaved changes. */
  async goTo(page: 'list'): Promise<boolean>;
  async goTo(page: 'entity', entity: string): Promise<boolean>;
  async goTo(page: DomainAppPage, entity: string | null = null): Promise<boolean> {
    return this._goTo(page, entity);
  }

  /** Whether a row is addressed by a path segment (`${base}/${key}`) rather than by `?entity=`:
   * the platform's own `/domains/…` routes are, an app mounted at `/apps/…` is not (ruling R2). */
  get entityPath(): boolean {
    return this.base.startsWith('/domains/');
  }

  /** How a row is spelled in a `/domains` path: the business key when it is unambiguous, the id
   * otherwise — no key declared, a null component, or a composite key whose values carry a '-'
   * (the TS half of `domain_row_meta.dart` `deepLink`). */
  keyOf(row: RowView): string {
    const parts: string[] = [];
    for (const column of this.table.info.businessKey) {
      const value = row[column];
      if (value === null || value === undefined)
        return row.id;
      parts.push(String(value));
    }
    if (parts.length === 0 || (parts.length > 1 && parts.some((part) => part.includes('-'))))
      return row.id;
    return parts.join('-');
  }

  /** The ribbon `appView` places: New, Save, Discard and the ⋯ menu; the search box and the
   * filters on the list page. Built once, owned by the app. */
  ribbon(): (Control | HTMLElement)[][] {
    if (this._ribbon !== undefined)
      return this._ribbon;
    const add = new Control(button('New', () => void this.goTo('entity', DomainApp.NEW)));
    add.root.dataset.u2 = 'new-button';
    add.effect(() => add.root.hidden = !this.listSource.access.value.can('insert'));
    const save = domains.saveButton(this.session);
    const discard = domains.discardButton(this.session);
    const more = this._actionsMenu();
    const search = domains.search(this.listSource);
    const filters = domains.filters(this.listSource);
    // the rows answer the filter that was applied, not the one being written: while the box is
    // invalid they are shown as stale and the count stands down
    this.effect(() => this.list.stale.value = filters.problem.value !== null);
    // nothing to save in the trash: its rows are read-only until they are restored
    this.effect(() => {
      const trash = this._trash.value;
      save.root.hidden = trash;
      discard.root.hidden = trash;
    });
    for (const control of [add, save, discard, more, search, filters])
      this.own(() => control.dispose());
    this._listOnly(search, filters);
    return this._ribbon = [[add, save, discard, more], [search, filters]];
  }

  /** The table-wide actions of the ⋯ menu — Trash today, Import and Bulk edit as they land. */
  menuActions(): Action[] {
    const trash = this._trash.peek();
    return [{name: trash ? 'Exit trash' : 'Trash', icon: 'trash-restore', requires: 'delete',
      run: () => void this.setTrash(!trash)}];
  }

  /** Enters or leaves trash mode: the list's source answers the deleted rows instead of the live
   * ones, under the same query and search. Through the gate — pending changes of the live list
   * are not dropped silently — and resolves to whether the move was made. */
  async setTrash(on: boolean): Promise<boolean> {
    if (on === this._trash.peek())
      return true;
    if (!await confirmDiscard(this.session, {action: on ? 'open the trash' : 'leave the trash'}))
      return false;
    // the trash is a list view: an entity page open over a live row has nothing to show in it
    if (on && this._page.peek() === 'entity')
      this._show('list', null);
    this._trash.value = on;
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
  private _wireTrash(): void {
    this.effect(() => {
      const on = this._trash.value;
      this.listSource.deleted.value = on ? 'only' : 'exclude';
      this.root.classList.toggle('u2-domain-app-trash', on);
    });
  }

  private _withTrash(path: string): string {
    return this._trash.value ? `${path}${path.includes('?') ? '&' : '?'}trash=1` : path;
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
   * base — the row segment, `?q=`, `?search=` — and its view comes to the front. */
  static openAt(base: string | undefined, address: string): boolean {
    const app = DomainApp._at(base);
    if (app === undefined)
      return false;
    void app.open(address);
    DomainApp._front(app);
    return true;
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
    this._trash.value = params.get('trash') === '1';
    const segment = this._segmentOf(at < 0 ? address : address.slice(0, at));
    const entity = segment ?? params.get('entity');
    if (entity !== null && entity !== '')
      return this._goTo('entity', entity, segment === null ? null : this._keyQuery(segment));
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

  /** The address bar back into the app — only while the shell shows this app's view, and only for
   * an address under one of its routes: another view's URL is not the app's to read. */
  private _restoreFromUrl(): void {
    if (this._isCurrentView() && this._restOf(location.pathname) !== null)
      void this.open(`${location.pathname}${location.search}`);
  }

  /** What an address carries under one of the app's routes: '' for the base itself, '/<row>…'
   * below it; null when the address is another view's. Both routes count: a `/domains/…` link
   * still reaches an app that has rebased onto `/apps/…`. */
  private _restOf(path: string): string | null {
    const here = path.toLowerCase();
    for (const base of [this.base, this._options.base]) {
      if (!here.startsWith(base.toLowerCase()))
        continue;
      const rest = path.slice(base.length);
      if (rest === '' || rest.startsWith('/'))
        return rest;
    }
    return null;
  }

  /** The trailing segment of an address under one of the app's routes — the row it names; null
   * when the address is the app's base itself, or another view's. */
  private _segmentOf(path: string): string | null {
    const rest = this._restOf(path);
    if (rest === null || rest === '')
      return null;
    const segment = rest.slice(1).split('/')[0];
    return segment === '' ? null : decodeURIComponent(segment);
  }

  /** The query finding the row an address segment names, null where the segment IS the id: a uuid
   * is one, and so is anything under a table with no business key; one key column takes the whole
   * segment, a composite key splits on '-' and must match its arity (`domain_entity_view.dart`
   * `_resolve`). A miss shows the form's not-found, as the Dart view does. */
  private _keyQuery(segment: string): FilterGroup | null {
    const key = this.table.info.businessKey;
    if (key.length === 0 || ID.test(segment))
      return null;
    const parts = key.length === 1 ? [segment] : segment.split('-');
    if (parts.length !== key.length)
      return null;
    const schema = this.listSource.schema;
    return Filters.group('and', key.map((column, at) => {
      const property = Filters.property(schema, column);
      const kind = property === null ? Filters.KIND.STRING : Filters.kindOf(property);
      const numeric = kind === Filters.KIND.INT || kind === Filters.KIND.FLOAT;
      return Filters.cond(column, '=', numeric ? Number(parts[at]) : parts[at]);
    }));
  }

  private _state(): AppState {
    const q = this.listSource.query.peek();
    return {page: this._page.peek(), entity: this._entity.peek(), search: this.listSource.search.peek(),
      query: typeof q === 'string' ? q : Filters.format(q)};
  }

  /** A move rather than a rewrite of the page in place: the business key replacing the id the page
   * was addressed by, and the draft that became a row, are the same page. */
  private static _moved(from: AppState, to: AppState): boolean {
    return from.page !== to.page || from.query !== to.query || from.search !== to.search ||
      (from.entity !== to.entity && from.entity !== DomainApp.NEW);
  }

  private _show(page: DomainAppPage, entity: string | null, query: FilterGroup | null = null): void {
    this._close();
    if (page === 'entity' && entity !== null) {
      const draft = entity === DomainApp.NEW;
      const source = draft ? this._source({draft: true}) :
        this._source({query: query ?? Filters.group('and', [Filters.cond('id', '=', entity)]), pageSize: 1});
      // a source of one has no list to make its row current
      source.effect(() => {
        const rows = source.rows.items.value;
        if (source.state.value === 'ready' && rows.length > 0 && source.currentRow.peek() === null)
          source.currentRow.value = rows[0];
      });
      const form = this.runInScope(() => new DomainForm(source, {system: 'footer', include: this._options.include,
        empty: `${this.table.info.singularName || 'Row'} "${entity}" was not found.`}));
      this._formHost.replaceChildren(form.root);
      const {children, history} = this._options;
      this._panes = this.runInScope(() => [
        ...(children === false ? [] : [domains.children(source, children === true ? undefined : children)]),
        ...(history === false ? [] : [domains.history(source)]),
      ]);
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
