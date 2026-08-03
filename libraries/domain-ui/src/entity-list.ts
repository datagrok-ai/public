/**
 * {@link EntityListWidget} — the "list of entities" every domain app opens on:
 * campaigns, profiles, studies, issues. Cards, a brief list or an editable grid
 * over one domain table, with search, the platform's own item rendering and
 * per-item commands, all permission-gated.
 *
 * @module entity-list
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';

import {DomainFrameEditor} from './frame-editor';
import {DomainGrid} from './domain-grid';
import {domainHandler} from './handler';
import {confirmDiscardChanges} from './unsaved';
import {applyDomainUiStyles} from './styles';

/** How {@link EntityListWidget} renders its items — the platform's own three
 * gallery modes ('Card' / 'Brief' / 'Grid' in the built-in Domain View). */
export type EntityListMode = 'cards' | 'brief' | 'grid';

/**
 * Rows a list loads at most, mirroring the built-in Domain View's own cap
 * (`DomainView.rowCap`, `core/client/xamgle/lib/src/views/domain_view.dart`): a
 * gallery stays responsive, and everything past it is a filtering job, not a
 * scrolling one. Beyond the cap the list says so and invites refining the
 * filter, in the platform's words. Full-subset operations (export, Open in Table
 * View) bypass it by running the query themselves.
 */
export const ENTITY_LIST_ROW_CAP = 1000;

/** Options of {@link EntityListWidget.create}. */
export interface EntityListOptions {
  /** What to list; the search box ANDs its term into this filter. */
  query?: DG.DomainQuerySpec;
  /** Initial render mode (default `'cards'`). */
  mode?: EntityListMode;
  /** Show the toolbar — search, mode switch, New (default true). */
  toolbar?: boolean;
  /** Offer in-grid editing in `'grid'` mode; ANDed with the table's `canEdit`
   * (default true — the built-in Domain View's grid mode edits too; pass false
   * for a strictly read-only browse list). */
  editable?: boolean;
  /** Row cap (default {@link ENTITY_LIST_ROW_CAP}). */
  limit?: number;
  /** What opening an item does (default: the platform's Entity View for it). */
  onOpen?: (row: DG.DomainRow) => void;
  /** Asked before anything that would rebuild an editable grid; the hosting
   * {@link AppView} passes its own gate so the prompt covers every editor of the
   * page. Defaults to prompting for this list's own grid. */
  confirmDiscard?: () => Promise<boolean>;
}

/**
 * A list of one domain table's rows.
 *
 * ```ts
 * const list = await EntityListWidget.create(grok.dapi.domains.table('grit.issue'));
 * grok.shell.newView('Issues', [list.root]);
 * ```
 *
 * What it composes:
 * - items rendered by the platform — `ui.renderCard` / `ui.render` dispatch to
 *   the {@link DG.ObjectHandler} registered for the table, so a plugin's cards
 *   show up here without this widget knowing about them, and a table without a
 *   handler gets the reflective default;
 * - commands from {@link DG.DomainObjectHandler.getRibbonActions} — Open, Edit,
 *   Clone, Delete, Share, Watch, History, Copy link — resolved per item when its
 *   context menu opens, so a list of a thousand rows costs no permission probes
 *   until one is asked for;
 * - `'grid'` mode is a {@link DomainGrid}: batch editing, one-transaction save,
 *   the standard conflict flow;
 * - permission gating from {@link DG.DomainTableCapabilities}: no New without
 *   `canInsert`, no editing without `canEdit`, and each item offers only what
 *   its own row permissions allow.
 *
 * Capabilities are SNAPSHOT when the list is built (see {@link DomainGrid}) —
 * rebuild it after a grant change.
 */
export class EntityListWidget {
  /** Mount this: toolbar + list. */
  readonly root: HTMLElement;
  readonly client: DG.DomainTableClient;
  /** The handler every item command goes through (a plugin's, when registered). */
  readonly handler: DG.DomainObjectHandler;
  readonly capabilities: DG.DomainTableCapabilities;
  readonly info: DG.DomainTableInfo;

  /** Holds either the gallery (cards / brief) or the grid: a flex column that
   * takes the space the toolbar leaves, so the platform's own gallery styling
   * applies to the gallery alone and a grid fills the rest of the page. */
  private readonly _body = ui.div([], 'domain-ui-list-body');
  private readonly _banner = ui.divText('', 'grok-domains-cap-banner');
  private readonly _counts = ui.divText('', 'grok-items-view-counts');
  private readonly _toolbar: HTMLElement;
  private readonly _modeIcons: {[mode: string]: HTMLElement} = {};
  private readonly _newButton: HTMLElement | null;
  private readonly _onItemSelected = new rxjs.Subject<DG.DomainRow | null>();
  private readonly _onItemOpened = new rxjs.Subject<DG.DomainRow>();
  private readonly _onRefreshed = new rxjs.Subject<EntityListWidget>();
  private readonly _options: EntityListOptions;
  private readonly _searchColumns: string[];
  private _mode: EntityListMode;
  private _query: DG.DomainQuerySpec;
  private _search = '';
  private _searchTimer: any = null;
  private _rows: DG.DomainRow[] = [];
  private _elements: HTMLElement[] = [];
  private _current: DG.DomainRow | null = null;
  private _grid: DomainGrid | null = null;
  /** Guards against a slow response of an abandoned query overwriting a newer
   * one — every keystroke in the search box starts one. */
  private _generation = 0;

  private constructor(client: DG.DomainTableClient, capabilities: DG.DomainTableCapabilities,
    info: DG.DomainTableInfo, options: EntityListOptions) {
    applyDomainUiStyles();
    this.client = client;
    this.capabilities = capabilities;
    this.info = info;
    this.handler = domainHandler(this.table);
    this._options = options;
    this._mode = options.mode ?? 'cards';
    this._query = options.query ?? {};
    // The search box ANDs a condition into the query's filter, which a raw
    // smart-filter string cannot take (it is parsed server-side, as one whole) —
    // such a list is not searchable, and says so by not offering the box.
    this._searchColumns = typeof this._query.filter === 'string' ? [] : this._searchNames();

    this._newButton = capabilities.canInsert
      ? ui.button(`New ${info.singularName}...`, () => this.createRow(),
        `Create a new ${info.singularName}`) : null;
    for (const mode of ['cards', 'brief', 'grid'] as EntityListMode[])
      this._modeIcons[mode] = ui.iconFA(EntityListWidget._MODE_ICONS[mode],
        () => this.setMode(mode), `Switch to ${mode === 'cards' ? 'card' : mode} view`);
    this._toolbar = ui.divH([
      this._newButton,
      ui.div(Object.values(this._modeIcons), 'grok-items-view-commands,grok-items-view-toggle'),
      this._counts,
      this._searchColumns.length === 0 ? null : ui.div([this._searchInput().root], 'd4-search-bar'),
    ], 'grok-gallery-search-bar');
    ui.setDisplay(this._toolbar, options.toolbar ?? true);
    this.root = ui.divV([this._toolbar, this._banner, this._body], 'domain-ui-entity-list');
    this._syncModeIcons();
  }

  /** Icons of the mode switch — the platform's own (`DataSourceCardView`). */
  private static readonly _MODE_ICONS: {[mode: string]: string} =
    {cards: 'grip-horizontal', brief: 'grip-lines', grid: 'table'};

  /** Builds the list: probes the table's capabilities and metadata, then loads
   * the first page. */
  static async create(client: DG.DomainTableClient,
    options?: EntityListOptions): Promise<EntityListWidget> {
    const address = `${client.schema}.${client.table}`;
    const [capabilities, info] = await Promise.all([
      client.capabilities(), grok.dapi.domains.registry.tableInfo(address)]);
    const list = new EntityListWidget(client, capabilities, info, options ?? {});
    await list.refresh();
    return list;
  }

  /** `'<schema>.<table>'`. */
  get table(): string { return `${this.client.schema}.${this.client.table}`; }

  /** The rows currently listed (capped at {@link ENTITY_LIST_ROW_CAP}). */
  get rows(): DG.DomainRow[] { return this._rows; }

  /** The item the user last clicked, or null. */
  get current(): DG.DomainRow | null { return this._current; }

  /** Fires when the current item changes (null when the list is rebuilt). */
  get onItemSelected(): rxjs.Observable<DG.DomainRow | null> { return this._onItemSelected; }
  /** Fires when an item is opened (double click, or the Open command). */
  get onItemOpened(): rxjs.Observable<DG.DomainRow> { return this._onItemOpened; }
  /** Fires after every (re)load. */
  get onRefreshed(): rxjs.Observable<EntityListWidget> { return this._onRefreshed; }

  /** Render mode; setting it reloads (a grid needs a DataFrame, cards do not). */
  get mode(): EntityListMode { return this._mode; }

  /** The query the list shows, WITHOUT the search term. Setting it reloads —
   * through the unsaved-changes gate, so pending grid edits are never dropped
   * silently. */
  get query(): DG.DomainQuerySpec { return this._query; }

  /** The search term ANDed into {@link query} (empty = no search). */
  get search(): string { return this._search; }

  /** The `'grid'` mode grid, or null in the other modes — where an app finds the
   * {@link DomainFrameEditor} to bind a save bar or a dirty check to. */
  get grid(): DomainGrid | null { return this._grid; }

  /** Every editor this list owns — what a caller's refresh policy gates on. */
  get editors(): DomainFrameEditor[] { return this._grid == null ? [] : [this._grid.editor]; }

  /** Whether anything in the list has unsaved changes. */
  get isDirty(): boolean { return this.editors.some((e) => e.isDirty); }

  /** Switches the render mode (through the gate — a grid is rebuilt). */
  async setMode(mode: EntityListMode): Promise<boolean> {
    if (mode === this._mode)
      return true;
    if (!(await this._confirm(`switch to ${mode}`)))
      return false;
    this._mode = mode;
    this._syncModeIcons();
    await this._reload();
    return true;
  }

  /** Replaces the query the list shows (through the gate). */
  async setQuery(query: DG.DomainQuerySpec): Promise<boolean> {
    if (!(await this._confirm('change the query')))
      return false;
    this._query = query ?? {};
    await this._reload();
    return true;
  }

  /** Sets the search term (through the gate — the list re-queries). */
  async setSearch(text: string): Promise<boolean> {
    if (!(await this._confirm('search')))
      return false;
    this._search = `${text ?? ''}`;
    await this._reload();
    return true;
  }

  /** Opens the platform's create dialog and reloads if a row was saved — the
   * whole "New ..." flow, permission-gated by the button's own presence. */
  async createRow(): Promise<boolean> {
    const saved = await this.handler.editRow();
    if (saved)
      await this.refresh();
    return saved;
  }

  /**
   * Reloads from the server.
   *
   * Every entry point that rebuilds the list goes through the unsaved-changes
   * gate FIRST ({@link EntityListOptions.confirmDiscard}, the hosting AppView's
   * by default): a rebuilt grid is a new {@link DomainFrameEditor} state, and
   * the editors never reconcile — refreshing IS discarding, so nothing here may
   * happen behind the user's back.
   */
  async refresh(): Promise<boolean> {
    if (!(await this._confirm('reload the list')))
      return false;
    await this._reload();
    return true;
  }

  /** The rebuild itself — every caller has already passed the gate. */
  private async _reload(): Promise<void> {
    const generation = ++this._generation;
    const spec = this._spec();
    this._body.innerHTML = '';
    this._body.appendChild(ui.div([ui.div()], 'grok-items-view-loading'));
    try {
      if (this._mode === 'grid')
        await this._showGrid(spec, generation);
      else
        await this._showItems(spec, generation);
    } catch (e: any) {
      if (generation !== this._generation)
        return;
      this._body.innerHTML = '';
      this._body.appendChild(ui.divText(`Cannot list ${this.info.pluralName}: ${e?.message ?? e}`));
      grok.log.error(`${this.table}: ${e?.message ?? e}`);
    }
    if (generation === this._generation)
      this._onRefreshed.next(this);
  }

  /** Releases the list's subscriptions (and its grid's, with its editor's). */
  detach(): void {
    if (this._searchTimer != null)
      clearTimeout(this._searchTimer);
    this._searchTimer = null;
    this._detachGrid();
  }

  // ─────────────────────── internals ─────────────────────────

  /** The columns the search box matches on — the row's identity as the platform
   * defines it: the display-name column and the business key (what the built-in
   * gallery searches). Empty means nothing identifies a row by text, and the box
   * is not offered. */
  private _searchNames(): string[] {
    const columns = this.info.nameColumn == null ? [] : [this.info.nameColumn];
    for (const key of this.info.businessKey)
      if (!columns.includes(key))
        columns.push(key);
    return columns;
  }

  /** The query as run: the caller's, with the row cap and the search term. A
   * query asking for MORE rows than the cap is clamped to it; asking for fewer
   * is honoured (a deep link's `limit` means what it says). */
  private _spec(): DG.DomainQuerySpec {
    const spec: DG.DomainQuerySpec = Object.assign({}, this._query);
    const cap = this._options.limit ?? ENTITY_LIST_ROW_CAP;
    spec.limit = Math.min(this._query.limit ?? cap, cap);
    const text = this._search.trim();
    if (text !== '' && this._searchColumns.length > 0) {
      // Any identity column matching is a hit; the caller's own filter still has
      // to hold, so the two are ANDed.
      const match = DG.or(...this._searchColumns.map((c) => DG.cond(c, 'like', `%${text}%`)));
      spec.filter = spec.filter == null ? match : DG.and(spec.filter as any, match);
    }
    return spec;
  }

  private _searchInput(): DG.InputBase<string> {
    const input = ui.input.search('', {
      placeholder: `Search ${this.info.pluralName} by ${this._searchColumns.join(' or ')}`,
      onValueChanged: (value) => {
        // Debounced: every keystroke would otherwise be a query, and the last
        // one to arrive would not necessarily be the last one typed.
        if (this._searchTimer != null)
          clearTimeout(this._searchTimer);
        this._searchTimer = setTimeout(async () => {
          // Refused (pending changes the user chose to keep): put the box back
          // where the list actually is, instead of showing a term nothing ran.
          if (!(await this.setSearch(value)))
            input.value = this._search;
        }, 300);
      },
    });
    return input;
  }

  private _syncModeIcons(): void {
    for (const mode of Object.keys(this._modeIcons))
      this._modeIcons[mode].classList.toggle('d4-current', mode === this._mode);
  }

  /** The caller's gate when it supplied one (an app prompts for the whole page),
   * this list's own otherwise. */
  private _confirm(action: string): Promise<boolean> {
    return this._options.confirmDiscard != null ? this._options.confirmDiscard()
      : confirmDiscardChanges(this.editors, {action: action, subject: this.table});
  }

  private _detachGrid(): void {
    this._grid?.detach();
    this._grid = null;
  }

  private async _showGrid(spec: DG.DomainQuerySpec, generation: number): Promise<void> {
    // A live grid re-runs its own query (rebuilding the frame and its editing
    // state); only the first one is built from scratch.
    if (this._grid != null)
      await this._grid.editor.refresh(spec);
    else {
      const grid = await DomainGrid.create(this.client, {query: spec,
        editable: this._options.editable ?? true, capabilities: this.capabilities});
      if (generation !== this._generation) {
        grid.detach();
        return;
      }
      this._grid = grid;
    }
    if (generation !== this._generation)
      return;
    this._rows = [];
    this._elements = [];
    this._setCurrent(null);
    this._body.innerHTML = '';
    this._body.appendChild(this._grid!.root);
    this._setCounts(this._grid!.dataFrame.rowCount, null);
  }

  private async _showItems(spec: DG.DomainQuerySpec, generation: number): Promise<void> {
    this._detachGrid();
    const [values, total] = await Promise.all([
      this.client.query(spec) as Promise<{[key: string]: any}[]>,
      this.client.count(spec.filter),
    ]);
    if (generation !== this._generation)
      return;
    // Local materialization — the rows are already here; asking the server for
    // each one again (getById) would be a round trip per card.
    this._rows = values.map((v) => this.handler.rowFrom(v));
    this._elements = this._rows.map((row) => this._renderItem(row));
    // The platform's own gallery container: `mode` is what its stylesheet keys
    // the card wall and the brief list on.
    const gallery = ui.div([], 'grok-gallery-grid');
    gallery.setAttribute('mode', this._mode === 'cards' ? 'Card' : 'Brief');
    if (this._rows.length === 0)
      gallery.appendChild(ui.div([
        ui.iconFA('empty-set'),
        ui.divText(this._search.trim() === '' && this._query.filter == null
          ? `No ${this.info.pluralName} yet.`
          : `No ${this.info.pluralName} match the current filter.`),
      ], 'grok-items-view-empty'));
    else
      for (const element of this._elements)
        gallery.appendChild(element);
    this._body.innerHTML = '';
    this._body.appendChild(gallery);
    this._setCurrent(null);
    this._setCounts(this._rows.length, total);
  }

  /** One item, rendered by the PLATFORM (so a plugin's card wins) and wired to
   * the list's interactions. */
  private _renderItem(row: DG.DomainRow): HTMLElement {
    const element = this._mode === 'cards' ? ui.renderCard(row) : ui.render(row);
    element.onclick = () => this._setCurrent(row);
    element.ondblclick = () => this.open(row);
    element.oncontextmenu = (e) => {
      e.preventDefault();
      this._setCurrent(row);
      this._showMenu(row, e);
    };
    return element;
  }

  /** Opens [row] — the caller's handler when it supplied one, the platform's
   * Entity View otherwise. */
  open(row: DG.DomainRow): void {
    this._onItemOpened.next(row);
    if (this._options.onOpen != null)
      this._options.onOpen(row);
    else
      this.handler.openRow(row as any);
  }

  private _setCurrent(row: DG.DomainRow | null): void {
    this._current = row;
    const index = row == null ? -1 : this._rows.indexOf(row);
    for (let i = 0; i < this._elements.length; i++)
      this._elements[i].classList.toggle('d4-current', i === index);
    this._onItemSelected.next(row);
  }

  /** The item's commands, resolved when the menu opens: the actions are
   * permission-probed per row, so building them for every listed item would be
   * a round trip per item. */
  private _showMenu(row: DG.DomainRow, e: MouseEvent): void {
    const pending = DG.Menu.popup().item('Loading...', () => {});
    pending.show({causedBy: e});
    this.handler.getRibbonActions(row as any).then((actions) => {
      pending.hide();
      const menu = DG.Menu.popup();
      for (const action of actions)
        menu.item(action.name, () => this._run(action));
      menu.show({causedBy: e});
    }).catch((x) => {
      pending.hide();
      grok.shell.error(`${x?.message ?? x}`);
    });
  }

  /** Runs an item command and reloads when it may have changed the list. */
  private async _run(action: DG.DomainAction): Promise<void> {
    const changed = await action.run();
    if (changed !== false && action.name !== 'Open' && action.name !== 'Copy link')
      await this.refresh();
  }

  private _setCounts(loaded: number, total: number | null): void {
    this._counts.textContent = total == null || total === loaded ? `${loaded}` : `${loaded} / ${total}`;
    const cap = this._options.limit ?? ENTITY_LIST_ROW_CAP;
    const capped = total != null && loaded >= cap && total > cap;
    ui.setDisplay(this._banner, capped);
    if (capped)
      this._banner.textContent =
        `Showing the first ${cap} of ${total} rows — refine the filter to narrow the results.`;
  }
}
