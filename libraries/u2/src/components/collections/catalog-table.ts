/* The catalog table: the records a person comes back to (campaigns, profiles, templates) in a real
   table whose columns fit their content. A header gesture sorts, a grouping folds the rows into
   sections with counts, a menu chooses the columns, each row carries its actions as icons and as a
   context menu, and every choice is kept through a persistence adapter. BasicTable renders it; this
   control owns the order, the groups and the choices (`state`) and reports them (`getState`). */
import {signal, computed, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {Input} from '../../core/input-base.js';
import type {LiveOption} from '../../core/input-base.js';
import {AsyncValue} from '../../core/async-value.js';
import {Tooltip} from '../../core/tooltip.js';
import {button, effect, link} from '../../core/elements.js';
import type {IItemsState, IItemState} from '../../core/widget-like.js';
import {Access} from '../../core/access.js';
import {Format} from '../../core/format.js';
import {BasicTable} from './table.js';
import type {TableColumn, TableSection} from './table.js';
import {SortOrder} from './sort-order.js';
import {Menu} from '../navigation/menu.js';
import {iconButton} from '../actions/buttons.js';
import {rowActions, allowedActions} from '../actions/actions.js';
import type {Action} from '../actions/actions.js';
import {icon} from '../display/icon.js';
import {skeleton, emptyState} from '../display/async-view.js';
import {notify} from '../display/notify.js';

export type CatalogValueType = 'auto' | 'string' | 'number' | 'date';

export interface CatalogColumn<T> {
  /** Its identity in `state`, the kept choices and `data-column`: 'Code', 'campaignFields.Target'. */
  name: string;
  /** The header text; `name` by default. */
  header?: string;
  /** What the column sorts by and, through `Format.value`, what its cell shows without a `render`. */
  value: (item: T) => unknown;
  /** The cell: a string is its text, an element its content (a status dot, a chip). */
  render?: (item: T) => HTMLElement | string;
  /** The cell is a link that activates its row (`onActivate`); a string is the link's tooltip. */
  link?: boolean | string;
  /** How the values compare. 'auto' (default) looks at the column: numbers when every value is one,
   * dates when every value is a date (a `Date`, or text like 2026/3/30 or 2026-03-30), else text in
   * natural order (KHL-2 before KHL-12). Empty values go last either way. */
  type?: CatalogValueType;
  /** The column's own ascending comparator, in place of `type`. */
  compare?: (a: T, b: T) => number;
  /** True by default. */
  sortable?: boolean;
  align?: 'left' | 'right' | 'center';
  /** Always shown, first, never in the columns menu (a record's name). */
  pinned?: boolean;
  /** The columns menu lists it in this submenu ('Campaign Fields'). */
  menuGroup?: string;
  /** Grey text after it in the columns menu (the templates that define a field). */
  menuDetail?: string;
  /** The header's hover text, before how to sort by it. */
  description?: string;
  /** A CSS length: a longer cell ends in an ellipsis and shows its whole text on hover. */
  maxWidth?: string;
  width?: string;
}

export interface CatalogGrouping<T> {
  /** Its identity in `state` and the kept choices. */
  name: string;
  /** The group menu's text; `name` by default. */
  label?: string;
  /** The group an item falls in, by its text (`Format.value`, as a cell reads the same value); a key with no text
   * (null, undefined, '', a record) goes to `emptyGroupLabel`. */
  key: (item: T) => unknown;
  /** A group header's text: the key's text by default (a template's name for its key). */
  groupLabel?: (key: string, first: T) => string;
  /** The order of the groups; by default the order they first appear in. */
  compare?: (a: string, b: string) => number;
  menuGroup?: string;
  menuDetail?: string;
}

export interface CatalogState {
  /** The columns shown after the pinned ones, in order. A name no column has (yet) is kept: columns
   * that arrive with the items show once they do. */
  columns: string[];
  /** The order, as every u2 table keeps it ({@link SortOrder}): 'Code', '!Created', '' for the default order. */
  sort: string;
  groupBy: string | null;
  /** Per grouping, the keys of its folded groups. */
  collapsed: {[groupBy: string]: string[]};
}

/** Where the choices live: `load` once when the table is built, `save` after every change. */
export interface CatalogPersistence {
  load(): Partial<CatalogState>;
  save(state: CatalogState): void;
}

/** What the user may do with a row: its edit actions show with 'edit', none with 'view', and with
 * 'pending' (still asked) the cell keeps their room empty. */
export type CatalogAccess = 'edit' | 'view' | 'pending';

export interface CatalogChromeOptions {
  /** The columns icon's tooltip: 'Choose columns'. */
  columns?: string;
  /** The grouping icon's tooltip for the current grouping's menu label ('None'). */
  group?: (current: string) => string;
  /** The refresh icon's tooltip: 'Refresh'. */
  refresh?: string;
}

export interface CatalogTableOptions<T> {
  name?: string;
  /** The table's accessible name. */
  label?: string;
  /** Live. From an `AsyncValue` the table shows a skeleton on the first load, keeps the rows while a
   * newer answer is coming (`aria-busy`), and its refresh is the refresh icon's and Retry's. */
  items: readonly T[] | ReadonlySignal<readonly T[]> | AsyncValue<readonly T[]>;
  keyOf: (item: T) => string;
  columns: LiveOption<CatalogColumn<T>[]>;
  groupings?: LiveOption<CatalogGrouping<T>[]>;
  /** The columns shown before the user chose any: every column the table starts with by default. */
  defaultColumns?: string[];
  /** The order while no column is sorted, and the order a sort breaks its ties by; the items' own by default. */
  defaultOrder?: (a: T, b: T) => number;
  /** Sorts by a column in place of the column's comparator — an app whose order probes the whole
   * column its own way. Gets the items in the default order. */
  sorter?: (items: T[], column: CatalogColumn<T>, ascending: boolean) => T[];
  /** What sorts by a header: a double click (default, as the platform grid) or a click. */
  sortGesture?: 'dblclick' | 'click';
  /** The group of the items with no key: 'No information'. */
  emptyGroupLabel?: string;
  /** It comes after the other groups (default true). */
  emptyGroupLast?: boolean;
  /** A grouping with one group still shows its header (default false: the rows alone). */
  singleGroupHeader?: boolean;
  /** A row's actions: the ones with an icon in its last cell, every one in its context menu. */
  actions?: (item: T) => Action[];
  /** What the user may do with each row (see {@link CatalogAccess}); 'edit' for every row by default. */
  access?: (item: T) => CatalogAccess | ReadonlySignal<CatalogAccess>;
  /** The icons show on every row ('always') or on the row under the pointer ('hover', default). */
  actionsVisible?: 'always' | 'hover';
  /** A view-only row shows this icon where the actions would be. */
  viewOnly?: {icon: string, tooltip: string};
  /** A link column's click, or Enter on it. */
  onActivate?: (item: T) => void;
  /** What shows when there are no items: a sentence, or an element. */
  empty?: string | (() => HTMLElement);
  /** The line under the table ('37 campaigns'). */
  footer?: (shown: readonly T[]) => string;
  /** Row and header height: 'control' (28px, default) or 'inline' (22px, as `ui.table`). */
  rowHeight?: 'control' | 'inline';
  /** Loading while items come (the first load shows a skeleton); with an `AsyncValue` its own state. */
  loading?: ReadonlySignal<boolean>;
  /** Why the items could not come; with no items shown, it shows instead of the table with Retry. */
  error?: ReadonlySignal<string | null>;
  /** What the refresh icon and Retry run; an `AsyncValue`'s refresh by default. */
  refresh?: () => Promise<unknown>;
  /** A localStorage key (`u2.catalog.<key>`, the whole state as JSON), or an adapter of the app's own. */
  persist?: string | CatalogPersistence;
  /** The choices to start from where nothing was kept. */
  state?: Partial<CatalogState>;
}

interface Group<T> {
  key: string;
  label: string;
  items: T[];
  folded: boolean;
}

interface View<T> {
  /** Every item, in the order shown (group by group when grouped). */
  ordered: readonly T[];
  /** What the BasicTable shows: section headers and rows. */
  rows: readonly (T | TableSection)[];
  groups: readonly Group<T>[];
  /** Whether the groups have headers. */
  sectioned: boolean;
  /** A group header's text by its key. */
  labels: ReadonlyMap<string, string>;
}

const NONE: readonly string[] = [];
const NUMBER = /^[-+]?(\d+\.?\d*|\.\d+)(e[-+]?\d+)?$/i;
const DATE = /^\d{4}[-/.]\d{1,2}[-/.]\d{1,2}([T ]\d{1,2}:\d{2}(:\d{2}(\.\d+)?)?(Z|[+-]\d{2}:?\d{2})?)?$/;
const COLLATOR = new Intl.Collator(undefined, {numeric: true});

export class CatalogTable<T> extends Control {
  /** The choices: shown columns, sort, grouping, folded groups. Code (and an assistant) writes them
   * here; every write is kept through `persist`. */
  readonly state: Signal<CatalogState>;
  /** The items in the order shown. */
  readonly shown: ReadonlySignal<readonly T[]>;
  /** Items are loading or refreshing. */
  readonly loading: ReadonlySignal<boolean>;

  private readonly _keyOf: (item: T) => string;
  private readonly _items: ReadonlySignal<readonly T[]>;
  private readonly _error: ReadonlySignal<string | null>;
  private readonly _refreshFn: (() => Promise<unknown>) | undefined;
  private readonly _refreshing = signal(false);
  private readonly _columns: Signal<ReadonlySignal<readonly CatalogColumn<T>[]>>;
  private readonly _groupings: Signal<ReadonlySignal<readonly CatalogGrouping<T>[]>>;
  private readonly _o: CatalogTableOptions<T>;
  private readonly _emptyLabel: string;
  private readonly _view: ReadonlySignal<View<T>>;
  private readonly _shownColumns: ReadonlySignal<readonly CatalogColumn<T>[]>;
  private readonly _table: BasicTable<T>;
  private readonly _scroll: HTMLElement;
  private readonly _host: HTMLElement;
  private readonly _footer: HTMLElement;
  private _modeScope: Scope | undefined;

  constructor(options: CatalogTableOptions<T>) {
    super(document.createElement('div'));
    this._o = options;
    this._keyOf = options.keyOf;
    this._emptyLabel = options.emptyGroupLabel ?? 'No information';
    this.root.classList.add('u2-catalog');
    this.root.dataset.u2 = 'catalog-table';
    if (options.name !== undefined)
      this.name = options.name;
    if (options.label !== undefined)
      this.root.setAttribute('aria-label', options.label);
    if (options.actionsVisible === 'always')
      this.root.classList.add('u2-catalog-actions-always');

    const items = options.items;
    if (items instanceof AsyncValue) {
      const last = signal<readonly T[]>([]);
      this.effect(() => {
        const s = items.state.value;
        if (s.kind === 'ready')
          last.value = s.value;
      });
      this._items = last;
      this.loading = computed(() => items.state.value.kind === 'loading' || this._refreshing.value);
      this._error = computed(() => {
        const s = items.state.value;
        return s.kind === 'error' ? `Could not load: ${s.message}` : null;
      });
      this._refreshFn = options.refresh ?? (() => items.refresh());
    }
    else {
      this._items = items instanceof Signal ? items : signal(items);
      const loading = options.loading;
      this.loading = computed(() => (loading?.value ?? false) || this._refreshing.value);
      this._error = options.error ?? signal(null);
      this._refreshFn = options.refresh;
    }

    const columns = options.columns;
    this._columns = signal(columns instanceof Signal ? columns : signal(columns));
    const groupings = options.groupings ?? [];
    this._groupings = signal(groupings instanceof Signal ? groupings : signal(groupings));

    const persist = typeof options.persist === 'string' ? CatalogTable.storage(options.persist) : options.persist;
    const initial: CatalogState = {
      columns: options.defaultColumns ?? this._columns.peek().peek().filter((c) => !c.pinned).map((c) => c.name),
      sort: '', groupBy: null, collapsed: {},
    };
    this.state = signal(CatalogTable._merge(CatalogTable._merge(initial, options.state ?? {}), persist?.load() ?? {}));

    // the columns a recipe computes are new objects on every items change: the same shape is the same column, so
    // the table keeps its header and rows (and the focus in them); a cell reads its column by name
    this._shownColumns = CatalogTable._same(() => {
      const all = this._columns.value.value;
      const chosen = this.state.value.columns;
      const pinned = all.filter((c) => c.pinned);
      const rest = chosen.map((n) => all.find((c) => c.name === n && !c.pinned))
        .filter((c): c is CatalogColumn<T> => c !== undefined);
      return [...pinned, ...rest];
    }, CatalogTable._shape);
    this._view = computed(() => this._compute(this._items.value, this._columns.value.value,
      this._groupings.value.value, this.state.value));
    this.shown = computed(() => this._view.value.ordered);

    this._table = this.build(() => new BasicTable<T>({
      columns: computed(() => this._tableColumns(this._shownColumns.value)),
      items: computed(() => this._view.value.rows),
      keyOf: options.keyOf,
      row: (_item, _index, tr) => tr.classList.add('u2-row'),
      sort: computed(() => this.state.value.sort),
      onSort: (sort) => this._patchState({sort}),
      sortGesture: options.sortGesture ?? 'dblclick',
      collapsed: computed(() => {
        const s = this.state.value;
        return s.groupBy === null ? NONE : s.collapsed[s.groupBy] ?? NONE;
      }),
      onCollapse: (keys) => {
        const s = this.state.peek();
        if (s.groupBy !== null)
          this._patchState({collapsed: {...s.collapsed, [s.groupBy]: [...keys]}});
      },
      contextActions: options.actions === undefined ? undefined : (item) => allowedActions(options.actions!(item),
        {access: this._accessNow(item) === 'edit' ? Access.full : Access.readOnly}),
      sectionLabel: (key) => this._view.peek().labels.get(key) ?? key,
      sectionCounts: true,
      rowHeight: options.rowHeight ?? 'control',
    }));
    this._scroll = document.createElement('div');
    this._scroll.className = 'u2-catalog-scroll';
    this._scroll.append(this._table.root);
    this._host = document.createElement('div');
    this._host.className = 'u2-catalog-body';
    this._footer = document.createElement('div');
    this._footer.className = 'u2-catalog-footer';
    this._footer.dataset.u2Part = 'footer';
    this.root.append(this._host, this._footer);

    if (persist !== undefined) {
      let first = true;
      this.effect(() => {
        const s = this.state.value;
        if (first)
          first = false;
        else
          untracked(() => persist.save(s));
      });
    }
    const mode = computed(() => this._items.value.length > 0 ? 'table' : this._error.value !== null ? 'error' :
      this.loading.value ? 'loading' : 'empty');
    this.effect(() => {
      const m = mode.value;
      untracked(() => this._showMode(m));
    });
    this.effect(() => {
      const busy = this.loading.value;
      if (busy)
        this.root.setAttribute('aria-busy', 'true');
      else
        this.root.removeAttribute('aria-busy');
    });
    this.effect(() => {
      const shown = this._view.value.ordered;
      const text = mode.value === 'table' ? options.footer?.(shown) ?? '' : '';
      this._footer.textContent = text;
      this._footer.hidden = text === '';
    });
    this._watchOverflow();
  }

  /** The heading's icon buttons — columns, grouping, and refresh when there is one — each opening its
   * menu under itself. Built in this table's scope: they go with it. */
  chrome(options?: CatalogChromeOptions): HTMLElement[] {
    return this.build(() => {
      const disabled = computed(() => this._items.value.length === 0 && this.loading.value ? 'Loading…' : null);
      const groupTip = options?.group ?? ((g: string) => `Group by. Current: ${g}`);
      const columns = iconButton('eye', () => {}, {tooltip: options?.columns ?? 'Choose columns', disabled,
        menu: () => this.build(() => this.columnsMenu())});
      const group = iconButton('layer-group', () => {}, {disabled,
        tooltip: computed(() => groupTip(this._groupingLabel(this.state.value.groupBy))),
        menu: () => this.build(() => this.groupMenu())});
      const out: HTMLElement[] = [columns, group];
      const name = this.name;
      const parts = ['columns', 'group'];
      if (this._refreshFn !== undefined) {
        const sync = iconButton('sync', () => void this.refresh(), {tooltip: options?.refresh ?? 'Refresh', disabled});
        effect(() => {
          const busy = this.loading.value;
          sync.classList.toggle('u2-catalog-refreshing', busy);
          if (busy)
            sync.setAttribute('aria-busy', 'true');
          else
            sync.removeAttribute('aria-busy');
        });
        out.push(sync);
        parts.push('refresh');
      }
      if (name !== undefined)
        out.forEach((el, i) => el.dataset.u2Name = `${name}-${parts[i]}`);
      return out;
    });
  }

  /** Runs the refresh; a failure over rows already shown keeps them and says why in a balloon. */
  async refresh(): Promise<void> {
    const fn = this._refreshFn;
    if (fn === undefined || this._refreshing.peek())
      return;
    this._refreshing.value = true;
    let failure: string | null = null;
    try {
      await fn();
      failure = this._error.peek();
    }
    catch (e) {
      failure = e instanceof Error ? e.message : String(e);
    }
    finally {
      if (!this.scope.isDisposed)
        this._refreshing.value = false;
    }
    if (failure !== null && this._items.peek().length > 0)
      notify.error(`Could not refresh: ${failure}`);
  }

  /** The columns menu: each column a check item, a column's `menuGroup` a submenu. */
  columnsMenu(): Menu {
    const menu = Menu.popup();
    const shown = this.state.peek().columns;
    const groups = new Map<string, CatalogColumn<T>[]>();
    for (const c of this._columns.peek().peek()) {
      if (c.pinned)
        continue;
      if (c.menuGroup !== undefined)
        groups.set(c.menuGroup, [...groups.get(c.menuGroup) ?? [], c]);
      else
        this._columnItem(menu, c, shown);
    }
    for (const [label, cols] of groups) {
      menu.separator().group(label);
      for (const c of cols)
        this._columnItem(menu, c, shown);
      menu.endGroup();
    }
    return menu;
  }

  /** The group menu: None, then each grouping, a grouping's `menuGroup` a submenu. */
  groupMenu(): Menu {
    const menu = Menu.popup();
    const current = this.state.peek().groupBy;
    menu.item('None', () => this.groupBy(null), {check: current === null});
    const groups = new Map<string, CatalogGrouping<T>[]>();
    for (const g of this._groupings.peek().peek()) {
      if (g.menuGroup !== undefined)
        groups.set(g.menuGroup, [...groups.get(g.menuGroup) ?? [], g]);
      else
        menu.item(g.label ?? g.name, () => this.groupBy(g.name), {check: current === g.name, detail: g.menuDetail});
    }
    for (const [label, gs] of groups) {
      menu.separator().group(label);
      for (const g of gs)
        menu.item(g.label ?? g.name, () => this.groupBy(g.name), {check: current === g.name, detail: g.menuDetail});
      menu.endGroup();
    }
    return menu;
  }

  /** Sorts by a column, or back to the default order with null. */
  sortBy(column: string | null, ascending = true): void {
    this._patchState({sort: column === null ? '' : SortOrder.of(column, !ascending)});
  }

  /** Groups by a grouping, or not at all with null. */
  groupBy(grouping: string | null): void {
    this._patchState({groupBy: grouping});
  }

  /** Shows a column at the end of the shown ones, or hides it. */
  showColumn(name: string, show = true): void {
    const columns = this.state.peek().columns.filter((c) => c !== name);
    this._patchState({columns: show ? [...columns, name] : columns});
  }

  /** Folds or unfolds a group of the current grouping, by its key. */
  foldGroup(key: string, folded = true): void {
    const by = this.state.peek().groupBy;
    if (by === null)
      return;
    const now = this.state.peek().collapsed[by] ?? [];
    const next = folded ? [...now.filter((k) => k !== key), key] : now.filter((k) => k !== key);
    this._patchState({collapsed: {...this.state.peek().collapsed, [by]: next}});
  }

  /** Group headers (`section:<key>`, `expanded`) and rows (by `keyOf`, the first column's text), in the
   * order shown; a row of a folded group is not rendered. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const view = this._view.peek();
    const first = this._shownColumns.peek()[0];
    const list: IItemState[] = [];
    const row = (item: T, rendered: boolean) => list.push({key: this._keyOf(item), index: list.length,
      label: first === undefined ? '' : Format.value(first.value(item)), rendered});
    if (!view.sectioned) {
      for (const item of view.ordered)
        row(item, true);
    }
    else {
      for (const g of view.groups) {
        list.push({key: `section:${g.key}`, index: list.length, label: g.label, expanded: !g.folded, rendered: true});
        for (const item of g.items)
          row(item, !g.folded);
      }
    }
    return Input.itemsState(list, opts?.all);
  }

  /** `extra`: `sort` ('Molecules desc'), `groupBy`, `columns` and `hidden` (comma lists of names),
   * `collapsed` (the folded group keys), `groups` ('In Progress (8), Submitted (3)'), `rows` (rendered)
   * and `total`. */
  getState(): ReturnType<Control['getState']> {
    const s = this.state.peek();
    const view = this._view.peek();
    const shown = this._shownColumns.peek();
    const hidden = this._columns.peek().peek().filter((c) => !c.pinned && !shown.some((x) => x.name === c.name));
    const rows = view.groups.reduce((n, g) => n + (g.folded && view.sectioned ? 0 : g.items.length), 0);
    const sort = SortOrder.parse(s.sort);
    return {...super.getState(), extra: {
      sort: sort === null ? '' : `${sort.column} ${sort.descending ? 'desc' : 'asc'}`,
      groupBy: s.groupBy ?? '',
      columns: shown.map((c) => c.name).join(','),
      hidden: hidden.map((c) => c.name).join(','),
      collapsed: view.groups.filter((g) => g.folded).map((g) => g.key).join(','),
      groups: view.sectioned ? view.groups.map((g) => `${g.label} (${g.items.length})`).join(', ') : '',
      rows,
      total: view.ordered.length,
    }};
  }

  /** A row of a folded group unfolds it first. Exact key first, then the label; two labels reject. */
  async revealItem(keyOrLabel: string): Promise<HTMLElement | null> {
    const items = this.getItems({all: true}).list;
    let match = items.find((i) => i.key === keyOrLabel);
    if (match === undefined) {
      const wanted = keyOrLabel.trim().toLowerCase();
      const labelled = items.filter((i) => i.label.toLowerCase() === wanted);
      if (labelled.length > 1)
        throw new Error(`u2: ${labelled.length} items match "${keyOrLabel}": ${labelled.map((i) => i.key).join(', ')}`);
      match = labelled[0];
    }
    if (match === undefined)
      return null;
    if (!match.rendered) {
      const group = this._view.peek().groups.find((g) => g.items.some((item) => this._keyOf(item) === match!.key));
      if (group !== undefined)
        this.foldGroup(group.key, false);
    }
    for (const el of this._table.root.querySelectorAll<HTMLElement>('[data-u2-item]')) {
      if (el.dataset.u2Item === match.key)
        return el;
    }
    return null;
  }

  /** A localStorage adapter: the whole state as JSON under `u2.catalog.<key>`. */
  static storage(key: string): CatalogPersistence {
    const id = `u2.catalog.${key}`;
    return {
      load: () => {
        try {
          const text = globalThis.localStorage?.getItem(id);
          return text ? JSON.parse(text) as Partial<CatalogState> : {};
        }
        catch {
          return {};
        }
      },
      save: (state) => {
        try {
          globalThis.localStorage?.setItem(id, JSON.stringify(state));
        }
        catch {
        }
      },
    };
  }

  /** How 'auto' reads a column: numbers when every value that is not empty is one, then dates, else text. */
  private static _typeOf(values: Iterable<unknown>): 'string' | 'number' | 'date' {
    let numbers = true;
    let dates = true;
    let any = false;
    for (const v of values) {
      if (CatalogTable._empty(v))
        continue;
      any = true;
      numbers &&= typeof v === 'number' ? Number.isFinite(v) : typeof v === 'string' && NUMBER.test(v.trim());
      dates &&= v instanceof Date ? !isNaN(v.getTime()) :
        typeof v === 'string' && DATE.test(v.trim()) && !isNaN(Date.parse(v));
      if (!numbers && !dates)
        break;
    }
    return !any ? 'string' : numbers ? 'number' : dates ? 'date' : 'string';
  }

  /** `items` sorted by the column — stable, so equal values keep the order they came in; empty values
   * last both ways. */
  static sortItems<T>(items: T[], column: CatalogColumn<T>, ascending: boolean): T[] {
    const sign = ascending ? 1 : -1;
    if (column.compare !== undefined)
      return items.sort((a, b) => sign * column.compare!(a, b));
    const values = new Map(items.map((item) => [item, column.value(item)]));
    const type = column.type === undefined || column.type === 'auto' ? CatalogTable._typeOf(values.values()) :
      column.type;
    const keyOf = (v: unknown) => type === 'number' ? Number(v) :
      type === 'date' ? (v instanceof Date ? v.getTime() : Date.parse(String(v))) : Format.value(v);
    const keys = new Map<T, number | string>();
    for (const [item, v] of values) {
      if (!CatalogTable._empty(v))
        keys.set(item, keyOf(v));
    }
    return items.sort((a, b) => {
      const ka = keys.get(a);
      const kb = keys.get(b);
      if (ka === undefined || kb === undefined)
        return ka === kb ? 0 : ka === undefined ? 1 : -1;
      return sign * (typeof ka === 'number' ? ka - (kb as number) : COLLATOR.compare(ka, kb as string));
    });
  }

  private _compute(items: readonly T[], columns: readonly CatalogColumn<T>[],
    groupings: readonly CatalogGrouping<T>[], s: CatalogState): View<T> {
    let ordered = [...items];
    if (this._o.defaultOrder !== undefined)
      ordered.sort(this._o.defaultOrder);
    const sort = SortOrder.parse(s.sort);
    const column = sort === null ? undefined : columns.find((c) => c.name === sort.column);
    if (sort !== null && column !== undefined && column.sortable !== false) {
      ordered = this._o.sorter !== undefined ? this._o.sorter(ordered, column, !sort.descending) :
        CatalogTable.sortItems(ordered, column, !sort.descending);
    }
    const grouping = s.groupBy === null ? undefined : groupings.find((g) => g.name === s.groupBy);
    if (grouping === undefined) {
      return {ordered, rows: ordered, groups: [{key: '', label: '', items: ordered, folded: false}],
        sectioned: false, labels: new Map()};
    }

    const by = new Map<string, T[]>();
    for (const item of ordered) {
      const text = Format.value(grouping.key(item));
      const key = text === '' ? this._emptyLabel : text;
      const members = by.get(key);
      if (members === undefined)
        by.set(key, [item]);
      else
        members.push(item);
    }
    const emptyLast = (this._o.emptyGroupLast ?? true) && by.has(this._emptyLabel);
    const keys = [...by.keys()].filter((k) => !emptyLast || k !== this._emptyLabel);
    if (grouping.compare !== undefined)
      keys.sort(grouping.compare);
    if (emptyLast)
      keys.push(this._emptyLabel);
    const folded = s.collapsed[grouping.name] ?? [];
    const groups = keys.map((key): Group<T> => {
      const members = by.get(key)!;
      const label = key === this._emptyLabel ? key : grouping.groupLabel?.(key, members[0]) ?? key;
      return {key, label, items: members, folded: folded.includes(key)};
    });
    const labels = new Map(groups.map((g) => [g.key, g.label]));
    const sectioned = groups.length > 1 || (this._o.singleGroupHeader ?? false) && groups.length === 1;
    if (!sectioned)
      return {ordered, rows: ordered, groups: groups.map((g) => ({...g, folded: false})), sectioned, labels};
    const rows: (T | TableSection)[] = [];
    for (const g of groups)
      rows.push({section: g.key}, ...g.items);
    return {ordered: groups.flatMap((g) => g.items), rows, groups, sectioned, labels};
  }

  private _tableColumns(columns: readonly CatalogColumn<T>[]): TableColumn<T>[] {
    const out: TableColumn<T>[] = columns.map((c, i) => ({
      header: c.header ?? c.name,
      name: c.name,
      render: (item: T) => this._cell(c.name, item),
      align: c.align,
      width: c.width,
      maxWidth: c.maxWidth,
      sortable: c.sortable ?? true,
      tooltip: c.description,
      pinned: i === 0 && c.pinned ? 'start' as const : undefined,
    }));
    if (this._o.actions !== undefined)
      out.push({header: '', render: (item: T) => this._actionsCell(item), pinned: 'end'});
    return out;
  }

  /** The cell of the column named so, as the columns are now: one of the same shape may be a new object. */
  private _cell(name: string, item: T): HTMLElement | string {
    const column = this._columns.peek().peek().find((c) => c.name === name)!;
    if (column.render !== undefined)
      return column.render(item);
    const text = Format.value(column.value(item));
    if (!column.link)
      return text;
    const el = link(text, () => this._o.onActivate?.(item));
    if (typeof column.link === 'string')
      Tooltip.bind(el, column.link);
    return el;
  }

  /** The row's icon actions for what its access allows, the view-only mark, or — while the access is
   * asked — their room kept empty. Built in the row's scope. */
  private _actionsCell(item: T): HTMLElement {
    const cell = document.createElement('div');
    cell.className = 'u2-catalog-actions';
    const access = this._o.access?.(item) ?? 'edit';
    const fill = (a: CatalogAccess) => {
      // an action outlives the row's icons: what it opens (a confirmation) is the table's, not the row's
      const actions = this._o.actions!(item).map((action) =>
        ({...action, run: () => this.build(() => action.run())}));
      const el = rowActions(actions, {access: a === 'view' ? Access.readOnly : Access.full});
      if (a === 'pending') {
        el.classList.add('u2-catalog-actions-pending');
        el.setAttribute('aria-hidden', 'true');
      }
      else if (a === 'view' && this._o.viewOnly !== undefined)
        el.append(icon(this._o.viewOnly.icon, {tooltip: this._o.viewOnly.tooltip, cls: 'u2-catalog-view-only'}));
      cell.replaceChildren(el);
    };
    if (typeof access === 'string')
      fill(access);
    else {
      effect(() => {
        const a = access.value;
        untracked(() => fill(a));
      });
    }
    return cell;
  }

  private _accessNow(item: T): CatalogAccess {
    const access = this._o.access?.(item) ?? 'edit';
    return typeof access === 'string' ? access : access.peek();
  }

  /** The table, its first-load skeleton, its empty state or the reason it could not load — each built
   * in a scope of its own, the table kept. */
  private _showMode(mode: string): void {
    this._modeScope?.dispose();
    this._modeScope = undefined;
    if (mode === 'table') {
      this._host.replaceChildren(this._scroll);
      return;
    }
    const scope = new Scope(this.scope, {label: 'catalog-state'});
    this._modeScope = scope;
    this._host.replaceChildren(scope.build(() => this._stateView(mode)));
  }

  private _stateView(mode: string): HTMLElement {
    if (mode === 'loading')
      return skeleton(8);
    if (mode === 'error') {
      const retry = this._refreshFn === undefined ? [] : [button('Retry', () => void this.refresh())];
      return emptyState(this._error.peek() ?? '', {error: true, actions: retry});
    }
    const empty = this._o.empty;
    return typeof empty === 'function' ? empty() : emptyState(empty ?? 'Nothing here yet', {icon: 'folder-open'});
  }

  /** `u2-catalog-overflow` while the table is wider than its host: the pinned end column casts an edge. */
  private _watchOverflow(): void {
    const scroll = this._scroll;
    const observer = new ResizeObserver(() =>
      this.root.classList.toggle('u2-catalog-overflow', scroll.scrollWidth > scroll.clientWidth + 1));
    observer.observe(scroll);
    observer.observe(this._table.root);
    this.own(() => observer.disconnect());
  }

  /** A column's check item: the menu stays open, so several columns are toggled in one visit. */
  private _columnItem(menu: Menu, c: CatalogColumn<T>, shown: readonly string[]): void {
    menu.item(c.header ?? c.name, () => this.showColumn(c.name, !this.state.peek().columns.includes(c.name)),
      {check: shown.includes(c.name), detail: c.menuDetail, keepOpen: true});
  }

  private _groupingLabel(by: string | null): string {
    const g = by === null ? undefined : this._groupings.peek().peek().find((x) => x.name === by);
    return g === undefined ? 'None' : g.label ?? g.name;
  }

  private _patchState(patch: Partial<CatalogState>): void {
    this.state.value = {...this.state.peek(), ...patch};
  }

  /** `base` with what `patch` holds of the right shape: what an app or a browser kept may be anything — a sort
   * kept before it was a string (`{column, ascending}`, or null for none) included. */
  private static _merge(base: CatalogState, patch: Partial<CatalogState>): CatalogState {
    const p = (patch ?? {}) as Record<string, unknown>;
    const strings = (x: unknown): x is string[] => Array.isArray(x) && x.every((s) => typeof s === 'string');
    const kept = p.sort as {column?: unknown, ascending?: unknown} | string | null | undefined;
    const sort = typeof kept === 'string' ? kept : kept === null ? '' :
      typeof kept?.column === 'string' && typeof kept.ascending === 'boolean' ?
        SortOrder.of(kept.column, !kept.ascending) : base.sort;
    const collapsed = p.collapsed;
    return {
      columns: strings(p.columns) ? [...p.columns] : base.columns,
      sort,
      groupBy: p.groupBy === null || typeof p.groupBy === 'string' ? p.groupBy as string | null : base.groupBy,
      collapsed: typeof collapsed === 'object' && collapsed !== null && !Array.isArray(collapsed) ?
        Object.fromEntries(Object.entries(collapsed).filter(([, v]) => strings(v))) as CatalogState['collapsed'] :
        base.collapsed,
    };
  }

  private static _empty(v: unknown): boolean {
    return v === null || v === undefined || v === '' || (typeof v === 'number' && isNaN(v));
  }

  /** What makes a column look the way it does in the header and the rows; its `value` and `render` are read by
   * name when a cell is built. */
  private static _shape<X>(c: CatalogColumn<X>): string {
    return JSON.stringify([c.name, c.header, c.width, c.align, c.maxWidth, c.pinned, c.sortable, c.description,
      c.link]);
  }

  /** A computed list that keeps its last array while its elements have the same keys in the same order, so what
   * follows it (the table's columns) is not built again for nothing. */
  private static _same<X>(fn: () => readonly X[], keyOf: (x: X) => string): ReadonlySignal<readonly X[]> {
    let last: readonly X[] = [];
    return computed(() => {
      const next = fn();
      if (next.length !== last.length || next.some((x, i) => keyOf(x) !== keyOf(last[i])))
        last = next;
      return last;
    });
  }
}
