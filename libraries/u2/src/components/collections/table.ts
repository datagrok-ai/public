/* Basic (non-virtualized) table — the `ui.table`/`ui.tableFromMap` parity control for small data:
   real table semantics, compact platform density, cells rendered by the caller. Rows all live in the
   DOM, so anything big belongs in VirtualList. Rows are keyed: an items change keeps the rows whose
   key stays and builds only the items that changed, each in a scope of its own. Columns are live (a
   change builds the header and the rows again), sortable headers set an order the owner applies, and
   sections can fold. */
import {signal, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {each} from '../../core/region.js';
import {Input} from '../../core/input-base.js';
import type {LiveOption} from '../../core/input-base.js';
import {Format} from '../../core/format.js';
import {Tooltip} from '../../core/tooltip.js';
import {icon} from '../display/icon.js';
import {badge} from '../display/badge.js';
import {actionsMenu} from '../actions/actions.js';
import type {Action} from '../actions/actions.js';
import type {IItemsState} from '../../core/widget-like.js';
import {SortOrder} from './sort-order.js';

export interface TableColumn<T> {
  header: string;
  /** A string becomes the cell text; an element is appended as-is. */
  render: (item: T, index: number) => HTMLElement | string;
  width?: string;
  align?: 'left' | 'right' | 'center';
  /** The column's identity: `data-column` on its header and cells, and the name `sort` holds. */
  name?: string;
  /** Its header sets `sort` (see {@link BasicTableOptions.sortGesture}): ascending, descending, unsorted.
   * Needs `name`. */
  sortable?: boolean;
  /** The header's hover text; a sortable header adds how to sort by it. */
  tooltip?: string;
  /** Kept at the start or the end of the row while the table scrolls sideways inside its host. */
  pinned?: 'start' | 'end';
  /** Caps the column (a CSS length): a longer cell ends in an ellipsis and shows its whole text on hover. */
  maxWidth?: string;
}

/** An item that is exactly `{section: 'Title'}`: a header row spanning every column over the rows
 * after it. Not a row of the selection — clicks and the keyboard pass it by. */
export interface TableSection {
  section: string;
}

export interface BasicTableOptions<T> {
  /** Live: a change builds the header and every row again (the rows keep their keys). */
  columns: LiveOption<TableColumn<T>[]>;
  items?: readonly (T | TableSection)[] | ReadonlySignal<readonly (T | TableSection)[]>;
  /** Identity of an item: its row is built once and kept across items changes, the selection
   * follows it, and it is the row's `data-u2-item`. By default a string or number item is its own
   * key, else its string `name` or `id`; anything else is keyed by position. */
  keyOf?: (item: T) => string;
  /** Runs once per row when it is built, in the row's scope: name, class or title the `<tr>`
   * (`tr.dataset.u2Name`, `tr.classList`, `tr.title`), or bind them to signals. */
  row?: (item: T, index: number, tr: HTMLTableRowElement) => void;
  onRowClick?: (item: T, index: number) => void;
  selectable?: boolean;
  headerVisible?: boolean;
  /** The table takes its host's width and the columns share it (`width` pins one), so a long cell
   * ends in an ellipsis rather than widening the table. Off: the table sizes to its content. */
  fixedLayout?: boolean;
  /** The order the sortable headers show and set, as DataTable's ({@link SortOrder}): `'name'` ascending, `'!name'`
   * descending, `''` none. The table does not reorder its items — the owner sorts them. A header writes the next
   * order here, or hands it to `onSort`. */
  sort?: ReadonlySignal<string>;
  /** Takes the next order a header asks for, in place of writing `sort`: an order kept in state of the owner's. */
  onSort?: (sort: string) => void;
  /** What sorts by a header: a click (default, as DataTable) or a double click (as the platform grid).
   * Enter or Space on a focused header sorts either way. */
  sortGesture?: 'click' | 'dblclick';
  /** The folded sections, by title. With it, a section header folds and unfolds the rows under it (a
   * click, or Enter or Space on the focused header) and says so with `aria-expanded`; the rows of a
   * folded section are not rendered. A fold writes the titles here, or hands them to `onCollapse`. */
  collapsed?: ReadonlySignal<readonly string[]>;
  /** Takes the folded titles a section header asks for, in place of writing `collapsed`. */
  onCollapse?: (collapsed: readonly string[]) => void;
  /** A row's actions as its context menu — a right-click, the context-menu key or Shift+F10 — as VirtualList's. */
  contextActions?: (item: T, index: number) => Action[];
  /** The text a section header shows for its title: the title by default. */
  sectionLabel?: (title: string) => string;
  /** A section header counts the rows under it, in a badge after its text. */
  sectionCounts?: boolean;
  /** Row and header height: `'inline'` (default, 22px, as `ui.table`) or `'control'` (28px, a form row). */
  rowHeight?: 'inline' | 'control';
}

interface Entry<T> {
  item: T | TableSection;
  key: string;
  index: number;
}

export class BasicTable<T> extends Control {
  /** Index into the items (sections count) of the selected row, -1 when nothing is selected; only
   * meaningful for a `selectable` table. An items change keeps it on the same key. */
  readonly selectedIndex: Signal<number> = signal(-1);

  private _cols: readonly TableColumn<T>[] = [];
  private readonly _columnsSource: Signal<ReadonlySignal<readonly TableColumn<T>[]>>;
  private readonly _keyOf: ((item: T) => string) | undefined;
  private readonly _rowHook: BasicTableOptions<T>['row'];
  private readonly _onRowClick: ((item: T, index: number) => void) | undefined;
  private readonly _selectable: boolean;
  private readonly _headerVisible: boolean;
  private readonly _sort: ReadonlySignal<string> | undefined;
  private readonly _onSort: (sort: string) => void;
  private readonly _gesture: 'click' | 'dblclick';
  private readonly _collapsed: ReadonlySignal<readonly string[]> | undefined;
  private readonly _onCollapse: (collapsed: readonly string[]) => void;
  private readonly _contextActions: BasicTableOptions<T>['contextActions'];
  private readonly _sectionLabel: ((title: string) => string) | undefined;
  private readonly _sectionCounts: boolean;
  private readonly _body = document.createElement('tbody');
  private _head: HTMLElement | undefined;
  private _colgroup: HTMLElement | undefined;
  private readonly _source = signal<ReadonlySignal<ReadonlyArray<T | TableSection>>>(signal([]));
  private readonly _entries = signal<readonly Entry<T>[]>([]);
  /** The item each kept row shows, so an items change tells a moved row from a changed one. */
  private readonly _shown = new WeakMap<HTMLElement, T | TableSection>();
  /** Bumped by a columns change: a row built for older columns is built again. */
  private _shape = 0;
  private readonly _rowShape = new WeakMap<HTMLElement, number>();
  private _items: readonly (T | TableSection)[] = [];
  /** Indices of the rows under a folded section. */
  private _hidden = new Set<number>();

  constructor(options: BasicTableOptions<T>) {
    super(document.createElement('table'));
    const columns = options.columns;
    this._columnsSource = signal(columns instanceof Signal ? columns : signal(columns));
    this._keyOf = options.keyOf;
    this._rowHook = options.row;
    this._onRowClick = options.onRowClick;
    this._selectable = options.selectable ?? false;
    this._headerVisible = options.headerVisible ?? true;
    const sort = this._sort = options.sort;
    this._onSort = options.onSort ?? ((x) => {
      (sort as Signal<string>).value = x;
    });
    this._gesture = options.sortGesture ?? 'click';
    const collapsed = this._collapsed = options.collapsed;
    this._onCollapse = options.onCollapse ?? ((x) => {
      (collapsed as Signal<readonly string[]>).value = x;
    });
    this._contextActions = options.contextActions;
    this._sectionLabel = options.sectionLabel;
    this._sectionCounts = options.sectionCounts ?? false;

    this.root.classList.add('u2-table');
    this.root.dataset.u2 = 'table';
    if (this._selectable || this._onRowClick)
      this.root.classList.add('u2-table-interactive');
    if (options.fixedLayout)
      this.root.classList.add('u2-table-fixed');
    if (options.rowHeight === 'control')
      this.root.classList.add('u2-table-control');
    if (this._selectable)
      this.root.tabIndex = 0;
    this.root.append(this._body);

    this.scope.listen(this.root, this._gesture, (e) => this._onHeaderGesture(e));
    this.scope.listen(this._body, 'click', (e) => this._onClick(e));
    this.scope.listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this.scope.listen(this._body, 'pointerover', (e) => {
      const cell = (e.target as Element).closest('.u2-table-clip') as HTMLElement | null;
      if (cell !== null)
        Tooltip.titleIfCut(cell);
    });
    if (this._contextActions !== undefined)
      this.scope.listen(this._body, 'contextmenu', (e) => this._onContextMenu(e as MouseEvent));

    if (options.items)
      this.setItems(options.items);
    this.effect(() => {
      const cols = this._columnsSource.value.value;
      untracked(() => this._layout(cols));
    });
    this.effect(() => {
      const items = this._source.value.value;
      this._collapsed?.value;
      untracked(() => this._apply(items));
    });
    this.build(() => each(this._body, this._entries, (e) => e.key, (e) => this._row(e),
      {update: (el, e) => this._patch(el, e)}));
    this.effect(() => this._select(this.selectedIndex.value));
    if (sort !== undefined)
      this.effect(() => this._showSort(sort.value));
  }

  /** The columns shown, in order. A change builds the header and the rows again; the rows keep their
   * keys and the selection stays. */
  get columns(): readonly TableColumn<T>[] {
    return this._cols;
  }

  set columns(x: LiveOption<TableColumn<T>[]>) {
    this._columnsSource.value = x instanceof Signal ? x : signal(x);
  }

  setItems(items: readonly (T | TableSection)[] | ReadonlySignal<readonly (T | TableSection)[]>): void {
    this._source.value = items instanceof Signal ? items : signal(items);
  }

  /** One item per row (sections are not rows): its key, the first cell's text and whether it is
   * selected. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const lead = this.selectedIndex.peek();
    const list = [];
    for (const row of this._body.querySelectorAll<HTMLElement>('.u2-table-row')) {
      const index = Number(row.dataset.index);
      list.push({key: row.dataset.u2Item!, index: list.length,
        label: Control.text(row.querySelector('[data-u2-part="item-label"]') ?? row),
        selected: index === lead, rendered: true});
    }
    return Input.itemsState(list, opts?.all);
  }

  /** A row under a folded section unfolds it first, found by its key or by its first cell's text. */
  async revealItem(keyOrLabel: string): Promise<HTMLElement | null> {
    const found = await super.revealItem(keyOrLabel);
    if (found !== null)
      return found.classList.contains('u2-table-section') ? null : found;
    const collapsed = this._collapsed;
    const wanted = keyOrLabel.trim().toLowerCase();
    const at = this._items.findIndex((item, i) => this._hidden.has(i) &&
      (this._identity(item as T) === keyOrLabel || this._labelOf(item as T, i).toLowerCase() === wanted));
    if (collapsed === undefined || at < 0)
      return null;
    const title = this._sectionOf(at);
    this._onCollapse(collapsed.peek().filter((t) => t !== title));
    return super.revealItem(keyOrLabel);
  }

  /** Whether `x` is a section header: exactly `{section: string}`. */
  static isSection(x: unknown): x is TableSection {
    return typeof x === 'object' && x !== null && typeof (x as TableSection).section === 'string' &&
      Object.keys(x).length === 1;
  }

  /** Folds `title` when it is open, else opens it. */
  toggleSection(title: string): void {
    const collapsed = this._collapsed;
    if (collapsed === undefined)
      return;
    const now = collapsed.peek();
    this._onCollapse(now.includes(title) ? now.filter((t) => t !== title) : [...now, title]);
  }

  /** The text a row's first cell shows — built in a scope of its own, for a row not rendered. */
  private _labelOf(item: T, index: number): string {
    const first = this._cols[0];
    if (first === undefined)
      return '';
    const scope = new Scope(undefined, {label: 'table-label', lazy: true});
    try {
      const content = scope.build(() => first.render(item, index));
      return (typeof content === 'string' ? content : Control.text(content)).trim();
    }
    finally {
      scope.dispose();
    }
  }

  /** New items: the selection follows its key (or keeps its index where items have none), then
   * the rows follow theirs; the rows of a folded section are left out. */
  private _apply(items: readonly (T | TableSection)[]): void {
    const lead = this.selectedIndex.peek();
    const leadKey = lead >= 0 && lead < this._items.length && !BasicTable.isSection(this._items[lead]) ?
      this._identity(this._items[lead] as T) : undefined;
    this._items = items;
    const folded = this._collapsed?.peek() ?? [];
    const entries: Entry<T>[] = [];
    const hidden = new Set<number>();
    let inFolded = false;
    for (const [index, item] of items.entries()) {
      if (BasicTable.isSection(item)) {
        inFolded = folded.includes(item.section);
        entries.push({item, index, key: `section:${item.section}`});
      }
      else if (inFolded)
        hidden.add(index);
      else
        entries.push({item, index, key: this._identity(item) ?? `#${index}`});
    }
    this._hidden = hidden;
    if (lead >= 0) {
      const next = leadKey === undefined ? (lead < items.length && !BasicTable.isSection(items[lead]) ? lead : -1) :
        items.findIndex((item) => !BasicTable.isSection(item) && this._identity(item) === leadKey);
      if (next !== lead)
        this.selectedIndex.value = next;
    }
    this._entries.value = entries;
  }

  private _identity(item: T): string | undefined {
    if (this._keyOf !== undefined)
      return this._keyOf(item);
    if (typeof item === 'string' || typeof item === 'number')
      return String(item);
    if (item === null || typeof item !== 'object')
      return undefined;
    const {name, id} = item as {name?: unknown, id?: unknown};
    return typeof name === 'string' ? name : typeof id === 'string' ? id : undefined;
  }

  /** The header and the colgroup for `cols`; the rows are built again on the next apply. */
  private _layout(cols: readonly TableColumn<T>[]): void {
    this._cols = cols;
    this._shape++;
    this._colgroup?.remove();
    this._colgroup = undefined;
    if (cols.some((c) => c.width)) {
      this._colgroup = this._colgroupOf(cols);
      this.root.prepend(this._colgroup);
    }
    this._head?.remove();
    this._head = undefined;
    if (this._headerVisible) {
      this._head = this._headOf(cols);
      this.root.insertBefore(this._head, this._body);
      if (this._sort !== undefined)
        this._showSort(this._sort.peek());
    }
    if (this._entries.peek().length > 0 || this._items.length > 0)
      this._apply(this._items);
  }

  private _colgroupOf(cols: readonly TableColumn<T>[]): HTMLElement {
    const group = document.createElement('colgroup');
    for (const column of cols) {
      const col = document.createElement('col');
      if (column.width)
        col.style.width = column.width;
      group.append(col);
    }
    return group;
  }

  private _headOf(cols: readonly TableColumn<T>[]): HTMLElement {
    const head = document.createElement('thead');
    const row = document.createElement('tr');
    row.className = 'u2-table-header';
    const how = this._gesture === 'dblclick' ? 'Double-click' : 'Click';
    for (const column of cols) {
      const cell = document.createElement('th');
      cell.dataset.u2 = 'column-header';
      cell.setAttribute('role', 'columnheader');
      cell.textContent = column.header;
      BasicTable._decorate(cell, column);
      const sortable = column.sortable === true && column.name !== undefined && this._sort !== undefined;
      if (sortable) {
        cell.classList.add('u2-table-sortable');
        cell.tabIndex = 0;
      }
      const hint = sortable ? `${how} to sort by ${column.header}` : undefined;
      const tip = [column.tooltip, hint].filter((t) => t).join('\n');
      if (tip !== '')
        BasicTable._tip(cell, tip);
      row.append(cell);
    }
    head.append(row);
    return head;
  }

  /** A row built in its own scope (the region's), so signal-bound cells live as long as it shows
   * its item. */
  private _row(entry: Entry<T>): HTMLElement {
    const row = document.createElement('tr');
    this._shown.set(row, entry.item);
    this._rowShape.set(row, this._shape);
    row.dataset.index = String(entry.index);
    if (BasicTable.isSection(entry.item)) {
      row.className = 'u2-table-section';
      row.dataset.u2Part = 'section';
      const cell = document.createElement('th');
      cell.setAttribute('role', 'rowheader');
      row.append(cell);
      this._fillSection(row, entry);
      return row;
    }
    const item = entry.item as T;
    row.className = 'u2-table-row';
    row.setAttribute('role', 'row');
    for (const column of this._cols) {
      const cell = document.createElement('td');
      BasicTable._decorate(cell, column);
      const content = column.render(item, entry.index);
      if (typeof content === 'string')
        cell.textContent = content;
      else
        cell.append(content);
      if (row.children.length === 0)
        cell.dataset.u2Part = 'item-label';
      row.append(cell);
    }
    this._rowHook?.(item, entry.index, row as HTMLTableRowElement);
    this._mark(row, entry.index === this.selectedIndex.peek());
    return row;
  }

  /** A section header's content: the chevron of a foldable one, its text and its count. */
  private _fillSection(row: HTMLElement, entry: Entry<T>): void {
    const title = (entry.item as TableSection).section;
    const cell = row.children[0] as HTMLTableCellElement;
    cell.colSpan = Math.max(1, this._cols.length);
    const label = document.createElement('span');
    label.className = 'u2-table-section-label';
    label.dataset.u2Part = 'item-label';
    label.textContent = this._sectionLabel?.(title) ?? title;
    const parts: HTMLElement[] = [label];
    if (this._sectionCounts)
      parts.push(badge(String(this._rowsUnder(entry.index))));
    const collapsed = this._collapsed;
    if (collapsed !== undefined) {
      const folded = collapsed.peek().includes(title);
      parts.unshift(icon(folded ? 'angle-right' : 'angle-down', {cls: 'u2-table-chevron'}));
      row.classList.add('u2-table-section-foldable');
      row.setAttribute('aria-expanded', String(!folded));
      cell.tabIndex = 0;
    }
    cell.replaceChildren(...parts);
  }

  private _rowsUnder(index: number): number {
    let n = 0;
    for (let i = index + 1; i < this._items.length && !BasicTable.isSection(this._items[i]); i++)
      n++;
    return n;
  }

  private _sectionOf(index: number): string | undefined {
    for (let i = index; i >= 0; i--) {
      const item = this._items[i];
      if (BasicTable.isSection(item))
        return item.section;
    }
    return undefined;
  }

  /** A kept key: a moved row only learns its index, a section header its new count and fold; a changed
   * item, or a row built for older columns, is built again. */
  private _patch(row: HTMLElement, entry: Entry<T>): boolean {
    if (this._rowShape.get(row) !== this._shape)
      return false;
    const section = BasicTable.isSection(entry.item);
    if (!section && this._shown.get(row) !== entry.item)
      return false;
    this._shown.set(row, entry.item);
    row.dataset.index = String(entry.index);
    if (section)
      this._fillSection(row, entry);
    else
      this._mark(row, entry.index === this.selectedIndex.peek());
    return true;
  }

  private _select(index: number): void {
    for (const row of this._body.querySelectorAll<HTMLElement>('.u2-table-row'))
      this._mark(row, Number(row.dataset.index) === index);
  }

  private _mark(row: HTMLElement, selected: boolean): void {
    row.classList.toggle('u2-table-row-selected', selected);
    if (!this._selectable)
      return;
    if (selected)
      row.setAttribute('aria-selected', 'true');
    else
      row.removeAttribute('aria-selected');
  }

  private _onClick(e: Event): void {
    const target = e.target as Element;
    const section = target.closest('.u2-table-section-foldable') as HTMLElement | null;
    if (section !== null && this._body.contains(section)) {
      this._toggleRow(section);
      return;
    }
    const row = target.closest('.u2-table-row') as HTMLElement | null;
    if (!row)
      return;
    const index = Number(row.dataset.index);
    if (this._selectable)
      this.selectedIndex.value = index;
    this._fire(index);
  }

  private _toggleRow(section: HTMLElement): void {
    const item = this._items[Number(section.dataset.index)];
    if (BasicTable.isSection(item))
      this.toggleSection(item.section);
  }

  private _onHeaderGesture(e: Event): void {
    const head = (e.target as Element).closest('th.u2-table-sortable') as HTMLElement | null;
    if (head === null || this._head === undefined || !this._head.contains(head))
      return;
    e.preventDefault();
    this._sortBy(head);
  }

  private _sortBy(head: HTMLElement): void {
    const sort = this._sort;
    if (sort !== undefined && head.dataset.column !== undefined)
      this._onSort(SortOrder.next(sort.peek(), head.dataset.column));
  }

  /** The row's context menu at the pointer. */
  private _onContextMenu(e: MouseEvent): void {
    if (!this._showContextMenu(e.target as Element, {x: e.clientX, y: e.clientY}))
      return;
    e.preventDefault();
    // the row menu is the only menu, as VirtualList's: an ancestor's contextmenu hook must not add its own
    e.stopPropagation();
  }

  /** The context menu of the row `target` is in, when it has actions: at a point, or under the row. */
  private _showContextMenu(target: Element, at?: {x: number, y: number}): boolean {
    const row = target.closest('.u2-table-row') as HTMLElement | null;
    const index = row === null ? -1 : Number(row.dataset.index);
    const item = this._items[index];
    if (row === null || item === undefined || BasicTable.isSection(item))
      return false;
    const actions = this._contextActions!(item as T, index);
    if (actions.length === 0)
      return false;
    if (this._selectable)
      this.selectedIndex.value = index;
    actionsMenu(actions).show(at ?? {anchor: row});
    return true;
  }

  /** Enter or Space on a header sorts and on a section header folds; the context-menu key or Shift+F10 opens a row's
   * menu; the movement keys go over the rows of a selectable table. */
  private _onKeyDown(e: KeyboardEvent): void {
    const target = e.target as Element;
    if (this._contextActions !== undefined && (e.key === 'ContextMenu' || e.shiftKey && e.key === 'F10')) {
      if (this._showContextMenu(target))
        e.preventDefault();
      return;
    }
    if (e.key === 'Enter' || e.key === ' ') {
      const head = target.closest('th.u2-table-sortable') as HTMLElement | null;
      if (head !== null && this._head?.contains(head)) {
        e.preventDefault();
        this._sortBy(head);
        return;
      }
      const section = target.closest('.u2-table-section-foldable') as HTMLElement | null;
      if (section !== null && this._body.contains(section)) {
        e.preventDefault();
        this._toggleRow(section);
        return;
      }
    }
    if (this._selectable)
      this._onRowKey(e);
  }

  /** The movement keys over the rows shown, sections passed by. */
  private _onRowKey(e: KeyboardEvent): void {
    const rows: number[] = [];
    for (const [i, item] of this._items.entries()) {
      if (!BasicTable.isSection(item) && !this._hidden.has(i))
        rows.push(i);
    }
    if (rows.length === 0)
      return;
    const at = rows.indexOf(this.selectedIndex.peek());
    let next: number;
    switch (e.key) {
      case 'ArrowDown': next = at + 1; break;
      case 'ArrowUp': next = at < 0 ? rows.length - 1 : at - 1; break;
      case 'Home': next = 0; break;
      case 'End': next = rows.length - 1; break;
      case 'Enter':
        e.preventDefault();
        this._fire(this.selectedIndex.peek());
        return;
      default: return;
    }
    e.preventDefault();
    this.selectedIndex.value = rows[Math.max(0, Math.min(rows.length - 1, next))];
  }

  private _fire(index: number): void {
    const item = this._items[index];
    if (this._onRowClick && index >= 0 && index < this._items.length && !BasicTable.isSection(item))
      this._onRowClick(item as T, index);
  }

  private _showSort(sort: string): void {
    if (this._head === undefined)
      return;
    const key = SortOrder.parse(sort);
    for (const cell of this._head.querySelectorAll<HTMLElement>('th.u2-table-sortable')) {
      const sorted = key !== null && cell.dataset.column === key.column;
      cell.setAttribute('aria-sort', !sorted ? 'none' : key.descending ? 'descending' : 'ascending');
      cell.classList.toggle('u2-table-sorted', sorted);
      cell.querySelector('.u2-table-sort')?.remove();
      if (sorted)
        cell.append(icon(key.descending ? 'arrow-down' : 'arrow-up', {cls: 'u2-table-sort'}));
    }
  }

  /** The column's identity, alignment, pin and cap on one of its cells. */
  private static _decorate<X>(cell: HTMLElement, column: TableColumn<X>): void {
    if (column.name !== undefined)
      cell.dataset.column = column.name;
    if (column.align)
      cell.classList.add(`u2-table-align-${column.align}`);
    if (column.pinned)
      cell.classList.add(`u2-table-pin-${column.pinned}`);
    if (column.maxWidth) {
      cell.classList.add('u2-table-clip');
      cell.style.maxWidth = column.maxWidth;
    }
  }

  private static _tip(el: HTMLElement, text: string): void {
    const scope = Scope.ambient;
    if (scope === undefined)
      el.title = text;
    else
      Tooltip.bind(el, text, {owner: scope});
  }
}

/** Two-column key/value table with muted keys — the `ui.tableFromMap` parity helper.
 * Element values are appended as-is, everything else reads as `Format.value` gives it. */
export function tableFromMap(map: Record<string, unknown>): HTMLElement {
  const table = document.createElement('table');
  table.className = 'u2-table u2-table-map';
  table.dataset.u2 = 'table';
  const body = document.createElement('tbody');
  for (const key of Object.keys(map)) {
    const row = document.createElement('tr');
    row.className = 'u2-table-row';
    row.dataset.u2Item = key;
    const name = document.createElement('td');
    name.className = 'u2-table-key';
    name.dataset.u2Part = 'item-label';
    name.textContent = key;
    const value = document.createElement('td');
    const content = map[key];
    if (content instanceof HTMLElement)
      value.append(content);
    else
      value.textContent = Format.value(content);
    row.append(name, value);
    body.append(row);
  }
  table.append(body);
  return table;
}
