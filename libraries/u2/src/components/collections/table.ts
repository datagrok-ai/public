/* Basic (non-virtualized) table — the `ui.table`/`ui.tableFromMap` parity control for small data:
   real table semantics, compact platform density, cells rendered by the caller. Rows all live in the
   DOM, so anything big belongs in VirtualList. Rows are keyed: an items change keeps the rows whose
   key stays and builds only the items that changed, each in a scope of its own. */
import {signal, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {each} from '../../core/region.js';
import {Input} from '../../core/input-base.js';
import type {IItemsState} from '../../core/widget-like.js';

export interface TableColumn<T> {
  header: string;
  /** A string becomes the cell text; an element is appended as-is. */
  render: (item: T, index: number) => HTMLElement | string;
  width?: string;
  align?: 'left' | 'right' | 'center';
}

/** An item that is exactly `{section: 'Title'}`: a header row spanning every column over the rows
 * after it. Not a row of the selection — clicks and the keyboard pass it by. */
export interface TableSection {
  section: string;
}

export interface BasicTableOptions<T> {
  columns: TableColumn<T>[];
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

  private readonly _columns: TableColumn<T>[];
  private readonly _keyOf: ((item: T) => string) | undefined;
  private readonly _rowHook: BasicTableOptions<T>['row'];
  private readonly _onRowClick: ((item: T, index: number) => void) | undefined;
  private readonly _selectable: boolean;
  private readonly _body = document.createElement('tbody');
  private readonly _source = signal<ReadonlySignal<ReadonlyArray<T | TableSection>>>(signal([]));
  private readonly _entries = signal<readonly Entry<T>[]>([]);
  /** The item each kept row shows, so an items change tells a moved row from a changed one. */
  private readonly _shown = new WeakMap<HTMLElement, T | TableSection>();
  private _items: readonly (T | TableSection)[] = [];

  constructor(options: BasicTableOptions<T>) {
    super(document.createElement('table'));
    this._columns = options.columns;
    this._keyOf = options.keyOf;
    this._rowHook = options.row;
    this._onRowClick = options.onRowClick;
    this._selectable = options.selectable ?? false;

    this.root.classList.add('u2-table');
    this.root.dataset.u2 = 'table';
    if (this._selectable || this._onRowClick)
      this.root.classList.add('u2-table-interactive');
    if (options.fixedLayout)
      this.root.classList.add('u2-table-fixed');
    if (this._selectable)
      this.root.tabIndex = 0;

    if (this._columns.some((c) => c.width))
      this.root.append(this._colgroup());
    if (options.headerVisible ?? true)
      this.root.append(this._head());
    this.root.append(this._body);

    this.scope.listen(this._body, 'click', (e) => this._onClick(e));
    if (this._selectable)
      this.scope.listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));

    if (options.items)
      this.setItems(options.items);
    this.effect(() => {
      const items = this._source.value.value;
      untracked(() => this._apply(items));
    });
    this.runInScope(() => each(this._body, this._entries, (e) => e.key, (e) => this._row(e),
      {update: (el, e) => this._patch(el, e)}));
    this.effect(() => this._select(this.selectedIndex.value));
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

  async revealItem(keyOrLabel: string): Promise<HTMLElement | null> {
    const found = await super.revealItem(keyOrLabel);
    return found?.classList.contains('u2-table-section') ? null : found;
  }

  /** Whether `x` is a section header: exactly `{section: string}`. */
  static isSection(x: unknown): x is TableSection {
    return typeof x === 'object' && x !== null && typeof (x as TableSection).section === 'string' &&
      Object.keys(x).length === 1;
  }

  /** New items: the selection follows its key (or keeps its index where items have none), then
   * the rows follow theirs. */
  private _apply(items: readonly (T | TableSection)[]): void {
    const lead = this.selectedIndex.peek();
    const leadKey = lead >= 0 && lead < this._items.length && !BasicTable.isSection(this._items[lead]) ?
      this._identity(this._items[lead] as T) : undefined;
    this._items = items;
    const entries = items.map((item, index): Entry<T> => ({item, index,
      key: BasicTable.isSection(item) ? `section:${item.section}` : this._identity(item) ?? `#${index}`}));
    if (lead >= 0) {
      const next = leadKey === undefined ? (lead < items.length && !BasicTable.isSection(items[lead]) ? lead : -1) :
        entries.findIndex((e) => e.key === leadKey);
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

  private _colgroup(): HTMLElement {
    const group = document.createElement('colgroup');
    for (const column of this._columns) {
      const col = document.createElement('col');
      if (column.width)
        col.style.width = column.width;
      group.append(col);
    }
    return group;
  }

  private _head(): HTMLElement {
    const head = document.createElement('thead');
    const row = document.createElement('tr');
    row.className = 'u2-table-header';
    for (const column of this._columns) {
      const cell = document.createElement('th');
      cell.textContent = column.header;
      BasicTable._align(cell, column.align);
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
    row.dataset.index = String(entry.index);
    if (BasicTable.isSection(entry.item)) {
      row.className = 'u2-table-section';
      row.dataset.u2Part = 'section';
      const cell = document.createElement('th');
      cell.colSpan = this._columns.length;
      cell.setAttribute('role', 'rowheader');
      cell.textContent = entry.item.section;
      row.append(cell);
      return row;
    }
    const item = entry.item as T;
    row.className = 'u2-table-row';
    row.setAttribute('role', 'row');
    for (const column of this._columns) {
      const cell = document.createElement('td');
      BasicTable._align(cell, column.align);
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

  /** A kept key: a moved row only learns its index; a changed item is built again. */
  private _patch(row: HTMLElement, entry: Entry<T>): boolean {
    if (this._shown.get(row) !== entry.item)
      return false;
    row.dataset.index = String(entry.index);
    if (!BasicTable.isSection(entry.item))
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
    const row = (e.target as Element).closest('.u2-table-row') as HTMLElement | null;
    if (!row)
      return;
    const index = Number(row.dataset.index);
    if (this._selectable)
      this.selectedIndex.value = index;
    this._fire(index);
  }

  /** The movement keys over the rows, sections passed by. */
  private _onKeyDown(e: KeyboardEvent): void {
    const rows: number[] = [];
    for (const [i, item] of this._items.entries()) {
      if (!BasicTable.isSection(item))
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

  private static _align(cell: HTMLElement, align?: string): void {
    if (align)
      cell.classList.add(`u2-table-align-${align}`);
  }
}

/** Two-column key/value table with muted keys — the `ui.tableFromMap` parity helper.
 * Element values are appended as-is, everything else is stringified. */
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
      value.textContent = content === null || content === undefined ? '' : String(content);
    row.append(name, value);
    body.append(row);
  }
  table.append(body);
  return table;
}
