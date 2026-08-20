/* Basic (non-virtualized) table — the `ui.table`/`ui.tableFromMap` parity control for small data:
   real table semantics, compact platform density, cells rendered by the caller. Rows all live in the
   DOM, so anything big belongs in VirtualList. */
import {signal, Signal, ReadonlySignal} from '../core/signals.js';
import {Control} from '../core/component.js';
import {Scope} from '../core/scope.js';

export interface TableColumn<T> {
  header: string;
  /** A string becomes the cell text; an element is appended as-is. */
  render: (item: T, index: number) => HTMLElement | string;
  width?: string;
  align?: 'left' | 'right' | 'center';
}

export interface BasicTableOptions<T> {
  columns: TableColumn<T>[];
  items?: T[] | ReadonlySignal<T[]>;
  onRowClick?: (item: T, index: number) => void;
  selectable?: boolean;
  headerVisible?: boolean;
}

export class BasicTable<T> extends Control {
  /** -1 when nothing is selected; only meaningful for a `selectable` table. */
  readonly selectedIndex: Signal<number> = signal(-1);

  private readonly _columns: TableColumn<T>[];
  private readonly _onRowClick: ((item: T, index: number) => void) | undefined;
  private readonly _selectable: boolean;
  private readonly _body = document.createElement('tbody');
  private readonly _source = signal<ReadonlySignal<T[]>>(signal<T[]>([]));
  private _items: T[] = [];
  private _rows: HTMLElement[] = [];
  private _content: Scope | undefined;

  constructor(options: BasicTableOptions<T>) {
    super(document.createElement('table'));
    this._columns = options.columns;
    this._onRowClick = options.onRowClick;
    this._selectable = options.selectable ?? false;

    this.root.classList.add('u2-table');
    this.root.dataset.u2 = 'table';
    if (this._selectable || this._onRowClick)
      this.root.classList.add('u2-table-interactive');
    if (this._selectable)
      this.root.tabIndex = 0;

    if (this._columns.some((c) => c.width))
      this.root.append(this._colgroup());
    if (options.headerVisible ?? true)
      this.root.append(this._head());
    this.root.append(this._body);

    this._listen(this._body, 'click', (e) => this._onClick(e));
    if (this._selectable)
      this._listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this.own(() => this._release());

    if (options.items)
      this.setItems(options.items);
    this.effect(() => {
      this._items = this._source.value.value;
      this._render();
    });
    this.effect(() => this._select(this.selectedIndex.value));
  }

  setItems(items: T[] | ReadonlySignal<T[]>): void {
    this._source.value = Array.isArray(items) ? signal(items) : items;
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

  /** Cell renderers run under a per-render scope, so signal-bound cells are released when the
   * rows they live in are replaced. */
  private _render(): void {
    this._release();
    const scope = new Scope();
    this._content = scope;
    this._rows = Scope.runWith(scope, () => this._items.map((item, i) => this._row(item, i)));
    this._body.replaceChildren(...this._rows);
    if (this.selectedIndex.peek() >= this._items.length)
      this.selectedIndex.value = -1;
    this._select(this.selectedIndex.peek());
  }

  private _row(item: T, index: number): HTMLElement {
    const row = document.createElement('tr');
    row.className = 'u2-table-row';
    row.dataset.index = String(index);
    for (const column of this._columns) {
      const cell = document.createElement('td');
      BasicTable._align(cell, column.align);
      const content = column.render(item, index);
      if (typeof content === 'string')
        cell.textContent = content;
      else
        cell.append(content);
      row.append(cell);
    }
    return row;
  }

  private _select(index: number): void {
    for (let i = 0; i < this._rows.length; i++) {
      const row = this._rows[i];
      const selected = i === index;
      row.classList.toggle('u2-table-row-selected', selected);
      if (!this._selectable)
        continue;
      if (selected)
        row.setAttribute('aria-selected', 'true');
      else
        row.removeAttribute('aria-selected');
    }
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

  private _onKeyDown(e: KeyboardEvent): void {
    const last = this._items.length - 1;
    if (last < 0)
      return;
    const current = this.selectedIndex.peek();
    let next: number;
    switch (e.key) {
      case 'ArrowDown': next = current + 1; break;
      case 'ArrowUp': next = current < 0 ? last : current - 1; break;
      case 'Home': next = 0; break;
      case 'End': next = last; break;
      case 'Enter':
        e.preventDefault();
        this._fire(current);
        return;
      default: return;
    }
    e.preventDefault();
    this.selectedIndex.value = Math.max(0, Math.min(last, next));
  }

  private _fire(index: number): void {
    if (this._onRowClick && index >= 0 && index < this._items.length)
      this._onRowClick(this._items[index], index);
  }

  private _release(): void {
    if (this._content)
      this._content.dispose();
    this._content = undefined;
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
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
    const name = document.createElement('td');
    name.className = 'u2-table-key';
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
