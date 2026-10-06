/* Virtualized grid of fixed-size cells — VirtualList's render-range diff and cell pool in two
   dimensions. The column count follows the root's width, so one grid serves any cell size; the
   index space stays flat (row-major), which is what selection, keys and `keyOf` address. */

import {Control} from '../../core/component.js';
import {computed, signal, Signal, ReadonlySignal} from '../../core/signals.js';
import {complement, IndexRange} from './list.js';

export interface VirtualGridOptions<T> {
  cellWidth: number;
  cellHeight: number;
  /** `cell` is the pooled cell element, cleared before each call — set `title` and classes on it. */
  render: (item: T, index: number, cell: HTMLElement) => HTMLElement;
  /** Identity of an item, used to keep the selection on the same item across `setItems`. */
  keyOf?: (item: T) => string;
  /** A click, Enter or Space on the selected cell. */
  onActivate?: (item: T, index: number) => void;
}

const OVERSCAN_ROWS = 2;
let gridCount = 0;

export class VirtualGrid<T> extends Control {
  readonly selectedIndex: Signal<number> = signal(-1);
  /** Cells per row — what the width holds, at least one. */
  readonly columns: ReadonlySignal<number>;

  private readonly _options: VirtualGridOptions<T>;
  private readonly _width = signal(0);
  private readonly _content = document.createElement('div');
  private readonly _source = signal<ReadonlySignal<T[]>>(signal<T[]>([]));
  private readonly _cells = new Map<number, HTMLElement>();
  private readonly _pool: HTMLElement[] = [];
  private readonly _idPrefix = `u2-grid-${++gridCount}-cell-`;
  private _items: T[] = [];
  private _rendered: IndexRange = {start: 0, end: 0};

  constructor(options: VirtualGridOptions<T>) {
    super();
    this._options = options;
    this.columns = computed(() => Math.max(1, Math.floor(this._width.value / options.cellWidth)));

    this.root.classList.add('u2-grid');
    this.root.tabIndex = 0;
    this.root.setAttribute('data-u2', 'grid');
    this.root.setAttribute('role', 'listbox');
    this._content.className = 'u2-grid-content';
    this.root.append(this._content);

    this._on(this.root, 'scroll', () => {
      this._measure();
      this._renderVisible();
    });
    this._on(this.root, 'keydown', (e) => this.onKeyDown(e));
    this._on(this._content, 'click', (e) => {
      const cell = (e.target as Element).closest('.u2-grid-cell') as HTMLElement | null;
      if (cell)
        this._activate(Number(cell.dataset.index));
    });

    const resize = new ResizeObserver(() => this._measure());
    resize.observe(this.root);
    this.own(() => resize.disconnect());
    this._measure();

    this.effect(() => {
      const keyOf = options.keyOf;
      const selected = keyOf ? this._items[this.selectedIndex.peek()] : undefined;
      this._items = this._source.value.value;
      this._layout();
      if (selected !== undefined)
        this.selectedIndex.value = this._items.findIndex((item) => keyOf!(item) === keyOf!(selected));
    });

    this.effect(() => {
      const index = this.selectedIndex.value;
      if (index < 0)
        this.root.removeAttribute('aria-activedescendant');
      else
        this.root.setAttribute('aria-activedescendant', this._idPrefix + index);
      for (const [i, cell] of this._cells)
        this._setSelected(cell, i === index);
    });
  }

  setItems(items: readonly T[] | ReadonlySignal<T[]>): void {
    this._source.value = Array.isArray(items) ? signal(items as T[]) : items as ReadonlySignal<T[]>;
  }

  scrollToIndex(index: number): void {
    const top = Math.floor(index / this.columns.peek()) * this._options.cellHeight;
    const bottom = top + this._options.cellHeight - this.root.clientHeight;
    if (top < this.root.scrollTop)
      this.root.scrollTop = top;
    else if (bottom > this.root.scrollTop)
      this.root.scrollTop = bottom;
    this._renderVisible();
  }

  get renderedCount(): number {
    return this._cells.size;
  }

  /** The grid's own key handling, public so a host can forward keys typed elsewhere (a search box
   * above the grid): arrows move by cell and row, Page keys by viewport, Enter and Space activate. */
  onKeyDown(e: KeyboardEvent): void {
    const last = this._items.length - 1;
    if (last < 0)
      return;
    const cols = this.columns.peek();
    const page = Math.max(1, Math.floor(this.root.clientHeight / this._options.cellHeight) - 1) * cols;
    const current = this.selectedIndex.peek();
    let next: number;
    switch (e.key) {
      case 'ArrowRight': next = current + 1; break;
      case 'ArrowLeft': next = current < 0 ? last : current - 1; break;
      case 'ArrowDown': next = current + cols; break;
      case 'ArrowUp': next = current < 0 ? last : current - cols; break;
      case 'PageDown': next = current + page; break;
      case 'PageUp': next = current < 0 ? last : current - page; break;
      case 'Home': next = 0; break;
      case 'End': next = last; break;
      case 'Enter': case ' ':
        if (current < 0)
          return;
        e.preventDefault();
        this._activate(current);
        return;
      default: return;
    }
    e.preventDefault();
    this.selectedIndex.value = Math.max(0, Math.min(last, next));
    this.scrollToIndex(this.selectedIndex.value);
  }

  private _on<K extends keyof HTMLElementEventMap>(target: HTMLElement, type: K,
    handler: (e: HTMLElementEventMap[K]) => void): void {
    target.addEventListener(type, handler as EventListener);
    this.own(() => target.removeEventListener(type, handler as EventListener));
  }

  /** Outside any effect: the width feeds `columns`, which the items effect reads. */
  private _measure(): void {
    this._width.value = this.root.clientWidth;
  }

  private _activate(index: number): void {
    this.selectedIndex.value = index;
    this._options.onActivate?.(this._items[index], index);
  }

  private _layout(): void {
    for (const cell of this._cells.values()) {
      cell.remove();
      this._pool.push(cell);
    }
    this._cells.clear();
    this._rendered = {start: 0, end: 0};
    const rows = Math.ceil(this._items.length / this.columns.value);
    this._content.style.height = `${rows * this._options.cellHeight}px`;
    if (this.selectedIndex.peek() >= this._items.length)
      this.selectedIndex.value = -1;
    this._renderVisible();
  }

  private _renderVisible(): void {
    const {cellHeight} = this._options;
    const cols = this.columns.peek();
    const count = this._items.length;
    const firstRow = Math.floor(this.root.scrollTop / cellHeight);
    const visibleRows = Math.ceil(this.root.clientHeight / cellHeight);
    const clamp = (i: number) => Math.max(0, Math.min(count, i));
    const range: IndexRange = {start: clamp((firstRow - OVERSCAN_ROWS) * cols),
      end: clamp((firstRow + visibleRows + OVERSCAN_ROWS) * cols)};

    for (const r of complement(this._rendered, range)) {
      for (let i = r.start; i < r.end; i++)
        this._releaseCell(i);
    }
    for (const r of complement(range, this._rendered)) {
      for (let i = r.start; i < r.end; i++)
        this._insertCell(i);
    }
    this._rendered = range;
  }

  private _insertCell(index: number): void {
    const {cellWidth, cellHeight, render} = this._options;
    const cols = this.columns.peek();
    const cell = this._pool.pop() ?? this._createCell();
    cell.className = 'u2-grid-cell';
    cell.removeAttribute('title');
    cell.style.left = `${(index % cols) * cellWidth}px`;
    cell.style.top = `${Math.floor(index / cols) * cellHeight}px`;
    cell.id = this._idPrefix + index;
    cell.dataset.index = String(index);
    cell.setAttribute('aria-posinset', String(index + 1));
    cell.setAttribute('aria-setsize', String(this._items.length));
    cell.replaceChildren(render(this._items[index], index, cell));
    this._setSelected(cell, this.selectedIndex.peek() === index);
    this._cells.set(index, cell);
    this._content.append(cell);
  }

  private _releaseCell(index: number): void {
    const cell = this._cells.get(index)!;
    this._cells.delete(index);
    cell.remove();
    this._pool.push(cell);
  }

  private _createCell(): HTMLElement {
    const cell = document.createElement('div');
    cell.setAttribute('role', 'option');
    cell.style.width = `${this._options.cellWidth}px`;
    cell.style.height = `${this._options.cellHeight}px`;
    return cell;
  }

  private _setSelected(cell: HTMLElement, selected: boolean): void {
    cell.classList.toggle('u2-grid-cell-selected', selected);
    cell.setAttribute('aria-selected', String(selected));
  }
}
