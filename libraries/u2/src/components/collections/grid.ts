/* Virtualized grid of fixed-size cells — VirtualList's render-range diff and cell pool in two
   dimensions. The column count follows the root's width, so one grid serves any cell size; the
   index space stays flat (row-major), which is what selection, keys and `keyOf` address. Cells are
   keyed like VirtualList's rows: an items change keeps the cells whose key stays in the window and
   renders only items that changed, each in a scope of its own. */

import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {Input} from '../../core/input-base.js';
import {computed, signal, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import type {IItemsState} from '../../core/widget-like.js';
import {complement, IndexRange} from './list.js';

export interface VirtualGridOptions<T> {
  /** The narrowest a cell may be: as many columns as fit, their cells stretched to share the row. */
  cellWidth: number;
  cellHeight: number;
  /** `cell` is the pooled cell element, cleared before each call — set `title` and classes on it.
   * Runs once per item shown, with the cell's scope ambient: what it binds lives as long as the
   * cell shows the item. */
  render: (item: T, index: number, cell: HTMLElement) => HTMLElement;
  /** Identity of an item: its cell survives items changes, the selection follows it, and it is the
   * cell's `data-u2-item`. By default a string or number item is its own key, else its string
   * `name` or `id`; anything else is keyed by position. */
  keyOf?: (item: T) => string;
  /** The text an item is found by (`revealItem`, the state block): the cell's text or title while
   * it is rendered, else the key, by default. */
  labelOf?: (item: T) => string;
  /** A click, Enter or Space on the selected cell. */
  onActivate?: (item: T, index: number) => void;
}

/** A rendered cell: the pooled element, the item and key it shows, and the scope its render owns. */
interface Cell<T> {
  el: HTMLElement;
  key: string;
  item: T;
  scope: Scope;
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
  private readonly _cells = new Map<number, Cell<T>>();
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

    this.scope.listen(this.root, 'scroll', () => {
      this._measure();
      this._renderVisible();
    });
    this.scope.listen(this.root, 'keydown', (e) => this.onKeyDown(e as KeyboardEvent));
    this.scope.listen(this._content, 'click', (e) => {
      const cell = (e.target as Element).closest('.u2-grid-cell') as HTMLElement | null;
      if (cell)
        this._activate(Number(cell.dataset.index));
    });

    const resize = new ResizeObserver(() => {
      this._measure();
      this._renderVisible();
    });
    resize.observe(this.root);
    this.own(() => resize.disconnect());
    this._measure();

    this.effect(() => {
      const items = this._source.value.value;
      const columns = this.columns.value;
      untracked(() => {
        const lead = this.selectedIndex.peek();
        const selected = lead >= 0 && lead < this._items.length ? this._identity(this._items[lead]) : undefined;
        this._items = items;
        this._layout(columns);
        if (selected !== undefined)
          this.selectedIndex.value = this._items.findIndex((item) => this._identity(item) === selected);
      });
    });

    this.effect(() => {
      const index = this.selectedIndex.value;
      if (index < 0)
        this.root.removeAttribute('aria-activedescendant');
      else
        this.root.setAttribute('aria-activedescendant', this._idPrefix + index);
      for (const [i, cell] of this._cells)
        this._setSelected(cell.el, i === index);
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

  /** The element showing `items[index]`; null while it is not rendered. */
  cellAt(index: number): HTMLElement | null {
    return this._cells.get(index)?.el ?? null;
  }

  /** Renders the cell of `key` again, or every cell on screen — for an item mutated in place. */
  refresh(key?: string): void {
    for (const [index, cell] of this._cells) {
      if (key === undefined || cell.key === key)
        this._fill(cell, index);
    }
  }

  /** Every item (not the rendered window), capped at 500 unless `all`: its key, its label, whether
   * it is selected and rendered. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const lead = this.selectedIndex.peek();
    return Input.itemsState(this._items.map((_, index) => ({key: this._keyAt(index), label: this._label(index),
      index, selected: index === lead, rendered: this._cells.has(index)})), opts?.all);
  }

  /** Scrolls the item — by key, else by label — into the window and resolves to its cell once it
   * is rendered; null when no item matches, a rejection naming the candidates when several do. */
  async revealItem(keyOrLabel: string): Promise<HTMLElement | null> {
    const index = this._find(keyOrLabel);
    if (index < 0)
      return null;
    this.scrollToIndex(index);
    // a frame for the layout to settle — bounded, since a background tab may never paint one
    await new Promise<void>((resolve) => {
      setTimeout(resolve, 100);
      if (typeof requestAnimationFrame === 'function')
        requestAnimationFrame(() => resolve());
    });
    if (!this._cells.has(index))
      this._renderVisible();
    return this.cellAt(index);
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

  /** Outside any effect: the width feeds `columns`, which the items effect reads. */
  private _measure(): void {
    this._width.value = this.root.clientWidth;
  }

  private _activate(index: number): void {
    this.selectedIndex.value = index;
    this._options.onActivate?.(this._items[index], index);
  }

  /** New items or a new column count: the height follows, and the window is fitted to them by key
   * — a cell whose key is still in it keeps its element and scope at its new place, rendered again
   * only when its item changed; the rest are released, and keys new to the window take pooled
   * elements. */
  private _layout(columns: number): void {
    const rows = Math.ceil(this._items.length / columns);
    this._content.style.height = `${rows * this._options.cellHeight}px`;
    if (this.selectedIndex.peek() >= this._items.length)
      this.selectedIndex.value = -1;
    const old = new Map<string, Cell<T>>();
    for (const cell of this._cells.values()) {
      if (old.has(cell.key))
        this._release(cell);
      else
        old.set(cell.key, cell);
    }
    this._cells.clear();
    const range = this._window();
    for (let i = range.start; i < range.end; i++) {
      const cell = old.get(this._keyAt(i));
      if (cell === undefined)
        continue;
      old.delete(cell.key);
      this._cells.set(i, cell);
      if (cell.item === this._items[i])
        this._place(cell.el, i);
      else
        this._fill(cell, i);
    }
    for (const cell of old.values())
      this._release(cell);
    for (let i = range.start; i < range.end; i++) {
      if (!this._cells.has(i))
        this._insertCell(i);
    }
    this._rendered = range;
  }

  private _renderVisible(): void {
    const range = this._window();
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

  /** The cells on screen plus the overscan rows; none while the grid has no width yet (built
   * before it is attached) — the first measure lays it out. */
  private _window(): IndexRange {
    if (this.root.clientWidth === 0)
      return {start: 0, end: 0};
    const {cellHeight} = this._options;
    const cols = this.columns.peek();
    const count = this._items.length;
    const firstRow = Math.floor(this.root.scrollTop / cellHeight);
    const visibleRows = Math.ceil(this.root.clientHeight / cellHeight);
    const clamp = (i: number) => Math.max(0, Math.min(count, i));
    return {start: clamp((firstRow - OVERSCAN_ROWS) * cols),
      end: clamp((firstRow + visibleRows + OVERSCAN_ROWS) * cols)};
  }

  /** The item's own identity, undefined for one keyed by position. */
  private _identity(item: T): string | undefined {
    if (this._options.keyOf !== undefined)
      return this._options.keyOf(item);
    if (typeof item === 'string' || typeof item === 'number')
      return String(item);
    if (item === null || typeof item !== 'object')
      return undefined;
    const {name, id} = item as {name?: unknown, id?: unknown};
    return typeof name === 'string' ? name : typeof id === 'string' ? id : undefined;
  }

  private _keyAt(index: number): string {
    return this._identity(this._items[index]) ?? `#${index}`;
  }

  private _label(index: number): string {
    if (this._options.labelOf !== undefined)
      return this._options.labelOf(this._items[index]);
    const el = this.cellAt(index);
    return (el === null ? '' : Control.text(el) || el.title) || this._keyAt(index);
  }

  /** The index of the item `keyOrLabel` names — its key, else its label (trimmed, any case); -1
   * for none; several label matches throw, naming their keys. */
  private _find(keyOrLabel: string): number {
    const count = this._items.length;
    for (let i = 0; i < count; i++) {
      if (this._keyAt(i) === keyOrLabel)
        return i;
    }
    const wanted = keyOrLabel.trim().toLowerCase();
    const matches: number[] = [];
    for (let i = 0; i < count; i++) {
      if (this._label(i).toLowerCase() === wanted)
        matches.push(i);
    }
    if (matches.length > 1) {
      throw new Error(`u2: ${matches.length} items match "${keyOrLabel}": ` +
        matches.map((i) => this._keyAt(i)).join(', '));
    }
    return matches[0] ?? -1;
  }

  private _insertCell(index: number): void {
    const cell: Cell<T> = {el: this._pool.pop() ?? this._createCell(), key: '', item: this._items[index],
      scope: this._cellScope()};
    this._cells.set(index, cell);
    this._render(cell, index);
    this._content.append(cell.el);
  }

  /** A fresh render: the old render's scope dies BEFORE the new one runs, so its effects never
   * write into the item that replaced it. */
  private _fill(cell: Cell<T>, index: number): void {
    cell.scope.dispose();
    cell.scope = this._cellScope();
    this._render(cell, index);
  }

  private _render(cell: Cell<T>, index: number): void {
    cell.item = this._items[index];
    cell.key = this._keyAt(index);
    const el = cell.el;
    el.className = 'u2-grid-cell';
    el.removeAttribute('title');
    el.dataset.u2Item = cell.key;
    this._place(el, index);
    const content = untracked(() => cell.scope.run(() => this._options.render(cell.item, index, el)));
    el.replaceChildren(...content === undefined ? [] : [content]);
  }

  private _place(el: HTMLElement, index: number): void {
    const cols = this.columns.peek();
    el.style.left = `${(index % cols) * 100 / cols}%`;
    el.style.width = `${100 / cols}%`;
    el.style.top = `${Math.floor(index / cols) * this._options.cellHeight}px`;
    el.id = this._idPrefix + index;
    el.dataset.index = String(index);
    el.setAttribute('aria-posinset', String(index + 1));
    el.setAttribute('aria-setsize', String(this._items.length));
    this._setSelected(el, this.selectedIndex.peek() === index);
  }

  /** Lazy: a cell whose render binds nothing costs no live scope. */
  private _cellScope(): Scope {
    return new Scope(this.scope, {label: 'row', lazy: true});
  }

  private _releaseCell(index: number): void {
    const cell = this._cells.get(index)!;
    this._cells.delete(index);
    this._release(cell);
  }

  private _release(cell: Cell<T>): void {
    cell.scope.dispose();
    cell.el.remove();
    this._pool.push(cell.el);
  }

  private _createCell(): HTMLElement {
    const cell = document.createElement('div');
    cell.setAttribute('role', 'option');
    cell.style.height = `${this._options.cellHeight}px`;
    return cell;
  }

  private _setSelected(cell: HTMLElement, selected: boolean): void {
    cell.classList.toggle('u2-grid-cell-selected', selected);
    cell.setAttribute('aria-selected', String(selected));
  }
}
