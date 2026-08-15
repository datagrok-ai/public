/* Virtualized list — port-and-adapt of VS Code's `src/vs/base/browser/ui/list/{listView,rangeMap}.ts`
   (github.com/microsoft/vscode, main @ 5616258b86991889c797f410af2ef6e281ef8406, fetched 2026-08-13), MIT.
   Kept: the render-range diff (relative complement of the previously rendered range against the
   current one, so scroll work is proportional to the viewport) and row recycling through a free pool.
   Dropped: variable item heights — with a fixed height RangeMap collapses to arithmetic — plus
   drag-and-drop, touch, custom scrollables, multi-template renderers, and the mouse-controller and
   accessibility-provider indirection layers. */

import {Component} from '../core/component.js';
import {signal, Signal, ReadonlySignal} from '../core/signals.js';

export interface VirtualListOptions<T> {
  itemHeight?: number;
  /** ARIA role of the pooled row elements; defaults to `option`. */
  rowRole?: string;
  /** Identity of an item, used to keep the selection on the same item across `setItems`. */
  keyOf?: (item: T) => string;
  /** `row` is the pooled row element — set per-row attributes on it; it is reset before each call. */
  render: (item: T, index: number, row: HTMLElement) => HTMLElement;
}

interface IndexRange {
  start: number;
  end: number;
}

function intersect(a: IndexRange, b: IndexRange): IndexRange {
  const start = Math.max(a.start, b.start);
  return {start, end: Math.max(start, Math.min(a.end, b.end))};
}

/** The parts of `a` not covered by `b` — at most two ranges, in ascending order. */
function complement(a: IndexRange, b: IndexRange): IndexRange[] {
  if (a.end <= a.start)
    return [];
  const shared = intersect(a, b);
  if (shared.end <= shared.start)
    return [a];
  const result: IndexRange[] = [];
  if (a.start < shared.start)
    result.push({start: a.start, end: shared.start});
  if (shared.end < a.end)
    result.push({start: shared.end, end: a.end});
  return result;
}

const OVERSCAN = 3;
let listCount = 0;

export class VirtualList<T> extends Component {
  readonly selectedIndex: Signal<number> = signal(-1);

  private readonly _itemHeight: number;
  private readonly _rowRole: string;
  private readonly _keyOf: ((item: T) => string) | undefined;
  private readonly _renderItem: (item: T, index: number, row: HTMLElement) => HTMLElement;
  private readonly _content = document.createElement('div');
  private readonly _source = signal<ReadonlySignal<T[]>>(signal<T[]>([]));
  private readonly _rows = new Map<number, HTMLElement>();
  private readonly _pool: HTMLElement[] = [];
  private readonly _baseline = new Map<string, string>();
  private readonly _idPrefix = `u2-list-${++listCount}-row-`;
  private _items: T[] = [];
  private _rendered: IndexRange = {start: 0, end: 0};

  constructor(options: VirtualListOptions<T>) {
    super();
    this._itemHeight = options.itemHeight ?? 22;
    this._rowRole = options.rowRole ?? 'option';
    this._keyOf = options.keyOf;
    this._renderItem = options.render;

    this.root.classList.add('u2-list');
    this.root.tabIndex = 0;
    this.root.setAttribute('data-u2', 'list');
    this.root.setAttribute('role', 'listbox');
    this._content.className = 'u2-list-content';
    this.root.append(this._content);

    this._on(this.root, 'scroll', () => this._renderVisible());
    this._on(this.root, 'keydown', (e) => this._onKeyDown(e));
    this._on(this._content, 'click', (e) => {
      const row = (e.target as Element).closest('.u2-list-row') as HTMLElement | null;
      if (row)
        this.selectedIndex.value = Number(row.dataset.index);
    });

    const resize = new ResizeObserver(() => this._renderVisible());
    resize.observe(this.root);
    this.own(() => resize.disconnect());

    this.effect(() => {
      const keyOf = this._keyOf;
      let key: string | undefined;
      if (keyOf) {
        const selected = this._items[this.selectedIndex.peek()];
        if (selected !== undefined)
          key = keyOf(selected);
      }
      this._items = this._source.value.value;
      this._reset();
      if (keyOf && key !== undefined)
        this.selectedIndex.value = this._items.findIndex((item) => keyOf(item) === key);
    });

    this.effect(() => {
      const index = this.selectedIndex.value;
      if (index < 0)
        this.root.removeAttribute('aria-activedescendant');
      else
        this.root.setAttribute('aria-activedescendant', this._idPrefix + index);
      for (const [i, row] of this._rows)
        this._setSelected(row, i === index);
    });
  }

  setItems(items: T[] | ReadonlySignal<T[]>): void {
    this._source.value = Array.isArray(items) ? signal(items) : items;
  }

  scrollToIndex(index: number): void {
    const top = index * this._itemHeight;
    const bottom = top + this._itemHeight - this.root.clientHeight;
    if (top < this.root.scrollTop)
      this.root.scrollTop = top;
    else if (bottom > this.root.scrollTop)
      this.root.scrollTop = bottom;
    this._renderVisible();
  }

  get renderedCount(): number {
    return this._rows.size;
  }

  private _on<K extends keyof HTMLElementEventMap>(target: HTMLElement, type: K,
    handler: (e: HTMLElementEventMap[K]) => void): void {
    target.addEventListener(type, handler as EventListener);
    this.own(() => target.removeEventListener(type, handler as EventListener));
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const last = this._items.length - 1;
    if (last < 0)
      return;
    const page = Math.max(1, Math.floor(this.root.clientHeight / this._itemHeight) - 1);
    const current = this.selectedIndex.value;
    let next: number;
    switch (e.key) {
      case 'ArrowDown': next = current + 1; break;
      case 'ArrowUp': next = current < 0 ? last : current - 1; break;
      case 'PageDown': next = current + page; break;
      case 'PageUp': next = current < 0 ? last : current - page; break;
      case 'Home': next = 0; break;
      case 'End': next = last; break;
      default: return;
    }
    e.preventDefault();
    this.selectedIndex.value = Math.max(0, Math.min(last, next));
    this.scrollToIndex(this.selectedIndex.value);
  }

  private _reset(): void {
    for (const row of this._rows.values()) {
      row.remove();
      this._pool.push(row);
    }
    this._rows.clear();
    this._rendered = {start: 0, end: 0};
    this._content.style.height = `${this._items.length * this._itemHeight}px`;
    if (this.selectedIndex.peek() >= this._items.length)
      this.selectedIndex.value = -1;
    this._renderVisible();
  }

  private _renderVisible(): void {
    const count = this._items.length;
    const first = Math.floor(this.root.scrollTop / this._itemHeight);
    const visible = Math.ceil(this.root.clientHeight / this._itemHeight);
    const clamp = (i: number) => Math.max(0, Math.min(count, i));
    const range: IndexRange = {start: clamp(first - OVERSCAN), end: clamp(first + visible + OVERSCAN)};

    for (const r of complement(this._rendered, range)) {
      for (let i = r.start; i < r.end; i++)
        this._releaseRow(i);
    }
    for (const r of complement(range, this._rendered)) {
      for (let i = r.start; i < r.end; i++)
        this._insertRow(i);
    }
    this._rendered = range;
  }

  private _insertRow(index: number): void {
    const row = this._pool.pop() ?? this._createRow();
    this._resetRow(row);
    row.style.top = `${index * this._itemHeight}px`;
    row.id = this._idPrefix + index;
    row.dataset.index = String(index);
    row.setAttribute('aria-posinset', String(index + 1));
    row.setAttribute('aria-setsize', String(this._items.length));
    row.replaceChildren(this._renderItem(this._items[index], index, row));
    this._setSelected(row, index === this.selectedIndex.peek());
    this._rows.set(index, row);
    this._content.append(row);
  }

  private _releaseRow(index: number): void {
    const row = this._rows.get(index)!;
    this._rows.delete(index);
    row.remove();
    this._pool.push(row);
  }

  private _createRow(): HTMLElement {
    const row = document.createElement('div');
    row.className = 'u2-list-row';
    row.setAttribute('role', this._rowRole);
    row.style.height = `${this._itemHeight}px`;
    if (this._baseline.size === 0) {
      for (const name of row.getAttributeNames())
        this._baseline.set(name, row.getAttribute(name)!);
    }
    return row;
  }

  /** Renderers write attributes on the row they are handed, so a recycled row is restored to the
   * attribute set a freshly created one has before it is rendered again. */
  private _resetRow(row: HTMLElement): void {
    for (const name of row.getAttributeNames()) {
      if (!this._baseline.has(name))
        row.removeAttribute(name);
    }
    for (const [name, value] of this._baseline)
      row.setAttribute(name, value);
  }

  private _setSelected(row: HTMLElement, selected: boolean): void {
    row.classList.toggle('u2-list-row-selected', selected);
    row.setAttribute('aria-selected', String(selected));
  }
}
