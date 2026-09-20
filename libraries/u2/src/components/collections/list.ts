/* Virtualized rows — port-and-adapt of VS Code's `src/vs/base/browser/ui/list/{listView,rangeMap}.ts`
   (github.com/microsoft/vscode, main @ 5616258b86991889c797f410af2ef6e281ef8406, fetched 2026-08-13), MIT.
   Kept: the render-range diff (relative complement of the previously rendered range against the
   current one, so scroll work is proportional to the viewport) and row recycling through a free pool.
   Dropped: variable item heights — with a fixed height RangeMap collapses to arithmetic — plus
   drag-and-drop, touch, custom scrollables, multi-template renderers, and the mouse-controller and
   accessibility-provider indirection layers.

   `VirtualRows` is that scroller and the selection model over it; `VirtualList` is one row shape on
   top of it (a flex line), `DataTable` another (a grid of pooled cells). */

import {Control} from '../../core/component.js';
import {batch, signal, Signal, ReadonlySignal} from '../../core/signals.js';
import {Action, ActionGroup, actionsMenu} from '../actions/actions.js';

export interface IndexRange {
  start: number;
  end: number;
}

export function intersect(a: IndexRange, b: IndexRange): IndexRange {
  const start = Math.max(a.start, b.start);
  return {start, end: Math.max(start, Math.min(a.end, b.end))};
}

/** The parts of `a` not covered by `b` — at most two ranges, in ascending order. */
export function complement(a: IndexRange, b: IndexRange): IndexRange[] {
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
let rowsCount = 0;
const NONE: ReadonlySet<number> = new Set<number>();

export interface VirtualRowsOptions<T> {
  itemHeight?: number;
  /** Identity of an item, used to keep the selection on the same item across `setItems`. */
  keyOf?: (item: T) => string;
  /** The item's FULL action list: right-click selects the row and opens it as a menu at the
   * cursor. The hover block (`rowActions`) shows the icon-bearing subset of the same list. */
  contextActions?: (item: T, index: number) => Action[];
  /** Submenus under those actions, asked once the right-click has settled the selection: what a
   * multi-selection offers as one ("3 issues"). */
  contextGroups?: (item: T, index: number) => ActionGroup[];
}

/** The chrome a row shape names itself by: `host` is the scroller's class, and `<host>-content`,
 * `<host>-row` and `<host>-row-selected` follow from it. */
interface RowsChrome {
  host: string;
  u2: string;
  role: string;
}

/** The scroller both virtualized collections are: the render-range diff over a fixed row height,
 * a recycling pool whose rows are restored to the attribute set a fresh one carries, and the
 * selection model — a lead plus the multi-selection (F5) — with its gestures and movement keys.
 * A subclass owns the shape of a row and nothing else. */
export abstract class VirtualRows<T> extends Control {
  /** The lead of the selection — where the keyboard and `aria-activedescendant` sit. */
  readonly selectedIndex: Signal<number> = signal(-1);
  /** The full multi-selection (F5): Ctrl/Cmd+click toggles, Shift+click ranges from the anchor
   * (the last plain click); a plain click, a keyboard move, `setItems` and any programmatic
   * {@link selectedIndex} write collapse it to the lead alone. */
  readonly selectedIndices: ReadonlySignal<ReadonlySet<number>>;

  protected readonly itemHeight: number;
  protected readonly content = document.createElement('div');
  protected items: readonly T[] = [];
  /** Identity of an item; a subclass whose items bring their own key replaces it in `setItems`. */
  protected keyOf: ((item: T) => string) | undefined;

  private readonly _chrome: RowsChrome;
  private readonly _selected = signal<ReadonlySet<number>>(NONE);
  private readonly _source = signal<ReadonlySignal<readonly T[]>>(signal<readonly T[]>([]));
  private readonly _rows = new Map<number, HTMLElement>();
  private readonly _pool: HTMLElement[] = [];
  private readonly _baseline = new Map<string, string>();
  private readonly _idPrefix: string;
  private _anchor = -1;
  /** The lead the last multi-gesture wrote; the collapse effect consumes it instead of collapsing.
   * A boolean guard would already be down when an outer `batch()` defers the flush past the write. */
  private _multiWrite: number | undefined;
  private _rendered: IndexRange = {start: 0, end: 0};
  /** The frame a zero-height render asked for, and whether this round already asked for one —
   * see {@link _renderVisible}. */
  private _retry = 0;
  private _retried = false;

  protected constructor(options: VirtualRowsOptions<T>, chrome: RowsChrome) {
    super();
    this.selectedIndices = this._selected;
    this._chrome = chrome;
    this.itemHeight = options.itemHeight ?? 22;
    this.keyOf = options.keyOf;
    this._idPrefix = `${chrome.host}-${++rowsCount}-row-`;

    this.root.classList.add(chrome.host);
    this.root.tabIndex = 0;
    this.root.setAttribute('data-u2', chrome.u2);
    this.root.setAttribute('role', chrome.role);
    this.content.className = `${chrome.host}-content`;
    this.root.append(this.content);

    this.listen(this.root, 'scroll', () => this._renderVisible());
    this.listen(this.root, 'keydown', (e) => this._onKeyDown(e));
    this.listen(this.content, 'click', (e) => {
      const index = this.indexOf(e);
      if (index < 0)
        return;
      if (e.ctrlKey || e.metaKey)
        this._toggle(index);
      else if (e.shiftKey)
        this._range(index);
      else
        this._single(index);
    });
    this.listen(this.content, 'dblclick', (e) => {
      const index = this.indexOf(e);
      if (index >= 0)
        this.activate(index);
    });
    const contextActions = options.contextActions;
    if (contextActions) {
      this.listen(this.content, 'contextmenu', (e) => {
        const index = this.indexOf(e);
        if (index < 0)
          return;
        // the Explorer convention: a right-click inside the selection keeps it, outside collapses.
        // Settled FIRST, so what the menu is built from is the selection the user sees
        if (this._selected.peek().has(index))
          this._write(this._selected.peek(), index);
        else
          this._single(index);
        const actions = contextActions(this.items[index], index);
        const groups = options.contextGroups?.(this.items[index], index);
        if (actions.length === 0 && (groups ?? []).length === 0)
          return;
        e.preventDefault();
        // the row menu is the only menu: an ancestor's contextmenu hook must not add its own
        e.stopPropagation();
        actionsMenu(actions, {groups}).show({x: e.clientX ?? 0, y: e.clientY ?? 0});
      });
    }

    const resize = new ResizeObserver(() => this._renderVisible());
    resize.observe(this.root);
    this.own(() => {
      resize.disconnect();
      if (this._retry !== 0)
        cancelAnimationFrame(this._retry);
    });

    // no items yet, so nothing renders here: a subclass's fields are assigned after `super()`,
    // and it calls `setItems` once they are
    this.effect(() => {
      const keyOf = this.keyOf;
      let key: string | undefined;
      if (keyOf) {
        const selected = this.items[this.selectedIndex.peek()];
        if (selected !== undefined)
          key = keyOf(selected);
      }
      this.items = this._source.value.value;
      this._reset();
      if (keyOf && key !== undefined)
        this.selectedIndex.value = this.items.findIndex((item) => keyOf(item) === key);
    });

    // a lead write that is not one of the multi gestures above — keyboard, setItems' keyed
    // re-select, a consumer's programmatic assignment — collapses the selection to the lead
    this.effect(() => {
      const index = this.selectedIndex.value;
      const gesture = this._multiWrite;
      this._multiWrite = undefined;
      if (gesture === index)
        return;
      this._anchor = index;
      this._selected.value = index < 0 ? NONE : new Set([index]);
    });

    this.effect(() => {
      const index = this.selectedIndex.value;
      const selected = this._selected.value;
      if (index < 0)
        this.root.removeAttribute('aria-activedescendant');
      else
        this.root.setAttribute('aria-activedescendant', this._idPrefix + index);
      for (const [i, row] of this._rows)
        this.setSelected(row, selected.has(i));
    });
  }

  /** Selects the rows the keys name, keeping whichever of them leads today (the first otherwise);
   * a key no longer among the items is dropped, and an empty result clears the selection. What
   * survives a collection being read again — the keys are the identity, the indices are not. */
  selectKeys(keys: Iterable<string>): void {
    const keyOf = this.keyOf;
    if (keyOf === undefined)
      return;
    const wanted = new Set(keys);
    const set = new Set<number>();
    for (const [i, item] of this.items.entries()) {
      if (wanted.has(keyOf(item)))
        set.add(i);
    }
    const lead = this.selectedIndex.peek();
    this._write(set, set.size === 0 ? -1 : set.has(lead) ? lead : [...set][0]);
  }

  /** The keys of the selected rows — the form a selection travels in. */
  selectedKeys(): string[] {
    const keyOf = this.keyOf;
    return keyOf === undefined ? [] :
      [...this._selected.peek()].map((i) => this.items[i]).filter((x) => x !== undefined).map(keyOf);
  }

  setItems(items: readonly T[] | ReadonlySignal<T[]>): void {
    this._source.value = Array.isArray(items) ?
      signal(items as readonly T[]) : items as ReadonlySignal<readonly T[]>;
  }

  scrollToIndex(index: number): void {
    const top = index * this.itemHeight;
    const bottom = top + this.itemHeight - this.viewport;
    if (top < this.root.scrollTop)
      this.root.scrollTop = top;
    else if (bottom > this.root.scrollTop)
      this.root.scrollTop = bottom;
    this._renderVisible();
  }

  get renderedCount(): number {
    return this._rows.size;
  }

  /** Renders the rows on screen again — what a change the items signal did not carry needs
   * (a per-cell verdict raised over the window already rendered). */
  refresh(): void {
    for (const [index, row] of this._rows)
      this.fillRow(row, index);
  }

  /** A pooled row element with whatever fixed structure the shape gives it (its cells, say);
   * called once per element, then recycled. */
  protected abstract createRow(): HTMLElement;

  /** The row's content for `items[index]` — the element is restored to the attribute set a fresh
   * one carries before every call, so a renderer's writes never leak into the next item. */
  protected abstract fillRow(row: HTMLElement, index: number): void;

  /** Enter or Delete on the selected row, and anything else the movement keys leave alone:
   * `true` when it was handled (and `preventDefault` called). */
  protected abstract rowKey(e: KeyboardEvent, index: number): boolean;

  /** A double-click on a row. */
  protected abstract activate(index: number): void;

  /** How the shape marks the selection; an override adds to it, it does not replace it. */
  protected setSelected(row: HTMLElement, selected: boolean): void {
    row.classList.toggle(`${this._chrome.host}-row-selected`, selected);
    row.setAttribute('aria-selected', String(selected));
  }

  /** What the rows have to themselves — the scroller, minus whatever chrome sits over them. */
  protected get viewport(): number {
    return this.root.clientHeight;
  }

  /** The index of the row an event landed in, or -1. */
  protected indexOf(e: Event): number {
    const row = (e.target as Element).closest(`.${this._chrome.host}-row`) as HTMLElement | null;
    return row === null ? -1 : Number(row.dataset.index);
  }

  /** A listener owned by the component scope. */
  protected listen<K extends keyof HTMLElementEventMap>(target: HTMLElement, type: K,
    handler: (e: HTMLElementEventMap[K]) => void): void {
    target.addEventListener(type, handler as EventListener);
    this.own(() => target.removeEventListener(type, handler as EventListener));
  }

  /** The attribute set a freshly created element carries, captured once from the first one; what
   * {@link restore} puts a recycled one back to. */
  protected static capture(el: HTMLElement, into: Map<string, string>): void {
    if (into.size > 0)
      return;
    for (const name of el.getAttributeNames())
      into.set(name, el.getAttribute(name)!);
  }

  protected static restore(el: HTMLElement, baseline: Map<string, string>): void {
    for (const name of el.getAttributeNames()) {
      if (!baseline.has(name))
        el.removeAttribute(name);
    }
    for (const [name, value] of baseline)
      el.setAttribute(name, value);
  }

  private _single(index: number): void {
    this._anchor = index;
    this._write(new Set([index]), index);
  }

  private _toggle(index: number): void {
    const set = new Set(this._selected.peek());
    if (set.has(index)) {
      // a selection must exist: toggling the sole member off is a no-op (the canvas rule)
      if (set.size === 1)
        return;
      set.delete(index);
      const lead = this.selectedIndex.peek();
      this._write(set, lead === index ? [...set][set.size - 1] : lead);
    } else {
      set.add(index);
      this._write(set, index);
    }
  }

  private _range(index: number): void {
    const lead = this.selectedIndex.peek();
    const last = this.items.length - 1;
    const from = this._anchor >= 0 ? this._anchor : lead >= 0 ? lead : index;
    const anchor = Math.max(0, Math.min(last, from));
    const step = index >= anchor ? 1 : -1;
    const set = new Set<number>();
    for (let i = anchor; i !== index + step; i += step)
      set.add(i);
    this._write(set, index);
  }

  /** One gesture writes both signals coherently, past the collapse effect's guard. */
  private _write(set: ReadonlySet<number>, lead: number): void {
    this._multiWrite = lead;
    batch(() => {
      this.selectedIndex.value = lead;
      this._selected.value = set;
    });
  }

  /** The scroller's own keys only: an event from inside a row — a rename box, a cell editor, an
   * action button — keeps its caret, its Home and its End. */
  private _onKeyDown(e: KeyboardEvent): void {
    const last = this.items.length - 1;
    if (last < 0 || e.target !== this.root)
      return;
    const current = this.selectedIndex.value;
    if (current >= 0 && this.rowKey(e, current))
      return;
    const page = Math.max(1, Math.floor(this.viewport / this.itemHeight) - 1);
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
    this.content.style.height = `${this.items.length * this.itemHeight}px`;
    if (this.selectedIndex.peek() >= this.items.length)
      this.selectedIndex.value = -1;
    // even when the lead keeps its index: new items mean the old multi-selection is meaningless,
    // and so is the anchor a later Shift+click would range from
    const lead = this.selectedIndex.peek();
    this._anchor = lead;
    this._selected.value = lead < 0 ? NONE : new Set([lead]);
    this._renderVisible();
  }

  private _renderVisible(): void {
    const count = this.items.length;
    const height = this.viewport;
    // A scroller measured mid-swap answers 0 and the window it yields is the overscan alone —
    // three rows under a status line saying "50" — and since the height afterwards is the height
    // it had BEFORE, no resize ever fires to put the rest back. One frame later the layout has
    // settled: what it renders then is the real window. Once per zero, so a list that is truly
    // off screen settles instead of asking for frames forever.
    if (height <= 0 && count > 0 && !this._retried) {
      this._retried = true;
      this._later();
    } else if (height > 0)
      this._retried = false;
    const first = Math.floor(this.root.scrollTop / this.itemHeight);
    const visible = Math.ceil(height / this.itemHeight);
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

  /** One frame later, once — by then the layout has settled, and a scroller that is still empty
   * is genuinely not on screen: its resize observer is what wakes it. */
  private _later(): void {
    if (this._retry !== 0 || typeof requestAnimationFrame !== 'function')
      return;
    this._retry = requestAnimationFrame(() => {
      this._retry = 0;
      if (!this.scope.isDisposed)
        this._renderVisible();
    });
  }

  private _insertRow(index: number): void {
    let row = this._pool.pop();
    if (row === undefined) {
      row = this.createRow();
      VirtualRows.capture(row, this._baseline);
    } else
      VirtualRows.restore(row, this._baseline);
    row.style.top = `${index * this.itemHeight}px`;
    row.id = this._idPrefix + index;
    row.dataset.index = String(index);
    this.fillRow(row, index);
    this.setSelected(row, this._selected.peek().has(index));
    this._rows.set(index, row);
    this.content.append(row);
  }

  private _releaseRow(index: number): void {
    const row = this._rows.get(index)!;
    this._rows.delete(index);
    row.remove();
    this._pool.push(row);
  }
}

export interface VirtualListOptions<T> extends VirtualRowsOptions<T> {
  /** ARIA role of the pooled row elements; defaults to `option`. */
  rowRole?: string;
  /** `row` is the pooled row element — set per-row attributes on it; it is reset before each call. */
  render: (item: T, index: number, row: HTMLElement) => HTMLElement;
  /** Enter on the selected row, or a double-click on it — its default action. */
  onEnter?: (item: T, index: number) => void;
  /** Delete on the selected row. */
  onDelete?: (item: T, index: number) => void;
}

export class VirtualList<T> extends VirtualRows<T> {
  private readonly _rowRole: string;
  private readonly _renderItem: (item: T, index: number, row: HTMLElement) => HTMLElement;
  private readonly _onEnter: ((item: T, index: number) => void) | undefined;
  private readonly _onDelete: ((item: T, index: number) => void) | undefined;

  constructor(options: VirtualListOptions<T>) {
    super(options, {host: 'u2-list', u2: 'list', role: 'listbox'});
    this._rowRole = options.rowRole ?? 'option';
    this._renderItem = options.render;
    this._onEnter = options.onEnter;
    this._onDelete = options.onDelete;
  }

  protected createRow(): HTMLElement {
    const row = document.createElement('div');
    row.className = 'u2-list-row';
    row.setAttribute('role', this._rowRole);
    row.style.height = `${this.itemHeight}px`;
    return row;
  }

  protected fillRow(row: HTMLElement, index: number): void {
    row.setAttribute('aria-posinset', String(index + 1));
    row.setAttribute('aria-setsize', String(this.items.length));
    row.replaceChildren(this._renderItem(this.items[index], index, row));
  }

  protected rowKey(e: KeyboardEvent, index: number): boolean {
    if (e.key !== 'Enter' && e.key !== 'Delete')
      return false;
    const handler = e.key === 'Enter' ? this._onEnter : this._onDelete;
    if (handler)
      e.preventDefault();
    // Enter and Delete are never movement keys: handled or not, they stop here
    handler?.(this.items[index], index);
    return true;
  }

  protected activate(index: number): void {
    this._onEnter?.(this.items[index], index);
  }

  /** Roving tabindex: only the selected row's action buttons are in the tab order, so Tab from the
   * list lands on its actions and not on every row's. */
  protected setSelected(row: HTMLElement, selected: boolean): void {
    super.setSelected(row, selected);
    for (const b of Array.from(row.querySelectorAll<HTMLElement>('.u2-row-actions button')))
      b.tabIndex = selected ? 0 : -1;
  }
}
