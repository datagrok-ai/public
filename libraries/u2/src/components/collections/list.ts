/* Virtualized rows — port-and-adapt of VS Code's `src/vs/base/browser/ui/list/{listView,rangeMap}.ts`
   (github.com/microsoft/vscode, main @ 5616258b86991889c797f410af2ef6e281ef8406, fetched 2026-08-13), MIT.
   Kept: the render-range diff (relative complement of the previously rendered range against the
   current one, so scroll work is proportional to the viewport) and row recycling through a free pool.
   Dropped: variable item heights — with a fixed height RangeMap collapses to arithmetic — plus
   drag-and-drop, touch, custom scrollables, multi-template renderers, and the mouse-controller and
   accessibility-provider indirection layers.

   `VirtualRows` is that scroller and the selection model over it; `VirtualList` is one row shape on
   top of it (a flex line), `DataTable` another (a grid of pooled cells). Rows are keyed: an items
   change keeps the rows whose key stays in the window and renders only what changed, and every
   rendered row owns a scope its render runs in. */

import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {batch, signal, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import {Action, ActionGroup, actionsMenu} from '../actions/actions.js';
import {ITEMS_CAP} from '../../core/widget-like.js';
import type {IItemState, IItemsState} from '../../core/widget-like.js';

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

/** A rendered row: the pooled element, the item and key it shows, and the scope its render owns. */
interface Row<T> {
  el: HTMLElement;
  key: string;
  item: T;
  scope: Scope;
}

export interface VirtualRowsOptions<T> {
  itemHeight?: number;
  /** Identity of an item: its row survives items changes, its selection too, and it is the row's
   * `data-u2-item`. By default a string or number item is its own key, else its string `name` or
   * `id`; anything else is keyed by position. */
  keyOf?: (item: T) => string;
  /** The text an item is found by (`revealItem`, the state block) while it is not rendered; the
   * rendered row's text, else the key, by default. */
  labelOf?: (item: T) => string;
  /** Items are live views (a frame row proxy) whose content changes without a new object: every
   * items change renders the window again. Off, a row whose key and item are unchanged keeps its
   * DOM, and an item mutated in place needs {@link VirtualRows.refresh}. Default: off for
   * `VirtualList`; `DataTable` declares its own, decided by what feeds it. */
  liveItems?: boolean;
  /** The item's FULL action list: right-click selects the row and opens it as a menu at the
   * cursor. The hover block (`rowActions`) shows the icon-bearing subset of the same list. */
  contextActions?: (item: T, index: number) => Action[];
  /** Submenus under those actions, asked once the right-click has settled the selection: what a
   * multi-selection offers as one ("3 issues"). */
  contextGroups?: (item: T, index: number) => ActionGroup[];
}

/** The chrome a row shape names itself by: `host` is the scroller's class, and `<host>-content`,
 * `<host>-row`, `<host>-row-selected` and `<host>-row-current` follow from it. */
interface RowsChrome {
  host: string;
  u2: string;
  role: string;
}

/** The scroller both virtualized collections are: the render-range diff over a fixed row height,
 * keyed rows — each with its own scope — over a recycling pool whose elements are restored to the
 * attribute set a fresh one carries, and the selection model — a lead plus the multi-selection
 * (F5) — with its gestures and movement keys. A subclass owns the shape of a row and nothing else. */
export abstract class VirtualRows<T> extends Control {
  /** The lead of the selection — where the keyboard and `aria-activedescendant` sit. */
  readonly selectedIndex: Signal<number> = signal(-1);
  /** The full multi-selection (F5): Ctrl/Cmd+click toggles, Shift+click ranges from the anchor
   * (the last plain click); a plain click, a keyboard move and any programmatic
   * {@link selectedIndex} write collapse it to the lead alone. An items change keeps every selected
   * item still present when items have keys (see {@link VirtualRowsOptions.keyOf}), else the lead. */
  readonly selectedIndices: ReadonlySignal<ReadonlySet<number>>;

  readonly itemHeight: number;
  protected readonly content = document.createElement('div');
  protected items: readonly T[] = [];
  /** Identity of an item; a subclass whose items bring their own key replaces it in `setItems`. */
  protected keyOf: ((item: T) => string) | undefined;

  private readonly _chrome: RowsChrome;
  private readonly _labelOf: ((item: T) => string) | undefined;
  private readonly _live: boolean;
  private readonly _selected = signal<ReadonlySet<number>>(NONE);
  private readonly _source = signal<ReadonlySignal<readonly T[]>>(signal<readonly T[]>([]));
  private readonly _rows = new Map<number, Row<T>>();
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
    this._labelOf = options.labelOf;
    this._live = options.liveItems ?? true;
    this._idPrefix = `${chrome.host}-${++rowsCount}-row-`;

    this.root.classList.add(chrome.host);
    this.root.tabIndex = 0;
    this.root.setAttribute('data-u2', chrome.u2);
    this.root.setAttribute('role', chrome.role);
    this.content.className = `${chrome.host}-content`;
    this.root.append(this.content);

    const scope = this.scope;
    scope.listen(this.root, 'scroll', () => this._renderVisible());
    scope.listen(this.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    scope.listen(this.content, 'click', (e) => {
      const index = this.indexOf(e);
      if (index < 0)
        return;
      const mouse = e as MouseEvent;
      if (mouse.ctrlKey || mouse.metaKey)
        this._toggle(index);
      else if (mouse.shiftKey)
        this._range(index);
      else
        this._single(index);
    });
    scope.listen(this.content, 'dblclick', (e) => {
      const index = this.indexOf(e);
      if (index >= 0)
        this.activate(index);
    });
    const contextActions = options.contextActions;
    if (contextActions) {
      scope.listen(this.content, 'contextmenu', (event) => {
        const e = event as MouseEvent;
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
      const items = this._source.value.value;
      untracked(() => this._apply(items));
    });

    // a lead write that is not one of the multi gestures above — keyboard, a consumer's
    // programmatic assignment — collapses the selection to the lead
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
        this._mark(row.el, i, selected, index);
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

  /** The element showing `items[index]`; null while it is not rendered. */
  rowAt(index: number): HTMLElement | null {
    return this._rows.get(index)?.el ?? null;
  }

  /** Renders the row of `key` again, or every row on screen — what a change the items signal did
   * not carry needs: an item mutated in place, a per-cell verdict raised over the rendered window. */
  refresh(key?: string): void {
    for (const [index, row] of this._rows) {
      if (key === undefined || row.key === key)
        this._rerender(row, index);
    }
  }

  /** Every item (not the rendered window), capped at 500 unless `all`: its key, its label, whether
   * it is selected and rendered. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const count = this.items.length;
    const end = opts?.all ? count : Math.min(count, ITEMS_CAP);
    const selected = this._selected.peek();
    const list: IItemState[] = [];
    for (let i = 0; i < end; i++) {
      const el = this.rowAt(i);
      list.push({key: this.keyAt(i), label: this._label(i), index: i, selected: selected.has(i),
        disabled: el?.getAttribute('aria-disabled') === 'true', rendered: el !== null});
    }
    return {count, list, truncated: end < count};
  }

  /** Scrolls the item — by key, else by label — into the window and resolves to its row once it
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
    if (!this._rows.has(index))
      this._renderVisible();
    return this.rowAt(index);
  }

  /** The base state, with the selection as the keys of every selected item, in selection order. */
  getState(): ReturnType<Control['getState']> {
    return {...super.getState(), selection: [...this._selected.peek()].map((i) => this.keyAt(i))};
  }

  /** The key of `items[index]` — see {@link VirtualRowsOptions.keyOf}. */
  protected keyAt(index: number): string {
    return this._identity(this.items[index]) ?? `#${index}`;
  }

  /** A pooled row element with whatever fixed structure the shape gives it (its cells, say);
   * called once per element, then recycled. */
  protected abstract createRow(): HTMLElement;

  /** The row's content for `items[index]` — the element is restored to the attribute set a fresh
   * one carries before every call, so a renderer's writes never leak into the next item. Runs with
   * the row's own scope ambient: what it binds or builds lives as long as the row shows the item. */
  protected abstract fillRow(row: HTMLElement, index: number): void;

  /** A kept row whose item was replaced: patches it in place, `true` when it did — else the row
   * is rendered again. Runs in the row's scope. */
  protected patchRow(_row: HTMLElement, _index: number): boolean {
    return false;
  }

  /** A kept row after an items change, at `index` now: writes what its position says, `true` when
   * that is all it needs — else the row is rendered again (the default, for a shape whose content
   * may depend on the index). */
  protected placeRow(_row: HTMLElement, _index: number): boolean {
    return false;
  }

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

  /** New items: the selection follows its keys (or keeps the lead's index where items have none),
   * then the window is fitted to them by key. */
  private _apply(items: readonly T[]): void {
    const previous = this.items;
    const lead = this.selectedIndex.peek();
    const leadKey = lead >= 0 && lead < previous.length ? this._identity(previous[lead]) : undefined;
    const keys: string[] = [];
    let keyed = lead < 0 || leadKey !== undefined;
    for (const i of this._selected.peek()) {
      const key = keyed && i < previous.length ? this._identity(previous[i]) : undefined;
      if (key === undefined) {
        keyed = false;
        break;
      }
      keys.push(key);
    }
    this.items = items;
    this.content.style.height = `${items.length * this.itemHeight}px`;

    let set: Set<number>;
    let next: number;
    if (keyed && keys.length > 0) {
      const at = this._indicesOf(new Set(keys));
      set = new Set(keys.map((key) => at.get(key)).filter((i): i is number => i !== undefined));
      next = leadKey === undefined ? -1 : at.get(leadKey) ?? -1;
      if (next < 0 && set.size > 0)
        next = [...set][0];
    } else {
      next = lead < items.length ? lead : -1;
      set = next < 0 ? new Set() : new Set([next]);
    }
    // the anchor a later Shift+click ranges from follows the lead, never a row that may be gone
    this._anchor = next;
    this._write(set, next);
    this._reconcile();
  }

  /** The first index of every key in `wanted`, found by one walk that stops once all are. */
  private _indicesOf(wanted: ReadonlySet<string>): Map<string, number> {
    const at = new Map<string, number>();
    for (let i = 0; i < this.items.length && at.size < wanted.size; i++) {
      const key = this._identity(this.items[i]);
      if (key !== undefined && wanted.has(key) && !at.has(key))
        at.set(key, i);
    }
    return at;
  }

  /** Fits the window to new items by key: a row whose key is still in it keeps its element and
   * scope at its new place, rendered again only when its item changed (or items are live); the
   * rest are released, and the keys new to the window take pooled elements. */
  private _reconcile(): void {
    const old = new Map<string, Row<T>>();
    for (const row of this._rows.values()) {
      if (old.has(row.key))
        this._release(row);
      else
        old.set(row.key, row);
    }
    this._rows.clear();
    const range = this._window();
    for (let i = range.start; i < range.end; i++) {
      const row = old.get(this.keyAt(i));
      if (row === undefined)
        continue;
      old.delete(row.key);
      this._rows.set(i, row);
      this._keep(row, i);
    }
    for (const row of old.values())
      this._release(row);
    for (let i = range.start; i < range.end; i++) {
      if (!this._rows.has(i))
        this._insertRow(i);
    }
    this._rendered = range;
    this._order();
  }

  /** Puts the rows in index order in the DOM — what a screen reader and an ordinal ("the second
   * item") read — moving only the rows out of place, so a kept row keeps its focus. */
  private _order(): void {
    let next = this.content.firstChild;
    for (let i = this._rendered.start; i < this._rendered.end; i++) {
      const el = this._rows.get(i)?.el;
      if (el === undefined)
        continue;
      if (el === next)
        next = el.nextSibling;
      else
        this.content.insertBefore(el, next);
    }
  }

  /** A row kept across an items change, now at `index`: moved, then patched or rendered again
   * when its item is not the one it shows. */
  private _keep(row: Row<T>, index: number): void {
    const same = !this._live && this.items[index] === row.item;
    row.item = this.items[index];
    this._place(row.el, index);
    const kept = this._hook(row, () => this.placeRow(row.el, index)) &&
      (same || !this._live && this._hook(row, () => this.patchRow(row.el, index)));
    if (kept)
      this._mark(row.el, index);
    else
      this._rerender(row, index);
  }

  /** Runs a shape hook in the row's scope, untracked; a throw is reported and answers false. */
  private _hook(row: Row<T>, fn: () => boolean): boolean {
    return untracked(() => row.scope.run(fn)) === true;
  }

  private _renderVisible(): void {
    const range = this._window();
    for (const r of complement(this._rendered, range)) {
      for (let i = r.start; i < r.end; i++)
        this._releaseRow(i);
    }
    const before = this._rendered;
    for (const r of complement(range, before)) {
      for (let i = r.start; i < r.end; i++)
        this._insertRow(i);
    }
    this._rendered = range;
    // rows scrolled in below are appended in place already; above, they go before the kept ones
    if (range.start < before.start && before.end > before.start)
      this._order();
  }

  /** The indices on screen plus the overscan. */
  private _window(): IndexRange {
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
    return {start: clamp(first - OVERSCAN), end: clamp(first + visible + OVERSCAN)};
  }

  /** The item's own identity, undefined for one keyed by position. */
  private _identity(item: T): string | undefined {
    if (this.keyOf !== undefined)
      return this.keyOf(item);
    if (typeof item === 'string' || typeof item === 'number')
      return String(item);
    if (item === null || typeof item !== 'object')
      return undefined;
    const {name, id} = item as {name?: unknown, id?: unknown};
    return typeof name === 'string' ? name : typeof id === 'string' ? id : undefined;
  }

  /** What the user reads on the item: `labelOf`, else the rendered row's text, else the key. */
  private _label(index: number): string {
    if (this._labelOf !== undefined)
      return this._labelOf(this.items[index]);
    const el = this.rowAt(index);
    return el === null ? this.keyAt(index) : Control.text(el.querySelector('[data-u2-part="item-label"]') ?? el);
  }

  /** The index of the item `keyOrLabel` names — its key, else its label (trimmed, any case); -1
   * for none; several label matches throw, naming their keys. */
  private _find(keyOrLabel: string): number {
    const count = this.items.length;
    for (let i = 0; i < count; i++) {
      if (this.keyAt(i) === keyOrLabel)
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
        matches.map((i) => this.keyAt(i)).join(', '));
    }
    return matches[0] ?? -1;
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
    let el = this._pool.pop();
    if (el === undefined) {
      el = this.createRow();
      VirtualRows.capture(el, this._baseline);
    } else
      VirtualRows.restore(el, this._baseline);
    const row: Row<T> = {el, key: this.keyAt(index), item: this.items[index], scope: this._rowScope()};
    this._rows.set(index, row);
    this._fill(row, index);
    this.content.append(el);
  }

  /** A fresh render of a kept row: the old render's scope dies BEFORE the new one runs, so its
   * effects never write into the item that replaced it. */
  private _rerender(row: Row<T>, index: number): void {
    row.scope.dispose();
    row.scope = this._rowScope();
    row.item = this.items[index];
    VirtualRows.restore(row.el, this._baseline);
    this._fill(row, index);
  }

  private _fill(row: Row<T>, index: number): void {
    this._place(row.el, index);
    row.el.dataset.u2Item = row.key;
    untracked(() => row.scope.run(() => this.fillRow(row.el, index)));
    this._mark(row.el, index);
  }

  /** The selection, and apart from it the lead: inside a multi-selection the row Enter and the
   * keys act on has to be found. */
  private _mark(el: HTMLElement, index: number, selected = this._selected.peek(),
    lead = this.selectedIndex.peek()): void {
    this.setSelected(el, selected.has(index));
    el.classList.toggle(`${this._chrome.host}-row-current`, index === lead);
  }

  private _place(el: HTMLElement, index: number): void {
    el.style.top = `${index * this.itemHeight}px`;
    el.id = this._idPrefix + index;
    el.dataset.index = String(index);
  }

  /** Lazy: a row whose render binds nothing costs no live scope. */
  private _rowScope(): Scope {
    return new Scope(this.scope, {label: 'row', lazy: true});
  }

  private _releaseRow(index: number): void {
    const row = this._rows.get(index)!;
    this._rows.delete(index);
    this._release(row);
  }

  private _release(row: Row<T>): void {
    row.scope.dispose();
    row.el.remove();
    this._pool.push(row.el);
  }
}

export interface VirtualListOptions<T> extends VirtualRowsOptions<T> {
  /** ARIA role of the pooled row elements; defaults to `option`. */
  rowRole?: string;
  /** `row` is the pooled row element — set per-row attributes on it; it is reset before each call.
   * Runs once per item shown, with the row's scope ambient: `span(signal)` follows the signal and
   * stops with the row. A row is rendered again only when its item changes (a new object under the
   * same key) — an item mutated in place needs `refresh(key)`. */
  render: (item: T, index: number, row: HTMLElement) => HTMLElement;
  /** A kept row's item was replaced by another under the same key: patch `el` (what `render`
   * returned) for `item`, in the row's scope; answering `false` renders the row again instead. */
  update?: (item: T, index: number, row: HTMLElement, el: HTMLElement) => boolean | void;
  /** Enter on the selected row, or a double-click on it — its default action. */
  onEnter?: (item: T, index: number) => void;
  /** Delete on the selected row. */
  onDelete?: (item: T, index: number) => void;
}

export class VirtualList<T> extends VirtualRows<T> {
  private readonly _rowRole: string;
  private readonly _renderItem: (item: T, index: number, row: HTMLElement) => HTMLElement;
  private readonly _update: VirtualListOptions<T>['update'];
  private readonly _onEnter: ((item: T, index: number) => void) | undefined;
  private readonly _onDelete: ((item: T, index: number) => void) | undefined;

  constructor(options: VirtualListOptions<T>) {
    super({...options, liveItems: options.liveItems ?? false}, {host: 'u2-list', u2: 'list', role: 'listbox'});
    this._rowRole = options.rowRole ?? 'option';
    this._renderItem = options.render;
    this._update = options.update;
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
    this.placeRow(row, index);
    row.replaceChildren(this._renderItem(this.items[index], index, row));
  }

  protected placeRow(row: HTMLElement, index: number): boolean {
    row.setAttribute('aria-posinset', String(index + 1));
    row.setAttribute('aria-setsize', String(this.items.length));
    return true;
  }

  protected patchRow(row: HTMLElement, index: number): boolean {
    const el = row.children[0] as HTMLElement | undefined;
    return this._update !== undefined && el !== undefined &&
      this._update(this.items[index], index, row, el) !== false;
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
