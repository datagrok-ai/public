import {signal, untracked, Signal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {Control} from '../../core/component.js';
import {div, effect, on, span, Text} from '../../core/elements.js';
import {each} from '../../core/region.js';
import {Reason} from '../../core/reason.js';
import {Ownership} from '../../core/owner.js';
import {Input} from '../../core/input-base.js';
import {Tooltip} from '../../core/tooltip.js';
import type {IItemsState} from '../../core/widget-like.js';
import {Splitter} from '../containers/splitter.js';

export interface MasterDetailRow {
  label: Text;
  /** A check before the label, two-way: the item is on or off (a function to run, a library to load). */
  check?: Signal<boolean>;
  /** The check's tooltip: what checking the item does. */
  checkTooltip?: string;
  /** Grey text after the label: what the detail holds ("3 selected"). */
  hint?: Text;
}

export interface MasterDetailOptions<T> {
  items: readonly T[] | ReadonlySignal<readonly T[]>;
  keyOf: (item: T) => string;
  row: (item: T) => MasterDetailRow;
  /** The detail of an item: built the first time the item is current and kept while it stays, so what was typed
   * in it survives a move to another item. */
  detail: (item: T) => Control | HTMLElement;
  /** Why an item's detail is off: it stays in view, dimmed and out of reach, the reason over it; `null` while on. */
  disabled?: (item: T) => Reason;
  /** The current item's key, adopted: the first item's by default. */
  current?: Signal<string | null>;
  /** The list's share of the width, 0–1; 0.33 by default. The user drags the sash between them. */
  split?: number;
  /** Shown as the detail while there are no items. */
  emptyText?: string;
  name?: string;
  /** The list's accessible name. */
  label?: string;
}

interface Row {
  el: HTMLElement;
  label: HTMLElement;
  check: Signal<boolean> | undefined;
}

interface Pane {
  el: HTMLElement;
  scope: Scope;
}

/** A list beside the detail of its current item, split by a sash: rows with a check and a hint, the arrow keys
 * moving the current one and Space checking it, each detail built when first shown and kept. For lists up to a few
 * hundred items: it renders every row. A longer list is a `VirtualList` (`new VirtualList({keyOf, render})`, then
 * `setItems(items)`; its `selectedIndex` is the current row) beside a detail of your own, in a `Splitter`. */
export class MasterDetail<T> extends Control {
  /** The current item's key. */
  readonly current: Signal<string | null>;

  private readonly _options: MasterDetailOptions<T>;
  private readonly _items: ReadonlySignal<readonly T[]>;
  private readonly _list: HTMLElement;
  private readonly _rows = new Map<string, Row>();
  private readonly _panes = new Map<string, Pane>();
  private readonly _ids = new Map<string, string>();

  constructor(options: MasterDetailOptions<T>) {
    super();
    this._options = options;
    const items = options.items;
    this._items = items instanceof Signal ? items as ReadonlySignal<readonly T[]> : signal(items as readonly T[]);
    this.current = options.current ?? signal<string | null>(null);
    const root = this.root;
    root.classList.add('u2-master-detail');
    root.dataset.u2 = 'master-detail';
    if (options.name !== undefined)
      this.name = options.name;

    const list = this._list = div([], 'u2-master-detail-list');
    list.dataset.u2Part = 'list';
    list.setAttribute('role', 'listbox');
    if (options.label !== undefined)
      list.setAttribute('aria-label', options.label);
    list.tabIndex = 0;
    const empty = span(options.emptyText ?? 'Nothing to show', 'u2-master-detail-empty');
    const detail = div([empty], 'u2-master-detail-detail');
    detail.dataset.u2Part = 'detail';
    const split = options.split ?? 0.33;
    const splitter = this.build(() => {
      each(list, this._items, (item) => options.keyOf(item), (item) => this._row(item));
      return new Splitter([list, detail], {direction: 'horizontal', sizes: [split, 1 - split], minSize: 120});
    });
    root.append(splitter.root);

    this.effect(() => {
      const all = this._items.value;
      const key = this.current.value;
      untracked(() => this._show(all, key, detail, empty));
    });
    this.scope.listen(list, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
  }

  /** One item per row, keyed by its key: the current one selected, a checked one checked. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const current = this.current.peek();
    return Input.itemsState(this._items.peek().map((item, index) => {
      const key = this._options.keyOf(item);
      const row = this._rows.get(key);
      return {key, label: row === undefined ? key : Control.text(row.label), index, selected: key === current,
        ...row?.check === undefined ? {} : {checked: row.check.peek()}, rendered: row !== undefined};
    }), opts?.all);
  }

  /** The base state plus the current item's key. */
  getState(): ReturnType<Control['getState']> {
    return {...super.getState(), current: this.current.peek() ?? undefined};
  }

  private _show(items: readonly T[], key: string | null, detail: HTMLElement, empty: HTMLElement): void {
    const keys = items.map(this._options.keyOf);
    for (const [k, pane] of this._panes) {
      if (keys.includes(k))
        continue;
      pane.scope.dispose();
      pane.el.remove();
      this._panes.delete(k);
    }
    empty.hidden = keys.length > 0;
    if (key === null || !keys.includes(key)) {
      this.current.value = keys[0] ?? null;
      if (keys.length > 0)
        return;
    }
    for (const [k, pane] of this._panes)
      pane.el.hidden = k !== key;
    if (key !== null && keys.includes(key) && !this._panes.has(key))
      detail.append(this._pane(items[keys.indexOf(key)], key));
    const row = key === null ? undefined : this._rows.get(key);
    if (row === undefined)
      this._list.removeAttribute('aria-activedescendant');
    else
      this._list.setAttribute('aria-activedescendant', row.el.id);
  }

  private _row(item: T): HTMLElement {
    const key = this._options.keyOf(item);
    const spec = this._options.row(item);
    const label = span(spec.label, 'u2-master-detail-label');
    label.dataset.u2Part = 'item-label';
    const el = div([], 'u2-master-detail-row');
    el.setAttribute('role', 'option');
    el.id = this._idOf(key);
    label.id = `${el.id}-label`;
    const check = spec.check;
    if (check !== undefined) {
      // the row is the option and says it is checked; the box is what the pointer presses
      const box = document.createElement('input');
      box.type = 'checkbox';
      box.className = 'u2-input-checkbox u2-master-detail-check';
      box.tabIndex = -1;
      box.setAttribute('aria-hidden', 'true');
      on(box, 'change', () => check.value = box.checked);
      if (spec.checkTooltip !== undefined)
        Tooltip.bind(box, spec.checkTooltip);
      effect(() => {
        box.checked = check.value;
        el.setAttribute('aria-checked', String(check.value));
      });
      el.append(box);
    }
    el.append(label);
    if (spec.hint !== undefined)
      el.append(span(spec.hint, 'u2-master-detail-hint'));
    on(el, 'click', () => this.current.value = key);
    effect(() => {
      const current = this.current.value === key;
      el.classList.toggle('u2-master-detail-row-current', current);
      el.setAttribute('aria-selected', String(current));
    });
    const row = {el, label, check};
    this._rows.set(key, row);
    Ownership.scopeOf(undefined, 'a MasterDetail row').own(() => {
      if (this._rows.get(key) === row)
        this._rows.delete(key);
    });
    return el;
  }

  private _pane(item: T, key: string): HTMLElement {
    const scope = this.scope.child('master-detail-pane');
    const el = scope.build(() => {
      const built = this._options.detail(item);
      const content = div([Control.is(built) ? built.root : built], 'u2-master-detail-content');
      const reason = span('', 'u2-master-detail-reason');
      reason.dataset.u2Part = 'reason';
      const pane = div([reason, content], 'u2-master-detail-pane');
      pane.setAttribute('role', 'region');
      pane.setAttribute('aria-labelledby', `${this._idOf(key)}-label`);
      Reason.follow(this._options.disabled?.(item) ?? null, (why) => {
        reason.textContent = why ?? '';
        reason.hidden = why === null;
        Reason.region(content, why);
      }, scope);
      return pane;
    });
    this._panes.set(key, {el, scope});
    return el;
  }

  // a rebuilt row keeps its key's id, so its pane stays labelled by it
  private _idOf(key: string): string {
    let id = this._ids.get(key);
    if (id === undefined)
      this._ids.set(key, id = `${this.rootId}-item-${this._ids.size}`);
    return id;
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const keys = this._items.peek().map(this._options.keyOf);
    const at = keys.indexOf(this.current.peek() ?? '');
    const to = e.key === 'ArrowDown' ? Math.min(at + 1, keys.length - 1) : e.key === 'ArrowUp' ? Math.max(at - 1, 0) :
      e.key === 'Home' ? 0 : e.key === 'End' ? keys.length - 1 : -1;
    if (to >= 0 && keys.length > 0) {
      e.preventDefault();
      this.current.value = keys[to];
      return;
    }
    const check = at < 0 ? undefined : this._rows.get(keys[at])?.check;
    if (e.key !== ' ' || check === undefined)
      return;
    e.preventDefault();
    check.value = !check.peek();
  }
}
