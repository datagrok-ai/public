/* ChecklistInput: checkboxes the user can also put in order — the compute functions to run and the order to run
   them in, the columns to show. The value is the checked items in the order shown; `order` is every item's. */
import {Input, InputOptions} from '../../core/input-base.js';
import {signal, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import {div, span, effect, on} from '../../core/elements.js';
import {each} from '../../core/region.js';
import {Reorder} from './reorder.js';
import type {ChoiceItem} from './choice-input.js';
import type {IItemsState} from '../../core/widget-like.js';

export interface ChecklistInputOptions extends InputOptions<string[]> {
  /** The items, a value each; a signal is followed, and a checked item that leaves the list leaves the value. */
  items: ChoiceItem[] | ReadonlySignal<ChoiceItem[]>;
  /** A drag handle per row and Alt+↑/↓: the order is part of the pick. Checked items given as the value come
   * first, in the value's order. */
  orderable?: boolean;
  /** What stands in for an empty item list: 'No items' by default, a row high, as MultiChoiceInput's. */
  emptyText?: string;
}

export class ChecklistInput extends Input<string[], ChecklistInputOptions> {
  // built by createEditor(), which the base constructor calls before subclass initializers would run
  private declare _items: Signal<ChoiceItem[]>;
  private declare _order: Signal<string[]>;
  private declare _list: HTMLElement;

  constructor(options: ChecklistInputOptions) {
    super(options, []);
    this.root.dataset.u2 = 'checklist-input';
  }

  /** Every item's value in the order shown, checked or not — what an orderable list persists. */
  get order(): ReadonlySignal<string[]> {
    return this._order;
  }

  get items(): ChoiceItem[] {
    return this._items.peek();
  }

  set items(x: ChoiceItem[]) {
    this.setItems(x);
  }

  /** Replaces the items. They show in their own order, as MultiChoiceInput's; an orderable list keeps the user's
   * order for the ones still there and puts the new ones last. A checked item that left takes its check with it
   * (the input's own write, not the user's). */
  setItems(items: ChoiceItem[]): void {
    const values = items.map(ChecklistInput._value);
    const kept = this._order.peek().filter((v) => values.includes(v));
    this._items.value = items;
    this._order.value = this.options.orderable ? [...kept, ...values.filter((v) => !kept.includes(v))] : values;
    const value = this.peek() ?? [];
    const pruned = value.filter((v) => values.includes(v));
    if (pruned.length !== value.length)
      this.systemWrite(pruned);
  }

  /** Moves the item shown at `from` to `to`. */
  move(from: number, to: number): void {
    this._place(from, to);
    this._commit();
  }

  /** One item per row in the order shown; the checked ones are the selection. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const checked = new Set(this.peek() ?? []);
    return Input.itemsState(this._order.peek().map((value, index) => ({key: value, label: this._label(value), index,
      checked: checked.has(value), selected: checked.has(value), rendered: true})), opts?.all);
  }

  /** The checked captions, in order. */
  protected displayText(): string {
    return (this.peek() ?? []).map((v) => this._label(v)).join(', ');
  }

  protected createEditor(): HTMLElement {
    const items = this.options.items;
    const list = this._list = div([], 'u2-checklist');
    list.setAttribute('role', 'list');
    this._items = signal(items instanceof Signal ? items.peek() : items);
    const values = this._items.peek().map(ChecklistInput._value);
    const checked = (this.peek() ?? []).filter((v) => values.includes(v));
    this._order = signal(this.options.orderable ? [...checked, ...values.filter((v) => !checked.includes(v))] : values);
    each(list, this._order, (value) => value, (value) => this._row(value),
      {empty: () => span(this.options.emptyText ?? 'No items', 'u2-checklist-empty')});
    if (items instanceof Signal) {
      let first = true;
      this.effect(() => {
        const next = items.value;
        if (!first)
          untracked(() => this.setItems(next));
        first = false;
      });
    }
    // a value written from elsewhere puts its checked items first, in its order; our own commits already agree
    if (this.options.orderable) {
      this.effect(() => {
        const value = (this.value ?? []).filter((v) => this._order.peek().includes(v));
        const shown = this._order.peek().filter((v) => value.includes(v));
        if (shown.some((v, i) => v !== value[i]))
          this._order.value = [...value, ...this._order.peek().filter((v) => !value.includes(v))];
      });
    }
    return list;
  }

  private _row(value: string): HTMLElement {
    const row = div([], 'u2-checklist-item');
    row.setAttribute('role', 'listitem');
    if (this.options.orderable) {
      const handle = Reorder.handle();
      row.append(handle);
      Reorder.wire(row, handle, () => Array.from(this._list.children) as HTMLElement[],
        (from, to) => this._place(from, to), () => this._commit());
    }
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.className = 'u2-input-checkbox';
    box.value = value;
    const label = document.createElement('label');
    label.className = 'u2-checklist-label';
    const text = span(this._label(value));
    text.dataset.u2Part = 'item-label';
    label.append(box, text);
    row.append(label);
    effect(() => {
      box.checked = (this.value ?? []).includes(value);
    });
    on(box, 'change', () => this._commit());
    box.disabled = this.disabled !== null;
    return row;
  }

  private _place(from: number, to: number): void {
    const order = this._order.peek().slice();
    const [value] = order.splice(from, 1);
    order.splice(to, 0, value);
    this._order.value = order;
  }

  /** The value from the boxes, in the order shown; written only when it differs. */
  private _commit(): void {
    const boxes = Array.from(this._list.querySelectorAll<HTMLInputElement>('.u2-input-checkbox'));
    const checked = new Set(boxes.filter((b) => b.checked).map((b) => b.value));
    const next = this._order.peek().filter((v) => checked.has(v));
    const value = this.peek() ?? [];
    if (next.length !== value.length || next.some((v, i) => v !== value[i]))
      this.value = next;
  }

  private _label(value: string): string {
    const item = this._items.peek().find((i) => ChecklistInput._value(i) === value);
    return item === undefined ? value : typeof item === 'string' ? item : item.label;
  }

  private static _value(item: ChoiceItem): string {
    return typeof item === 'string' ? item : item.value;
  }
}
