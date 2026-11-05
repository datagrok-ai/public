import {signal, untracked, Signal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {Usage} from '../../core/scope.js';
import {Input, InputOptions} from '../../core/input-base.js';
import {iconOf} from '../display/icon.js';
import type {IItemsState} from '../../core/widget-like.js';

/** `group` puts the item under that heading; items of one group are listed together. */
export type ChoiceItem = string | {value: string, label: string, group?: string};

function itemValue(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.value;
}

function itemLabel(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.label;
}

/** The items less any null one — `grok.shell.tableNames` lists a nameless table as `undefined` (SEAM S-23) —
 * which is reported and skipped, so the rest still render. */
function present(items: ChoiceItem[]): ChoiceItem[] {
  if (items.every((item) => item != null))
    return items;
  Usage.report('item-null', 'a choice list holds a null item, which is skipped: pass strings or {value, label} ' +
    'items (for open tables, `Platform.tableNames()`)');
  return items.filter((item) => item != null);
}

export interface ChoiceInputOptions extends InputOptions<string | null> {
  /** The choices — a list, or a signal followed live (each change is {@link ChoiceInput.setItems}). The
   * value may arrive before or after the item it names: it stays selected once both are there. */
  items: ChoiceItem[] | ReadonlySignal<ChoiceItem[]>;
  /** Shown as a disabled placeholder option while the item list is empty, so a picker over a
   * live collection reads as empty rather than broken. */
  emptyText?: string;
}

/** Single choice over a styled native `<select>`, as `ui.input.choice` does — the platform look
 * comes free and the picker is the OS one. Combobox covers the custom-popup case. */
export class ChoiceInput extends Input<string | null, ChoiceInputOptions> {
  private declare _select: HTMLSelectElement;
  private declare _items: ChoiceItem[];

  constructor(options: ChoiceInputOptions) {
    super(options, null);
    this.root.dataset.u2 = 'choice-input';
    const items = options.items;
    if (items instanceof Signal) {
      let first = true;
      this.effect(() => {
        const list = items.value;
        if (!first)
          untracked(() => this.setItems(list));
        first = false;
      });
    }
  }

  get items(): ChoiceItem[] {
    return this._items;
  }

  set items(x: ChoiceItem[]) {
    this.setItems(x);
  }

  /** Replaces the item list, keeping the current value if it is still one of the items. Dropping
   * a value whose item vanished is a {@link Input.system} write, not a user edit. */
  setItems(items: ChoiceItem[]): void {
    this._items = present(items);
    this._fill();
    const value = this.peek();
    this.systemWrite(value !== null && this._items.some((i) => itemValue(i) === value) ? value : null);
    this._select.value = this.peek() ?? '';
  }

  /** One item per choice (the empty option is none), keyed by its value, the chosen one selected. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const value = this.peek();
    return Input.itemsState(this._items.map((item, index) => ({key: itemValue(item), label: itemLabel(item),
      index, selected: itemValue(item) === value, rendered: true})), opts?.all);
  }

  protected createEditor(): HTMLElement {
    const select = document.createElement('select');
    this._select = select;
    const items = this.options.items;
    this._items = present(items instanceof Signal ? items.peek() : items);
    this._fill();
    this.scope.listen(select, 'change', () => this.value = select.value === '' ? null : select.value);
    this.effect(() => select.value = this.value ?? '');
    return select;
  }

  private _fill(): void {
    const select = this._select;
    select.textContent = '';
    if (this._items.length === 0 && this.options.emptyText !== undefined) {
      const placeholder = new Option(this.options.emptyText, '');
      placeholder.disabled = true;
      placeholder.selected = true;
      select.append(placeholder);
      return;
    }
    if (this.nullable)
      select.append(new Option('', ''));
    let parent: HTMLElement = select;
    let heading: string | undefined;
    for (const item of this._items) {
      const name = typeof item === 'string' ? undefined : item.group;
      if (name !== heading) {
        heading = name;
        parent = name === undefined ? select : Object.assign(document.createElement('optgroup'), {label: name});
        if (parent !== select)
          select.append(parent);
      }
      parent.append(new Option(itemLabel(item), itemValue(item)));
    }
  }
}

export interface MultiChoiceInputOptions extends InputOptions<string[]> {
  items: ChoiceItem[];
  /** Stands in for the checkbox list while there is nothing to check, so a picker over a live
   * collection reads as empty rather than broken. 'No items' by default. */
  emptyText?: string;
  /** Summary row atop the list: a tri-state checkbox plus "<N> of <M>". Clicking the checkbox
   * selects all / none; clicking the text (or its chevron) collapses/expands the list. */
  showSummaryCheckbox?: boolean;
}

let summaryCount = 0;

/** Compact checkbox list; the value holds the checked items in item order. */
export class MultiChoiceInput extends Input<string[], MultiChoiceInputOptions> {
  // built by createEditor(), which the base constructor calls before subclass initializers
  // would run — `declare` keeps them out of the emitted class fields
  private declare _itemsHost: HTMLElement;
  private declare _boxes: HTMLInputElement[];
  private declare _items: ChoiceItem[];
  private declare _summaryBox: HTMLInputElement | undefined;
  private declare _summaryText: HTMLElement | undefined;
  private declare _summaryToggle: HTMLElement | undefined;
  private declare _expanded: Signal<boolean> | undefined;

  /** What stands in for an empty item list. */
  get emptyText(): string {
    return this.options.emptyText ?? 'No items';
  }

  set emptyText(x: string) {
    this.options.emptyText = x ?? undefined;
    const empty = this._itemsHost.querySelector('.u2-multi-choice-empty');
    if (empty !== null)
      empty.textContent = this.emptyText;
  }

  get disabled(): string | null {
    return super.disabled;
  }

  // the base sweep disables the native controls; the collapse toggle is a span and must give up
  // its tab stop and gestures itself
  set disabled(reason: string | null) {
    super.disabled = reason;
    if (this._summaryToggle !== undefined)
      this._summaryToggle.tabIndex = this.disabled === null ? 0 : -1;
  }

  constructor(options: MultiChoiceInputOptions) {
    super(options, []);
    this.root.dataset.u2 = 'multi-choice-input';
  }

  get items(): ChoiceItem[] {
    return this._items;
  }

  set items(x: ChoiceItem[]) {
    this.setItems(x);
  }

  /** One item per choice, keyed by its value; the checked ones are the selection. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const checked = new Set(this.peek() ?? []);
    return Input.itemsState(this._items.map((item, index) => {
      const on = checked.has(itemValue(item));
      return {key: itemValue(item), label: itemLabel(item), index, checked: on, selected: on, rendered: true};
    }), opts?.all);
  }

  /** The checked captions, in item order. */
  protected displayText(): string {
    const checked = new Set(this.peek() ?? []);
    return this._items.filter((i) => checked.has(itemValue(i))).map(itemLabel).join(', ');
  }

  /** Collapse state of the item list; present only with {@link MultiChoiceInputOptions.showSummaryCheckbox}. */
  get expanded(): Signal<boolean> | undefined {
    return this._expanded;
  }

  /** Replaces the item list, keeping the items that are still there checked
   * (`multi_choice_input.dart:19-37`). The value is rewritten only when something actually
   * vanished, and that rewrite is a {@link Input.system} write rather than a user edit — Dart
   * instead RETAINS a checked item the new list no longer offers, which leaves a value the user
   * can neither see nor uncheck. */
  setItems(items: ChoiceItem[]): void {
    this._fill(items);
    const current = new Set(this.peek() ?? []);
    const pruned = this._items.filter((i) => current.has(itemValue(i))).map(itemValue);
    if (pruned.length !== current.size)
      this.systemWrite(pruned);
  }

  protected createEditor(): HTMLElement {
    const list = document.createElement('div');
    list.className = 'u2-multi-choice';
    list.setAttribute('role', 'group');
    this._itemsHost = list;
    if (this.options.showSummaryCheckbox)
      this._buildSummary(list);
    this._fill(this.options.items);
    this.scope.listen(list, 'change', (e) => {
      if (this._summaryBox && e.target === this._summaryBox) {
        if (this._items.length === 0)
          this._updateSummary();
        else
          this.value = this._summaryBox.checked ? this._items.map(itemValue) : [];
        return;
      }
      this.value = this._boxes.filter((b) => b.checked).map((b) => b.value);
    });
    this.effect(() => {
      const selected = new Set(this.value ?? []);
      for (const box of this._boxes)
        box.checked = selected.has(box.value);
      this._updateSummary();
    });
    return list;
  }

  private _buildSummary(list: HTMLElement): void {
    const row = document.createElement('div');
    row.className = 'u2-multi-choice-summary';
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.className = 'u2-input-checkbox';
    this._summaryBox = box;
    const text = document.createElement('span');
    text.className = 'u2-multi-choice-summary-text';
    text.id = `u2-multi-choice-summary-${++summaryCount}`;
    box.setAttribute('aria-labelledby', text.id);
    this._summaryText = text;
    const chevron = document.createElement('span');
    chevron.className = 'u2-multi-choice-chevron';
    chevron.append(iconOf('chevron-down'));
    const toggle = document.createElement('span');
    toggle.className = 'u2-multi-choice-summary-toggle';
    toggle.tabIndex = 0;
    toggle.setAttribute('role', 'button');
    // its own name — the content-derived "N of M" would collide with the checkbox's label
    toggle.setAttribute('aria-label', 'Collapse list');
    toggle.append(text, chevron);
    row.append(box, toggle);
    const items = document.createElement('div');
    items.className = 'u2-multi-choice-items';
    list.append(row, items);
    this._itemsHost = items;
    this._summaryToggle = toggle;
    const expanded = signal(true);
    this._expanded = expanded;
    const flip = () => {
      if (this.disabled === null)
        expanded.value = !expanded.peek();
    };
    const onKeyDown = (e: Event) => {
      const key = (e as KeyboardEvent).key;
      if (key === 'Enter' || key === ' ') {
        e.preventDefault();
        flip();
      }
    };
    this.scope.listen(toggle, 'click', flip);
    this.scope.listen(toggle, 'keydown', onKeyDown);
    this.effect(() => {
      const on = expanded.value;
      toggle.setAttribute('aria-expanded', String(on));
      items.style.display = on ? '' : 'none';
    });
  }

  private _updateSummary(): void {
    const box = this._summaryBox;
    if (!box)
      return;
    const total = this._items.length;
    const checked = this._boxes.filter((b) => b.checked).length;
    box.checked = checked === total && total > 0;
    box.indeterminate = checked > 0 && checked < total;
    this._summaryText!.textContent = `${checked} of ${total}`;
  }

  private _fill(list: ChoiceItem[]): void {
    const items = present(list);
    this._items = items;
    const selected = new Set(this.peek() ?? []);
    this._itemsHost.textContent = '';
    this._boxes = [];
    if (items.length === 0) {
      const empty = document.createElement('div');
      empty.className = 'u2-multi-choice-empty';
      empty.textContent = this.options.emptyText ?? 'No items';
      this._itemsHost.append(empty);
      this._updateSummary();
      return;
    }
    for (const item of items) {
      const row = document.createElement('label');
      row.className = 'u2-multi-choice-item';
      row.dataset.u2Item = itemValue(item);
      row.dataset.value = itemValue(item);
      const box = document.createElement('input');
      box.type = 'checkbox';
      box.className = 'u2-input-checkbox';
      box.value = itemValue(item);
      box.checked = selected.has(box.value);
      const text = document.createElement('span');
      text.dataset.u2Part = 'item-label';
      text.textContent = itemLabel(item);
      row.append(box, text);
      this._itemsHost.append(row);
      this._boxes.push(box);
    }
    this._updateSummary();
    this.refreshDisabled();
  }
}
