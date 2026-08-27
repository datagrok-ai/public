import {Input, InputOptions} from '../../core/input-base.js';

export type ChoiceItem = string | {value: string, label: string};

function itemValue(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.value;
}

function itemLabel(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.label;
}

export interface ChoiceInputOptions extends InputOptions<string | null> {
  items: ChoiceItem[];
  /** Shown as a disabled placeholder option while the item list is empty, so a picker over a
   * live collection reads as empty rather than broken. */
  emptyText?: string;
}

/** Single choice over a styled native `<select>`, as `ui.input.choice` does — the platform look
 * comes free and the picker is the OS one. Combobox covers the custom-popup case. */
export class ChoiceInput extends Input<string | null, ChoiceInputOptions> {
  private _select!: HTMLSelectElement;
  private _items!: ChoiceItem[];
  private _echo = false;

  constructor(options: ChoiceInputOptions) {
    super(options, null);
    this.root.dataset.u2 = 'choice-input';
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
    this._items = items;
    this._fill();
    const value = this.value.peek();
    Input.system(() => this.value.value =
      value !== null && items.some((i) => itemValue(i) === value) ? value : null);
    this._select.value = this.value.peek() ?? '';
  }

  protected createEditor(): HTMLElement {
    const select = document.createElement('select');
    this._select = select;
    this._items = this.options.items;
    this._fill();
    const onChange = () => {
      this._echo = true;
      try {
        this.value.value = select.value === '' ? null : select.value;
      } finally {
        this._echo = false;
      }
    };
    select.addEventListener('change', onChange);
    this.own(() => select.removeEventListener('change', onChange));
    this.effect(() => {
      const value = this.value.value;
      if (!this._echo)
        select.value = value ?? '';
    });
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
    for (const item of this._items)
      select.append(new Option(itemLabel(item), itemValue(item)));
  }
}

export interface MultiChoiceInputOptions extends InputOptions<string[]> {
  items: ChoiceItem[];
  /** Stands in for the checkbox list while there is nothing to check, so a picker over a live
   * collection reads as empty rather than broken. 'No items' by default. */
  emptyText?: string;
}

/** Compact checkbox list; the value holds the checked items in item order. */
export class MultiChoiceInput extends Input<string[], MultiChoiceInputOptions> {
  private _list!: HTMLElement;
  private _boxes!: HTMLInputElement[];
  private _items!: ChoiceItem[];

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

  /** Replaces the item list, keeping the items that are still there checked
   * (`multi_choice_input.dart:19-37`). The value is rewritten only when something actually
   * vanished, and that rewrite is a {@link Input.system} write rather than a user edit — Dart
   * instead RETAINS a checked item the new list no longer offers, which leaves a value the user
   * can neither see nor uncheck. */
  setItems(items: ChoiceItem[]): void {
    this._fill(items);
    const current = new Set(this.value.peek());
    const pruned = items.filter((i) => current.has(itemValue(i))).map(itemValue);
    if (pruned.length !== current.size)
      Input.system(() => this.value.value = pruned);
  }

  protected createEditor(): HTMLElement {
    const list = document.createElement('div');
    this._list = list;
    list.className = 'u2-multi-choice';
    list.setAttribute('role', 'group');
    this._fill(this.options.items);
    const onChange = () => {
      this.value.value = this._boxes.filter((b) => b.checked).map((b) => b.value);
    };
    list.addEventListener('change', onChange);
    this.own(() => list.removeEventListener('change', onChange));
    this.effect(() => {
      const selected = new Set(this.value.value);
      for (const box of this._boxes)
        box.checked = selected.has(box.value);
    });
    return list;
  }

  private _fill(items: ChoiceItem[]): void {
    this._items = items;
    const selected = new Set(this.value.peek());
    this._list.textContent = '';
    this._boxes = [];
    if (items.length === 0) {
      const empty = document.createElement('div');
      empty.className = 'u2-multi-choice-empty';
      empty.textContent = this.options.emptyText ?? 'No items';
      this._list.append(empty);
      return;
    }
    for (const item of items) {
      const row = document.createElement('label');
      row.className = 'u2-multi-choice-item';
      const box = document.createElement('input');
      box.type = 'checkbox';
      box.className = 'u2-input-checkbox';
      box.value = itemValue(item);
      box.checked = selected.has(box.value);
      const text = document.createElement('span');
      text.textContent = itemLabel(item);
      row.append(box, text);
      this._list.append(row);
      this._boxes.push(box);
    }
    this.refreshEnabled();
  }
}
