import {Input, InputOptions} from '../core/input-base.js';

export type ChoiceItem = string | {value: string, label: string};

function itemValue(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.value;
}

function itemLabel(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.label;
}

export interface ChoiceInputOptions extends InputOptions<string | null> {
  items: ChoiceItem[];
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

  /** Replaces the item list, keeping the current value if it is still one of the items. */
  setItems(items: ChoiceItem[]): void {
    this._items = items;
    this._fill();
    const value = this.value.peek();
    this.value.value = value !== null && items.some((i) => itemValue(i) === value) ? value : null;
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
    if (this.nullable)
      select.append(new Option('', ''));
    for (const item of this._items)
      select.append(new Option(itemLabel(item), itemValue(item)));
  }
}

export interface MultiChoiceInputOptions extends InputOptions<string[]> {
  items: ChoiceItem[];
}

/** Compact checkbox list; the value holds the checked items in item order. */
export class MultiChoiceInput extends Input<string[], MultiChoiceInputOptions> {
  constructor(options: MultiChoiceInputOptions) {
    super(options, []);
    this.root.dataset.u2 = 'multi-choice-input';
  }

  protected createEditor(): HTMLElement {
    const list = document.createElement('div');
    list.className = 'u2-multi-choice';
    list.setAttribute('role', 'group');
    const boxes: HTMLInputElement[] = [];
    for (const item of this.options.items) {
      const row = document.createElement('label');
      row.className = 'u2-multi-choice-item';
      const box = document.createElement('input');
      box.type = 'checkbox';
      box.className = 'u2-input-checkbox';
      box.value = itemValue(item);
      const text = document.createElement('span');
      text.textContent = itemLabel(item);
      row.append(box, text);
      list.append(row);
      boxes.push(box);
    }
    const onChange = () => {
      this.value.value = boxes.filter((b) => b.checked).map((b) => b.value);
    };
    list.addEventListener('change', onChange);
    this.own(() => list.removeEventListener('change', onChange));
    this.effect(() => {
      const selected = new Set(this.value.value);
      for (const box of boxes)
        box.checked = selected.has(box.value);
    });
    return list;
  }
}
