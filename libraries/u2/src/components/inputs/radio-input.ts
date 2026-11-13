import {Input, InputOptions} from '../../core/input-base.js';
import type {IItemsState} from '../../core/widget-like.js';

export interface RadioInputOptions extends InputOptions<string | null> {
  items: string[];
  /** Raised-button variant, as Dart's `RadioInput.buttonStyle`. */
  buttons?: boolean;
  itemTooltips?: Record<string, string>;
}

/** Radio group over native `input[type=radio]` sharing one `name`, so arrow-key navigation,
 * label clicks and screen readers work without machinery of our own. Each row is an addressable
 * item (`data-u2-item`); the radio inside it carries the role and the checked state, and its
 * `data-value`, as the platform's radios do. */
export class RadioInput extends Input<string | null, RadioInputOptions> {
  private static _seq = 0;

  private declare _group: HTMLElement;
  private declare _name: string;
  private declare _items: string[];
  private declare _radios: HTMLInputElement[];

  constructor(options: RadioInputOptions) {
    super(options, null);
    this.root.dataset.u2 = 'radio-input';
  }

  get items(): string[] {
    return this._items;
  }

  set items(x: string[]) {
    this.setItems(x ?? []);
  }

  get itemTooltips(): Record<string, string> | undefined {
    return this.options.itemTooltips;
  }

  set itemTooltips(x: Record<string, string> | undefined) {
    this.options.itemTooltips = x ?? undefined;
    this._fill();
    this._apply(this.peek());
  }

  /** Replaces the item list, keeping the current value if it is still one of the items. */
  setItems(items: string[]): void {
    this._items = items;
    this._fill();
    const value = this.peek();
    this.value = value != null && items.includes(value) ? value : null;
    this._apply(this.peek());
  }

  /** One item per radio, keyed by its value, the chosen one checked. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const value = this.peek();
    return Input.itemsState(this._items.map((item, index) =>
      ({key: item, label: item, index, checked: item === value, rendered: true})), opts?.all);
  }

  /** The chosen item. */
  protected displayText(): string {
    return this.peek() ?? '';
  }

  protected createEditor(): HTMLElement {
    this._name = `u2-radio-${++RadioInput._seq}`;
    this._items = this.options.items;
    const group = document.createElement('div');
    this._group = group;
    group.className = 'u2-radio';
    group.setAttribute('role', 'radiogroup');
    if (this.options.buttons)
      group.classList.add('u2-radio-buttons');
    this._fill();

    this.scope.listen(group, 'change', () => {
      const checked = this._radios.find((r) => r.checked);
      this.value = checked ? checked.value : null;
    });
    this.effect(() => this._apply(this.value));
    return group;
  }

  private _fill(): void {
    const tooltips = this.options.itemTooltips;
    this._group.textContent = '';
    this._radios = [];
    for (const item of this._items) {
      const row = document.createElement('label');
      row.className = 'u2-radio-item';
      row.dataset.u2Item = item;
      const radio = document.createElement('input');
      radio.type = 'radio';
      radio.name = this._name;
      radio.value = item;
      radio.dataset.value = item;
      radio.className = 'u2-input-radio';
      const text = document.createElement('span');
      text.dataset.u2Part = 'item-label';
      text.textContent = item;
      row.append(radio, text);
      const tooltip = tooltips ? tooltips[item] : undefined;
      if (tooltip)
        row.title = tooltip;
      this._group.append(row);
      this._radios.push(radio);
    }
    this.refreshDisabled();
  }

  private _apply(value: string | null): void {
    for (const radio of this._radios) {
      radio.checked = radio.value === value;
      (radio.parentElement as HTMLElement).classList.toggle('u2-radio-pressed', radio.checked);
    }
  }
}
