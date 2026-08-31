import {Input, InputOptions} from '../../core/input-base.js';

export interface RadioInputOptions extends InputOptions<string | null> {
  items: string[];
  /** Raised-button variant, as Dart's `RadioInput.buttonStyle`. */
  buttons?: boolean;
  itemTooltips?: Record<string, string>;
}

/** Radio group over native `input[type=radio]` sharing one `name`, so arrow-key navigation,
 * label clicks and screen readers work without machinery of our own. */
export class RadioInput extends Input<string | null, RadioInputOptions> {
  private static _seq = 0;

  private _group!: HTMLElement;
  private _name!: string;
  private _items!: string[];
  private _radios!: HTMLInputElement[];

  constructor(options: RadioInputOptions) {
    super(options, null);
    this.root.dataset.u2 = 'radio-input';
  }

  /** Replaces the item list, keeping the current value if it is still one of the items. */
  setItems(items: string[]): void {
    this._items = items;
    this._fill();
    const value = this.value.peek();
    this.value.value = value !== null && items.includes(value) ? value : null;
    this._apply(this.value.peek());
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

    const onChange = () => {
      const checked = this._radios.find((r) => r.checked);
      this.value.value = checked ? checked.value : null;
    };
    group.addEventListener('change', onChange);
    this.own(() => group.removeEventListener('change', onChange));
    this.effect(() => this._apply(this.value.value));
    return group;
  }

  private _fill(): void {
    const tooltips = this.options.itemTooltips;
    this._group.textContent = '';
    this._radios = [];
    for (const item of this._items) {
      const row = document.createElement('label');
      row.className = 'u2-radio-item';
      const radio = document.createElement('input');
      radio.type = 'radio';
      radio.name = this._name;
      radio.value = item;
      radio.className = 'u2-input-radio';
      const text = document.createElement('span');
      text.textContent = item;
      row.append(radio, text);
      const tooltip = tooltips ? tooltips[item] : undefined;
      if (tooltip)
        row.title = tooltip;
      this._group.append(row);
      this._radios.push(radio);
    }
    this.refreshEnabled();
  }

  private _apply(value: string | null): void {
    for (const radio of this._radios) {
      radio.checked = radio.value === value;
      (radio.parentElement as HTMLElement).classList.toggle('u2-radio-pressed', radio.checked);
    }
  }
}
