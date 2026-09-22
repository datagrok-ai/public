/* Toggle chips: every item is a pressable chip, the value holds the pressed ones in item order —
   the inline choice for a handful of names where a checkbox list is too tall and a popup
   (`MultiSelect`) hides what is chosen. `tag()` (badge.ts) stays the removable, non-toggling chip.
   Not a `ButtonGroup({toggle: 'multi'})`: that is a segmented Control with fixed items, this is an
   `Input` (label, enabled, validation, `setItems`) whose chips wrap. */
import {Input, InputOptions} from '../../core/input-base.js';
import {ChoiceItem} from './choice-input.js';

export interface ChipsInputOptions extends InputOptions<string[]> {
  items: ChoiceItem[];
  /** Stands in for the chips while there is nothing to offer ('No items' by default). */
  emptyText?: string;
}

function itemValue(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.value;
}

function itemLabel(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.label;
}

export class ChipsInput extends Input<string[], ChipsInputOptions> {
  private _host!: HTMLElement;
  private _items!: ChoiceItem[];

  constructor(options: ChipsInputOptions) {
    super(options, []);
    this.root.dataset.u2 = 'chips-input';
  }

  get items(): ChoiceItem[] {
    return this._items;
  }

  /** Replaces the chips, dropping picks whose item vanished — a {@link Input.system} write. */
  setItems(items: ChoiceItem[]): void {
    this._items = items;
    const keep = new Set(items.map(itemValue));
    const value = this.value.peek();
    const pruned = value.filter((v) => keep.has(v));
    if (pruned.length !== value.length)
      Input.system(() => this.value.value = pruned);
    this._render();
  }

  protected createEditor(): HTMLElement {
    this._items = this.options.items;
    const host = document.createElement('div');
    this._host = host;
    host.className = 'u2-chips';
    host.setAttribute('role', 'group');
    const onClick = (e: Event) => {
      const chip = (e.target as Element).closest('.u2-chip') as HTMLElement | null;
      if (chip !== null && this.enabled)
        this._toggle(chip.dataset.value!);
    };
    host.addEventListener('click', onClick);
    this.own(() => host.removeEventListener('click', onClick));
    this.effect(() => {
      this.value.value;
      this._render();
    });
    return host;
  }

  private _toggle(value: string): void {
    const picked = new Set(this.value.peek());
    if (picked.has(value))
      picked.delete(value);
    else
      picked.add(value);
    this.value.value = this._items.map(itemValue).filter((v) => picked.has(v));
  }

  private _render(): void {
    const picked = new Set(this.value.peek());
    this._host.textContent = '';
    if (this._items.length === 0) {
      const empty = document.createElement('span');
      empty.className = 'u2-chips-empty';
      empty.textContent = this.options.emptyText ?? 'No items';
      this._host.append(empty);
      return;
    }
    for (const item of this._items) {
      const chip = document.createElement('button');
      chip.type = 'button';
      chip.className = 'u2-chip';
      chip.dataset.value = itemValue(item);
      chip.textContent = itemLabel(item);
      const on = picked.has(itemValue(item));
      chip.classList.toggle('u2-chip-on', on);
      chip.setAttribute('aria-pressed', String(on));
      this._host.append(chip);
    }
    this.refreshEnabled();
  }
}
