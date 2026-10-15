/* Toggle chips: every item is a pressable chip, the value holds the pressed ones in item order —
   the inline choice for a handful of names where a checkbox list is too tall and a popup
   (`MultiSelect`) hides what is chosen. `tag()` (badge.ts) stays the removable, non-toggling chip.
   Not a `ButtonGroup({toggle: 'multi'})`: that is a segmented Control with fixed items, this is an
   `Input` (label, enabled, validation, `setItems`) whose chips wrap. */
import {signal, Signal} from '../../core/signals.js';
import {Input, InputOptions} from '../../core/input-base.js';
import {each} from '../../core/region.js';
import {ChoiceItem} from './choice-input.js';
import type {IItemsState} from '../../core/widget-like.js';

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
  private declare _items: Signal<ChoiceItem[]>;

  constructor(options: ChipsInputOptions) {
    super(options, []);
    this.root.dataset.u2 = 'chips-input';
  }

  get items(): ChoiceItem[] {
    return this._items.peek();
  }

  /** What stands in for an empty item list. */
  get emptyText(): string {
    return this.options.emptyText ?? 'No items';
  }

  set emptyText(x: string) {
    this.options.emptyText = x ?? undefined;
    const empty = this.root.querySelector('.u2-chips-empty');
    if (empty !== null)
      empty.textContent = this.emptyText;
  }

  set items(x: ChoiceItem[]) {
    this.setItems(x);
  }

  /** Replaces the chips, dropping picks whose item vanished — a {@link Input.system} write. Chips whose
   * item stays keep their element (and focus). */
  setItems(items: ChoiceItem[]): void {
    this._items.value = items;
    const keep = new Set(items.map(itemValue));
    const value = this.value.peek() ?? [];
    const pruned = value.filter((v) => keep.has(v));
    if (pruned.length !== value.length)
      Input.system(() => this.value.value = pruned);
  }

  /** One item per chip, keyed by its value; the pressed ones are the selection. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const picked = new Set(this.value.peek());
    return Input.itemsState(this._items.peek().map((item, index) => ({key: itemValue(item), label: itemLabel(item),
      index, selected: picked.has(itemValue(item)), rendered: true})), opts?.all);
  }

  /** The pressed captions, in item order. */
  protected displayText(): string {
    const picked = new Set(this.value.peek());
    return this._items.peek().filter((i) => picked.has(itemValue(i))).map(itemLabel).join(', ');
  }

  protected createEditor(): HTMLElement {
    this._items = signal(this.options.items);
    const host = document.createElement('div');
    host.className = 'u2-chips';
    host.setAttribute('role', 'group');
    this.scope.listen(host, 'click', (e) => {
      const chip = (e.target as Element).closest('.u2-chips-chip') as HTMLElement | null;
      if (chip !== null && this.enabled)
        this._toggle(chip.dataset.value!);
    });
    each(host, this._items, (item) => itemValue(item), (item) => this._chip(item), {
      update: (chip, item) => {
        chip.querySelector('[data-u2-part="item-label"]')!.textContent = itemLabel(item);
      },
      empty: () => {
        const empty = document.createElement('span');
        empty.className = 'u2-chips-empty';
        empty.textContent = this.options.emptyText ?? 'No items';
        return empty;
      },
    });
    this.effect(() => {
      const picked = new Set(this.value.value);
      for (const chip of host.querySelectorAll<HTMLElement>('.u2-chips-chip'))
        ChipsInput._press(chip, picked.has(chip.dataset.value!));
    });
    return host;
  }

  private _chip(item: ChoiceItem): HTMLElement {
    const chip = document.createElement('button');
    chip.type = 'button';
    chip.className = 'u2-chips-chip';
    chip.setAttribute('role', 'button');
    chip.dataset.value = itemValue(item);
    const text = document.createElement('span');
    text.dataset.u2Part = 'item-label';
    text.textContent = itemLabel(item);
    chip.append(text);
    chip.disabled = !this.enabled;
    ChipsInput._press(chip, new Set(this.value.peek()).has(itemValue(item)));
    return chip;
  }

  private _toggle(value: string): void {
    const picked = new Set(this.value.peek());
    if (picked.has(value))
      picked.delete(value);
    else
      picked.add(value);
    this.value.value = this._items.peek().map(itemValue).filter((v) => picked.has(v));
  }

  private static _press(chip: HTMLElement, on: boolean): void {
    chip.classList.toggle('u2-chips-chip-on', on);
    chip.setAttribute('aria-pressed', String(on));
  }
}
