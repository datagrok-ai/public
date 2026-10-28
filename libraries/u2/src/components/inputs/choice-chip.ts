import {signal, computed, Signal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Input, LiveOption} from '../../core/input-base.js';
import {Reason} from '../../core/reason.js';
import {Tooltip} from '../../core/tooltip.js';
import type {IItemsState} from '../../core/widget-like.js';
import {icon} from '../display/icon.js';
import {Menu} from '../navigation/menu.js';
import type {Action} from '../actions/actions.js';

/** A choice of a {@link ChoiceChip}: a value, or a value with the label its menu item shows and what it stands
 * for (the item's grey detail and the chip's tooltip). */
export type ChoiceChipItem = string | {value: string, label: string, description?: string};

export interface ChoiceChipOptions {
  /** The small caption before the value ("R1"); also what the chip is called. */
  label?: string;
  /** The automation id (`data-u2-name`). */
  name?: string;
  value?: string | null;
  /** Two-way binding; without it the chip owns its signal. */
  bind?: Signal<string | null>;
  /** The choices — a list, or a signal of one; the menu reads them each time it opens. */
  items: ChoiceChipItem[] | ReadonlySignal<ChoiceChipItem[]>;
  /** Runs after the user picked an item in the menu, never for a write from code. */
  onInput?: (value: string | null, chip: ChoiceChip) => void;
  /** Why the chip is off: greyed, its menu closed, the reason as its tooltip. `null` while on. */
  disabled?: Reason;
  /** Why the current choice cannot stand: the chip turns the error color, from the start, and says why as its
   * tooltip. `null` while it is fine. */
  invalid?: Reason;
  /** What the chip shows after its label: the current item's label by default, '–' while the value is none of
   * the items. */
  valueText?: LiveOption<string>;
  /** What the current choice stands for, said by the tooltip and to a screen reader: the current item's
   * description by default. */
  description?: LiveOption<string>;
  /** Entries after the choices, under a separator — what else the chip offers ("Edit definition…"); a function
   * is asked on every open. */
  actions?: Action[] | (() => Action[]);
  /** Shows the chip pressed while something it opened elsewhere is showing (its fields under the chips). */
  active?: LiveOption<boolean>;
}

const isSignal = <T>(x: LiveOption<T> | undefined): x is ReadonlySignal<T> => x instanceof Signal;

/** A compact choice inline among others, reading "label value ▾": a press opens a menu of the choices, each with
 * what it stands for, the current one checked. Under the pointer it opens down when it fits, else up. */
export class ChoiceChip extends Control {
  /** The value's signal, for wiring. */
  readonly signal: Signal<string | null>;
  /** Why the current choice cannot stand, or `null`. */
  readonly validity: ReadonlySignal<string | null>;

  private readonly _options: ChoiceChipOptions;
  private readonly _items: ReadonlySignal<ChoiceChipItem[]>;
  private readonly _text: ReadonlySignal<string>;
  private readonly _onInput: ((value: string | null, chip: ChoiceChip) => void) | undefined;
  private _menu: Menu | undefined;

  constructor(options: ChoiceChipOptions) {
    super(Reason.focusable(document.createElement('button')));
    this._options = options;
    const root = this.root as HTMLButtonElement;
    root.type = 'button';
    root.classList.add('u2-choice-chip');
    root.dataset.u2 = 'choice-chip';
    root.setAttribute('aria-haspopup', 'menu');
    root.setAttribute('aria-expanded', 'false');
    if (options.name !== undefined)
      this.name = options.name;
    this.signal = options.bind ?? signal(options.value ?? null);
    const items = options.items;
    this._items = items instanceof Signal ? items : signal(items);
    this.validity = computed(() => Reason.read(options.invalid));
    const current = computed(() => {
      const value = this.signal.value;
      return this._items.value.map(ChoiceChip._item).find((item) => item.value === value) ?? null;
    });
    this._text = computed(() => (isSignal(options.valueText) ? options.valueText.value : options.valueText) ??
      current.value?.label ?? '–');
    const description = computed(() =>
      (isSignal(options.description) ? options.description.value : options.description) ??
      current.value?.description ?? null);
    this._onInput = options.onInput === undefined ? undefined :
      Input.bindOwned(options.onInput, {owner: this.scope}) as (value: string | null, chip: ChoiceChip) => void;

    const value = document.createElement('span');
    value.className = 'u2-choice-chip-value';
    value.dataset.u2Part = 'value';
    if (options.label) {
      const label = document.createElement('span');
      label.className = 'u2-choice-chip-label';
      label.dataset.u2Part = 'label';
      label.textContent = options.label;
      // the space reads "R1 OH" to a screen reader; the flex row does not render it
      root.append(label, ' ');
    }
    const caret = icon('chevron-down', {cls: 'u2-choice-chip-caret'});
    caret.setAttribute('aria-hidden', 'true');
    root.append(value, caret);

    const said = computed(() => this.validity.value ?? description.value);
    this.effect(() => {
      value.textContent = this._text.value;
      const off = Reason.read(options.disabled);
      const problem = this.validity.value;
      Reason.apply(root, off);
      if (off === null) {
        const text = said.value;
        if (text === null)
          root.removeAttribute('aria-description');
        else
          root.setAttribute('aria-description', text);
      }
      root.classList.toggle('u2-choice-chip-invalid', problem !== null);
      root.setAttribute('aria-invalid', String(problem !== null));
      if (off !== null)
        this._menu?.close();
    });
    const active = options.active;
    if (active !== undefined)
      this.effect(() => root.classList.toggle('u2-choice-chip-active', isSignal(active) ? active.value : active));
    Tooltip.bind(root, () => {
      const off = Reason.now(options.disabled);
      if (off !== null)
        return off;
      const text = said.peek() ?? this._text.peek();
      return options.label ? `${options.label}: ${text}` : text;
    }, {owner: this.scope});

    this.scope.listen(root, 'click', () => this._toggle());
    this.scope.listen(root, 'keydown', (e) => {
      if ((e as KeyboardEvent).key !== 'ArrowDown' || this.isOpen)
        return;
      e.preventDefault();
      this._toggle();
    });
    this.own(() => this._menu?.close());
  }

  /** The value: read inside an `effect` or a `computed`, it is followed; set it to change it. */
  get value(): string | null {
    return this.signal.value;
  }

  set value(x: string | null) {
    this.signal.value = x;
  }

  /** The value, not followed. */
  peek(): string | null {
    return this.signal.peek();
  }

  /** Keeps the value on `read` (the model the chip shows), re-read whenever a signal it reads changes. */
  follow(read: () => string | null): this {
    return this.effect(() => {
      const value = read();
      if (!Object.is(value, this.signal.peek()))
        this.signal.value = value;
    });
  }

  /** Why the chip is off — `null` while it is on. */
  get disabled(): string | null {
    return Reason.now(this._options.disabled);
  }

  /** Whether its menu is open. */
  get isOpen(): boolean {
    return this._menu?.isOpen.peek() === true;
  }

  /** One item per choice, keyed by its value; the current one selected and checked. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const value = this.signal.peek();
    return Input.itemsState(this._items.peek().map(ChoiceChip._item).map((item, index) => ({key: item.value,
      label: item.label, index, selected: item.value === value, checked: item.value === value, rendered: true})),
    opts?.all);
  }

  /** The base state plus the caption, the value and the text shown for it, the verdict and whether the menu is
   * open. */
  getState(): ReturnType<Control['getState']> {
    const error = this.validity.peek();
    return {...super.getState(), value: this.signal.peek(), text: this._text.peek(), valid: error === null, error,
      open: this.isOpen};
  }

  protected accessibleLabel(): string | undefined {
    return this._options.label || super.accessibleLabel();
  }

  private _toggle(): void {
    if (this.isOpen) {
      this._menu!.close();
      return;
    }
    if (this.disabled !== null)
      return;
    this._menu = this.build(() => this._buildMenu());
    this._menu.show({anchor: this.root});
  }

  private _buildMenu(): Menu {
    const menu = Menu.popup();
    const value = this.signal.peek();
    for (const item of this._items.peek().map(ChoiceChip._item)) {
      menu.item(item.label, () => this._pick(item.value),
        {check: item.value === value, ...item.description === undefined ? {} : {detail: item.description}});
    }
    const actions = this._options.actions;
    const extra = typeof actions === 'function' ? actions() : actions ?? [];
    if (extra.length > 0)
      menu.separator();
    for (const action of extra) {
      menu.item(action.text, action.run, {disabled: action.disabled, shortcut: action.shortcut, detail: action.detail,
        check: action.check});
    }
    return menu;
  }

  private _pick(value: string): void {
    this.signal.value = value;
    this._onInput?.(value, this);
  }

  private static _item(item: ChoiceChipItem): {value: string, label: string, description?: string} {
    return typeof item === 'string' ? {value: item, label: item} : item;
  }
}
