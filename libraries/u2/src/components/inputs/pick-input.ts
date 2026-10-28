/* A reference picker: the value is the id a reference holds, the box shows the item's caption,
   over any `TypeAhead` (a domain table's rows, the platform's users or groups). */
import {Input, InputOptions} from '../../core/input-base.js';
import {div, span} from '../../core/elements.js';
import {TypeAhead} from './typeahead.js';
import {iconButton} from '../actions/buttons.js';

export interface PickInputOptions<T> extends InputOptions<string | null> {
  /** Builds the type-ahead the input is a face of. */
  typeAhead: () => TypeAhead<T>;
  idOf: (item: T) => string;
  /** The item behind a value written from outside — the id a form loaded; null when unknown. */
  resolve: (id: string) => Promise<T | null>;
  /** Follows the reference: an Open icon beside the clear while a value is held. */
  open?: (id: string) => void;
}

export class PickInput<T> extends Input<string | null, PickInputOptions<T>> {
  private declare _typeAhead: TypeAhead<T>;
  private declare _hint: HTMLElement;

  constructor(options: PickInputOptions<T>) {
    super(options, null);
    this.root.dataset.u2 = 'pick-input';
    this.root.append(this._hint);
  }

  get typeAhead(): TypeAhead<T> {
    return this._typeAhead;
  }

  /** No typing: the box reads as a dropdown and a click opens the list — for a short lookup. */
  asDropdown(): void {
    const input = this._typeAhead.root.querySelector('input')!;
    input.readOnly = true;
    this.root.classList.add('u2-pick-dropdown');
    this.scope.listen(input, 'click', () => {
      if (!this._typeAhead.isOpen.peek())
        this._typeAhead.open();
    });
  }

  // runs from the base constructor, so the coupling state lives here rather than in fields the
  // subclass initializers would reset afterwards
  protected createEditor(): HTMLElement {
    const {typeAhead, idOf, resolve} = this.options;
    const picker = typeAhead();
    this._typeAhead = picker;
    // the value is held through typing and clearing the text: only ✕ (or a write from outside) drops it
    let held: T | null = null;
    let gen = 0;
    const lookUp = async (id: string) => {
      const mine = ++gen;
      try {
        const item = await resolve(id);
        if (mine === gen && this.peek() === id && item !== null)
          picker.selected.value = item;
      } catch {
        // an unknown id keeps the value and shows nothing
      }
    };
    // a form seeds a string field with '' for an empty cell: empty either way, and never rewritten
    // to null on the input's own account — that write would dirty the row on render
    const empty = (id: string | null) => id === null || id === '';
    this.effect(() => {
      const id = this.value;
      const current = picker.selected.peek();
      if ((current === null ? null : idOf(current)) === id)
        return;
      held = null;
      if (empty(id))
        picker.selected.value = null;
      else
        void lookUp(id!);
    });
    this.effect(() => {
      const item = picker.selected.value;
      if (item === null)
        return;
      held = item;
      const id = idOf(item);
      if (id !== this.peek())
        this.value = id;
    });
    const input = picker.root.querySelector('input')!;
    const hint = span('', 'u2-pick-hint');
    hint.dataset.u2Part = 'hint';
    hint.hidden = true;
    this._hint = hint;
    // text that is not a pick never survives a blur: the held pick comes back (with a hint if
    // something was typed), else the box is cleared
    const onBlur = () => {
      if (picker.isPickPending || picker.selected.peek() !== null)
        return;
      const typed = picker.text.peek();
      if (held !== null) {
        picker.selected.value = held;
        hint.textContent = typed === '' ? '' : `No match for "${typed}" — kept ${picker.text.peek()}`;
        hint.hidden = typed === '';
      } else if (typed !== '')
        picker.resync();
    };
    const onFocus = () => hint.hidden = true;
    input.addEventListener('blur', onBlur);
    input.addEventListener('focus', onFocus);
    this.own(() => {
      input.removeEventListener('blur', onBlur);
      input.removeEventListener('focus', onFocus);
    });
    const clear = this.createClearButton(() => {
      held = null;
      picker.selected.value = null;
      picker.resync();
      this.value = null;
      input.focus();
    });
    this.effect(() => clear.hidden = empty(this.value) && picker.text.value === '');
    const follow = this.options.open;
    if (follow === undefined)
      return div([picker.root, clear], 'u2-pick-box');
    const open = iconButton('external-link', () => {
      const id = this.peek();
      if (!empty(id))
        follow(id!);
    }, {tooltip: 'Open'});
    open.classList.add('u2-pick-open');
    this.effect(() => open.hidden = empty(this.value));
    return div([picker.root, open, clear], 'u2-pick-box u2-pick-box-open');
  }
}
