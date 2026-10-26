/* A reference picker with an `Input<string | null>` face: the value is the id a reference holds,
   the box shows the item's caption. The face over any `TypeAhead` — a domain table's rows, the
   platform's users or groups. */
import {Input, InputOptions} from '../../core/input-base.js';
import {div, span} from '../../core/elements.js';
import {TypeAhead} from './typeahead.js';

export interface PickInputOptions<T> extends InputOptions<string | null> {
  /** Builds the type-ahead the input is a face of. */
  typeAhead: () => TypeAhead<T>;
  idOf: (item: T) => string;
  /** The item behind a value written from outside — the id a form loaded; null when unknown. */
  resolve: (id: string) => Promise<T | null>;
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

  // runs from the base constructor, so the coupling state lives here rather than in fields the
  // subclass initializers would reset afterwards
  protected createEditor(): HTMLElement {
    const {typeAhead, idOf, resolve} = this.options;
    const picker = typeAhead();
    this._typeAhead = picker;
    // the last pick is held through typing, so blur can put its name back
    let held: T | null = null;
    let resolving: string | null = null;
    let gen = 0;
    const lookUp = async (id: string) => {
      const mine = ++gen;
      resolving = id;
      try {
        const item = await resolve(id);
        if (mine === gen && this.value.peek() === id && item !== null)
          picker.selected.value = item;
      } catch {
        // an unknown id keeps the value and shows nothing
      } finally {
        if (mine === gen)
          resolving = null;
      }
    };
    // a form seeds a string field with '' for an empty cell: empty either way, and never rewritten
    // to null on the input's own account — that write would dirty the row on render
    const empty = (id: string | null) => id === null || id === '';
    this.effect(() => {
      const id = this.value.value;
      const current = picker.selected.peek();
      if ((current === null ? null : idOf(current)) === id)
        return;
      if (empty(id))
        picker.selected.value = null;
      else
        void lookUp(id!);
    });
    this.effect(() => {
      const item = picker.selected.value;
      // a value still awaiting its item is not a clear
      if (item === null && resolving !== null)
        return;
      const id = item === null ? null : idOf(item);
      const current = this.value.peek();
      // a selection dropped by typing keeps the value; one dropped by clearing the text does not
      if (id !== null)
        held = item;
      else if (picker.text.peek() !== '')
        return;
      else
        held = null;
      if (id === null ? !empty(current) : id !== current)
        this.value.value = id;
    });
    const input = picker.root.querySelector('input')!;
    const hint = span('', 'u2-pick-hint');
    hint.dataset.u2Part = 'hint';
    hint.hidden = true;
    this._hint = hint;
    // text that is not a pick never survives a blur: the held pick's name comes back (and a line
    // says so until the next focus), else the box is cleared — text that looks chosen is worse
    // than none (an Enter still awaiting its candidates is left to land)
    const onBlur = () => {
      if (picker.isPickPending || picker.selected.peek() !== null || picker.text.peek() === '')
        return;
      if (held !== null) {
        const typed = picker.text.peek();
        picker.selected.value = held;
        hint.textContent = `No match for "${typed}" — kept ${picker.text.peek()}`;
        hint.hidden = false;
      } else
        picker.resync();
    };
    const onFocus = () => hint.hidden = true;
    input.addEventListener('blur', onBlur);
    input.addEventListener('focus', onFocus);
    this.own(() => {
      input.removeEventListener('blur', onBlur);
      input.removeEventListener('focus', onFocus);
    });
    const clear = document.createElement('button');
    clear.type = 'button';
    clear.className = 'u2-input-clear';
    clear.textContent = '✕';
    clear.setAttribute('aria-label', 'Clear');
    const onClear = () => {
      held = null;
      picker.selected.value = null;
      picker.resync();
      this.value.value = null;
      input.focus();
    };
    clear.addEventListener('click', onClear);
    this.own(() => clear.removeEventListener('click', onClear));
    this.effect(() => clear.hidden = empty(this.value.value) && picker.text.value === '');
    return div([picker.root, clear], 'u2-pick-box');
  }
}
