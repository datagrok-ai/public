import {Input, InputOptions} from '../../core/input-base.js';
import {signal, Signal} from '../../core/signals.js';
import {text} from '../../core/text.js';

const INTEGER = /^[+-]?\d+$/;

export interface BigIntInputOptions extends InputOptions<bigint | null> {
  placeholder?: string;
}

/** Editor for an integer of any size, the port of Dart's `BigIntInput` (`big_int_input.dart:3-31`).
 * The text is the storage — a `bigint` is only ever made from it and printed back — so digits are
 * never routed through a double and nothing is lost, however long the number is.
 *
 * Dart reports null for text that is not an integer (`big_int_input.dart:29` throws on it, in
 * fact); u2 keeps the text, marks the input invalid and leaves the value at the last good one. */
export class BigIntInput extends Input<bigint | null, BigIntInputOptions> {
  // assigned in createEditor, which the base constructor runs before field initializers
  private declare _input: HTMLInputElement;
  private declare _bad: Signal<boolean>;

  constructor(options: BigIntInputOptions = {}) {
    super(options, null);
    this.root.dataset.u2 = 'bigint-input';
  }

  /** The integer in `text`, null when it is blank, undefined when it is not an integer. */
  static parse(text: string): bigint | null | undefined {
    const trimmed = text.trim();
    if (trimmed === '')
      return null;
    return INTEGER.test(trimmed) ? BigInt(trimmed) : undefined;
  }

  protected createEditor(): HTMLElement {
    this._bad = signal(false);
    const input = document.createElement('input');
    this._input = input;
    input.type = 'text';
    input.autocomplete = 'off';
    input.inputMode = 'numeric';
    input.placeholder = this.options.placeholder ?? '';

    // set by an edit that changed the value, consumed by the effect run it causes — after the
    // handler, when its batch ends — so the text being typed is never rewritten
    let echo = false;
    this.effect(() => {
      const value = this.value;
      if (echo)
        echo = false;
      else
        input.value = text(value);
    });
    this.scope.listen(input, 'input', () => {
      const parsed = BigIntInput.parse(input.value);
      this._bad.value = parsed === undefined;
      if (parsed === undefined || this.options.commitOn === 'change')
        return;
      const before = this.peek();
      this.value = parsed;
      echo = !Object.is(this.peek(), before);
    });
    this.scope.listen(input, 'blur', () => this._commit());
    this.scope.listen(input, 'keydown', (e) => {
      if ((e as KeyboardEvent).key === 'Enter')
        this._commit();
    });

    // registered here, before the base builds its validation effect, so the effect tracks `_bad`
    this.addValidator(() => this._bad.value ? 'Not an integer' : null);
    return input;
  }

  get unreadable(): boolean {
    return this._bad.value;
  }

  /** Writes the text on blur or Enter (`commitOn: 'change'` waits for it) and redraws the value,
   * so `+007` settles as `7`. */
  private _commit(): void {
    const parsed = BigIntInput.parse(this._input.value);
    if (parsed === undefined)
      return;
    this.value = parsed;
    this._input.value = text(this.peek());
  }
}
