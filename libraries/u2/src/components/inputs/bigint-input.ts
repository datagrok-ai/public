import {Input, InputOptions} from '../../core/input-base.js';
import {signal, Signal} from '../../core/signals.js';

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
  private _input!: HTMLInputElement;
  // assigned in createEditor, which the base constructor runs before field initializers
  private _bad!: Signal<boolean>;

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

    let echo = false;
    this.effect(() => {
      const value = this.value.value;
      if (echo)
        return;
      input.value = value === null ? '' : String(value);
    });
    this._listen(input, 'input', () => {
      const parsed = BigIntInput.parse(input.value);
      this._bad.value = parsed === undefined;
      if (parsed === undefined)
        return;
      echo = true;
      try {
        this.value.value = parsed;
      } finally {
        echo = false;
      }
    });
    this._listen(input, 'blur', () => this._commit());

    // registered here, before the base builds its validation effect, so the effect tracks `_bad`
    this.addValidator(() => this._bad.value ? 'Not an integer' : null);
    return input;
  }

  /** Redraws the committed value on blur, so `+007` settles as `7`. */
  private _commit(): void {
    if (!this._bad.peek()) {
      const value = this.value.peek();
      this._input.value = value === null ? '' : String(value);
    }
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
