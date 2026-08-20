import {Input, InputOptions} from '../../core/input-base.js';
import {signal, Signal} from '../../core/signals.js';
import {QNum} from '../../core/qnum.js';

export interface QNumInputOptions extends InputOptions<number | null> {
  placeholder?: string;
  /** Display format for the numeric part; the qualifier glyph is prepended to whatever it returns. */
  format?: (value: number) => string;
}

/** Editor for a qualified number — `<5.2`, `10`, `>1e3` — the port of Dart's `QNumInput`
 * (`qnum_input.dart:3-37`). The value signal holds the packed double `QNum` produces, so it can be
 * handed to a qnum column or a qnum property as it stands.
 *
 * Text that carries no recognisable number stays on screen and marks the input invalid, as the
 * rest of u2 does (the plan's ruling 1); Dart's editor reports null for it and shows nothing. */
export class QNumInput extends Input<number | null, QNumInputOptions> {
  private _input!: HTMLInputElement;
  // assigned in createEditor, which the base constructor runs before field initializers
  private _bad!: Signal<boolean>;

  constructor(options: QNumInputOptions = {}) {
    super(options, null);
    this.root.dataset.u2 = 'qnum-input';
  }

  protected createEditor(): HTMLElement {
    this._bad = signal(false);
    const input = document.createElement('input');
    this._input = input;
    input.type = 'text';
    input.autocomplete = 'off';
    input.placeholder = this.options.placeholder ?? '';

    let echo = false;
    this.effect(() => {
      const value = this.value.value;
      if (echo)
        return;
      input.value = QNum.format(value, this.options.format);
    });
    this._listen(input, 'input', () => {
      const text = input.value.trim();
      const parsed = text === '' ? null : QNum.parse(text);
      this._bad.value = text !== '' && parsed === null;
      if (this._bad.peek())
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
    this.addValidator(() => this._bad.value ? 'Not a qualified number' : null);
    return input;
  }

  /** Redraws the committed value on blur, so `> 10.` settles as `>10` — Dart reformats on every
   * `value` write (`qnum_input.dart:33-36`). */
  private _commit(): void {
    if (!this._bad.peek())
      this._input.value = QNum.format(this.value.peek(), this.options.format);
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
