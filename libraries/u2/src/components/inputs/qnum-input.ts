import {Input, InputOptions, LiveOption} from '../../core/input-base.js';
import {signal, Signal} from '../../core/signals.js';
import {QNum} from '../../core/qnum.js';

export interface QNumInputOptions extends InputOptions<number | null> {
  placeholder?: LiveOption<string>;
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
  // assigned in createEditor, which the base constructor runs before field initializers
  private declare _input: HTMLInputElement;
  private declare _bad: Signal<boolean>;

  constructor(options: QNumInputOptions = {}) {
    super(options, null);
    this.root.dataset.u2 = 'qnum-input';
    this.liveOption('placeholder', options.placeholder, (x) => this._input.placeholder = x);
  }

  protected createEditor(): HTMLElement {
    this._bad = signal(false);
    const input = document.createElement('input');
    this._input = input;
    input.type = 'text';
    input.autocomplete = 'off';

    // set by an edit that changed the value, consumed by the effect run it causes — after the
    // handler, when its batch ends — so the text being typed is never rewritten
    let echo = false;
    this.effect(() => {
      const value = this.value.value;
      if (echo)
        echo = false;
      else
        input.value = QNum.format(value ?? null, this.options.format);
    });
    this.scope.listen(input, 'input', () => {
      const parsed = this._parse();
      if (parsed === undefined || this.options.commitOn === 'change')
        return;
      const before = this.value.peek();
      this.value.value = parsed;
      echo = !Object.is(this.value.peek(), before);
    });
    this.scope.listen(input, 'blur', () => this._commit());
    this.scope.listen(input, 'keydown', (e) => {
      if ((e as KeyboardEvent).key === 'Enter')
        this._commit();
    });

    // registered here, before the base builds its validation effect, so the effect tracks `_bad`
    this.addValidator(() => this._bad.value ? 'Not a qualified number' : null);
    return input;
  }

  /** Redraws the committed value on blur, so `> 10.` settles as `>10` — Dart reformats on every
   * `value` write (`qnum_input.dart:33-36`). */
  private _commit(): void {
    const parsed = this._parse();
    if (parsed === undefined)
      return;
    this.value.value = parsed;
    this._input.value = QNum.format(this.value.peek() ?? null, this.options.format);
  }

  /** The packed number in the field, null when it is blank, undefined (and invalid) when it carries
   * no recognisable number. */
  private _parse(): number | null | undefined {
    const text = this._input.value.trim();
    const parsed = text === '' ? null : QNum.parse(text);
    this._bad.value = text !== '' && parsed === null;
    return this._bad.peek() ? undefined : parsed;
  }
}
