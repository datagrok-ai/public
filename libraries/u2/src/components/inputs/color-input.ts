import {Input, InputOptions} from '../../core/input-base.js';
import {signal, Signal} from '../../core/signals.js';

const HEX = /^#([0-9a-f]{3}|[0-9a-f]{6})$/i;

export interface ColorInputOptions extends InputOptions<string> {
  /** Swatch only, no hex field — Dart's `showOnlyColorBox`. */
  swatchOnly?: boolean;
}

/** Hex field plus a swatch that opens the browser's own picker: the platform palette lives in
 * Dart, and L0-L2 stay platform-free. Text that is not a hex color stays on screen and only marks
 * the input invalid — the value signal and the swatch keep the last good color.
 *
 * Layout is the platform's (`color_input.dart:6-9,49-50`): a 100px hex field with the swatch on
 * the options rail right of it, and `swatchOnly` hides the field and drops the rail underline, as
 * `ui-input-color-box-only` does (`ui.css:256-260,898-901`). */
export class ColorInput extends Input<string, ColorInputOptions> {
  private declare _text: HTMLInputElement;
  private declare _picker: HTMLInputElement;
  private declare _swatch: HTMLElement;
  private declare _bad: Signal<boolean>;

  constructor(options: ColorInputOptions = {}) {
    super(options, '#000000');
    this.root.dataset.u2 = 'color-input';
  }

  protected createEditor(): HTMLElement {
    this._bad = signal(false);
    const text = document.createElement('input');
    this._text = text;
    text.type = 'text';
    text.autocomplete = 'off';
    text.spellcheck = false;
    text.className = 'u2-color-text';
    text.placeholder = '#rrggbb';

    const swatch = document.createElement('button');
    this._swatch = swatch;
    swatch.type = 'button';
    swatch.className = 'u2-color-swatch';
    swatch.setAttribute('aria-label', typeof this.options.label === 'string' ? `Pick ${this.options.label}` : 'Pick a color');

    // the native picker is the popup; it is never the tab target, the swatch is
    const picker = document.createElement('input');
    this._picker = picker;
    picker.type = 'color';
    picker.className = 'u2-color-picker';
    picker.tabIndex = -1;
    picker.setAttribute('aria-hidden', 'true');

    // set by a typed write, dropped by the effect run it causes (handlers run batched)
    let echo = false;
    this.effect(() => {
      const value = this.value.value;
      if (echo)
        echo = false;
      else {
        // a value set from outside replaces what was typed, its "not a color" included
        text.value = value;
        this._bad.value = false;
      }
      this._paint(value);
    });
    this.scope.listen(text, 'input', () => {
      const typed = text.value.trim();
      this._bad.value = !HEX.test(typed);
      if (this._bad.peek() || typed === this.value.peek())
        return;
      echo = true;
      this.value.value = typed;
    });
    this.scope.listen(picker, 'input', () => {
      this._bad.value = false;
      this.value.value = picker.value;
    });
    this.scope.listen(swatch, 'click', () => picker.click());

    const box = document.createElement('span');
    box.className = 'u2-color-box';
    box.append(swatch, picker);
    this.addOptions(box);
    if (this.options.swatchOnly) {
      this.root.classList.add('u2-color-box-only');
      text.hidden = true;
    }
    // registered here, before the base builds its validation effect, so the effect tracks `_bad`
    this.addValidator(() => this._bad.value ? 'Not a color' : null);
    return text;
  }

  protected resetEditor(): void {
    this._bad.value = false;
    this._text.value = this.value.peek();
    this._paint(this.value.peek());
  }

  private _paint(color: string): void {
    this._swatch.style.background = color;
    this._picker.value = ColorInput.toHex6(color);
  }

  /** `#abc` → `#aabbcc`; the native picker accepts nothing else. */
  static toHex6(color: string): string {
    const short = /^#([0-9a-f])([0-9a-f])([0-9a-f])$/i.exec(color);
    return short ? `#${short[1]}${short[1]}${short[2]}${short[2]}${short[3]}${short[3]}`.toLowerCase() :
      color.toLowerCase();
  }
}
