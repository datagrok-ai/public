import {Input, InputOptions} from '../../core/input-base.js';
import {signal, Signal} from '../../core/signals.js';
import {icon} from '../display/icon.js';

export interface NumberInputOptions extends InputOptions<number | null> {
  /** `int` accepts an optional minus and digits; `float` (default) also decimals and exponents. */
  mode?: 'int' | 'float';
  min?: number;
  max?: number;
  /** Spinner, clicker and arrow-key increment; 1 by default. */
  step?: number;
  /** Range control under the text box (Dart `showSlider`); renders only within real bounds. */
  slider?: boolean;
  /** − / + buttons on the options rail (Dart `showPlusMinus`). */
  clicker?: boolean;
  /** u2's own hover spinner inside the field; off by default, the platform has no counterpart. */
  spinner?: boolean;
  /** Display format for committed values; text being typed is never reformatted. */
  format?: (value: number) => string;
}

const INT = /^-?\d+$/;
const FLOAT = /^-?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/;

/** Numeric editor over a text field: `input type=number` brings a native spinner and its own
 * validation, both of which fight the token chrome. Text that does not parse stays visible and
 * only marks the input invalid; the value signal keeps the last good number.
 *
 * Out-of-range text behaves as the platform's `NumberInput` does (`number_input.dart:40-61`): the
 * typed number reaches the value signal and the min/max validators flag it, nothing is clamped
 * behind the user's back. Only the bounded controls — slider and stepper — stay inside the range,
 * and both work off the value signal, never off the displayed text a `format` may have rewritten
 * into something the parser rejects. On an empty field the stepper does nothing, as the platform's
 * clicker does (`number_input.dart:66-68` steps `null` to `null`).
 *
 * The chrome is the platform's, ported as it stands (`number_input.dart:65-92`, `d4.css:75-134`):
 * at rest a bare field, and hovering anywhere on the input reveals the − + pair on the options
 * rail (20×20, blue, left of the units) and the slider, absolutely positioned under the field at
 * its full width. u2's own hover spinner has no Dart counterpart and is opt-in (`spinner: true`).
 *
 * Ported quirks — the platform's behavior, deliberately not improved here:
 * - the revealed slider overlaps whatever is below it by 9px (`d4.css:81-89`: `bottom: -9px;
 *   z-index: 2`), so in a form it sits over the next row;
 * - at rest nothing hints that a slider or a clicker exists — they are `visibility: hidden` until
 *   the pointer enters the input (`d4.css:130-134`);
 * - the platform's clicker steps past min/max (`number_input.dart:66-68` adds `step` unguarded);
 *   u2 clamps it, which is the one divergence kept from wave 1. */
export class NumberInput extends Input<number | null, NumberInputOptions> {
  private _input!: HTMLInputElement;
  // assigned in createEditor, which the base constructor runs before field initializers
  private _slider!: HTMLInputElement | null;
  private _bad!: Signal<boolean>;

  constructor(options: NumberInputOptions = {}) {
    super(options, null);
    this.root.dataset.u2 = 'number-input';
  }

  protected createEditor(): HTMLElement {
    this._bad = signal(false);
    this._slider = null;
    const input = document.createElement('input');
    this._input = input;
    input.type = 'text';
    input.autocomplete = 'off';
    input.inputMode = this._int ? 'numeric' : 'decimal';

    let echo = false;
    this.effect(() => {
      const value = this.value.value;
      if (echo)
        return;
      input.value = this._format(value);
      this._syncSlider();
    });
    this._listen(input, 'input', () => {
      const parsed = this._parse(input.value);
      this._bad.value = parsed === undefined;
      if (parsed === undefined)
        return;
      echo = true;
      try {
        this.value.value = parsed;
      } finally {
        echo = false;
      }
      this._syncSlider();
    });
    this._listen(input, 'blur', () => this._commit());
    this._listen(input, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this._listen(input, 'wheel', (e) => {
      if (document.activeElement !== input)
        return;
      e.preventDefault();
      this._step((e as WheelEvent).deltaY < 0 ? 1 : -1);
    }, {passive: false});

    const box = document.createElement('span');
    box.className = 'u2-number-box';
    box.append(input);
    if (this.options.spinner)
      box.append(this._createSpinner());

    const editor = document.createElement('div');
    editor.className = 'u2-number';
    editor.append(box);
    // the pair goes on the shared rail, as `showPlusMinus` does (`number_input.dart:70-71`)
    if (this.options.clicker) {
      this.addOptions(this._createClicker('minus', -1, 'Decrease'));
      this.addOptions(this._createClicker('plus', 1, 'Increase'));
    }
    const slider = this._createSlider();
    if (slider !== null)
      editor.append(slider);

    // registered here, before the base builds its validation effect, so the effect tracks `_bad`
    this.addValidator(() => this._bad.value ? 'Not a number' : null);
    this.addValidator((value) => this._outOfRange(value));
    return editor;
  }

  private get _int(): boolean {
    return this.options.mode === 'int';
  }

  private _createSpinner(): HTMLElement {
    const spinner = document.createElement('span');
    spinner.className = 'u2-number-spinner';
    spinner.append(this._createStep('▴', 1, 'Increase'), this._createStep('▾', -1, 'Decrease'));
    this._listen(spinner, 'mousedown', (e) => e.preventDefault());
    return spinner;
  }

  private _createStep(glyph: string, direction: number, label: string): HTMLElement {
    const step = document.createElement('button');
    step.type = 'button';
    step.className = 'u2-number-spin';
    step.textContent = glyph;
    step.setAttribute('aria-label', label);
    this._listen(step, 'click', () => this._step(direction));
    return step;
  }

  /** The platform draws its pair with the same two icons, not with text glyphs
   * (`number_input.dart:66-69`). */
  private _createClicker(glyph: string, direction: number, label: string): HTMLElement {
    const button = this._createStep('', direction, label);
    button.className = 'u2-number-click';
    button.append(icon(glyph));
    return button;
  }

  /** The slider materializes only within real bounds, as the platform's does
   * (`number_input.dart:81`); floats default its step to a hundredth of the range (`:83`). */
  private _createSlider(): HTMLInputElement | null {
    const {min, max, step} = this.options;
    if (this.options.slider !== true || min === undefined || max === undefined || max <= min)
      return null;
    const slider = document.createElement('input');
    slider.type = 'range';
    slider.className = 'u2-number-slider';
    slider.min = String(min);
    slider.max = String(max);
    slider.step = String(step ?? (this._int ? 1 : (max - min) / 100));
    this._listen(slider, 'input', () => {
      const value = Number(slider.value);
      this._set(this._int ? Math.round(value) : value);
    });
    this._slider = slider;
    this._syncSlider();
    return slider;
  }

  private _syncSlider(): void {
    if (this._slider === null)
      return;
    const value = this.value.peek();
    this._slider.value = value === null ? '' : String(value);
  }

  private _onKeyDown(e: KeyboardEvent): void {
    if (e.key === 'Enter') {
      this._commit();
      return;
    }
    const direction = e.key === 'ArrowUp' ? 1 : e.key === 'ArrowDown' ? -1 : 0;
    if (!direction)
      return;
    e.preventDefault();
    this._step(direction * (e.shiftKey ? 10 : 1));
  }

  private _step(multiplier: number): void {
    const current = this.value.peek();
    if (current === null)
      return;
    const next = current + multiplier * (this.options.step ?? 1);
    this._set(this._clamp(this._int ? Math.round(next) : Number(next.toPrecision(12))));
  }

  private _commit(): void {
    const parsed = this._parse(this._input.value);
    if (parsed !== undefined)
      this._set(parsed);
  }

  private _set(value: number | null): void {
    this._bad.value = false;
    this.value.value = value;
    this._input.value = this._format(this.value.peek());
    this._syncSlider();
  }

  private _parse(text: string): number | null | undefined {
    const trimmed = text.trim();
    if (!trimmed)
      return null;
    return (this._int ? INT : FLOAT).test(trimmed) ? Number(trimmed) : undefined;
  }

  private _clamp(value: number): number {
    const {min, max} = this.options;
    if (min !== undefined && value < min)
      return min;
    if (max !== undefined && value > max)
      return max;
    return value;
  }

  private _outOfRange(value: number | null): string | null {
    const {min, max} = this.options;
    if (value === null)
      return null;
    if (min !== undefined && value < min)
      return `Value must be at least ${min}`;
    if (max !== undefined && value > max)
      return `Value must be at most ${max}`;
    return null;
  }

  private _format(value: number | null): string {
    if (value === null)
      return '';
    const format = this.options.format;
    if (format !== undefined)
      return format(value);
    if (this._int)
      return String(value);
    const decimals = this._decimals(value);
    if (decimals > 0)
      return value.toFixed(decimals);
    const rounded = Number(value.toFixed(6));
    return String(rounded === 0 && value !== 0 ? value : rounded);
  }

  /** Dart's dynamic precision: as many decimals as the most precise of value/min/max/step carries
   * (`float_input.dart:21-24`). */
  private _decimals(value: number): number {
    const {min, max, step} = this.options;
    let decimals = NumberInput._decimalsOf(value);
    for (const bound of [min, max, step])
      decimals = Math.max(decimals, NumberInput._decimalsOf(bound));
    return decimals;
  }

  private static _decimalsOf(value: number | undefined): number {
    if (value === undefined || value === null)
      return 0;
    const text = String(value);
    // '1.5e-7' carries no decimals to count; the caller's toFixed(6) trim renders it instead
    if (text.includes('e') || text.includes('E'))
      return 0;
    const dot = text.lastIndexOf('.');
    return dot < 0 ? 0 : text.length - dot - 1;
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void,
    options?: AddEventListenerOptions): void {
    el.addEventListener(type, handler, options);
    this.own(() => el.removeEventListener(type, handler, options));
  }
}
