import {Input, InputOptions} from '../core/input-base.js';
import {signal, Signal} from '../core/signals.js';

export interface NumberInputOptions extends InputOptions<number | null> {
  /** `int` accepts an optional minus and digits; `float` (default) also decimals and exponents. */
  mode?: 'int' | 'float';
  min?: number;
  max?: number;
  /** Spinner and arrow-key increment; 1 by default. */
  step?: number;
}

const INT = /^-?\d+$/;
const FLOAT = /^-?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$/;

/** Numeric editor over a text field: `input type=number` brings a native spinner and its own
 * validation, both of which fight the token chrome. Text that does not parse stays visible and
 * only marks the input invalid; the value signal keeps the last good number. */
export class NumberInput extends Input<number | null, NumberInputOptions> {
  private _input!: HTMLInputElement;
  private _bad!: Signal<boolean>;

  constructor(options: NumberInputOptions = {}) {
    super(options, null);
    this.root.dataset.u2 = 'number-input';
  }

  protected createEditor(): HTMLElement {
    this._bad = signal(false);
    const input = document.createElement('input');
    this._input = input;
    input.type = 'text';
    input.autocomplete = 'off';
    input.inputMode = this._int ? 'numeric' : 'decimal';

    let echo = false;
    this.effect(() => {
      const value = this.value.value;
      if (!echo)
        input.value = this._format(value);
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
    });
    this._listen(input, 'blur', () => this._commit());
    this._listen(input, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this._listen(input, 'wheel', (e) => {
      if (document.activeElement !== input)
        return;
      e.preventDefault();
      this._step((e as WheelEvent).deltaY < 0 ? 1 : -1);
    }, {passive: false});

    const editor = document.createElement('div');
    editor.className = 'u2-number';
    editor.append(input, this._createSpinner());
    // registered here, before the base builds its validation effect, so the effect tracks `_bad`
    this.addValidator(() => this._bad.value ? 'Not a number' : null);
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
    const parsed = this._parse(this._input.value);
    const current = parsed === undefined || parsed === null ? 0 : parsed;
    const next = current + multiplier * (this.options.step ?? 1);
    this._set(this._int ? Math.round(next) : Number(next.toPrecision(12)));
  }

  private _commit(): void {
    const parsed = this._parse(this._input.value);
    if (parsed !== undefined)
      this._set(parsed);
  }

  private _set(value: number | null): void {
    this._bad.value = false;
    this.value.value = value === null ? null : this._clamp(value);
    this._input.value = this._format(this.value.peek());
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

  private _format(value: number | null): string {
    if (value === null)
      return '';
    if (this._int)
      return String(value);
    const rounded = Number(value.toFixed(6));
    return String(rounded === 0 && value !== 0 ? value : rounded);
  }

  private _listen(el: HTMLElement, type: string, handler: (e: Event) => void,
    options?: AddEventListenerOptions): void {
    el.addEventListener(type, handler, options);
    this.own(() => el.removeEventListener(type, handler, options));
  }
}
