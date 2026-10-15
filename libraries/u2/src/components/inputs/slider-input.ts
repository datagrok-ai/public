import {Input, InputOptions} from '../../core/input-base.js';
import {Tooltip} from '../../core/tooltip.js';

export interface SliderInputOptions extends InputOptions<number | null> {
  min?: number;
  max?: number;
  /** Track increment; `(max - min) / 100` by default, as Dart's `SliderInput` does. */
  step?: number;
  /** Tooltip text; the plain number by default. */
  format?: (value: number) => string;
}

/** Bare range control — the counterpart of Dart's `SliderInput`, which has no textbox of its own.
 * Ported quirk (`slider_input.dart:7-23`, `ui.css:562-567`): a naked track shows no number, so the
 * value only surfaces in the tooltip, on hover and while dragging. Where the exact number matters
 * as much as the range, `NumberInput` with `slider: true` is the control. */
export class SliderInput extends Input<number | null, SliderInputOptions> {
  private declare _range: HTMLInputElement;

  constructor(options: SliderInputOptions = {}) {
    super(options, null);
    this.root.dataset.u2 = 'slider-input';
  }

  get min(): number | undefined {
    return this.options.min;
  }

  set min(x: number | undefined) {
    this.options.min = x ?? undefined;
    this._bounds();
  }

  get max(): number | undefined {
    return this.options.max;
  }

  set max(x: number | undefined) {
    this.options.max = x ?? undefined;
    this._bounds();
  }

  get step(): number | undefined {
    return this.options.step;
  }

  set step(x: number | undefined) {
    this.options.step = x ?? undefined;
    this._bounds();
  }

  protected createEditor(): HTMLElement {
    const range = document.createElement('input');
    this._range = range;
    range.type = 'range';
    range.className = 'u2-slider-track';
    this._bounds();

    // set by a drag that changed the value, consumed by the effect run it causes
    let echo = false;
    this.effect(() => {
      const value = this.value.value;
      if (echo)
        echo = false;
      else
        range.value = String(value ?? this.options.min ?? 0);
    });
    const commit = () => {
      const before = this.value.peek();
      this.value.value = Number(range.value);
      echo = !Object.is(this.value.peek(), before);
    };
    // `commitOn: 'change'` writes on release; the tooltip follows the thumb either way
    this.scope.listen(range, 'input', () => {
      if (this.options.commitOn !== 'change')
        commit();
      Tooltip.show(this._format(Number(range.value)), range);
    });
    if (this.options.commitOn === 'change')
      this.scope.listen(range, 'change', commit);
    Tooltip.bind(range, () => this.text, this.scope);

    const editor = document.createElement('div');
    editor.className = 'u2-slider';
    editor.append(range);
    return editor;
  }

  private _bounds(): void {
    const min = this.options.min ?? 0;
    const max = this.options.max ?? 100;
    const range = this._range;
    range.min = String(min);
    range.max = String(max);
    range.step = String(this.options.step ?? (max > min ? (max - min) / 100 : 1));
    range.value = String(this.value.peek() ?? min);
  }

  /** What the tooltip shows — the formatted value, empty while there is none. */
  get text(): string {
    const value = this.value.peek();
    return value == null ? '' : this._format(value);
  }

  private _format(value: number): string {
    const format = this.options.format;
    return format ? format(value) : String(Math.round(value * 100) / 100);
  }
}
