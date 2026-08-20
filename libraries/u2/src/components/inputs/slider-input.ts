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
  private _range!: HTMLInputElement;

  constructor(options: SliderInputOptions = {}) {
    super(options, null);
    this.root.dataset.u2 = 'slider-input';
  }

  protected createEditor(): HTMLElement {
    const min = this.options.min ?? 0;
    const max = this.options.max ?? 100;
    const range = document.createElement('input');
    this._range = range;
    range.type = 'range';
    range.className = 'u2-slider-track';
    range.min = String(min);
    range.max = String(max);
    range.step = String(this.options.step ?? (max > min ? (max - min) / 100 : 1));

    let echo = false;
    this.effect(() => {
      const value = this.value.value;
      if (!echo)
        range.value = String(value ?? min);
    });
    const onInput = () => {
      echo = true;
      try {
        this.value.value = Number(range.value);
      } finally {
        echo = false;
      }
      Tooltip.show(this.text, range);
    };
    range.addEventListener('input', onInput);
    this.own(() => range.removeEventListener('input', onInput));
    Tooltip.bind(range, () => this.text, this.scope);

    const editor = document.createElement('div');
    editor.className = 'u2-slider';
    editor.append(range);
    return editor;
  }

  /** What the tooltip shows — the formatted value, empty while there is none. */
  get text(): string {
    const value = this.value.peek();
    if (value === null)
      return '';
    const format = this.options.format;
    return format ? format(value) : String(value);
  }
}
