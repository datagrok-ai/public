/* Progress bar: a thin track whose fill follows a 0..1 signal, plus the indeterminate variant for
   work with no known extent. */
import {signal, Signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {span} from '../../core/elements.js';
import type {BindSource} from '../../core/widget-like.js';

export interface ProgressBarOptions {
  indeterminate?: boolean;
  showPercent?: boolean;
  description?: string | ReadonlySignal<unknown>;
}

export class ProgressBar extends Control {
  /** The fraction done (0..1) as a signal, for wiring. */
  readonly signal: Signal<number>;

  private readonly _track = document.createElement('div');
  private readonly _fill = document.createElement('div');
  private readonly _percent?: HTMLElement;
  private _indeterminate: boolean;

  constructor(options?: ProgressBarOptions) {
    super();
    this.signal = signal(0);
    this._indeterminate = options?.indeterminate ?? false;
    this.root.classList.add('u2-progress');
    this.root.dataset.u2 = 'progress-bar';

    this._track.className = 'u2-progress-track';
    this._track.dataset.u2Part = 'track';
    this._track.setAttribute('role', 'progressbar');
    this._track.setAttribute('aria-valuemin', '0');
    this._track.setAttribute('aria-valuemax', '1');
    this._fill.className = 'u2-progress-fill';
    this._track.append(this._fill);

    const row = document.createElement('div');
    row.className = 'u2-progress-row';
    row.append(this._track);
    if (options?.showPercent) {
      this._percent = span('', 'u2-progress-percent');
      this._percent.dataset.u2Part = 'percent';
      row.append(this._percent);
    }
    this.root.append(row);

    const description = options?.description;
    if (description !== undefined) {
      const label = this.build(() => span(description, 'u2-progress-description'));
      label.dataset.u2Part = 'label';
      this.root.append(label);
    }

    this.root.classList.toggle('u2-progress-indeterminate', this._indeterminate);
    this.root.setAttribute('aria-busy', String(this._indeterminate));
    this.effect(() => this._apply(this.value));
  }

  /** The fraction done, 0..1: followed when read inside an `effect` or a `computed`; set it to move the bar. */
  get value(): number {
    return this.signal.value;
  }

  set value(x: number) {
    this.signal.value = x;
  }

  /** The fraction done, not followed. */
  peek(): number {
    return this.signal.peek();
  }

  bindStep(name: string): Signal<unknown> | BindSource | null {
    return name === '' ? this.signal as Signal<unknown> : super.bindStep(name);
  }

  get indeterminate(): boolean {
    return this._indeterminate;
  }

  set indeterminate(x: boolean) {
    this._indeterminate = x;
    this.root.classList.toggle('u2-progress-indeterminate', x);
    this.root.setAttribute('aria-busy', String(x));
    this._apply(this.peek());
  }

  /** The base state plus the fraction done, the percent shown, and busy while indeterminate. */
  getState(): ReturnType<Control['getState']> {
    const percent = Math.round(this.peek() * 100);
    return {...super.getState(), value: this.peek(), text: this._indeterminate ? '' : `${percent}%`};
  }

  private _apply(v: number): void {
    const clamped = Number.isFinite(v) ? Math.min(1, Math.max(0, v)) : 0;
    if (clamped !== v) {
      this.value = clamped;
      return;
    }
    const percent = Math.round(clamped * 1000) / 10;
    if (this._percent)
      this._percent.textContent = this._indeterminate ? '' : `${Math.round(percent)}%`;
    if (this._indeterminate) {
      this._fill.style.width = '';
      this._track.removeAttribute('aria-valuenow');
      return;
    }
    this._fill.style.width = `${percent}%`;
    this._track.setAttribute('aria-valuenow', String(percent / 100));
  }
}
