/* Progress bar: a thin track whose fill follows a 0..1 signal, plus the indeterminate variant for
   work with no known extent. */
import {signal, Signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {span} from '../../core/elements.js';

export interface ProgressBarOptions {
  indeterminate?: boolean;
  showPercent?: boolean;
  description?: string | ReadonlySignal<unknown>;
}

export class ProgressBar extends Control {
  readonly value: Signal<number>;

  private readonly _track = document.createElement('div');
  private readonly _fill = document.createElement('div');
  private readonly _percent?: HTMLElement;
  private _indeterminate: boolean;

  constructor(options?: ProgressBarOptions) {
    super();
    this.value = signal(0);
    this._indeterminate = options?.indeterminate ?? false;
    this.root.classList.add('u2-progress');
    this.root.dataset.u2 = 'progress-bar';

    this._track.className = 'u2-progress-track';
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
      row.append(this._percent);
    }
    this.root.append(row);

    const description = options?.description;
    if (description !== undefined)
      this.root.append(this.run(() => span(description, 'u2-progress-description')));

    this.root.classList.toggle('u2-progress-indeterminate', this._indeterminate);
    this.effect(() => this._apply(this.value.value));
  }

  get indeterminate(): boolean {
    return this._indeterminate;
  }

  set indeterminate(x: boolean) {
    this._indeterminate = x;
    this.root.classList.toggle('u2-progress-indeterminate', x);
    this._apply(this.value.peek());
  }

  private _apply(v: number): void {
    const clamped = Number.isFinite(v) ? Math.min(1, Math.max(0, v)) : 0;
    if (clamped !== v) {
      this.value.value = clamped;
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
