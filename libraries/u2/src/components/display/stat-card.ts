/* KPI card: a large formatted value, its label, an optional signed delta and icon. Shares the
   card surface through a double class instead of extending Card — the two have the surface in
   common and nothing else. A `source` takes over the value area, so loading and error states
   come from the AsyncValue contract rather than from a hand-rolled spinner. */
import {Signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {div, span, Text} from '../../core/elements.js';
import {AsyncValue} from '../../core/async-value.js';
import {icon} from './icon.js';
import {skeleton} from './async-view.js';

export type StatValue = number | string;

const EM_DASH = '—';

export interface StatCardOptions {
  label: Text;
  value?: StatValue | ReadonlySignal<StatValue | undefined>;
  /** Wins over {@link value}: loading renders a skeleton, error the message, ready the value. */
  source?: AsyncValue<StatValue>;
  /** Default: numbers through `toLocaleString()`, strings as they are. */
  format?: (x: StatValue) => string;
  /** A fraction: 0.12 reads "+12%". */
  delta?: number | ReadonlySignal<number | undefined>;
  deltaFormat?: (d: number) => string;
  /** Negative is the good direction (error rates, latency). */
  deltaInverted?: boolean;
  icon?: string;
}

export class StatCard extends Control {
  private readonly _options: StatCardOptions;
  private readonly _valueEl: HTMLElement;
  private readonly _deltaEl: HTMLElement | undefined;

  constructor(options: StatCardOptions) {
    super();
    this._options = options;
    this.root.classList.add('u2-card', 'u2-stat-card');
    this.root.dataset.u2 = 'stat-card';
    this._valueEl = div(undefined, 'u2-stat-value');
    this._deltaEl = options.delta === undefined ? undefined : div(undefined, 'u2-stat-delta');
    this.runInScope(() => {
      if (options.icon !== undefined)
        this.root.append(icon(options.icon, {cls: 'u2-stat-icon'}));
      this.root.append(this._valueEl, span(options.label, 'u2-stat-label'));
      if (this._deltaEl)
        this.root.append(this._deltaEl);
    });
    this.effect(() => this._renderValue());
    if (this._deltaEl)
      this.effect(() => this._renderDelta());
  }

  private _renderValue(): void {
    const source = this._options.source;
    if (source === undefined) {
      const value = StatCard._current(this._options.value);
      this._valueEl.textContent = value === undefined ? EM_DASH : this._format(value);
      return;
    }
    const state = source.state.value;
    if (state.kind === 'ready')
      this._valueEl.textContent = this._format(state.value);
    else if (state.kind === 'error')
      this._valueEl.replaceChildren(div([state.message], 'u2-stat-error'));
    else
      this._valueEl.replaceChildren(skeleton(1));
  }

  private _renderDelta(): void {
    const el = this._deltaEl as HTMLElement;
    const delta = StatCard._current(this._options.delta);
    el.classList.remove('u2-stat-delta-up', 'u2-stat-delta-down');
    el.textContent = '';
    el.hidden = delta === undefined || !Number.isFinite(delta);
    if (el.hidden)
      return;
    const d = delta as number;
    if (d !== 0) {
      el.classList.add(d > 0 !== !!this._options.deltaInverted ? 'u2-stat-delta-up' : 'u2-stat-delta-down');
      el.append(icon(d > 0 ? 'arrow-up' : 'arrow-down', {cls: 'u2-stat-delta-icon'}));
    }
    const format = this._options.deltaFormat;
    el.append(span(format ? format(d) : `${d > 0 ? '+' : ''}${Math.round(d * 100)}%`));
  }

  private _format(x: StatValue): string {
    const format = this._options.format;
    if (format)
      return format(x);
    return typeof x === 'number' ? x.toLocaleString() : x;
  }

  private static _current<T>(x: T | ReadonlySignal<T> | undefined): T | undefined {
    return x instanceof Signal ? (x as ReadonlySignal<T>).value : x as T | undefined;
  }
}
