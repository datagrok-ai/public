/* Two-handle range selector — the counterpart of Dart's `RangeSlider` (d4 range_selector), DOM
   instead of SVG. Drag a handle to move one end, the band to move the whole window, click the
   track to pull the nearest handle; each handle is a keyboard slider. Vertical runs bottom-to-top. */
import {signal, batch, Signal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {div} from '../../core/elements.js';
import {Tooltip} from '../../core/tooltip.js';

export interface RangeSliderOptions {
  min?: number;
  max?: number;
  /** Snap increment for dragging; continuous by default. The keyboard falls back to
   * `(max - min) / 100`, as `SliderInput` does. */
  step?: number;
  /** Smallest allowed `hi - lo`; 0 by default. */
  minRange?: number;
  /** Bottom-to-top track. */
  vertical?: boolean;
  /** Initial value, or an external two-way signal the slider adopts. */
  lo?: number | Signal<number>;
  hi?: number | Signal<number>;
  onChanged?: (lo: number, hi: number) => void;
  /** Tooltip and `aria-valuetext`; the plain number by default. */
  format?: (value: number) => string;
}

type DragKind = 'lo' | 'hi' | 'band';

export class RangeSlider extends Control {
  readonly lo: Signal<number>;
  readonly hi: Signal<number>;

  private readonly options: RangeSliderOptions;
  /** Bumped when a bound (min/max/step/minRange) changes, so the render effect re-runs. */
  private readonly _layout = signal(0);
  private readonly _track = div([], 'u2-range-track');
  private readonly _band = div([], 'u2-range-band');
  private readonly _loHandle: HTMLElement;
  private readonly _hiHandle: HTMLElement;

  constructor(options: RangeSliderOptions = {}) {
    super();
    this.options = options;
    this.lo = options.lo instanceof Signal ? options.lo : signal(options.lo ?? this.min);
    this.hi = options.hi instanceof Signal ? options.hi : signal(options.hi ?? this.max);

    this.root.classList.add('u2-range-slider');
    this.root.classList.toggle('u2-range-slider-vertical', this.vertical);
    this.root.dataset.u2 = 'range-slider';
    this._loHandle = this._createHandle('lo', 'Minimum');
    this._hiHandle = this._createHandle('hi', 'Maximum');
    this.root.append(this._track, this._band, this._loHandle, this._hiHandle);

    this.effect(() => this._render());
    const onChanged = options.onChanged;
    if (onChanged) {
      let initial = true;
      this.effect(() => {
        const lo = this.lo.value;
        const hi = this.hi.value;
        if (initial)
          initial = false;
        else
          onChanged(lo, hi);
      });
    }

    const onPointerDown = (e: Event) => this._onPointerDown(e as PointerEvent);
    this.root.addEventListener('pointerdown', onPointerDown);
    this.own(() => this.root.removeEventListener('pointerdown', onPointerDown));
  }

  get min(): number {
    return this.options.min ?? 0;
  }

  set min(x: number) {
    this.options.min = x;
    this._relayout();
  }

  get max(): number {
    return this.options.max ?? 100;
  }

  set max(x: number) {
    this.options.max = x;
    this._relayout();
  }

  get step(): number | undefined {
    return this.options.step;
  }

  set step(x: number | undefined) {
    this.options.step = x;
    this._relayout();
  }

  get minRange(): number {
    return this.options.minRange ?? 0;
  }

  set minRange(x: number) {
    this.options.minRange = x;
    this._relayout();
  }

  get vertical(): boolean {
    return this.options.vertical ?? false;
  }

  private _relayout(): void {
    this._layout.value = this._layout.peek() + 1;
  }

  /** Writes both ends at once, clamped to the bounds with the order and `minRange` kept. */
  setValues(lo: number, hi: number): void {
    const clampedLo = Math.min(Math.max(lo, this.min), this.max);
    const clampedHi = Math.min(Math.max(Math.max(hi, clampedLo + this.minRange), this.min), this.max);
    batch(() => {
      this.lo.value = clampedLo;
      this.hi.value = clampedHi;
    });
  }

  private _format(value: number): string {
    return this.options.format ? this.options.format(value) : String(Math.round(value * 100) / 100);
  }

  private _createHandle(kind: DragKind, label: string): HTMLElement {
    const handle = div([], `u2-range-handle u2-range-handle-${kind}`);
    handle.tabIndex = 0;
    handle.setAttribute('role', 'slider');
    handle.setAttribute('aria-label', label);
    if (this.vertical)
      handle.setAttribute('aria-orientation', 'vertical');
    const onKeyDown = (e: Event) => this._onKeyDown(e as KeyboardEvent, kind);
    handle.addEventListener('keydown', onKeyDown);
    this.own(() => handle.removeEventListener('keydown', onKeyDown));
    Tooltip.bind(handle, () => this._format((kind === 'lo' ? this.lo : this.hi).peek()), this.scope);
    return handle;
  }

  /** Where a fraction of the track lands as a value; 0 is the left/bottom end. */
  private _valueAt(e: PointerEvent): number {
    const rect = this._track.getBoundingClientRect();
    const length = this.vertical ? rect.height : rect.width;
    if (length === 0)
      return this.min;
    const fraction = this.vertical ?
      (rect.bottom - e.clientY) / length : (e.clientX - rect.left) / length;
    return this.min + Math.min(Math.max(fraction, 0), 1) * (this.max - this.min);
  }

  private _snap(value: number): number {
    const step = this.options.step;
    return step && step > 0 ? this.min + Math.round((value - this.min) / step) * step : value;
  }

  private _keyStep(): number {
    const step = this.options.step;
    if (step && step > 0)
      return step;
    return this.max > this.min ? (this.max - this.min) / 100 : 1;
  }

  private _onKeyDown(e: KeyboardEvent, kind: DragKind): void {
    const step = this._keyStep();
    const current = kind === 'lo' ? this.lo.peek() : this.hi.peek();
    let value: number;
    switch (e.key) {
      case 'ArrowRight': case 'ArrowUp': value = current + step; break;
      case 'ArrowLeft': case 'ArrowDown': value = current - step; break;
      case 'PageUp': value = current + step * 10; break;
      case 'PageDown': value = current - step * 10; break;
      case 'Home': value = this.min; break;
      case 'End': value = this.max; break;
      default: return;
    }
    e.preventDefault();
    this._writeEnd(kind, value);
    this._fire('input');
    this._fire('change');
  }

  private _writeEnd(kind: DragKind, value: number): void {
    if (kind === 'lo')
      this.lo.value = Math.min(Math.max(value, this.min), this.hi.peek() - this.minRange);
    else
      this.hi.value = Math.max(Math.min(value, this.max), this.lo.peek() + this.minRange);
  }

  /** A grabbed handle drags relative to the grab point (no teleport, no snap jump); the whole
   * root height between the handles pans the band — the 4px band line is too thin a target to
   * demand a direct hit; only a press outside the selection jumps the nearer handle to it. */
  private _onPointerDown(e: PointerEvent): void {
    const target = e.target as HTMLElement;
    const startValue = this._valueAt(e);
    const startLo = this.lo.peek();
    const startHi = this.hi.peek();
    let kind: DragKind;
    let grabbed = false;
    if (this._loHandle.contains(target)) {
      kind = 'lo';
      grabbed = true;
    } else if (this._hiHandle.contains(target)) {
      kind = 'hi';
      grabbed = true;
    } else if (target === this._band || (startValue > startLo && startValue < startHi))
      kind = 'band';
    else
      kind = startValue <= startLo ? 'lo' : 'hi';
    const offset = grabbed ? (kind === 'lo' ? startLo : startHi) - startValue : 0;
    e.preventDefault();

    const apply = (ev: PointerEvent) => {
      const value = this._valueAt(ev);
      if (kind === 'band') {
        const span = startHi - startLo;
        const lo = Math.min(Math.max(this._snap(startLo + value - startValue), this.min),
          this.max - span);
        batch(() => {
          this.lo.value = lo;
          this.hi.value = lo + span;
        });
        Tooltip.show(`${this._format(lo)} – ${this._format(lo + span)}`, this._band);
      } else {
        this._writeEnd(kind, this._snap(value + offset));
        const handle = kind === 'lo' ? this._loHandle : this._hiHandle;
        Tooltip.show(this._format(kind === 'lo' ? this.lo.peek() : this.hi.peek()), handle);
      }
      this._fire('input');
    };
    const onMove = (ev: Event) => apply(ev as PointerEvent);
    const stop = () => {
      window.removeEventListener('pointermove', onMove);
      window.removeEventListener('pointerup', onUp);
      this.scope.disown(stop);
    };
    const onUp = () => {
      stop();
      Tooltip.hide();
      this._fire('change');
    };
    window.addEventListener('pointermove', onMove);
    window.addEventListener('pointerup', onUp);
    this.scope.own(stop);
    if (kind !== 'band' && !grabbed)
      apply(e);
  }

  private _fire(type: string): void {
    this.root.dispatchEvent(new CustomEvent(type));
  }

  private _render(): void {
    this._layout.value;
    const span = this.max - this.min;
    const fraction = (value: number) =>
      span <= 0 ? 0 : Math.min(Math.max((value - this.min) / span, 0), 1);
    const lo = this.lo.value;
    const hi = this.hi.value;
    const fLo = fraction(lo);
    const fHi = Math.max(fraction(hi), fLo);
    const pct = (f: number) => `${f * 100}%`;
    if (this.vertical) {
      this._band.style.bottom = pct(fLo);
      this._band.style.height = pct(fHi - fLo);
      this._loHandle.style.bottom = pct(fLo);
      this._hiHandle.style.bottom = pct(fHi);
    } else {
      this._band.style.left = pct(fLo);
      this._band.style.width = pct(fHi - fLo);
      this._loHandle.style.left = pct(fLo);
      this._hiHandle.style.left = pct(fHi);
    }
    this._aria(this._loHandle, lo, this.min, hi - this.minRange);
    this._aria(this._hiHandle, hi, lo + this.minRange, this.max);
  }

  private _aria(handle: HTMLElement, value: number, min: number, max: number): void {
    handle.setAttribute('aria-valuenow', String(value));
    handle.setAttribute('aria-valuemin', String(min));
    handle.setAttribute('aria-valuemax', String(max));
    handle.setAttribute('aria-valuetext', this._format(value));
  }
}
