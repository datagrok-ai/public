/* Splitter — the interaction model is a port-and-adapt of VS Code's
   `src/vs/base/browser/ui/sash/sash.ts` (github.com/microsoft/vscode, main, fetched 2026-08-13), MIT.
   Kept: the sash as an absolutely positioned hit area wider than its visible line, pointer capture
   for the whole drag, the delayed `hover` class, the `active` class while dragging, and
   double-click as a reset gesture. Dropped: touch/gesture factories, linked and orthogonal sashes,
   sash states (disabled/minimum/maximum), and the layout-provider indirection — sizes are a signal
   here, so panels lay themselves out. Added: keyboard resizing, which upstream has no equivalent of. */

import {Control} from '../../core/component.js';
import {signal, Signal} from '../../core/signals.js';

export interface SplitterOptions {
  direction: 'horizontal' | 'vertical';
  sizes?: number[];
  minSize?: number;
}

const HOVER_DELAY_MS = 150;
const SASH_SIZE = 4;

interface Drag {
  index: number;
  startSizes: number[];
  startPos: number;
  total: number;
}

export class Splitter extends Control {
  readonly sizes: Signal<number[]>;

  private readonly _horizontal: boolean;
  private readonly _minSize: number;
  private readonly _panels: HTMLElement[] = [];
  private readonly _sashes: HTMLElement[] = [];
  private _drag: Drag | undefined;
  private _pointerPos = 0;
  private _frame = 0;
  private _hoverTimer = 0;

  constructor(panels: HTMLElement[], options: SplitterOptions) {
    super();
    this._horizontal = options.direction === 'horizontal';
    this._minSize = options.minSize ?? 60;

    const initial = options.sizes?.length === panels.length ?
      options.sizes : panels.map(() => 1 / panels.length);
    const sum = initial.reduce((a, b) => a + b, 0);
    this.sizes = signal(initial.map((s) => s / sum));

    this.root.classList.add('u2-splitter', this._horizontal ? 'u2-splitter-horizontal' : 'u2-splitter-vertical');
    this.root.dataset.u2 = 'splitter';

    for (const content of panels) {
      const panel = document.createElement('div');
      panel.className = 'u2-splitter-panel';
      panel.append(content);
      this._panels.push(panel);
      this.root.append(panel);
    }
    for (let i = 0; i < panels.length - 1; i++) {
      const sash = this._createSash(i);
      this._sashes.push(sash);
      this.root.append(sash);
    }

    this.own(() => {
      window.clearTimeout(this._hoverTimer);
      if (this._frame)
        cancelAnimationFrame(this._frame);
    });
    this.effect(() => this._layout());
  }

  private _createSash(index: number): HTMLElement {
    const sash = document.createElement('div');
    sash.className = 'u2-sash';
    sash.tabIndex = 0;
    sash.setAttribute('role', 'separator');
    sash.setAttribute('aria-orientation', this._horizontal ? 'vertical' : 'horizontal');
    sash.setAttribute('aria-valuemin', '0');
    sash.setAttribute('aria-valuemax', '100');
    this._listen(sash, 'pointerdown', (e) => this._onPointerDown(e as PointerEvent, index));
    this._listen(sash, 'pointermove', (e) => this._onPointerMove(e as PointerEvent));
    this._listen(sash, 'pointerup', () => this._endDrag());
    this._listen(sash, 'pointercancel', () => this._endDrag());
    this._listen(sash, 'pointerenter', () => this._hover(sash, true));
    this._listen(sash, 'pointerleave', () => this._hover(sash, false));
    this._listen(sash, 'dblclick', () => this._resetPair(index));
    this._listen(sash, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent, index));
    return sash;
  }

  private _hover(sash: HTMLElement, on: boolean): void {
    window.clearTimeout(this._hoverTimer);
    if (on)
      this._hoverTimer = window.setTimeout(() => sash.classList.add('u2-sash-hover'), HOVER_DELAY_MS);
    else
      sash.classList.remove('u2-sash-hover');
  }

  private _onPointerDown(e: PointerEvent, index: number): void {
    if (e.button !== 0)
      return;
    e.preventDefault();
    const sash = this._sashes[index];
    sash.setPointerCapture(e.pointerId);
    sash.focus();
    sash.classList.add('u2-sash-active');
    this.root.classList.add('u2-splitter-dragging');
    this._drag = {index, startSizes: this.sizes.peek(), startPos: this._coord(e), total: this._total()};
  }

  private _onPointerMove(e: PointerEvent): void {
    if (!this._drag)
      return;
    this._pointerPos = this._coord(e);
    if (this._frame)
      return;
    this._frame = requestAnimationFrame(() => {
      this._frame = 0;
      this._applyDrag();
    });
  }

  private _endDrag(): void {
    if (!this._drag)
      return;
    if (this._frame) {
      cancelAnimationFrame(this._frame);
      this._frame = 0;
      this._applyDrag();
    }
    this._sashes[this._drag.index].classList.remove('u2-sash-active');
    this.root.classList.remove('u2-splitter-dragging');
    this._drag = undefined;
  }

  private _applyDrag(): void {
    const drag = this._drag!;
    this._resize(drag.index, drag.startSizes, (this._pointerPos - drag.startPos) / drag.total, drag.total);
  }

  private _onKeyDown(e: KeyboardEvent, index: number): void {
    const back = this._horizontal ? 'ArrowLeft' : 'ArrowUp';
    const forward = this._horizontal ? 'ArrowRight' : 'ArrowDown';
    if (e.key !== back && e.key !== forward)
      return;
    e.preventDefault();
    const step = (e.shiftKey ? 0.1 : 0.02) * (e.key === forward ? 1 : -1);
    this._resize(index, this.sizes.peek(), step, this._total());
  }

  private _resetPair(index: number): void {
    const sizes = this.sizes.peek().slice();
    const half = (sizes[index] + sizes[index + 1]) / 2;
    sizes[index] = half;
    sizes[index + 1] = half;
    this.sizes.value = sizes;
  }

  private _resize(index: number, from: number[], delta: number, total: number): void {
    const min = this._minSize / total;
    const lowest = min - from[index];
    const highest = from[index + 1] - min;
    if (lowest > highest)
      return;
    const applied = Math.max(lowest, Math.min(highest, delta));
    const sizes = from.slice();
    sizes[index] += applied;
    sizes[index + 1] -= applied;
    this.sizes.value = sizes;
  }

  private _layout(): void {
    const sizes = this.sizes.value;
    let offset = 0;
    for (let i = 0; i < this._panels.length; i++) {
      const percent = sizes[i] * 100;
      this._panels[i].style.flexBasis = `${percent}%`;
      offset += percent;
      const sash = this._sashes[i];
      if (!sash)
        continue;
      sash.style.setProperty(this._horizontal ? 'left' : 'top', `calc(${offset}% - ${SASH_SIZE / 2}px)`);
      sash.setAttribute('aria-valuenow', String(Math.round(offset)));
    }
  }

  private _total(): number {
    return this._horizontal ? this.root.clientWidth : this.root.clientHeight;
  }

  private _coord(e: PointerEvent): number {
    return this._horizontal ? e.clientX : e.clientY;
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
