/* Tour — a guided walkthrough over UI that already exists. One fixed layer dims the page and an
   SVG mask punches the spotlight holes, so any number of holes is one declarative element that no
   stacking context can occlude. The layer is click-through: a step may ask the user to actually use
   the control it points at (`advanceOn`). The popup is anchored by floating-ui's autoUpdate, the
   same convention Overlay.show uses — nothing polls. */
import {autoUpdate, computePosition, flip, offset, shift} from '@floating-ui/dom';
import {Overlay} from '../../core/overlay.js';
import {Scope} from '../../core/scope.js';
import {signal, ReadonlySignal} from '../../core/signals.js';
import {button, div, span} from '../../core/elements.js';

export interface TourStep {
  /** A string is always a spec name, resolved as `[data-u2-name="…"]` inside `root`; pass an
   * element or a resolver for anything else (a CSS selector: `() => document.querySelector(…)`). */
  target: string | HTMLElement | (() => HTMLElement | null);
  content: string | HTMLElement | (() => HTMLElement);
  position?: 'top' | 'bottom' | 'left' | 'right';
  /** A second spotlight hole. */
  extra?: string | HTMLElement;
  /** Auto-advance: a signal advances the tour when it turns truthy, a promise when it resolves.
   * NEXT stays available either way. */
  advanceOn?: ReadonlySignal<unknown> | Promise<unknown>;
}

export interface TourOptions {
  steps: TourStep[];
  /** Query scope for `data-u2-name` targets; the document by default. */
  root?: ParentNode;
  onFinish?: (result: 'done' | 'skipped') => void;
  buttons?: {next?: string, skip?: string, done?: string};
}

const SVG_NS = 'http://www.w3.org/2000/svg';
const HOLE_PADDING = 4;
// the CSS property, not the SVG presentation attribute: only the property resolves var()
const DIM = 'var(--dg-tour-dim, rgba(0, 0, 0, 0.6))';
const FOCUSABLE = 'a[href], button:not([disabled]), input:not([disabled]), select:not([disabled]), ' +
  'textarea:not([disabled]), [tabindex]:not([tabindex="-1"])';

export class Tour {
  private static _current: Tour | undefined;
  private static _seq = 0;

  readonly step: ReadonlySignal<number>;
  readonly active: ReadonlySignal<boolean>;
  readonly overlay: HTMLElement;
  readonly popup: HTMLElement;

  private readonly _options: TourOptions;
  private readonly _steps: TourStep[];
  private readonly _scope = new Scope();
  private readonly _step = signal(0);
  private readonly _active = signal(false);
  private readonly _hole: SVGElement;
  private readonly _extraHole: SVGElement;
  private readonly _content = div([], 'u2-tour-content');
  private readonly _counter = span('', 'u2-tour-counter');
  private readonly _next: HTMLButtonElement;
  private readonly _nextLabel: string;
  private readonly _doneLabel: string;
  private readonly _warned = new Set<string>();
  private _stepScope: Scope | undefined;
  private _target: HTMLElement | null = null;
  private _extra: HTMLElement | null = null;
  private _stopAutoUpdate: () => void = () => {};
  private _restore: HTMLElement | undefined;

  private constructor(options: TourOptions) {
    this._options = options;
    this._steps = options.steps;
    this.step = this._step;
    this.active = this._active;

    const maskId = `u2-tour-mask-${++Tour._seq}`;
    this._hole = Tour._svg('rect', {class: 'u2-tour-hole', fill: 'black', rx: String(HOLE_PADDING)});
    this._extraHole = Tour._svg('rect', {class: 'u2-tour-hole', fill: 'black', rx: String(HOLE_PADDING)});
    const mask = Tour._svg('mask', {id: maskId});
    mask.append(Tour._svg('rect', {x: '0', y: '0', width: '100%', height: '100%', fill: 'white'}),
      this._hole, this._extraHole);
    const defs = Tour._svg('defs', {});
    defs.append(mask);
    const dim = Tour._svg('rect',
      {x: '0', y: '0', width: '100%', height: '100%', mask: `url(#${maskId})`});
    dim.style.fill = DIM;
    const svg = Tour._svg('svg', {});
    svg.append(defs, dim);

    this.overlay = div([], 'u2-tour-overlay');
    this.overlay.dataset.u2 = 'tour';
    this.overlay.append(svg);

    const labels = options.buttons ?? {};
    this._nextLabel = labels.next ?? 'NEXT';
    this._doneLabel = labels.done ?? 'DONE';
    const [skip, next] = Scope.runWith(this._scope, () => [
      button(labels.skip ?? 'SKIP', () => this.skip()),
      button(this._nextLabel, () => this.next(), {primary: true}),
    ]);
    skip.classList.add('u2-tour-skip');
    next.classList.add('u2-tour-next');
    this._next = next;
    this.popup = div([this._content, div([skip, this._counter, next], 'u2-tour-footer')], 'u2-tour-popup');
    this.popup.setAttribute('role', 'dialog');
    this.popup.tabIndex = -1;
  }

  /** Constructs the tour and shows its first renderable step. */
  static run(options: TourOptions): Tour {
    const tour = new Tour(options);
    tour._start();
    return tour;
  }

  next(): void {
    if (this._active.value)
      this._show(this._step.value + 1);
  }

  skip(): void {
    this.finish('skipped');
  }

  finish(result: 'done' | 'skipped' = 'done'): void {
    if (!this._active.value)
      return;
    this._active.value = false;
    if (Tour._current === this)
      Tour._current = undefined;
    this._stepScope?.dispose();
    this._stepScope = undefined;
    this._scope.dispose();
    this.overlay.remove();
    this.popup.remove();
    if (this._restore?.isConnected)
      this._restore.focus();
    this._options.onFinish?.(result);
  }

  private _start(): void {
    // two tours would fight over the dim layer and the keyboard; the running one ends first
    Tour._current?.finish('skipped');
    Tour._current = this;
    this._restore = (document.activeElement as HTMLElement | null) ?? undefined;
    this.overlay.style.zIndex = String(Overlay.nextLayer());
    this.popup.style.zIndex = String(Overlay.nextLayer());
    Overlay.host.append(this.overlay, this.popup);
    this._active.value = true;

    const onKeyDown = (e: Event) => this._onKeyDown(e as KeyboardEvent);
    document.addEventListener('keydown', onKeyDown, true);
    this._scope.own(() => document.removeEventListener('keydown', onKeyDown, true));
    this._scope.own(() => this._stopAutoUpdate());
    this._show(0);
  }

  private _show(index: number): void {
    this._stepScope?.dispose();
    this._stepScope = undefined;
    this._stopAutoUpdate();
    this._stopAutoUpdate = () => {};

    let at = index;
    while (at < this._steps.length && !this._enterStep(at))
      at++;
    if (at >= this._steps.length) {
      this.finish('done');
      return;
    }

    const step = this._steps[at];
    this._step.value = at;
    this._content.replaceChildren();
    if (typeof step.content === 'string')
      this._content.textContent = step.content;
    else
      this._content.append(typeof step.content === 'function' ? step.content() : step.content);
    this._counter.textContent = `${at + 1} / ${this._steps.length}`;
    this._next.textContent = at === this._steps.length - 1 ? this._doneLabel : this._nextLabel;

    const placement = step.position ?? 'bottom';
    this._place(placement);
    this._stopAutoUpdate = autoUpdate(this._target!, this.popup, () => this._place(placement),
      {animationFrame: true});

    const scope = new Scope();
    this._stepScope = scope;
    const advanceOn = step.advanceOn;
    if (advanceOn instanceof Promise)
      void advanceOn.then(() => this._advance(scope));
    else if (advanceOn) {
      scope.effect(() => {
        if (advanceOn.value)
          this._advance(scope);
      });
    }
    this._next.focus();
  }

  /** Auto-advance fires only while its own step is the current one. */
  private _advance(scope: Scope): void {
    if (this._stepScope === scope)
      this.next();
  }

  /** Makes `index` the step being shown — answers false, and points at nothing, when its target is
   * not on screen. */
  private _enterStep(index: number): boolean {
    const step = this._steps[index];
    const target = this._element(step.target);
    if (!Tour._renderable(target)) {
      const key = typeof step.target === 'string' ? step.target : `step ${index + 1}`;
      if (!this._warned.has(key)) {
        this._warned.add(key);
        console.warn(`u2: tour target "${key}" not found — step skipped`);
      }
      return false;
    }
    this._target = target;
    this._extra = step.extra === undefined ? null : this._element(step.extra);
    return true;
  }

  private _element(target: TourStep['target']): HTMLElement | null {
    if (typeof target === 'function')
      return target();
    if (typeof target !== 'string')
      return target;
    return (this._options.root ?? document).querySelector(`[data-u2-name="${target}"]`);
  }

  private _place(placement: NonNullable<TourStep['position']>): void {
    const target = this._target;
    if (!target || !target.isConnected) {
      this.next();
      return;
    }
    Tour._setHole(this._hole, target);
    Tour._setHole(this._extraHole, this._extra);
    void computePosition(target, this.popup, {
      strategy: 'fixed',
      placement,
      middleware: [offset(8), flip({padding: 4}), shift({padding: 4})],
    }).then(({x, y}) => {
      this.popup.style.left = `${x}px`;
      this.popup.style.top = `${y}px`;
    });
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const target = e.target as HTMLElement | null;
    if (e.key === 'Escape') {
      e.preventDefault();
      this.skip();
    } else if (e.key === 'Enter' && !Tour._isTextEntry(target)) {
      e.preventDefault();
      this.next();
    } else if (e.key === 'Tab' && target && this.popup.contains(target))
      this._trapTab(e);
  }

  private _trapTab(e: KeyboardEvent): void {
    const els = Tour._focusable(this.popup);
    if (els.length === 0)
      return;
    const active = document.activeElement;
    if (e.shiftKey && active === els[0]) {
      e.preventDefault();
      els[els.length - 1].focus();
    } else if (!e.shiftKey && active === els[els.length - 1]) {
      e.preventDefault();
      els[0].focus();
    }
  }

  private static _focusable(el: HTMLElement): HTMLElement[] {
    const els = Array.from(el.querySelectorAll(FOCUSABLE)) as HTMLElement[];
    return els.filter((x) => x.offsetParent !== null);
  }

  private static _isTextEntry(el: HTMLElement | null): boolean {
    return !!el && (/^(INPUT|TEXTAREA|SELECT)$/.test(el.tagName) || el.isContentEditable);
  }

  private static _renderable(el: HTMLElement | null): el is HTMLElement {
    if (!el || !el.isConnected)
      return false;
    const rect = el.getBoundingClientRect();
    return rect.width > 0 && rect.height > 0;
  }

  private static _setHole(hole: SVGElement, el: HTMLElement | null): void {
    const rect = el ? el.getBoundingClientRect() : new DOMRect();
    const pad = el ? HOLE_PADDING : 0;
    hole.setAttribute('x', String(rect.left - pad));
    hole.setAttribute('y', String(rect.top - pad));
    hole.setAttribute('width', String(rect.width + 2 * pad));
    hole.setAttribute('height', String(rect.height + 2 * pad));
  }

  private static _svg(tag: string, attrs: Record<string, string>): SVGElement {
    const el = document.createElementNS(SVG_NS, tag);
    for (const name in attrs)
      el.setAttribute(name, attrs[name]);
    return el;
  }
}
