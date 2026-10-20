/* Fluent element factories (L3). Every text/child slot takes `T | ReadonlySignal<T>`; a signal
   is bound through the ambient Scope, so live values never outlive their owner. */
import {Signal, ReadonlySignal, batch} from './signals.js';
import {Scope} from './scope.js';
import {Control} from './component.js';
import {bindText} from './bind.js';
import {Dates} from './dates.js';
import {Tooltip} from './tooltip.js';

export type Child = HTMLElement | Control | string | ReadonlySignal<unknown>;
export type Text = string | ReadonlySignal<unknown>;

const NO_OWNER = 'u2: signal binding needs an owner — ' +
  'wrap the code in Control.build(...) or component.runInScope(...)';

function isSignal(x: unknown): x is ReadonlySignal<unknown> {
  return x instanceof Signal;
}

function owner(): Scope {
  const scope = Scope.ambient;
  if (!scope)
    throw new Error(NO_OWNER);
  return scope;
}

function element<K extends keyof HTMLElementTagNameMap>(tag: K,
  ...classes: (string | undefined)[]): HTMLElementTagNameMap[K] {
  const el = document.createElement(tag);
  const cls = classes.filter((c) => c).join(' ');
  if (cls)
    el.className = cls;
  return el;
}

function setText<T extends HTMLElement>(el: T, text: Text): T {
  if (isSignal(text))
    bindText(owner(), el, text);
  else
    el.textContent = text;
  return el;
}

function node(child: Child): HTMLElement | string {
  if (isSignal(child))
    return setText(element('span'), child);
  return Control.is(child) ? child.root : child;
}

function fill<T extends HTMLElement>(el: T, children?: Child[]): T {
  for (const child of children ?? [])
    el.append(node(child));
  return el;
}

/** Reports the throws of handlers on elements built with no owner. Never owns anything, so it
 * is never counted. */
const UNOWNED = new Scope(undefined, {label: 'unowned', lazy: true});

/** A factory-made element's handler: owned by whoever is building (removed on its dispose), and
 * like every u2 handler batched and contained — a throw is reported, never thrown at the browser. */
function listen<E extends Event>(el: HTMLElement, type: string, handler: (e: E) => void): void {
  const owner = Scope.ambient;
  if (owner) {
    owner.listen(el, type, handler as (e: Event) => void);
    return;
  }
  el.addEventListener(type, (e) => {
    try {
      batch(() => handler(e as E));
    } catch (error) {
      UNOWNED.reportError(error);
    }
  });
}

export function div(children?: Child[], cls?: string): HTMLDivElement {
  return fill(element('div', 'u2-div', cls), children);
}

export function divV(children?: Child[], cls?: string): HTMLDivElement {
  return fill(element('div', 'u2-div', 'u2-div-v', cls), children);
}

export function divH(children?: Child[], cls?: string): HTMLDivElement {
  return fill(element('div', 'u2-div', 'u2-div-h', cls), children);
}

export function panel(children?: Child[]): HTMLDivElement {
  return fill(element('div', 'u2-div', 'u2-panel'), children);
}

export function empty(): HTMLDivElement {
  return element('div', 'u2-div');
}

export function span(text: Text, cls?: string): HTMLSpanElement {
  return setText(element('span', cls), text);
}

export interface DateLike { toDate(): Date }

/** Timestamp, full truth on hover (docs/recipes/list-item-rendering.md): the locale date with its
 * year — the form a grid cell shows — and the time of day where the value carries one, with the
 * full date-time in the tooltip. Accepts anything date-shaped (`dayjs` satisfies {@link DateLike});
 * an invalid date renders empty. */
export function timestamp(value: Date | number | string | DateLike, cls?: string,
  options: {utcDates?: boolean} = {}): HTMLSpanElement {
  const date = value instanceof Date ? value :
    typeof value === 'object' && value !== null ? value.toDate() :
      new Date(value);
  const el = element('span', 'u2-timestamp', cls);
  if (!(date instanceof Date) || isNaN(date.getTime()))
    return el;
  const day: Intl.DateTimeFormatOptions = {month: 'short', day: 'numeric', year: 'numeric'};
  const full: Intl.DateTimeFormatOptions = {...day, hour: '2-digit', minute: '2-digit'};
  // a value that falls on midnight UTC is a DATE, not a moment: read in the local zone it would
  // be the day before for half the world (`2026-10-02T00:00:00Z` → "Oct 1, 2026, 08:00 PM"), so
  // it is both decided and printed in UTC. Everything else is a real moment and stays local.
  if (options.utcDates === true && Dates.isDateOnly(date)) {
    el.textContent = el.title = date.toLocaleDateString(undefined, {...day, timeZone: 'UTC'});
    return el;
  }
  const midnight = date.getHours() === 0 && date.getMinutes() === 0;
  el.textContent = midnight ? date.toLocaleDateString(undefined, day) : date.toLocaleString(undefined, full);
  el.title = date.toLocaleString(undefined, full);
  return el;
}

export function label(text: Text, cls?: string): HTMLLabelElement {
  return setText(element('label', 'u2-label', cls), text);
}

export function h1(text: Text): HTMLHeadingElement {
  return setText(element('h1'), text);
}

export function h2(text: Text): HTMLHeadingElement {
  return setText(element('h2'), text);
}

export function h3(text: Text): HTMLHeadingElement {
  return setText(element('h3'), text);
}

/** `onClick` gets the click, or the Enter keydown that activated the link. */
export function link(text: Text,
  onClickOrUrl: ((e: MouseEvent | KeyboardEvent) => void) | string): HTMLAnchorElement {
  const el = setText(element('a', 'u2-link'), text);
  if (typeof onClickOrUrl === 'string') {
    el.classList.add('u2-link-external');
    el.href = onClickOrUrl;
    el.target = '_blank';
    return el;
  }
  el.tabIndex = 0;
  listen<MouseEvent>(el, 'click', (e) => onClickOrUrl(e));
  listen<KeyboardEvent>(el, 'keydown', (e) => {
    if (e.key === 'Enter')
      onClickOrUrl(e);
  });
  return el;
}

/** Adds `handler` to `el`, owned by whoever is building — a row of `each()`, a card — and removed
 * with it; batched, and a throw reported rather than thrown. What a built element listens with:
 * `scope.listen` on a longer-lived scope would keep every element ever built alive. */
export function on<E extends Event = Event>(el: HTMLElement, type: string, handler: (e: E) => void): HTMLElement {
  listen<E>(el, type, handler);
  return el;
}

/** Runs `fn` now and again whenever a signal it read changes, owned by whoever is building — an `each()`
 * row, a card — and stopped with it: what a plain element built there follows its model with, as
 * {@link on} is what it listens with. */
export function effect(fn: () => void): void {
  owner().effect(fn);
}

export interface ButtonOptions {
  primary?: boolean;
  /** Hover text; a function is asked at hover time. */
  tooltip?: string | (() => string);
  /** `true`, or WHY the button is disabled — the reason becomes its tooltip, so the user learns what
   * to fix; `false`, `null` or `''` enable it. A signal is followed live. */
  disabled?: boolean | string | null | ReadonlySignal<boolean | string | null>;
}

export function button(text: Text | HTMLElement, onClick: (e: MouseEvent) => void,
  options?: ButtonOptions): HTMLButtonElement {
  const el = element('button', 'u2-btn', options?.primary ? 'u2-btn-primary' : undefined);
  el.type = 'button';
  el.dataset.u2 = 'button';
  if (text instanceof HTMLElement)
    el.append(text);
  else
    setText(el, text);
  listen<MouseEvent>(el, 'click', (e) => onClick(e));
  const disabled = options?.disabled;
  const apply = (x: boolean | string | null | undefined) => el.disabled = x === true || (typeof x === 'string' && x !== '');
  if (isSignal(disabled))
    owner().effect(() => apply(disabled.value as boolean | string | null));
  else
    apply(disabled);
  if (disabled !== undefined || options?.tooltip !== undefined) {
    const text = () => {
      const x = isSignal(disabled) ? disabled.peek() : disabled;
      const tip = options?.tooltip;
      return typeof x === 'string' && x !== '' ? x : typeof tip === 'function' ? tip() : tip ?? '';
    };
    // no owner to release the hover listeners with: the browser's own title, as iconButton does
    const scope = Scope.ambient;
    if (scope === undefined)
      el.title = text();
    else
      Tooltip.bind(el, text, scope);
  }
  return el;
}

export function bigButton(text: Text, onClick: (e: MouseEvent) => void): HTMLButtonElement {
  return button(text, onClick, {primary: true});
}
