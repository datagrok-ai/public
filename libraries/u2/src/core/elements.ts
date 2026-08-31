/* Fluent element factories (L3). Every text/child slot takes `T | ReadonlySignal<T>`; a signal
   is bound through the ambient Scope, so live values never outlive their owner. */
import {Signal, ReadonlySignal} from './signals.js';
import {Scope} from './scope.js';
import {Control} from './component.js';
import {bindText} from './bind.js';

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

/** Listeners on a factory-made element die with the element; owning them is still right when
 * someone is building inside a scope, so disposal is provably listener-free. */
function listen(el: HTMLElement, type: string, handler: (e: Event) => void): void {
  el.addEventListener(type, handler);
  const scope = Scope.ambient;
  if (scope)
    scope.own(() => el.removeEventListener(type, handler));
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

/** Compact timestamp, full truth on hover (docs/recipes/list-item-rendering.md): renders the
 * locale short date — the year appears only when it is not the current one — with the full
 * date-time in the tooltip. Accepts anything date-shaped (`dayjs` satisfies {@link DateLike});
 * an invalid date renders empty. */
export function timestamp(value: Date | number | string | DateLike, cls?: string): HTMLSpanElement {
  const date = value instanceof Date ? value :
    typeof value === 'object' && value !== null ? value.toDate() :
      new Date(value);
  const el = element('span', 'u2-timestamp', cls);
  if (!(date instanceof Date) || isNaN(date.getTime()))
    return el;
  const withYear = date.getFullYear() !== new Date().getFullYear();
  el.textContent = date.toLocaleDateString(undefined,
    {month: 'short', day: 'numeric', ...(withYear ? {year: 'numeric'} : {})});
  el.title = date.toLocaleString(undefined,
    {month: 'short', day: 'numeric', year: 'numeric', hour: '2-digit', minute: '2-digit'});
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

export function link(text: Text, onClickOrUrl: (() => void) | string): HTMLAnchorElement {
  const el = setText(element('a', 'u2-link'), text);
  if (typeof onClickOrUrl === 'string') {
    el.classList.add('u2-link-external');
    el.href = onClickOrUrl;
    el.target = '_blank';
    return el;
  }
  el.tabIndex = 0;
  listen(el, 'click', () => onClickOrUrl());
  listen(el, 'keydown', (e) => {
    if ((e as KeyboardEvent).key === 'Enter')
      onClickOrUrl();
  });
  return el;
}

export function button(text: Text | HTMLElement, onClick: () => void,
  options?: {primary?: boolean}): HTMLButtonElement {
  const el = element('button', 'u2-btn', options?.primary ? 'u2-btn-primary' : undefined);
  el.type = 'button';
  if (text instanceof HTMLElement)
    el.append(text);
  else
    setText(el, text);
  listen(el, 'click', () => onClick());
  return el;
}

export function bigButton(text: Text, onClick: () => void): HTMLButtonElement {
  return button(text, onClick, {primary: true});
}
