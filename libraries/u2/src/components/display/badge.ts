/* Badges: status pills, count bubbles, status dots and removable tags. Text and count slots take
   `T | Signal<T>` like the element factories — a signal binds through the ambient Scope, so live
   values never outlive their owner. */
import {ReadonlySignal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {span, effect, Text} from '../../core/elements.js';
import {icon} from './icon.js';

export type BadgeVariant = 'default' | 'accent' | 'success' | 'warning' | 'error';

function owner(): Scope {
  const scope = Scope.ambient;
  if (!scope) {
    throw new Error('u2: signal binding needs an owner — ' +
      'wrap the code in Control.build(...) or component.runInScope(...)');
  }
  return scope;
}

export function badge(text: Text, options?: {variant?: BadgeVariant}): HTMLElement {
  const el = span(text, 'u2-badge');
  el.classList.add(`u2-badge-${options?.variant ?? 'default'}`);
  el.dataset.u2 = 'badge';
  return el;
}

/** Round counter for unread-style counts: anything above `max` reads `max+`, and zero hides the
 * badge instead of rendering an empty circle. A signal count is followed in a scope of the badge's
 * own, a child of the ambient one. */
export function countBadge(count: number | ReadonlySignal<number>, options?: {max?: number}): HTMLElement {
  const max = options?.max ?? 99;
  const el = span('', 'u2-count-badge');
  el.dataset.u2 = 'count-badge';
  const apply = (n: number) => {
    el.textContent = n > max ? `${max}+` : String(n);
    el.style.display = n === 0 ? 'none' : '';
  };
  if (typeof count === 'number')
    apply(count);
  else
    new Scope(owner(), {label: 'count-badge', lazy: true}).effect(() => apply(count.value));
  return el;
}

export function dot(variant?: BadgeVariant): HTMLElement {
  const el = span('', 'u2-dot');
  el.classList.add(`u2-dot-${variant ?? 'default'}`);
  el.dataset.u2 = 'dot';
  el.setAttribute('aria-hidden', 'true');
  return el;
}

export function tag(text: string, options?: {onRemove?: () => void}): HTMLElement {
  const el = span('', 'u2-tag');
  el.dataset.u2 = 'tag';
  el.append(span(text, 'u2-tag-text'));
  const onRemove = options?.onRemove;
  if (onRemove) {
    const remove = document.createElement('button');
    remove.type = 'button';
    remove.className = 'u2-tag-remove';
    remove.textContent = '✕';
    remove.setAttribute('aria-label', `Remove ${text}`);
    remove.addEventListener('click', () => onRemove());
    el.append(remove);
  }
  return el;
}

/** The line beside an action that says what it will do or why it cannot yet: `problem` while there is
 * one — in the error color, after a warning mark — else `otherwise`, quiet ("Saves into Peptides").
 * Two lines at most; the whole text is its tooltip. Followed in a scope of the line's own, a child of
 * the ambient one; screen readers hear it change (`role="status"`). */
export function statusLine(problem: ReadonlySignal<string | null>, otherwise?: Text): HTMLElement {
  const text = span('', 'u2-status-line-text');
  const mark = icon('exclamation-triangle', {cls: 'u2-status-line-mark'});
  mark.setAttribute('aria-hidden', 'true');
  const el = span('', 'u2-status-line');
  el.dataset.u2 = 'status-line';
  el.setAttribute('role', 'status');
  el.append(mark, text);
  effect(() => {
    const message = problem.value;
    const quiet = otherwise === undefined ? '' : typeof otherwise === 'string' ? otherwise : String(otherwise.value ?? '');
    el.classList.toggle('u2-status-line-problem', message !== null);
    mark.hidden = message === null;
    text.textContent = message ?? quiet;
    el.title = text.textContent;
  });
  return el;
}
