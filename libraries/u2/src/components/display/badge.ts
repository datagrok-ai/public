/* Badges: status pills, count bubbles, status dots and removable tags. Text and count slots take
   `T | Signal<T>` like the element factories — a signal binds through the ambient Scope, so live
   values never outlive their owner. */
import {ReadonlySignal} from '../../core/signals.js';
import {Scope} from '../../core/scope.js';
import {span, Text} from '../../core/elements.js';

export type BadgeVariant = 'default' | 'accent' | 'success' | 'warning' | 'error';

function owner(): Scope {
  const scope = Scope.ambient;
  if (!scope) {
    throw new Error('u2: signal binding needs an owner — ' +
      'wrap the code in Control.build(...) or component.run(...)');
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
 * badge instead of rendering an empty circle. */
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
    owner().effect(() => apply(count.value));
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
