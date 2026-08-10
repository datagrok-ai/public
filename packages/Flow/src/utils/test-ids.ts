/** Stable `data-testid` values for Flow UI surfaces: `ff-<area>-<thing>[-<dynamic-slug>]`.
 *  Build with `tid(...)`, stamp plain DOM elements with `setTid(el, ...)`. */

export const TID_PREFIX = 'ff';

/** Slugify a dynamic label into a stable, selector-safe token. */
export function tidSlug(s: string | number): string {
  return String(s ?? '')
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-+|-+$/g, '') || 'x';
}

/** Compose a namespaced data-testid value, slugifying every part. */
export function tid(...parts: Array<string | number>): string {
  return [TID_PREFIX, ...parts.map(tidSlug)].join('-');
}

/** Stamp a plain DOM element's `data-testid`; returns the element for chaining. */
export function setTid<T extends HTMLElement>(el: T, ...parts: Array<string | number>): T {
  el.dataset.testid = tid(...parts);
  return el;
}
