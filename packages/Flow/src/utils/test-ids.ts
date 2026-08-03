/** Stable `data-testid` attributes for every Flow UI surface, so the canvas,
 *  panels, ribbon, and dynamic lists can be addressed deterministically from
 *  UI tests (Playwright's `getByTestId()` reads `data-testid` by default).
 *
 *  Convention: values are `ff-<area>-<thing>[-<dynamic-slug>]`. Static parts are
 *  written literally; dynamic parts (a function name, a param name, a node type)
 *  are slugified so the id is selector-safe and reproducible. Examples:
 *    ff-node · ff-browser-search · ff-browser-item-<funcSlug> ·
 *    ff-prop-input-<paramSlug> · ff-socket-input-<keySlug>
 *
 *  Use `tid(...)` to build a value and `setTid(el, ...)` to stamp a plain DOM
 *  element; React components take the result of `tid(...)` as a `data-testid` prop. */

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

/** Stamp a plain DOM element's `data-testid` (camelCased dataset key maps to
 *  the `data-testid` attribute). Returns the element for chaining. */
export function setTid<T extends HTMLElement>(el: T, ...parts: Array<string | number>): T {
  el.dataset.testid = tid(...parts);
  return el;
}
