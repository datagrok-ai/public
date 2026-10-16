/* Focus across a structural rebuild (the designer's patches). Controls keep their elements with
   keyed regions (`each`, core/region.ts) and need none of this. `focusMark` names the focused
   element by the keys on the way down from the host (`data-u2-key`, or an input's `data-u2-name`,
   one per level of keyed nesting); `refocus` focuses the element the same keys name after the
   rebuild: the element itself where it is focusable, else its first focusable descendant. A keyed
   element gone from the rebuilt DOM (a removed row) hands focus to the neighbour that followed it,
   else the one before, else the nearest surviving level's first focusable. */
const FOCUSABLE = 'input, select, textarea, button, [tabindex]';

/** Where focus stood before a rebuild: the keys down to it, and at every level the keys of the
 * nearest keyed elements in order, the neighbours to fall back on. */
export interface FocusMark {
  path: string[];
  levels: string[][];
}

function keyOf(el: Element): string | undefined {
  const data = (el as HTMLElement).dataset;
  return data.u2Key ?? data.u2Name;
}

/** The keyed descendants nearest to `el`: a keyed element hides what is keyed below it. */
function nearestKeyed(el: Element): Element[] {
  const found: Element[] = [];
  for (const child of el.children) {
    if (keyOf(child) !== undefined)
      found.push(child);
    else
      found.push(...nearestKeyed(child));
  }
  return found;
}

function focusable(el: Element): boolean {
  return el.matches(FOCUSABLE) && (el as HTMLInputElement).disabled !== true;
}

function firstFocusable(el: Element): HTMLElement | undefined {
  return (focusable(el) ? el : [...el.querySelectorAll(FOCUSABLE)].find(focusable)) as HTMLElement | undefined;
}

/** The keys from `host` down to the focused element; null while focus is elsewhere, empty for
 * a focused element nothing keys. */
export function focusPath(host: Element): string[] | null {
  return focusMark(host)?.path ?? null;
}

export function focusMark(host: Element): FocusMark | null {
  const active = document.activeElement;
  if (active === null || active === host || !host.contains(active))
    return null;
  const keyed: Element[] = [];
  for (let el: Element | null = active; el !== null && el !== host; el = el.parentElement) {
    if (keyOf(el) !== undefined)
      keyed.unshift(el);
  }
  const levels: string[][] = [];
  let parent: Element = host;
  for (const el of keyed) {
    levels.push(nearestKeyed(parent).map((c) => keyOf(c)!));
    parent = el;
  }
  return {path: keyed.map((el) => keyOf(el)!), levels};
}

/** The element `path` names under `el`, from level `from` on; undefined where a key is missing. */
function resolve(el: Element, path: string[], from: number): HTMLElement | undefined {
  for (let i = from; i < path.length; i++) {
    const next = nearestKeyed(el).find((c) => keyOf(c) === path[i]);
    if (next === undefined)
      return undefined;
    el = next;
  }
  return firstFocusable(el);
}

/** Focuses the element the mark (or a bare path) names under `host`; a key no longer there
 * falls back to its neighbours, then to the first focusable of the nearest surviving level. */
export function refocus(host: Element, mark: FocusMark | string[] | null): void {
  if (mark === null)
    return;
  const path = Array.isArray(mark) ? mark : mark.path;
  const levels = Array.isArray(mark) ? [] : mark.levels;
  if (path.length === 0)
    return;
  const survivors: Element[] = [host];
  for (let i = 0; i < path.length; i++) {
    const el = survivors[survivors.length - 1];
    const keyed = nearestKeyed(el);
    const next = keyed.find((c) => keyOf(c) === path[i]);
    if (next !== undefined) {
      survivors.push(next);
      continue;
    }
    const old = levels[i] ?? [];
    const at = old.indexOf(path[i]);
    const neighbours = [...old.slice(at + 1), ...old.slice(0, Math.max(at, 0)).reverse()];
    for (const key of neighbours) {
      const candidate = keyed.find((c) => keyOf(c) === key);
      const target = candidate === undefined ? undefined : resolve(candidate, path, i + 1) ?? firstFocusable(candidate);
      if (target !== undefined)
        return target.focus();
    }
    break;
  }
  for (const el of survivors.reverse()) {
    const target = firstFocusable(el);
    if (target !== undefined)
      return target.focus();
  }
}
