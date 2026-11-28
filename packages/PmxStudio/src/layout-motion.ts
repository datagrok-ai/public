// How a page's panels move while it is customized (panels.ts): where a dragged panel goes (the nearest slot to the
// pointer, wherever it is released), the width a resized edge snaps to, and the motion itself: the others slide out
// of the way on a spring (a slight overshoot), the dragged one lifts and follows the pointer a touch behind, and on
// release settles into its slot from where it was (FLIP). Only transform and opacity move, frame by frame; with
// reduced motion everything moves at once.

/** A panel's box on screen (a DOMRect's corners). */
export interface Box {left: number, top: number, width: number, height: number}

/** The overshooting spring every move eases with, and how long it takes. */
export const SPRING = 'cubic-bezier(0.34, 1.4, 0.64, 1)';
export const SPRING_MS = 240;

/** The user asked for less motion: every move is instant. */
export const reducedMotion = () => typeof matchMedia === 'function' &&
  matchMedia('(prefers-reduced-motion: reduce)').matches;

/** Where a panel dragged to (x, y) goes among the others (`boxes`, in the page's order without it): before or after
 *  the nearest one, by which side of its centre the pointer is on (across a row: left or right; else above or
 *  below). Never needs the pointer over a panel: past the last, it goes last; above the first, first. */
export function insertionIndex(boxes: Box[], x: number, y: number): number {
  const n = nearestBox(boxes, x, y);
  return n.i < 0 ? 0 : n.i + (n.after ? 1 : 0);
}

/** The panel nearest (x, y) among `boxes` (-1: none) and whether the pointer is past its centre (insertionIndex). */
export function nearestBox(boxes: Box[], x: number, y: number): {i: number, after: boolean} {
  if (!boxes.length)
    return {i: -1, after: true};
  let best = 0;
  let bestD = Infinity;
  boxes.forEach((b, i) => {
    // distance to the box itself (0 inside), so a tall panel is near wherever along it the pointer is
    const dx = Math.max(b.left - x, 0, x - (b.left + b.width));
    const dy = Math.max(b.top - y, 0, y - (b.top + b.height));
    const d = dx * dx + dy * dy;
    if (d < bestD) {
      bestD = d;
      best = i;
    }
  });
  const b = boxes[best];
  const inRow = y >= b.top && y <= b.top + b.height;
  return {i: best, after: inRow ? x > b.left + b.width / 2 : y > b.top + b.height / 2};
}

/** The widths a panel takes, in twelfths. */
export const SPAN_STOPS = [3, 4, 6, 8, 12] as const;
export type SpanStop = typeof SPAN_STOPS[number];

/** The width an edge dragged to `cols` twelfths snaps to: the nearest of a quarter, a third, a half, two thirds and
 *  the row. */
export function snapSpan(cols: number): SpanStop {
  return SPAN_STOPS.reduce((a, b) => Math.abs(b - cols) < Math.abs(a - cols) ? b : a);
}

/** The stop `by` stops from `span` (snapped first): ← and → on a panel's edge, held at a quarter and the row. */
export function stepSpan(span: number, by: number): SpanStop {
  const i = SPAN_STOPS.indexOf(snapSpan(span)) + by;
  return SPAN_STOPS[Math.max(0, Math.min(SPAN_STOPS.length - 1, i))];
}

/** The width a panel dragged into a row takes (`row`: the widths of the others there, in twelfths): its own
 *  when it fits beside them; else, when the row has a gap, the widest stop that fills it (a whole-row panel beside a
 *  half is a half); else its own (a row of its own). */
export function fitSpan(row: number[], own: number): number {
  const free = 12 - row.reduce((a, b) => a + b, 0);
  return own <= free || free < SPAN_STOPS[0] ? own : [...SPAN_STOPS].reverse().find((s) => s <= free)!;
}

/** Moves `els` from where they were (`first`, their boxes before the page changed) to where they are now, on the
 *  spring: each starts offset back to its old place and slides home. An element new on the page (not in `first`)
 *  grows in. */
export function flip(els: HTMLElement[], first: Map<HTMLElement, Box>): void {
  if (reducedMotion())
    return;
  for (const el of els) {
    const was = first.get(el);
    const now = el.getBoundingClientRect();
    if (!was) {
      el.animate([{transform: 'scale(0.92)', opacity: 0}, {transform: 'none', opacity: 1}],
        {duration: SPRING_MS, easing: SPRING});
      continue;
    }
    const [dx, dy] = [was.left - now.left, was.top - now.top];
    const [sx, sy] = [now.width ? was.width / now.width : 1, now.height ? was.height / now.height : 1];
    if (Math.abs(dx) < 0.5 && Math.abs(dy) < 0.5 && Math.abs(sx - 1) < 0.01 && Math.abs(sy - 1) < 0.01)
      continue;
    el.animate([{transformOrigin: 'top left', transform: `translate(${dx}px, ${dy}px) scale(${sx}, ${sy})`},
      {transformOrigin: 'top left', transform: 'none'}], {duration: SPRING_MS, easing: SPRING});
  }
}

/** The boxes of `els` now, to flip them from after the page changes. */
export const boxesOf = (els: HTMLElement[]) => new Map(els.map((el) => [el, el.getBoundingClientRect() as Box]));

/** Shrinks and fades `el` out where it was (`box`, inside `parent`) after the page has already let it go: the page's
 *  state changes at once, so clicks in quick succession all count (QA 0.1.3: ten × hid five), and the fade runs on
 *  its own, out of the way of clicks. Nothing of it stays on `el` afterwards: shown again (Ctrl+Z, show), mid-fade or
 *  later, it is fully visible (it stayed at opacity 0 behind a forwards fill). */
export function fadeOut(el: HTMLElement, parent: HTMLElement, box: Box): void {
  if (reducedMotion() || typeof el.animate !== 'function')
    return;
  const was = el.getAttribute('style');
  Object.assign(el.style, {position: 'absolute', left: `${box.left}px`, top: `${box.top}px`, width: `${box.width}px`,
    height: `${box.height}px`, margin: '0', pointerEvents: 'none', zIndex: '4'});
  parent.append(el);
  const a = el.animate([{transform: 'none', opacity: 1}, {transform: 'scale(0.92)', opacity: 0}],
    {duration: 160, easing: 'ease-in', fill: 'forwards'});
  const end = () => {
    a.cancel();
    if (was === null)
      el.removeAttribute('style');
    else
      el.setAttribute('style', was);
    // (still the ghost: not shown again meanwhile)
    if (el.parentElement === parent)
      el.remove();
  };
  a.finished.then(end, end);
}

/** The element that scrolls `el`'s page (the window's when none of its ancestors does). */
export function scrollerOf(el: HTMLElement): HTMLElement {
  for (let p = el.parentElement; p; p = p.parentElement) {
    const o = getComputedStyle(p).overflowY;
    if ((o === 'auto' || o === 'scroll') && p.scrollHeight > p.clientHeight)
      return p;
  }
  return document.scrollingElement as HTMLElement ?? document.documentElement;
}

/** How far a drag near the top or bottom of the scroller (within `edge` pixels) scrolls it this frame: faster the
 *  closer it is; 0 away from both edges. */
export function autoScrollStep(y: number, top: number, bottom: number, edge = 56, max = 18): number {
  if (y < top + edge)
    return -Math.ceil(max * Math.min(1, (top + edge - y) / edge));
  if (y > bottom - edge)
    return Math.ceil(max * Math.min(1, (y - (bottom - edge)) / edge));
  return 0;
}
