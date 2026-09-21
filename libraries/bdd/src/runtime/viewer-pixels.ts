/* What a viewer's pixels say: repainted, ink, colors, the selection highlight, the value range
   and the color scale — every "than before" against the snapshot the last change took. A throw
   inside an `expect.poll` callback ends the poll, so a read the viewer may be between layouts of
   returns `false` and keeps the reason for the failure that follows. */
import {Page} from '@playwright/test';
import {expect, pollMs} from './patience.js';
import type {ElementRef} from './args.js';
import {reasonOf} from './failure.js';
import {AreaChange, AreaColor, AreaColors, CanvasChange, Range, RangeChange, ScaleChange, ScaleRange} from './viewer-runtime.js';
import {hitArea, onViewer} from './viewers.js';

const COLOR_MIN_PX = 10;
const SIGNIFICANT_PX = 30;
const HUE_TOLERANCE = 20;
const GREY_SATURATION = 0.15;
const LIGHTNESS_TOLERANCE = 0.15;
const HIGHLIGHT_FLOOR = 200;
const HIGHLIGHT_PER_ROW = 2;

export function canvasChange(page: Page, target: ElementRef): Promise<CanvasChange> {
  return onViewer(page, target, (el) => (window as any).__bdd.change(el), undefined);
}

/** Waits until the canvas differs from the last snapshot by at least `minPx` pixels. */
export async function expectRepainted(page: Page, target: ElementRef, minPx = 1): Promise<void> {
  await expect.poll(async () => (await canvasChange(page, target)).delta,
    {timeout: pollMs(10000), message: minPx > 1 ? `${target.phrase} did not repaint by ${minPx} pixels` : `${target.phrase} did not repaint`}).toBeGreaterThanOrEqual(minPx);
}

/** A hit area's own repaint: the pixels inside its rectangle that differ from the snapshot. */
export async function expectAreaRepainted(page: Page, target: ElementRef, area: string, minPx = 1): Promise<void> {
  await expect.poll(() => onViewer(page, target, (el, a) => (window as any).__bdd.areaDelta(el, a), area),
    {timeout: pollMs(10000), message: `the "${area}" area of ${target.phrase} did not repaint`}).toBeGreaterThanOrEqual(minPx);
}

/** The canvas is what it was at the snapshot, read after the frame a repaint would have landed
 * on (a viewer may run a render pass that draws the same picture — the mouse-over row does). */
export async function expectNotRepainted(page: Page, target: ElementRef): Promise<void> {
  const still: {renders: number; delta: number; where: string} = await onViewer(page, target, (el) => (window as any).__bdd.stillness(el), undefined);
  expect(still.delta, `${target.phrase} repainted: ${still.delta} px changed in ${still.renders} render(s)${still.where ? ', ' + still.where : ''}`).toBe(0);
}

export async function expectInk(page: Page, target: ElementRef, compare: 'less' | 'more' | 'some'): Promise<void> {
  if (compare === 'some') {
    let why = '';
    try {
      await expect.poll(() => onViewer(page, target, (el) => (window as any).__bdd.ink(el), undefined).catch((e: Error) => { why = reasonOf(e); return 0; }),
        {timeout: pollMs(10000)}).toBeGreaterThan(0);
    }
    catch {
      throw new Error(`${target.phrase} is blank${why === '' ? '' : `: ${why}`}`);
    }
    return;
  }
  await expect.poll(async () => {
    const c = await canvasChange(page, target);
    return compare === 'less' ? c.ink < c.inkBefore : c.ink > c.inkBefore;
  }, {timeout: pollMs(10000), message: `${target.phrase} does not have ${compare} ink than before`}).toBe(true);
}

/** A hit area's ink against the snapshot before the last change. */
export async function expectAreaInk(page: Page, target: ElementRef, area: string, compare: 'less' | 'more'): Promise<void> {
  let last: AreaChange | undefined;
  let failure = '';
  const holds = async (): Promise<boolean> => {
    try {
      const c: AreaChange = last = await onViewer(page, target, (el, a) => (window as any).__bdd.areaChange(el, a), area);
      return compare === 'less' ? c.ink < c.inkBefore : c.ink > c.inkBefore;
    }
    catch (e) {
      failure = reasonOf(e);
      return false;
    }
  };
  try {
    await expect.poll(holds, {timeout: pollMs(10000)}).toBe(true);
  }
  catch {
    throw new Error(`the "${area}" area of ${target.phrase} does not have ${compare} ink than before` +
      (last ? ` (${last.inkBefore} px before, ${last.ink} now, in ${JSON.stringify(last.rect)} of the bitmap)` : failure ? `: ${failure}` : ''));
  }
}

export async function expectAreaPainted(page: Page, target: ElementRef, area: string): Promise<void> {
  await hitArea(page, target, area);
  await expect.poll(() => onViewer(page, target, (el, a) => (window as any).__bdd.areaInk(el, a), area),
    {timeout: pollMs(5000), message: `the "${area}" area of ${target.phrase} is blank`}).toBeGreaterThan(0);
}

export async function expectPalette(page: Page, target: ElementRef, min: number): Promise<void> {
  await expect.poll(() => onViewer(page, target, (el) => (window as any).__bdd.palette(el, 500), undefined),
    {timeout: pollMs(5000), message: `${target.phrase} is painted in fewer than ${min} colors`}).toBeGreaterThanOrEqual(min);
}

// --- colors ----------------------------------------------------------------------------------------

export function hsl(color: string): {h: number; s: number; l: number} {
  const c = parseInt(color.slice(1), 16);
  const r = ((c >> 16) & 255) / 255;
  const g = ((c >> 8) & 255) / 255;
  const b = (c & 255) / 255;
  const max = Math.max(r, g, b);
  const min = Math.min(r, g, b);
  const l = (max + min) / 2;
  const d = max - min;
  const s = d === 0 ? 0 : d / (1 - Math.abs(2 * l - 1));
  const h = d === 0 ? 0 : max === r ? 60 * (((g - b) / d) % 6) : max === g ? 60 * ((b - r) / d + 2) : 60 * ((r - g) / d + 4);
  return {h: (h + 360) % 360, s, l};
}

/** Two colors are the same paint: markers are drawn with alpha and edges anti-aliased, so a
 * color lands on the canvas as its hue at less saturation — the hue is what survives. Greys
 * (no hue) compare by lightness. */
export function near(a: string, b: string): boolean {
  const x = hsl(a);
  const y = hsl(b);
  if (x.s < GREY_SATURATION || y.s < GREY_SATURATION)
    return x.s < GREY_SATURATION && y.s < GREY_SATURATION && Math.abs(x.l - y.l) <= LIGHTNESS_TOLERANCE;
  const dh = Math.abs(x.h - y.h);
  return Math.min(dh, 360 - dh) <= HUE_TOLERANCE;
}

/** `#rrggbb`, upper-cased, from a color a feature wrote. */
export function parseHex(color: string): string {
  const want = '#' + color.replace(/^#/, '').toUpperCase();
  if (!/^#[0-9A-F]{6}$/.test(want))
    throw new Error(`"${color}" is not a #rrggbb color`);
  return want;
}

function areaColors(page: Page, target: ElementRef, area: string): Promise<AreaColors> {
  return onViewer(page, target, (el, a) => (window as any).__bdd.areaColors(el, a), area);
}

const describeColors = (read: AreaColors | undefined): string => read ?
  `its colors: ${read.colors.slice(0, 5).map((c) => `${c.hex} (${c.count})`).join(', ') || 'none'} in ${JSON.stringify(read.rect)} of a ${read.bitmap.join('x')} bitmap` : '';

const pixelsNear = (colors: AreaColor[], want: string): number => colors.filter((c) => near(c.hex, want)).reduce((n, c) => n + c.count, 0);

/** The color (within a shade of anti-aliasing) covers some pixels of the area. */
export async function expectAreaColor(page: Page, target: ElementRef, area: string, color: string): Promise<void> {
  const want = parseHex(color);
  await hitArea(page, target, area);
  let last: AreaColors | undefined;
  try {
    await expect.poll(async () => pixelsNear((last = await areaColors(page, target, area)).colors, want), {timeout: pollMs(5000)}).toBeGreaterThanOrEqual(COLOR_MIN_PX);
  }
  catch {
    throw new Error(`the "${area}" area of ${target.phrase} is not painted in ${want}; ${describeColors(last)}`);
  }
}

/** No pixel of the color (nor a shade of it) inside the area, read once. */
export async function expectAreaNotColor(page: Page, target: ElementRef, area: string, color: string): Promise<void> {
  const want = parseHex(color);
  await hitArea(page, target, area);
  const read = await areaColors(page, target, area);
  const count = pixelsNear(read.colors, want);
  expect(count, `the "${area}" area of ${target.phrase} is painted in ${want} (${count} px); ${describeColors(read)}`).toBeLessThan(COLOR_MIN_PX);
}

/** The area is painted in at least `count` hues: colors covering some pixels each, grouped by
 * `near` (a hue and its anti-aliased shades are one), greys and white aside. */
export async function expectAreaColors(page: Page, target: ElementRef, area: string, count: number): Promise<void> {
  await hitArea(page, target, area);
  let last: AreaColors | undefined;
  const hues = async (): Promise<number> => {
    const read: AreaColors = last = await areaColors(page, target, area);
    const groups: string[] = [];
    for (const c of read.colors.filter((x) => x.count >= COLOR_MIN_PX && hsl(x.hex).s >= GREY_SATURATION))
      if (!groups.some((g) => near(g, c.hex)))
        groups.push(c.hex);
    return groups.length;
  };
  try {
    await expect.poll(hues, {timeout: pollMs(5000)}).toBeGreaterThanOrEqual(count);
  }
  catch {
    throw new Error(`the "${area}" area of ${target.phrase} is painted in fewer than ${count} hues; ${describeColors(last)}`);
  }
}

/** One area has a significant color the other has nothing near. */
const hasOwnColor = (mine: AreaColor[], theirs: AreaColor[]): boolean =>
  mine.some((c) => c.count >= SIGNIFICANT_PX && !theirs.some((d) => d.count >= COLOR_MIN_PX && near(c.hex, d.hex)));

/** Two areas painted in different colors (`same` = false), or in the same ones: every significant
 * color of either has a near color in the other (a linked color coding, a category and its swatch). */
export async function expectAreaColorsAlike(page: Page, target: ElementRef, a: string, b: string, same: boolean): Promise<void> {
  await hitArea(page, target, a);
  await hitArea(page, target, b);
  let shownA: AreaColors | undefined;
  let shownB: AreaColors | undefined;
  const holds = async (): Promise<boolean> => {
    const ra: AreaColors = shownA = await areaColors(page, target, a);
    const rb: AreaColors = shownB = await areaColors(page, target, b);
    return (hasOwnColor(ra.colors, rb.colors) || hasOwnColor(rb.colors, ra.colors)) !== same;
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the "${a}" and "${b}" areas of ${target.phrase} are painted in ${same ? 'different' : 'the same'} colors; "${a}" ${describeColors(shownA)}; "${b}" ${describeColors(shownB)}`);
  }
}

// --- the selection highlight ---------------------------------------------------------------------------

/** The highlight a selection must add or drop, in device pixels: a floor for a handful of rows,
 * a share per selected row for more, capped at a quarter of the view (markers overlap). */
export function highlightMargin(c: CanvasChange): number {
  return Math.min(Math.max(HIGHLIGHT_FLOOR, HIGHLIGHT_PER_ROW * c.selected) * c.dpr * c.dpr, c.viewPx / 4);
}

/** Pixels in the selection hue against the snapshot: `more`/`less` than before by the margin the
 * selection warrants, `some`, or `none`. */
export async function expectHighlight(page: Page, target: ElementRef, compare: 'more' | 'less' | 'some' | 'none'): Promise<void> {
  const wrong = {none: 'shows a selection highlight', some: 'shows no selection highlight',
    more: 'does not show more selection highlight than before', less: 'does not show less selection highlight than before'};
  let last: CanvasChange | undefined;
  const holds = async (): Promise<boolean> => {
    const c = last = await canvasChange(page, target);
    const margin = highlightMargin(c);
    return compare === 'none' ? c.hue === 0 : compare === 'some' ? c.hue > 0 : compare === 'more' ? c.hue >= c.hueBefore + margin : c.hue <= c.hueBefore - margin;
  };
  try {
    await expect.poll(holds, {timeout: pollMs(10000)}).toBe(true);
  }
  catch {
    throw new Error(`${target.phrase} ${wrong[compare]}` +
      (last ? ` (${last.hueBefore} px in the selection color before, ${last.hue} now, ${last.selected} rows selected, margin ${Math.round(highlightMargin(last))})` : ''));
  }
}

// --- the value range and the color scale ------------------------------------------------------------

/** Both axes of a viewport, when it has a horizontal one: a viewer that restores the Y window and
 * loses the X window is not showing the same range. */
function sameRange(a: Range, b: Range, eps: number): boolean {
  if (Math.abs(a.top - b.top) >= eps || Math.abs(a.bottom - b.bottom) >= eps)
    return false;
  const horizontal = [a.left, a.right, b.left, b.right].every((x) => typeof x === 'number' && isFinite(x));
  return !horizontal || (Math.abs(a.left - b.left) < eps && Math.abs(a.right - b.right) < eps);
}

function rangeArea(r: Range): number {
  const width = typeof r.left === 'number' && typeof r.right === 'number' && isFinite(r.left) && isFinite(r.right)
    ? Math.abs(r.right - r.left) : 1;
  return Math.abs(r.height) * (width || 1);
}

export function rememberRange(page: Page, target: ElementRef): Promise<void> {
  return onViewer(page, target, (el) => { (window as any).__bdd.rememberRange(el); }, undefined);
}

/** The value range equals the one remembered for this viewer type — across a close and a reopen. */
export async function expectRememberedRange(page: Page, target: ElementRef): Promise<void> {
  let last: RangeChange = {};
  const holds = async (): Promise<boolean> => {
    last = await onViewer(page, target, (el) => (window as any).__bdd.rememberedRange(el), undefined);
    return !!last.before && !!last.now && sameRange(last.before, last.now, 0.5);
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`${target.phrase} does not show the remembered value range (remembered ${JSON.stringify(last.before)}, now ${JSON.stringify(last.now)})`);
  }
}

export function rangeChange(page: Page, target: ElementRef): Promise<RangeChange> {
  return onViewer(page, target, (el) => (window as any).__bdd.rangeChange(el), undefined);
}

/** The value range (viewport) against the snapshot before the last change; `same` is read once the
 * viewer is quiet, so a reset that lands a tick later fails it rather than slipping past. */
export async function expectValueRange(page: Page, target: ElementRef, compare: 'narrower' | 'same' | 'wider'): Promise<void> {
  let last: RangeChange = {};
  const holds = async (): Promise<boolean | string> => {
    last = compare === 'same' ? await onViewer(page, target, (el) => (window as any).__bdd.quietRangeChange(el), undefined) : await rangeChange(page, target);
    if (!last.before || !last.now)
      return 'no range';
    if (compare === 'same')
      return sameRange(last.before, last.now, 1e-6);
    return compare === 'narrower' ? rangeArea(last.now) < rangeArea(last.before) * 0.95 : rangeArea(last.now) > rangeArea(last.before) * 1.05;
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`${target.phrase} does not show ${compare === 'same' ? 'the same' : `a ${compare}`} value range` +
      ` (before ${JSON.stringify(last.before)}, now ${JSON.stringify(last.now)})`);
  }
}

/** The value range lies within the column's values: no empty space beyond the data. */
export async function expectValueRangeWithin(page: Page, target: ElementRef, column: string): Promise<void> {
  const r = await onViewer(page, target, (el, c) => {
    const b = (window as any).__bdd;
    const v = b.viewerOf(el);
    const col = b.col(c, v.dataFrame);
    return {range: b.rangeChange(el).now as Range | undefined, min: col.stats.min as number, max: col.stats.max as number};
  }, column);
  expect(r.range, `${target.phrase} reports no value range`).toBeTruthy();
  expect(r.range!.top, `${target.phrase}'s range starts below "${column}" (${r.range!.top} < ${r.min})`).toBeGreaterThanOrEqual(r.min * 0.9);
  expect(r.range!.bottom, `${target.phrase}'s range ends above "${column}" (${r.range!.bottom} > ${r.max})`).toBeLessThanOrEqual(r.max * 1.1);
}

/** The range the color scale labels against the snapshot before the last change. */
export async function expectScaleRange(page: Page, target: ElementRef, compare: 'narrower' | 'wider'): Promise<void> {
  const span = (r: ScaleRange) => r.max - r.min;
  let last: ScaleChange = {};
  const holds = async (): Promise<boolean | string> => {
    last = await onViewer(page, target, (el) => (window as any).__bdd.scaleChange(el), undefined);
    if (!last.before || !last.now)
      return `no color scale ${!last.before ? 'at the snapshot' : 'now'}`;
    return compare === 'narrower' ? span(last.now) < span(last.before) : span(last.now) > span(last.before);
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the color scale of ${target.phrase} does not cover a ${compare} range (before ${JSON.stringify(last.before)}, now ${JSON.stringify(last.now)})`);
  }
}
