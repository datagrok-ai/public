/* Viewer runtime: what the `viewers` tier's steps are made of — an in-page helper set installed
   once per page (`window.__bdd`) and the Playwright-side readers over it.
   The principle: this is our platform, not a black box. A viewer that waits out a debounce, a menu
   that builds on a timer, a canvas that paints on the next frame — the platform can tell us when
   (its events) or stop waiting altogether (`immediateRendering`, armed here on every viewer the
   page will ever hold), so nothing in this module sleeps. When a signal or a name is missing, it
   is added to the core, never faked here. */
import {Locator, Page} from '@playwright/test';
import {expect, pollMs} from './patience.js';
import type {ElementRef} from './args.js';
import {reasonOf} from './failure.js';
import {typeVerified} from './gestures.js';
import {exactText, locate} from './locate.js';

declare const grok: any;
declare const DG: any;

export interface Box {
  x: number;
  y: number;
  width: number;
  height: number;
}

/** A viewer's canvas after a change: `delta` is the number of pixels that differ from the snapshot
 * taken before the change (the histogram distance when the bitmap changed size), `ink` the painted
 * pixels now, `inkBefore` those of the snapshot, `hue` and `hueBefore` the pixels in the selection
 * hue; `selected` rows and the `view` area in device pixels size the highlight a selection must
 * paint. */
export interface CanvasChange {
  delta: number;
  ink: number;
  inkBefore: number;
  hue: number;
  hueBefore: number;
  selected: number;
  viewPx: number;
  dpr: number;
}

/** A named hit area's painted pixels now and at the snapshot before the last change, and the
 * bitmap rectangle read (for the report when a reading is off). */
export interface AreaChange {
  ink: number;
  inkBefore: number;
  rect: {x: number; y: number; w: number; h: number};
}

/** A color drawn in an area and how many pixels of it. */
export interface AreaColor {
  hex: string;
  count: number;
}

export interface ScaleRange {
  min: number;
  max: number;
}

/** The range a viewer's color scale labels now and at the snapshot before the last change. */
export interface ScaleChange {
  before?: ScaleRange;
  now?: ScaleRange;
}

/** A viewer's value range (its viewport) now and at the snapshot before the last change. */
export interface RangeChange {
  before?: Range;
  now?: Range;
}

export interface Range {
  top: number;
  bottom: number;
  height: number;
  left: number;
  right: number;
}

export interface Balloon {
  type: string;
  message: string;
}

const POPUP = '.d4-menu-popup';

/** The sentence of an in-page error, without Playwright's evaluate prefix and stack. */
function reasonText(e: unknown): string {
  return String(e instanceof Error ? e.message : e).replace(/^[\w.]+: /, '').split('\n')[0];
}
const SEP = /\s*[>|]\s*/;
const MENU_ITEM = 'xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " d4-menu-item ")][1]';
// a hidden tooltip keeps its last content, so every read matches the visible one only
const TOOLTIP_ROWS = '.d4-tooltip:visible table.d4-row-tooltip-table tr';
const TOOLTIP_COLUMNS = TOOLTIP_ROWS + ' td:first-child';

/** Everything the in-page side needs, on `window.__bdd`. Self-contained: it runs in the browser. */
function install(): void {
  const w = window as any;
  if (w.__bdd)
    return;
  const renders = new WeakMap<Element, {count: number; last: number; sub?: any}>();
  const snapshots = new WeakMap<Element, {colors: Map<number, number>; ink: number; hue: number; renders: number; range?: Range;
    areas: Record<string, number>; rects: Record<string, Box>; scale?: ScaleRange; values: Record<string, unknown>; bitmap: ImageData;
    legend?: {mode: string; slot: string; items: number; keys: string[]; selected: string[]; width: number; height: number}}>();
  const rememberedValues: Record<string, unknown> = {};
  // how long a repaint the viewer says is pending may take before that is a platform failure
  const PENDING_CAP = 10000;
  const sizes = new WeakMap<Element, {width: string; height: string}>();
  const forced = new WeakMap<HTMLElement, MutationObserver>();
  const listeners = new WeakMap<Element, Record<string, {count: number; sub: any}>>();
  const armed: Record<string, Promise<unknown>> = {};
  const balloons: Balloon[] = [];
  const remembered: Record<string, Range | undefined> = {};
  let layout: any;
  let tokens = 0;
  const norm = (s: unknown) => String(s ?? '').toLowerCase().replace(/[^a-z0-9]/g, '');
  // the platform's selected-rows orange, as the pixels of a marker or a box drawn in it
  const isHue = (r: number, g: number, b: number) => r >= 150 && g >= 100 && g <= 200 && b <= 110;
  const isBlank = (r: number, g: number, b: number, a: number) => a === 0 || (r >= 250 && g >= 250 && b >= 250);

  const viewers = (): any[] => {
    const all: any[] = [];
    for (const view of Array.from(grok.shell.tableViews ?? []) as any[]) {
      for (const v of Array.from(view.viewers ?? []) as any[])
        all.push(v);
    }
    return all;
  };
  const viewerOf = (el: Element): any => {
    const v = viewers().find((x) => x.root === el || x.root.contains(el) || el.contains(x.root));
    if (!v)
      throw new Error('the element is not a viewer of an open table view');
    return v;
  };
  const arm = (v: any): void => {
    if (!v || renders.has(v.root))
      return;
    const stamp: {count: number; last: number; sub?: any} = {count: 0, last: 0};
    renders.set(v.root, stamp);
    try {
      v.immediateRendering = true;
    } catch { /* a JS viewer without the flag */ }
    // a JS viewer announces its own renders on `onRendered` (the host's onViewerRendered never fires for it)
    try {
      stamp.sub = (v.onRendered ?? v.onViewerRendered).subscribe(() => {
        stamp.count++;
        stamp.last = Date.now();
      });
    } catch { /* a viewer without the event */ }
  };
  const stampAll = (): void => {
    for (const v of viewers())
      arm(v);
  };
  /** Whether the viewer says a change is on its way to its canvas; undefined for a viewer without
   * the signal (a JS viewer). */
  const pending = (v: any): boolean | undefined => {
    try {
      const p = v.isRenderPending;
      return typeof p === 'boolean' ? p : undefined;
    }
    catch {
      return undefined;
    }
  };
  /** Resolves with the renders since the call once the viewer has nothing on its way: the platform
   * says whether a refresh or a repaint is pending, so a property that paints nothing returns at
   * once and a repaint is waited for as long as it takes — through every render it announces (a
   * legend change re-anchors in up to four passes 100 ms apart; a settle that ended on the first
   * render left the later ones to land on the next step's baseline). A pending one that never
   * lands is a platform failure, reported after `PENDING_CAP`. A viewer without the signal falls
   * back to a render event within `capMs`. */
  const settle = (el: Element, capMs: number): Promise<number> => {
    const v = viewerOf(el);
    arm(v);
    const stamp = renders.get(v.root)!;
    const before = stamp.count;
    const t0 = Date.now();
    return new Promise((resolve, reject) => {
      const tick = () => {
        const p = pending(v);
        if (stamp.count > before && p !== true)
          return resolve(stamp.count - before);
        if (p === false)
          return resolve(0);
        if (p === undefined) {
          if (Date.now() - t0 >= capMs)
            return resolve(0);
          return setTimeout(tick, 10);
        }
        if (Date.now() - t0 >= PENDING_CAP)
          return reject(new Error(`${v.type}: a repaint has been pending for ${PENDING_CAP} ms`));
        setTimeout(tick, 0);
      };
      setTimeout(tick, 0);
    });
  };
  /** After the frame a change would land on, and until the viewer has nothing pending. */
  const quiet = async (v: any): Promise<void> => {
    await new Promise((resolve) => requestAnimationFrame(() => setTimeout(resolve, 0)));
    const t0 = Date.now();
    while (pending(v) === true) {
      if (Date.now() - t0 >= PENDING_CAP)
        throw new Error(`${v.type}: a repaint has been pending for ${PENDING_CAP} ms`);
      await new Promise((resolve) => setTimeout(resolve, 0));
    }
  };
  const partsOf = (v: any): Record<string, Element> => v.getWidgetStatus?.()?.parts ?? {};
  const canvasOf = (v: any): HTMLCanvasElement => {
    const part = partsOf(v).canvas as HTMLCanvasElement | undefined;
    const cv = part ?? v.root.querySelector('canvas[name="canvas"]') ?? v.root.querySelector('canvas');
    if (!cv)
      throw new Error(`${v.type} has no canvas`);
    return cv;
  };
  /** What a viewer's hit areas are relative to: its canvas, else the root it reports, else its
   * root (a form, a filter panel — DOM viewers with no canvas). */
  const anchorOf = (v: any): Element => {
    const parts = partsOf(v);
    return parts.canvas ?? parts.root ?? v.root;
  };
  // a WebGL canvas (the 3D scatter plot) has no 2D context and no pixels to read: its picture is
  // the viewer's own `scene signature` reading
  const pixels = (cv: HTMLCanvasElement): ImageData => {
    const ctx = cv.getContext('2d');
    return ctx ? ctx.getImageData(0, 0, cv.width, cv.height) : new ImageData(1, 1);
  };
  /** The viewer's picture: its canvas with its `overlay` part composited on top when there is one
   * of the same size (the scatter plot draws regression lines, labels and stats on the overlay,
   * the grid its selection and current cell) — so paint on either layer counts. */
  const pixelsOf = (v: any): ImageData => {
    const cv = canvasOf(v);
    const img = pixels(cv);
    const over = partsOf(v).overlay as HTMLCanvasElement | undefined;
    if (!(over instanceof HTMLCanvasElement) || over === cv || over.width !== cv.width || over.height !== cv.height)
      return img;
    const top = pixels(over);
    const a = img.data;
    const b = top.data;
    for (let i = 0; i < b.length; i += 4) {
      const alpha = b[i + 3];
      if (alpha === 0)
        continue;
      const w = alpha / 255;
      a[i] = Math.round(b[i] * w + a[i] * (1 - w));
      a[i + 1] = Math.round(b[i + 1] * w + a[i + 1] * (1 - w));
      a[i + 2] = Math.round(b[i + 2] * w + a[i + 2] * (1 - w));
      a[i + 3] = Math.max(a[i + 3], alpha);
    }
    return img;
  };
  const histogram = (img: ImageData) => {
    const data = img.data;
    const colors = new Map<number, number>();
    let ink = 0;
    let hue = 0;
    for (let i = 0; i < data.length; i += 4) {
      const key = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
      colors.set(key, (colors.get(key) ?? 0) + 1);
      if (isBlank(data[i], data[i + 1], data[i + 2], data[i + 3]))
        continue;
      ink++;
      if (isHue(data[i], data[i + 1], data[i + 2]))
        hue++;
    }
    return {colors, ink, hue};
  };
  const areasOf = (v: any): Record<string, Box> => v.getWidgetStatus()?.hitAreas ?? {};
  const areaKey = (v: any, name: string): string => {
    const areas = areasOf(v);
    const key = Object.keys(areas).find((k) => norm(k) === norm(name));
    if (!key)
      throw new Error(`${v.type} has no "${name}" area right now; it has: ${Object.keys(areas).join(', ') || 'none'}`);
    return key;
  };
  /** A hit area (canvas coordinates) as device pixels of the canvas bitmap. */
  const deviceRect = (cv: HTMLCanvasElement, r: Box): {x: number; y: number; w: number; h: number} => {
    const scale = cv.width / cv.getBoundingClientRect().width;
    return {x: Math.floor(r.x * scale), y: Math.floor(r.y * scale), w: Math.max(1, Math.floor(r.width * scale)), h: Math.max(1, Math.floor(r.height * scale))};
  };
  /** The colors drawn inside a rectangle of the bitmap, by pixel count, blanks aside. */
  const colorsIn = (img: ImageData, r: {x: number; y: number; w: number; h: number}): Map<number, number> => {
    const colors = new Map<number, number>();
    const data = img.data;
    for (let y = r.y; y < Math.min(r.y + r.h, img.height); y++) {
      for (let x = r.x; x < Math.min(r.x + r.w, img.width); x++) {
        const i = (y * img.width + x) * 4;
        if (isBlank(data[i], data[i + 1], data[i + 2], data[i + 3]))
          continue;
        const key = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
        colors.set(key, (colors.get(key) ?? 0) + 1);
      }
    }
    return colors;
  };
  const inkIn = (img: ImageData, r: {x: number; y: number; w: number; h: number}): number => {
    let n = 0;
    for (const count of colorsIn(img, r).values())
      n += count;
    return n;
  };
  /** The colors drawn in at least `minPx` pixels, blanks and near-whites aside. */
  const palette = (el: Element, minPx: number): number => {
    const data = histogram(pixelsOf(viewerOf(el)));
    let n = 0;
    for (const [c, count] of data.colors) {
      const r = (c >> 16) & 255;
      const g = (c >> 8) & 255;
      const b = c & 255;
      if (count >= minPx && !(r >= 250 && g >= 250 && b >= 250) && c !== 0)
        n++;
    }
    return n;
  };
  const hex = (c: number): string => '#' + c.toString(16).padStart(6, '0').toUpperCase();
  /** The colors drawn inside a hit area, most pixels first, with the bitmap rectangle read. */
  const areaColors = (el: Element, name: string): {colors: AreaColor[]; rect: {x: number; y: number; w: number; h: number}; bitmap: number[]} => {
    const v = viewerOf(el);
    const cv = canvasOf(v);
    const rect = deviceRect(cv, areasOf(v)[areaKey(v, name)]);
    const colors = colorsIn(pixelsOf(v), rect);
    return {colors: [...colors.entries()].sort((a, b) => b[1] - a[1]).map(([c, count]) => ({hex: hex(c), count})), rect, bitmap: [cv.width, cv.height]};
  };
  /** The legend a viewer hosts, as the attributes it publishes on every commit: its mode and slot,
   * the item total, the rendered items' keys and the selected ones; undefined without a legend.
   * A legend shown in the tooltip is re-parented out of the viewer, so it is looked up there. */
  const legendState = (v: any): {mode: string; slot: string; items: number; keys: string[]; selected: string[]; width: number; height: number} | undefined => {
    const l: HTMLElement | null = v.root.querySelector('[name="legend"]') ?? document.querySelector('.d4-tooltip [name="legend"]');
    if (!l)
      return undefined;
    const items = Array.from(l.querySelectorAll('[name="legend-item"]'));
    const key = (i: Element) => i.getAttribute('data-item-key') ?? i.getAttribute('aria-label') ?? '';
    return {mode: l.dataset.legendMode ?? '', slot: l.dataset.legendSlot ?? '', items: Number(l.dataset.legendItems ?? items.length),
      keys: items.map(key), selected: items.filter((i) => i.getAttribute('aria-selected') === 'true' || i.getAttribute('data-item-selected') === 'true').map(key),
      width: l.offsetWidth, height: l.offsetHeight};
  };
  const legendChange = (el: Element): {before?: ReturnType<typeof legendState>; now?: ReturnType<typeof legendState>} => {
    const v = viewerOf(el);
    return {before: snapshots.get(v.root)?.legend, now: legendState(v)};
  };
  /** A reading kept by viewer type and name, for "as remembered" across a close and a reopen. */
  const rememberValue = (el: Element, name: string): void => {
    const v = viewerOf(el);
    const r = valueChange(el, name);
    if (r.now === undefined || r.now === null)
      throw new Error(`${v.type} has no "${name}" reading; it reports: ${r.has.join(', ') || 'no readings'}`);
    rememberedValues[`${v.type}|${norm(name)}`] = r.now;
  };
  const rememberedValue = (el: Element, name: string): {before?: unknown; now?: unknown; has: string[]} => {
    const v = viewerOf(el);
    const r = valueChange(el, name);
    return {before: rememberedValues[`${v.type}|${norm(name)}`], now: r.now, has: r.has};
  };
  /** A hit area's rectangle now and at the snapshot before the last change (anchor coordinates). */
  const areaRectChange = (el: Element, name: string): {before?: Box; now?: Box; has: string[]} => {
    const v = viewerOf(el);
    const areas = areasOf(v);
    const key = Object.keys(areas).find((k) => norm(k) === norm(name));
    const before = snapshots.get(v.root)?.rects ?? {};
    return {before: before[norm(name)], now: key === undefined ? undefined : areas[key], has: Object.keys(areas)};
  };
  /** The viewer's named readings (`getWidgetStatus().values`). */
  const valuesOf = (v: any): Record<string, unknown> => v.getWidgetStatus?.()?.values ?? {};
  const scaleRange = (v: any): ScaleRange | undefined => {
    const values = valuesOf(v);
    const min = values['color scale min'];
    const max = values['color scale max'];
    return typeof min === 'number' && typeof max === 'number' ? {min, max} : undefined;
  };
  /** A named reading now and at the snapshot before the last change. */
  const valueChange = (el: Element, name: string): {before?: unknown; now?: unknown; has: string[]} => {
    const v = viewerOf(el);
    const now = valuesOf(v);
    const key = Object.keys(now).find((k) => norm(k) === norm(name));
    const before = snapshots.get(v.root)?.values ?? {};
    const keyBefore = Object.keys(before).find((k) => norm(k) === norm(name));
    return {before: keyBefore === undefined ? undefined : before[keyBefore], now: key === undefined ? undefined : now[key], has: Object.keys(now)};
  };
  const tableOf = (el: Element): string => String(viewerOf(el).dataFrame?.name ?? '');
  const rangeOf = (v: any): Range | undefined => {
    try {
      const vp = v.viewport;
      return vp ? {top: vp.top, bottom: vp.bottom, height: vp.height, left: vp.left, right: vp.right} : undefined;
    }
    catch {
      return undefined;
    }
  };
  const property = (v: any, caption: string): any => {
    const want = norm(caption);
    const short = norm(caption.replace(/\s+columns?$/i, ''));
    const friendly = (name: string) => name.replace(/ColumnNames?$/, '');
    const props: any[] = v.getProperties();
    // "Marker Size" is the markerSize property; "Marker Size Column" is markerSizeColumnName
    const p = props.find((x) => norm(x.caption) === want) ?? props.find((x) => norm(x.name) === want) ??
      props.find((x) => norm(friendly(x.name)) === want) ??
      props.find((x) => /ColumnNames?$/.test(x.name) && norm(friendly(x.name)) === short);
    if (!p) {
      // a property's caption defaults to its name; the grid shows the name split into words
      const title = (x: any): string => x.caption && x.caption !== x.name ? x.caption :
        x.name.replace(/ColumnNames?$/, ' Column').replace(/[a-z\d](?=[A-Z])|[a-z](?=\d)|[A-Z](?=[A-Z][a-z])/g, '$& ')
          .replace(/^./, (c: string) => c.toUpperCase());
      const distance = (a: string, b: string): number => {
        const row = [...Array(b.length + 1).keys()];
        for (let i = 1; i <= a.length; i++) {
          let prev = row[0]++;
          for (let j = 1; j <= b.length; j++) {
            const cur = row[j];
            row[j] = Math.min(row[j] + 1, row[j - 1] + 1, prev + (a[i - 1] === b[j - 1] ? 0 : 1));
            prev = cur;
          }
        }
        return row[b.length];
      };
      const titles = props.map(title);
      const nearest = titles.map((t) => ({t, d: distance(norm(t), want)})).sort((x, y) => x.d - y.d).slice(0, 3).map((x) => x.t);
      throw new Error(`${v.type} has no "${caption}" property; nearest: ${nearest.join(', ')}. All: ${titles.join(', ')}`);
    }
    return p;
  };
  const HEX = /^#?([0-9a-f]{6})$/i;
  const convert = (p: any, text: string): unknown => {
    const s = text.replace(/\\n/g, '\n');
    const type = String(p.propertyType ?? '');
    if (type === 'bool')
      return /^(true|yes|on|checked|1)$/i.test(s);
    if (type === 'int' || type === 'double' || type === 'num' || type === 'bigint') {
      if (s === '')
        return null;
      const hex = HEX.exec(s);
      if (hex && type === 'int')
        return (0xFF000000 | parseInt(hex[1], 16)) >>> 0;
      return Number(s);
    }
    if (type === 'string_list' || type === 'list' || type === 'column_list')
      return s === '' ? [] : s.split(/\s*,\s*/);
    return s;
  };
  /** The property as text: "" for none, lists comma-joined, an int read against a `#rrggbb`
   * expectation as `#rrggbb` (colors are ints in the property bag). */
  const readProperty = (el: Element, caption: string, expected = ''): string => {
    const v = viewerOf(el);
    const p = property(v, caption);
    const value = v.props[p.name];
    if (value == null)
      return '';
    if (Array.isArray(value))
      return value.join(', ');
    if (typeof value === 'number' && String(p.propertyType) === 'int' && HEX.test(expected))
      return '#' + (value & 0xFFFFFF).toString(16).padStart(6, '0').toUpperCase();
    return String(value);
  };
  /** The canvas now — its histogram, the ink of every hit area, the value range, the color scale's
   * range — as the baseline the "than before" checks compare with. */
  const snapshot = (el: Element): number => {
    const v = viewerOf(el);
    arm(v);
    const hit = areasOf(v);
    const rects: Record<string, Box> = {};
    for (const key of Object.keys(hit))
      rects[norm(key)] = {...hit[key]};
    const legend = legendState(v);
    let cv: HTMLCanvasElement | undefined;
    try {
      cv = canvasOf(v);
    } catch { /* a DOM viewer: its areas, readings and legend are the snapshot */ }
    if (!cv) {
      snapshots.set(v.root, {colors: new Map(), ink: 0, hue: 0, renders: renders.get(v.root)!.count, range: rangeOf(v), areas: {}, rects,
        scale: scaleRange(v), values: {...valuesOf(v)}, bitmap: new ImageData(1, 1), legend});
      return 0;
    }
    const img = pixelsOf(v);
    const shot = histogram(img);
    const areas: Record<string, number> = {};
    for (const key of Object.keys(hit))
      areas[norm(key)] = inkIn(img, deviceRect(cv, hit[key]));
    snapshots.set(v.root, {...shot, renders: renders.get(v.root)!.count, range: rangeOf(v), areas, rects, scale: scaleRange(v), values: {...valuesOf(v)}, bitmap: img, legend});
    return shot.ink;
  };
  /** Pixels that differ between two bitmaps of the same size — a reorder of equal bars keeps the
   * color histogram and is still a repaint; bitmaps of different sizes differ by their histogram
   * and their size. */
  const pixelDelta = (before: ImageData, now: ImageData, colorsBefore: Map<number, number>, colorsNow: Map<number, number>): number => {
    if (before.width !== now.width || before.height !== now.height) {
      let delta = Math.abs(before.width * before.height - now.width * now.height);
      for (const [c, n] of colorsNow)
        delta += Math.abs(n - (colorsBefore.get(c) ?? 0));
      for (const [c, n] of colorsBefore) {
        if (!colorsNow.has(c))
          delta += n;
      }
      return delta;
    }
    const a = before.data;
    const b = now.data;
    let delta = 0;
    for (let i = 0; i < a.length; i += 4) {
      if (a[i] !== b[i] || a[i + 1] !== b[i + 1] || a[i + 2] !== b[i + 2] || a[i + 3] !== b[i + 3])
        delta++;
    }
    return delta;
  };
  /** Where the picture differs from the snapshot, for the failure to say: the changed pixels'
   * box in CSS px and the hit areas it touches, or the two canvas sizes when they differ. */
  const changedWhere = (v: any, cv: HTMLCanvasElement, before: ImageData, now: ImageData): string => {
    if (before.width !== now.width || before.height !== now.height)
      return `canvas ${before.width}x${before.height} → ${now.width}x${now.height}`;
    const a = before.data;
    const b = now.data;
    let x0 = now.width, y0 = now.height, x1 = -1, y1 = -1;
    for (let i = 0, p = 0; i < a.length; i += 4, p++) {
      if (a[i] !== b[i] || a[i + 1] !== b[i + 1] || a[i + 2] !== b[i + 2] || a[i + 3] !== b[i + 3]) {
        const x = p % now.width;
        const y = (p - x) / now.width;
        if (x < x0) x0 = x;
        if (x > x1) x1 = x;
        if (y < y0) y0 = y;
        if (y > y1) y1 = y;
      }
    }
    if (x1 < 0)
      return '';
    const s = cv.width / cv.getBoundingClientRect().width;
    const box = {x: x0 / s, y: y0 / s, width: (x1 - x0 + 1) / s, height: (y1 - y0 + 1) / s};
    const touched = Object.entries(areasOf(v)).filter(([n, r]) => n !== 'view' && r.x < box.x + box.width && box.x < r.x + r.width && r.y < box.y + box.height && box.y < r.y + r.height).map(([n]) => n);
    return `within ${Math.round(box.x)},${Math.round(box.y)} ${Math.round(box.width)}x${Math.round(box.height)} of the canvas`
      + (touched.length ? `, over: ${touched.join(', ')}` : '');
  };
  /** The painted pixels now, without moving the snapshot. */
  const ink = (el: Element): number => histogram(pixelsOf(viewerOf(el))).ink;
  // the baseline "should have repainted" compares with: the canvas before the change
  const baseline = (el: Element): void => {
    try {
      snapshot(el);
    } catch { /* not a viewer of an open view */ }
  };
  /** Every viewer's baseline at once — before a change that reaches them all (a filter, a
   * selection, a column's colors). */
  const baselineAll = (): void => {
    for (const v of viewers())
      baseline(v.root);
  };
  const writeProperties = async (el: Element, entries: [string, string][], capMs: number): Promise<number> => {
    const v = viewerOf(el);
    arm(v);
    // the state a "than before" claim compares with is the one the viewer had FINISHED, not
    // whatever it happened to be showing mid-render: a word cloud between two layouts reports no
    // words at all, and the claim after the write then has nothing to compare with
    await quiet(v).catch(() => undefined);
    baseline(el);
    const settled = settle(el, capMs);
    for (const [caption, text] of entries) {
      const p = property(v, caption);
      v.props[p.name] = convert(p, text);
    }
    return settled;
  };
  const canvasBox = (v: any): Box => {
    const r = anchorOf(v).getBoundingClientRect();
    return {x: r.x, y: r.y, width: r.width, height: r.height};
  };
  /** The named hit area in client coordinates, or the names the viewer reports instead. */
  const findArea = (el: Element, name: string, beforeChange = false): {box?: Box; has: string[]} => {
    const v = viewerOf(el);
    if (beforeChange)
      baseline(el);
    const areas: Record<string, Box> = v.getWidgetStatus()?.hitAreas ?? {};
    const has = Object.keys(areas);
    const key = has.find((k) => norm(k) === norm(name));
    if (!key)
      return {has};
    const r = areas[key];
    const cv = canvasBox(v);
    return {box: {x: cv.x + r.x, y: cv.y + r.y, width: r.width, height: r.height}, has};
  };
  const hitArea = (el: Element, name: string, beforeChange = false): Box => {
    const found = findArea(el, name, beforeChange);
    if (!found.box)
      throw new Error(`${viewerOf(el).type} has no "${name}" area right now; it has: ${found.has.join(', ') || 'none'}`);
    return found.box;
  };
  /** Every hit area the viewer reports, in client coordinates. */
  const areas = (el: Element): Record<string, Box> => {
    const v = viewerOf(el);
    const cv = canvasBox(v);
    const result: Record<string, Box> = {};
    for (const [key, r] of Object.entries(areasOf(v)))
      result[key] = {x: cv.x + r.x, y: cv.y + r.y, width: r.width, height: r.height};
    return result;
  };
  /** Painted pixels inside a hit area (the canvas may be scaled to the device). */
  const areaInk = (el: Element, name: string): number => {
    const v = viewerOf(el);
    const cv = canvasOf(v);
    return inkIn(pixelsOf(v), deviceRect(cv, areasOf(v)[areaKey(v, name)]));
  };
  /** A hit area's ink now against the snapshot's (the area must have been reported then too). */
  const areaChange = (el: Element, name: string): AreaChange => {
    const v = viewerOf(el);
    const before = snapshots.get(v.root);
    if (!before)
      throw new Error(`${v.type}: no snapshot to compare with`);
    const inkBefore = before.areas[norm(name)];
    if (inkBefore === undefined)
      throw new Error(`${v.type} reported no "${name}" area at the snapshot; it had: ${Object.keys(before.areas).join(', ') || 'none'}`);
    const cv = canvasOf(v);
    return {ink: areaInk(el, name), inkBefore, rect: deviceRect(cv, areasOf(v)[areaKey(v, name)])};
  };
  /** The pixels that changed inside one hit area since the snapshot — the area's own repaint,
   * as opposed to the whole canvas'. */
  const areaDelta = (el: Element, name: string): number => {
    const v = viewerOf(el);
    const before = snapshots.get(v.root);
    if (!before)
      throw new Error(`${v.type}: no snapshot to compare with`);
    const cv = canvasOf(v);
    const now = pixelsOf(v);
    if (before.bitmap.width !== now.width || before.bitmap.height !== now.height)
      throw new Error(`${v.type}: the canvas resized since the snapshot`);
    const r = deviceRect(cv, areasOf(v)[areaKey(v, name)]);
    const a = before.bitmap.data;
    const b = now.data;
    let delta = 0;
    for (let y = r.y; y < Math.min(r.y + r.h, now.height); y++) {
      for (let x = r.x; x < Math.min(r.x + r.w, now.width); x++) {
        const i = (y * now.width + x) * 4;
        if (a[i] !== b[i] || a[i + 1] !== b[i + 1] || a[i + 2] !== b[i + 2] || a[i + 3] !== b[i + 3])
          delta++;
      }
    }
    return delta;
  };
  const change = (el: Element): CanvasChange => {
    const v = viewerOf(el);
    const cv = canvasOf(v);
    const img = pixelsOf(v);
    const now = histogram(img);
    const before = snapshots.get(v.root);
    if (!before)
      throw new Error(`${v.type}: no snapshot to compare with`);
    const delta = pixelDelta(before.bitmap, img, before.colors, now.colors);
    const view = areasOf(v)['view'];
    const viewPx = view ? deviceRect(cv, view).w * deviceRect(cv, view).h : cv.width * cv.height;
    let selected = 0;
    try {
      selected = v.dataFrame.selection.trueCount;
    } catch { /* a viewer without a table */ }
    return {delta, ink: now.ink, inkBefore: before.ink, hue: now.hue, hueBefore: before.hue, selected, viewPx, dpr: window.devicePixelRatio};
  };
  const rangeChange = (el: Element): RangeChange => {
    const v = viewerOf(el);
    return {before: snapshots.get(v.root)?.range, now: rangeOf(v)};
  };
  const scaleChange = (el: Element): ScaleChange => {
    const v = viewerOf(el);
    return {before: snapshots.get(v.root)?.scale, now: scaleRange(v)};
  };
  /** A range kept by the viewer's type, so it survives the viewer being closed and reopened (a
   * project round-trip). */
  const rememberRange = (el: Element): void => {
    const v = viewerOf(el);
    remembered[String(v.type)] = rangeOf(v);
  };
  const rememberedRange = (el: Element): RangeChange => {
    const v = viewerOf(el);
    return {before: remembered[String(v.type)], now: rangeOf(v)};
  };
  /** Whether the viewer painted since its snapshot, read once it is quiet: after the frame a
   * change would land on, and after whatever it says is still pending. */
  const stillness = async (el: Element): Promise<{renders: number; delta: number; where: string}> => {
    const v = viewerOf(el);
    const before = snapshots.get(v.root);
    if (!before)
      throw new Error(`${v.type}: no snapshot to compare with`);
    await quiet(v);
    const cv = canvasOf(v);
    const delta = change(el).delta;
    return {renders: renders.get(v.root)!.count - before.renders, delta, where: delta ? changedWhere(v, cv, before.bitmap, pixelsOf(v)) : ''};
  };
  /** The value range against the snapshot's, read once the viewer is quiet (a reset that lands a
   * tick after the change is read, not missed). */
  const quietRangeChange = async (el: Element): Promise<RangeChange> => {
    await quiet(viewerOf(el));
    return rangeChange(el);
  };
  const quietValueChange = async (el: Element, name: string): Promise<{before?: unknown; now?: unknown; has: string[]}> => {
    await quiet(viewerOf(el));
    return valueChange(el, name);
  };
  /** One subscription per viewer and event, alive from "listens for" until `unlisten` (the
   * "should have fired" read, or the viewer closing). */
  const unlisten = (v: any, event?: string): void => {
    const all = listeners.get(v.root);
    if (!all)
      return;
    for (const e of event === undefined ? Object.keys(all) : [event]) {
      all[e]?.sub.unsubscribe();
      delete all[e];
    }
  };
  const listen = (el: Element, event: string): void => {
    const v = viewerOf(el);
    unlisten(v, event);
    const all = listeners.get(v.root) ?? {};
    listeners.set(v.root, all);
    const entry = {count: 0, sub: undefined as any};
    entry.sub = v.onEvent(event).subscribe(() => { entry.count++; });
    all[event] = entry;
  };
  const firedCount = (el: Element, event: string): number => listeners.get(viewerOf(el).root)?.[event]?.count ?? -1;
  const forget = (v: any): void => {
    unlisten(v);
    renders.get(v.root)?.sub?.unsubscribe();
    renders.delete(v.root);
  };
  /** Resolves once the element's box has been the same for two frames, or the cap runs out. The
   * dock manager sizes the element it hosts, so a size written while it is still laying a freshly
   * docked viewer out is overwritten by the pass that follows — and the viewer then reads at
   * whatever width the dock gave it, which is a fixture silently gone. */
  /** Resolves once [read] says the same thing on [frames] consecutive frames after the first — what
   * it describes has stopped moving — or after [capMs]. */
  const stable = (read: () => string, capMs: number, frames: number): Promise<void> => new Promise((done) => {
    const t0 = Date.now();
    let last = '';
    let same = 0;
    const tick = (): void => {
      const now = read();
      same = now === last ? same + 1 : 0;
      last = now;
      if (same >= frames || Date.now() - t0 >= capMs)
        return done();
      requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
  });
  const stableBox = (root: HTMLElement, capMs: number): Promise<void> => stable(() => {
    const r = root.getBoundingClientRect();
    return `${Math.round(r.width)}x${Math.round(r.height)}`;
  }, capMs, 2);
  /** A gesture aims at where the viewer has finished putting the thing, and a finished render is not
   * the end of that: a layout pass on the next frame can still move the whole viewer (a title just
   * set moved the pivot's grid down a row, and the right-click meant for a header landed on a cell).
   * The areas are relative to the anchor, so it is the anchor's box that must agree on two
   * consecutive frames — a cheap read, against the whole status per frame. */
  const stableArea = (el: Element, _name: string, capMs: number): Promise<void> => stable(() => {
    const r = anchorOf(viewerOf(el)).getBoundingClientRect();
    return [r.x, r.y, r.width, r.height].map(Math.round).join(',');
  }, capMs, 1);
  // the repaint a resize causes lands on the next task, so the settle is armed before the event
  const resize = async (el: Element, width: number | null, height: number | null, capMs: number): Promise<number> => {
    const root = viewerOf(el).root as HTMLElement;
    if (!sizes.has(root))
      sizes.set(root, {width: root.style.width, height: root.style.height});
    await stableBox(root, 1000);
    await quiet(viewerOf(el)).catch(() => undefined);
    baseline(el);
    const settled = settle(el, capMs);
    // and the dock keeps sizing it afterwards, whenever anything else in the view is laid out: the
    // size a step asked for is the fixture every claim after it is made against, so it is held
    // until the step that gives it back
    forced.get(root)?.disconnect();
    const want = {
      width: width === null ? null : `${width}px`,
      height: height === null ? null : `${height}px`,
    };
    const apply = (): void => {
      if (want.width !== null && root.style.width !== want.width)
        root.style.width = want.width;
      if (want.height !== null && root.style.height !== want.height)
        root.style.height = want.height;
    };
    const observer = new MutationObserver(apply);
    observer.observe(root, {attributes: true, attributeFilter: ['style']});
    forced.set(root, observer);
    apply();
    window.dispatchEvent(new Event('resize'));
    return settled;
  };
  const restoreSize = (el: Element, capMs: number): Promise<number> => {
    const root = viewerOf(el).root as HTMLElement;
    const size = sizes.get(root);
    if (!size)
      return Promise.resolve(0);
    forced.get(root)?.disconnect();
    forced.delete(root);
    baseline(el);
    const settled = settle(el, capMs);
    root.style.width = size.width;
    root.style.height = size.height;
    sizes.delete(root);
    window.dispatchEvent(new Event('resize'));
    return settled;
  };
  /** Subscribes to a `grok.events` stream before a gesture; `waitArmed` collects the outcome after. */
  const armEvent = (name: string, capMs: number): string => {
    const token = `t${++tokens}`;
    armed[token] = new Promise((resolve) => {
      let sub: any;
      try {
        sub = grok.events[name].subscribe((args: unknown) => {
          sub.unsubscribe();
          resolve(args ?? true);
        });
      } catch {
        resolve(undefined);
        return;
      }
      setTimeout(() => {
        try {
          sub?.unsubscribe();
        } catch { /* already gone */ }
        resolve(undefined);
      }, capMs);
    });
    return token;
  };
  const waitArmed = async (token: string): Promise<boolean> => {
    const result = await armed[token];
    delete armed[token];
    return result !== undefined;
  };
  /** Closes the open popup menu the platform's way (a click outside) and waits for it to be gone. */
  const closeMenu = async (): Promise<void> => {
    if (!document.querySelector('.d4-menu-popup'))
      return;
    const closed = armEvent('onContextMenuClosed', 1500);
    document.body.click();
    await waitArmed(closed);
  };
  /** Ready for a right-click: no popup open, `onContextMenuShown` armed (`onContextMenu` fires
   * before the menu exists). */
  const openMenu = async (capMs: number): Promise<string> => {
    await closeMenu();
    return armEvent('onContextMenuShown', capMs);
  };
  /** Where to right-click an element for its context menu — a named hit area, else the viewer's
   * `view` area, else the element's centre — with the canvas baseline taken and the menu armed. */
  const menuPoint = async (el: Element, area: string | null, capMs: number): Promise<{x: number; y: number; token: string}> => {
    // the point the menu opens at must be where the viewer has finished putting the thing: a title
    // just set moves the grid under it, and the right-click lands a row off
    await settle(el, 300).catch(() => undefined);
    await stableArea(el, area ?? 'view', 1000);
    let box: Box | undefined;
    try {
      box = hitArea(el, area ?? 'view', true);
    } catch (e) {
      if (area !== null)
        throw e;
    }
    if (!box) {
      const r = el.getBoundingClientRect();
      box = {x: r.x, y: r.y, width: r.width, height: r.height};
    }
    return {x: box.x + box.width / 2, y: box.y + box.height / 2, token: await openMenu(capMs)};
  };
  const addViewer = (type: string): void => {
    const types: string[] = DG.Viewer.getViewerTypes();
    const exact = types.find((t) => norm(t) === norm(type));
    if (!exact)
      throw new Error(`no viewer type "${type}"; the platform has: ${types.join(', ')}`);
    arm(grok.shell.tv.addViewer(exact));
  };
  /** The balloons shown since the last read, and clears them. */
  const takeBalloons = (): Balloon[] => balloons.splice(0, balloons.length);
  const saveLayout = (): void => { layout = grok.shell.tv.saveLayout(); };
  /** Saves through the server and keeps only the id: "loads the saved layout" then fetches what
   * the server stored, so the round-trip covers the serialization too. */
  const saveLayoutToServer = async (): Promise<string> => {
    const l = grok.shell.tv.saveLayout();
    await grok.dapi.layouts.save(l);
    layout = {serverId: l.id};
    return l.id;
  };
  const loadLayout = async (): Promise<void> => {
    if (!layout)
      throw new Error('no layout saved in this feature');
    const saved = layout.serverId ? await grok.dapi.layouts.find(layout.serverId) : layout;
    if (!saved)
      throw new Error(`the server has no layout ${layout.serverId}`);
    grok.shell.tv.loadLayout(saved);
  };
  const deleteLayout = async (id: string): Promise<void> => {
    const saved = await grok.dapi.layouts.find(id);
    if (saved)
      await grok.dapi.layouts.delete(saved);
  };

  /** The function call a menu command starts, and the current table's columns before it: the
   * platform announces every call (`onBeforeRunAction` / `onAfterRunAction`), and the one whose
   * function is registered under the picked menu path is the command — its end is when the
   * command is done, whatever dialog it showed in between. */
  let command: {name: string; done: Promise<void>; started: number} | undefined;
  let columnsBefore: {table: any; names: string[]} | undefined;
  let commandArm: any;
  const armCommand = (path: string): void => {
    const t = grok.shell.t;
    columnsBefore = t ? {table: t.dart, names: t.columns.names()} : undefined;
    command = undefined;
    commandArm?.unsubscribe();
    const want = norm(path);
    const sub = grok.functions.onBeforeRunAction.subscribe((fc: any) => {
      const menu = fc?.func?.topMenu;
      if (!menu || norm(menu) !== want)
        return;
      sub.unsubscribe();
      if (commandArm === sub)
        commandArm = undefined;
      let resolve!: () => void;
      const done = new Promise<void>((r) => { resolve = r; });
      const after = grok.functions.onAfterRunAction.subscribe((ended: any) => {
        if (ended?.dart === fc.dart || ended?.id === fc.id) {
          after.unsubscribe();
          resolve();
        }
      });
      command = {name: String(fc.func?.nqName ?? fc.func?.name ?? path), done, started: Date.now()};
    });
    commandArm = sub;
    // a pick that started no call within 5 s is disarmed — this arm only, never a later pick's
    setTimeout(() => {
      sub.unsubscribe();
      if (commandArm === sub)
        commandArm = undefined;
    }, 5000);
  };
  const waitCommand = async (capMs: number): Promise<string> => {
    // the click's handler may start the call a task or two later
    const t0 = Date.now();
    while (!command && commandArm && Date.now() - t0 < 5000)
      await new Promise((r) => setTimeout(r, 10));
    if (!command)
      throw new Error('no menu command has started a function call in this scenario (a Dart command, or the menu item ran nothing)');
    const c = command;
    const timeout = new Promise<'timeout'>((r) => setTimeout(() => r('timeout'), capMs));
    if (await Promise.race([c.done.then(() => 'done'), timeout]) === 'timeout')
      throw new Error(`${c.name} has been running for ${Math.round((Date.now() - c.started) / 1000)} s`);
    return c.name;
  };
  /** The current table's columns now and before the last menu command; `same` says whether the
   * table is still the one the command started on. */
  const columnsSince = (): {before: string[] | null; now: string[]; same: boolean} => {
    const t = grok.shell.t;
    return {before: columnsBefore?.names ?? null, now: t ? t.columns.names() : [], same: !!t && columnsBefore?.table === t.dart};
  };

  /** Custom platform events (`grok.events.fireCustomEvent`) by id, counted from "listens for"
   * until the page resets; a read takes the count and the last arguments and zeroes them. */
  const customEvents = new Map<string, {count: number; last: unknown; sub: any}>();
  const listenCustom = (id: string): void => {
    if (customEvents.has(id))
      return;
    const entry = {count: 0, last: undefined as unknown, sub: undefined as any};
    entry.sub = grok.events.onCustomEvent(id).subscribe((args: unknown) => { entry.count++; entry.last = args; });
    customEvents.set(id, entry);
  };
  const customFired = (id: string, take: boolean): {count: number; last: unknown} => {
    const entry = customEvents.get(id);
    if (!entry)
      return {count: -1, last: undefined};
    let last: unknown;
    try {
      last = JSON.parse(JSON.stringify(entry.last ?? null));
    }
    catch {
      last = String(entry.last);
    }
    const read = {count: entry.count, last};
    if (take) {
      entry.count = 0;
      entry.last = undefined;
    }
    return read;
  };

  w.__bdd = {viewerOf, arm, stampAll, settle, quiet, readProperty, writeProperties, findArea, hitArea, areas, areaInk, areaChange, areaDelta, areaColors,
    areaRectChange, legendState: (el: Element) => legendState(viewerOf(el)), legendChange, rememberValue, rememberedValue,
    snapshot, baselineAll, change, rangeChange, quietRangeChange, scaleChange, valueChange, quietValueChange, rememberRange, rememberedRange, stillness,
    palette, tableOf, listen, unlisten, firedCount, resize, restoreSize, armEvent, waitArmed, closeMenu, openMenu, menuPoint, stableArea, addViewer,
    takeBalloons, saveLayout, saveLayoutToServer, loadLayout, deleteLayout, ink, armCommand, waitCommand, columnsSince, listenCustom, customFired};
  stampAll();
  grok.events.onViewerAdded.subscribe((a: any) => arm(a?.args?.viewer));
  grok.events.onViewerClosed.subscribe((a: any) => a?.args?.viewer && forget(a.args.viewer));
  grok.events.onEvent('d4-balloon-shown').subscribe((a: any) => balloons.push({type: String(a?.args?.type ?? ''), message: String(a?.args?.message ?? '')}));
}

const installed = new WeakSet<Page>();
const watched = new WeakSet<Page>();

/** Installs the in-page side once the shell is up. Cheap on every step: the page is remembered,
 * and forgotten when its main frame navigates (a reload drops the in-page side). */
export async function installViewerRuntime(page: Page): Promise<void> {
  if (installed.has(page))
    return;
  await page.evaluate(install);
  installed.add(page);
  if (!watched.has(page)) {
    watched.add(page);
    page.on('framenavigated', (frame) => {
      if (frame === page.mainFrame())
        installed.delete(page);
    });
  }
}

/** The single visible element of a viewer phrase (a closed view can leave a zero-size twin). */
export async function viewerLocator(page: Page, target: ElementRef): Promise<Locator> {
  const loc = await locate(page, target);
  return loc.filter({visible: true}).first();
}

/** `hitArea` in client coordinates — the rectangle the viewer reports for that named region;
 * `beforeChange` takes the canvas baseline in the same call (a click or a double-click follows).
 * Waits for the area the way a locator waits for its element: a viewer that renders twice on a
 * change (a category switch relays out after its first paint) reports the area after the second. */
export async function hitArea(page: Page, target: ElementRef, name: string, beforeChange = false): Promise<Box> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  const find = (): Promise<{box?: Box; has: string[]}> =>
    loc.evaluate((el, [n, b]) => (window as any).__bdd.findArea(el, n, b), [name, beforeChange] as [string, boolean]);
  // a gesture aims at where the viewer has finished putting the thing: a network diagram still
  // running its physics moves the node between the read and the click. A viewer that says nothing
  // is pending answers at once; one that never stops is the claim's problem, not the gesture's
  if (beforeChange)
    await loc.evaluate((el) => (window as any).__bdd.settle(el, 300)).catch(() => undefined);
  let found = await find();
  if (!found.box) {
    await expect.poll(async () => (found = await find()).box !== undefined,
      {timeout: pollMs(5000), message: `${target.phrase} reports no "${name}" area; it has: ${found.has.join(', ') || 'none'}`}).toBe(true);
  }
  return found.box!;
}

export async function expectHasArea(page: Page, target: ElementRef, name: string, negate = false): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  const find = (): Promise<{box?: Box; has: string[]}> => loc.evaluate((el, n) => (window as any).__bdd.findArea(el, n, false), name);
  let found = await find();
  const poll = expect.poll(async () => (found = await find()).box !== undefined,
    {timeout: pollMs(5000), message: `${target.phrase} ${negate ? 'still reports' : 'reports no'} "${name}" area; it has: ${found.has.join(', ') || 'none'}`});
  await (negate ? poll.not : poll).toBe(true);
}

export function centerOf(box: Box): {x: number; y: number} {
  return {x: box.x + box.width / 2, y: box.y + box.height / 2};
}

/** Every hit area the viewer reports right now, in client coordinates — for a step that reasons
 * over them together (the order of the bars, a spot no bar covers). */
export async function hitAreas(page: Page, target: ElementRef): Promise<Record<string, Box>> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  return loc.evaluate((el) => (window as any).__bdd.areas(el));
}

export async function addViewer(page: Page, type: string): Promise<void> {
  await installViewerRuntime(page);
  await page.evaluate((t) => { (window as any).__bdd.addViewer(t); }, type);
  await page.locator(`[name="viewer-${type.replace(/\s+/g, '-')}" i]`).filter({visible: true}).first().waitFor();
}

/** Sets properties by caption in one go: one settle for the group (the sets coalesce into a
 * single repaint under immediate rendering; the viewer says whether one is pending, the cap only
 * matters for a viewer without that signal). The canvas is snapshotted first, so `should have
 * repainted` compares with the state before the change. One roundtrip: `locator.evaluate` waits
 * for the element itself. */
export async function setProperties(page: Page, target: ElementRef, entries: [string, string][], capMs = 300): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await loc.evaluate((el, [e, cap]) => (window as any).__bdd.writeProperties(el, e, cap), [entries, capMs] as [[string, string][], number]);
}

export async function readProperty(page: Page, target: ElementRef, caption: string, expected = ''): Promise<string> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  return loc.evaluate((el, [c, x]) => (window as any).__bdd.readProperty(el, c, x), [caption, expected] as [string, string]);
}

export async function expectProperty(page: Page, target: ElementRef, caption: string, value: string, negate = false): Promise<void> {
  const poll = expect.poll(() => readProperty(page, target, caption, value), {timeout: pollMs(5000), message: `"${caption}" of ${target.phrase}`});
  await (negate ? poll.not : poll).toBe(value.replace(/\\n/g, '\n'));
}

export async function snapshot(page: Page, target: ElementRef): Promise<number> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  return loc.evaluate((el) => (window as any).__bdd.snapshot(el));
}

/** Every viewer of every open table view gets its baseline — before a change that is not one
 * viewer's own (a filter, a selection, a column's colors). */
export async function baselineAll(page: Page): Promise<void> {
  await installViewerRuntime(page);
  await page.evaluate(() => { (window as any).__bdd.baselineAll(); });
}

/** Waits until no viewer of any open table view has a refresh or a repaint pending — after a
 * change that reaches them all, so the next step's baseline is the state after it. */
export async function settleAll(page: Page): Promise<void> {
  await installViewerRuntime(page);
  await page.evaluate(async () => {
    const b = (window as any).__bdd;
    for (const view of Array.from(grok.shell.tableViews ?? []) as any[]) {
      for (const v of Array.from(view.viewers ?? []) as any[])
        await b.quiet(v);
    }
  });
}

export async function canvasChange(page: Page, target: ElementRef): Promise<CanvasChange> {
  const loc = await viewerLocator(page, target);
  return loc.evaluate((el) => (window as any).__bdd.change(el));
}

/** Waits until the canvas differs from the last snapshot by at least `minPx` pixels, then makes
 * the new state the snapshot. */
export async function expectRepainted(page: Page, target: ElementRef, minPx = 1): Promise<void> {
  await expect.poll(async () => (await canvasChange(page, target)).delta,
    {timeout: pollMs(10000), message: minPx > 1 ? `${target.phrase} did not repaint by ${minPx} pixels` : `${target.phrase} did not repaint`}).toBeGreaterThanOrEqual(minPx);
}

/** A hit area's own repaint: the pixels inside its rectangle that differ from the snapshot. */
export async function expectAreaRepainted(page: Page, target: ElementRef, area: string, minPx = 1): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await expect.poll(() => loc.evaluate((el, a) => (window as any).__bdd.areaDelta(el, a), area),
    {timeout: pollMs(10000), message: `the "${area}" area of ${target.phrase} did not repaint`}).toBeGreaterThanOrEqual(minPx);
}

/** The table the viewer draws (`viewer.dataFrame`), not the property it was asked to bind. */
export async function expectBoundTable(page: Page, target: ElementRef, name: string): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await expect.poll(() => loc.evaluate((el) => (window as any).__bdd.tableOf(el)), {timeout: pollMs(5000), message: `the table ${target.phrase} is bound to`}).toBe(name);
}

/** A hit area's ink against the snapshot before the last change. */
export async function expectAreaInk(page: Page, target: ElementRef, area: string, compare: 'less' | 'more'): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  let last: AreaChange | undefined;
  let failure = '';
  const holds = async (): Promise<boolean> => {
    try {
      const c: AreaChange = last = await loc.evaluate((el, a) => (window as any).__bdd.areaChange(el, a), area);
      return compare === 'less' ? c.ink < c.inkBefore : c.ink > c.inkBefore;
    }
    catch (e) {
      failure = reasonText(e);
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

const COLOR_MIN_PX = 10;
const HUE_TOLERANCE = 20;
const GREY_SATURATION = 0.15;
const LIGHTNESS_TOLERANCE = 0.15;

function hsl(color: string): {h: number; s: number; l: number} {
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
function near(a: string, b: string): boolean {
  const x = hsl(a);
  const y = hsl(b);
  if (x.s < GREY_SATURATION || y.s < GREY_SATURATION)
    return x.s < GREY_SATURATION && y.s < GREY_SATURATION && Math.abs(x.l - y.l) <= LIGHTNESS_TOLERANCE;
  const dh = Math.abs(x.h - y.h);
  return Math.min(dh, 360 - dh) <= HUE_TOLERANCE;
}

interface AreaColors {
  colors: AreaColor[];
  rect: {x: number; y: number; w: number; h: number};
  bitmap: number[];
}

async function areaColors(page: Page, target: ElementRef, area: string): Promise<AreaColors> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  return loc.evaluate((el, a) => (window as any).__bdd.areaColors(el, a), area);
}

const describeColors = (read: AreaColors | undefined): string => read ?
  `its colors: ${read.colors.slice(0, 5).map((c) => `${c.hex} (${c.count})`).join(', ') || 'none'} in ${JSON.stringify(read.rect)} of a ${read.bitmap.join('x')} bitmap` : '';

/** The color (within a shade of anti-aliasing) covers some pixels of the area. */
export async function expectAreaColor(page: Page, target: ElementRef, area: string, color: string): Promise<void> {
  const want = '#' + color.replace(/^#/, '').toUpperCase();
  if (!/^#[0-9A-F]{6}$/.test(want))
    throw new Error(`"${color}" is not a #rrggbb color`);
  await hitArea(page, target, area);
  let last: AreaColors | undefined;
  const count = async (): Promise<number> => {
    const read: AreaColors = last = await areaColors(page, target, area);
    return read.colors.filter((c) => near(c.hex, want)).reduce((n, c) => n + c.count, 0);
  };
  try {
    await expect.poll(count, {timeout: pollMs(5000)}).toBeGreaterThanOrEqual(COLOR_MIN_PX);
  }
  catch {
    throw new Error(`the "${area}" area of ${target.phrase} is not painted in ${want}; ${describeColors(last)}`);
  }
}

const SIGNIFICANT_PX = 30;

/** The area is painted in at least `count` hues: colors covering some pixels each, grouped by
 * `near` (a hue and its anti-aliased shades are one), greys and white aside — a cell whose letters
 * take their colors from the data, not a text cell. */
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

/** Two areas are painted in different colors: one of them has a color (covering some pixels) the
 * other has nothing near. */
export async function expectAreasDiffer(page: Page, target: ElementRef, a: string, b: string): Promise<void> {
  await hitArea(page, target, a);
  await hitArea(page, target, b);
  let shownA: AreaColors | undefined;
  let shownB: AreaColors | undefined;
  const own = (mine: AreaColor[], theirs: AreaColor[]): boolean =>
    mine.some((c) => c.count >= SIGNIFICANT_PX && !theirs.some((d) => d.count >= COLOR_MIN_PX && near(c.hex, d.hex)));
  const differ = async (): Promise<boolean> => {
    const ra: AreaColors = shownA = await areaColors(page, target, a);
    const rb: AreaColors = shownB = await areaColors(page, target, b);
    return own(ra.colors, rb.colors) || own(rb.colors, ra.colors);
  };
  try {
    await expect.poll(differ, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the "${a}" and "${b}" areas of ${target.phrase} are painted in the same colors; "${a}" ${describeColors(shownA)}; "${b}" ${describeColors(shownB)}`);
  }
}

export type ReadingCompare = 'equal' | 'lower' | 'higher' | 'differ' | 'same';

/** A reading of the viewer as it is now; a name the viewer does not report fails naming the
 * readings it does. */
export async function readValue(page: Page, target: ElementRef, name: string): Promise<unknown> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  const r: {now?: unknown; has: string[]} = await loc.evaluate((el, n) => (window as any).__bdd.valueChange(el, n), name);
  if (r.now === undefined || r.now === null)
    throw new Error(`${target.phrase} has no "${name}" reading; it reports: ${r.has.join(', ') || 'no readings'}`);
  return r.now;
}

/** A reading of the viewer (`getWidgetStatus().values`: "rows shown", "bars", "scene signature")
 * equals a value, is lower/higher than at the snapshot before the last change, differs from it, or
 * is the same — the negative read once the viewer is quiet, and read once. */
export async function expectReading(page: Page, target: ElementRef, name: string, compare: ReadingCompare, value?: number): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  let last: {before?: unknown; now?: unknown; has: string[]} = {has: []};
  const holds = async (): Promise<boolean | string> => {
    last = await loc.evaluate((el, [n, q]) => {
      const b = (window as any).__bdd;
      return q ? b.quietValueChange(el, n) : b.valueChange(el, n);
    }, [name, compare === 'same'] as [string, boolean]);
    if (last.now === undefined || last.now === null)
      return `no "${name}" reading`;
    if (compare === 'equal')
      return last.now === value;
    if (last.before === undefined || last.before === null)
      return `no "${name}" reading at the snapshot`;
    if (compare === 'differ' || compare === 'same')
      return (last.now !== last.before) === (compare === 'differ');
    if (typeof last.now !== 'number' || typeof last.before !== 'number')
      return `"${name}" is not a number`;
    return compare === 'lower' ? last.now < last.before : last.now > last.before;
  };
  const what = {equal: `${value}`, lower: 'lower than before', higher: 'higher than before', differ: 'different from before', same: 'the same as before'}[compare];
  const report = (): never => {
    throw new Error(`"${name}" of ${target.phrase} is ${String(last.now)}, not ${what}` +
      (compare === 'equal' ? '' : ` (${String(last.before)})`) +
      (last.now === undefined || last.now === null ? `; the viewer reports: ${last.has.join(', ') || 'no readings'}` : ''));
  };
  if (compare === 'same') {
    if (await holds() !== true)
      report();
    return;
  }
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    report();
  }
}

/** The legend sits on the side the viewer's own slot class names (`d4-legend-left`, …). */
export async function expectLegendSide(page: Page, target: ElementRef, side: string): Promise<void> {
  const want = side.toLowerCase();
  const legend = (await viewerLocator(page, target)).locator('[name="legend"]').filter({visible: true}).first();
  await expect(legend, `the legend of ${target.phrase}`).toBeVisible();
  const classes = (await legend.getAttribute('class') ?? '').split(/\s+/);
  const sides = classes.filter((c) => /^d4-legend-(left|right|top|bottom)$/.test(c)).map((c) => c.slice('d4-legend-'.length));
  expect(sides, `the legend of ${target.phrase} is on the ${sides.join(', ') || 'unknown'} side, not the ${want}`).toContain(want);
}

/** The range the color scale labels against the snapshot before the last change. */
export async function expectScaleRange(page: Page, target: ElementRef, compare: 'narrower' | 'same' | 'wider'): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  const span = (r: ScaleRange) => r.max - r.min;
  let last: ScaleChange = {};
  const holds = async (): Promise<boolean | string> => {
    last = await loc.evaluate((el) => (window as any).__bdd.scaleChange(el));
    if (!last.before || !last.now)
      return `no color scale ${!last.before ? 'at the snapshot' : 'now'}`;
    if (compare === 'same')
      return Math.abs(span(last.before) - span(last.now)) < 1e-6;
    return compare === 'narrower' ? span(last.now) < span(last.before) : span(last.now) > span(last.before);
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the color scale of ${target.phrase} does not cover ${compare === 'same' ? 'the same' : `a ${compare}`} range` +
      ` (before ${JSON.stringify(last.before)}, now ${JSON.stringify(last.now)})`);
  }
}

/** The canvas is what it was at the snapshot, read after the frame a repaint would have landed
 * on (a viewer may run a render pass that draws the same picture — the mouse-over row does). */
export async function expectNotRepainted(page: Page, target: ElementRef): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  const still: {renders: number; delta: number; where: string} = await loc.evaluate((el) => (window as any).__bdd.stillness(el));
  expect(still.delta, `${target.phrase} repainted: ${still.delta} px changed in ${still.renders} render(s)${still.where ? ', ' + still.where : ''}`).toBe(0);
}

export async function expectInk(page: Page, target: ElementRef, compare: 'less' | 'more' | 'some'): Promise<void> {
  if (compare === 'some') {
    await installViewerRuntime(page);
    const loc = await viewerLocator(page, target);
    // a viewer between two layouts has no canvas to read, and a throw from the callback would end
    // the poll on the spot: the reason is kept and told at the end instead
    let why = '';
    try {
      await expect.poll(() => loc.evaluate((el) => (window as any).__bdd.ink(el)).catch((e: Error) => { why = reasonOf(e); return 0; }),
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

const HIGHLIGHT_FLOOR = 200;
const HIGHLIGHT_PER_ROW = 2;

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

export async function rememberRange(page: Page, target: ElementRef): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await loc.evaluate((el) => { (window as any).__bdd.rememberRange(el); });
}

/** Both axes of a viewport, when it has a horizontal one: a viewer that restores the Y window and
 * loses the X window is not showing the same range. */
function sameRange(a: Range, b: Range, eps: number): boolean {
  if (Math.abs(a.top - b.top) >= eps || Math.abs(a.bottom - b.bottom) >= eps)
    return false;
  const horizontal = [a.left, a.right, b.left, b.right].every((x) => typeof x === 'number' && isFinite(x));
  return !horizontal || (Math.abs(a.left! - b.left!) < eps && Math.abs(a.right! - b.right!) < eps);
}

function rangeArea(r: Range): number {
  const width = typeof r.left === 'number' && typeof r.right === 'number' && isFinite(r.left) && isFinite(r.right)
    ? Math.abs(r.right - r.left) : 1;
  return Math.abs(r.height) * (width || 1);
}

/** The value range equals the one remembered for this viewer type — across a close and a reopen. */
export async function expectRememberedRange(page: Page, target: ElementRef): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  let last: RangeChange = {};
  const holds = async (): Promise<boolean> => {
    last = await loc.evaluate((el) => (window as any).__bdd.rememberedRange(el));
    return !!last.before && !!last.now && sameRange(last.before, last.now, 0.5);
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`${target.phrase} does not show the remembered value range (remembered ${JSON.stringify(last.before)}, now ${JSON.stringify(last.now)})`);
  }
}

export async function expectPalette(page: Page, target: ElementRef, min: number): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await expect.poll(() => loc.evaluate((el) => (window as any).__bdd.palette(el, 500)),
    {timeout: pollMs(5000), message: `${target.phrase} is painted in fewer than ${min} colors`}).toBeGreaterThanOrEqual(min);
}

export async function expectAreaPainted(page: Page, target: ElementRef, area: string): Promise<void> {
  await hitArea(page, target, area);
  const loc = await viewerLocator(page, target);
  await expect.poll(() => loc.evaluate((el, a) => (window as any).__bdd.areaInk(el, a), area),
    {timeout: pollMs(5000), message: `the "${area}" area of ${target.phrase} is blank`}).toBeGreaterThan(0);
}

export async function rangeChange(page: Page, target: ElementRef): Promise<RangeChange> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  return loc.evaluate((el) => (window as any).__bdd.rangeChange(el));
}

/** The value range (viewport) against the snapshot before the last change; `same` is read once the
 * viewer is quiet, so a reset that lands a tick later fails it rather than slipping past. */
export async function expectValueRange(page: Page, target: ElementRef, compare: 'narrower' | 'same' | 'wider'): Promise<void> {
  const same = (a: Range, b: Range) => sameRange(a, b, 1e-6);
  const loc = await viewerLocator(page, target);
  let last: RangeChange = {};
  const holds = async (): Promise<boolean | string> => {
    last = compare === 'same' ? await loc.evaluate((el) => (window as any).__bdd.quietRangeChange(el)) : await rangeChange(page, target);
    if (!last.before || !last.now)
      return 'no range';
    if (compare === 'same')
      return same(last.before, last.now);
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
  const loc = await viewerLocator(page, target);
  const r = await loc.evaluate((el, c) => {
    const b = (window as any).__bdd;
    const v = b.viewerOf(el);
    const col = v.dataFrame.col(c);
    if (!col)
      throw new Error(`no "${c}" column in ${v.dataFrame.name}`);
    return {range: b.rangeChange(el).now, min: col.stats.min, max: col.stats.max};
  }, column);
  expect(r.range, `${target.phrase} reports no value range`).toBeTruthy();
  expect(r.range!.top, `${target.phrase}'s range starts below "${column}" (${r.range!.top} < ${r.min})`).toBeGreaterThanOrEqual(r.min * 0.9);
  expect(r.range!.bottom, `${target.phrase}'s range ends above "${column}" (${r.range!.bottom} > ${r.max})`).toBeLessThanOrEqual(r.max * 1.1);
}

export async function listenFor(page: Page, target: ElementRef, event: string): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await loc.evaluate((el, e) => (window as any).__bdd.listen(el, e), event);
}

export async function expectFired(page: Page, target: ElementRef, event: string): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await expect.poll(() => loc.evaluate((el, e) => (window as any).__bdd.firedCount(el, e), event),
    {timeout: pollMs(5000), message: `"${event}" did not fire on ${target.phrase} (listen for it before the gesture)`}).toBeGreaterThan(0);
  await loc.evaluate((el, e) => { const b = (window as any).__bdd; b.unlisten(b.viewerOf(el), e); }, event);
}

/** Not fired so far; the subscription stays for a later "should have fired". */
export async function expectNotFired(page: Page, target: ElementRef, event: string): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  const count: number = await loc.evaluate((el, e) => (window as any).__bdd.firedCount(el, e), event);
  expect(count, count < 0 ? `"${event}" is not listened for on ${target.phrase}` : `"${event}" fired ${count} time(s) on ${target.phrase}`).toBe(0);
}

export async function resize(page: Page, target: ElementRef, width: number | null, height: number | null): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await loc.evaluate((el, [w, h]) => (window as any).__bdd.resize(el, w, h, 500), [width, height] as [number | null, number | null]);
}

export async function restoreSize(page: Page, target: ElementRef): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await loc.evaluate((el) => (window as any).__bdd.restoreSize(el, 500));
}

export async function saveLayout(page: Page): Promise<void> {
  await installViewerRuntime(page);
  await page.evaluate(() => { (window as any).__bdd.saveLayout(); });
}

/** Saves the layout through the server and returns its id; the caller registers the deletion. */
export async function saveLayoutToServer(page: Page): Promise<string> {
  await installViewerRuntime(page);
  return page.evaluate(() => (window as any).__bdd.saveLayoutToServer());
}

export async function deleteLayout(page: Page, id: string): Promise<void> {
  await page.evaluate((i) => (window as any).__bdd.deleteLayout(i), id);
}

export async function loadLayout(page: Page): Promise<void> {
  await installViewerRuntime(page);
  await page.evaluate(() => (window as any).__bdd.loadLayout());
}

/** The balloons (info, warning, error) shown since the last read; reading clears them. */
export async function takeBalloons(page: Page): Promise<Balloon[]> {
  await installViewerRuntime(page);
  return page.evaluate(() => (window as any).__bdd.takeBalloons());
}

// --- context menus ----------------------------------------------------------------------------------

/** Arms a `grok.events.<name>` subscription before a gesture; the returned function waits for it. */
export async function armEvent(page: Page, name: string, capMs = 3000): Promise<() => Promise<boolean>> {
  await installViewerRuntime(page);
  const token: string = await page.evaluate(([n, cap]) => (window as any).__bdd.armEvent(n, cap), [name, capMs] as [string, number]);
  return () => page.evaluate((t) => (window as any).__bdd.waitArmed(t), token);
}

export async function closeContextMenu(page: Page): Promise<void> {
  await installViewerRuntime(page);
  await page.evaluate(() => (window as any).__bdd.closeMenu());
  await expect(page.locator(POPUP)).toHaveCount(0);
}

/** A real right-click at an armed point; resolves once the platform says the popup is in the DOM
 * (`onContextMenuShown` — `onContextMenu` fires before the menu exists). */
async function rightClickArmed(page: Page, x: number, y: number, token: string): Promise<Locator> {
  await page.mouse.click(x, y, {button: 'right'});
  if (!await page.evaluate((t) => (window as any).__bdd.waitArmed(t), token))
    throw new Error(`no context menu opened at (${Math.round(x)}, ${Math.round(y)})`);
  return page.locator(POPUP).last();
}

export async function openContextMenuAt(page: Page, x: number, y: number): Promise<Locator> {
  await installViewerRuntime(page);
  const token: string = await page.evaluate((cap) => (window as any).__bdd.openMenu(cap), 3000);
  return rightClickArmed(page, x, y, token);
}

/** The context menu of an element: at a named hit area, else its `view` area when it reports
 * one, else its centre. Three roundtrips: the point and the arming in one, the click, the wait. */
export async function openContextMenuOf(page: Page, target: ElementRef, area?: string): Promise<Locator> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  const {x, y, token} = await loc.evaluate((el, [a, cap]) => (window as any).__bdd.menuPoint(el, a, cap), [area ?? null, 3000] as [string | null, number]);
  return rightClickArmed(page, x, y, token);
}

/** Every label of the open menus reading exactly that — a viewer menu can hold two groups of the
 * same name (the axis "Annotations" and the viewer's own), and a group's children are in the same
 * popup, hidden until it opens. */
export function menuLabelsMatching(page: Page, label: string): Locator {
  return page.locator(POPUP).locator('.d4-menu-item-label', {hasText: exactText(label)});
}

/** Every item of the open menus carrying that label. */
export function menuItems(page: Page, label: string): Locator {
  return menuLabelsMatching(page, label).locator(MENU_ITEM);
}

/** Whether that label is on screen right now — several items can carry it and all but one be a
 * closed group's child, so the question is whether any match is visible, not whether the first is. */
async function menuShows(page: Page, label: string): Promise<boolean> {
  return (await menuLabelsMatching(page, label).filter({visible: true}).count()) > 0;
}

/** A menu item of the open popup by its own label (a group item also contains its children's). */
export function menuItem(page: Page, label: string): Locator {
  return menuItems(page, label).first();
}

/** A group item opens on a pointer move over it: the move enters from the left, since a pointer
 * already resting on the item (a hover before the menu was reopened) would move nowhere. The
 * group's items live in the same popup, hidden until it opens, so the wait is on how many labels
 * are visible; when several items share the label, each is tried until [wanted] shows up. */
async function openGroup(page: Page, label: string, wanted?: string): Promise<void> {
  // a group already open, or one the menu renders inline, needs no hover
  if (wanted !== undefined && await menuShows(page, wanted))
    return;
  const candidates = menuItems(page, label);
  await candidates.filter({visible: true}).first().waitFor({state: 'visible', timeout: 5000});
  const count = await candidates.count();
  for (let i = 0; i < count; i++) {
    const item = candidates.nth(i);
    if (!(await item.isVisible()))
      continue;
    const box = await item.boundingBox();
    if (!box)
      continue;
    const cy = box.y + box.height / 2;
    // enter from the left: a pointer already resting on the item would move nowhere
    await page.mouse.move(Math.max(0, box.x - 8), cy);
    await page.mouse.move(box.x + box.width / 2, cy);
    if (wanted === undefined)
      return;
    const opened = await expect.poll(() => menuShows(page, wanted), {timeout: pollMs(2000)}).toBe(true).then(() => true, () => false);
    if (opened) {
      // step into the flyout along the group's own row: the menu hides a submenu when the pointer
      // leaves the group item at a steep angle, which would take the item away mid-click
      const target = await menuLabelsMatching(page, wanted).filter({visible: true}).first().boundingBox();
      if (target != null) {
        await page.mouse.move(target.x + 4, cy);
        await page.mouse.move(target.x + 4, target.y + target.height / 2);
      }
      return;
    }
  }
  const visible = await page.locator(POPUP).locator('.d4-menu-item-label').filter({visible: true}).allTextContents();
  throw new Error(`the "${label}" group did not show "${wanted}"; the menu shows: ${visible.map((t) => t.trim()).filter(Boolean).join(' | ')}`);
}

/** `Misc > Show Inside Values`: hovers the groups, clicks the leaf. */
export async function pickMenuPath(page: Page, path: string): Promise<void> {
  const segments = path.split(SEP).filter((s) => s.length > 0);
  for (let i = 0; i < segments.length; i++) {
    const last = i === segments.length - 1;
    // the item ancestor is what the platform wires the click to, but an inline group's rows are
    // laid out so that only the label itself has a box
    // the item ancestor is the usual click target, but a group's rows lay the box out on the
    // label itself, so the ancestor is there and not clickable
    const act = last
      ? menuItems(page, segments[i]).filter({visible: true}).first().click({timeout: 3000})
        .catch(() => menuLabelsMatching(page, segments[i]).filter({visible: true}).first().click({timeout: 3000}))
      : openGroup(page, segments[i], segments[i + 1]);
    await act.catch(async () => {
      const visible = await page.locator(POPUP).locator('.d4-menu-item-label').filter({visible: true}).allTextContents();
      throw new Error(`no "${segments[i]}" in the menu; it shows: ${visible.map((s) => s.trim()).filter(Boolean).join(' | ')}`);
    });
  }
}

/** The labels of the open menu, opening every group of [path] first ("" reads the top level). */
export async function menuLabels(page: Page, path: string): Promise<string[]> {
  const segments = path.split(SEP).filter((s) => s.length > 0);
  for (let i = 0; i < segments.length; i++)
    await openGroup(page, segments[i], segments[i + 1]);
  const labels = await page.locator(POPUP).locator('.d4-menu-item-label').filter({visible: true}).allTextContents();
  return labels.map((t) => t.trim()).filter(Boolean);
}

// --- area gestures and geometry ----------------------------------------------------------------------

/** A plain drag from the centre of one hit area to the centre of another (a column header to a
 * new place, a range handle to a bin); the baseline is taken before the drag. */
export async function dragArea(page: Page, target: ElementRef, from: string, to: string): Promise<void> {
  const a = centerOf(await hitArea(page, target, from, true));
  const b = centerOf(await hitArea(page, target, to));
  await page.mouse.move(a.x, a.y);
  await page.mouse.down();
  await page.mouse.move(b.x, b.y, {steps: 3});
  await page.mouse.up();
}

/** A drag of the area's centre by a distance in a direction (a resizer, a splitter). */
export async function dragAreaBy(page: Page, target: ElementRef, area: string, px: number, direction: string): Promise<void> {
  const d = direction.toLowerCase();
  if (!['left', 'right', 'up', 'down'].includes(d))
    throw new Error(`a drag goes left, right, up or down, not "${direction}"`);
  const c = centerOf(await hitArea(page, target, area, true));
  const dx = d === 'left' ? -px : d === 'right' ? px : 0;
  const dy = d === 'up' ? -px : d === 'down' ? px : 0;
  await page.mouse.move(c.x, c.y);
  await page.mouse.down();
  await page.mouse.move(c.x + dx / 2, c.y + dy / 2);
  await page.mouse.move(c.x + dx, c.y + dy);
  await page.mouse.up();
}

/** A drag with a modifier across the inner 80% of an area: Control+Shift removes from the
 * selection, Alt zooms. */
export async function dragBoxOverArea(page: Page, target: ElementRef, area: string, keys: string[]): Promise<void> {
  const b = await hitArea(page, target, area, true);
  for (const k of keys)
    await page.keyboard.down(k);
  try {
    await page.mouse.move(b.x + b.width * 0.1, b.y + b.height * 0.1);
    await page.mouse.down();
    await page.mouse.move(b.x + b.width * 0.9, b.y + b.height * 0.9, {steps: 3});
    await page.mouse.up();
  }
  finally {
    for (const k of [...keys].reverse())
      await page.keyboard.up(k);
  }
}

/** Types into a hit area that holds an editor (a range input, a form field): a click on its
 * centre, select all, the text, Enter. */
/** A click on the area, the text typed over what the editor there holds, Enter. The click must have
 * put the focus into an editor inside the viewer — a histogram's range input took the click and not
 * the focus once in twenty runs, and the text then opened a cell editor on the grid, unseen — so it
 * is repeated, at the area's current place, until one did. */
export async function typeIntoArea(page: Page, target: ElementRef, area: string, text: string): Promise<void> {
  const loc = await viewerLocator(page, target);
  // the editor is pinned by a mark of its own rather than by `:focus`: the focus can leave it while
  // the text is being read back (a table view focuses its grid a second after it opens, whatever
  // the user is doing), and a locator on `:focus` would then wait on nothing
  const mark = `bdd-editor-${Date.now()}`;
  let where = '';
  try {
    await expect.poll(async () => {
      const c = centerOf(await hitArea(page, target, area, true));
      await page.mouse.click(c.x, c.y);
      where = await loc.evaluate((el, m) => {
        const a = document.activeElement;
        if (!a)
          return 'nothing';
        if (el.contains(a) && (a.tagName === 'INPUT' || a.tagName === 'TEXTAREA' || (a as HTMLElement).isContentEditable)) {
          a.setAttribute('data-bdd-editor', m);
          return '';
        }
        const name = a.getAttribute('name');
        const cls = String(a.className ?? '').trim();
        return a.tagName.toLowerCase() + (name ? `[name="${name}"]` : '') + (cls ? '.' + cls.split(/\s+/).join('.') : '');
      }, mark);
      return where === '';
    }, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`a click on the "${area}" area of ${target.phrase} did not focus an editor there; the focus is on ${where}`);
  }
  const editor = loc.locator(`[data-bdd-editor="${mark}"]`);
  try {
    await typeVerified(editor, text.replace(/\\n/g, '\n'), `the "${area}" area of ${target.phrase}`);
    await editor.press('Enter');
  }
  finally {
    await editor.evaluate((e) => e.removeAttribute('data-bdd-editor')).catch(() => undefined);
  }
}

export async function expectAreaSize(page: Page, target: ElementRef, area: string, dimension: 'tall' | 'wide', min: number): Promise<void> {
  const box = await hitArea(page, target, area);
  const size = dimension === 'tall' ? box.height : box.width;
  expect(size, `the "${area}" area of ${target.phrase} is ${Math.round(size)} px ${dimension}, not at least ${min}`).toBeGreaterThanOrEqual(min);
}

/** The area's rectangle against the snapshot's: taller or wider than before. */
export async function expectAreaGrew(page: Page, target: ElementRef, area: string, dimension: 'taller' | 'wider' | 'shorter' | 'narrower'): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  let last: {before?: Box; now?: Box; has: string[]} = {has: []};
  const holds = async (): Promise<boolean | string> => {
    last = await loc.evaluate((el, a) => (window as any).__bdd.areaRectChange(el, a), area);
    if (!last.now)
      return `no "${area}" area now`;
    if (!last.before)
      return `no "${area}" area at the snapshot`;
    const now = dimension === 'taller' || dimension === 'shorter' ? last.now.height : last.now.width;
    const before = dimension === 'taller' || dimension === 'shorter' ? last.before.height : last.before.width;
    return dimension === 'taller' || dimension === 'wider' ? now > before + 0.5 : now < before - 0.5;
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the "${area}" area of ${target.phrase} is not ${dimension} than before (before ${JSON.stringify(last.before)}, now ${JSON.stringify(last.now)}` +
      (last.now ? ')' : `; it has: ${last.has.join(', ') || 'none'})`));
  }
}

/** Two areas are painted in the same colors: every significant color of either has a near color
 * in the other (a linked color coding, a category and its swatch). */
export async function expectSameColors(page: Page, target: ElementRef, a: string, b: string): Promise<void> {
  await hitArea(page, target, a);
  await hitArea(page, target, b);
  let shownA: AreaColors | undefined;
  let shownB: AreaColors | undefined;
  const own = (mine: AreaColor[], theirs: AreaColor[]): boolean =>
    mine.some((c) => c.count >= SIGNIFICANT_PX && !theirs.some((d) => d.count >= COLOR_MIN_PX && near(c.hex, d.hex)));
  const same = async (): Promise<boolean> => {
    const ra: AreaColors = shownA = await areaColors(page, target, a);
    const rb: AreaColors = shownB = await areaColors(page, target, b);
    return !own(ra.colors, rb.colors) && !own(rb.colors, ra.colors);
  };
  try {
    await expect.poll(same, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the "${a}" and "${b}" areas of ${target.phrase} are painted in different colors; "${a}" ${describeColors(shownA)}; "${b}" ${describeColors(shownB)}`);
  }
}

/** No pixel of the color (nor a shade of it) inside the area, read once. */
export async function expectAreaNotColor(page: Page, target: ElementRef, area: string, color: string): Promise<void> {
  const want = '#' + color.replace(/^#/, '').toUpperCase();
  if (!/^#[0-9A-F]{6}$/.test(want))
    throw new Error(`"${color}" is not a #rrggbb color`);
  await hitArea(page, target, area);
  const read = await areaColors(page, target, area);
  const count = read.colors.filter((c) => near(c.hex, want)).reduce((n, c) => n + c.count, 0);
  expect(count, `the "${area}" area of ${target.phrase} is painted in ${want} (${count} px); ${describeColors(read)}`).toBeLessThan(COLOR_MIN_PX);
}

export async function rememberReading(page: Page, target: ElementRef, name: string): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await loc.evaluate((el, n) => { (window as any).__bdd.rememberValue(el, n); }, name);
}

export async function expectRememberedReading(page: Page, target: ElementRef, name: string, not = false): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  let last: {before?: unknown; now?: unknown; has: string[]} = {has: []};
  const holds = async (): Promise<boolean | string> => {
    last = await loc.evaluate((el, n) => (window as any).__bdd.rememberedValue(el, n), name);
    if (last.before === undefined)
      return `"${name}" was not remembered`;
    if (last.now === undefined || last.now === null)
      return `no "${name}" reading`;
    return last.now === last.before;
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(!not);
  }
  catch {
    if (not)
      throw new Error(`"${name}" of ${target.phrase} is still the remembered ${String(last.before)}`);
    throw new Error(`"${name}" of ${target.phrase} is ${String(last.now)}, not the remembered ${String(last.before)}` +
      (last.now === undefined || last.now === null ? `; the viewer reports: ${last.has.join(', ') || 'no readings'}` : ''));
  }
}

// --- the legend --------------------------------------------------------------------------------------

interface LegendState {
  mode: string;
  slot: string;
  items: number;
  keys: string[];
  selected: string[];
  width: number;
  height: number;
}

/** The legend the viewer hosts, wherever it sits now (in the viewer, or re-parented into the
 * tooltip when collapsed to the mini icon). */
function legendRoot(page: Page, viewer: Locator): Locator {
  return viewer.locator('[name="legend"]').or(page.locator('.d4-tooltip [name="legend"]')).filter({visible: true}).first();
}

/** A legend item by the label it shows (`aria-label`, else the label text). */
export async function legendItem(page: Page, target: ElementRef, label: string): Promise<Locator> {
  const legend = legendRoot(page, await viewerLocator(page, target));
  const items = legend.locator('[name="legend-item"]');
  const byAria = items.filter({has: page.locator(`:scope[aria-label="${label.replace(/"/g, '\\"')}" i]`)});
  const byText = items.filter({has: page.locator('.d4-legend-value', {hasText: exactText(label)})});
  return byAria.or(byText).first();
}

export async function legendStateOf(page: Page, target: ElementRef): Promise<LegendState | undefined> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  return loc.evaluate((el) => (window as any).__bdd.legendState(el));
}

function describeLegend(s: LegendState | undefined): string {
  return s ? `mode ${s.mode || 'none'}, slot ${s.slot || 'none'}, ${s.items} items (${s.keys.length} rendered), ${s.width}x${s.height}` : 'no legend';
}

/** The legend's item total (`data-legend-items`, every section, rendered or not). */
export async function expectLegendItems(page: Page, target: ElementRef, count: number): Promise<void> {
  let last: LegendState | undefined;
  try {
    await expect.poll(async () => (last = await legendStateOf(page, target))?.items, {timeout: pollMs(5000)}).toBe(count);
  }
  catch {
    throw new Error(`the legend of ${target.phrase} lists ${last?.items ?? 'no'} items, not ${count} (${describeLegend(last)})`);
  }
}

export async function expectLegendItemsChange(page: Page, target: ElementRef, compare: 'fewer' | 'same'): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  let last: {before?: LegendState; now?: LegendState} = {};
  const holds = async (): Promise<boolean | string> => {
    last = await loc.evaluate((el) => (window as any).__bdd.legendChange(el));
    if (!last.before)
      return 'no legend at the snapshot';
    if (!last.now)
      return 'no legend now';
    return compare === 'fewer' ? last.now.items < last.before.items :
      last.now.items === last.before.items && last.now.keys.join('|') === last.before.keys.join('|');
  };
  try {
    await expect.poll(holds, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the legend of ${target.phrase} does not list ${compare === 'fewer' ? 'fewer items than' : 'the same items as'} before (before ${describeLegend(last.before)}; now ${describeLegend(last.now)})`);
  }
}

const LEGEND_MODES: Record<string, string> = {docked: 'docked', corner: 'corner', 'mini icon': 'miniIcon', tooltip: 'tooltip', hidden: 'hidden'};

export async function expectLegendMode(page: Page, target: ElementRef, mode: string): Promise<void> {
  const want = LEGEND_MODES[mode];
  if (!want)
    throw new Error(`a legend is docked, in a corner, collapsed to the mini icon, shown in the tooltip or hidden — not "${mode}"`);
  let last: LegendState | undefined;
  try {
    await expect.poll(async () => (last = await legendStateOf(page, target))?.mode, {timeout: pollMs(5000)}).toBe(want);
  }
  catch {
    throw new Error(`the legend of ${target.phrase} is not ${mode} (${describeLegend(last)})`);
  }
}

export async function expectLegendSlot(page: Page, target: ElementRef, slot: string, negate = false): Promise<void> {
  const norm = (s: string) => s.toLowerCase().replace(/[^a-z]/g, '');
  let last: LegendState | undefined;
  const poll = expect.poll(async () => norm((last = await legendStateOf(page, target))?.slot ?? ''), {timeout: pollMs(5000)});
  try {
    await (negate ? poll.not : poll).toBe(norm(slot));
  }
  catch {
    throw new Error(`the legend of ${target.phrase} is ${negate ? 'still' : 'not'} in the ${slot} slot (${describeLegend(last)})`);
  }
}

/** Mode, slot and size as at the snapshot before the last change, read once the viewer is quiet:
 * a legend that keeps its slot but jumps or resizes inside it has moved. */
export async function expectLegendPlacedAsBefore(page: Page, target: ElementRef): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  const r: {before?: LegendState; now?: LegendState} = await loc.evaluate(async (el) => {
    const b = (window as any).__bdd;
    await b.quiet(b.viewerOf(el));
    return b.legendChange(el);
  });
  if (!r.before || !r.now)
    throw new Error(`the legend of ${target.phrase}: ${!r.before ? 'no legend at the snapshot' : 'no legend now'}`);
  const where = (l: LegendState) => `${l.mode}/${l.slot}/${Math.round(l.width)}x${Math.round(l.height)}`;
  expect(where(r.now), `the legend of ${target.phrase} moved (before ${describeLegend(r.before)}; now ${describeLegend(r.now)})`).toBe(where(r.before));
}

/** A click on a legend item (the category filters the viewer): the baseline is taken first and
 * the viewer settles after, so the checks that follow read the state after the repaint. */
export async function clickLegendItem(page: Page, target: ElementRef, label: string, options: {key?: string; cross?: boolean} = {}): Promise<void> {
  await snapshot(page, target);
  const item = await legendItem(page, target, label);
  await item.waitFor({state: 'visible', timeout: 5000}).catch(async () => {
    const shown = await legendRoot(page, await viewerLocator(page, target)).locator('.d4-legend-value').allTextContents();
    throw new Error(`no "${label}" item in the legend of ${target.phrase}; it shows: ${shown.map((s) => s.trim()).filter(Boolean).join(' | ') || 'nothing'}`);
  });
  const what = options.cross ? item.locator('.d4-legend-cross').first() : item.locator('.d4-legend-value').first();
  if (options.cross)
    await item.hover();
  const keys = options.key ? options.key.split('+') : [];
  for (const k of keys)
    await page.keyboard.down(k);
  try {
    await what.click();
  }
  finally {
    for (const k of [...keys].reverse())
      await page.keyboard.up(k);
  }
  const loc = await viewerLocator(page, target);
  await loc.evaluate((el) => (window as any).__bdd.settle(el, 300));
}

function cssToHex(color: string): string {
  const m = /rgba?\((\d+),\s*(\d+),\s*(\d+)/.exec(color);
  if (m)
    return '#' + [m[1], m[2], m[3]].map((x) => Number(x).toString(16).padStart(2, '0')).join('').toUpperCase();
  const h = /^#([0-9a-f]{6})$/i.exec(color.trim());
  return h ? '#' + h[1].toUpperCase() : '';
}

/** The color a legend item is drawn in (its inline color, the category's color). */
export async function legendItemColor(page: Page, target: ElementRef, label: string): Promise<string> {
  const item = await legendItem(page, target, label);
  await item.waitFor({state: 'visible', timeout: 5000});
  return cssToHex(await item.evaluate((e) => (e as HTMLElement).style.color || getComputedStyle(e).color));
}

export async function expectLegendItemColor(page: Page, target: ElementRef, label: string, color: string, negate = false): Promise<void> {
  const want = '#' + color.replace(/^#/, '').toUpperCase();
  if (!/^#[0-9A-F]{6}$/.test(want))
    throw new Error(`"${color}" is not a #rrggbb color`);
  let last = '';
  try {
    await expect.poll(async () => near(last = await legendItemColor(page, target, label), want), {timeout: pollMs(5000)}).toBe(!negate);
  }
  catch {
    throw new Error(`the "${label}" item in the legend of ${target.phrase} is ${negate ? 'still' : 'not'} colored ${want}` +
      (negate ? '' : `; it is ${last || 'colored nothing'}`));
  }
}

export async function expectLegendItemsDiffer(page: Page, target: ElementRef, a: string, b: string): Promise<void> {
  let ca = '';
  let cb = '';
  try {
    await expect.poll(async () => {
      ca = await legendItemColor(page, target, a);
      cb = await legendItemColor(page, target, b);
      return ca !== '' && cb !== '' && !near(ca, cb);
    }, {timeout: pollMs(5000)}).toBe(true);
  }
  catch {
    throw new Error(`the "${a}" and "${b}" items in the legend of ${target.phrase} are colored alike (${ca || 'nothing'} and ${cb || 'nothing'})`);
  }
}

/** A real drag of the legend splitter: along the axis the bar resizes (a vertical bar moves
 * horizontally), positive towards the legend's far side; the viewer settles after. */
export async function dragLegendSplitter(page: Page, target: ElementRef, px: number): Promise<void> {
  await snapshot(page, target);
  const loc = await viewerLocator(page, target);
  const bar = loc.locator('[name="legend-splitter"]').filter({visible: true}).first();
  await bar.waitFor({state: 'visible', timeout: 5000}).catch(() => {
    throw new Error(`${target.phrase} shows no legend splitter (a docked legend has one)`);
  });
  const vertical = (await bar.getAttribute('class') ?? '').includes('vertical');
  const box = await bar.boundingBox();
  if (!box)
    throw new Error('the legend splitter has no box');
  const c = {x: box.x + box.width / 2, y: box.y + box.height / 2};
  await page.mouse.move(c.x, c.y);
  await page.mouse.down();
  await page.mouse.move(c.x + (vertical ? px / 2 : 0), c.y + (vertical ? 0 : px / 2));
  await page.mouse.move(c.x + (vertical ? px : 0), c.y + (vertical ? 0 : px));
  await page.mouse.up();
  await loc.evaluate((el) => (window as any).__bdd.settle(el, 300));
}

// --- tooltips ---------------------------------------------------------------------------------------

/** The value the row tooltip shows for a column name. */
export async function tooltipValue(page: Page, column: string): Promise<string | undefined> {
  const rows = await page.locator(TOOLTIP_ROWS).evaluateAll((trs) =>
    trs.map((tr) => Array.from(tr.querySelectorAll('td')).map((td) => (td.textContent ?? '').trim())));
  const hit = rows.find((cells) => cells[0]?.toUpperCase() === column.toUpperCase());
  return hit ? hit.slice(1).join(' ').trim() : undefined;
}

export async function expectTooltipValue(page: Page, column: string, value: string): Promise<void> {
  await expect.poll(() => tooltipValue(page, column), {timeout: pollMs(5000), message: `"${column}" in the row tooltip`}).toBe(value);
}

export async function tooltipColumns(page: Page): Promise<string[]> {
  const cells = await page.locator(TOOLTIP_COLUMNS).allTextContents();
  return [...new Set(cells.map((c) => c.trim().toUpperCase()).filter((c) => c.length > 0))].sort();
}

export async function expectTooltipColumns(page: Page, list: string, negate = false): Promise<void> {
  const want = [...new Set(list.split(/\s*,\s*/).map((c) => c.trim().toUpperCase()).filter((c) => c.length > 0))].sort();
  const poll = expect.poll(() => tooltipColumns(page), {timeout: pollMs(5000), message: negate ? 'the tooltip shows exactly these columns' : 'the tooltip does not show these columns'});
  await (negate ? poll.not : poll).toEqual(want);
}
