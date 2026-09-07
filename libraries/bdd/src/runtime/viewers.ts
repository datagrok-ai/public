/* Viewer runtime: what the `viewers` tier's steps are made of — an in-page helper set installed
   once per page (`window.__bdd`) and the Playwright-side readers over it.
   The principle: this is our platform, not a black box. A viewer that waits out a debounce, a menu
   that builds on a timer, a canvas that paints on the next frame — the platform can tell us when
   (its events) or stop waiting altogether (`immediateRendering`, armed here on every viewer the
   page will ever hold), so nothing in this module sleeps. When a signal or a name is missing, it
   is added to the core, never faked here. */
import {expect, Locator, Page} from '@playwright/test';
import type {ElementRef} from './args.js';
import {exactText, locate} from './locate.js';

declare const grok: any;
declare const DG: any;

export interface Box {
  x: number;
  y: number;
  width: number;
  height: number;
}

/** A viewer's canvas after a change: `delta` is the histogram distance from the snapshot taken
 * before the change, `ink` the painted pixels now, `inkBefore` those of the snapshot, `hue` and
 * `hueBefore` the pixels in the selection hue; `selected` rows and the `view` area in device
 * pixels size the highlight a selection must paint. */
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
const MENU_ITEM = 'xpath=ancestor::*[contains(concat(" ", normalize-space(@class), " "), " d4-menu-item ")][1]';
const TOOLTIP_COLUMNS = '.d4-tooltip table.d4-row-tooltip-table tr td:first-child';

/** Everything the in-page side needs, on `window.__bdd`. Self-contained: it runs in the browser. */
function install(): void {
  const w = window as any;
  if (w.__bdd)
    return;
  const renders = new WeakMap<Element, {count: number; last: number; sub?: any}>();
  const snapshots = new WeakMap<Element, {colors: Map<number, number>; ink: number; hue: number; renders: number; range?: Range;
    areas: Record<string, number>; scale?: ScaleRange; values: Record<string, unknown>}>();
  // how long a repaint the viewer says is pending may take before that is a platform failure
  const PENDING_CAP = 10000;
  const sizes = new WeakMap<Element, {width: string; height: string}>();
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
    try {
      stamp.sub = v.onViewerRendered.subscribe(() => {
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
   * once and a repaint is waited for as long as it takes (a pending one that never lands is a
   * platform failure, reported after `PENDING_CAP`). A viewer without the signal falls back to a
   * render event within `capMs`. */
  const settle = (el: Element, capMs: number): Promise<number> => {
    const v = viewerOf(el);
    arm(v);
    const stamp = renders.get(v.root)!;
    const before = stamp.count;
    const t0 = Date.now();
    return new Promise((resolve, reject) => {
      const tick = () => {
        if (stamp.count > before)
          return resolve(stamp.count - before);
        const p = pending(v);
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
  const canvasOf = (v: any): HTMLCanvasElement => {
    const part = v.getWidgetStatus?.()?.parts?.canvas;
    const cv = part ?? v.root.querySelector('canvas[name="canvas"]') ?? v.root.querySelector('canvas');
    if (!cv)
      throw new Error(`${v.type} has no canvas`);
    return cv;
  };
  const pixels = (cv: HTMLCanvasElement): ImageData => cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height);
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
    const data = histogram(pixels(canvasOf(viewerOf(el))));
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
    const colors = colorsIn(pixels(cv), rect);
    return {colors: [...colors.entries()].sort((a, b) => b[1] - a[1]).map(([c, count]) => ({hex: hex(c), count})), rect, bitmap: [cv.width, cv.height]};
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
    const cv = canvasOf(v);
    const img = pixels(cv);
    const shot = histogram(img);
    const areas: Record<string, number> = {};
    const hit = areasOf(v);
    for (const key of Object.keys(hit))
      areas[norm(key)] = inkIn(img, deviceRect(cv, hit[key]));
    snapshots.set(v.root, {...shot, renders: renders.get(v.root)!.count, range: rangeOf(v), areas, scale: scaleRange(v), values: {...valuesOf(v)}});
    return shot.ink;
  };
  // the baseline "should have repainted" compares with: the canvas before the change
  const baseline = (el: Element): void => {
    try {
      snapshot(el);
    } catch { /* a viewer without a canvas */ }
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
    baseline(el);
    const settled = settle(el, capMs);
    for (const [caption, text] of entries) {
      const p = property(v, caption);
      v.props[p.name] = convert(p, text);
    }
    return settled;
  };
  const canvasBox = (v: any): Box => {
    const r = canvasOf(v).getBoundingClientRect();
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
  /** Painted pixels inside a hit area (the canvas may be scaled to the device). */
  const areaInk = (el: Element, name: string): number => {
    const v = viewerOf(el);
    const cv = canvasOf(v);
    return inkIn(pixels(cv), deviceRect(cv, areasOf(v)[areaKey(v, name)]));
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
  const change = (el: Element): CanvasChange => {
    const v = viewerOf(el);
    const cv = canvasOf(v);
    const now = histogram(pixels(cv));
    const before = snapshots.get(v.root);
    if (!before)
      throw new Error(`${v.type}: no snapshot to compare with`);
    let delta = 0;
    for (const [c, n] of now.colors)
      delta += Math.abs(n - (before.colors.get(c) ?? 0));
    for (const [c, n] of before.colors) {
      if (!now.colors.has(c))
        delta += n;
    }
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
  const stillness = async (el: Element): Promise<{renders: number; delta: number}> => {
    const v = viewerOf(el);
    const before = snapshots.get(v.root);
    if (!before)
      throw new Error(`${v.type}: no snapshot to compare with`);
    await quiet(v);
    return {renders: renders.get(v.root)!.count - before.renders, delta: change(el).delta};
  };
  /** The value range against the snapshot's, read once the viewer is quiet (a reset that lands a
   * tick after the change is read, not missed). */
  const quietRangeChange = async (el: Element): Promise<RangeChange> => {
    await quiet(viewerOf(el));
    return rangeChange(el);
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
  // the repaint a resize causes lands on the next task, so the settle is armed before the event
  const resize = (el: Element, width: number | null, height: number | null, capMs: number): Promise<number> => {
    const root = viewerOf(el).root as HTMLElement;
    if (!sizes.has(root))
      sizes.set(root, {width: root.style.width, height: root.style.height});
    baseline(el);
    const settled = settle(el, capMs);
    if (width !== null)
      root.style.width = `${width}px`;
    if (height !== null)
      root.style.height = `${height}px`;
    window.dispatchEvent(new Event('resize'));
    return settled;
  };
  const restoreSize = (el: Element, capMs: number): Promise<number> => {
    const root = viewerOf(el).root as HTMLElement;
    const size = sizes.get(root);
    if (!size)
      return Promise.resolve(0);
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
  const loadLayout = (): void => {
    if (!layout)
      throw new Error('no layout saved in this feature');
    grok.shell.tv.loadLayout(layout);
  };

  w.__bdd = {viewerOf, arm, stampAll, settle, quiet, readProperty, writeProperties, findArea, hitArea, areaInk, areaChange, areaColors,
    snapshot, baselineAll, change, rangeChange, quietRangeChange, scaleChange, valueChange, rememberRange, rememberedRange, stillness,
    palette, tableOf, listen, unlisten, firedCount, resize, restoreSize, armEvent, waitArmed, closeMenu, openMenu, menuPoint, addViewer,
    takeBalloons, saveLayout, loadLayout};
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
  let found = await find();
  if (!found.box) {
    await expect.poll(async () => (found = await find()).box !== undefined,
      {timeout: 5000, message: `${target.phrase} reports no "${name}" area; it has: ${found.has.join(', ') || 'none'}`}).toBe(true);
  }
  return found.box!;
}

export async function expectHasArea(page: Page, target: ElementRef, name: string, negate = false): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  const find = (): Promise<{box?: Box; has: string[]}> => loc.evaluate((el, n) => (window as any).__bdd.findArea(el, n, false), name);
  let found = await find();
  const poll = expect.poll(async () => (found = await find()).box !== undefined,
    {timeout: 5000, message: `${target.phrase} ${negate ? 'still reports' : 'reports no'} "${name}" area; it has: ${found.has.join(', ') || 'none'}`});
  await (negate ? poll.not : poll).toBe(true);
}

export function centerOf(box: Box): {x: number; y: number} {
  return {x: box.x + box.width / 2, y: box.y + box.height / 2};
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
  const poll = expect.poll(() => readProperty(page, target, caption, value), {timeout: 5000, message: `"${caption}" of ${target.phrase}`});
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
    {timeout: 10000, message: minPx > 1 ? `${target.phrase} did not repaint by ${minPx} pixels` : `${target.phrase} did not repaint`}).toBeGreaterThanOrEqual(minPx);
}

/** The table the viewer draws (`viewer.dataFrame`), not the property it was asked to bind. */
export async function expectBoundTable(page: Page, target: ElementRef, name: string): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await expect.poll(() => loc.evaluate((el) => (window as any).__bdd.tableOf(el)), {timeout: 5000, message: `the table ${target.phrase} is bound to`}).toBe(name);
}

/** A hit area's ink against the snapshot before the last change. */
export async function expectAreaInk(page: Page, target: ElementRef, area: string, compare: 'less' | 'more'): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  let last: AreaChange | undefined;
  const holds = async (): Promise<boolean> => {
    const c: AreaChange = last = await loc.evaluate((el, a) => (window as any).__bdd.areaChange(el, a), area);
    return compare === 'less' ? c.ink < c.inkBefore : c.ink > c.inkBefore;
  };
  try {
    await expect.poll(holds, {timeout: 10000}).toBe(true);
  }
  catch {
    throw new Error(`the "${area}" area of ${target.phrase} does not have ${compare} ink than before` +
      (last ? ` (${last.inkBefore} px before, ${last.ink} now, in ${JSON.stringify(last.rect)} of the bitmap)` : ''));
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
    await expect.poll(count, {timeout: 5000}).toBeGreaterThanOrEqual(COLOR_MIN_PX);
  }
  catch {
    throw new Error(`the "${area}" area of ${target.phrase} is not painted in ${want}; ${describeColors(last)}`);
  }
}

const SIGNIFICANT_PX = 30;

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
    await expect.poll(differ, {timeout: 5000}).toBe(true);
  }
  catch {
    throw new Error(`the "${a}" and "${b}" areas of ${target.phrase} are painted in the same colors; "${a}" ${describeColors(shownA)}; "${b}" ${describeColors(shownB)}`);
  }
}

/** A numeric reading of the viewer (`getWidgetStatus().values`, "rows shown") equals a value, or
 * is lower/higher than at the snapshot before the last change. */
export async function expectReading(page: Page, target: ElementRef, name: string, compare: 'equal' | 'lower' | 'higher', value?: number): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  let last: {before?: unknown; now?: unknown; has: string[]} = {has: []};
  const holds = async (): Promise<boolean | string> => {
    last = await loc.evaluate((el, n) => (window as any).__bdd.valueChange(el, n), name);
    if (typeof last.now !== 'number')
      return `no "${name}" reading`;
    if (compare === 'equal')
      return last.now === value;
    if (typeof last.before !== 'number')
      return `no "${name}" reading at the snapshot`;
    return compare === 'lower' ? last.now < last.before : last.now > last.before;
  };
  try {
    await expect.poll(holds, {timeout: 5000}).toBe(true);
  }
  catch {
    const what = compare === 'equal' ? `${value}` : `${compare} than before (${String(last.before)})`;
    throw new Error(`"${name}" of ${target.phrase} is ${String(last.now)}, not ${what}` +
      (typeof last.now !== 'number' ? `; the viewer reports: ${last.has.join(', ') || 'no readings'}` : ''));
  }
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
    await expect.poll(holds, {timeout: 5000}).toBe(true);
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
  const still: {renders: number; delta: number} = await loc.evaluate((el) => (window as any).__bdd.stillness(el));
  expect(still.delta, `${target.phrase} repainted: ${still.delta} px changed in ${still.renders} render(s)`).toBe(0);
}

export async function expectInk(page: Page, target: ElementRef, compare: 'less' | 'more' | 'some'): Promise<void> {
  if (compare === 'some') {
    await expect.poll(() => snapshot(page, target), {timeout: 10000, message: `${target.phrase} is blank`}).toBeGreaterThan(0);
    return;
  }
  await expect.poll(async () => {
    const c = await canvasChange(page, target);
    return compare === 'less' ? c.ink < c.inkBefore : c.ink > c.inkBefore;
  }, {timeout: 10000, message: `${target.phrase} does not have ${compare} ink than before`}).toBe(true);
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
    await expect.poll(holds, {timeout: 10000}).toBe(true);
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

/** The value range equals the one remembered for this viewer type — across a close and a reopen. */
export async function expectRememberedRange(page: Page, target: ElementRef): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  let last: RangeChange = {};
  const holds = async (): Promise<boolean> => {
    last = await loc.evaluate((el) => (window as any).__bdd.rememberedRange(el));
    return !!last.before && !!last.now && Math.abs(last.before.top - last.now.top) < 0.5 && Math.abs(last.before.bottom - last.now.bottom) < 0.5;
  };
  try {
    await expect.poll(holds, {timeout: 5000}).toBe(true);
  }
  catch {
    throw new Error(`${target.phrase} does not show the remembered value range (remembered ${JSON.stringify(last.before)}, now ${JSON.stringify(last.now)})`);
  }
}

export async function expectPalette(page: Page, target: ElementRef, min: number): Promise<void> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  await expect.poll(() => loc.evaluate((el) => (window as any).__bdd.palette(el, 500)),
    {timeout: 5000, message: `${target.phrase} is painted in fewer than ${min} colors`}).toBeGreaterThanOrEqual(min);
}

export async function expectAreaPainted(page: Page, target: ElementRef, area: string): Promise<void> {
  await hitArea(page, target, area);
  const loc = await viewerLocator(page, target);
  await expect.poll(() => loc.evaluate((el, a) => (window as any).__bdd.areaInk(el, a), area),
    {timeout: 5000, message: `the "${area}" area of ${target.phrase} is blank`}).toBeGreaterThan(0);
}

export async function rangeChange(page: Page, target: ElementRef): Promise<RangeChange> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  return loc.evaluate((el) => (window as any).__bdd.rangeChange(el));
}

/** The value range (viewport) against the snapshot before the last change; `same` is read once the
 * viewer is quiet, so a reset that lands a tick later fails it rather than slipping past. */
export async function expectValueRange(page: Page, target: ElementRef, compare: 'narrower' | 'same' | 'wider'): Promise<void> {
  const same = (a: Range, b: Range) => Math.abs(a.top - b.top) < 1e-6 && Math.abs(a.bottom - b.bottom) < 1e-6;
  const loc = await viewerLocator(page, target);
  let last: RangeChange = {};
  const holds = async (): Promise<boolean | string> => {
    last = compare === 'same' ? await loc.evaluate((el) => (window as any).__bdd.quietRangeChange(el)) : await rangeChange(page, target);
    if (!last.before || !last.now)
      return 'no range';
    if (compare === 'same')
      return same(last.before, last.now);
    return compare === 'narrower' ? last.now.height < last.before.height * 0.95 : last.now.height > last.before.height * 1.05;
  };
  try {
    await expect.poll(holds, {timeout: 5000}).toBe(true);
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
    {timeout: 5000, message: `"${event}" did not fire on ${target.phrase} (listen for it before the gesture)`}).toBeGreaterThan(0);
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

export async function loadLayout(page: Page): Promise<void> {
  await installViewerRuntime(page);
  await page.evaluate(() => { (window as any).__bdd.loadLayout(); });
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

/** A menu item of the open popup by its own label (a group item also contains its children's). */
export function menuItem(page: Page, label: string): Locator {
  return page.locator(POPUP).last().locator('.d4-menu-item-label', {hasText: exactText(label)}).first().locator(MENU_ITEM);
}

/** `Misc > Show Inside Values`: hovers the groups, clicks the leaf. */
export async function pickMenuPath(page: Page, path: string): Promise<void> {
  const segments = path.split(/\s*[>|]\s*/).filter((s) => s.length > 0);
  for (let i = 0; i < segments.length; i++) {
    const item = menuItem(page, segments[i]);
    const act = i < segments.length - 1 ? item.hover({timeout: 5000}) : item.click({timeout: 5000});
    await act.catch(async (e: Error) => {
      const visible = await page.locator(POPUP).last().locator('.d4-menu-item-label').allTextContents();
      throw new Error(`no "${segments[i]}" in the menu; it shows: ${visible.map((s) => s.trim()).filter(Boolean).join(' | ')}\n${e.message}`);
    });
  }
}

// --- tooltips ---------------------------------------------------------------------------------------

export async function tooltipColumns(page: Page): Promise<string[]> {
  const cells = await page.locator(TOOLTIP_COLUMNS).allTextContents();
  return [...new Set(cells.map((c) => c.trim().toUpperCase()).filter((c) => c.length > 0))].sort();
}

export async function expectTooltipColumns(page: Page, list: string, negate = false): Promise<void> {
  const want = [...new Set(list.split(/\s*,\s*/).map((c) => c.trim().toUpperCase()).filter((c) => c.length > 0))].sort();
  const poll = expect.poll(() => tooltipColumns(page), {timeout: 5000, message: negate ? 'the tooltip shows exactly these columns' : 'the tooltip does not show these columns'});
  await (negate ? poll.not : poll).toEqual(want);
}
