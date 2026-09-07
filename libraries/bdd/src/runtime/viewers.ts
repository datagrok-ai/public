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
 * `hueBefore` the pixels in the selection hue. */
export interface CanvasChange {
  delta: number;
  ink: number;
  inkBefore: number;
  hue: number;
  hueBefore: number;
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
  const snapshots = new WeakMap<Element, {colors: Map<number, number>; ink: number; hue: number; renders: number; range?: Range}>();
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
  const settle = (el: Element, capMs: number): Promise<number> => {
    const v = viewerOf(el);
    arm(v);
    const stamp = renders.get(v.root)!;
    const before = stamp.count;
    const t0 = Date.now();
    return new Promise((resolve) => {
      const tick = () => {
        if (stamp.count > before || Date.now() - t0 >= capMs)
          resolve(stamp.count - before);
        else
          setTimeout(tick, 10);
      };
      setTimeout(tick, 0);
    });
  };
  const canvasOf = (v: any): HTMLCanvasElement => {
    const part = v.getWidgetStatus?.()?.parts?.canvas;
    const cv = part ?? v.root.querySelector('canvas[name="canvas"]') ?? v.root.querySelector('canvas');
    if (!cv)
      throw new Error(`${v.type} has no canvas`);
    return cv;
  };
  const histogram = (cv: HTMLCanvasElement) => {
    const data = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
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
  /** The colors drawn in at least `minPx` pixels, blanks and near-whites aside. */
  const palette = (el: Element, minPx: number): number => {
    const data = canvasOf(viewerOf(el)).getContext('2d')!.getImageData(0, 0, 1, 1) && histogram(canvasOf(viewerOf(el)));
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
  const snapshot = (el: Element): number => {
    const v = viewerOf(el);
    arm(v);
    const shot = histogram(canvasOf(v));
    snapshots.set(v.root, {...shot, renders: renders.get(v.root)!.count, range: rangeOf(v)});
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
    const areas: Record<string, Box> = v.getWidgetStatus()?.hitAreas ?? {};
    const key = Object.keys(areas).find((k) => norm(k) === norm(name));
    if (!key)
      throw new Error(`${v.type} has no "${name}" area right now; it has: ${Object.keys(areas).join(', ') || 'none'}`);
    const r = areas[key];
    const cv = canvasOf(v);
    const scale = cv.width / cv.getBoundingClientRect().width;
    const data = cv.getContext('2d')!.getImageData(Math.floor(r.x * scale), Math.floor(r.y * scale),
      Math.max(1, Math.floor(r.width * scale)), Math.max(1, Math.floor(r.height * scale))).data;
    let n = 0;
    for (let i = 0; i < data.length; i += 4) {
      if (!isBlank(data[i], data[i + 1], data[i + 2], data[i + 3]))
        n++;
    }
    return n;
  };
  const change = (el: Element): CanvasChange => {
    const v = viewerOf(el);
    const now = histogram(canvasOf(v));
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
    return {delta, ink: now.ink, inkBefore: before.ink, hue: now.hue, hueBefore: before.hue};
  };
  const rangeChange = (el: Element): RangeChange => {
    const v = viewerOf(el);
    return {before: snapshots.get(v.root)?.range, now: rangeOf(v)};
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
  /** Whether the viewer painted since its snapshot, read after the frame a change would land on:
   * a repaint is scheduled on the next task, so one animation frame and one task later there is
   * nothing left to wait for. */
  const stillness = async (el: Element): Promise<{renders: number; delta: number}> => {
    const v = viewerOf(el);
    const before = snapshots.get(v.root);
    if (!before)
      throw new Error(`${v.type}: no snapshot to compare with`);
    await new Promise((resolve) => requestAnimationFrame(() => setTimeout(resolve, 0)));
    return {renders: renders.get(v.root)!.count - before.renders, delta: change(el).delta};
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

  w.__bdd = {viewerOf, arm, stampAll, settle, readProperty, writeProperties, findArea, hitArea, areaInk, snapshot, baselineAll,
    change, rangeChange, rememberRange, rememberedRange, stillness, palette, listen, unlisten, firedCount, resize, restoreSize,
    armEvent, waitArmed, closeMenu, openMenu, menuPoint, addViewer, takeBalloons, saveLayout, loadLayout};
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
 * single repaint under immediate rendering; the cap only matters for a property that paints
 * nothing). The canvas is snapshotted first, so `should have repainted` compares with the state
 * before the change. One roundtrip: `locator.evaluate` waits for the element itself. */
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

export async function canvasChange(page: Page, target: ElementRef): Promise<CanvasChange> {
  const loc = await viewerLocator(page, target);
  return loc.evaluate((el) => (window as any).__bdd.change(el));
}

/** Waits until the canvas differs from the last snapshot, then makes the new state the snapshot. */
export async function expectRepainted(page: Page, target: ElementRef): Promise<void> {
  await expect.poll(async () => (await canvasChange(page, target)).delta, {timeout: 10000, message: `${target.phrase} did not repaint`}).toBeGreaterThan(0);
  await snapshot(page, target);
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
  await snapshot(page, target);
}

/** Pixels in the selection hue against the snapshot: `more`/`less` than before, `some`, or `none`. */
export async function expectHighlight(page: Page, target: ElementRef, compare: 'more' | 'less' | 'some' | 'none'): Promise<void> {
  const wrong = {none: 'shows a selection highlight', some: 'shows no selection highlight',
    more: 'does not show more selection highlight than before', less: 'does not show less selection highlight than before'};
  await expect.poll(async () => {
    const c = await canvasChange(page, target);
    return compare === 'none' ? c.hue === 0 : compare === 'some' ? c.hue > 0 : compare === 'more' ? c.hue > c.hueBefore : c.hue < c.hueBefore;
  }, {timeout: 10000, message: `${target.phrase} ${wrong[compare]}`}).toBe(true);
  await snapshot(page, target);
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

/** The value range (viewport) against the snapshot before the last change. */
export async function expectValueRange(page: Page, target: ElementRef, compare: 'narrower' | 'same' | 'wider'): Promise<void> {
  const same = (a: Range, b: Range) => Math.abs(a.top - b.top) < 1e-6 && Math.abs(a.bottom - b.bottom) < 1e-6;
  let last: RangeChange = {};
  const holds = async (): Promise<boolean | string> => {
    last = await rangeChange(page, target);
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
  await snapshot(page, target);
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
