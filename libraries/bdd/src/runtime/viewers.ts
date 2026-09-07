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
 * before the change, `ink` the painted pixels now, `inkBefore` those of the snapshot. */
export interface CanvasChange {
  delta: number;
  ink: number;
  inkBefore: number;
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
  const snapshots = new WeakMap<Element, {colors: Map<number, number>; ink: number}>();
  const sizes = new WeakMap<Element, {width: string; height: string}>();
  const listeners = new WeakMap<Element, Record<string, {count: number; sub: any}>>();
  const armed: Record<string, Promise<unknown>> = {};
  let tokens = 0;
  const norm = (s: unknown) => String(s ?? '').toLowerCase().replace(/[^a-z0-9]/g, '');

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
    for (let i = 0; i < data.length; i += 4) {
      const key = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
      colors.set(key, (colors.get(key) ?? 0) + 1);
      if (data[i + 3] !== 0 && !(data[i] >= 250 && data[i + 1] >= 250 && data[i + 2] >= 250))
        ink++;
    }
    return {colors, ink};
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
    if (!p)
      throw new Error(`${v.type} has no "${caption}" property; it has: ${props.map((x) => x.caption ?? x.name).join(', ')}`);
    return p;
  };
  const convert = (p: any, text: string): unknown => {
    const s = text.replace(/\\n/g, '\n');
    const type = String(p.propertyType ?? '');
    if (type === 'bool')
      return /^(true|yes|on|checked|1)$/i.test(s);
    if (type === 'int' || type === 'double' || type === 'num' || type === 'bigint') {
      if (s === '')
        return null;
      const hex = /^#?([0-9a-f]{6})$/i.exec(s);
      if (hex && type === 'int')
        return (0xFF000000 | parseInt(hex[1], 16)) >>> 0;
      return Number(s);
    }
    if (type === 'string_list' || type === 'list' || type === 'column_list')
      return s === '' ? [] : s.split(/\s*,\s*/);
    return s;
  };
  const readProperty = (el: Element, caption: string): string => {
    const v = viewerOf(el);
    const value = v.props[property(v, caption).name];
    if (value == null)
      return '';
    return Array.isArray(value) ? value.join(', ') : String(value);
  };
  const snapshot = (el: Element): number => {
    const v = viewerOf(el);
    const shot = histogram(canvasOf(v));
    snapshots.set(v.root, shot);
    return shot.ink;
  };
  // the baseline "should have repainted" compares with: the canvas before the change
  const baseline = (el: Element): void => {
    try {
      snapshot(el);
    } catch { /* a viewer without a canvas */ }
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
    const cv = canvasOf(v).getBoundingClientRect();
    return {box: {x: cv.x + r.x, y: cv.y + r.y, width: r.width, height: r.height}, has};
  };
  const hitArea = (el: Element, name: string, beforeChange = false): Box => {
    const found = findArea(el, name, beforeChange);
    if (!found.box)
      throw new Error(`${viewerOf(el).type} has no "${name}" area right now; it has: ${found.has.join(', ') || 'none'}`);
    return found.box;
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
    return {delta, ink: now.ink, inkBefore: before.ink};
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

  w.__bdd = {viewerOf, arm, stampAll, settle, readProperty, writeProperties, findArea, hitArea, snapshot, change,
    listen, unlisten, firedCount, resize, restoreSize, armEvent, waitArmed, closeMenu, openMenu, menuPoint, addViewer};
  stampAll();
  grok.events.onViewerAdded.subscribe((a: any) => arm(a?.args?.viewer));
  grok.events.onViewerClosed.subscribe((a: any) => a?.args?.viewer && forget(a.args.viewer));
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

export async function readProperty(page: Page, target: ElementRef, caption: string): Promise<string> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  return loc.evaluate((el, c) => (window as any).__bdd.readProperty(el, c), caption);
}

export async function expectProperty(page: Page, target: ElementRef, caption: string, value: string, negate = false): Promise<void> {
  const poll = expect.poll(() => readProperty(page, target, caption), {timeout: 5000, message: `"${caption}" of ${target.phrase}`});
  await (negate ? poll.not : poll).toBe(value.replace(/\\n/g, '\n'));
}

export async function snapshot(page: Page, target: ElementRef): Promise<number> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  return loc.evaluate((el) => (window as any).__bdd.snapshot(el));
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
