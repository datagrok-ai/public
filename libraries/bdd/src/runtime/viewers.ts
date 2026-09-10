/* Viewer readers on the Playwright side, over the in-page runtime (`viewer-runtime.ts`): a
   viewer's element, its hit areas, its properties, its size, its layout and the events it fires.
   The pixel claims are in `viewer-pixels.ts`, the context menus in `viewer-menus.ts`, the legend
   and the tooltip in `viewer-legend.ts`; all of them are re-exported here, so `import * as v`
   reaches everything.
   Roundtrips are the cost, not the page: a step is locate + one action, and the baseline a
   "than before" claim compares with is taken inside the same in-page call as the change. */
import {Locator, Page} from '@playwright/test';
import {expect, pollMs} from './patience.js';
import type {ElementRef} from './args.js';
import {typeVerified, withKeys} from './gestures.js';
import {locate} from './locate.js';
import {Balloon, Box, evaluate, installViewerRuntime, Reading} from './viewer-runtime.js';

export * from './viewer-runtime.js';
export * from './viewer-pixels.js';
export * from './viewer-menus.js';
export * from './viewer-legend.js';

/** The single visible element of a viewer phrase (a closed view can leave a zero-size twin). */
export async function viewerLocator(page: Page, target: ElementRef): Promise<Locator> {
  const loc = await locate(page, target);
  return loc.filter({visible: true}).first();
}

/** `locator.evaluate` on the viewer's element with the in-page runtime installed — one roundtrip,
 * and the locator waits for the element itself. */
export async function onViewer<R>(page: Page, target: ElementRef, fn: (el: Element, arg: any) => R | Promise<R>, arg?: unknown): Promise<R> {
  await installViewerRuntime(page);
  const loc = await viewerLocator(page, target);
  return loc.evaluate(fn as (el: SVGElement | HTMLElement, arg: unknown) => R, arg);
}

export function centerOf(box: Box): {x: number; y: number} {
  return {x: box.x + box.width / 2, y: box.y + box.height / 2};
}

/** Waits for the area the way a locator waits for its element: a viewer that renders twice on a
 * change (a category switch relays out after its first paint) reports the area after the second. */
async function awaitArea(page: Page, target: ElementRef, name: string, negate: boolean, beforeChange: boolean): Promise<Box | undefined> {
  const find = (): Promise<{box?: Box; has: string[]}> =>
    onViewer(page, target, (el, [n, b]) => (window as any).__bdd.findArea(el, n, b), [name, beforeChange] as [string, boolean]);
  let found = await find();
  if ((found.box !== undefined) === !negate)
    return found.box;
  const poll = expect.poll(async () => (found = await find()).box !== undefined,
    {timeout: pollMs(5000), message: `${target.phrase} ${negate ? 'still reports' : 'reports no'} "${name}" area; it has: ${found.has.join(', ') || 'none'}`});
  await (negate ? poll.not : poll).toBe(true);
  return found.box;
}

/** `hitArea` in client coordinates — the rectangle the viewer reports for that named region;
 * `beforeChange` takes the canvas baseline in the same call (a click or a double-click follows),
 * after the viewer has finished putting the thing where it is: a network diagram still running
 * its physics moves the node between the read and the click. */
export async function hitArea(page: Page, target: ElementRef, name: string, beforeChange = false): Promise<Box> {
  if (beforeChange)
    await onViewer(page, target, (el) => (window as any).__bdd.settle(el, 300), undefined).catch(() => undefined);
  return (await awaitArea(page, target, name, false, beforeChange))!;
}

export async function expectHasArea(page: Page, target: ElementRef, name: string, negate = false): Promise<void> {
  await awaitArea(page, target, name, negate, false);
}

/** Every hit area the viewer reports right now, in client coordinates — for a step that reasons
 * over them together (the order of the bars, a spot no bar covers). */
export function hitAreas(page: Page, target: ElementRef): Promise<Record<string, Box>> {
  return onViewer(page, target, (el) => (window as any).__bdd.areas(el), undefined);
}

export async function addViewer(page: Page, type: string): Promise<void> {
  await evaluate(page, (t) => { (window as any).__bdd.addViewer(t); }, type);
  await page.locator(`[name="viewer-${type.replace(/\s+/g, '-')}" i]`).filter({visible: true}).first().waitFor();
}

/** Sets properties by caption in one go: one settle for the group (the sets coalesce into a
 * single repaint under immediate rendering; the cap only matters for a viewer without the pending
 * signal). The canvas is snapshotted first, so `should have repainted` compares with the state
 * before the change. */
export async function setProperties(page: Page, target: ElementRef, entries: [string, string][], capMs = 300): Promise<void> {
  await onViewer(page, target, (el, [e, cap]) => (window as any).__bdd.writeProperties(el, e, cap), [entries, capMs] as [[string, string][], number]);
}

export function readProperty(page: Page, target: ElementRef, caption: string, expected = ''): Promise<string> {
  return onViewer(page, target, (el, [c, x]) => (window as any).__bdd.readProperty(el, c, x), [caption, expected] as [string, string]);
}

export async function expectProperty(page: Page, target: ElementRef, caption: string, value: string, negate = false): Promise<void> {
  const poll = expect.poll(() => readProperty(page, target, caption, value), {timeout: pollMs(5000), message: `"${caption}" of ${target.phrase}`});
  await (negate ? poll.not : poll).toBe(value.replace(/\\n/g, '\n'));
}

export function snapshot(page: Page, target: ElementRef): Promise<void> {
  return onViewer(page, target, (el) => { (window as any).__bdd.snapshot(el); }, undefined);
}

/** Every viewer of every open table view gets its baseline — before a change that is not one
 * viewer's own (a filter, a selection, a column's colors). */
export function baselineAll(page: Page): Promise<void> {
  return evaluate(page, () => { (window as any).__bdd.baselineAll(); }, undefined);
}

/** Waits until no viewer of any open table view has a refresh or a repaint pending — after a
 * change that reaches them all, so the next step's baseline is the state after it. */
export function settleAll(page: Page): Promise<void> {
  return evaluate(page, async () => {
    const b = (window as any).__bdd;
    for (const view of Array.from((window as any).grok.shell.tableViews ?? []) as any[]) {
      for (const v of Array.from(view.viewers ?? []) as any[])
        await b.quiet(v);
    }
  }, undefined);
}

/** Resolves once the viewer has nothing pending (`capMs` for a viewer without the signal). */
export function settle(page: Page, target: ElementRef, capMs = 300): Promise<number> {
  return onViewer(page, target, (el, ms) => (window as any).__bdd.settle(el, ms), capMs);
}

/** The table the viewer draws (`viewer.dataFrame`), not the property it was asked to bind. */
export async function expectBoundTable(page: Page, target: ElementRef, name: string): Promise<void> {
  await expect.poll(() => onViewer(page, target, (el) => (window as any).__bdd.tableOf(el), undefined),
    {timeout: pollMs(5000), message: `the table ${target.phrase} is bound to`}).toBe(name);
}

/** A reading of the viewer as it is now; a name the viewer does not report fails naming the
 * readings it does. */
export async function readValue(page: Page, target: ElementRef, name: string): Promise<unknown> {
  const r: Reading = await onViewer(page, target, (el, n) => (window as any).__bdd.valueChange(el, n), name);
  if (r.now === undefined || r.now === null)
    throw new Error(`${target.phrase} has no "${name}" reading; it reports: ${r.has.join(', ') || 'no readings'}`);
  return r.now;
}

export type ReadingCompare = 'equal' | 'lower' | 'higher' | 'differ' | 'same';

/** A reading of the viewer (`getWidgetStatus().values`: "rows shown", "bars", "scene signature")
 * equals a value, is lower/higher than at the snapshot before the last change, differs from it, or
 * is the same — the negative read once the viewer is quiet, and read once. */
export async function expectReading(page: Page, target: ElementRef, name: string, compare: ReadingCompare, value?: number): Promise<void> {
  let last: Reading = {has: []};
  const holds = async (): Promise<boolean | string> => {
    last = await onViewer(page, target, (el, [n, q]) => {
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

export function rememberReading(page: Page, target: ElementRef, name: string): Promise<void> {
  return onViewer(page, target, (el, n) => { (window as any).__bdd.rememberValue(el, n); }, name);
}

export async function expectRememberedReading(page: Page, target: ElementRef, name: string, not = false): Promise<void> {
  let last: Reading = {has: []};
  const holds = async (): Promise<boolean | string> => {
    last = await onViewer(page, target, (el, n) => (window as any).__bdd.rememberedValue(el, n), name);
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

/** The area's rectangle against the snapshot's: taller or wider than before. */
export async function expectAreaGrew(page: Page, target: ElementRef, area: string, dimension: 'taller' | 'wider' | 'shorter' | 'narrower'): Promise<void> {
  let last: {before?: Box; now?: Box; has: string[]} = {has: []};
  const holds = async (): Promise<boolean | string> => {
    last = await onViewer(page, target, (el, a) => (window as any).__bdd.areaRectChange(el, a), area);
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

export async function expectAreaSize(page: Page, target: ElementRef, area: string, dimension: 'tall' | 'wide', min: number): Promise<void> {
  const box = await hitArea(page, target, area);
  const size = dimension === 'tall' ? box.height : box.width;
  expect(size, `the "${area}" area of ${target.phrase} is ${Math.round(size)} px ${dimension}, not at least ${min}`).toBeGreaterThanOrEqual(min);
}

// --- events ------------------------------------------------------------------------------------------

export function listenFor(page: Page, target: ElementRef, event: string): Promise<void> {
  return onViewer(page, target, (el, e) => (window as any).__bdd.listen(el, e), event);
}

export async function expectFired(page: Page, target: ElementRef, event: string): Promise<void> {
  await expect.poll(() => onViewer(page, target, (el, e) => (window as any).__bdd.firedCount(el, e), event),
    {timeout: pollMs(5000), message: `"${event}" did not fire on ${target.phrase} (listen for it before the gesture)`}).toBeGreaterThan(0);
  await onViewer(page, target, (el, e) => { const b = (window as any).__bdd; b.unlisten(b.viewerOf(el), e); }, event);
}

/** Not fired so far; the subscription stays for a later "should have fired". */
export async function expectNotFired(page: Page, target: ElementRef, event: string): Promise<void> {
  const count: number = await onViewer(page, target, (el, e) => (window as any).__bdd.firedCount(el, e), event);
  expect(count, count < 0 ? `"${event}" is not listened for on ${target.phrase}` : `"${event}" fired ${count} time(s) on ${target.phrase}`).toBe(0);
}

/** The balloons (info, warning, error) shown since the last read; reading clears them. */
export function takeBalloons(page: Page): Promise<Balloon[]> {
  return evaluate(page, () => (window as any).__bdd.takeBalloons(), undefined);
}

// --- size and layout -----------------------------------------------------------------------------------

export async function resize(page: Page, target: ElementRef, width: number | null, height: number | null): Promise<void> {
  await onViewer(page, target, (el, [w, h]) => (window as any).__bdd.resize(el, w, h, 500), [width, height] as [number | null, number | null]);
}

export async function restoreSize(page: Page, target: ElementRef): Promise<void> {
  await onViewer(page, target, (el) => (window as any).__bdd.restoreSize(el, 500), undefined);
}

export function saveLayout(page: Page): Promise<void> {
  return evaluate(page, () => { (window as any).__bdd.saveLayout(); }, undefined);
}

/** Saves the layout through the server and returns its id; the caller registers the deletion. */
export function saveLayoutToServer(page: Page): Promise<string> {
  return evaluate(page, () => (window as any).__bdd.saveLayoutToServer(), undefined);
}

export function deleteLayout(page: Page, id: string): Promise<void> {
  return page.evaluate((i) => (window as any).__bdd.deleteLayout(i), id);
}

export function loadLayout(page: Page): Promise<void> {
  return evaluate(page, () => (window as any).__bdd.loadLayout(), undefined);
}

// --- area gestures ------------------------------------------------------------------------------------

export {withKeys};

/** A plain drag from the centre of one hit area to the centre of another (a column header to a
 * new place, a range handle to a bin), with keys held; the baseline is taken before the drag. */
export async function dragArea(page: Page, target: ElementRef, from: string, to: string, keys: string[] = []): Promise<void> {
  const a = centerOf(await hitArea(page, target, from, true));
  const b = centerOf(await hitArea(page, target, to));
  await withKeys(page, keys, async () => {
    await page.mouse.move(a.x, a.y);
    await page.mouse.down();
    await page.mouse.move(b.x, b.y, {steps: 3});
    await page.mouse.up();
  });
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

/** A drag across the inner 80% of an area with keys held: Shift selects, Control+Shift removes
 * from the selection, Alt zooms. */
export async function dragBoxOverArea(page: Page, target: ElementRef, area: string, keys: string[]): Promise<void> {
  const b = await hitArea(page, target, area, true);
  await withKeys(page, keys, async () => {
    await page.mouse.move(b.x + b.width * 0.1, b.y + b.height * 0.1);
    await page.mouse.down();
    await page.mouse.move(b.x + b.width * 0.9, b.y + b.height * 0.9, {steps: 3});
    await page.mouse.up();
  });
}

/** A click on the area, the text typed over what the editor there holds, Enter. The click must have
 * put the focus into an editor inside the viewer — a histogram's range input took the click and not
 * the focus once in twenty runs, and the text then opened a cell editor on the grid, unseen — so it
 * is repeated, at the area's current place, until one did. The editor is pinned by a mark of its
 * own rather than by `:focus`: the focus can leave it while the text is being read back. */
export async function typeIntoArea(page: Page, target: ElementRef, area: string, text: string): Promise<void> {
  const loc = await viewerLocator(page, target);
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
