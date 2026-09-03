import {expect, Page} from '@playwright/test';
import * as v from '../../helpers/viewers';
import {isLocalBootNoise} from '../../spec-login';

declare const grok: any;

export const SP_TYPE = 'Scatter plot';

export interface Rect {x: number; y: number; width: number; height: number}

export const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Failed to connect to Claude runtime/.test(text) ||
  /powerPreference option is currently ignored/.test(text) ||
  /willReadFrequently/.test(text) || isLocalBootNoise(text);

export interface ErrorTracker { count: () => number; all: () => string[]; }

export function trackErrors(page: Page, isBenign: (text: string) => boolean = isAmbientError): ErrorTracker {
  const errors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenign(String(e))) errors.push(String(e)); });
  page.on('console', (m) => { if (m.type() === 'error' && !isBenign(m.text())) errors.push(m.text()); });
  return {count: () => errors.length, all: () => errors};
}

export const canvasRect = (page: Page): Promise<Rect> => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const r = sp.root.querySelector('canvas[name="canvas"]').getBoundingClientRect();
  return {x: r.x, y: r.y, width: r.width, height: r.height};
});

export async function addScatterPlot(page: Page): Promise<void> {
  await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
  await page.waitForFunction(() => {
    const sp = grok.shell.tv?.viewers?.find((x: any) => x.type === 'Scatter plot') as any;
    const c = sp?.root?.querySelector('canvas[name="canvas"]');
    return !!c && c.getBoundingClientRect().width > 0 && c.getBoundingClientRect().height > 0;
  }, null, {timeout: 20_000});
  await v.installEventWaits(page);
  // installs the render stamp the later quiet-waits read; a plot painted before this point costs the cap
  await v.waitForViewerQuiet(page, SP_TYPE, {gapMs: 200, capMs: 600});
}

export const readProp = (page: Page, name: string) => page.evaluate((n: string) => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return sp ? sp.props[n] ?? null : null;
}, name);

export const propIs = (page: Page, propName: string, value: unknown, timeoutMs = 2500) =>
  v.pollValue(() => readProp(page, propName), (x) => x === value, timeoutMs, 50);

export const pickOnViewer = (page: Page, role: string, column: string) =>
  v.pickColumnViaSelectorTrusted(page, {role, columnName: column});

const waitBackdrop = (page: Page, timeout = 5000) =>
  page.waitForFunction(() => !!document.querySelector('.d4-column-selector-backdrop'), null, {timeout})
    .then(() => true).catch(() => false);

// The popup's search box appears once a letter reaches the combobox root, which the platform
// focuses a tick after the popup opens; a letter can land before or twice, so the box's text is
// replaced wholesale rather than appended to.
async function openSearchBox(page: Page): Promise<void> {
  const searchFocused = () => page.evaluate(() =>
    document.activeElement?.classList.contains('d4-column-selector-search-input') ?? false);
  await v.pollValue(() => page.evaluate(() =>
    document.activeElement?.classList.contains('d4-column-selector') ?? false), (ok) => ok, 500, 25);
  for (let i = 0; i < 2 && !await searchFocused(); i++) {
    await page.keyboard.press('a');
    await v.pollValue(searchFocused, (ok) => ok, 1000, 25);
  }
  if (!await searchFocused()) throw new Error('the column popup search box did not open');
  await page.keyboard.press('Control+A');
}

async function commitColumn(page: Page, column: string): Promise<void> {
  await openSearchBox(page);
  await page.keyboard.type(column.toLowerCase());
  await page.keyboard.press('Enter');
}

// The popup's empty row is the first row of a canvas grid, so there is nothing to click by
// name; an empty search text committed with Enter resolves to that row (ColumnGrid.currentColumnName).
async function commitEmptyColumn(page: Page): Promise<void> {
  await openSearchBox(page);
  await page.keyboard.press('Backspace');
  await page.keyboard.press('Enter');
}

export async function clearOnViewer(page: Page, role: string, propName = `${role}ColumnName`): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.move(r.x + r.width / 2, r.y + r.height / 2);
  const pt = await v.pollValue(() => page.evaluate((rl: string) => {
    const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
    const sel = sp.root.querySelector(`[name="div-column-combobox-${rl}"]`) as HTMLElement | null;
    const el = sel?.querySelector('.d4-column-selector-column') ?? sel;
    const b = el?.getBoundingClientRect();
    return b && b.width > 0 && b.height > 0 ? {x: b.x + b.width / 2, y: b.y + b.height / 2} : null;
  }, role), (p) => p !== null, 1000, 50);
  if (!pt) throw new Error(`on-viewer ${role} selector has no clickable text`);
  await page.mouse.click(pt.x, pt.y);
  if (!await waitBackdrop(page)) throw new Error(`${role} column popup did not open`);
  await commitEmptyColumn(page);
  const left = await v.pollValue(() => readProp(page, propName), (p) => p === '' || p === null, 2500, 50);
  if (left !== '' && left !== null) throw new Error(`${propName} was left as ${left} by the ${role} selector`);
}

export const openSettings = (page: Page) => v.openViewerSettings(page, SP_TYPE);

const editorSized = (page: Page, selector: string) => page.evaluate((sel: string) => {
  const el = document.querySelector(sel) as HTMLElement | null;
  if (!el || !el.offsetParent) return false;
  const b = el.getBoundingClientRect();
  return b.width > 0 && b.height > 0;
}, selector);

export async function revealPropEditor(page: Page, editorSelector: string, category: string): Promise<void> {
  const ready = () => editorSized(page, editorSelector);
  for (let i = 0; i < 4; i++) {
    if (await ready()) return;
    const header = page.locator(`[name="prop-category-${category}"]`);
    if (await header.count() > 0 && await header.isVisible()) await header.click();
    if (await v.pollValue(ready, (ok) => ok, 800, 50)) return;
  }
  // the named category is a hint: property rows move between categories, so expanding
  // every header reaches the row wherever it actually lives
  for (const h of await page.locator('[name^="prop-category-"]').all())
    if (await h.isVisible()) await h.click();
  if (await v.pollValue(ready, (ok) => ok, 2000, 50)) return;
  throw new Error(`property editor ${editorSelector} never became reachable`);
}

export const rowProp = (rowName: string) =>
  rowName.replace(/^prop-/, '').replace(/-([a-z])/g, (_, c: string) => c.toUpperCase());

export const propCellText = (page: Page, viewCell: string) => page.evaluate((n: string) =>
  (document.querySelector(`[name="${n}"]`)?.textContent ?? '').trim(), viewCell);

export const rowOpacity = (page: Page, rowName: string) => page.evaluate((n: string) =>
  (document.querySelector(`[name="${n}"]`) as HTMLElement | null)?.style.opacity ?? null, rowName);

export async function setChoiceProp(
  page: Page, rowName: string, category: string, value: string, propName = rowProp(rowName),
): Promise<void> {
  const viewCell = rowName.replace(/^prop-/, 'prop-view-');
  await openSettings(page);
  await revealPropEditor(page, `[name="${viewCell}"]`, category);
  const cell = page.locator(`[name="${viewCell}"]`);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  const select = page.locator(`[name="${rowName}"] select.property-grid-item-editor-spinner`);
  await select.waitFor({state: 'visible', timeout: 4000});
  await select.selectOption(value);
  const got = await v.pollValue(() => readProp(page, propName),
    (x) => String(x).toLowerCase() === value.toLowerCase(), 2500, 50);
  expect(String(got).toLowerCase()).toBe(value.toLowerCase());
  await v.waitForViewerRendered(page, SP_TYPE, 300);
}

export async function setCheckboxProp(
  page: Page, rowName: string, category: string, value: boolean, propName = rowProp(rowName),
): Promise<void> {
  await openSettings(page);
  const box = `[name="${rowName}"] input.property-grid-item-editor-checkbox`;
  await revealPropEditor(page, box, category);
  for (let i = 0; i < 3 && await readProp(page, propName) !== value; i++) {
    const locator = page.locator(box);
    await locator.scrollIntoViewIfNeeded();
    await locator.click();
    await propIs(page, propName, value, 2000);
  }
  if (await readProp(page, propName) !== value)
    throw new Error(`${propName} did not reach ${value} from the ${rowName} row`);
  await v.waitForViewerRendered(page, SP_TYPE, 300);
}

export async function setNumericProp(
  page: Page, rowName: string, category: string, value: number, propName = rowProp(rowName),
): Promise<void> {
  await openSettings(page);
  const selector = `[name="${rowName}"] input.property-grid-slider-textbox, [name="${rowName}"] input`;
  await revealPropEditor(page, selector, category);
  const locator = page.locator(selector).first();
  await locator.scrollIntoViewIfNeeded();
  await locator.click();
  await v.pollValue(() => page.evaluate(() => document.activeElement instanceof HTMLInputElement), (f) => f, 1000, 25);
  await page.keyboard.press('Control+A');
  await page.keyboard.type(String(value));
  await page.keyboard.press('Enter');
  if (await propIs(page, propName, value, 2500) !== value)
    throw new Error(`${propName} did not reach ${value} from the ${rowName} row`);
  await v.waitForViewerRendered(page, SP_TYPE, 300);
}

export async function setTextProp(
  page: Page, rowName: string, category: string, value: string, propName = rowProp(rowName),
): Promise<void> {
  const viewCell = rowName.replace(/^prop-/, 'prop-view-');
  await openSettings(page);
  await revealPropEditor(page, `[name="${viewCell}"]`, category);
  const cell = page.locator(`[name="${viewCell}"]`);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  const editor = page.locator(
    `[name="${rowName}"] input.property-grid-item-editor-textbox, ` +
    `[name="${rowName}"] input.property-grid-ellipsis-editor-input`).first();
  await editor.waitFor({state: 'visible', timeout: 2000});
  await editor.click();
  await page.keyboard.press('Control+A');
  if (value === '') await page.keyboard.press('Delete');
  else await page.keyboard.type(value);
  await page.keyboard.press('Enter');
  if (await propIs(page, propName, value, 2500) !== value)
    throw new Error(`${propName} did not reach "${value}" from the ${rowName} row`);
}

async function openPanelColumnPopup(page: Page, rowName: string, comboName: string, category: string): Promise<void> {
  await openSettings(page);
  const combo = `[name="${rowName}"] [name="${comboName}"]`;
  await revealPropEditor(page, combo, category);
  const sel = page.locator(combo);
  await sel.scrollIntoViewIfNeeded();
  await sel.click();
  if (!await waitBackdrop(page)) throw new Error(`${comboName} popup did not open`);
}

export async function pickPanelColumn(
  page: Page, rowName: string, comboName: string, category: string, column: string, propName = rowProp(rowName) + 'ColumnName',
): Promise<void> {
  await openPanelColumnPopup(page, rowName, comboName, category);
  await commitColumn(page, column);
  if (await propIs(page, propName, column, 2500) !== column)
    throw new Error(`${propName} did not reach ${column} from the ${rowName} row`);
  await v.waitForViewerRendered(page, SP_TYPE, 300);
}

export async function clearPanelColumn(
  page: Page, rowName: string, comboName: string, category: string, propName = rowProp(rowName) + 'ColumnName',
): Promise<void> {
  await openPanelColumnPopup(page, rowName, comboName, category);
  await commitEmptyColumn(page);
  const left = await v.pollValue(() => readProp(page, propName), (p) => p === '' || p === null, 2500, 50);
  if (left !== '' && left !== null) throw new Error(`${propName} was left as ${left} by the ${rowName} row`);
  await v.waitForViewerRendered(page, SP_TYPE, 300);
}

export const filterCount = (page: Page) =>
  page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount as number);

/** The filtered row count once it has left `from` and gone quiet. */
export const filterMoved = (page: Page, from: number, capMs = 5000): Promise<number> =>
  page.evaluate(({from, cap}) =>
    (window as any).__moved(() => grok.shell.tv.dataFrame.filter.trueCount, from, cap), {from, cap: capMs});

/** A hold proving the filtered row count does NOT move: a capped poll for a change. */
export const filterHeld = (page: Page, capMs = 1000): Promise<number> =>
  page.evaluate((cap) => {
    const read = () => grok.shell.tv.dataFrame.filter.trueCount as number;
    const first = read();
    return (window as any).__poll(read, (c: number) => c !== first, cap, 50);
  }, capMs);

export const menuLeafNames = (page: Page) => page.evaluate(() =>
  [...document.querySelectorAll('.d4-menu-popup [name]')]
    .map((e) => e.getAttribute('name')!).filter((n) => !!n));

export async function openPlotContextMenu(page: Page): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.click(r.x + r.width / 2, r.y + r.height / 2, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});
  await v.pollValue(() => menuLeafNames(page), (names) => names.includes('div-Properties...'), 3000, 50);
}

export async function dismissMenu(page: Page): Promise<void> {
  await page.keyboard.press('Escape');
  await v.pollValue(() => page.locator('.d4-menu-popup:visible').count(), (n) => n === 0, 1000, 50);
}

/** Drives an open popup along `names` (item name attributes, groups first, the leaf last) and clicks the leaf. */
export async function navigatePopup(page: Page, names: string[]): Promise<void> {
  await page.evaluate(async (path: string[]) => {
    const w = window as any;
    const find = (n: string) => {
      const popup = [...document.querySelectorAll('.d4-menu-popup')].pop();
      const el = popup?.querySelector(`[name="${n}"]`) as HTMLElement | null;
      return el ? (el.closest('.d4-menu-item') as HTMLElement | null) ?? el : null;
    };
    for (let i = 0; i < path.length; i++) {
      const item = await w.__poll(() => find(path[i]), (el: Element | null) => el !== null, 5000, 50);
      if (!item) {
        const popup = [...document.querySelectorAll('.d4-menu-popup')].pop();
        throw new Error(`menu: "${path[i]}" not found; visible: ` +
          [...(popup?.querySelectorAll('[name]') ?? [])].map((x) => x.getAttribute('name')).join(' | '));
      }
      if (i === path.length - 1) { item.click(); return; }
      const b = item.getBoundingClientRect();
      for (const type of ['mouseover', 'mousemove'])
        item.dispatchEvent(new MouseEvent(type, {bubbles: true, clientX: b.x + 5, clientY: b.y + 5}));
    }
  }, names);
}

export async function clickContextMenuLeaf(page: Page, names: string[]): Promise<void> {
  await openPlotContextMenu(page);
  await navigatePopup(page, names);
  await v.waitForViewerRendered(page, SP_TYPE, 1000);
  await v.pollValue(() => page.locator('.d4-menu-popup:visible').count(), (n) => n === 0, 1000, 50);
}

export async function resetView(page: Page): Promise<void> {
  await openPlotContextMenu(page);
  await page.locator('.d4-menu-popup [name="div-Reset-View"]').last().click();
  await v.waitForViewerRendered(page, SP_TYPE, 800);
  await dismissMenu(page);
}

export const viewport = (page: Page): Promise<Rect> => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const vp = sp.viewport;
  return {x: vp.x, y: vp.y, width: vp.width, height: vp.height};
});

/** sp.viewport once it has left `from` and gone quiet. */
export const viewportMoved = (page: Page, from: Rect, capMs = 4000): Promise<Rect> =>
  page.evaluate(async ({from, cap}) => {
    const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
    const read = () => { const vp = sp.viewport; return `${vp.x}|${vp.y}|${vp.width}|${vp.height}`; };
    await (window as any).__moved(read, `${from.x}|${from.y}|${from.width}|${from.height}`, cap);
    const vp = sp.viewport;
    return {x: vp.x, y: vp.y, width: vp.width, height: vp.height};
  }, {from, cap: capMs});

export async function parkPointer(page: Page, capMs = 1500): Promise<void> {
  await page.mouse.move(4, 4);
  await v.waitForViewerQuiet(page, SP_TYPE, {gapMs: 200, capMs});
}

/** Sampled non-white pixel count on one scatter plot canvas layer. */
export const canvasInk = (page: Page, layer: 'canvas' | 'overlay') => page.evaluate((name: string) => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const c = sp?.root.querySelector(`canvas[name="${name}"]`) as HTMLCanvasElement | null;
  const ctx = c?.getContext('2d');
  if (!c || !ctx) return -1;
  let data: Uint8ClampedArray;
  try { data = ctx.getImageData(0, 0, c.width, c.height).data; } catch (_) { return -1; }
  let n = 0;
  for (let k = 0; k < data.length; k += 16)
    if (data[k + 3] !== 0 && !(data[k] >= 250 && data[k + 1] >= 250 && data[k + 2] >= 250)) n++;
  return n;
}, layer);

/**
 * The ink of a layer once it has settled: with `from`, after it has first moved away from that
 * value by more than `tolerance`; without, after two reads agree within `tolerance`.
 */
export async function settledInk(
  page: Page, layer: 'canvas' | 'overlay', tolerance: number, from?: number, capMs = 4000,
): Promise<number> {
  await parkPointer(page);
  if (from !== undefined)
    await v.pollValue(() => canvasInk(page, layer), (cur) => cur >= 0 && Math.abs(cur - from) > tolerance, capMs, 100);
  return v.pollStable(() => canvasInk(page, layer), (a, b) => a >= 0 && Math.abs(a - b) <= tolerance, capMs, 150);
}

export const selectionCount = (page: Page) =>
  page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount as number);

export const selectionMoved = (page: Page, from: number, capMs = 4000): Promise<number> =>
  page.evaluate(({from, cap}) =>
    (window as any).__moved(() => grok.shell.tv.dataFrame.selection.trueCount, from, cap), {from, cap: capMs});

export const selectionHeld = (page: Page, capMs = 600): Promise<number> =>
  page.evaluate((cap) => {
    const read = () => grok.shell.tv.dataFrame.selection.trueCount as number;
    const first = read();
    return (window as any).__poll(read, (c: number) => c !== first, cap, 50);
  }, capMs);

export interface Frac {fx: number; fy: number}

export const at = (r: Rect, p: Frac) => ({x: r.x + r.width * p.fx, y: r.y + r.height * p.fy});

export async function dragCanvas(page: Page, from: Frac, to: Frac, mods: string[] = []): Promise<void> {
  const r = await canvasRect(page);
  const p1 = at(r, from);
  const p2 = at(r, to);
  const shown1 = await v.armEvent(page, 'grok.events.onTooltipShown', 100);
  await page.mouse.move(p1.x, p1.y);
  await shown1();
  for (const m of mods) await page.keyboard.down(m);
  await page.mouse.down();
  await page.mouse.move((p1.x + p2.x) / 2, (p1.y + p2.y) / 2, {steps: 8});
  const shown2 = await v.armEvent(page, 'grok.events.onTooltipShown', 150);
  await page.mouse.move(p2.x, p2.y, {steps: 8});
  await shown2();
  await page.mouse.up();
  for (const m of [...mods].reverse()) await page.keyboard.up(m);
  await v.waitForViewerRendered(page, SP_TYPE, 300);
}

export async function clickCanvas(page: Page, p: Frac): Promise<void> {
  const r = await canvasRect(page);
  const pt = at(r, p);
  await page.mouse.click(pt.x, pt.y);
  await v.waitForViewerRendered(page, SP_TYPE, 250);
}

// A Filters viewer left open with every card removed never shows a card again, and
// openFilterPanel waits on one; close that empty panel first so the group rebuilds.
export async function openFilterPanel(page: Page): Promise<void> {
  const emptyPanel = await page.evaluate(() => {
    const f = grok.shell.tv.viewers.find((x: any) => x.type === 'Filters');
    if (!f || f.root.querySelector('.d4-filter')) return false;
    f.close();
    return true;
  });
  if (emptyPanel)
    await v.pollValue(() => page.evaluate(() =>
      !grok.shell.tv.viewers.find((x: any) => x.type === 'Filters')), (gone) => gone, 1000, 50);
  await v.openFilterPanel(page);
}

export async function closeFilterPanel(page: Page): Promise<void> {
  await page.evaluate(() => grok.shell.tv.viewers.find((x: any) => x.type === 'Filters')?.close());
  await v.pollValue(() => page.evaluate(() =>
    !grok.shell.tv.viewers.find((x: any) => x.type === 'Filters')), (gone) => gone, 1000, 50);
}

/** Opens a second table view on `path` next to the current one (openTable would close it). */
export async function addTableView(page: Page, path: string, probeColumn: string): Promise<void> {
  await page.evaluate(async ({p, col}) => {
    const df = await (window as any).__readCsv(p);
    grok.shell.addTableView(df);
    await (window as any).__poll(() => !!grok.shell.tv?.dataFrame?.col(col), (ok: boolean) => ok, 10_000, 50);
  }, {p: path, col: probeColumn});
}
