import {Page} from '@playwright/test';
import * as v from '../../helpers/viewers';

export const PIVOT = '[name="viewer-Pivot-table"]';
export const INNER_CANVAS = `${PIVOT} .grok-pivot-grid [name="viewer-Grid"] canvas[name="canvas"]`;
export const DEMOG = 'System:DemoFiles/demog.csv';

export async function openPivot(page: Page, path = DEMOG): Promise<void> {
  await v.openTable(page, {path, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'pivot-table', 'Pivot-table', 15000);
  await page.locator(`${PIVOT} .grok-pivot-column-tags-title[d4-name="Group by"]`).waitFor({timeout: 15000});
  // installs the render stamp so later waits see renders that happened before they were called
  await v.waitForViewerRendered(page, 'Pivot table', 50);
}

export async function rowChips(page: Page, rowTitle: string): Promise<string[]> {
  return page.evaluate((title) => {
    const root = document.querySelector('[name="viewer-Pivot-table"]');
    const panel = Array.from(root?.querySelectorAll('.grok-pivot-column-panel') ?? [])
      .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === title);
    if (!panel) return [];
    return Array.from(panel.querySelectorAll('.d4-tag'))
      .map((t) => (t.querySelector('span')?.textContent ?? t.textContent ?? '').trim());
  }, rowTitle);
}

export async function pivotProps(page: Page) {
  return page.evaluate(() => {
    const pv = Array.from((window as any).grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
    return {
      groupBy: pv.props.groupByColumnNames as string[],
      pivot: pv.props.pivotColumnNames as string[],
      agg: pv.props.aggregateColumnNames as string[],
      aggTypes: pv.props.aggregateAggTypes as string[],
      rowSource: pv.props.rowSource as string,
      filteringEnabled: pv.props.filteringEnabled as boolean,
    };
  });
}

export async function ensurePivotPlusClickable(page: Page, plusName: string): Promise<void> {
  const clear = () => page.evaluate((s) => {
    const plus = document.querySelector(s) as HTMLElement | null;
    if (!plus) return false;
    const r = plus.getBoundingClientRect();
    if (r.width === 0 || r.height === 0) return false;
    const top = document.elementFromPoint(r.x + r.width / 2, r.y + r.height / 2);
    return !!top && plus.contains(top);
  }, `${PIVOT} [name="${plusName}"]`);
  if (await clear()) return;
  await page.keyboard.press('Escape');
  if (!await v.pollValue(clear, (c) => c, 4000, 100))
    throw new Error(`pivot ${plusName} + icon stayed obscured (overlay never cleared)`);
}

const SEARCH_INPUT = 'input.d4-column-selector-search-input';

export async function addColumnViaPlus(page: Page, plusName: string, columnName: string): Promise<void> {
  await ensurePivotPlusClickable(page, plusName);
  const plus = `${PIVOT} [name="${plusName}"]`;
  await page.locator(plus).click();
  await page.waitForSelector('.d4-column-selector-backdrop', {timeout: 6000});
  // Enter commits the row under the pointer over the typed text (column_combo_box.dart:312-323),
  // and the click left the pointer over the popup's first rows
  await page.mouse.move(0, 0);
  // the popup focuses the plus icon on a timer, and the keydown that creates the search box is on that icon
  await v.pollValue(() => page.evaluate((s) => {
    const origin = document.querySelector(s)?.firstElementChild as HTMLElement | null;
    if (!origin) return false;
    if (document.activeElement !== origin) origin.focus();
    return document.activeElement === origin;
  }, plus), (focused) => focused, 2000, 50);
  const searchFocused = () => page.evaluate((sel) => document.activeElement?.matches(sel) === true, SEARCH_INPUT);
  await page.keyboard.press(columnName[0]);
  if (!await v.pollValue(searchFocused, (f) => f, 1500, 50)) {
    await page.keyboard.press(columnName[0]);
    await page.waitForFunction((sel) => document.activeElement?.matches(sel) === true, SEARCH_INPUT, {timeout: 3000});
  }
  if (columnName.length > 1) await page.keyboard.type(columnName.slice(1));
  await page.waitForFunction(({sel, text}) =>
    (document.querySelector(sel) as HTMLInputElement | null)?.value === text, {sel: SEARCH_INPUT, text: columnName}, {timeout: 3000});
  await page.keyboard.press('Enter');
  await page.locator('.d4-column-selector-search').waitFor({state: 'detached', timeout: 3000}).catch(() => {});
}

async function leafBox(page: Page, leafName: string) {
  return page.evaluate((name) => {
    const el = document.querySelector(`.d4-menu-popup[name="pivot-tag"] [name="${name}"]`) as HTMLElement | null;
    if (!el) return null;
    const r = el.getBoundingClientRect();
    return r.width > 0 && r.height > 0 ? {x: r.x + r.width / 2, y: r.y + r.height / 2} : null;
  }, leafName);
}

export async function expandSubmenu(page: Page, parent: 'Aggregation' | 'Column', leaf: string): Promise<void> {
  const box = await page.locator(`.d4-menu-popup[name="pivot-tag"] [name="div-${parent}"]`).boundingBox();
  if (!box) throw new Error(`pivot-tag ${parent} parent not found`);
  const px = box.x + box.width / 2, py = box.y + box.height / 2;
  let jitter = 0;
  const ready = await v.pollValue(async () => {
    await page.mouse.move(px + (jitter++ % 2), py);
    return (await leafBox(page, leaf)) !== null;
  }, (r) => r, 4000, 100);
  if (!ready) throw new Error(`pivot-tag ${parent} flyout did not reveal ${leaf}`);
}

export async function clickFlyoutLeaf(page: Page, parent: 'Aggregation' | 'Column', leafName: string): Promise<void> {
  for (let attempt = 0; attempt < 3; attempt++) {
    await expandSubmenu(page, parent, leafName);
    const box = await leafBox(page, leafName);
    if (!box) continue;
    await page.mouse.move(box.x, box.y);
    const onLeaf = await v.pollValue(() => page.evaluate(({x, y, name}) => {
      const top = document.elementFromPoint(x, y);
      const el = document.querySelector(`.d4-menu-popup[name="pivot-tag"] [name="${name}"]`);
      return !!(top && el && (el.contains(top) || el === top));
    }, {x: box.x, y: box.y, name: leafName}), (on) => on, 500, 25);
    if (!onLeaf) continue;
    await page.mouse.down();
    await page.mouse.up();
    return;
  }
  throw new Error(`could not click pivot-tag ${parent} leaf ${leafName}`);
}

export const pickAggregation = (page: Page, aggType: string) =>
  clickFlyoutLeaf(page, 'Aggregation', `div-Aggregation---${aggType}`);
export const pickColumn = (page: Page, colLeafName: string) =>
  clickFlyoutLeaf(page, 'Column', `div-Column---${colLeafName}`);

export async function openChipMenu(page: Page, rowTitle: string, containing?: string): Promise<void> {
  await page.evaluate(({title, text}) => {
    const root = document.querySelector('[name="viewer-Pivot-table"]');
    const panel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
      .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === title);
    const chips = Array.from(panel!.querySelectorAll('.d4-tag'));
    const chip = text ? chips.find((t) => (t.textContent ?? '').includes(text)) : chips[0];
    chip!.dispatchEvent(new MouseEvent('contextmenu', {bubbles: true, cancelable: true, button: 2}));
  }, {title: rowTitle, text: containing ?? null});
  await page.locator('.d4-menu-popup[name="pivot-tag"] [name="div-Aggregation"]').waitFor({timeout: 5000});
}

export async function checkedAggregation(page: Page): Promise<string[]> {
  return page.evaluate(() => {
    const menu = document.querySelector('.d4-menu-popup[name="pivot-tag"]');
    if (!menu) return [];
    return Array.from(menu.querySelectorAll('[name^="div-Aggregation---"]'))
      .filter((mi) => mi.querySelector('.d4-menu-item-check i')?.className.includes('fa-dot-circle'))
      .map((mi) => mi.getAttribute('d4-name') ?? '');
  });
}

export async function closePivotTagMenu(page: Page): Promise<void> {
  await page.keyboard.press('Escape');
  await page.locator('.d4-menu-popup[name="pivot-tag"]').waitFor({state: 'detached', timeout: 5000}).catch(() => {});
}

export async function removeChip(page: Page, rowTitle: string, containing?: string): Promise<void> {
  await page.evaluate(({title, text}) => {
    const root = document.querySelector('[name="viewer-Pivot-table"]');
    const panel = Array.from(root!.querySelectorAll('.grok-pivot-column-panel'))
      .find((p) => p.querySelector('.grok-pivot-column-tags-title')?.getAttribute('d4-name') === title);
    const chips = Array.from(panel!.querySelectorAll('.d4-tag'));
    const chip = text ? chips.find((t) => (t.textContent ?? '').includes(text)) : chips[0];
    (chip!.querySelector('i, .grok-icon') as HTMLElement | null)?.click();
  }, {title: rowTitle, text: containing ?? null});
}

export interface GridLookColumn {
  columnName: string; visible: boolean; colorCodingType: string; isColorCoded: boolean; width: number;
}

export async function gridLookColumn(page: Page, columnName: string): Promise<GridLookColumn | undefined> {
  const cols: GridLookColumn[] = await page.evaluate(() => {
    const pv = Array.from((window as any).grok.shell.tv.viewers).find((x: any) => x.type === 'Pivot table') as any;
    return (pv.getOptions(true).look.gridLook?.columns ?? []).map((c: any) => ({
      columnName: c.columnName, visible: c.visible,
      colorCodingType: c.colorCodingType, isColorCoded: c.isColorCoded, width: c.width,
    }));
  });
  return cols.find((c) => c.columnName === columnName);
}

export async function pivotTitleDom(page: Page): Promise<string[]> {
  return page.evaluate(() =>
    Array.from(document.querySelectorAll('.panel-titlebar-tabhost .panel-titlebar-text'))
      .map((n) => (n.textContent || '').trim())
      .filter(Boolean));
}

export const HEADER_Y = 12;
export const rowY = (r: number) => 24 + r * 24 + 12;

export async function innerRect(page: Page) {
  const canvas = page.locator(INNER_CANVAS).first();
  await canvas.waitFor({state: 'visible', timeout: 10000});
  await canvas.scrollIntoViewIfNeeded().catch(() => {});
  const vp = page.viewportSize() ?? {width: 1920, height: 1080};
  const box = await v.pollValue(() => canvas.boundingBox(), (b) => !!b && b.width > 20 && b.height > 20 &&
    b.x >= 0 && b.y >= 0 && b.x + 40 <= vp.width && b.y + 20 <= vp.height, 10000, 100);
  if (!box) throw new Error('inner pivot grid canvas not visible');
  return box;
}

export function headerPoint(box: {x: number; y: number; width: number; height: number},
    vp: {width: number; height: number}, xLocal: number): {px: number; py: number} {
  const clamp = (val: number, lo: number, hi: number) => Math.max(lo, Math.min(hi, val));
  const px = clamp(clamp(box.x + xLocal, box.x + 2, box.x + box.width - 2), 1, vp.width - 1);
  const py = clamp(clamp(box.y + HEADER_Y, box.y + 2, box.y + box.height - 2), 1, vp.height - 1);
  return {px, py};
}

export async function rightClickHeader(page: Page, xLocal: number): Promise<void> {
  const vp = page.viewportSize() ?? {width: 1920, height: 1080};
  let lastErr: unknown = null;
  for (let attempt = 0; attempt < 3; attempt++) {
    try {
      const box = await innerRect(page);
      const {px, py} = headerPoint(box, vp, xLocal);
      await page.evaluate(({sel, cx, cy}) => {
        const cv = (document.querySelector(sel.replace('canvas[name="canvas"]', 'canvas[name="overlay"]'))
          ?? document.querySelector(sel)) as HTMLCanvasElement | null;
        if (!cv) throw new Error('inner-grid overlay canvas not found for contextmenu dispatch');
        cv.dispatchEvent(new MouseEvent('contextmenu',
          {bubbles: true, cancelable: true, clientX: cx, clientY: cy, button: 2, buttons: 2}));
      }, {sel: INNER_CANVAS, cx: px, cy: py});
      await page.locator('.d4-menu-popup').first().waitFor({timeout: 4000});
      return;
    } catch (e) {
      lastErr = e;
      await page.keyboard.press('Escape').catch(() => {});
      await page.locator('.d4-menu-popup').first().waitFor({state: 'detached', timeout: 400}).catch(() => {});
    }
  }
  throw new Error(`rightClickHeader: could not open the inner-grid header menu at xLocal=${xLocal} — ${String(lastErr)}`);
}

async function gridMenuItemBox(page: Page, name: string) {
  return page.evaluate((n) => {
    const el = document.querySelector(`.d4-menu-popup [name="${n}"]`) as HTMLElement | null;
    if (!el) return null;
    const r = el.getBoundingClientRect();
    return r.width > 0 && r.height > 0 ? {x: r.x + r.width / 2, y: r.y + r.height / 2} : null;
  }, name);
}

export async function expandGroup(page: Page, groupName: string, childName: string): Promise<void> {
  const box = await gridMenuItemBox(page, groupName);
  if (!box) throw new Error(`grid-menu group ${groupName} not found`);
  let jitter = 0;
  const ready = await v.pollValue(async () => {
    await page.mouse.move(box.x + (jitter++ % 2), box.y);
    return (await gridMenuItemBox(page, childName)) !== null;
  }, (r) => r, 4000, 100);
  if (!ready) throw new Error(`grid-menu ${groupName} flyout did not reveal ${childName}`);
}

export async function clickLeaf(page: Page, leafName: string): Promise<void> {
  const box = await gridMenuItemBox(page, leafName);
  if (!box) throw new Error(`grid-menu leaf ${leafName} not laid out`);
  await page.mouse.move(box.x, box.y);
  await v.pollValue(() => page.evaluate(({x, y, name}) => {
    const top = document.elementFromPoint(x, y);
    const el = document.querySelector(`.d4-menu-popup [name="${name}"]`);
    return !!(top && el && (el.contains(top) || el === top));
  }, {x: box.x, y: box.y, name: leafName}), (on) => on, 500, 25);
  await page.mouse.down();
  await page.mouse.up();
}

export async function clickGridLeaf(page: Page, leafName: string): Promise<void> {
  await expandGroup(page, 'div-Grid', leafName);
  await clickLeaf(page, leafName);
}

export async function applyLinearColorCoding(page: Page, xLocal: number): Promise<void> {
  await rightClickHeader(page, xLocal);
  await expandGroup(page, 'div-Grid', 'div-Grid---Color-Coding');
  await expandGroup(page, 'div-Grid---Color-Coding', 'div-Grid---Color-Coding---Linear');
  await clickLeaf(page, 'div-Grid---Color-Coding---Linear');
  await page.keyboard.press('Escape');
}
