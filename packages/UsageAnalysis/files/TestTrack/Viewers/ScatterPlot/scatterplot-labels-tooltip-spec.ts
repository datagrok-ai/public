/* ---
realizes: [scatterplot.cp.labels-tooltip, viewers.scatter-plot]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const SETUP_X = 'WEIGHT';
const SETUP_Y = 'HEIGHT';

const REFERENCE_ROW = 10;
const CUSTOM_TOOLTIP_COLUMNS = ['AGE', 'SEX'];
const LABEL_COLUMN = 'AGE';
const SPGI_X = 'Whole blood assay 2 Date';
const SPGI_Y = 'Stereo Category';
const SPGI_LABEL_COLUMN = 'Id';

const INK_SETTLE = 60;
const LABEL_INK_DELTA = 500;
const OVERLAY_RESTORE_TOLERANCE = 150;

const MARKER_AIM_TOLERANCE = 0.05;

const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Failed to connect to Claude runtime/.test(text) ||
  /powerPreference option is currently ignored/.test(text) ||
  /willReadFrequently/.test(text);

const isBenignError = (text: string) => isAmbientError(text);

const CRASH_SIGNATURE = /Infinity\.ceil/i;

interface Rect {x: number; y: number; width: number; height: number}

const canvasRect = (page: Page): Promise<Rect> => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const r = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
  return {x: r.x, y: r.y, width: r.width, height: r.height};
});

const pickOnViewer = (page: Page, role: string, column: string) =>
  v.pickColumnViaSelectorTrusted(page, {role, columnName: column});

/** [rebind] re-opens the gear even when a panel is already built: after openTable adds a
 *  NEW view + viewer, the existing panel still edits the PREVIOUS view's viewer, so writes
 *  silently land on a viewer nothing is asserting against. */
async function openSettings(page: Page, _rebind = false): Promise<void> {
  await v.openViewerSettings(page, 'Scatter plot');
}

async function revealPropEditor(page: Page, editorSelector: string, category: string): Promise<void> {
  for (let i = 0; i < 8; i++) {
    const ready = await v.pollValue(() => page.evaluate((sel: string) => {
      const el = document.querySelector(sel) as HTMLElement | null;
      if (!el || !el.offsetParent) return false;
      const b = el.getBoundingClientRect();
      return b.width > 0 && b.height > 0;
    }, editorSelector), (ok) => ok, i === 0 ? 300 : 900, 100);
    if (ready) return;
    const header = page.locator(`[name="prop-category-${category}"]`);
    if (await header.count() > 0 && await header.isVisible()) await header.click();
  }
  throw new Error(`property editor ${editorSelector} never became reachable`);
}

async function setChoiceProp(
  page: Page, rowName: string, viewCell: string, category: string, value: string,
): Promise<void> {
  await openSettings(page);
  await revealPropEditor(page, `[name="${viewCell}"]`, category);
  const cell = page.locator(`[name="${viewCell}"]`);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  const editor = page.locator(`[name="${rowName}"] select.property-grid-item-editor-spinner`);
  await editor.waitFor({state: 'visible', timeout: 4000});
  await editor.selectOption(value);
  const prop = rowName.replace(/^prop-/, '').replace(/-([a-z])/g, (_, c: string) => c.toUpperCase());
  await v.pollValue(() => readProp(page, prop), (x) => x === value, 2500, 100);
  await v.waitForViewerRendered(page, 'Scatter plot', 250);
}

const readProp = (page: Page, name: string) => page.evaluate((n: string) => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return sp.props[n];
}, name);

async function setListProp(page: Page, name: string, value: string | string[]): Promise<void> {
  await v.setViewerProps(page, 'Scatter plot', [{set: {[name]: value}, wait: 400}]);
}

const canvasInk = (page: Page, layer: 'canvas' | 'overlay') => page.evaluate((name: string) => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'));
  const c = root?.querySelector(`canvas[name="${name}"]`) as HTMLCanvasElement | null;
  const ctx = c?.getContext('2d');
  if (!c || !ctx) return -1;
  let data: Uint8ClampedArray;
  try { data = ctx.getImageData(0, 0, c.width, c.height).data; } catch (_) { return -1; }
  let n = 0;
  for (let k = 0; k < data.length; k += 16)
    if (data[k + 3] !== 0 && !(data[k] >= 250 && data[k + 1] >= 250 && data[k + 2] >= 250)) n++;
  return n;
}, layer);

async function parkPointer(page: Page): Promise<void> {
  await page.mouse.move(4, 4);
  await v.waitForViewerRendered(page, 'Scatter plot', 250);
}

async function settledInk(page: Page, layer: 'canvas' | 'overlay'): Promise<number> {
  await parkPointer(page);
  let prev = await canvasInk(page, layer);
  let cur = prev;
  for (let i = 0; i < 12; i++) {

    await page.waitForTimeout(250);
    cur = await canvasInk(page, layer);
    if (cur >= 0 && Math.abs(cur - prev) <= INK_SETTLE) break;
    prev = cur;
  }
  return cur;
}

interface TooltipEntry {name: string; value: string}

const tooltipEntries = (page: Page): Promise<TooltipEntry[]> => page.evaluate(() => {
  const tt = document.querySelector('.d4-tooltip') as HTMLElement | null;
  if (!tt) return [];
  return [...tt.querySelectorAll('tr')].map((tr) => {
    const value = (tr.querySelector('.d4-tooltip-text-value-cell')?.textContent ?? '').trim();
    const cells = [...tr.querySelectorAll('td')];
    const first = cells.length > 0 ? (cells[0].textContent ?? '').trim() : '';
    return {name: first === value ? '' : first, value};
  }).filter((e) => e.value.length > 0 || e.name.length > 0);
});

const tooltipText = (page: Page): Promise<string> => page.evaluate(() => {
  const tt = document.querySelector('.d4-tooltip') as HTMLElement | null;
  if (!tt || tt.offsetParent === null) return '';
  return (tt.textContent ?? '').trim();
});

interface HoverResult {
  hitRow: number;
  entries: TooltipEntry[];
  aimOffsetX: number;
  aimOffsetY: number;
  viewportWidth: number;
  viewportHeight: number;
}

async function hoverReferenceMarker(
  page: Page, row: number, expectTooltip = true,
): Promise<HoverResult> {
  const rect = await canvasRect(page);
  await page.mouse.move(rect.x + rect.width * 0.03, rect.y + rect.height * 0.03);
  const parked = await v.pollValue(() => tooltipText(page), (t) => t.length === 0, 800, 100);
  const aim = await page.evaluate((i: number) => {
    const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
      .find((e) => !e.closest('.d4-dialog'))!;
    const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
    const df = grok.shell.tv.dataFrame;
    const s = sp.worldToScreen(df.col(sp.props.xColumnName).get(i), df.col(sp.props.yColumnName).get(i));
    const b = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
    return {x: b.x + s.x, y: b.y + s.y};
  }, row);
  await page.mouse.move(aim.x, aim.y, {steps: 8});
  if (expectTooltip) {
    await v.pollValue(() => tooltipText(page), (t) => t.length > 0 && t !== parked, 3000, 100);

    await v.pollValue(() => tooltipEntries(page), (e) => e.length > 0, 200, 50);
  } else

    await v.waitForViewerRendered(page, 'Scatter plot', 1600);
  const state = await page.evaluate((i: number) => {
    const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
    const df = grok.shell.tv.dataFrame;
    const hit = df.mouseOverRowIdx as number;
    const xc = df.col(sp.props.xColumnName);
    const yc = df.col(sp.props.yColumnName);
    const vp = sp.viewport;
    return {
      hitRow: hit,
      aimOffsetX: hit >= 0 ? Math.abs(xc.get(hit) - xc.get(i)) : -1,
      aimOffsetY: hit >= 0 ? Math.abs(yc.get(hit) - yc.get(i)) : -1,
      viewportWidth: vp.width as number,
      viewportHeight: vp.height as number,
    };
  }, row);
  return {...state, entries: await tooltipEntries(page)};
}

const rowValues = (page: Page, row: number, columns: string[]) =>
  page.evaluate(({i, cols}: {i: number; cols: string[]}) => {
    const df = grok.shell.tv.dataFrame;
    const out: Record<string, string | null> = {};
    for (const c of cols) {
      const col = df.col(c);
      out[c] = col ? String(col.getString(i)) : null;
    }
    return out;
  }, {i: row, cols: columns});

interface Frac {fx: number; fy: number}

const at = (r: Rect, p: Frac) => ({x: r.x + r.width * p.fx, y: r.y + r.height * p.fy});

async function dragCanvas(page: Page, from: Frac, to: Frac, mods: string[] = []): Promise<void> {
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
  await v.waitForViewerRendered(page, 'Scatter plot', 300);
}

async function clickCanvas(page: Page, p: Frac): Promise<void> {
  const r = await canvasRect(page);
  const pt = at(r, p);
  await page.mouse.click(pt.x, pt.y);
  await v.waitForViewerRendered(page, 'Scatter plot', 250);
}

const liveCount = (page: Page, kind: 'selection' | 'filter') =>
  page.evaluate((k: string) =>
    (k === 'selection' ? grok.shell.tv.dataFrame.selection.trueCount
      : grok.shell.tv.dataFrame.filter.trueCount) as number, kind);

async function settledCount(page: Page, kind: 'selection' | 'filter'): Promise<number> {
  let last = -1;
  let stable = 0;
  for (let i = 0; i < 30; i++) {
    const c = await liveCount(page, kind);
    if (c === last) {
      stable++;
      if (stable >= 2) return c;
    } else {
      stable = 0;
      last = c;
    }

    await page.waitForTimeout(200);
  }
  return last;
}

async function settledCountAfterChange(
  page: Page, kind: 'selection' | 'filter', from: number, timeoutMs = 7000,
): Promise<number> {
  const deadline = Date.now() + timeoutMs;
  let prev = from;
  while (Date.now() < deadline) {
    const c = await liveCount(page, kind);
    if (c !== from) {
      if (c === prev) return c;
      prev = c;
    }

    await page.waitForTimeout(120);
  }
  return await liveCount(page, kind);
}

const rowCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.rowCount as number);

async function filterCardIndex(page: Page, column: string, timeoutMs = 45_000): Promise<number> {
  return v.pollValue(() => page.evaluate((col: string) =>
    [...document.querySelectorAll('[name="viewer-Filters"] .d4-filter')]
      .findIndex((c) => (c.querySelector('.d4-filter-column-name')?.textContent ?? '').trim() === col),
  column), (idx) => idx >= 0, timeoutMs, 300);
}

interface CategoryState {filtered: number; survivors: string[]}

const categoryState = (page: Page, column: string): Promise<CategoryState> =>
  page.evaluate((col: string) => {
    const df = grok.shell.tv.dataFrame;
    const c = df.col(col);
    const counts: Record<string, number> = {};
    for (const x of c.categories) counts[x] = 0;
    for (let i = 0; i < df.rowCount; i++) if (df.filter.get(i)) counts[c.get(i)]++;
    return {
      filtered: df.filter.trueCount as number,
      survivors: Object.keys(counts).filter((k) => counts[k] > 0),
    };
  }, column);

async function resetFilterPanel(page: Page, from?: number): Promise<void> {
  const reset = page.locator('[name="viewer-Filters"] [name="icon-arrow-rotate-left"]').first();
  await reset.scrollIntoViewIfNeeded();
  await reset.click();

  if (from === undefined) await v.waitForViewerRendered(page, 'Scatter plot', 800);
  else await settledCountAfterChange(page, 'filter', from, 8000);
}

async function narrowToOneCategory(
  page: Page, cardIndex: number, column: string, fullRows: number,
): Promise<CategoryState | null> {
  const card = page.locator('[name="viewer-Filters"] .d4-filter').nth(cardIndex);
  const canvas = card.locator('canvas').last();
  for (const fy of [0.5, 0.3, 0.7, 0.15, 0.85]) {
    await canvas.scrollIntoViewIfNeeded();
    const box = await canvas.boundingBox();
    if (!box || box.width === 0 || box.height === 0) break;
    await canvas.click({position: {x: box.width / 2, y: box.height * fy}});
    await settledCountAfterChange(page, 'filter', fullRows, 3500);
    const state = await categoryState(page, column);
    if (state.survivors.length === 1 && state.filtered > 0 && state.filtered < fullRows) return state;
    await resetFilterPanel(page, state.filtered === fullRows ? undefined : state.filtered);
  }
  return null;
}

test('Scatter Plot — Marker Labels and Tooltip', async ({page}: {page: Page}) => {
  test.setTimeout(1_500_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  const allMessages: string[] = [];
  page.on('pageerror', (e) => {
    allMessages.push(String(e));
    if (!isBenignError(String(e))) pageErrors.push(String(e));
  });
  page.on('console', (m) => {
    if (m.type() !== 'error') return;
    allMessages.push(m.text());
    if (!isBenignError(m.text()) && !isLocalBootNoise(m.text())) consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;
  const crashCount = () => allMessages.filter((t) => CRASH_SIGNATURE.test(t)).length;

  await openDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
  await v.waitForViewerRendered(page, 'Scatter plot', 500);
  await pickOnViewer(page, 'x', SETUP_X);
  await pickOnViewer(page, 'y', SETUP_Y);
  expect(await readProp(page, 'xColumnName')).toBe(SETUP_X);
  expect(await readProp(page, 'yColumnName')).toBe(SETUP_Y);
  expect(await readProp(page, 'showTooltip')).toBe('inherit from table');
  expect(await readProp(page, 'labelColumnNames')).toEqual([]);

  await softStep('Tooltip inherited from the table', async () => {
    const errBefore = errCount();
    const hover = await hoverReferenceMarker(page, REFERENCE_ROW);
    expect(hover.hitRow).toBeGreaterThanOrEqual(0);
    expect(hover.aimOffsetX).toBeLessThanOrEqual(hover.viewportWidth * MARKER_AIM_TOLERANCE);
    expect(hover.aimOffsetY).toBeLessThanOrEqual(hover.viewportHeight * MARKER_AIM_TOLERANCE);

    const entries = hover.entries;
    expect(entries.length).toBeGreaterThan(1);
    const named = entries.filter((e) => e.name.length > 0);
    expect(named.length).toBe(entries.length);
    const table = await rowValues(page, hover.hitRow, named.map((e) => e.name));
    const unknown = named.filter((e) => table[e.name] === null).map((e) => e.name);
    const mismatched = named
      .filter((e) => table[e.name] !== null && table[e.name] !== e.value)
      .map((e) => `${e.name}: tooltip ${e.value} vs table ${table[e.name]}`);
    console.log(`inherited tooltip: row=${hover.hitRow} columns=${named.map((e) => e.name).join(',')}`);

    expect(unknown).toEqual([]);
    expect(mismatched).toEqual([]);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Custom tooltip column list', async () => {
    const errBefore = errCount();
    await setChoiceProp(page, 'prop-show-tooltip', 'prop-view-show-tooltip', 'tooltip', 'show custom tooltip');
    expect(await readProp(page, 'showTooltip')).toBe('show custom tooltip');
    await setListProp(page, 'rowTooltip', CUSTOM_TOOLTIP_COLUMNS.join('\n'));
    expect(await readProp(page, 'rowTooltip')).toBe(CUSTOM_TOOLTIP_COLUMNS.join('\n'));
    await setChoiceProp(page, 'prop-data-values', 'prop-view-data-values', 'tooltip', 'Do not add');
    expect(await readProp(page, 'dataValues')).toBe('Do not add');

    const custom = await hoverReferenceMarker(page, REFERENCE_ROW);
    expect(custom.hitRow).toBeGreaterThanOrEqual(0);
    const configured = await rowValues(page, custom.hitRow, CUSTOM_TOOLTIP_COLUMNS);
    console.log(`custom tooltip: row=${custom.hitRow} values=${custom.entries.map((e) => e.value).join('|')}`);

    expect(custom.entries.map((e) => e.value))
      .toEqual(CUSTOM_TOOLTIP_COLUMNS.map((c) => configured[c]));

    await setChoiceProp(page, 'prop-show-tooltip', 'prop-view-show-tooltip', 'tooltip', 'do not show');
    expect(await readProp(page, 'showTooltip')).toBe('do not show');
    const silent = await hoverReferenceMarker(page, REFERENCE_ROW, false);
    expect(silent.hitRow).toBeGreaterThanOrEqual(0);
    expect(silent.entries).toEqual([]);

    await setChoiceProp(page, 'prop-show-tooltip', 'prop-view-show-tooltip', 'tooltip', 'inherit from table');
    await setListProp(page, 'rowTooltip', '');
    await setChoiceProp(page, 'prop-data-values', 'prop-view-data-values', 'tooltip', 'Merge');
    expect(await readProp(page, 'showTooltip')).toBe('inherit from table');
    expect(await readProp(page, 'rowTooltip')).toBe('');

    const restored = await hoverReferenceMarker(page, REFERENCE_ROW);
    expect(restored.hitRow).toBeGreaterThanOrEqual(0);
    expect(restored.entries.length).toBeGreaterThan(CUSTOM_TOOLTIP_COLUMNS.length);
    const namedAgain = restored.entries.filter((e) => e.name.length > 0);
    expect(namedAgain.length).toBe(restored.entries.length);
    const tableAgain = await rowValues(page, restored.hitRow, namedAgain.map((e) => e.name));
    expect(namedAgain.filter((e) => tableAgain[e.name] !== e.value)).toEqual([]);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Labels for the selected rows', async () => {
    const errBefore = errCount();
    // a corner click is not a reliable "clear selection": earlier steps can leave a row
    // selected and (0.02, 0.02) may still land on a marker. Clear it outright, then confirm.
    await page.evaluate(() => grok.shell.tv.dataFrame.selection.setAll(false));
    await v.waitForViewerRendered(page, 'Scatter plot', 300);
    expect(await settledCount(page, 'selection')).toBe(0);
    const baseline = await settledInk(page, 'overlay');

    expect(baseline).toBeGreaterThanOrEqual(0);

    await setListProp(page, 'labelColumnNames', [LABEL_COLUMN]);
    expect(await readProp(page, 'labelColumnNames')).toEqual([LABEL_COLUMN]);
    await setChoiceProp(page, 'prop-show-labels-for', 'prop-view-show-labels-for', 'labels', 'Selected');
    expect(await readProp(page, 'showLabelsFor')).toBe('Selected');

    const configured = await settledInk(page, 'overlay');
    expect(Math.abs(configured - baseline)).toBeLessThan(OVERLAY_RESTORE_TOLERANCE);

    await dragCanvas(page, {fx: 0.3, fy: 0.3}, {fx: 0.7, fy: 0.7}, ['Shift']);
    const selected = await settledCountAfterChange(page, 'selection', 0);
    expect(selected).toBeGreaterThan(0);
    const labelled = await settledInk(page, 'overlay');
    console.log(`label overlay ink: baseline=${baseline} configured=${configured} ` +
      `labelled=${labelled} selected=${selected}`);

    expect(labelled - baseline).toBeGreaterThanOrEqual(LABEL_INK_DELTA);
    expect(errCount()).toBe(errBefore);

    await clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await settledCountAfterChange(page, 'selection', selected)).toBe(0);
    const cleared = await settledInk(page, 'overlay');
    console.log(`label overlay ink: cleared=${cleared}`);
    expect(Math.abs(cleared - baseline)).toBeLessThan(OVERLAY_RESTORE_TOLERANCE);

    await setListProp(page, 'labelColumnNames', []);
    await setChoiceProp(page, 'prop-show-labels-for', 'prop-view-show-labels-for', 'labels', 'All');
    expect(await readProp(page, 'labelColumnNames')).toEqual([]);
    expect(await readProp(page, 'showLabelsFor')).toBe('All');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Labels plus filtering on a datetime axis', async () => {
    await v.openTable(page, {path: spgiPath});
    const fullRows = await rowCount(page);
    expect(fullRows).toBeGreaterThan(0);
    await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
    await v.waitForViewerRendered(page, 'Scatter plot', 600);
    // this is a NEW view + viewer; rebind the settings panel or later prop writes go to
    // the previous view's scatter plot and this one keeps its defaults
    await openSettings(page, true);
    await pickOnViewer(page, 'x', SPGI_X);
    await pickOnViewer(page, 'y', SPGI_Y);
    expect(await readProp(page, 'xColumnName')).toBe(SPGI_X);
    expect(await readProp(page, 'yColumnName')).toBe(SPGI_Y);

    expect(await page.evaluate((c: string) =>
      String(grok.shell.tv.dataFrame.col(c).type), SPGI_X)).toBe('datetime');

    await setListProp(page, 'labelColumnNames', [SPGI_LABEL_COLUMN]);
    await setChoiceProp(page, 'prop-show-labels-for', 'prop-view-show-labels-for', 'labels', 'Selected');
    expect(await readProp(page, 'labelColumnNames')).toEqual([SPGI_LABEL_COLUMN]);
    expect(await readProp(page, 'showLabelsFor')).toBe('Selected');

    const errBefore = errCount();
    const crashBefore = crashCount();
    await v.openFilterPanel(page);
    const cardIndex = await filterCardIndex(page, SPGI_Y);
    expect(cardIndex).toBeGreaterThanOrEqual(0);
    const narrowed = await narrowToOneCategory(page, cardIndex, SPGI_Y, fullRows);
    expect(narrowed).not.toBeNull();

    console.log(`stereo category filter: kept=${narrowed!.survivors[0]} ` +
      `rows=${narrowed!.filtered}/${fullRows}`);
    expect(narrowed!.survivors.length).toBe(1);
    expect(narrowed!.filtered).toBeGreaterThan(0);
    expect(narrowed!.filtered).toBeLessThan(fullRows);

    await dragCanvas(page, {fx: 0.2, fy: 0.2}, {fx: 0.8, fy: 0.8}, ['Shift']);
    const selected = await settledCountAfterChange(page, 'selection', 0);
    expect(selected).toBeGreaterThan(0);

    await v.waitForViewerRendered(page, 'Scatter plot', 2000);
    expect(crashCount()).toBe(crashBefore);
    expect(errCount()).toBe(errBefore);

    await resetFilterPanel(page, narrowed!.filtered);
    await clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await settledCountAfterChange(page, 'selection', selected)).toBe(0);
    expect(await settledCount(page, 'filter')).toBe(fullRows);
    expect(crashCount()).toBe(crashBefore);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
