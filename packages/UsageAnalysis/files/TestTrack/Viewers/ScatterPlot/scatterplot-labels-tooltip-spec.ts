/* ---
realizes: [scatterplot.cp.labels-tooltip, viewers.scatter-plot]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as sp from './scatterplot-shared';

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

const CRASH_SIGNATURE = /Infinity\.ceil/i;

const setTooltipChoice = (page: Page, rowName: string, value: string) =>
  sp.setChoiceProp(page, rowName, 'tooltip', value);

const setListProp = (page: Page, name: string, value: string | string[]) =>
  v.setViewerProps(page, sp.SP_TYPE, [{set: {[name]: value}, wait: 400}]);

const settledOverlay = (page: Page) => sp.settledInk(page, 'overlay', INK_SETTLE);

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

async function hoverReferenceMarker(page: Page, row: number, expectTooltip = true): Promise<HoverResult> {
  const rect = await sp.canvasRect(page);
  await page.mouse.move(rect.x + rect.width * 0.03, rect.y + rect.height * 0.03);
  const parked = await v.pollValue(() => tooltipText(page), (t) => t.length === 0, 800, 50);
  const aim = await page.evaluate((i: number) => {
    const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
    const df = grok.shell.tv.dataFrame;
    const pt = s.worldToScreen(df.col(s.props.xColumnName).get(i), df.col(s.props.yColumnName).get(i));
    const b = s.root.querySelector('canvas[name="canvas"]').getBoundingClientRect();
    return {x: b.x + pt.x, y: b.y + pt.y};
  }, row);
  const before = await page.evaluate(() => grok.shell.tv.dataFrame.mouseOverRowIdx as number);
  await page.mouse.move(aim.x, aim.y, {steps: 8});
  if (expectTooltip) {
    await v.pollValue(() => tooltipText(page), (t) => t.length > 0 && t !== parked, 3000, 50);
    await v.pollValue(() => tooltipEntries(page), (e) => e.length > 0, 200, 50);
  } else
    await v.pollValue(() => page.evaluate(() => grok.shell.tv.dataFrame.mouseOverRowIdx as number),
      (r) => r >= 0 && r !== before, 1600, 50);
  const state = await page.evaluate((i: number) => {
    const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
    const df = grok.shell.tv.dataFrame;
    const hit = df.mouseOverRowIdx as number;
    const xc = df.col(s.props.xColumnName);
    const yc = df.col(s.props.yColumnName);
    const vp = s.viewport;
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

const rowCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.rowCount as number);

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

async function resetFilterPanel(page: Page, from: number): Promise<void> {
  const reset = page.locator('[name="viewer-Filters"] [name="icon-arrow-rotate-left"]').first();
  await reset.scrollIntoViewIfNeeded();
  await reset.click();
  await sp.filterMoved(page, from, 8000);
}

test('Scatter Plot — Marker Labels and Tooltip', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const allMessages: string[] = [];
  page.on('pageerror', (e) => allMessages.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') allMessages.push(m.text()); });
  const errors = sp.trackErrors(page);
  const errCount = errors.count;
  const crashCount = () => allMessages.filter((t) => CRASH_SIGNATURE.test(t)).length;

  await openDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await sp.addScatterPlot(page);
  await sp.pickOnViewer(page, 'x', SETUP_X);
  await sp.pickOnViewer(page, 'y', SETUP_Y);
  expect(await sp.readProp(page, 'xColumnName')).toBe(SETUP_X);
  expect(await sp.readProp(page, 'yColumnName')).toBe(SETUP_Y);
  expect(await sp.readProp(page, 'showTooltip')).toBe('inherit from table');
  expect(await sp.readProp(page, 'labelColumnNames')).toEqual([]);

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
    await setTooltipChoice(page, 'prop-show-tooltip', 'show custom tooltip');
    expect(await sp.readProp(page, 'showTooltip')).toBe('show custom tooltip');
    await setListProp(page, 'rowTooltip', CUSTOM_TOOLTIP_COLUMNS.join('\n'));
    expect(await sp.readProp(page, 'rowTooltip')).toBe(CUSTOM_TOOLTIP_COLUMNS.join('\n'));
    await setTooltipChoice(page, 'prop-data-values', 'Do not add');
    expect(await sp.readProp(page, 'dataValues')).toBe('Do not add');

    const custom = await hoverReferenceMarker(page, REFERENCE_ROW);
    expect(custom.hitRow).toBeGreaterThanOrEqual(0);
    const configured = await rowValues(page, custom.hitRow, CUSTOM_TOOLTIP_COLUMNS);
    console.log(`custom tooltip: row=${custom.hitRow} values=${custom.entries.map((e) => e.value).join('|')}`);

    expect(custom.entries.map((e) => e.value))
      .toEqual(CUSTOM_TOOLTIP_COLUMNS.map((c) => configured[c]));

    await setTooltipChoice(page, 'prop-show-tooltip', 'do not show');
    expect(await sp.readProp(page, 'showTooltip')).toBe('do not show');
    const silent = await hoverReferenceMarker(page, REFERENCE_ROW, false);
    expect(silent.hitRow).toBeGreaterThanOrEqual(0);
    expect(silent.entries).toEqual([]);

    await setTooltipChoice(page, 'prop-show-tooltip', 'inherit from table');
    await setListProp(page, 'rowTooltip', '');
    await setTooltipChoice(page, 'prop-data-values', 'Merge');
    expect(await sp.readProp(page, 'showTooltip')).toBe('inherit from table');
    expect(await sp.readProp(page, 'rowTooltip')).toBe('');

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
    await v.waitForViewerRendered(page, sp.SP_TYPE, 300);
    expect(await sp.selectionHeld(page)).toBe(0);
    const baseline = await settledOverlay(page);
    expect(baseline).toBeGreaterThanOrEqual(0);

    await setListProp(page, 'labelColumnNames', [LABEL_COLUMN]);
    expect(await sp.readProp(page, 'labelColumnNames')).toEqual([LABEL_COLUMN]);
    // Show Labels For is still at its default, All, so every row is labelled here; Selected
    // takes those labels away again, and that is the change waited for
    const allLabelled = await settledOverlay(page);
    await sp.setChoiceProp(page, 'prop-show-labels-for', 'labels', 'Selected');
    expect(await sp.readProp(page, 'showLabelsFor')).toBe('Selected');

    const configured = await sp.settledInk(page, 'overlay', INK_SETTLE, allLabelled);
    expect(Math.abs(configured - baseline)).toBeLessThan(OVERLAY_RESTORE_TOLERANCE);

    await sp.dragCanvas(page, {fx: 0.3, fy: 0.3}, {fx: 0.7, fy: 0.7}, ['Shift']);
    const selected = await sp.selectionMoved(page, 0);
    expect(selected).toBeGreaterThan(0);
    const labelled = await sp.settledInk(page, 'overlay', INK_SETTLE, baseline);
    console.log(`label overlay ink: baseline=${baseline} configured=${configured} ` +
      `labelled=${labelled} selected=${selected}`);

    expect(labelled - baseline).toBeGreaterThanOrEqual(LABEL_INK_DELTA);
    expect(errCount()).toBe(errBefore);

    await sp.clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await sp.selectionMoved(page, selected)).toBe(0);
    const cleared = await sp.settledInk(page, 'overlay', INK_SETTLE, labelled);
    console.log(`label overlay ink: cleared=${cleared}`);
    expect(Math.abs(cleared - baseline)).toBeLessThan(OVERLAY_RESTORE_TOLERANCE);

    await setListProp(page, 'labelColumnNames', []);
    await sp.setChoiceProp(page, 'prop-show-labels-for', 'labels', 'All');
    expect(await sp.readProp(page, 'labelColumnNames')).toEqual([]);
    expect(await sp.readProp(page, 'showLabelsFor')).toBe('All');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Labels plus filtering on a datetime axis', async () => {
    await v.openTable(page, {path: spgiPath});
    const fullRows = await rowCount(page);
    expect(fullRows).toBeGreaterThan(0);
    await sp.addScatterPlot(page);
    await sp.pickOnViewer(page, 'x', SPGI_X);
    await sp.pickOnViewer(page, 'y', SPGI_Y);
    expect(await sp.readProp(page, 'xColumnName')).toBe(SPGI_X);
    expect(await sp.readProp(page, 'yColumnName')).toBe(SPGI_Y);

    expect(await page.evaluate((c: string) =>
      String(grok.shell.tv.dataFrame.col(c).type), SPGI_X)).toBe('datetime');

    await setListProp(page, 'labelColumnNames', [SPGI_LABEL_COLUMN]);
    await sp.setChoiceProp(page, 'prop-show-labels-for', 'labels', 'Selected');
    expect(await sp.readProp(page, 'labelColumnNames')).toEqual([SPGI_LABEL_COLUMN]);
    expect(await sp.readProp(page, 'showLabelsFor')).toBe('Selected');

    const errBefore = errCount();
    const crashBefore = crashCount();
    await sp.openFilterPanel(page);
    const firstCategory = await page.evaluate((col: string) =>
      grok.shell.tv.dataFrame.col(col).categories[0] as string, SPGI_Y);
    await v.applyCategoricalFilter(page, SPGI_Y, [firstCategory], 600);
    const narrowedCount = await sp.filterMoved(page, fullRows);
    const narrowed = await categoryState(page, SPGI_Y);

    console.log(`stereo category filter: kept=${narrowed.survivors[0]} rows=${narrowed.filtered}/${fullRows}`);
    expect(narrowed.filtered).toBe(narrowedCount);
    expect(narrowed.survivors).toEqual([firstCategory]);
    expect(narrowed.filtered).toBeGreaterThan(0);
    expect(narrowed.filtered).toBeLessThan(fullRows);

    await sp.dragCanvas(page, {fx: 0.03, fy: 0.03}, {fx: 0.97, fy: 0.97}, ['Shift']);
    const selected = await sp.selectionMoved(page, 0);
    expect(selected).toBeGreaterThan(0);

    await v.waitForViewerQuiet(page, sp.SP_TYPE, {gapMs: 300, capMs: 2000});
    expect(crashCount()).toBe(crashBefore);
    expect(errCount()).toBe(errBefore);

    await resetFilterPanel(page, narrowed.filtered);
    await sp.clickCanvas(page, {fx: 0.02, fy: 0.02});
    expect(await sp.selectionMoved(page, selected)).toBe(0);
    expect(await sp.filterHeld(page, 500)).toBe(fullRows);
    expect(crashCount()).toBe(crashBefore);
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
