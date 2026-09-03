/* ---
realizes: []
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as sp from './scatterplot-shared';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const SETUP_X = 'WEIGHT';
const SETUP_Y = 'HEIGHT';
const WHISKER_X_MIN = 'AGE';
const WHISKER_X_MAX = 'WEIGHT';
const WHISKER_Y_MIN = 'HEIGHT';
const WHISKER_Y_MAX = 'WEIGHT';
const HISTOGRAM_BINS = 20;
const DEFAULT_HISTOGRAM_BINS = 10;
const TITLE_TEXT = 'Test Plot';
const DESCRIPTION_TEXT = 'Test description';
const LINES_ORDER_COLUMN = 'AGE';
const COLOR_COLUMN = 'RACE';
const LINES_BY_COLUMN = 'SEX';

const CANVAS_SETTLE_TOLERANCE = 400;
const CANVAS_CHANGE_MIN = 2000;
const CANVAS_RESTORE_MAX = 1500;

const captureCanvas = (page: Page, key: string) => page.evaluate((k: string) => {
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const c = s?.root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement | null;
  const ctx = c?.getContext('2d');
  if (!c || !ctx) return false;
  let data: Uint8ClampedArray;
  try { data = ctx.getImageData(0, 0, c.width, c.height).data; } catch (_) { return false; }
  const colors = new Map<number, number>();
  for (let i = 0; i < data.length; i += 4) {
    const rgb = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
    colors.set(rgb, (colors.get(rgb) ?? 0) + 1);
  }
  const w = window as any;
  w.__spCanvasSnap = w.__spCanvasSnap || {};
  w.__spCanvasSnap[k] = colors;
  return true;
}, key);

const diffCanvas = (page: Page, key: string) => page.evaluate((k: string) => {
  const w = window as any;
  const prev = w.__spCanvasSnap?.[k] as Map<number, number> | undefined;
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const c = s?.root.querySelector('canvas[name="canvas"]') as HTMLCanvasElement | null;
  const ctx = c?.getContext('2d');
  if (!prev || !c || !ctx) return -1;
  let data: Uint8ClampedArray;
  try { data = ctx.getImageData(0, 0, c.width, c.height).data; } catch (_) { return -1; }
  const colors = new Map<number, number>();
  for (let i = 0; i < data.length; i += 4) {
    const rgb = (data[i] << 16) | (data[i + 1] << 8) | data[i + 2];
    colors.set(rgb, (colors.get(rgb) ?? 0) + 1);
  }
  let delta = 0;
  for (const [rgb, n] of colors) delta += Math.abs(n - (prev.get(rgb) ?? 0));
  for (const [rgb, n] of prev) if (!colors.has(rgb)) delta += n;
  return delta;
}, key);

// The canvas delta against `key` once the plot is quiet; with expectChange, once it has first
// moved away from `key` by more than the settle tolerance.
async function settledCanvasDiff(page: Page, key: string, expectChange = false): Promise<number> {
  await sp.parkPointer(page);
  if (expectChange)
    await v.pollValue(() => diffCanvas(page, key), (d) => d > CANVAS_SETTLE_TOLERANCE, 4000, 100);
  await v.waitForViewerQuiet(page, sp.SP_TYPE, {gapMs: 250, capMs: 1500});
  const delta = await v.pollStable(() => diffCanvas(page, key), (a, b) => Math.abs(a - b) <= CANVAS_SETTLE_TOLERANCE, 3000, 100);
  expect(delta).toBeGreaterThanOrEqual(0);
  return delta;
}

async function captureBaseline(page: Page, key: string): Promise<void> {
  await sp.parkPointer(page);
  await v.waitForViewerQuiet(page, sp.SP_TYPE, {gapMs: 250, capMs: 1500});
  expect(await captureCanvas(page, key)).toBe(true);
}

const categoryCount = (page: Page, column: string) => page.evaluate((c: string) =>
  grok.shell.tv.dataFrame.col(c).categories.length as number, column);

const axisSelectorState = (page: Page) => page.evaluate(() => {
  const root = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot').root as HTMLElement;
  const read = (role: string) => {
    const el = root.querySelector(`[name="div-column-combobox-${role}"]`) as HTMLElement | null;
    if (!el) return null;
    const b = el.getBoundingClientRect();
    return {
      visibility: getComputedStyle(el).visibility,
      inDom: true,
      hasOffsetParent: !!el.offsetParent,
      width: b.width,
      height: b.height,
    };
  };
  return {x: read('x'), y: read('y')};
});

const viewerAlive = (page: Page) => page.evaluate(() => {
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return !!s && document.body.contains(s.root);
});

test('Scatter Plot — Secondary Settings Surface', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const errors = sp.trackErrors(page);
  const errCount = errors.count;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await sp.addScatterPlot(page);

  await sp.pickOnViewer(page, 'x', SETUP_X);
  await sp.pickOnViewer(page, 'y', SETUP_Y);
  await sp.openSettings(page);

  await softStep('Scenario 1 — Axis histograms', async () => {
    const errBefore = errCount();
    await captureBaseline(page, 'hist-base');

    await sp.setCheckboxProp(page, 'prop-show-x-histogram', 'x-axis', true);
    expect(await settledCanvasDiff(page, 'hist-base', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    expect(await captureCanvas(page, 'hist-x')).toBe(true);
    await sp.setCheckboxProp(page, 'prop-show-y-histogram', 'x-axis', true);
    expect(await settledCanvasDiff(page, 'hist-x', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    expect(await captureCanvas(page, 'hist-xy')).toBe(true);
    await sp.setNumericProp(page, 'prop-histogram-bins', 'x-axis', HISTOGRAM_BINS);
    expect(await settledCanvasDiff(page, 'hist-xy', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    await sp.setCheckboxProp(page, 'prop-show-x-histogram', 'x-axis', false);
    await sp.setCheckboxProp(page, 'prop-show-y-histogram', 'x-axis', false);
    await sp.setNumericProp(page, 'prop-histogram-bins', 'x-axis', DEFAULT_HISTOGRAM_BINS);
    expect(await settledCanvasDiff(page, 'hist-base')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 2 — Grid lines, axes and selector visibility (GROK-13533)', async () => {
    const errBefore = errCount();
    await captureBaseline(page, 'vis-base');
    const before = await axisSelectorState(page);
    expect(before.x?.visibility).toBe('visible');
    expect(before.y?.visibility).toBe('visible');

    await sp.setCheckboxProp(page, 'prop-show-vertical-grid-lines', 'x', false);
    await sp.setCheckboxProp(page, 'prop-show-horizontal-grid-lines', 'y', false);
    expect(await settledCanvasDiff(page, 'vis-base', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    expect(await captureCanvas(page, 'vis-nogrid')).toBe(true);
    await sp.setCheckboxProp(page, 'prop-show-x-axis', 'x', false);
    await sp.setCheckboxProp(page, 'prop-show-y-axis', 'y', false);
    expect(await settledCanvasDiff(page, 'vis-nogrid', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    await sp.setCheckboxProp(page, 'prop-show-x-selector', 'x', false);
    await sp.setCheckboxProp(page, 'prop-show-y-selector', 'y', false);
    await sp.parkPointer(page);
    const hidden = await v.pollValue(() => axisSelectorState(page),
      (s) => s.x?.visibility === 'hidden' && s.y?.visibility === 'hidden', 2000, 50);

    expect(hidden.x?.visibility).toBe('hidden');
    expect(hidden.y?.visibility).toBe('hidden');
    expect(hidden.x?.inDom).toBe(true);
    expect(hidden.y?.inDom).toBe(true);
    expect(hidden.x?.hasOffsetParent).toBe(true);
    expect(hidden.y?.hasOffsetParent).toBe(true);
    expect(hidden.x!.width).toBeGreaterThan(0);
    expect(hidden.x!.height).toBeGreaterThan(0);
    expect(hidden.y!.width).toBeGreaterThan(0);
    expect(hidden.y!.height).toBeGreaterThan(0);

    await sp.setCheckboxProp(page, 'prop-show-x-selector', 'x', true);
    await sp.setCheckboxProp(page, 'prop-show-y-selector', 'y', true);
    await sp.setCheckboxProp(page, 'prop-show-x-axis', 'x', true);
    await sp.setCheckboxProp(page, 'prop-show-y-axis', 'y', true);
    await sp.setCheckboxProp(page, 'prop-show-vertical-grid-lines', 'x', true);
    await sp.setCheckboxProp(page, 'prop-show-horizontal-grid-lines', 'y', true);
    await sp.parkPointer(page);
    const restored = await v.pollValue(() => axisSelectorState(page),
      (s) => s.x?.visibility === 'visible' && s.y?.visibility === 'visible', 2000, 50);
    expect(restored.x?.visibility).toBe('visible');
    expect(restored.y?.visibility).toBe('visible');
    expect(await settledCanvasDiff(page, 'vis-base')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 3 — Whiskers', async () => {
    const errBefore = errCount();
    await captureBaseline(page, 'whisker-base');

    await sp.pickPanelColumn(page, 'prop-x-whisker-min', 'div-column-combobox-x--whisker--min', 'x', WHISKER_X_MIN);
    await sp.pickPanelColumn(page, 'prop-x-whisker-max', 'div-column-combobox-x--whisker--max', 'x', WHISKER_X_MAX);
    await sp.pickPanelColumn(page, 'prop-y-whisker-min', 'div-column-combobox-y--whisker--min', 'y', WHISKER_Y_MIN);
    await sp.pickPanelColumn(page, 'prop-y-whisker-max', 'div-column-combobox-y--whisker--max', 'y', WHISKER_Y_MAX);
    expect(await settledCanvasDiff(page, 'whisker-base', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);
    expect(errCount()).toBe(errBefore);

    await sp.clearPanelColumn(page, 'prop-x-whisker-min', 'div-column-combobox-x--whisker--min', 'x');
    await sp.clearPanelColumn(page, 'prop-x-whisker-max', 'div-column-combobox-x--whisker--max', 'x');
    await sp.clearPanelColumn(page, 'prop-y-whisker-min', 'div-column-combobox-y--whisker--min', 'y');
    await sp.clearPanelColumn(page, 'prop-y-whisker-max', 'div-column-combobox-y--whisker--max', 'y');
    expect(await settledCanvasDiff(page, 'whisker-base')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 4 — Context menu', async () => {
    const errBefore = errCount();
    await sp.openPlotContextMenu(page);

    const entries = await sp.menuLeafNames(page);
    expect(entries).toContain('div-Reset-View');
    expect(entries).toContain('div-Lasso-Tool');
    expect(entries).toContain('div-Tools');
    expect(entries).toContain('div-Properties...');

    await page.keyboard.press('Escape');
    await page.locator('.d4-menu-popup').first().waitFor({state: 'detached', timeout: 8000});
    expect(await page.locator('.d4-menu-popup').count()).toBe(0);
    expect(await viewerAlive(page)).toBe(true);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 5 — Title and description', async () => {
    const errBefore = errCount();

    const viewerText = () => page.evaluate(() => {
      const root = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot').root as HTMLElement;
      const panel = root.closest('.panel-base') as HTMLElement;
      return {
        inViewer: (root.innerText ?? '').replace(/\s+/g, ' ').trim(),
        onViewer: (panel.innerText ?? '').replace(/\s+/g, ' ').trim(),
      };
    });

    const before = await viewerText();
    expect(before.onViewer).not.toContain(TITLE_TEXT);
    expect(before.inViewer).not.toContain(DESCRIPTION_TEXT);

    await sp.setTextProp(page, 'prop-title', 'description', TITLE_TEXT);
    expect((await v.pollValue(viewerText, (t) => t.onViewer.includes(TITLE_TEXT), 2000, 50)).onViewer).toContain(TITLE_TEXT);

    await sp.setTextProp(page, 'prop-description', 'description', DESCRIPTION_TEXT);
    const withBoth = await v.pollValue(viewerText, (t) => t.inViewer.includes(DESCRIPTION_TEXT), 2000, 50);
    expect(withBoth.onViewer).toContain(TITLE_TEXT);
    expect(withBoth.inViewer).toContain(DESCRIPTION_TEXT);

    await sp.setTextProp(page, 'prop-title', 'description', '');
    await sp.setTextProp(page, 'prop-description', 'description', '');
    const cleared = await v.pollValue(viewerText,
      (t) => !t.onViewer.includes(TITLE_TEXT) && !t.inViewer.includes(DESCRIPTION_TEXT), 2000, 50);
    expect(cleared.onViewer).not.toContain(TITLE_TEXT);
    expect(cleared.inViewer).not.toContain(DESCRIPTION_TEXT);

    expect(errCount()).toBe(errBefore);
  });

  await softStep('Lines By overrides the color column when splitting connecting lines', async () => {
    const errBefore = errCount();

    await sp.openSettings(page);
    await sp.revealPropEditor(page, '[name="prop-lines-by"]', 'data');
    expect(await sp.rowOpacity(page, 'prop-lines-by')).toBe('0.5');
    await captureBaseline(page, 'lines-base');

    await sp.pickPanelColumn(page, 'prop-color', 'div-column-combobox-color', 'color', COLOR_COLUMN);

    const colorCategories = await categoryCount(page, COLOR_COLUMN);
    const linesByCategories = await categoryCount(page, LINES_BY_COLUMN);
    expect(linesByCategories).toBeGreaterThan(1);
    expect(colorCategories).toBeGreaterThan(linesByCategories);
    await settledCanvasDiff(page, 'lines-base', true);
    expect(await captureCanvas(page, 'lines-nolines')).toBe(true);

    await sp.pickPanelColumn(page, 'prop-lines-order', 'div-column-combobox-lines--order', 'data', LINES_ORDER_COLUMN);
    expect(await v.pollValue(() => sp.rowOpacity(page, 'prop-lines-by'), (o) => o === '1', 2000, 50)).toBe('1');
    expect(await settledCanvasDiff(page, 'lines-nolines', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    expect(await captureCanvas(page, 'lines-color-split')).toBe(true);

    await sp.pickPanelColumn(page, 'prop-lines-by', 'div-column-combobox-lines--by', 'data', COLOR_COLUMN);
    expect(await settledCanvasDiff(page, 'lines-color-split')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    await sp.pickPanelColumn(page, 'prop-lines-by', 'div-column-combobox-lines--by', 'data', LINES_BY_COLUMN);
    expect(await settledCanvasDiff(page, 'lines-color-split', true)).toBeGreaterThanOrEqual(CANVAS_CHANGE_MIN);

    await sp.clearPanelColumn(page, 'prop-lines-by', 'div-column-combobox-lines--by', 'data');
    expect(await settledCanvasDiff(page, 'lines-color-split')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    await sp.clearPanelColumn(page, 'prop-lines-order', 'div-column-combobox-lines--order', 'data');
    await sp.clearPanelColumn(page, 'prop-color', 'div-column-combobox-color', 'color');
    expect(await settledCanvasDiff(page, 'lines-base')).toBeLessThanOrEqual(CANVAS_RESTORE_MAX);

    expect(errCount()).toBe(errBefore);
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
