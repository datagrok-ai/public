/* ---
realizes: [scatterplot.cp.zoom-filter-sync, viewers.scatter-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const PROBE_X = 'ZZ_LOG_X_PROBE';
const PROBE_Y = 'ZZ_LOG_Y_PROBE';
const SPGI_X = 'CAST Idea ID';
const SPGI_Y = 'Idea ID';
const JITTER_X = 16;
const JITTER_Y = 17;

const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Failed to connect to Claude runtime/.test(text) ||
  /powerPreference option is currently ignored/.test(text) ||
  /willReadFrequently/.test(text);

const isBenignError = (text: string) => isAmbientError(text);

const canvasRect = (page: Page) => page.evaluate(() => {
  const root = [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .find((e) => !e.closest('.d4-dialog'))!;
  const r = root.querySelector('canvas[name="canvas"]')!.getBoundingClientRect();
  return {x: r.x, y: r.y, width: r.width, height: r.height};
});

const waitPlotCanvas = (page: Page) => page.waitForFunction(() => {
  return [...document.querySelectorAll('[name="viewer-Scatter-plot"]')]
    .filter((e) => !e.closest('.d4-dialog'))
    .some((e) => {
      const c = e.querySelector('canvas[name="canvas"]');
      return !!c && c.getBoundingClientRect().width > 0 && c.getBoundingClientRect().height > 0;
    });
}, null, {timeout: 20000});

const pickOnViewer = (page: Page, role: string, column: string) =>
  v.pickColumnViaSelectorTrusted(page, {role, columnName: column});

async function openSettings(page: Page): Promise<void> {
  for (let i = 0; i < 4; i++) {
    // Never return early on an already-built panel: once a second view is open it can
    // still be showing the FIRST view's scatter plot, and every edit below would land
    // there while the assertions read grok.shell.tv — the jitter step failed that way.
    await v.openViewerGear(page, 'Scatter plot');

    await v.pollValue(() => page.evaluate(() => document.querySelectorAll('[name^="prop-"]').length),
      (n) => n > 0, 1000, 200);
    if (await page.evaluate(() => !!document.querySelector('[name="prop-category-data"]'))) return;
  }
  throw new Error('the scatter plot settings panel did not build');
}

async function revealPropEditor(
  page: Page, editorSelector: string, category: string,
): Promise<void> {
  for (let i = 0; i < 8; i++) {
    const ready = await page.evaluate((sel: string) => {
      const el = document.querySelector(sel) as HTMLElement | null;
      if (!el || !el.offsetParent) return false;
      const b = el.getBoundingClientRect();
      return b.width > 0 && b.height > 0;
    }, editorSelector);
    if (ready) return;
    const header = page.locator(`[name="prop-category-${category}"]`);
    if (await header.count() > 0 && await header.isVisible()) await header.click();

    await v.pollValue(() => page.evaluate(() => document.querySelectorAll('[name^="prop-"]').length),
      (n) => n > 0, 800, 150);
  }
  throw new Error(`property editor ${editorSelector} never became reachable`);
}

const rowProp = (rowName: string) =>
  rowName.replace(/^prop-/, '').replace(/-([a-z])/g, (_, c: string) => c.toUpperCase());

const viewerProp = (page: Page, prop: string) => page.evaluate((p: string) => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return sp ? sp.props[p] : null;
}, prop);

async function setChoiceProp(
  page: Page, rowName: string, viewCell: string, category: string, value: string,
): Promise<void> {
  await openSettings(page);
  await revealPropEditor(page, `[name="${viewCell}"]`, category);
  const cell = page.locator(`[name="${viewCell}"]`);
  await cell.scrollIntoViewIfNeeded();
  await cell.click();
  await page.locator(`[name="${rowName}"] select.property-grid-item-editor-spinner`).selectOption(value);
  await v.pollValue(() => propCellText(page, viewCell), (t) => t === value, 4000, 100);
  await v.waitForViewerRendered(page, 'Scatter plot', 200);
}

async function setSliderProp(
  page: Page, rowName: string, category: string, value: number,
): Promise<void> {
  await openSettings(page);
  const box = `[name="${rowName}"] input.property-grid-slider-textbox`;
  await revealPropEditor(page, box, category);
  const locator = page.locator(box);
  await locator.scrollIntoViewIfNeeded();
  await locator.click();
  await page.keyboard.press('Control+A');
  await page.keyboard.type(String(value));
  await page.keyboard.press('Enter');
  await v.pollValue(() => viewerProp(page, rowProp(rowName)), (x) => x === value, 4000, 100);
  await v.waitForViewerRendered(page, 'Scatter plot', 200);
}

const propertyRowNames = (page: Page) => page.evaluate(() =>
  [...document.querySelectorAll('tr.property-grid-item[name]')]
    .map((e) => e.getAttribute('name') as string));

async function setFilterOutInvalid(page: Page, on: boolean): Promise<void> {
  await openSettings(page);
  await revealPropEditor(page, '[name="prop-view-filter-out-invalid"]', 'data');
  const box = page.locator('input[name="prop-view-filter-out-invalid"]');
  await box.scrollIntoViewIfNeeded();
  if (await box.isChecked() !== on) {
    await box.click();
    await v.pollValue(() => box.isChecked(), (c) => c === on, 3000, 100);
  }
  expect(await box.isChecked()).toBe(on);
}

const propCellText = (page: Page, viewCell: string) => page.evaluate((n: string) =>
  (document.querySelector(`[name="${n}"]`)?.textContent ?? '').trim(), viewCell);

const rowOpacity = (page: Page, rowName: string) => page.evaluate((n: string) =>
  (document.querySelector(`[name="${n}"]`) as HTMLElement | null)?.style.opacity ?? null, rowName);

async function openPlotContextMenu(page: Page): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.click(r.x + r.width / 2, r.y + r.height / 2, {button: 'right'});
  await page.locator('.d4-menu-popup').last().waitFor({timeout: 8000});

  await v.pollStable(() => page.evaluate(() => document.querySelectorAll('.d4-menu-popup').length),
    (a, b) => a === b, 500, 100);
}

const dismissPlotMenu = (page: Page, capMs: number) =>
  v.pollValue(() => page.locator('.d4-menu-popup').count(), (n) => n === 0, capMs, 50);

async function clickContextMenuLeaf(page: Page, groups: string[], leaf: string): Promise<void> {
  await openPlotContextMenu(page);
  for (let i = 0; i < groups.length; i++) {
    await page.locator(`.d4-menu-popup [name="${groups[i]}"]`).last().hover();
    const next = i + 1 < groups.length ? groups[i + 1] : leaf;
    await page.locator(`.d4-menu-popup [name="${next}"]`).last().waitFor({timeout: 8000});
  }
  await page.locator(`.d4-menu-popup [name="${leaf}"]`).last().click();
  await v.waitForViewerRendered(page, 'Scatter plot', 1000);
  await page.keyboard.press('Escape');
  await dismissPlotMenu(page, 500);
}

const X_AXIS_TYPE_MENU = ['div-Properties...', 'div-Properties...---X-Axis',
  'div-Properties...---X-Axis---X-Axis-Type'];
const xAxisTypeLeaf = (choice: 'Linear' | 'Logarithmic') =>
  `div-Properties...---X-Axis---X-Axis-Type---${choice}`;

async function wheelZoomIn(page: Page, steps = 1): Promise<void> {
  const r = await canvasRect(page);
  await page.mouse.move(r.x + r.width / 2, r.y + r.height / 2);
  for (let i = 0; i < steps; i++) {
    await page.mouse.wheel(0, -300);

    await v.waitForViewerRendered(page, 'Scatter plot', 600);
  }
}

async function resetView(page: Page): Promise<void> {
  await openPlotContextMenu(page);
  await page.locator('.d4-menu-popup [name="div-Reset-View"]').last().click();
  await v.waitForViewerRendered(page, 'Scatter plot', 800);
  await page.keyboard.press('Escape');
  await dismissPlotMenu(page, 400);
}

const filterCount = (page: Page) =>
  page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount as number);

async function settledFilterCountAfterChange(page: Page, from: number): Promise<number> {
  let prev = -1;
  let stable = 0;
  for (let i = 0; i < 70; i++) {
    const c = await filterCount(page);
    if (c === prev) {
      stable++;
      if (c !== from && stable >= 4) return c;
    } else {
      stable = 0;
      prev = c;
    }

    await page.waitForTimeout(250);
  }
  return prev;
}

async function settledFilterCountUnchanged(page: Page, floorMs = 2000): Promise<number> {

  await page.waitForTimeout(floorMs);
  let last = -1;
  let stable = 0;
  for (let i = 0; i < 40; i++) {
    const c = await filterCount(page);
    if (c === last) {
      stable++;
      if (stable >= 3) return c;
    } else {
      stable = 0;
      last = c;
    }

    await page.waitForTimeout(400);
  }
  return last;
}

const rowCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.rowCount as number);

interface Rect {x: number; y: number; width: number; height: number}

const viewport = (page: Page): Promise<Rect> => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  const vp = sp.viewport;
  return {x: vp.x, y: vp.y, width: vp.width, height: vp.height};
});

const VIEWPORT_TOLERANCE = 0.02;

function expectSameViewport(actual: Rect, reference: Rect): void {
  expect(Math.abs(actual.width - reference.width)).toBeLessThan(reference.width * VIEWPORT_TOLERANCE);
  expect(Math.abs(actual.height - reference.height)).toBeLessThan(reference.height * VIEWPORT_TOLERANCE);
  expect(Math.abs(actual.x - reference.x)).toBeLessThan(reference.width * VIEWPORT_TOLERANCE);
  expect(Math.abs(actual.y - reference.y)).toBeLessThan(reference.height * VIEWPORT_TOLERANCE);
}

async function applyRangeFilter(
  page: Page, column: string, loFraction: number, hiFraction: number,
): Promise<{min: number; max: number; from: number; to: number}> {
  const band = await page.evaluate(({col, lo, hi}: {col: string; lo: number; hi: number}) => {
    const stats = grok.shell.tv.dataFrame.col(col).stats;
    const min = stats.min as number;
    const max = stats.max as number;
    return {min, max, from: min + (max - min) * lo, to: min + (max - min) * hi};
  }, {col: column, lo: loFraction, hi: hiFraction});

  await v.applyNumericFilter(page, column, band.from, band.to, 1500);
  return band;
}

const viewerProps = (page: Page) => page.evaluate(() => {
  const sp = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return {
    x: sp.props.xColumnName, y: sp.props.yColumnName,
    xAxisType: sp.props.xAxisType, yAxisType: sp.props.yAxisType,
    zoomAndFilter: sp.props.zoomAndFilter,
    jitter: sp.props.jitterSize, jitterY: sp.props.jitterSizeY,
  };
});

const filteredIndicator = (page: Page) => page.evaluate(() => {
  const el = document.querySelector('[name="span-filtered"]') as HTMLElement | null;
  return {present: !!el, text: el ? (el.innerText ?? '').trim() : null};
});

// The status-bar indicator is written a frame after the row count settles, so a one-shot
// read races it — as `present` and as text, since it is built once with the count in it.
const settledIndicator = (page: Page, present: boolean) =>
  v.pollValue(() => filteredIndicator(page), (r) => r.present === present, 2000, 100);

async function closeFilterPanel(page: Page): Promise<void> {
  await page.evaluate(() => {
    const f = grok.shell.tv.viewers.find((x: any) => x.type === 'Filters');
    if (f) f.close();
  });
  await v.pollValue(() => page.evaluate(() =>
    !grok.shell.tv.viewers.find((x: any) => x.type === 'Filters')), (gone) => gone, 1000, 100);
}

test('Scatter Plot — Zoom and Filter Synchronization', async ({page}: {page: Page}) => {
  test.setTimeout(1_500_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await loginToDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
  await waitPlotCanvas(page);
  await pickOnViewer(page, 'x', 'WEIGHT');
  await pickOnViewer(page, 'y', 'HEIGHT');
  const fullRowCount = await rowCount(page);
  expect(fullRowCount).toBeGreaterThan(0);

  await softStep('Zoom drives the table filter, reset restores it', async () => {
    const errBefore = errCount();

    await openSettings(page);
    await revealPropEditor(page, '[name="prop-view-zoom-and-filter"]', 'data');
    if (await propCellText(page, 'prop-view-zoom-and-filter') !== 'filter by zoom')
      await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data', 'filter by zoom');
    expect(await propCellText(page, 'prop-view-zoom-and-filter')).toBe('filter by zoom');

    expect(await settledFilterCountUnchanged(page, 500)).toBe(fullRowCount);

    await wheelZoomIn(page);
    const firstStep = await settledFilterCountAfterChange(page, fullRowCount);
    expect(firstStep).toBeLessThan(fullRowCount);

    await wheelZoomIn(page);
    const secondStep = await settledFilterCountAfterChange(page, firstStep);
    expect(secondStep).toBeLessThan(firstStep);

    await resetView(page);
    expect(await settledFilterCountAfterChange(page, secondStep)).toBe(fullRowCount);

    expect((await viewerProps(page)).zoomAndFilter).toBe('filter by zoom');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Pack and Zoom keeps the zoom filtering resettable', async () => {
    const errBefore = errCount();

    await wheelZoomIn(page, 2);
    const zoomed = await settledFilterCountAfterChange(page, fullRowCount);
    expect(zoomed).toBeLessThan(fullRowCount);

    await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data',
      'pack and zoom by filter');
    expect(await propCellText(page, 'prop-view-zoom-and-filter')).toBe('pack and zoom by filter');

    const beforeReset = await filterCount(page);
    await resetView(page);
    const afterReset = beforeReset === fullRowCount
      ? await settledFilterCountUnchanged(page)
      : await settledFilterCountAfterChange(page, beforeReset);
    expect(afterReset).toBe(fullRowCount);
    expect(errCount()).toBe(errBefore);

    await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data', 'filter by zoom');
    await resetView(page);
    expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);
    expect((await viewerProps(page)).zoomAndFilter).toBe('filter by zoom');
  });

  await softStep('External filtering drives the viewport in zoom by filter mode', async () => {
    const errBefore = errCount();

    await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data',
      'zoom by filter');
    expect((await viewerProps(page)).zoomAndFilter).toBe('zoom by filter');

    const baseline = await viewport(page);
    expect(baseline.width).toBeGreaterThan(0);
    expect(baseline.height).toBeGreaterThan(0);
    expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);

    await v.openFilterPanel(page);
    let narrowedBand = fullRowCount;
    try {

      const band = await applyRangeFilter(page, 'WEIGHT', 0.4, 0.6);
      expect(band.from).toBeGreaterThan(band.min);
      expect(band.to).toBeLessThan(band.max);

      const narrowed = await settledFilterCountAfterChange(page, fullRowCount);
      narrowedBand = narrowed;
      expect(narrowed).toBeGreaterThan(0);
      expect(narrowed).toBeLessThan(fullRowCount);

      const fitted = await viewport(page);
      expect(fitted.width).toBeGreaterThan(0);

      expect(fitted.width).toBeLessThan(baseline.width * (1 - VIEWPORT_TOLERANCE));
      expect(errCount()).toBe(errBefore);
    } finally {
      await v.resetFilters(page);
    }

    expect(await settledFilterCountAfterChange(page, narrowedBand)).toBe(fullRowCount);
    expectSameViewport(await viewport(page), baseline);

    await setChoiceProp(page, 'prop-zoom-and-filter', 'prop-view-zoom-and-filter', 'data',
      'filter by zoom');
    expect((await viewerProps(page)).zoomAndFilter).toBe('filter by zoom');
    await closeFilterPanel(page);
  });

  await softStep('Filter Panel reset clears the scatter plot\'s contribution', async () => {
    const errBefore = errCount();
    await v.openFilterPanel(page);
    // opening the panel clears the plot's own zoom contribution a beat AFTER the panel
    // itself appears (4891 -> 5850 observed). Zooming into that window makes the count
    // dip and bounce back: the dip satisfies the "changed" wait, then the reset lands
    // and takes the indicator with it, so the assertion below sees no indicator at all
    expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);

    await wheelZoomIn(page, 2);
    const zoomed = await settledFilterCountAfterChange(page, fullRowCount);
    expect(zoomed).toBeLessThan(fullRowCount);

    const reported = await settledIndicator(page, true);
    // one consolidated read: when the indicator is missing, the mode and the count are
    // what say whether the zoom actually filtered or the step inherited the wrong state
    expect({
      present: reported.present,
      mode: (await viewerProps(page)).zoomAndFilter,
      filtered: zoomed < fullRowCount,
    }).toEqual({present: true, mode: 'filter by zoom', filtered: true});
    expect(reported.text).toContain(String(zoomed));

    await page.locator('[name="viewer-Filters"] [name="icon-arrow-rotate-left"]').click();
    expect(await settledFilterCountAfterChange(page, zoomed)).toBe(fullRowCount);

    expect((await settledIndicator(page, false)).present).toBe(false);
    expect(errCount()).toBe(errBefore);

    await closeFilterPanel(page);
    await resetView(page);
    expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);
  });

  await softStep('Axis type switch on a datetime axis keeps the applied filter', async () => {
    const errBefore = errCount();
    await pickOnViewer(page, 'x', 'STARTED');
    expect((await viewerProps(page)).x).toBe('STARTED');

    await v.openFilterPanel(page);
    await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
    const baseline = await settledFilterCountAfterChange(page, fullRowCount);
    expect(baseline).toBeLessThan(fullRowCount);
    expect(baseline).toBeGreaterThan(0);

    await openSettings(page);
    await revealPropEditor(page, '[name="prop-view-x-axis-type"]', 'x-axis');
    expect(await rowOpacity(page, 'prop-x-axis-type')).toBe('0.5');

    expect((await viewerProps(page)).xAxisType).toBe('linear');
    await clickContextMenuLeaf(page, X_AXIS_TYPE_MENU, xAxisTypeLeaf('Logarithmic'));

    expect((await viewerProps(page)).xAxisType).toBe('logarithmic');

    expect(await settledFilterCountUnchanged(page)).toBe(baseline);
    expect(errCount()).toBe(errBefore);

    await clickContextMenuLeaf(page, X_AXIS_TYPE_MENU, xAxisTypeLeaf('Linear'));
    expect((await viewerProps(page)).xAxisType).toBe('linear');
    await v.resetFilters(page);
    await closeFilterPanel(page);
    await pickOnViewer(page, 'x', 'WEIGHT');
    await resetView(page);
    const reverted = await viewerProps(page);
    expect(reverted.x).toBe('WEIGHT');
    expect(await settledFilterCountAfterChange(page, baseline)).toBe(fullRowCount);
  });

  await softStep('Filter Out Invalid removes the rows a logarithmic axis cannot draw', async () => {
    const errBefore = errCount();

    const fixture = await page.evaluate(({xName, yName}) => {
      const df = grok.shell.tv.dataFrame;
      for (const n of [xName, yName]) if (df.col(n)) df.columns.remove(n);
      const xc = df.columns.addNewFloat(xName);
      const yc = df.columns.addNewFloat(yName);
      for (let i = 0; i < df.rowCount; i++) {
        xc.set(i, i + 1);
        yc.set(i, i % 4 === 0 ? -(i + 1) : (i + 1));
      }
      let nonPositive = 0;
      let xUndrawable = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (!(yc.get(i) > 0)) nonPositive++;
        if (xc.isNone(i) || !isFinite(xc.get(i))) xUndrawable++;
      }
      return {nonPositive, xUndrawable, rowCount: df.rowCount as number};
    }, {xName: PROBE_X, yName: PROBE_Y});
    await v.waitForViewerRendered(page, 'Scatter plot', 1500);
    expect(fixture.rowCount).toBe(fullRowCount);
    expect(fixture.xUndrawable).toBe(0);
    expect(fixture.nonPositive).toBeGreaterThan(0);
    expect(fixture.nonPositive).toBeLessThan(fullRowCount);

    try {
      await pickOnViewer(page, 'x', PROBE_X);
      await pickOnViewer(page, 'y', PROBE_Y);
      const axes = await viewerProps(page);
      expect(axes.x).toBe(PROBE_X);
      expect(axes.y).toBe(PROBE_Y);

      await openSettings(page);
      await revealPropEditor(page, '[name="prop-view-filter-out-invalid"]', 'data');
      expect(await propertyRowNames(page)).toContain('prop-filter-out-invalid');

      await setFilterOutInvalid(page, false);
      expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);

      await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y-axis', 'logarithmic');
      expect((await viewerProps(page)).yAxisType).toBe('logarithmic');
      expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);

      await setFilterOutInvalid(page, true);
      expect(await settledFilterCountAfterChange(page, fullRowCount))
        .toBe(fullRowCount - fixture.nonPositive);

      await setFilterOutInvalid(page, false);
      expect(await settledFilterCountAfterChange(page, fullRowCount - fixture.nonPositive))
        .toBe(fullRowCount);

      await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y-axis', 'linear');
      await setFilterOutInvalid(page, true);
      expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);
      expect(errCount()).toBe(errBefore);
    } finally {
      await setFilterOutInvalid(page, false);
      await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y-axis', 'linear');
      await pickOnViewer(page, 'x', 'WEIGHT');
      await pickOnViewer(page, 'y', 'HEIGHT');
      await page.evaluate((names: string[]) => {
        const df = grok.shell.tv.dataFrame;
        for (const n of names) if (df.col(n)) df.columns.remove(n);
      }, [PROBE_X, PROBE_Y]);
      await v.waitForViewerRendered(page, 'Scatter plot', 1500);
    }
    const reverted = await viewerProps(page);
    expect(reverted.x).toBe('WEIGHT');
    expect(reverted.y).toBe('HEIGHT');
    expect(reverted.yAxisType).toBe('linear');
    expect(await settledFilterCountUnchanged(page)).toBe(fullRowCount);
  });

  await softStep('Large jitter with a logarithmic axis does not filter rows', async () => {
    const errBefore = errCount();
    await page.evaluate(async (path: string) => {
      const df = await grok.dapi.files.readCsv(path);
      grok.shell.addTableView(df);
      await new Promise((resolve) => {
        const s = df.onSemanticTypeDetected.subscribe(() => { s.unsubscribe(); resolve(null); });
        setTimeout(resolve, 6000);
      });
    }, spgiPath);
    await page.waitForFunction((col: string) =>
      !!grok.shell.tv?.dataFrame?.col(col), SPGI_X, {timeout: 30000});
    await v.addViewerByIcon(page, 'scatter-plot', 'Scatter-plot');
    await waitPlotCanvas(page);

    await pickOnViewer(page, 'x', SPGI_X);
    await pickOnViewer(page, 'y', SPGI_Y);
    const axes = await viewerProps(page);
    expect(axes.x).toBe(SPGI_X);
    expect(axes.y).toBe(SPGI_Y);
    const yMin = await page.evaluate((c: string) =>
      grok.shell.tv.dataFrame.col(c).stats.min as number, SPGI_Y);
    expect(yMin).toBeGreaterThan(0);

    const spgiRowCount = await rowCount(page);
    expect(spgiRowCount).toBeGreaterThan(0);
    expect(await settledFilterCountUnchanged(page, 500)).toBe(spgiRowCount);

    await setSliderProp(page, 'prop-jitter-size', 'marker', JITTER_X);
    await setSliderProp(page, 'prop-jitter-size-y', 'marker', JITTER_Y);
    const jittered = await viewerProps(page);
    expect(jittered.jitter).toBe(JITTER_X);
    expect(jittered.jitterY).toBe(JITTER_Y);

    await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y-axis', 'logarithmic');
    expect((await viewerProps(page)).yAxisType).toBe('logarithmic');

    expect(await settledFilterCountUnchanged(page)).toBe(spgiRowCount);
    expect(errCount()).toBe(errBefore);

    await setChoiceProp(page, 'prop-y-axis-type', 'prop-view-y-axis-type', 'y-axis', 'linear');
    await setSliderProp(page, 'prop-jitter-size', 'marker', 0);
    await setSliderProp(page, 'prop-jitter-size-y', 'marker', 0);
    const back = await viewerProps(page);
    expect(back.yAxisType).toBe('linear');
    expect(back.jitter).toBe(0);
    expect(back.jitterY).toBe(0);
    expect(await settledFilterCountUnchanged(page, 500)).toBe(spgiRowCount);
    await page.evaluate(() => grok.shell.v.close());
    await v.pollValue(() => page.evaluate(() => grok.shell.tv?.dataFrame?.rowCount ?? -1),
      (n) => n === fullRowCount, 1500, 100);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
