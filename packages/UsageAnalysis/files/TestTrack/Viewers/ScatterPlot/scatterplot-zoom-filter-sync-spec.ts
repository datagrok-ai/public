/* ---
realizes: [scatterplot.cp.zoom-filter-sync, viewers.scatter-plot]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import * as sp from './scatterplot-shared';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const demogPath = 'System:DemoFiles/demog.csv';
const spgiPath = 'System:AppData/Chem/tests/spgi-100.csv';
const PROBE_X = 'ZZ_LOG_X_PROBE';
const PROBE_Y = 'ZZ_LOG_Y_PROBE';
const SPGI_X = 'CAST Idea ID';
const SPGI_Y = 'Idea ID';
const JITTER_X = 16;
const JITTER_Y = 17;

const propertyRowNames = (page: Page) => page.evaluate(() =>
  [...document.querySelectorAll('tr.property-grid-item[name]')]
    .map((e) => e.getAttribute('name') as string));

async function setFilterOutInvalid(page: Page, on: boolean): Promise<void> {
  await sp.openSettings(page);
  await sp.revealPropEditor(page, '[name="prop-view-filter-out-invalid"]', 'data');
  const box = page.locator('input[name="prop-view-filter-out-invalid"]');
  await box.scrollIntoViewIfNeeded();
  if (await box.isChecked() !== on) {
    await box.click();
    await v.pollValue(() => box.isChecked(), (c) => c === on, 3000, 50);
  }
  expect(await box.isChecked()).toBe(on);
  await v.waitForViewerRendered(page, sp.SP_TYPE, 500);
}

const setZoomAndFilter = (page: Page, value: string) =>
  sp.setChoiceProp(page, 'prop-zoom-and-filter', 'data', value, 'zoomAndFilter');

const X_AXIS_TYPE_MENU = ['div-Properties...', 'div-Properties...---X-Axis',
  'div-Properties...---X-Axis---X-Axis-Type'];
const xAxisTypeLeaf = (choice: 'Linear' | 'Logarithmic') =>
  `div-Properties...---X-Axis---X-Axis-Type---${choice}`;

async function wheelZoomIn(page: Page, steps = 1): Promise<void> {
  const r = await sp.canvasRect(page);
  await page.mouse.move(r.x + r.width / 2, r.y + r.height / 2);
  for (let i = 0; i < steps; i++) {
    await page.mouse.wheel(0, -300);
    await v.waitForViewerRendered(page, sp.SP_TYPE, 600);
  }
}

const rowCount = (page: Page) => page.evaluate(() => grok.shell.tv.dataFrame.rowCount as number);

const VIEWPORT_TOLERANCE = 0.02;

function expectSameViewport(actual: sp.Rect, reference: sp.Rect): void {
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
  const s = grok.shell.tv.viewers.find((x: any) => x.type === 'Scatter plot') as any;
  return {
    x: s.props.xColumnName, y: s.props.yColumnName,
    xAxisType: s.props.xAxisType, yAxisType: s.props.yAxisType,
    zoomAndFilter: s.props.zoomAndFilter,
    jitter: s.props.jitterSize, jitterY: s.props.jitterSizeY,
  };
});

const filteredIndicator = (page: Page) => page.evaluate(() => {
  const el = document.querySelector('[name="span-filtered"]') as HTMLElement | null;
  return {present: !!el, text: el ? (el.innerText ?? '').trim() : null};
});

// The status-bar indicator is written a frame after the row count settles, so a one-shot
// read races it — as `present` and as text, since it is built once with the count in it.
const settledIndicator = (page: Page, present: boolean) =>
  v.pollValue(() => filteredIndicator(page), (r) => r.present === present, 2000, 50);

test('Scatter Plot — Zoom and Filter Synchronization', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const errors = sp.trackErrors(page);
  const errCount = errors.count;

  await openDatagrok(page);
  await v.openTable(page, {path: demogPath, semTypeTimeoutMs: 3000});
  await sp.addScatterPlot(page);
  await sp.pickOnViewer(page, 'x', 'WEIGHT');
  await sp.pickOnViewer(page, 'y', 'HEIGHT');
  const fullRowCount = await rowCount(page);
  expect(fullRowCount).toBeGreaterThan(0);

  await softStep('Zoom drives the table filter, reset restores it', async () => {
    const errBefore = errCount();

    await sp.openSettings(page);
    await sp.revealPropEditor(page, '[name="prop-view-zoom-and-filter"]', 'data');
    if (await sp.readProp(page, 'zoomAndFilter') !== 'filter by zoom')
      await setZoomAndFilter(page, 'filter by zoom');
    expect(await sp.propCellText(page, 'prop-view-zoom-and-filter')).toBe('filter by zoom');

    expect(await sp.filterHeld(page, 500)).toBe(fullRowCount);

    await wheelZoomIn(page);
    const firstStep = await sp.filterMoved(page, fullRowCount);
    expect(firstStep).toBeLessThan(fullRowCount);

    await wheelZoomIn(page);
    const secondStep = await sp.filterMoved(page, firstStep);
    expect(secondStep).toBeLessThan(firstStep);

    await sp.resetView(page);
    expect(await sp.filterMoved(page, secondStep)).toBe(fullRowCount);

    expect((await viewerProps(page)).zoomAndFilter).toBe('filter by zoom');
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Pack and Zoom keeps the zoom filtering resettable', async () => {
    const errBefore = errCount();

    await wheelZoomIn(page, 2);
    const zoomed = await sp.filterMoved(page, fullRowCount);
    expect(zoomed).toBeLessThan(fullRowCount);

    await setZoomAndFilter(page, 'pack and zoom by filter');
    expect(await sp.propCellText(page, 'prop-view-zoom-and-filter')).toBe('pack and zoom by filter');

    const beforeReset = await sp.filterCount(page);
    await sp.resetView(page);
    const afterReset = beforeReset === fullRowCount
      ? await sp.filterHeld(page)
      : await sp.filterMoved(page, beforeReset);
    expect(afterReset).toBe(fullRowCount);
    expect(errCount()).toBe(errBefore);

    await setZoomAndFilter(page, 'filter by zoom');
    await sp.resetView(page);
    expect(await sp.filterHeld(page)).toBe(fullRowCount);
    expect((await viewerProps(page)).zoomAndFilter).toBe('filter by zoom');
  });

  await softStep('External filtering drives the viewport in zoom by filter mode', async () => {
    const errBefore = errCount();

    await setZoomAndFilter(page, 'zoom by filter');
    expect((await viewerProps(page)).zoomAndFilter).toBe('zoom by filter');

    const baseline = await sp.viewport(page);
    expect(baseline.width).toBeGreaterThan(0);
    expect(baseline.height).toBeGreaterThan(0);
    expect(await sp.filterHeld(page)).toBe(fullRowCount);

    await sp.openFilterPanel(page);
    let narrowedBand = fullRowCount;
    let fitted = baseline;
    try {
      const band = await applyRangeFilter(page, 'WEIGHT', 0.4, 0.6);
      expect(band.from).toBeGreaterThan(band.min);
      expect(band.to).toBeLessThan(band.max);

      const narrowed = await sp.filterMoved(page, fullRowCount);
      narrowedBand = narrowed;
      expect(narrowed).toBeGreaterThan(0);
      expect(narrowed).toBeLessThan(fullRowCount);

      fitted = await sp.viewportMoved(page, baseline);
      expect(fitted.width).toBeGreaterThan(0);
      expect(fitted.width).toBeLessThan(baseline.width * (1 - VIEWPORT_TOLERANCE));
      expect(errCount()).toBe(errBefore);
    } finally {
      await v.resetFilters(page);
      await sp.closeFilterPanel(page);
    }

    expect(await sp.filterMoved(page, narrowedBand)).toBe(fullRowCount);
    expectSameViewport(await sp.viewportMoved(page, fitted), baseline);

    await setZoomAndFilter(page, 'filter by zoom');
    expect((await viewerProps(page)).zoomAndFilter).toBe('filter by zoom');
  });

  await softStep('Filter Panel reset clears the scatter plot\'s contribution', async () => {
    const errBefore = errCount();
    await sp.openFilterPanel(page);
    // opening the panel clears the plot's own zoom contribution a beat AFTER the panel itself
    // appears; zooming into that window makes the count dip and bounce back
    expect(await sp.filterHeld(page)).toBe(fullRowCount);

    await wheelZoomIn(page, 2);
    const zoomed = await sp.filterMoved(page, fullRowCount);
    expect(zoomed).toBeLessThan(fullRowCount);

    const reported = await settledIndicator(page, true);
    expect({
      present: reported.present,
      mode: (await viewerProps(page)).zoomAndFilter,
      filtered: zoomed < fullRowCount,
    }).toEqual({present: true, mode: 'filter by zoom', filtered: true});
    expect(reported.text).toContain(String(zoomed));

    await page.locator('[name="viewer-Filters"] [name="icon-arrow-rotate-left"]').click();
    expect(await sp.filterMoved(page, zoomed)).toBe(fullRowCount);

    expect((await settledIndicator(page, false)).present).toBe(false);
    expect(errCount()).toBe(errBefore);

    await sp.closeFilterPanel(page);
    await sp.resetView(page);
    expect(await sp.filterHeld(page)).toBe(fullRowCount);
  });

  await softStep('Axis type switch on a datetime axis keeps the applied filter', async () => {
    const errBefore = errCount();
    await sp.pickOnViewer(page, 'x', 'STARTED');
    expect((await viewerProps(page)).x).toBe('STARTED');

    await sp.openFilterPanel(page);
    await v.applyCategoricalFilter(page, 'SEX', ['F'], 600);
    const baseline = await sp.filterMoved(page, fullRowCount);
    expect(baseline).toBeLessThan(fullRowCount);
    expect(baseline).toBeGreaterThan(0);

    await sp.openSettings(page);
    await sp.revealPropEditor(page, '[name="prop-view-x-axis-type"]', 'x-axis');
    expect(await v.pollValue(() => sp.rowOpacity(page, 'prop-x-axis-type'), (o) => o === '0.5', 2000, 50)).toBe('0.5');

    expect((await viewerProps(page)).xAxisType).toBe('linear');
    await sp.clickContextMenuLeaf(page, [...X_AXIS_TYPE_MENU, xAxisTypeLeaf('Logarithmic')]);
    expect(await sp.propIs(page, 'xAxisType', 'logarithmic')).toBe('logarithmic');

    expect(await sp.filterHeld(page)).toBe(baseline);
    expect(errCount()).toBe(errBefore);

    await sp.clickContextMenuLeaf(page, [...X_AXIS_TYPE_MENU, xAxisTypeLeaf('Linear')]);
    expect(await sp.propIs(page, 'xAxisType', 'linear')).toBe('linear');
    await v.resetFilters(page);
    await sp.closeFilterPanel(page);
    await sp.pickOnViewer(page, 'x', 'WEIGHT');
    await sp.resetView(page);
    expect((await viewerProps(page)).x).toBe('WEIGHT');
    expect(await sp.filterMoved(page, baseline)).toBe(fullRowCount);
  });

  await softStep('Filter Out Invalid removes the rows a logarithmic axis cannot draw', async () => {
    const errBefore = errCount();

    const fixture = await page.evaluate(({xName, yName}) => {
      const df = grok.shell.tv.dataFrame;
      for (const n of [xName, yName]) if (df.col(n)) df.columns.remove(n);
      const xs: number[] = [];
      const ys: number[] = [];
      for (let i = 0; i < df.rowCount; i++) {
        xs.push(i + 1);
        ys.push(i % 4 === 0 ? -(i + 1) : (i + 1));
      }
      df.columns.add(DG.Column.fromList('double', xName, xs));
      df.columns.add(DG.Column.fromList('double', yName, ys));
      const xc = df.col(xName);
      const yc = df.col(yName);
      let nonPositive = 0;
      let xUndrawable = 0;
      for (let i = 0; i < df.rowCount; i++) {
        if (!(yc.get(i) > 0)) nonPositive++;
        if (xc.isNone(i) || !isFinite(xc.get(i))) xUndrawable++;
      }
      return {nonPositive, xUndrawable, rowCount: df.rowCount as number};
    }, {xName: PROBE_X, yName: PROBE_Y});
    await v.waitForViewerRendered(page, sp.SP_TYPE, 1500);
    expect(fixture.rowCount).toBe(fullRowCount);
    expect(fixture.xUndrawable).toBe(0);
    expect(fixture.nonPositive).toBeGreaterThan(0);
    expect(fixture.nonPositive).toBeLessThan(fullRowCount);

    try {
      await sp.pickOnViewer(page, 'x', PROBE_X);
      await sp.pickOnViewer(page, 'y', PROBE_Y);
      const axes = await viewerProps(page);
      expect(axes.x).toBe(PROBE_X);
      expect(axes.y).toBe(PROBE_Y);

      await sp.openSettings(page);
      await sp.revealPropEditor(page, '[name="prop-view-filter-out-invalid"]', 'data');
      expect(await propertyRowNames(page)).toContain('prop-filter-out-invalid');

      await setFilterOutInvalid(page, false);
      expect(await sp.filterHeld(page)).toBe(fullRowCount);

      await sp.setChoiceProp(page, 'prop-y-axis-type', 'y-axis', 'logarithmic');
      expect((await viewerProps(page)).yAxisType).toBe('logarithmic');
      expect(await sp.filterHeld(page)).toBe(fullRowCount);

      await setFilterOutInvalid(page, true);
      expect(await sp.filterMoved(page, fullRowCount, 8000)).toBe(fullRowCount - fixture.nonPositive);

      await setFilterOutInvalid(page, false);
      expect(await sp.filterMoved(page, fullRowCount - fixture.nonPositive, 8000)).toBe(fullRowCount);

      await sp.setChoiceProp(page, 'prop-y-axis-type', 'y-axis', 'linear');
      await setFilterOutInvalid(page, true);
      expect(await sp.filterHeld(page)).toBe(fullRowCount);
      expect(errCount()).toBe(errBefore);
    } finally {
      await setFilterOutInvalid(page, false);
      await sp.setChoiceProp(page, 'prop-y-axis-type', 'y-axis', 'linear');
      await sp.pickOnViewer(page, 'x', 'WEIGHT');
      await sp.pickOnViewer(page, 'y', 'HEIGHT');
      await page.evaluate((names: string[]) => {
        const df = grok.shell.tv.dataFrame;
        for (const n of names) if (df.col(n)) df.columns.remove(n);
      }, [PROBE_X, PROBE_Y]);
      await v.waitForViewerRendered(page, sp.SP_TYPE, 1500);
    }
    const reverted = await viewerProps(page);
    expect(reverted.x).toBe('WEIGHT');
    expect(reverted.y).toBe('HEIGHT');
    expect(reverted.yAxisType).toBe('linear');
    expect(await sp.filterHeld(page)).toBe(fullRowCount);
  });

  await softStep('Large jitter with a logarithmic axis does not filter rows', async () => {
    const errBefore = errCount();
    await sp.addTableView(page, spgiPath, SPGI_X);
    await sp.addScatterPlot(page);

    await sp.pickOnViewer(page, 'x', SPGI_X);
    await sp.pickOnViewer(page, 'y', SPGI_Y);
    const axes = await viewerProps(page);
    expect(axes.x).toBe(SPGI_X);
    expect(axes.y).toBe(SPGI_Y);
    const yMin = await page.evaluate((c: string) =>
      grok.shell.tv.dataFrame.col(c).stats.min as number, SPGI_Y);
    expect(yMin).toBeGreaterThan(0);

    const spgiRowCount = await rowCount(page);
    expect(spgiRowCount).toBeGreaterThan(0);
    expect(await sp.filterHeld(page, 500)).toBe(spgiRowCount);

    await sp.setNumericProp(page, 'prop-jitter-size', 'marker', JITTER_X);
    await sp.setNumericProp(page, 'prop-jitter-size-y', 'marker', JITTER_Y);
    const jittered = await viewerProps(page);
    expect(jittered.jitter).toBe(JITTER_X);
    expect(jittered.jitterY).toBe(JITTER_Y);

    await sp.setChoiceProp(page, 'prop-y-axis-type', 'y-axis', 'logarithmic');
    expect((await viewerProps(page)).yAxisType).toBe('logarithmic');

    expect(await sp.filterHeld(page)).toBe(spgiRowCount);
    expect(errCount()).toBe(errBefore);

    await sp.setChoiceProp(page, 'prop-y-axis-type', 'y-axis', 'linear');
    await sp.setNumericProp(page, 'prop-jitter-size', 'marker', 0);
    await sp.setNumericProp(page, 'prop-jitter-size-y', 'marker', 0);
    const back = await viewerProps(page);
    expect(back.yAxisType).toBe('linear');
    expect(back.jitter).toBe(0);
    expect(back.jitterY).toBe(0);
    expect(await sp.filterHeld(page, 500)).toBe(spgiRowCount);
    await page.evaluate(() => grok.shell.v.close());
    await v.pollValue(() => page.evaluate(() => grok.shell.tv?.dataFrame?.rowCount ?? -1),
      (n) => n === fullRowCount, 1500, 50);
  });

  await v.cleanupShell(page);
  v.finishSpec();
});
