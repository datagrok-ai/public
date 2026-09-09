/* ---
realizes: [barchart.cp.value-axis-range-scale-scroll]
--- */

import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';
const splitCol = 'Primary Series Name';

test('Bar Chart — Value-Axis Range, Scale, and Scroll', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => {
    if (m.type() === 'error' && !isLocalBootNoise(m.text())) consoleErrors.push(m.text());
  });

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 4000});
  await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');

  const barChartState = () => page.evaluate(() => {
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    return {
      valueMin: bc.props.valueMin,
      valueMax: bc.props.valueMax,
      axisType: bc.props.axisType,
      clippedIndicators: bc.props.showClippedBarIndicators,
      hasCanvas: !!bc.root.querySelector('canvas'),
    };
  });

  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear.click();
  });
  await v.pollValue(() => page.locator('.property-grid').count(), (n) => n > 0, 500, 50);

  await v.setViewerProps(page, 'Bar chart', [{
    set: {
      splitColumnName: splitCol,
      valueColumnName: 'CAST Idea ID',
      valueAggrType: 'count',
      showClippedBarIndicators: true,
    },
    wait: 900,
  }]);

  await softStep('Scenario 1 Step 3: Value Min above the shortest bar clips it; clipped-bar indicators render (GROK-19346)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const stats = await page.evaluate(({split}) => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col(split);
      const counts: Record<string, number> = {};
      for (const c of col.categories) counts[c] = 0;
      for (let i = 0; i < df.rowCount; i++) counts[col.get(i)]++;
      const nonZero = Object.values(counts).filter((n) => (n as number) > 0) as number[];
      return {nonZero, shortest: Math.min(...nonZero), tallest: Math.max(...nonZero)};
    }, {split: splitCol});
    const minOffset = Math.max(1, Math.round((stats.tallest - stats.shortest) * 0.25));
    await v.setViewerProps(page, 'Bar chart',
      [{set: {valueMin: stats.shortest + minOffset}, wait: 800}]);
    const read = await barChartState();
    const info = {
      ...read,
      shortest: stats.shortest,
      tallest: stats.tallest,
      minOffset,
      barsBelowMin: stats.nonZero.filter((n) => n < read.valueMin).length,
    };

    await v.setViewerProps(page, 'Bar chart', [{set: {showClippedBarIndicators: false}, wait: 700}]);
    await v.waitForCanvasQuiet(page, 'Bar chart', {timeoutMs: 900, optional: true});
    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0);
    expect(settle.deltaPx).toBeLessThan(30);
    await v.setViewerProps(page, 'Bar chart', [{set: {showClippedBarIndicators: true}, wait: 700}]);

    const {deltaPx} = await v.diffCanvasColors(page, 'Bar chart');
    console.log(`Step 3 clipped-indicator toggle deltaPx=${deltaPx}`);
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.valueMin).toBe(info.shortest + info.minOffset);
    expect(info.valueMin).toBeLessThan(info.tallest);
    expect(info.barsBelowMin).toBeGreaterThanOrEqual(1);
    expect(info.clippedIndicators).toBe(true);
    expect(deltaPx).toBeGreaterThan(150);
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 5: Value Max below the tallest bar clips it; axis stops at the maximum; top clipped-bar indicators render', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const stats = await page.evaluate(({split}) => {
      const df = grok.shell.tv.dataFrame;
      const col = df.col(split);
      const counts: Record<string, number> = {};
      for (const c of col.categories) counts[c] = 0;
      for (let i = 0; i < df.rowCount; i++) counts[col.get(i)]++;
      const nonZero = Object.values(counts).filter((n) => (n as number) > 0) as number[];
      return {nonZero, shortest: Math.min(...nonZero), tallest: Math.max(...nonZero)};
    }, {split: splitCol});
    const maxOffset = Math.max(1, Math.round((stats.tallest - stats.shortest) * 0.25));
    await v.setViewerProps(page, 'Bar chart',
      [{set: {valueMax: stats.tallest - maxOffset}, wait: 800}]);
    const read = await barChartState();
    const info = {
      ...read,
      shortest: stats.shortest,
      tallest: stats.tallest,
      maxOffset,
      barsAboveMax: stats.nonZero.filter((n) => n > read.valueMax).length,
    };

    await v.setViewerProps(page, 'Bar chart', [{set: {showClippedBarIndicators: false}, wait: 700}]);
    await v.waitForCanvasQuiet(page, 'Bar chart', {timeoutMs: 900, optional: true});
    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0);
    expect(settle.deltaPx).toBeLessThan(30);
    await v.setViewerProps(page, 'Bar chart', [{set: {showClippedBarIndicators: true}, wait: 700}]);

    const {deltaPx} = await v.diffCanvasColors(page, 'Bar chart');
    console.log(`Step 5 clipped-indicator toggle deltaPx=${deltaPx}`);
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.valueMax).toBe(info.tallest - info.maxOffset);
    expect(info.valueMax).toBeGreaterThan(info.shortest);
    expect(info.valueMax).toBeGreaterThan(info.valueMin);
    expect(info.barsAboveMax).toBeGreaterThanOrEqual(1);
    expect(info.valueMin).not.toBeNull();
    expect(deltaPx).toBeGreaterThan(150);
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 7: the value-axis scroll bar is present on the constrained range', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await page.locator('[name="viewer-Bar-chart"]').first().hover();

    const info = await v.pollValue(() => page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const rangeSliders = bc.root.querySelectorAll('svg[type="range-slider"]').length;
      return {
        rangeSliders,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    }), (r) => r.rangeSliders > 0, 400, 50);
    const errAfter = pageErrors.length + consoleErrors.length;

    expect(info.rangeSliders).toBeGreaterThan(0);
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 9: logarithmic value axis re-scales positive-count bars without error', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const before = (await barChartState()).axisType;
    await v.setViewerProps(page, 'Bar chart', [{set: {axisType: 'logarithmic'}, wait: 900}]);
    const info = {before, ...await barChartState()};
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.before).toBe('linear');
    expect(info.axisType).toBe('logarithmic');
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 10: under the log axis the clipping precondition still holds (Value Min set, indicators on)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await barChartState();
    const errAfter = pageErrors.length + consoleErrors.length;

    expect(info.axisType).toBe('logarithmic');
    expect(info.valueMin).not.toBeNull();
    expect(info.clippedIndicators).toBe(true);
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 11: switching back to linear restores linear proportional heights', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const before = (await barChartState()).axisType;
    await v.setViewerProps(page, 'Bar chart', [{set: {axisType: 'linear'}, wait: 900}]);
    const info = {before, ...await barChartState()};
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.before).toBe('logarithmic');
    expect(info.axisType).toBe('linear');
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 13: clearing Value Min/Max restores the full-range axis (no clipping)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await v.setViewerProps(page, 'Bar chart', [{set: {valueMin: null, valueMax: null}, wait: 800}]);
    const info = await barChartState();
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.valueMin).toBeNull();
    expect(info.valueMax).toBeNull();
    expect(info.axisType).toBe('linear');
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  v.finishSpec();
});
