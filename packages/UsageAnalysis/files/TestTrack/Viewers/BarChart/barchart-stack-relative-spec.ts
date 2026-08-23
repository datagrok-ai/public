/* ---
realizes: [barchart.cp.stack-relative-negatives]
--- */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';
const splitCol = 'Primary Series Name';
const stackCol = 'Scaffold Names';
const valueCol = 'Chemical Space X';

test('Bar Chart — Stacking, Relative Values, and Negative Aggregates', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 4000});
  await v.addViewerByIcon(page, 'bar-chart', 'Bar-chart');

  await page.locator('[name="viewer-Bar-chart"]').first().hover();
  await v.waitForViewerRendered(page, 'Bar chart', 300);

  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear.click();
  });
  await v.pollValue(() => page.evaluate(() => !!document.querySelector('.property-grid')),
    (up) => up, 500, 100);

  await v.setViewerProps(page, 'Bar chart',
    [{set: {splitColumnName: splitCol, valueColumnName: valueCol, valueAggrType: 'sum'}, wait: 900}]);

  const legendItemCount = () => page.evaluate(() => {
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    return bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length as number;
  });
  const legendGone = () => page.evaluate(() => {
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    const el = bc.root.querySelector('[name="legend"]') as HTMLElement | null;
    const laidOut = !!el && getComputedStyle(el).display !== 'none' &&
      el.getBoundingClientRect().height > 0 && el.offsetParent !== null;
    const items = bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
    return !laidOut && items === 0;
  });
  const setStackColumn = (stack: string | null) => page.evaluate((s: string | null) => {
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    bc.props.stackColumnName = s;
  }, stack);

  await softStep('Scenario 1 Step 4: Relative Values + Stack normalizes bars to equal width (canvas delta), chart not blank (github-2659)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;

    await setStackColumn(stackCol);

    await v.waitForCanvasQuiet(page, 'Bar chart', {timeoutMs: 900, optional: true});
    const pre = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      return {stackDefaultRel: bc.props.relativeValues};
    });

    await v.waitForCanvasQuiet(page, 'Bar chart', {timeoutMs: 900, optional: true});

    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0); 
    expect(settle.deltaPx).toBeLessThan(500);
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.relativeValues = true;
    });

    const deltaPx = await v.waitForCanvasChange(page, 'Bar chart', {minDelta: 5001, timeoutMs: 1300});
    const info = await page.evaluate(({split, stack}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv ? cv.getBoundingClientRect() : {width: 0, height: 0};
      const legendHost = bc.root.querySelector('[name="legend"]');
      const legendItems = bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      const df = grok.shell.tv.dataFrame;
      return {
        stack: bc.props.stackColumnName,
        split: bc.props.splitColumnName,
        rel: bc.props.relativeValues,
        hasCanvas: !!cv,
        canvasW: rect.width,
        canvasH: rect.height,
        legendRendered: !!legendHost,
        legendItems,
        outerCats: df.col(split).categories.length,
        stackCats: df.col(stack).categories.length,
      };
    }, {split: splitCol, stack: stackCol});
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.stack).toBe(stackCol);
    expect(info.split).toBe(splitCol);
    expect(pre.stackDefaultRel).toBe(false);
    expect(info.rel).toBe(true);
    expect(deltaPx).toBeGreaterThan(5000);
    expect(info.hasCanvas).toBe(true);
    expect(info.canvasW).toBeGreaterThan(0);
    expect(info.canvasH).toBeGreaterThan(0);
    expect(info.outerCats).toBeGreaterThanOrEqual(2);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 5: at least two stacked segments visible (legend has >= 2 stack categories)', async () => {
    const info = await page.evaluate(({stack}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const legendItems = bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      return {
        legendRendered: !!bc.root.querySelector('[name="legend"]'),
        legendItems,
        stackCats: grok.shell.tv.dataFrame.col(stack).categories.length,
      };
    }, {stack: stackCol});
    expect(info.legendRendered).toBe(true);
    expect(info.legendItems).toBeGreaterThanOrEqual(2);
    expect(info.stackCats).toBeGreaterThanOrEqual(2);
  });

  await softStep('Scenario 1 Step 7: stacked bars render error-free with negative sum totals (GROK-19480)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(({value}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const df = grok.shell.tv.dataFrame;
      const col = df.col(value);
      let hasNeg = false;
      let minV = Infinity;
      for (let i = 0; i < df.rowCount; i++) {
        const val = col.get(i);
        if (val === null || val === undefined || Number.isNaN(val)) continue;
        if (val < minV) minV = val;
        if (val < 0) hasNeg = true;
      }
      return {
        aggr: bc.props.valueAggrType,
        value: bc.props.valueColumnName,
        hasNegativeValues: hasNeg,
        minValue: minV,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    }, {value: valueCol});
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.aggr).toBe('sum');
    expect(info.value).toBe(valueCol);
    expect(info.hasNegativeValues).toBe(true);
    expect(info.minValue).toBeLessThan(0);
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 9: disabling Relative Values reverts bars to absolute widths (canvas delta)', async () => {

    await v.waitForCanvasQuiet(page, 'Bar chart', {timeoutMs: 900, optional: true});

    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0); 
    expect(settle.deltaPx).toBeLessThan(500);
    await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.relativeValues = false;
    });
    const deltaPx = await v.waitForCanvasChange(page, 'Bar chart', {minDelta: 5001, timeoutMs: 1100});
    const info = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv ? cv.getBoundingClientRect() : {width: 0, height: 0};
      return {
        rel: bc.props.relativeValues,
        hasCanvas: !!cv,
        canvasW: rect.width,
        canvasH: rect.height,
        stack: bc.props.stackColumnName,
      };
    });
    expect(info.rel).toBe(false);
    expect(info.stack).toBe(stackCol);
    expect(info.hasCanvas).toBe(true);
    expect(info.canvasW).toBeGreaterThan(0);
    expect(deltaPx).toBeGreaterThan(5000);
  });

  await softStep('Scenario 1 Step 11: removing the Stack column collapses to single-segment bars, no legend', async () => {
    await setStackColumn(null);
    await v.pollValue(legendGone, (gone) => gone, 800, 100);
    const info = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const el = bc.root.querySelector('[name="legend"]') as HTMLElement | null;
      const legendLaidOut = !!el && getComputedStyle(el).display !== 'none' &&
        el.getBoundingClientRect().height > 0 && el.offsetParent !== null;
      return {
        stack: bc.props.stackColumnName,
        legendLaidOut,
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    });
    expect(info.stack).toBeNull();
    expect(info.legendLaidOut).toBe(false);
    expect(info.legendItems).toBe(0);
    expect(info.hasCanvas).toBe(true);
  });

  await softStep('Scenario 2 Step 3: Relative Values without a Stack column has no stacking effect, no error', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;

    const entry = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      return {
        stackAtEntry: bc.props.stackColumnName,
        legendAtBaseline: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
      };
    });
    await v.setViewerProps(page, 'Bar chart', [{set: {relativeValues: true}, wait: 800}]);
    const info = await page.evaluate((base: {stackAtEntry: string | null; legendAtBaseline: number}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv ? cv.getBoundingClientRect() : {width: 0, height: 0};
      const el = bc.root.querySelector('[name="legend"]') as HTMLElement | null;
      const legendLaidOut = !!el && getComputedStyle(el).display !== 'none' &&
        el.getBoundingClientRect().height > 0 && el.offsetParent !== null;
      return {
        stackAtEntry: base.stackAtEntry,
        legendAtBaseline: base.legendAtBaseline,
        rel: bc.props.relativeValues,
        stackAfter: bc.props.stackColumnName,
        legendLaidOut,
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
        hasCanvas: !!cv,
        canvasW: rect.width,
        canvasH: rect.height,
      };
    }, entry);
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.stackAtEntry).toBeNull();
    expect(info.legendAtBaseline).toBe(0);
    expect(info.rel).toBe(true);
    expect(info.stackAfter).toBeNull();
    expect(info.legendLaidOut).toBe(false);
    expect(info.legendItems).toBe(0);
    expect(info.hasCanvas).toBe(true);
    expect(info.canvasW).toBeGreaterThan(0);
    expect(info.canvasH).toBeGreaterThan(0);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 5-6: setting a Stack column activates stacking — bars normalize to equal width (canvas delta), legend renders', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;

    await v.waitForCanvasQuiet(page, 'Bar chart', {timeoutMs: 900, optional: true});

    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0);
    expect(settle.deltaPx).toBeLessThan(500);
    await setStackColumn(stackCol);
    const deltaPx = await v.waitForCanvasChange(page, 'Bar chart', {minDelta: 5001, timeoutMs: 1300});
    const info = await page.evaluate(({stack}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv ? cv.getBoundingClientRect() : {width: 0, height: 0};
      return {
        stack: bc.props.stackColumnName,
        rel: bc.props.relativeValues,
        hasCanvas: !!cv,
        canvasW: rect.width,
        legendRendered: !!bc.root.querySelector('[name="legend"]'),
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
        stackCats: grok.shell.tv.dataFrame.col(stack).categories.length,
      };
    }, {stack: stackCol});
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.stack).toBe(stackCol);
    expect(info.rel).toBe(true);
    expect(deltaPx).toBeGreaterThan(5000);
    expect(info.hasCanvas).toBe(true);
    expect(info.canvasW).toBeGreaterThan(0);
    expect(info.legendRendered).toBe(true);
    expect(info.legendItems).toBeGreaterThanOrEqual(2);
    expect(info.stackCats).toBeGreaterThanOrEqual(2);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 7-8: remove Stack reverts to inert; re-add re-activates (repeatable)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await setStackColumn(null);
    await v.pollValue(legendGone, (gone) => gone, 900, 100);
    const afterRemove = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const el = bc.root.querySelector('[name="legend"]') as HTMLElement | null;
      const legendLaidOut = !!el && getComputedStyle(el).display !== 'none' &&
        el.getBoundingClientRect().height > 0 && el.offsetParent !== null;
      return {
        stack: bc.props.stackColumnName,
        legendLaidOut,
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    });
    await setStackColumn(stackCol);
    await v.pollValue(legendItemCount, (n) => n >= 2, 900, 100);
    const info = await page.evaluate((removed: typeof afterRemove) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const el = bc.root.querySelector('[name="legend"]') as HTMLElement | null;
      const legendLaidOut = !!el && getComputedStyle(el).display !== 'none' &&
        el.getBoundingClientRect().height > 0 && el.offsetParent !== null;
      const afterReadd = {
        stack: bc.props.stackColumnName,
        legendLaidOut,
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
      };
      return {afterRemove: removed, afterReadd, rel: bc.props.relativeValues};
    }, afterRemove);
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.afterRemove.stack).toBeNull();
    expect(info.afterRemove.legendLaidOut).toBe(false);
    expect(info.afterRemove.legendItems).toBe(0);
    expect(info.afterRemove.hasCanvas).toBe(true);
    expect(info.afterReadd.stack).toBe(stackCol);
    expect(info.afterReadd.legendLaidOut).toBe(true);
    expect(info.afterReadd.legendItems).toBeGreaterThanOrEqual(2);
    expect(info.rel).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 9-10: disabling Relative Values (no Stack) restores the baseline', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await setStackColumn(null);
    await v.pollValue(legendGone, (gone) => gone, 700, 100);
    await v.setViewerProps(page, 'Bar chart', [{set: {relativeValues: false}, wait: 700}]);
    const info = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv ? cv.getBoundingClientRect() : {width: 0, height: 0};
      const el = bc.root.querySelector('[name="legend"]') as HTMLElement | null;
      const legendLaidOut = !!el && getComputedStyle(el).display !== 'none' &&
        el.getBoundingClientRect().height > 0 && el.offsetParent !== null;
      return {
        rel: bc.props.relativeValues,
        stack: bc.props.stackColumnName,
        hasCanvas: !!cv,
        canvasW: rect.width,
        canvasH: rect.height,
        legendLaidOut,
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.rel).toBe(false);
    expect(info.stack).toBeNull();
    expect(info.hasCanvas).toBe(true);
    expect(info.canvasW).toBeGreaterThan(0);
    expect(info.canvasH).toBeGreaterThan(0);
    expect(info.legendLaidOut).toBe(false);
    expect(info.legendItems).toBe(0);
    expect(errAfter).toBe(errBefore);
  });

  v.finishSpec();
});
