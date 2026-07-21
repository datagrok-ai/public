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
  await page.waitForTimeout(300);

  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear.click();
  });
  await page.waitForTimeout(500);

  await page.evaluate(async ({split, value}) => {
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    bc.props.splitColumnName = split;
    bc.props.valueColumnName = value;
    bc.props.valueAggrType = 'sum';
    await new Promise((r) => setTimeout(r, 900));
  }, {split: splitCol, value: valueCol});

  await softStep('Scenario 1 Step 4: Relative Values + Stack normalizes bars to equal width (canvas delta), chart not blank (github-2659)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    // Set the Stack column and let the absolute-width stacked layout settle.
    const pre = await page.evaluate(async ({stack}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.stackColumnName = stack;
      await new Promise((r) => setTimeout(r, 900));
      return {stackDefaultRel: bc.props.relativeValues};
    }, {stack: stackCol});
    // Snapshot the absolute-width stacked layout, then enable Relative Values:
    // each outer bar normalizes to equal width — the width normalization is a
    // large canvas color delta (not just a prop echo).
    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    await page.waitForTimeout(400);
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeLessThan(500);
    const info = await page.evaluate(async ({split, stack}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.relativeValues = true;
      await new Promise((r) => setTimeout(r, 1000));
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
    await page.waitForTimeout(300);
    const {deltaPx} = await v.diffCanvasColors(page, 'Bar chart');
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
    // Snapshot the normalized (equal-width) layout, disable Relative Values, and
    // measure the canvas color delta as bars revert to absolute widths.
    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    await page.waitForTimeout(400);
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeLessThan(500);
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.relativeValues = false;
      await new Promise((r) => setTimeout(r, 800));
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
    await page.waitForTimeout(300);
    const {deltaPx} = await v.diffCanvasColors(page, 'Bar chart');
    expect(info.rel).toBe(false);
    expect(info.stack).toBe(stackCol);
    expect(info.hasCanvas).toBe(true);
    expect(info.canvasW).toBeGreaterThan(0);
    expect(deltaPx).toBeGreaterThan(5000);
  });

  await softStep('Scenario 1 Step 11: removing the Stack column collapses to single-segment bars, no legend', async () => {
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.stackColumnName = null;
      await new Promise((r) => setTimeout(r, 800));
      return {
        stack: bc.props.stackColumnName,
        legendRendered: !!bc.root.querySelector('[name="legend"]'),
        legendItems: bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length,
        hasCanvas: !!bc.root.querySelector('canvas'),
      };
    });
    expect(info.stack).toBeNull();
    expect(info.legendRendered).toBe(false);
    expect(info.legendItems).toBe(0);
    expect(info.hasCanvas).toBe(true);
  });

  await softStep('Scenario 2 Step 3: Relative Values without a Stack column has no stacking effect, no error', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const stackAtEntry = bc.props.stackColumnName;
      bc.props.relativeValues = true;
      await new Promise((r) => setTimeout(r, 800));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const rect = cv ? cv.getBoundingClientRect() : {width: 0, height: 0};
      const legendItems = bc.root.querySelectorAll('[name="legend"] .d4-legend-item').length;
      bc.props.relativeValues = false;
      await new Promise((r) => setTimeout(r, 500));
      return {
        stackAtEntry,
        legendRendered: !!bc.root.querySelector('[name="legend"]'),
        legendItems,
        hasCanvas: !!cv,
        canvasW: rect.width,
        canvasH: rect.height,
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.stackAtEntry).toBeNull();
    expect(info.legendRendered).toBe(false);
    expect(info.legendItems).toBe(0);
    expect(info.hasCanvas).toBe(true);
    expect(info.canvasW).toBeGreaterThan(0);
    expect(errAfter).toBe(errBefore);
  });

  v.finishSpec();
});
