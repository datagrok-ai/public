/* ---
realizes: [barchart.cp.sort-and-orient]
--- */

import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';
const splitCol = 'Primary Series Name';
const valueCol = 'Chemical Space X';
const countCol = 'CAST Idea ID';

function categorySums(page: import('@playwright/test').Page, split: string, value: string) {
  return page.evaluate(({split, value}) => {
    const df = grok.shell.tv.dataFrame;
    const s = df.col(split), val = df.col(value);
    const sums: Record<string, number> = {};
    for (let i = 0; i < df.rowCount; i++) {
      if (!df.filter.get(i)) continue;
      const k = String(s.get(i));
      const val_i = val.get(i);
      if (val_i === null || val_i === undefined || Number.isNaN(val_i)) continue;
      sums[k] = (sums[k] || 0) + val_i;
    }
    return sums;
  }, {split, value});
}

test('Bar Chart — Sorting and Orientation', async ({page}) => {
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
    await new Promise((r) => setTimeout(r, 1000));
  }, {split: splitCol, value: valueCol});

  await softStep('Scenario 1 Step 3: vertical + descending by-value config took, error-free render', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.orientation = 'vertical';
      bc.props.barSortType = 'by value';
      bc.props.barSortOrder = 'desc';
      await new Promise((r) => setTimeout(r, 1200));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      return {
        orientation: bc.props.orientation,
        barSortType: bc.props.barSortType,
        barSortOrder: bc.props.barSortOrder,
        hasCanvas: !!cv && cv.getBoundingClientRect().width > 0,
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.orientation).toBe('vertical');
    expect(info.barSortType).toBe('by value');
    expect(info.barSortOrder).toBe('desc');
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 4: stacking holds on the negative-sum value column (GROK-19480)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.legendVisibility = 'Always';
      await new Promise((r) => setTimeout(r, 600));
    });
    // Legend appears as a RESULT of enabling stacking (with no stack column the
    // bar chart shows no legend even under legendVisibility Always), so the
    // before/after delta below is caused by the stack action itself. GROK-19480
    // broke stacking whenever the aggregated value column has negative sums,
    // which Chemical Space X does.
    const legendBefore = await v.readLegend(page, 'Bar chart');
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.stackColumnName = 'Stereo Category';
      await new Promise((r) => setTimeout(r, 900));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      return {
        splitColumnName: bc.props.splitColumnName,
        hasCanvas: !!cv && cv.getBoundingClientRect().width > 0,
      };
    });
    const legendAfter = await v.readLegend(page, 'Bar chart');
    await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.stackColumnName = null;
      await new Promise((r) => setTimeout(r, 600));
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(legendBefore.legendRendered).toBe(false);
    expect(legendAfter.legendRendered).toBe(true);
    expect(legendAfter.itemCount).toBeGreaterThanOrEqual(2);
    expect(info.splitColumnName).toBe(splitCol);
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 5: negative-sum category renders error-free below baseline (GROK-19480)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const sums = await categorySums(page, splitCol, valueCol);
    const negatives = Object.entries(sums).filter(([, val]) => val < 0);
    expect(negatives.length).toBeGreaterThanOrEqual(1);
    expect(Math.min(...Object.values(sums))).toBeLessThan(0);
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      await new Promise((r) => setTimeout(r, 600));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      return {
        aggr: bc.props.valueAggrType,
        value: bc.props.valueColumnName,
        hasCanvas: !!cv && cv.getBoundingClientRect().width > 0,
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.aggr).toBe('sum');
    expect(info.value).toBe(valueCol);
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 7: revert Bar Sort to Ascending; config took, error-free render', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const order = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.barSortOrder = 'asc';
      await new Promise((r) => setTimeout(r, 900));
      return bc.props.barSortOrder;
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(order).toBe('asc');
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 8: revert Orientation to Horizontal; ascending, error-free baseline', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.orientation = 'horizontal';
      await new Promise((r) => setTimeout(r, 900));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      return {
        order: bc.props.barSortOrder,
        hasCanvas: !!cv && cv.getBoundingClientRect().width > 0,
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.order).toBe('asc');
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 3: vertical with small-magnitude counts renders error-free (github-3417)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async ({value}) => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.valueColumnName = value;
      bc.props.valueAggrType = 'count';
      await new Promise((r) => setTimeout(r, 900));
      bc.props.orientation = 'vertical';
      await new Promise((r) => setTimeout(r, 1000));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      return {
        splitColumnName: bc.props.splitColumnName,
        barSortType: bc.props.barSortType,
        hasCanvas: !!cv && cv.getBoundingClientRect().width > 0,
      };
    }, {value: countCol});
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.splitColumnName).toBe(splitCol);
    expect(info.barSortType).toBe('by value');
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 4: intermediate revert to Horizontal renders error-free', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.orientation = 'horizontal';
      await new Promise((r) => setTimeout(r, 900));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      return {hasCanvas: !!cv && cv.getBoundingClientRect().width > 0};
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 5: reset Orientation to Auto returns to default layout, no error', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    const info = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.orientation = 'auto';
      await new Promise((r) => setTimeout(r, 900));
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      return {
        valueAggrType: bc.props.valueAggrType,
        valueColumnName: bc.props.valueColumnName,
        splitColumnName: bc.props.splitColumnName,
        hasCanvas: !!cv && cv.getBoundingClientRect().width > 0,
      };
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.valueAggrType).toBe('count');
    expect(info.valueColumnName).toBe(countCol);
    expect(info.splitColumnName).toBe(splitCol);
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  v.finishSpec();
});
