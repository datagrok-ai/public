/* ---
realizes: [barchart.cp.sort-and-orient]
--- */

import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
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

function barTopThirds(page: import('@playwright/test').Page) {
  return page.evaluate(() => {
    const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
    const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
    const data = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
    const third = cv.width / 3;
    let leftTopY = cv.height, rightTopY = cv.height, leftCount = 0, rightCount = 0;
    for (let y = 0; y < cv.height; y++)
      for (let x = 0; x < cv.width; x++) {
        const i = (y * cv.width + x) * 4;
        const r = data[i], g = data[i + 1], b = data[i + 2];
        if (r >= 120 && r <= 180 && g >= 190 && g <= 240 && b >= 120 && b <= 180) {
          if (x < third) { leftCount++; if (y < leftTopY) leftTopY = y; }
          else if (x >= cv.width - third) { rightCount++; if (y < rightTopY) rightTopY = y; }
        }
      }
    return {leftTopFrac: leftTopY / cv.height, rightTopFrac: rightTopY / cv.height, leftCount, rightCount};
  });
}

test('Bar Chart — Sorting and Orientation', async ({page}) => {
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

  await v.installEventWaits(page);

  await page.locator('[name="viewer-Bar-chart"]').first().hover();

  await page.locator('[name="viewer-Bar-chart"]')
    .locator('xpath=ancestor::div[contains(@class, "panel-base")]')
    .locator('[name="icon-font-icon-settings"]').first()
    .waitFor({state: 'attached', timeout: 15_000});

  await page.evaluate(() => {
    const bcEl = document.querySelector('[name="viewer-Bar-chart"]') as HTMLElement;
    const panelBase = bcEl.closest('.panel-base') as HTMLElement;
    const gear = panelBase.querySelector('[name="icon-font-icon-settings"]') as HTMLElement;
    gear.click();
  });

  await v.pollValue(() => page.locator('.property-grid').count(), (n) => n > 0, 500, 50);

  await v.setViewerProps(page, 'Bar chart',
    [{set: {splitColumnName: splitCol, valueColumnName: valueCol, valueAggrType: 'sum'}, wait: 1000}]);

  await softStep('Scenario 1 Step 3: vertical + descending by-value sort repaints the chart (canvas delta); tallest bars on the left', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;

    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);

    await v.waitForCanvasQuiet(page, 'Bar chart');
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0); 
    expect(settle.deltaPx).toBeLessThan(500);
    await v.setViewerProps(page, 'Bar chart',
      [{set: {orientation: 'vertical', barSortType: 'by value', barSortOrder: 'desc'}, wait: 1200}]);
    const info = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      return {
        orientation: bc.props.orientation,
        barSortType: bc.props.barSortType,
        barSortOrder: bc.props.barSortOrder,
        hasCanvas: !!cv && cv.getBoundingClientRect().width > 0,
      };
    });
    const {deltaPx} = await v.diffCanvasColors(page, 'Bar chart');

    const bars = await v.pollValue(() => barTopThirds(page),
      (b) => b.leftCount > 0 && b.rightCount > 0 && b.leftTopFrac < b.rightTopFrac, 1200, 100);
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.orientation).toBe('vertical');
    expect(info.barSortType).toBe('by value');
    expect(info.barSortOrder).toBe('desc');
    expect(deltaPx).toBeGreaterThan(1000);    expect(bars.leftCount).toBeGreaterThan(0);
    expect(bars.rightCount).toBeGreaterThan(0);
    expect(bars.leftTopFrac).toBeLessThan(bars.rightTopFrac);
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 4: stacking holds on the negative-sum value column (GROK-19480)', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await page.evaluate(async () => {
      const w = window as any;
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      await w.__settled('viewer:Bar chart.onViewerRendered', () => { bc.props.legendVisibility = 'Always'; }, 2000);
    });

    const legendBefore = await v.readLegend(page, 'Bar chart');
    await v.setViewerProps(page, 'Bar chart', [{set: {stackColumnName: 'Stereo Category'}, wait: 900}]);
    const info = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      return {
        splitColumnName: bc.props.splitColumnName,
        hasCanvas: !!cv && cv.getBoundingClientRect().width > 0,
      };
    });

    const legendAfter = await v.pollValue(() => v.readLegend(page, 'Bar chart'),
      (l) => l.legendRendered && l.itemCount >= 2, 900, 100);
    await v.setViewerProps(page, 'Bar chart', [{set: {stackColumnName: null}, wait: 600}]);
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
    const info = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
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

  await softStep('Scenario 1 Step 7: revert Bar Sort to Ascending re-renders the reorder — the tall-bar side swaps', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;

    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);

    await v.waitForCanvasQuiet(page, 'Bar chart');
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0); 
    expect(settle.deltaPx).toBeLessThan(500);
    const before = await barTopThirds(page);
    const order = (await v.setViewerProps(page, 'Bar chart',
      [{set: {barSortOrder: 'asc'}, wait: 900, read: 'barSortOrder'}]))[0];

    const after = await v.pollValue(() => barTopThirds(page),
      (b) => b.leftCount > 0 && b.rightCount > 0 && b.rightTopFrac < b.leftTopFrac, 900, 100);
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(order).toBe('asc');
    expect(before.leftTopFrac).toBeLessThan(before.rightTopFrac);
    expect(after.leftCount).toBeGreaterThan(0);
    expect(after.rightCount).toBeGreaterThan(0);
    expect(after.rightTopFrac).toBeLessThan(after.leftTopFrac);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 1 Step 8: revert Orientation to Horizontal; ascending, error-free baseline', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await v.setViewerProps(page, 'Bar chart', [{set: {orientation: 'horizontal'}, wait: 900}]);
    const info = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
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
    await v.setViewerProps(page, 'Bar chart', [
      {set: {valueColumnName: countCol, valueAggrType: 'count'}, wait: 900},
      {set: {orientation: 'vertical'}, wait: 1000},
    ]);
    const info = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      return {
        splitColumnName: bc.props.splitColumnName,
        hasCanvas: !!cv && cv.getBoundingClientRect().width > 0,
      };
    });

    const readGreenBars = () => page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      const data = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      let gMinY = cv.height, gMaxY = -1, gCount = 0;
      for (let y = 0; y < cv.height; y++)
        for (let x = 0; x < cv.width; x++) {
          const i = (y * cv.width + x) * 4;
          const r = data[i], g = data[i + 1], b = data[i + 2];
          if (r >= 120 && r <= 180 && g >= 190 && g <= 240 && b >= 120 && b <= 180) {
            gCount++;
            if (y < gMinY) gMinY = y; if (y > gMaxY) gMaxY = y;
          }
        }
      return {gCount, gTopFrac: gMinY / cv.height};
    });

    const bars = await v.pollValue(readGreenBars,
      (b) => b.gCount > 1000 && b.gTopFrac < 0.40, 1000, 100);
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.splitColumnName).toBe(splitCol);
    expect(info.hasCanvas).toBe(true);
    expect(bars.gCount).toBeGreaterThan(1000);
    expect(bars.gTopFrac).toBeLessThan(0.40);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 4: intermediate revert to Horizontal renders error-free', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await v.setViewerProps(page, 'Bar chart', [{set: {orientation: 'horizontal'}, wait: 900}]);
    const info = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      const cv = bc.root.querySelector('canvas') as HTMLCanvasElement;
      return {hasCanvas: !!cv && cv.getBoundingClientRect().width > 0};
    });
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.hasCanvas).toBe(true);
    expect(errAfter).toBe(errBefore);
  });

  await softStep('Scenario 2 Step 5: reset Orientation to Auto returns to default layout, no error', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    await v.setViewerProps(page, 'Bar chart', [{set: {orientation: 'auto'}, wait: 900}]);
    const info = await page.evaluate(() => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
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
