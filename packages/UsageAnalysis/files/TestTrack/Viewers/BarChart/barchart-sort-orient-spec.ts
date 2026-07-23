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

// Positional bar probe (vertical orientation): topmost bar-fill pixel row in
// the left and right thirds of the canvas, as fractions of the canvas height.
// The default bar fill (#96d794 ± tolerance) isolates bars from axis chrome.
// A smaller fraction = a taller bar in that third; a by-value sort puts the
// tallest bars on one side, and flipping the sort order swaps sides — the
// render signal for a reorder, which the color-histogram canvas diff cannot
// see (a pure permutation of equal-width bars keeps the histogram unchanged).
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

  await softStep('Scenario 1 Step 3: vertical + descending by-value sort repaints the chart (canvas delta); tallest bars on the left', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    // Snapshot the horizontal default-sort baseline; the settle precheck proves
    // the setup repaint has drained, so the delta below measures the
    // orientation + sort actuation, not a late setup tail.
    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    await page.waitForTimeout(400);
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0); // -1 = canvas fault
    expect(settle.deltaPx).toBeLessThan(500);
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
    const {deltaPx} = await v.diffCanvasColors(page, 'Bar chart');
    // Positional sort signal: descending by-value under vertical orientation
    // puts the tallest bar on the left — the left third's bar-fill top edge
    // sits above the right third's.
    const bars = await barTopThirds(page);
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
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.legendVisibility = 'Always';
      await new Promise((r) => setTimeout(r, 600));
    });
    // Legend appears as a RESULT of enabling stacking: with no stack column the
    // bar chart shows no legend even under legendVisibility Always, so a build
    // that renders an empty legend host under Always would fail the legendBefore
    // assert and need this step re-based. GROK-19480 broke stacking whenever the
    // aggregated value column has negative sums, which Chemical Space X does.
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

  await softStep('Scenario 1 Step 7: revert Bar Sort to Ascending re-renders the reorder — the tall-bar side swaps', async () => {
    const errBefore = pageErrors.length + consoleErrors.length;
    // A pure reorder permutes equal-width bars, which the color-histogram
    // canvas diff cannot see, so the settle diff here is a fault guard only;
    // the reorder signal is positional: under descending the tallest bars sit
    // in the left third (asserted at Step 3), and flipping to ascending must
    // swap the tall-bar side to the right third.
    expect(await v.snapshotCanvasColors(page, 'Bar chart')).toBe(true);
    await page.waitForTimeout(400);
    const settle = await v.diffCanvasColors(page, 'Bar chart');
    expect(settle.deltaPx).toBeGreaterThanOrEqual(0); // -1 = canvas fault
    expect(settle.deltaPx).toBeLessThan(500);
    const before = await barTopThirds(page);
    const order = await page.evaluate(async () => {
      const bc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Bar chart') as any;
      bc.props.barSortOrder = 'asc';
      await new Promise((r) => setTimeout(r, 900));
      return bc.props.barSortOrder;
    });
    const after = await barTopThirds(page);
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
        hasCanvas: !!cv && cv.getBoundingClientRect().width > 0,
      };
    }, {value: countCol});
    // github-3417 render signal: isolate the bar-fill (#96d794 ± tolerance) from
    // the axis chrome and measure how far the tallest bar's top sits from the
    // canvas top. The raw non-white bbox is useless here — axis labels/ticks
    // span nearly the whole height. The bug manifests as excess whitespace above
    // the tallest bar under vertical + small counts, so the guard asserts the
    // green bars reach the upper region: a regression that shrank the bars would
    // push gTopFrac well above the 0.40 ceiling.
    const bars = await page.evaluate(() => {
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
    const errAfter = pageErrors.length + consoleErrors.length;
    expect(info.splitColumnName).toBe(splitCol);
    expect(info.hasCanvas).toBe(true);
    expect(bars.gCount).toBeGreaterThan(1000);
    expect(bars.gTopFrac).toBeLessThan(0.40);
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
