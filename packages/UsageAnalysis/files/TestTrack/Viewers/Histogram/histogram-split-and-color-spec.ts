/* ---
realizes: [histogram.cp.split-and-color]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

// Computed opacity of the three Color property-grid rows: 1.0 = color-coding
// active, 0.5 = disabled because a split is active.
async function colorRowOpacities(page: import('@playwright/test').Page): Promise<string[]> {
  return page.evaluate(() => {
    const names = ['prop-color', 'prop-color-aggr-type', 'prop-invert-color-scheme'];
    return names.map((n) => {
      const tr = document.querySelector(`.property-grid tr[name="${n}"]`) as HTMLElement | null;
      return tr ? getComputedStyle(tr).opacity : 'absent';
    });
  });
}

test('Histogram — Split vs Color-coding Transition', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'histogram', 'Histogram');

  // grok.shell.warnings is not exposed to JS here, so uncaught page errors and
  // console errors are the no-error floor for the range-slider step.
  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') pageErrors.push(m.text()); });

  await v.setViewerProps(page, 'Histogram', [{set: {valueColumnName: 'AGE'}, wait: 400}]);
  await page.evaluate(() => {
    const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram');
    grok.shell.o = h;
  });

  // The Context Panel mounts asynchronously after grok.shell.o is assigned, so poll
  // for the Color row and re-issue the assignment on each attempt.
  await expect(async () => {
    const attached = await page.evaluate(() =>
      !!document.querySelector('.property-grid tr[name="prop-color"]'));
    if (!attached) {
      await page.evaluate(() => {
        const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram');
        grok.shell.o = h;
      });
      throw new Error('Context Panel Color property rows not attached yet');
    }
  }).toPass({timeout: 45000, intervals: [500, 1000, 1500, 2000]});
  await page.locator('.property-grid tr[name="prop-color"]').waitFor({state: 'attached', timeout: 15000});

  await softStep('S1: set Color SEX + avg + invert (color config first)', async () => {
    const result = await v.setViewerProps(page, 'Histogram', [
      {set: {colorColumnName: 'SEX'}, wait: 300, read: 'colorColumnName'},
      {set: {colorAggrType: 'avg'}, wait: 200, read: 'colorAggrType'},
      {set: {invertColorScheme: true}, wait: 300, read: 'invertColorScheme'},
    ]);
    expect(result).toEqual(['SEX', 'avg', true]);
  });

  await softStep('S1: Color property rows active (opacity 1.0) with no split', async () => {
    await page.locator('.property-grid tr[name="prop-color"]').waitFor({state: 'attached', timeout: 15000});
    const ops = await colorRowOpacities(page);
    expect(ops).toEqual(['1', '1', '1']);
  });

  const fullCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
  expect(fullCount).toBeGreaterThan(0);

  await softStep('S1: split RACE disables Color rows (opacity 0.5)', async () => {
    await v.setViewerProps(page, 'Histogram', [{set: {splitColumnName: 'RACE'}, wait: 600}]);
    await page.locator('.property-grid tr[name="prop-color"]').waitFor({state: 'attached', timeout: 15000});
    const ops = await colorRowOpacities(page);
    // GROK-19761: an active split disables the Color UI (opacity 0.5).
    for (const o of ops) expect(Number(o)).toBeLessThan(1);
  });

  await softStep('S1: range narrow under split keeps filter valid, no error', async () => {
    const errsBefore = pageErrors.length;
    const {filtered, rowCount} = await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      h.props.valueMin = 20;
      h.props.valueMax = 60;
      await new Promise((r) => setTimeout(r, 700));
      const df = grok.shell.tv.dataFrame;
      return {filtered: df.filter.trueCount, rowCount: df.rowCount};
    });
    // GROK-18399: a range-slider change under a split must not error and must keep
    // df.filter valid.
    expect(filtered).toBeGreaterThanOrEqual(0);
    expect(filtered).toBeLessThanOrEqual(rowCount);
    expect(filtered).toBeLessThanOrEqual(fullCount);
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  await page.evaluate(async () => {
    const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
    h.props.valueMin = null;
    h.props.valueMax = null;
    await new Promise((r) => setTimeout(r, 400));
  });

  await softStep('S1: clear split re-activates Color rows (opacity 1.0 round-trip)', async () => {
    await v.setViewerProps(page, 'Histogram', [{set: {splitColumnName: null}, wait: 600}]);
    await page.locator('.property-grid tr[name="prop-color"]').waitFor({state: 'attached', timeout: 15000});
    const ops = await colorRowOpacities(page);
    expect(ops).toEqual(['1', '1', '1']); // GROK-19761 regression guard (round-trip)
  });

  v.finishSpec();
});
