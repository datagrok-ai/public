/* ---
realizes: [histogram.cp.range-filter-bounds]
--- */
import {localTest as test, expect} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Histogram — Range filtering and bound validation', async ({page}) => {
  test.setTimeout(300_000);

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'histogram', 'Histogram');

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => {
    if (m.type() === 'error' && !isLocalBootNoise(m.text())) pageErrors.push(m.text());
  });

  const minSel = '[name="viewer-Histogram"] .d4-filter-input-min';
  const maxSel = '[name="viewer-Histogram"] .d4-filter-input-max';

  const waitForFilter = (ok: (n: number) => boolean, capMs: number) =>
    v.pollValue(() => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount), ok, capMs, 50);

  await v.setViewerProps(page, 'Histogram',
    [{set: {valueColumnName: 'AGE', showRangeInputs: true, filteringEnabled: true}, wait: 800}]);
  const rowCount = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount);
  expect(rowCount).toBeGreaterThan(0);
  await page.locator(minSel).first().waitFor({timeout: 15000});

  await softStep('S1: Min=30/Max=60 filters to AGE sub-range (< full row count)', async () => {
    await page.locator(minSel).first().fill('30');
    await page.locator(minSel).first().press('Enter');
    await page.locator(maxSel).first().fill('60');
    await page.locator(maxSel).first().press('Enter');
    await waitForFilter((n) => n < rowCount, 700);
    const {filtered, inRange, total} = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const age = df.getCol('AGE');
      let inRange = 0;
      for (let i = 0; i < df.rowCount; i++) {
        const val = age.get(i);
        if (val !== null && val >= 30 && val <= 60) inRange++;
      }
      return {filtered: df.filter.trueCount, inRange, total: df.rowCount};
    });
    expect(filtered).toBeLessThan(total);
    expect(Math.abs(filtered - inRange)).toBeLessThanOrEqual(1); 
  });

  await softStep('S1: filter.trueCount unchanged after Split Stack (GROK-18948)', async () => {
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await v.setViewerProps(page, 'Histogram', [{set: {splitColumnName: 'SEX', splitStack: true}, wait: 700}]);
    const after = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(after).toBe(before); 

    const teardownErrs = pageErrors.length;
    await v.setViewerProps(page, 'Histogram', [{set: {splitStack: false, splitColumnName: ''}, wait: 400}]);
    expect(pageErrors.slice(teardownErrs)).toEqual([]);
  });

  await softStep('S1: Normalise to Filter raises no error (github-2329)', async () => {
    const errsBefore = pageErrors.length;
    await v.setViewerProps(page, 'Histogram', [{set: {normalizeToFilter: true}, wait: 500}]);
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  await softStep('S1: Min=40 narrows filter further, no error (github-2329)', async () => {
    const errsBefore = pageErrors.length;
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await page.locator(minSel).first().fill('40');
    await page.locator(minSel).first().press('Enter');
    const after = await waitForFilter((n) => n < before, 700);
    expect(after).toBeLessThan(before); 
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  await softStep('S2: establish valid Min=40/Max=60 sub-range', async () => {
    await page.locator(minSel).first().fill('40');
    await page.locator(minSel).first().press('Enter');
    await page.locator(maxSel).first().fill('60');
    await page.locator(maxSel).first().press('Enter');

    await v.waitForViewerRendered(page, 'Histogram', 500);
    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBeGreaterThan(0);
    expect(filtered).toBeLessThan(rowCount);
  });

  await softStep('S2: Max below Min collapses filter without crash (GROK-19581, GROK-19760)', async () => {
    const errsBefore = pageErrors.length;
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await page.locator(maxSel).first().fill('20');
    await page.locator(maxSel).first().press('Enter');
    const after = await waitForFilter((n) => n < before, 600);
    expect(after).toBeLessThan(before); 
    expect(pageErrors.slice(errsBefore)).toEqual([]); 
  });

  await softStep('S2: Min below column min widens range without crash (GROK-19581)', async () => {
    const errsBefore = pageErrors.length;

    const collapsed = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await page.locator(maxSel).first().fill('60');
    await page.locator(maxSel).first().press('Enter');
    const before = await waitForFilter((n) => n > collapsed, 300);
    await page.locator(minSel).first().fill('-999');
    await page.locator(minSel).first().press('Enter');
    await waitForFilter((n) => n > before, 600);
    const {after, min, max} = await page.evaluate(({mn, mx}) => ({
      after: grok.shell.tv.dataFrame.filter.trueCount,
      min: parseFloat((document.querySelector(mn) as HTMLInputElement)?.value),
      max: parseFloat((document.querySelector(mx) as HTMLInputElement)?.value),
    }), {mn: minSel, mx: maxSel});
    expect(after).toBeGreaterThan(before); 
    expect(min).toBeLessThanOrEqual(max); 
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  await softStep('S2: Max above column max widens range without crash (GROK-19760)', async () => {
    const errsBefore = pageErrors.length;

    const widened = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await page.locator(minSel).first().fill('40');
    await page.locator(minSel).first().press('Enter');
    await page.locator(maxSel).first().fill('60');
    await page.locator(maxSel).first().press('Enter');
    const before = await waitForFilter((n) => n < widened, 300);
    await page.locator(maxSel).first().fill('999');
    await page.locator(maxSel).first().press('Enter');
    await waitForFilter((n) => n > before, 600);
    const {after, min, max} = await page.evaluate(({mn, mx}) => ({
      after: grok.shell.tv.dataFrame.filter.trueCount,
      min: parseFloat((document.querySelector(mn) as HTMLInputElement)?.value),
      max: parseFloat((document.querySelector(mx) as HTMLInputElement)?.value),
    }), {mn: minSel, mx: maxSel});
    expect(after).toBeGreaterThan(before); 
    expect(min).toBeLessThanOrEqual(max); 
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  await softStep('S3: restore full range returns filter to full row count (round-trip)', async () => {
    await v.setViewerProps(page, 'Histogram', [{set: {normalizeToFilter: false}, wait: 300}]);
    const {colMin, colMax} = await page.evaluate(() => {
      const age = grok.shell.tv.dataFrame.getCol('AGE');
      return {colMin: age.min, colMax: age.max};
    });
    await page.locator(minSel).first().fill(String(colMin));
    await page.locator(minSel).first().press('Enter');
    await page.locator(maxSel).first().fill(String(colMax));
    await page.locator(maxSel).first().press('Enter');
    const filtered = await waitForFilter((n) => n === rowCount, 700);
    expect(filtered).toBe(rowCount); 
  });

  await softStep('S3: disable Filter leaves viewer no-error (teardown floor)', async () => {
    const errsBefore = pageErrors.length;
    await v.setViewerProps(page, 'Histogram', [{set: {filteringEnabled: false}, wait: 500}]);
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  v.finishSpec();
});
