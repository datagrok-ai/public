/* ---
realizes: [histogram.cp.range-filter-bounds]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

// The Min/Max range text inputs (.d4-filter-input-min / .d4-filter-input-max) appear
// inside the viewer root when showRangeInputs=true; filter.trueCount is the readable
// downstream state, the canvas is not.

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

test('Histogram — Range filtering and bound validation', async ({page}) => {
  test.setTimeout(300_000);

  await loginToDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'histogram', 'Histogram');

  // grok.shell.warnings is undefined on this build, so uncaught page errors and
  // console errors are the no-throw floor for the canvas and validation steps.
  const pageErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error') pageErrors.push(m.text()); });

  const minSel = '[name="viewer-Histogram"] .d4-filter-input-min';
  const maxSel = '[name="viewer-Histogram"] .d4-filter-input-max';

  const rowCount = await page.evaluate(async () => {
    const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
    h.props.valueColumnName = 'AGE';
    h.props.showRangeInputs = true;
    h.props.filteringEnabled = true;
    await new Promise((r) => setTimeout(r, 800));
    return grok.shell.tv.dataFrame.rowCount;
  });
  expect(rowCount).toBeGreaterThan(0);
  await page.locator(minSel).first().waitFor({timeout: 15000});

  await softStep('S1: Min=30/Max=60 filters to AGE sub-range (< full row count)', async () => {
    await page.locator(minSel).first().fill('30');
    await page.locator(minSel).first().press('Enter');
    await page.locator(maxSel).first().fill('60');
    await page.locator(maxSel).first().press('Enter');
    await page.waitForTimeout(700);
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
    expect(Math.abs(filtered - inRange)).toBeLessThanOrEqual(1); // bin-edge tolerance
  });

  await softStep('S1: filter.trueCount unchanged after Split Stack (GROK-18948)', async () => {
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    const after = await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      h.props.splitColumnName = 'SEX';
      h.props.splitStack = true;
      await new Promise((r) => setTimeout(r, 700));
      return grok.shell.tv.dataFrame.filter.trueCount;
    });
    expect(after).toBe(before); // GROK-18948: range filter survives stacking
    // Disable-split teardown gets its own no-error slice — otherwise a fault here
    // would be mis-attributed to the next step's errsBefore baseline.
    const teardownErrs = pageErrors.length;
    await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      h.props.splitStack = false;
      h.props.splitColumnName = '';
      await new Promise((r) => setTimeout(r, 400));
    });
    expect(pageErrors.slice(teardownErrs)).toEqual([]);
  });

  await softStep('S1: Normalise to Filter raises no error (github-2329)', async () => {
    const errsBefore = pageErrors.length;
    await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      h.props.normalizeToFilter = true;
      await new Promise((r) => setTimeout(r, 500));
    });
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  await softStep('S1: Min=40 narrows filter further, no error (github-2329)', async () => {
    const errsBefore = pageErrors.length;
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await page.locator(minSel).first().fill('40');
    await page.locator(minSel).first().press('Enter');
    await page.waitForTimeout(700);
    const after = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(after).toBeLessThan(before); // narrower range excludes more rows
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  await softStep('S2: establish valid Min=40/Max=60 sub-range', async () => {
    await page.locator(minSel).first().fill('40');
    await page.locator(minSel).first().press('Enter');
    await page.locator(maxSel).first().fill('60');
    await page.locator(maxSel).first().press('Enter');
    await page.waitForTimeout(500);
    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBeGreaterThan(0);
    expect(filtered).toBeLessThan(rowCount);
  });

  // NOTE: on this build the range inputs are NOT rejected or clamped — an
  // out-of-range bound is applied verbatim (widening the effective range to the
  // data extent), and an inverted Max<Min collapses the filter rather than being
  // refused. The tickets (GROK-19581/GROK-19760) guaranteed only "no crash". The
  // real, non-tautological signal per step is the direction of the trueCount move
  // plus the no-error floor.

  await softStep('S2: Max below Min collapses filter without crash (GROK-19581, GROK-19760)', async () => {
    const errsBefore = pageErrors.length;
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await page.locator(maxSel).first().fill('20');
    await page.locator(maxSel).first().press('Enter');
    await page.waitForTimeout(600);
    const after = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(after).toBeLessThan(before); // inverted Max<Min collapses the sub-range
    expect(pageErrors.slice(errsBefore)).toEqual([]); // no crash on Min>=Max violation
  });

  await softStep('S2: Min below column min widens range without crash (GROK-19581)', async () => {
    const errsBefore = pageErrors.length;
    await page.locator(maxSel).first().fill('60');
    await page.locator(maxSel).first().press('Enter');
    await page.waitForTimeout(300);
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await page.locator(minSel).first().fill('-999');
    await page.locator(minSel).first().press('Enter');
    await page.waitForTimeout(600);
    const {after, min, max} = await page.evaluate(({mn, mx}) => ({
      after: grok.shell.tv.dataFrame.filter.trueCount,
      min: parseFloat((document.querySelector(mn) as HTMLInputElement)?.value),
      max: parseFloat((document.querySelector(mx) as HTMLInputElement)?.value),
    }), {mn: minSel, mx: maxSel});
    expect(after).toBeGreaterThan(before); // lower bound below the data extent lets in the rest
    expect(min).toBeLessThanOrEqual(max); // effective range not inverted
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  await softStep('S2: Max above column max widens range without crash (GROK-19760)', async () => {
    const errsBefore = pageErrors.length;
    await page.locator(minSel).first().fill('40');
    await page.locator(minSel).first().press('Enter');
    await page.locator(maxSel).first().fill('60');
    await page.locator(maxSel).first().press('Enter');
    await page.waitForTimeout(300);
    const before = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    await page.locator(maxSel).first().fill('999');
    await page.locator(maxSel).first().press('Enter');
    await page.waitForTimeout(600);
    const {after, min, max} = await page.evaluate(({mn, mx}) => ({
      after: grok.shell.tv.dataFrame.filter.trueCount,
      min: parseFloat((document.querySelector(mn) as HTMLInputElement)?.value),
      max: parseFloat((document.querySelector(mx) as HTMLInputElement)?.value),
    }), {mn: minSel, mx: maxSel});
    expect(after).toBeGreaterThan(before); // upper bound above the data extent lets in the rest
    expect(min).toBeLessThanOrEqual(max); // effective range not inverted
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  await softStep('S3: restore full range returns filter to full row count (round-trip)', async () => {
    const {colMin, colMax} = await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      h.props.normalizeToFilter = false;
      const age = grok.shell.tv.dataFrame.getCol('AGE');
      await new Promise((r) => setTimeout(r, 300));
      return {colMin: age.min, colMax: age.max};
    });
    await page.locator(minSel).first().fill(String(colMin));
    await page.locator(minSel).first().press('Enter');
    await page.locator(maxSel).first().fill(String(colMax));
    await page.locator(maxSel).first().press('Enter');
    await page.waitForTimeout(700);
    const filtered = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    expect(filtered).toBe(rowCount); // full range restores all rows
  });

  await softStep('S3: disable Filter leaves viewer no-error (teardown floor)', async () => {
    const errsBefore = pageErrors.length;
    await page.evaluate(async () => {
      const h = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Histogram') as any;
      h.props.filteringEnabled = false;
      await new Promise((r) => setTimeout(r, 500));
    });
    expect(pageErrors.slice(errsBefore)).toEqual([]);
  });

  v.finishSpec();
});
