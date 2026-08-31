/* ---
realizes: [linechart.cp.filter-follow-and-empty]
--- */
import {expect, type Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, stepErrors, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';
const fullRowCount = 100;

const pageErrors: string[] = [];
const consoleErrors: string[] = [];

function realErrors(): string[] {
  return [...pageErrors, ...consoleErrors];
}

async function setProps(page: Page, props: Record<string, any>) {
  await v.setViewerProps(page, 'Line chart', [{set: props, wait: 500}]);
}

async function getProps(page: Page, ...names: string[]): Promise<Record<string, any>> {
  return page.evaluate((ns) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    const out: Record<string, any> = {};
    for (const n of ns) out[n] = (lc.props as any)[n];
    return out;
  }, names);
}

async function filterSeries(
  page: Page, selected: string[] | null, settled: (count: number) => boolean,
): Promise<number> {
  await page.evaluate((sel) => {
    const tv = grok.shell.tv;
    const df = tv.dataFrame;
    const fg = tv.getFiltersGroup();
    const values = sel === null ? df.columns.byName('Series').categories : sel;
    fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Series', selected: values});
    df.rows.requestFilter();
  }, selected);
  return v.pollValue(() => trueCount(page), settled, 800, 100);
}

async function filterRange(
  page: Page, column: string, frac: [number, number] | null, settled: (count: number) => boolean,
): Promise<number> {
  await page.evaluate((args) => {
    const {col, f} = args;
    const tv = grok.shell.tv;
    const df = tv.dataFrame;
    const fg = tv.getFiltersGroup();
    const c = df.columns.byName(col);
    const lo = f === null ? c.min : c.min + (c.max - c.min) * f[0];
    const hi = f === null ? c.max : c.min + (c.max - c.min) * f[1];
    fg.updateOrAdd({type: 'histogram', column: col, min: lo, max: hi});
    df.rows.requestFilter();
  }, {col: column, f: frac});
  return v.pollValue(() => trueCount(page), settled, 800, 100);
}

async function trueCount(page: Page): Promise<number> {
  return page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
}

test('Line Chart — Filter Follow and Empty-Chart Resilience', async ({page}) => {
  test.setTimeout(300_000);
  stepErrors.length = 0;
  pageErrors.length = 0;
  consoleErrors.length = 0;

  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error' && !isLocalBootNoise(m.text())) consoleErrors.push(m.text()); });

  await openDatagrok(page);

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  await page.locator('[name="icon-line-chart"]').click();
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 15000});

  await setProps(page, {xColumnName: 'Chemical Space X', yColumnNames: ['Chemical Space Y'], axesFollowFilter: true});
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await v.pollValue(() => page.locator('[name="viewer-Filters"] .d4-filter').count(),
    (n) => n > 0, 1500, 100);

  const baseline = await trueCount(page);
  expect(baseline).toBe(fullRowCount);

  await softStep('S1: switch X axis to logarithmic', async () => {
    await setProps(page, {xAxisType: 'logarithmic'});
    expect((await getProps(page, 'xAxisType')).xAxisType).toBe('logarithmic');
  });

  await softStep('S1 Step 4: filter to zero rows, filter.trueCount === 0', async () => {
    const cnt = await filterSeries(page, [], (n) => n === 0);
    expect(cnt).toBe(0);
  });

  await softStep('S1 Step 5-6: hover empty log-axis chart raises no error (github-2574)', async () => {
    const before = realErrors().length;
    await page.locator('[name="viewer-Line-chart"]').hover();

    await page.waitForTimeout(800);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1 Step 7-8: remove filter restores baseline row count', async () => {
    const cnt = await filterSeries(page, null, (n) => n === baseline);
    expect(cnt).toBe(baseline);
  });

  await softStep('S1 Step 9: switch X axis back to linear', async () => {
    await setProps(page, {xAxisType: 'linear'});
    expect((await getProps(page, 'xAxisType')).xAxisType).toBe('linear');
  });

  await softStep('S2 Step 1-2: set X to a date column with Year-quarter time-split, no error', async () => {
    const before = realErrors().length;

    await setProps(page, {xColumnName: 'Competition assay Date', xMap: 'Year quarter'});
    const cfg = await getProps(page, 'xColumnName', 'xMap');
    expect(cfg.xColumnName).toBe('Competition assay Date');
    expect(cfg.xMap).toBe('Year quarter');
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2 Step 5: filter while time-split X keeps rows > 0 (GROK-18375)', async () => {
    const before = realErrors().length;
    const cnt = await filterSeries(page, ['Aminopiperidines', 'Pyrrolidines'],
      (n) => n > 0 && n < baseline);
    expect(cnt).toBeGreaterThan(0);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2 Step 7: remove filter, restore numeric X for next scenario', async () => {
    await filterSeries(page, null, (n) => n === baseline);
    await setProps(page, {xColumnName: 'Chemical Space X', xMap: ''});
    const cfg = await getProps(page, 'xColumnName', 'xMap');
    expect(cfg.xColumnName).toBe('Chemical Space X');
    expect(cfg.xMap).toBe(''); 
    expect(await trueCount(page)).toBe(baseline);
  });

  await softStep('S3 Step 1: confirm numeric X, linear scale', async () => {
    const cfg = await getProps(page, 'xColumnName', 'xAxisType');
    expect(cfg.xColumnName).toBe('Chemical Space X');
    expect(cfg.xAxisType).toBe('linear');
  });

  await softStep('S3 Step 4-5: narrow X range, filter.trueCount drops below baseline (GROK-20185)', async () => {
    const before = realErrors().length;
    const narrowed = await filterRange(page, 'Chemical Space X', [0.25, 0.75],
      (n) => n > 0 && n < baseline);
    expect(narrowed).toBeLessThan(baseline);
    expect(narrowed).toBeGreaterThan(0);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S3 Step 7-8: restore full range, filter.trueCount returns to baseline', async () => {
    const restored = await filterRange(page, 'Chemical Space X', null, (n) => n === baseline);
    expect(restored).toBe(baseline);
  });

  if (stepErrors.length > 0)
    throw new Error(`Line Chart filter-follow-and-empty failures:\n${stepErrors.join('\n')}`);
});
