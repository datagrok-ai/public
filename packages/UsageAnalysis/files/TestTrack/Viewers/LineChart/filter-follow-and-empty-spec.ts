/* ---
realizes: [linechart.cp.filter-follow-and-empty]
--- */
import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/SPGI.csv';
const fullRowCount = 3624;

const pageErrors: string[] = [];
const consoleErrors: string[] = [];

function realErrors(): string[] {
  return [...pageErrors, ...consoleErrors];
}

async function setProps(page: Page, props: Record<string, any>) {
  await page.evaluate((p) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    for (const [k, val] of Object.entries(p)) (lc.props as any)[k] = val;
  }, props);
  await page.waitForTimeout(500);
}

async function getProps(page: Page, ...names: string[]): Promise<Record<string, any>> {
  return page.evaluate((ns) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    const out: Record<string, any> = {};
    for (const n of ns) out[n] = (lc.props as any)[n];
    return out;
  }, names);
}

// The filter card is a canvas, so fg.updateOrAdd + requestFilter is the driving
// path. selected:[] yields 0 rows, while a histogram range clamps to >=1, so the
// empty-chart state has to come from the categorical filter.
async function filterSeries(page: Page, selected: string[] | null): Promise<number> {
  return page.evaluate((sel) => {
    const tv = grok.shell.tv;
    const df = tv.dataFrame;
    const fg = tv.getFiltersGroup();
    const values = sel === null ? df.columns.byName('Series').categories : sel;
    fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'Series', selected: values});
    df.rows.requestFilter();
    return new Promise<number>((resolve) => setTimeout(() => resolve(df.filter.trueCount), 800));
  }, selected);
}

async function filterRange(page: Page, column: string, frac: [number, number] | null): Promise<number> {
  return page.evaluate((args) => {
    const {col, f} = args;
    const tv = grok.shell.tv;
    const df = tv.dataFrame;
    const fg = tv.getFiltersGroup();
    const c = df.columns.byName(col);
    const lo = f === null ? c.min : c.min + (c.max - c.min) * f[0];
    const hi = f === null ? c.max : c.min + (c.max - c.min) * f[1];
    fg.updateOrAdd({type: 'histogram', column: col, min: lo, max: hi});
    df.rows.requestFilter();
    return new Promise<number>((resolve) => setTimeout(() => resolve(df.filter.trueCount), 800));
  }, {col: column, f: frac});
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
  page.on('console', (m) => { if (m.type() === 'error') consoleErrors.push(m.text()); });

  await loginToDatagrok(page);

  await page.evaluate(async (path) => {
    document.body.classList.add('selenium');
    grok.shell.settings.showFiltersIconsConstantly = true;
    grok.shell.windows.simpleMode = true;
    grok.shell.closeAll();
    const df = await grok.dapi.files.readCsv(path);
    grok.shell.addTableView(df);
    await new Promise((resolve) => {
      const sub = df.onSemanticTypeDetected.subscribe(() => { sub.unsubscribe(); resolve(undefined); });
      setTimeout(resolve, 3000);
    });
    for (let i = 0; i < 50; i++) {
      if (document.querySelector('[name="viewer-Grid"] canvas')) break;
      await new Promise((r) => setTimeout(r, 200));
    }
    await new Promise((r) => setTimeout(r, 5000));
  }, datasetPath);
  await page.locator('.d4-grid[name="viewer-Grid"]').waitFor({timeout: 60000});

  await page.locator('[name="icon-line-chart"]').click();
  await page.locator('[name="viewer-Line-chart"]').waitFor({timeout: 15000});

  await setProps(page, {xColumnName: 'Chemical Space X', yColumnNames: ['Chemical Space Y'], axesFollowFilter: true});
  await page.evaluate(() => grok.shell.tv.getFiltersGroup());
  await page.waitForTimeout(1500);

  const baseline = await trueCount(page);
  expect(baseline).toBe(fullRowCount);

  await softStep('S1: switch X axis to logarithmic', async () => {
    await setProps(page, {xAxisType: 'logarithmic'});
    expect((await getProps(page, 'xAxisType')).xAxisType).toBe('logarithmic');
  });

  await softStep('S1 Step 4: filter to zero rows, filter.trueCount === 0', async () => {
    const cnt = await filterSeries(page, []);
    expect(cnt).toBe(0);
  });

  await softStep('S1 Step 5: hover empty log-axis chart raises no error (github-2574)', async () => {
    const before = realErrors().length;
    await page.locator('[name="viewer-Line-chart"]').hover();
    await page.waitForTimeout(800);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1 Step 7-8: remove filter restores baseline row count', async () => {
    const cnt = await filterSeries(page, null);
    expect(cnt).toBe(baseline);
  });

  await softStep('S1 Step 9: switch X axis back to linear', async () => {
    await setProps(page, {xAxisType: 'linear'});
    expect((await getProps(page, 'xAxisType')).xAxisType).toBe('linear');
  });

  // GROK-18375: with a time-split X axis, applying a filter removed every row. The
  // Chem substructure filter cannot be driven headless, so a categorical filter on
  // the structure-family Series column stands in for it.
  await softStep('S2 Step 1-2: set X to a date column with Year-quarter time-split, no error', async () => {
    const before = realErrors().length;
    // xMap is the time-split the GROK-18375 repro requires — without it the
    // bug's configuration (time-split x filter) is not actually exercised.
    await setProps(page, {xColumnName: 'Competition assay Date', xMap: 'Year quarter'});
    const cfg = await getProps(page, 'xColumnName', 'xMap');
    expect(cfg.xColumnName).toBe('Competition assay Date');
    expect(cfg.xMap).toBe('Year quarter');
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2 Step 5: filter while time-split X keeps rows > 0 (GROK-18375)', async () => {
    const before = realErrors().length;
    const cnt = await filterSeries(page, ['Aminopiperidines', 'Pyrrolidines']);
    expect(cnt).toBeGreaterThan(0);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2 Step 7: remove filter, restore numeric X for next scenario', async () => {
    await filterSeries(page, null);
    await setProps(page, {xColumnName: 'Chemical Space X', xMap: ''});
    const cfg = await getProps(page, 'xColumnName', 'xMap');
    expect(cfg.xColumnName).toBe('Chemical Space X');
    expect(cfg.xMap).toBe(''); // time-split must not leak into Scenario 3
    expect(await trueCount(page)).toBe(baseline);
  });

  // The range-slider drag is gesture-uncontrollable-headless — MCP recon
  // 2026-07-21: the Filter Panel numeric filter is an embedded Histogram with a
  // canvas-drawn range strip (no DOM handle elements), canvas-strip drags at
  // two positions leave df.filter untouched, and the min/max d4-filter-inputs
  // exist but stay invisible (null boundingBox) even after hover. The range is
  // therefore driven through fg.updateOrAdd, which still verifies the
  // GROK-20185 observable: the live count update and its reversibility.
  await softStep('S3 Step 1: confirm numeric X, linear scale', async () => {
    const cfg = await getProps(page, 'xColumnName', 'xAxisType');
    expect(cfg.xColumnName).toBe('Chemical Space X');
    expect(cfg.xAxisType).toBe('linear');
  });

  await softStep('S3 Step 4-5: narrow X range, filter.trueCount drops below baseline (GROK-20185)', async () => {
    const before = realErrors().length;
    const narrowed = await filterRange(page, 'Chemical Space X', [0.25, 0.75]);
    expect(narrowed).toBeLessThan(baseline);
    expect(narrowed).toBeGreaterThan(0);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S3 Step 7-8: restore full range, filter.trueCount returns to baseline', async () => {
    const restored = await filterRange(page, 'Chemical Space X', null);
    expect(restored).toBe(baseline);
  });

  if (stepErrors.length > 0)
    throw new Error(`Line Chart filter-follow-and-empty failures:\n${stepErrors.join('\n')}`);
});
