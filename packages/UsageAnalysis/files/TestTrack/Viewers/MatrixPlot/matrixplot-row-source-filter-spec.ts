/* ---
realizes: [matrixplot.cp.row-source-filter, matrixplot.int.row-source-filter-narrow-cells, viewers.matrix-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const matrixInk = (page: Page) => page.evaluate(() => {
  const cells = document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer');
  let total = 0;
  for (const c of Array.from(cells) as HTMLCanvasElement[]) {
    const ctx = c.getContext('2d');
    if (!ctx) continue;
    let data: Uint8ClampedArray;
    try { data = ctx.getImageData(0, 0, c.width, c.height).data; } catch (_) { continue; }
    for (let k = 0; k < data.length; k += 16)
      if (data[k + 3] !== 0 && !(data[k] >= 250 && data[k + 1] >= 250 && data[k + 2] >= 250)) total++;
  }
  return total;
});

async function settledMatrixInk(page: Page): Promise<number> {
  let prev = await matrixInk(page);
  let cur = prev;
  for (let i = 0; i < 12; i++) {
    await page.waitForTimeout(300); 
    cur = await matrixInk(page);
    if (Math.abs(cur - prev) < 80) break;
    prev = cur;
  }
  return cur;
}

async function settledMatrixInkAfterChange(page: Page, from: number, capMs: number): Promise<number> {
  await v.pollValue(() => matrixInk(page), (ink) => Math.abs(ink - from) >= 80, capMs, 100);
  return settledMatrixInk(page);
}

test('Matrix Plot — Row Source and Filtering', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'matrix-plot', 'Matrix-plot');

  const filterTrueCount = () => page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount as number);

  await softStep('Scenario 1 — Row Source switched to Selected repaints over the selection', async () => {
    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      mp.props.rowSource = 'Filtered';
      mp.props.filter = '';
    });

    await v.waitForViewerRendered(page, 'Matrix plot', 700);
    const baseline = await settledMatrixInk(page);
    expect(baseline).toBeGreaterThan(0);

    await page.evaluate(() => {
      grok.shell.tv.dataFrame.selection.init((i: number) => i < 50);
    });
    const selCount = await v.pollValue(
      () => page.evaluate(() => grok.shell.tv.dataFrame.selection.trueCount as number),
      (n) => n === 50, 300, 100);
    expect(selCount).toBe(50);

    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.rowSource = 'Selected';
    });
    const selInk = await settledMatrixInkAfterChange(page, baseline, 800);
    console.log(`MatrixPlot rowSource ink: baseline=${baseline} selected=${selInk} delta=${Math.abs(selInk - baseline)}`);

    expect(Math.abs(selInk - baseline)).toBeGreaterThan(2000);

    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      const df = grok.shell.tv.dataFrame;
      mp.props.rowSource = 'Filtered';
      df.selection.setAll(false);
    });
    const restoreInk = await settledMatrixInkAfterChange(page, selInk, 800);

    expect(restoreInk).toBe(baseline);
  });

  await softStep('Scenario 2 — Filter Panel SEX = M narrows the shared filter', async () => {
    await page.evaluate(() => { grok.shell.tv.getFiltersGroup(); });
    await v.pollValue(() => page.evaluate(() => !!document.querySelector('[name="viewer-Filters"]')),
      (up) => up, 700, 100);
    const full = await filterTrueCount();
    await page.evaluate(() => {
      const fg = grok.shell.tv.getFiltersGroup();
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M']});
    });
    const filtered = await v.pollValue(filterTrueCount, (n) => n < full, 800, 100);
    await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const fg = grok.shell.tv.getFiltersGroup();
      const cats = df.col('SEX').categories;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: cats});
    });
    const restored = await v.pollValue(filterTrueCount, (n) => n === full, 800, 100);
    const result = {full, filtered, restored};

    expect(result.filtered).toBeLessThan(result.full);
    expect(result.restored).toBe(result.full);
  });

  await softStep('Scenario 3 — Viewer Filter property repaints without touching the shared filter', async () => {
    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.filter = '';
    });

    await v.waitForViewerRendered(page, 'Matrix plot', 500);
    const baseline = await settledMatrixInk(page);
    const fullFilter = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);

    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.filter = '${AGE} > 30';
    });
    const filterInk = await settledMatrixInkAfterChange(page, baseline, 900);
    const filterCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    console.log(`MatrixPlot viewer-filter ink: baseline=${baseline} filtered=${filterInk} delta=${Math.abs(filterInk - baseline)}`);

    expect(Math.abs(filterInk - baseline)).toBeGreaterThan(100);

    expect(filterCount).toBe(fullFilter);

    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.filter = '';
    });
    const restoreInk = await settledMatrixInkAfterChange(page, filterInk, 900);
    expect(restoreInk).toBe(baseline);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
