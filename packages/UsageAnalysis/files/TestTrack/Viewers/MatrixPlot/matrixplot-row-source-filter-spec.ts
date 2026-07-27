/* ---
realizes: [matrixplot.cp.row-source-filter, matrixplot-row-source-filter-narrow-cells, viewers.matrix-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

/** Summed non-white ink over every inner cell canvas — a deterministic
 * whole-matrix render signature at a fixed viewer size. */
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

/** Re-read the whole-matrix ink until two consecutive readings agree. */
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

test('Matrix Plot — Row Source and Filtering', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'matrix-plot', 'Matrix-plot');

  await softStep('Scenario 1 — Row Source switched to Selected repaints over the selection', async () => {
    await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      const df = grok.shell.tv.dataFrame;
      df.selection.setAll(false);
      mp.props.rowSource = 'Filtered';
      mp.props.filter = '';
      await new Promise((r) => setTimeout(r, 700));
    });
    const baseline = await settledMatrixInk(page);
    expect(baseline).toBeGreaterThan(0);

    const selCount = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      df.selection.init((i: number) => i < 50);
      await new Promise((r) => setTimeout(r, 300));
      return df.selection.trueCount;
    });
    expect(selCount).toBe(50);

    await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.rowSource = 'Selected';
      await new Promise((r) => setTimeout(r, 800));
    });
    const selInk = await settledMatrixInk(page);
    console.log(`MatrixPlot rowSource ink: baseline=${baseline} selected=${selInk} delta=${Math.abs(selInk - baseline)}`);
    // Plotting only 50 of 5850 rows is a dramatic repaint; the floor keeps a
    // wide margin below the observed delta.
    expect(Math.abs(selInk - baseline)).toBeGreaterThan(2000);

    await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      const df = grok.shell.tv.dataFrame;
      mp.props.rowSource = 'Filtered';
      df.selection.setAll(false);
      await new Promise((r) => setTimeout(r, 800));
    });
    const restoreInk = await settledMatrixInk(page);
    // Rendering is deterministic for identical state at a fixed viewer size —
    // the round-trip restores the exact baseline signature.
    expect(restoreInk).toBe(baseline);
  });

  await softStep('Scenario 2 — Filter Panel SEX = M narrows the shared filter', async () => {
    const result = await page.evaluate(async () => {
      const df = grok.shell.tv.dataFrame;
      const fg = grok.shell.tv.getFiltersGroup();
      await new Promise((r) => setTimeout(r, 700));
      const full = df.filter.trueCount;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: ['M']});
      await new Promise((r) => setTimeout(r, 800));
      const filtered = df.filter.trueCount;
      const cats = df.col('SEX').categories;
      fg.updateOrAdd({type: DG.FILTER_TYPE.CATEGORICAL, column: 'SEX', selected: cats});
      await new Promise((r) => setTimeout(r, 800));
      const restored = df.filter.trueCount;
      return {full, filtered, restored};
    });
    // The panel filter drives the shared table filter the matrix listens to.
    expect(result.filtered).toBeLessThan(result.full);
    expect(result.restored).toBe(result.full);
  });

  await softStep('Scenario 3 — Viewer Filter property repaints without touching the shared filter', async () => {
    await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.filter = '';
      await new Promise((r) => setTimeout(r, 500));
    });
    const baseline = await settledMatrixInk(page);
    const fullFilter = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);

    await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.filter = '${AGE} > 30';
      await new Promise((r) => setTimeout(r, 900));
    });
    const filterInk = await settledMatrixInk(page);
    const filterCount = await page.evaluate(() => grok.shell.tv.dataFrame.filter.trueCount);
    console.log(`MatrixPlot viewer-filter ink: baseline=${baseline} filtered=${filterInk} delta=${Math.abs(filterInk - baseline)}`);
    // The viewer-level formula narrows only this viewer's rows — a modest
    // repaint; the floor keeps the assert non-vacuous with margin.
    expect(Math.abs(filterInk - baseline)).toBeGreaterThan(100);
    // The viewer Filter prop does NOT drive df.filter — the shared count is unchanged.
    expect(filterCount).toBe(fullFilter);

    await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.filter = '';
      await new Promise((r) => setTimeout(r, 900));
    });
    const restoreInk = await settledMatrixInk(page);
    expect(restoreInk).toBe(baseline);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
