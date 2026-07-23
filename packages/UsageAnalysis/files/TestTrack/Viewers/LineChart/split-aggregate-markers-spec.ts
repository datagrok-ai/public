/* ---
realizes: [linechart.cp.setup-split-aggregate-markers]
--- */
import {test, expect, type Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep, stepErrors} from '../../spec-login';
import {countCanvasPixels, snapshotCanvasColors, diffCanvasColors} from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

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
  await page.waitForTimeout(400);
}

async function getProps(page: Page, ...names: string[]): Promise<Record<string, any>> {
  return page.evaluate((ns) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    const out: Record<string, any> = {};
    for (const n of ns) out[n] = (lc.props as any)[n];
    return out;
  }, names);
}

// Render assert (GROK-20255 / Step 4 trend line): a blank Line chart paints 0
// non-white px, while a rendered chart with a data line paints far more. The
// 28000 threshold sits far above any axes-only paint and far below the data
// render; -1 (no canvas / getImageData fault) fails the threshold too.
async function chartCanvasNonEmpty(page: Page): Promise<boolean> {
  return (await countCanvasPixels(page, 'Line chart')).total > 28000;
}

// Settle-gated canvas baseline: snapshot the per-color histogram, let any
// pending repaint drain, and prove the frame is stable (settle-precheck diff)
// so the next measured diff isolates the setter's own repaint. diffCanvasColors
// REPLACES the snapshot, so the settled frame becomes the diff baseline.
async function settleCanvasBaseline(page: Page) {
  expect(await snapshotCanvasColors(page, 'Line chart')).toBe(true);
  await page.waitForTimeout(600);
  const settle = await diffCanvasColors(page, 'Line chart');
  expect(settle.deltaPx).toBeGreaterThanOrEqual(0);
  expect(settle.deltaPx).toBeLessThan(200);
}

async function hoverChart(page: Page) {
  const box = await page.evaluate(() => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    const canvases = lc.root.querySelectorAll('canvas');
    let mc: HTMLCanvasElement | null = null, ma = 0;
    for (const c of canvases) {
      const r = (c as HTMLCanvasElement).getBoundingClientRect();
      if (r.width * r.height > ma) { ma = r.width * r.height; mc = c as HTMLCanvasElement; }
    }
    const rect = mc!.getBoundingClientRect();
    return {cx: rect.left + rect.width * 0.5, cy: rect.top + rect.height * 0.5};
  });
  await page.mouse.move(box.cx, box.cy);
  await page.waitForTimeout(300);
}

test('Line Chart — Setup, Split, Aggregate, Markers', async ({page}) => {
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

  await softStep('S1: set X and Y columns (connected trend)', async () => {
    const before = realErrors().length;
    await setProps(page, {xColumnName: 'CAST Idea ID', yColumnNames: ['Chemical Space X']});
    const props = await getProps(page, 'xColumnName', 'yColumnNames');
    expect(props.xColumnName).toBe('CAST Idea ID');
    expect(props.yColumnNames).toEqual(['Chemical Space X']);
    // The 28000-px floor assumes the canvas size this viewport produces — fail
    // loudly if a layout change resizes it rather than silently devaluing the
    // threshold.
    const cvDims = await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Line chart') as any;
      const cv = lc?.root?.querySelector('canvas') as HTMLCanvasElement | null;
      return cv ? {w: cv.width, h: cv.height} : {w: -1, h: -1};
    });
    expect(cvDims.w, 'canvas width changed — recalibrate the 28000-px floor in chartCanvasNonEmpty').toBeGreaterThan(800);
    expect(cvDims.h, 'canvas height changed — recalibrate the 28000-px floor in chartCanvasNonEmpty').toBeGreaterThan(350);
    // The canvas renders a connected trend line, not just a blank frame.
    expect(await chartCanvasNonEmpty(page)).toBe(true);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1: split by Stereo Category', async () => {
    const before = realErrors().length;
    await settleCanvasBaseline(page);
    await setProps(page, {splitColumnName: 'Stereo Category'});
    await page.waitForTimeout(600);
    // Splitting repaints the plot from one line to per-category lines — a
    // large share of the plot ink changes color/position.
    const {deltaPx} = await diffCanvasColors(page, 'Line chart');
    expect(deltaPx).toBeGreaterThan(1000);    expect((await getProps(page, 'splitColumnName')).splitColumnName).toBe('Stereo Category');
    expect(realErrors().length).toBe(before);
  });

  await softStep('S1: aggregation avg + whiskers std err', async () => {
    const before = realErrors().length;
    await settleCanvasBaseline(page);
    await setProps(page, {aggrType: 'avg', whiskersType: 'Avg | ±StError'});
    await page.waitForTimeout(600);
    // Aggregation collapses duplicate X values and the whiskers mode adds
    // error-bar ink on top of the averaged lines.
    const {deltaPx} = await diffCanvasColors(page, 'Line chart');
    expect(deltaPx).toBeGreaterThan(300);    const props = await getProps(page, 'aggrType', 'whiskersType');
    expect(props.aggrType).toBe('avg');
    expect(props.whiskersType).toBe('Avg | ±StError');
    expect(realErrors().length).toBe(before);
  });

  // GROK-19883: marker type together with a size-coding column.
  await softStep('S1: marker type + size-coding column', async () => {
    const before = realErrors().length;
    await settleCanvasBaseline(page);
    await setProps(page, {markerType: 'circle', markersSizeColumnName: 'Chemical Space Y'});
    await page.waitForTimeout(600);
    // Size-coded circle markers add ink at every data point.
    const {deltaPx} = await diffCanvasColors(page, 'Line chart');
    expect(deltaPx).toBeGreaterThan(300);    const props = await getProps(page, 'markerType', 'markersSizeColumnName');
    expect(props.markerType).toBe('circle');
    expect(props.markersSizeColumnName).toBe('Chemical Space Y');
    expect(realErrors().length).toBe(before);
  });

  // GROK-20255: a second Y column combined with a split.
  await softStep('S1: second Y column with split', async () => {
    const before = realErrors().length;
    await setProps(page, {yColumnNames: ['Chemical Space X', 'Chemical Space Y']});
    expect((await getProps(page, 'yColumnNames')).yColumnNames).toHaveLength(2);
    await hoverChart(page);
    // With 2+ split categories a per-category average line renders (GROK-20255):
    // the canvas stays non-empty, not a blank frame.
    expect(await chartCanvasNonEmpty(page)).toBe(true);
    expect(realErrors().length).toBe(before);
  });

  // GROK-17519: many split columns must not hang the page.
  await softStep('S1: add R1/R2/R3 splits, page stays responsive', async () => {
    const before = realErrors().length;
    for (const cols of [
      ['Stereo Category', 'R1'],
      ['Stereo Category', 'R1', 'R2'],
      ['Stereo Category', 'R1', 'R2', 'R3'],
    ]) {
      await page.evaluate((c) => {
        const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
        lc.props.splitColumnNames = c;
      }, cols);
      await page.waitForTimeout(500);
      // Page still responds: a follow-up evaluate resolves within timeout.
      const responsive = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount === 100);
      expect(responsive).toBe(true);
    }
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toHaveLength(4);
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: X axis logarithmic', async () => {
    const before = realErrors().length;
    await setProps(page, {xAxisType: 'logarithmic'});
    expect((await getProps(page, 'xAxisType')).xAxisType).toBe('logarithmic');
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: X axis back to linear', async () => {
    const before = realErrors().length;
    await setProps(page, {xAxisType: 'linear'});
    expect((await getProps(page, 'xAxisType')).xAxisType).toBe('linear');
    expect(realErrors().length).toBe(before);
  });

  await softStep('S2: clear all split columns', async () => {
    const before = realErrors().length;
    await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
      lc.props.splitColumnNames = [];
      lc.props.splitColumnName = '';
    });
    await page.waitForTimeout(400);
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toHaveLength(0);
    expect(realErrors().length).toBe(before);
  });

  if (stepErrors.length > 0)
    throw new Error(`Line Chart split-aggregate-markers failures:\n${stepErrors.join('\n')}`);
});
