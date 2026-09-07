/* ---
realizes: [linechart.cp.setup-split-aggregate-markers]
--- */
import {expect, type Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:AppData/Chem/tests/spgi-100.csv';

async function setProps(page: Page, props: Record<string, any>) {
  await v.setViewerProps(page, 'Line chart', [{set: props, wait: 400}]);
}

async function getProps(page: Page, ...names: string[]): Promise<Record<string, any>> {
  return page.evaluate((ns) => {
    const lc = Array.from(grok.shell.tv.viewers).find((v: any) => v.type === 'Line chart') as any;
    const out: Record<string, any> = {};
    for (const n of ns) out[n] = (lc.props as any)[n];
    return out;
  }, names);
}

async function chartCanvasNonEmpty(page: Page): Promise<boolean> {
  return (await v.countCanvasPixels(page, 'Line chart')).total > 28000;
}

// The previous step's repaint may still be in flight (a split repaints in a burst), so the
// baseline is the first interval in which the canvas holds still, not a snapshot taken on arrival.
async function settleCanvasBaseline(page: Page) {
  expect(await v.snapshotCanvasColors(page, 'Line chart')).toBe(true);
  const settle = await v.pollValue(() => v.diffCanvasColors(page, 'Line chart'), (d) => d.deltaPx < 200, 1500, 150);
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
  await v.waitForViewerRendered(page, 'Line chart', 300);
}

test('Line Chart — Setup, Split, Aggregate, Markers', async ({page}) => {
  test.setTimeout(300_000);

  const pageErrors: string[] = [];
  const consoleErrors: string[] = [];
  page.on('pageerror', (e) => pageErrors.push(String(e)));
  page.on('console', (m) => { if (m.type() === 'error' && !isLocalBootNoise(m.text())) consoleErrors.push(m.text()); });
  const errorCount = () => pageErrors.length + consoleErrors.length;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'line-chart', 'Line-chart', 15_000, 'Line chart');

  await softStep('S1: set X and Y columns (connected trend)', async () => {
    const before = errorCount();
    await setProps(page, {xColumnName: 'CAST Idea ID', yColumnNames: ['Chemical Space X']});
    const props = await getProps(page, 'xColumnName', 'yColumnNames');
    expect(props.xColumnName).toBe('CAST Idea ID');
    expect(props.yColumnNames).toEqual(['Chemical Space X']);

    const cvDims = await page.evaluate(() => {
      const lc = Array.from(grok.shell.tv.viewers).find((x: any) => x.type === 'Line chart') as any;
      const cv = lc?.root?.querySelector('canvas') as HTMLCanvasElement | null;
      return cv ? {w: cv.width, h: cv.height} : {w: -1, h: -1};
    });
    expect(cvDims.w, 'canvas width changed — recalibrate the 28000-px floor in chartCanvasNonEmpty').toBeGreaterThan(800);
    expect(cvDims.h, 'canvas height changed — recalibrate the 28000-px floor in chartCanvasNonEmpty').toBeGreaterThan(350);

    expect(await chartCanvasNonEmpty(page)).toBe(true);
    expect(errorCount()).toBe(before);
  });

  await softStep('S1: split by Stereo Category', async () => {
    const before = errorCount();
    await settleCanvasBaseline(page);
    await setProps(page, {splitColumnName: 'Stereo Category'});

    const deltaPx = await v.waitForCanvasChange(page, 'Line chart', {minDelta: 1000, timeoutMs: 600});
    expect(deltaPx).toBeGreaterThan(1000);
    expect((await getProps(page, 'splitColumnName')).splitColumnName).toBe('Stereo Category');
    expect(errorCount()).toBe(before);
  });

  await softStep('S1: aggregation avg + whiskers std err', async () => {
    const before = errorCount();
    await settleCanvasBaseline(page);
    await setProps(page, {aggrType: 'avg', whiskersType: 'Avg | ±StError'});

    const deltaPx = await v.waitForCanvasChange(page, 'Line chart', {minDelta: 300, timeoutMs: 600});
    expect(deltaPx).toBeGreaterThan(300);
    const props = await getProps(page, 'aggrType', 'whiskersType');
    expect(props.aggrType).toBe('avg');
    expect(props.whiskersType).toBe('Avg | ±StError');
    expect(errorCount()).toBe(before);
  });

  await softStep('S1: marker type + size-coding column', async () => {
    const before = errorCount();
    await settleCanvasBaseline(page);
    await setProps(page, {markerType: 'circle', markersSizeColumnName: 'Chemical Space Y'});

    const deltaPx = await v.waitForCanvasChange(page, 'Line chart', {minDelta: 300, timeoutMs: 600});
    expect(deltaPx).toBeGreaterThan(300);
    const props = await getProps(page, 'markerType', 'markersSizeColumnName');
    expect(props.markerType).toBe('circle');
    expect(props.markersSizeColumnName).toBe('Chemical Space Y');
    expect(errorCount()).toBe(before);
  });

  await softStep('S1: second Y column with split', async () => {
    const before = errorCount();
    await setProps(page, {yColumnNames: ['Chemical Space X', 'Chemical Space Y']});
    expect((await getProps(page, 'yColumnNames')).yColumnNames).toHaveLength(2);
    await hoverChart(page);

    expect(await chartCanvasNonEmpty(page)).toBe(true);
    expect(errorCount()).toBe(before);
  });

  await softStep('S1: add R1/R2/R3 splits, page stays responsive', async () => {
    const before = errorCount();
    for (const cols of [
      ['Stereo Category', 'R1'],
      ['Stereo Category', 'R1', 'R2'],
      ['Stereo Category', 'R1', 'R2', 'R3'],
    ]) {
      await v.setViewerProps(page, 'Line chart', [{set: {splitColumnNames: cols}, wait: 500}]);

      const responsive = await page.evaluate(() => grok.shell.tv.dataFrame.rowCount === 100);
      expect(responsive).toBe(true);
    }
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toHaveLength(4);
    expect(errorCount()).toBe(before);
  });

  await softStep('S2: X axis logarithmic', async () => {
    const before = errorCount();
    await setProps(page, {xAxisType: 'logarithmic'});
    expect((await getProps(page, 'xAxisType')).xAxisType).toBe('logarithmic');
    expect(errorCount()).toBe(before);
  });

  await softStep('S2: X axis back to linear', async () => {
    const before = errorCount();
    await setProps(page, {xAxisType: 'linear'});
    expect((await getProps(page, 'xAxisType')).xAxisType).toBe('linear');
    expect(errorCount()).toBe(before);
  });

  await softStep('S2: clear all split columns', async () => {
    const before = errorCount();
    await v.setViewerProps(page, 'Line chart',
      [{set: {splitColumnNames: [], splitColumnName: ''}, wait: 400}]);
    expect((await getProps(page, 'splitColumnNames')).splitColumnNames).toHaveLength(0);
    expect(errorCount()).toBe(before);
  });

  await v.closeAllAndWait(page);
  v.finishSpec();
});
