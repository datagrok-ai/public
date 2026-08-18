/* ---
realizes: []
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Stack trace [A-Za-z]+/.test(text);

const cellInk = (page: Page, idx: number) => page.evaluate((i: number) => {
  const cells = document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer');
  const c = cells[i] as HTMLCanvasElement | undefined;
  if (!c) return -1;
  const ctx = c.getContext('2d');
  if (!ctx) return -1;
  let data: Uint8ClampedArray;
  try { data = ctx.getImageData(0, 0, c.width, c.height).data; } catch (_) { return -2; }
  let n = 0;
  for (let k = 0; k < data.length; k += 16)
    if (data[k + 3] !== 0 && !(data[k] >= 250 && data[k + 1] >= 250 && data[k + 2] >= 250)) n++;
  return n;
}, idx);

async function settledCellInk(page: Page, idx: number): Promise<number> {
  let prev = await cellInk(page, idx);
  let cur = prev;
  for (let i = 0; i < 10; i++) {
    await page.waitForTimeout(300); 
    cur = await cellInk(page, idx);
    if (cur >= 0 && Math.abs(cur - prev) < 40) break;
    prev = cur;
  }
  return cur;
}

async function settledCellInkAfterChange(
  page: Page, idx: number, from: number, capMs: number,
): Promise<number> {
  await v.pollValue(() => cellInk(page, idx),
    (ink) => ink >= 0 && Math.abs(ink - from) >= 40, capMs, 100);
  return settledCellInk(page, idx);
}

const cellCount = (page: Page) => page.evaluate(() =>
  document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer').length);

test('Matrix plot tests', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text()); });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'matrix-plot', 'Matrix-plot');

  await softStep('Default State — auto-picked X/Y sets are numerical/datetime, at most 10 each', async () => {
    const sets = await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      const df = grok.shell.tv.dataFrame;
      const check = (names: string[]) => names.map((n) => {
        const c = df.col(n);
        return c.isNumerical || c.type === 'datetime';
      });
      return {x: mp.props.xColumnNames, y: mp.props.yColumnNames,
        xOk: check(mp.props.xColumnNames), yOk: check(mp.props.yColumnNames)};
    });
    expect(sets.x.length).toBeGreaterThan(0);
    expect(sets.x.length).toBeLessThanOrEqual(10);
    expect(sets.y.length).toBeLessThanOrEqual(10);
    expect(sets.xOk.every((b: boolean) => b)).toBe(true);
    expect(sets.yOk.every((b: boolean) => b)).toBe(true);
  });

  await softStep('Back Color — property round-trip only (paint not measured, open GROK-20439)', async () => {
    const errBefore = errCount();
    const [red, white] = await v.setViewerProps(page, 'Matrix plot', [
      {set: {backColor: 0xFFFF0000}, wait: 400, read: 'backColor'},
      {set: {backColor: 0xFFFFFFFF}, wait: 300, read: 'backColor'},
    ]);
    const result = {red: red >>> 0, white: white >>> 0};

    expect(result.red).toBe(0xFFFF0000);
    expect(result.white).toBe(0xFFFFFFFF);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Title and Description — rendered text tracks the properties', async () => {
    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.showTitle = true;
      mp.props.title = 'My Matrix';
      mp.props.description = 'Test description';
    });

    const set = await v.pollValue(() => page.evaluate(() => {

      const titleShown = [...document.querySelectorAll('.panel-titlebar-text')]
        .some((e) => e.textContent!.trim() === 'My Matrix');
      const root = document.querySelector('[name="viewer-Matrix-plot"]')!;
      const descShown = [...root.querySelectorAll('*')]
        .some((e) => e.children.length === 0 && e.textContent!.trim() === 'Test description');
      const titleTexts = [...document.querySelectorAll('.panel-titlebar-text')].map((e) => e.textContent!.trim());
      return {titleShown, descShown, titleTexts};
    }), (s) => s.titleShown && s.descShown, 900, 100);
    console.log(`MatrixPlot title/desc: shown=${set.titleShown} desc=${set.descShown} titles=${JSON.stringify(set.titleTexts)}`);
    expect(set.titleShown).toBe(true);
    expect(set.descShown).toBe(true);

    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;

      mp.props.showTitle = false;
      mp.props.title = '';
      mp.props.description = '';
    });
    const cleared = await v.pollValue(() => page.evaluate(() => {
      const titleGone = ![...document.querySelectorAll('.panel-titlebar-text')]
        .some((e) => e.textContent!.trim() === 'My Matrix');
      return {titleGone};
    }), (c) => c.titleGone, 800, 100);
    expect(cleared.titleGone).toBe(true);
  });

  await softStep('Inner Viewer Look — Cell Plot Type repaints; inner size is a prop round-trip + no-error floor', async () => {
    const errBefore = errCount();

    const densInk = await settledCellInk(page, 1);
    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.cellPlotType = 'Scatter plot';
    });
    const scatInk = await settledCellInkAfterChange(page, 1, densInk, 900);
    expect(Math.abs(scatInk - densInk)).toBeGreaterThan(500);

    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      const ivl = mp.props.innerViewerLook;
      ivl.markerSize = 10;
      mp.props.innerViewerLook = ivl;
    });

    await v.waitForViewerRendered(page, 'Matrix plot', 900);
    const marker = await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      return mp.props.innerViewerLook.markerSize;
    });
    expect(marker).toBe(10);
    expect(await cellCount(page)).toBe(16);

    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.cellPlotType = 'Density plot';
    });
    const backInk = await settledCellInkAfterChange(page, 1, scatInk, 900);
    expect(Math.abs(backInk - scatInk)).toBeGreaterThan(500);

    expect(Math.abs(backInk - densInk)).toBeLessThan(40);

    await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      const ivl = mp.props.innerViewerLook;
      if ('binCount' in ivl) ivl.binCount = 20; else ivl.bins = 20;
      mp.props.innerViewerLook = ivl;
    });
    await v.waitForViewerRendered(page, 'Matrix plot', 900);
    const bin = await page.evaluate(() => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      return mp.props.innerViewerLook.binCount ?? mp.props.innerViewerLook.bins;
    });
    expect(bin).toBe(20);
    expect(await cellCount(page)).toBe(16);
    expect(errCount()).toBe(errBefore);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
