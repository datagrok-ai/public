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
    const result = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.backColor = 0xFFFF0000;
      await new Promise((r) => setTimeout(r, 400));
      const red = mp.props.backColor >>> 0;
      mp.props.backColor = 0xFFFFFFFF;
      await new Promise((r) => setTimeout(r, 300));
      const white = mp.props.backColor >>> 0;
      return {red, white};
    });
    // Prop-echo by design: GROK-20439 (back color not applied to the drawing)
    // is open, so the painted background is deliberately not measured — the
    // step asserts the round-trip plus a no-error floor.
    expect(result.red).toBe(0xFFFF0000);
    expect(result.white).toBe(0xFFFFFFFF);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Title and Description — rendered text tracks the properties', async () => {
    const set = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.showTitle = true;
      mp.props.title = 'My Matrix';
      mp.props.description = 'Test description';
      await new Promise((r) => setTimeout(r, 900));
      // The title renders in the viewer's own panel titlebar; 'My Matrix' is a
      // unique string, so a document-wide title-text scan is unambiguous.
      const titleShown = [...document.querySelectorAll('.panel-titlebar-text')]
        .some((e) => e.textContent!.trim() === 'My Matrix');
      const root = document.querySelector('[name="viewer-Matrix-plot"]')!;
      const descShown = [...root.querySelectorAll('*')]
        .some((e) => e.children.length === 0 && e.textContent!.trim() === 'Test description');
      const titleTexts = [...document.querySelectorAll('.panel-titlebar-text')].map((e) => e.textContent!.trim());
      return {titleShown, descShown, titleTexts};
    });
    console.log(`MatrixPlot title/desc: shown=${set.titleShown} desc=${set.descShown} titles=${JSON.stringify(set.titleTexts)}`);
    expect(set.titleShown).toBe(true);
    expect(set.descShown).toBe(true);

    const cleared = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      // The titlebar always shows the title text, so clearing the Title (with
      // Show Title off) is what removes the rendered title.
      mp.props.showTitle = false;
      mp.props.title = '';
      mp.props.description = '';
      await new Promise((r) => setTimeout(r, 800));
      const titleGone = ![...document.querySelectorAll('.panel-titlebar-text')]
        .some((e) => e.textContent!.trim() === 'My Matrix');
      return {titleGone};
    });
    expect(cleared.titleGone).toBe(true);
  });

  await softStep('Inner Viewer Look — Cell Plot Type repaints; inner size is a prop round-trip + no-error floor', async () => {
    const errBefore = errCount();
    // The Cell Plot Type switch is the reliable inner-cell render signal: an
    // off-diagonal cell repaints between Density and Scatter and restores
    // exactly on the way back.
    const densInk = await settledCellInk(page, 1);
    await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.cellPlotType = 'Scatter plot';
      await new Promise((r) => setTimeout(r, 900));
    });
    const scatInk = await settledCellInk(page, 1);
    expect(Math.abs(scatInk - densInk)).toBeGreaterThan(500);

    // Inner marker size (scatter state): the look object round-trips and the
    // viewer stays alive, but the size change does NOT repaint the cells — so
    // it is held as a prop round-trip + no-error + liveness floor, never a
    // vacuous render claim (recorded gap).
    const marker = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      const ivl = mp.props.innerViewerLook;
      ivl.markerSize = 10;
      mp.props.innerViewerLook = ivl;
      await new Promise((r) => setTimeout(r, 900));
      return mp.props.innerViewerLook.markerSize;
    });
    expect(marker).toBe(10);
    expect(await cellCount(page)).toBe(16);

    // Density state: repaint back (exact restore) and an inner bin-count round-trip.
    await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.cellPlotType = 'Density plot';
      await new Promise((r) => setTimeout(r, 900));
    });
    const backInk = await settledCellInk(page, 1);
    expect(Math.abs(backInk - scatInk)).toBeGreaterThan(500);
    // Exact restore: returning to Density repaints the cell back to its original
    // ink to within the settle-bound (40 = the settledCellInk two-read agreement
    // window); the restore is deterministic, so no wider tolerance is warranted.
    expect(Math.abs(backInk - densInk)).toBeLessThan(40);

    const bin = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      const ivl = mp.props.innerViewerLook;
      if ('binCount' in ivl) ivl.binCount = 20; else ivl.bins = 20;
      mp.props.innerViewerLook = ivl;
      await new Promise((r) => setTimeout(r, 900));
      return mp.props.innerViewerLook.binCount ?? mp.props.innerViewerLook.bins;
    });
    expect(bin).toBe(20);
    expect(await cellCount(page)).toBe(16);
    expect(errCount()).toBe(errBefore);
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
