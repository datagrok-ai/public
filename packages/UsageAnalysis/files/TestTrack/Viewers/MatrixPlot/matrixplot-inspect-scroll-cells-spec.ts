/* ---
realizes: [matrixplot.cp.inspect-scroll-cells, matrixplot.int.large-matrix-scroll-viewport, matrixplot.int.cell-tooltip-then-open-fullscreen, viewers.matrix-plot]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {INNER_CELLS, cellCount, settledCellInks, stampRenders} from './matrix-helpers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const BASE_COLS = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];

const isBenignError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Stack trace [A-Za-z]+/.test(text);

// One drag of a viewport slider's max handle from `from` visible columns to `to`. Every
// mousemove that changes the integer viewport re-tiles the whole visible matrix (0.1s at
// 25 cells, 3s at 240), so the gesture is a single move to the computed position rather
// than a paced sweep; the value is proven by the cell count the caller asserts.
async function dragMaxTo(page: Page, slider: 'x-slider' | 'y-slider', from: number, to: number,
  expectRender: boolean): Promise<number> {
  return page.evaluate(async (a: {slider: string; from: number; to: number; expectRender: boolean; cells: string}) => {
    const w = window as any;
    const svg = document.querySelector(`[name="viewer-Matrix-plot"] svg[name="${a.slider}"]`)!;
    const min = svg.querySelector('[name="min-handle"]')!.getBoundingClientRect();
    const max = svg.querySelector('[name="max-handle"]') as HTMLElement;
    const mb = max.getBoundingClientRect();
    const vertical = a.slider === 'y-slider';
    const centre = (b: DOMRect) => vertical ? b.y + b.height / 2 : b.x + b.width / 2;
    // the handle centres sit one handle diameter apart beyond the value span
    const pxPerUnit = (centre(mb) - centre(min) - mb.width) / a.from;
    const d = (a.to - a.from) * pxPerUnit;
    const cx = mb.x + mb.width / 2, cy = mb.y + mb.height / 2;
    const act = () => w.__drag(max, {x: cx, y: cy}, {x: vertical ? cx : cx + d, y: vertical ? cy + d : cy},
      {steps: 1, stepMs: 0, holdMs: 0});
    if (a.expectRender) await w.__settled('viewer:Matrix plot.onViewerRendered', act, 5000);
    else await act();
    return document.querySelectorAll(a.cells).length;
  }, {slider, from, to, expectRender, cells: INNER_CELLS});
}

const topLabelSet = (page: Page) => page.evaluate(() => {
  const root = document.querySelector('[name="viewer-Matrix-plot"]')!;
  const names = new Set(['AGE', 'HEIGHT', 'WEIGHT', 'STARTED']);

  return [...new Set([...root.querySelectorAll('div')]
    .filter((e) => e.children.length === 0 && names.has(e.textContent!.trim())
      && getComputedStyle(e).writingMode === 'horizontal-tb')
    .map((e) => e.textContent!.trim()))];
});

async function setSets(page: Page, cols: string[], settleMs = 1200) {
  await v.setViewerProps(page, 'Matrix plot', [{set: {xColumnNames: cols, yColumnNames: cols}, wait: 300}]);
  await v.waitForViewerRendered(page, 'Matrix plot', settleMs);
}

test('Matrix Plot — Viewport Scrolling and Cell Inspection', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error' && !isBenignError(m.text()) && !isLocalBootNoise(m.text())) consoleErrors.push(m.text()); });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  const fixtureCols = await page.evaluate(() => {
    const df = grok.shell.tv.dataFrame;
    const names: string[] = [];
    for (let i = 1; i <= 12; i++) {
      const nm = 'MP_FX_' + i;
      if (!df.columns.names().includes(nm)) df.columns.addNewCalculated(nm, '${AGE} + ' + i);
      names.push(nm);
    }
    return names;
  });
  await v.pollValue(() => page.evaluate((names: string[]) => {
    const have = grok.shell.tv.dataFrame.columns.names();
    return names.every((nm) => have.includes(nm));
  }, fixtureCols), (ready) => ready, 800, 100);

  try {
    await v.addViewerByIcon(page, 'matrix-plot', 'Matrix-plot');
    await stampRenders(page);
    await setSets(page, BASE_COLS);

    await softStep('Scenario 1 — move the viewport with the sliders and hit the 250-cell cap', async () => {
      expect(await page.evaluate(() =>
        !!document.querySelector('[name="viewer-Matrix-plot"] svg[name="x-slider"]') &&
        !!document.querySelector('[name="viewer-Matrix-plot"] svg[name="y-slider"]'))).toBe(true);
      expect(await v.pollValue(() => cellCount(page), (n) => n === 16, 3000, 150)).toBe(16);

      const mid = await dragMaxTo(page, 'x-slider', 4, 2, true);
      expect(mid).toBeLessThan(16);
      expect(mid).toBeGreaterThanOrEqual(4);
      const midLabels = await topLabelSet(page);
      expect(midLabels.length).toBeLessThan(4);
      expect(midLabels.length).toBeGreaterThanOrEqual(1);

      const back = await dragMaxTo(page, 'x-slider', 2, 4, true);
      expect(back).toBe(16);
      expect((await topLabelSet(page)).length).toBe(4);

      const allCols = [...BASE_COLS, ...fixtureCols];
      await setSets(page, allCols, 1500);

      const initial = await v.pollValue(() => cellCount(page), (n) => n > 0, 3000, 150);
      expect(initial).toBeLessThan(80);

      const xFull = await dragMaxTo(page, 'x-slider', 5, 16, true);
      expect(xFull).toBeGreaterThan(initial);

      const yNearCap = await dragMaxTo(page, 'y-slider', 5, 15, true);
      // the 16th row would make 256 cells: the increment is rejected and nothing repaints
      await dragMaxTo(page, 'y-slider', 15, 16, false);
      const yFull = await v.pollValue(() => cellCount(page), (n) => n !== yNearCap, 1000, 100);
      const samples = [initial, xFull, yNearCap, yFull];
      const maxSeen = Math.max(...samples);
      console.log(`MatrixPlot cap: initial=${initial} xFull=${xFull} yNearCap=${yNearCap} yFull=${yFull} maxSeen=${maxSeen}`);

      expect(yFull).toBeGreaterThan(xFull);
      expect(yFull).toBe(yNearCap);
      expect(maxSeen).toBeLessThanOrEqual(250);
      expect(samples.includes(256)).toBe(false);

      await setSets(page, BASE_COLS, 1500);
      expect(await v.pollValue(() => cellCount(page), (n) => n === 16, 3000, 150)).toBe(16);
    });

    await softStep('Scenario 2 — hover a cell: tooltip identity and expand icon reveal', async () => {
      const visBefore = await page.evaluate(() => {
        const cells = document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer');
        const cell = cells[1] as HTMLElement;
        const iconBefore = cell.parentElement!.querySelector('[name="icon-expand-arrows"]') as HTMLElement | null;
        const visBefore = iconBefore ? getComputedStyle(iconBefore).visibility : null;
        const r = cell.getBoundingClientRect();
        const cx = r.x + r.width / 2, cy = r.y + r.height / 2;
        cell.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, clientX: cx, clientY: cy}));
        cell.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: cx, clientY: cy}));
        return visBefore;
      });
      const hovered = await v.pollValue(() => page.evaluate(() => {
        const cells = document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer');
        const parent = (cells[1] as HTMLElement).parentElement!;
        const tt = document.querySelector('.d4-tooltip') as HTMLElement | null;
        const ttWidth = tt ? tt.getBoundingClientRect().width : 0;
        const ttText = (tt && ttWidth > 0) ? tt.innerText : null;
        const iconAfter = parent.querySelector('[name="icon-expand-arrows"]') as HTMLElement | null;
        const visAfter = iconAfter ? getComputedStyle(iconAfter).visibility : null;
        return {visAfter, ttWidth, ttText};
      }), (h) => h.ttWidth > 0 && h.visAfter === 'visible', 600, 100);
      const result = {visBefore, ...hovered};

      expect(result.ttWidth).toBeGreaterThan(0);
      expect(result.ttText).toMatch(/HEIGHT/);
      expect(result.ttText).toMatch(/AGE/);

      expect(result.visBefore).toBe('hidden');
      expect(result.visAfter).toBe('visible');
    });

    await softStep('Scenario 3 — open a cell as a matching standalone viewer', async () => {
      const before = await page.evaluate(() => grok.shell.tv.viewers.map((vw: any) => vw.type));

      const viewerTypes = () => page.evaluate(() => grok.shell.tv.viewers.map((vw: any) => vw.type) as string[]);
      const expandCell = async (idx: number) => {
        const countBefore = (await viewerTypes()).length;
        await page.evaluate((i: number) => {
          const cells = document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer');
          const cell = cells[i] as HTMLElement;
          const r = cell.getBoundingClientRect();
          const cx = r.x + r.width / 2, cy = r.y + r.height / 2;
          cell.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, clientX: cx, clientY: cy}));
          cell.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: cx, clientY: cy}));
        }, idx);
        await v.pollValue(() => page.evaluate((i: number) => {
          const cells = document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer');
          const icon = (cells[i] as HTMLElement).parentElement!
            .querySelector('[name="icon-expand-arrows"]') as HTMLElement | null;
          return icon ? getComputedStyle(icon).visibility : null;
        }, idx), (vis) => vis === 'visible', 500, 100);
        await page.evaluate((i: number) => {
          const cells = document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer');
          const icon = (cells[i] as HTMLElement).parentElement!
            .querySelector('[name="icon-expand-arrows"]') as HTMLElement | null;
          if (icon) icon.click();
        }, idx);
        return v.pollValue(viewerTypes, (types) => types.length > countBefore, 1200, 100);
      };
      const closeAdded = async () => {
        await page.evaluate(() => {
          for (const vw of grok.shell.tv.viewers.slice())
            if (vw.type !== 'Grid' && vw.type !== 'Matrix plot') vw.close();
        });
        await v.pollValue(viewerTypes,
          (types) => types.every((t) => t === 'Grid' || t === 'Matrix plot'), 600, 100);
      };

      const afterOff = await expandCell(1);
      expect(afterOff).toContain('Density plot');
      await closeAdded();

      const afterDiag = await expandCell(0);
      expect(afterDiag).toContain('Histogram');
      await closeAdded();

      const finalSet = await page.evaluate(() => grok.shell.tv.viewers.map((vw: any) => vw.type));
      expect(finalSet.sort()).toEqual([...before].sort());
    });

    await softStep('Scenario 4 — wheel zoom is per-cell', async () => {
      const errBefore = errCount();
      const [targetBase, neighbourBase] = await settledCellInks(page, [1, 2]);
      const wheel = (idx: number, dy: number) => page.evaluate((args: {i: number; dy: number}) => {
        const cells = document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer');
        const c = cells[args.i] as HTMLElement;
        const r = c.getBoundingClientRect();
        c.dispatchEvent(new WheelEvent('wheel', {
          bubbles: true, cancelable: true,
          clientX: r.x + r.width / 2, clientY: r.y + r.height / 2, deltaY: args.dy,
        }));
      }, {i: idx, dy});

      await wheel(1, -300);
      await v.waitForViewerRendered(page, 'Matrix plot', 400);
      const [targetZoom, neighbourZoom] = await settledCellInks(page, [1, 2]);
      console.log(`MatrixPlot wheel zoom: target ${targetBase}->${targetZoom} neighbour ${neighbourBase}->${neighbourZoom}`);

      expect(Math.abs(targetZoom - targetBase)).toBeGreaterThan(300);
      expect(Math.abs(neighbourZoom - neighbourBase)).toBeLessThan(100);

      await wheel(1, 300);
      await v.waitForViewerRendered(page, 'Matrix plot', 400);
      await settledCellInks(page, [1]);
      expect(errCount()).toBe(errBefore);
    });
  } finally {
    await page.evaluate((names: string[]) => {
      const df = grok.shell.tv?.dataFrame;
      if (df) for (const nm of names) if (df.columns.names().includes(nm)) df.columns.remove(nm);
    }, fixtureCols);
    await v.closeAllAndWait(page);
  }

  v.finishSpec();
});
