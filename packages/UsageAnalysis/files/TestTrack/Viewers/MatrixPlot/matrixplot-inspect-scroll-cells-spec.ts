/* ---
realizes: [matrixplot.cp.inspect-scroll-cells, matrixplot-large-matrix-scroll-viewport, matrixplot-cell-tooltip-then-open-fullscreen, viewers.matrix-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const BASE_COLS = ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'];

const isBenignError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text) || /Stack trace [A-Za-z]+/.test(text);

const cellCount = (page: Page) => page.evaluate(() =>
  document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer').length);

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

/** Drag a viewport slider's max-handle to a fraction of its track (synthetic
 * mousedown -> document-mousemove steps -> mouseup drives it — same widget
 * family as the PC Plot axis sliders). Returns the final visible cell count. */
async function dragSliderMax(page: Page, slider: 'x-slider' | 'y-slider', axis: 'x' | 'y', fraction: number): Promise<number> {
  return page.evaluate(async (args: {slider: string; axis: string; fraction: number}) => {
    const root = document.querySelector('[name="viewer-Matrix-plot"]')!;
    const svg = root.querySelector(`svg[name="${args.slider}"]`)!;
    const sr = svg.getBoundingClientRect();
    const max = svg.querySelector('[name="max-handle"]')!;
    const b = max.getBoundingClientRect();
    let cx = b.x + b.width / 2, cy = b.y + b.height / 2;
    const mk = (x: number, y: number) => ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
    max.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
    await new Promise((r) => setTimeout(r, 50));
    const cellN = () => document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer').length;
    if (args.axis === 'x') {
      const target = sr.x + sr.width * args.fraction - 1;
      const step = target < cx ? -15 : 15;
      for (let x = cx; (step < 0 ? x >= target : x <= target); x += step) {
        document.dispatchEvent(new MouseEvent('mousemove', mk(x, cy)));
        svg.dispatchEvent(new MouseEvent('mousemove', mk(x, cy)));
        await new Promise((r) => setTimeout(r, 20));
      }
      document.dispatchEvent(new MouseEvent('mouseup', mk(target, cy)));
    } else {
      const target = sr.y + sr.height * args.fraction - 1;
      const step = target < cy ? -15 : 15;
      for (let y = cy; (step < 0 ? y >= target : y <= target); y += step) {
        document.dispatchEvent(new MouseEvent('mousemove', mk(cx, y)));
        svg.dispatchEvent(new MouseEvent('mousemove', mk(cx, y)));
        await new Promise((r) => setTimeout(r, 20));
      }
      document.dispatchEvent(new MouseEvent('mouseup', mk(cx, target)));
    }
    await new Promise((r) => setTimeout(r, 800));
    return cellN();
  }, {slider, axis, fraction});
}

/** Drag the max-handle fully open while sampling the visible cell count each
 * step — used to prove the 250-cell cap never lets a sample reach 256. */
async function dragFullSampled(page: Page, slider: 'x-slider' | 'y-slider', axis: 'x' | 'y'): Promise<{samples: number[]; final: number}> {
  return page.evaluate(async (args: {slider: string; axis: string}) => {
    const root = document.querySelector('[name="viewer-Matrix-plot"]')!;
    const svg = root.querySelector(`svg[name="${args.slider}"]`)!;
    const sr = svg.getBoundingClientRect();
    const max = svg.querySelector('[name="max-handle"]')!;
    const b = max.getBoundingClientRect();
    const cx = b.x + b.width / 2, cy = b.y + b.height / 2;
    const mk = (x: number, y: number) => ({bubbles: true, cancelable: true, clientX: x, clientY: y, button: 0});
    const cellN = () => document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer').length;
    const samples: number[] = [];
    max.dispatchEvent(new MouseEvent('mousedown', mk(cx, cy)));
    await new Promise((r) => setTimeout(r, 50));
    if (args.axis === 'x') {
      const end = sr.x + sr.width - 1;
      for (let x = cx; x <= end; x += 15) {
        document.dispatchEvent(new MouseEvent('mousemove', mk(x, cy)));
        svg.dispatchEvent(new MouseEvent('mousemove', mk(x, cy)));
        await new Promise((r) => setTimeout(r, 25));
        samples.push(cellN());
      }
      document.dispatchEvent(new MouseEvent('mouseup', mk(end, cy)));
    } else {
      const end = sr.y + sr.height - 1;
      for (let y = cy; y <= end; y += 15) {
        document.dispatchEvent(new MouseEvent('mousemove', mk(cx, y)));
        svg.dispatchEvent(new MouseEvent('mousemove', mk(cx, y)));
        await new Promise((r) => setTimeout(r, 25));
        samples.push(cellN());
      }
      document.dispatchEvent(new MouseEvent('mouseup', mk(cx, end)));
    }
    await new Promise((r) => setTimeout(r, 800));
    return {samples, final: cellN()};
  }, {slider, axis});
}

const topLabelSet = (page: Page) => page.evaluate(() => {
  const root = document.querySelector('[name="viewer-Matrix-plot"]')!;
  const names = new Set(['AGE', 'HEIGHT', 'WEIGHT', 'STARTED']);
  // Only the horizontal (top-strip) X labels; the vertical (vertical-rl) Y
  // labels on the left edge are always present and must be excluded.
  return [...new Set([...root.querySelectorAll('div')]
    .filter((e) => e.children.length === 0 && names.has(e.textContent!.trim())
      && getComputedStyle(e).writingMode === 'horizontal-tb')
    .map((e) => e.textContent!.trim()))];
});

async function setSets(page: Page, cols: string[], settleMs = 1200) {
  await page.evaluate(async (c: string[]) => {
    const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
    mp.props.xColumnNames = c;
    mp.props.yColumnNames = c;
    await new Promise((r) => setTimeout(r, 300));
  }, cols);
  await page.waitForTimeout(settleMs);
}

test('Matrix Plot — Viewport Scrolling and Cell Inspection', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => { if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text()); });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await loginToDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});

  // Derive 12 numeric fixture columns (AGE + i, i = 1..12) so the eligible
  // numeric set reaches 16 for the 250-cell cap probe. Removed in finally.
  const fixtureCols = await page.evaluate(async () => {
    const df = grok.shell.tv.dataFrame;
    const names: string[] = [];
    for (let i = 1; i <= 12; i++) {
      const nm = 'MP_FX_' + i;
      if (!df.columns.names().includes(nm)) df.columns.addNewCalculated(nm, '${AGE} + ' + i);
      names.push(nm);
    }
    await new Promise((r) => setTimeout(r, 800));
    return names;
  });

  try {
    await v.addViewerByIcon(page, 'matrix-plot', 'Matrix-plot');
    await setSets(page, BASE_COLS);

    await softStep('Scenario 1 — move the viewport with the sliders and hit the 250-cell cap', async () => {
      expect(await page.evaluate(() =>
        !!document.querySelector('[name="viewer-Matrix-plot"] svg[name="x-slider"]') &&
        !!document.querySelector('[name="viewer-Matrix-plot"] svg[name="y-slider"]'))).toBe(true);
      expect(await cellCount(page)).toBe(16);

      // Drag the x max-handle toward the slider midpoint: the visible X
      // columns narrow (the exact settle count is viewport-dependent, so the
      // assert is "fewer columns visible" rather than a hard-coded 2).
      const mid = await dragSliderMax(page, 'x-slider', 'x', 0.5);
      expect(mid).toBeLessThan(16);
      expect(mid).toBeGreaterThanOrEqual(4);
      const midLabels = await topLabelSet(page);
      expect(midLabels.length).toBeLessThan(4);
      expect(midLabels.length).toBeGreaterThanOrEqual(1);

      // Drag back to the right end: the full 4x4 viewport returns (round-trip).
      const back = await dragSliderMax(page, 'x-slider', 'x', 1.0);
      expect(back).toBe(16);
      expect((await topLabelSet(page)).length).toBe(4);

      // Build the 16x16 set and probe the cap.
      const allCols = [...BASE_COLS, ...fixtureCols];
      await setSets(page, allCols, 1500);
      const initial = await cellCount(page);
      // The re-tiled viewport shows only a small window of the 16x16 grid.
      expect(initial).toBeLessThan(80);

      const xFull = await dragFullSampled(page, 'x-slider', 'x');
      expect(xFull.final).toBeGreaterThan(initial);

      const yFull = await dragFullSampled(page, 'y-slider', 'y');
      const allSamples = [...xFull.samples, ...yFull.samples, xFull.final, yFull.final];
      const maxSeen = Math.max(...allSamples);
      console.log(`MatrixPlot cap: initial=${initial} xFull=${xFull.final} yFull=${yFull.final} maxSeen=${maxSeen}`);
      // The 250-cell cap rejects the final increments: opening both sliders
      // fully never exposes the full 16x16 = 256 grid; the count settles at or
      // below the cap. 256 > 250 makes this assert non-vacuous.
      expect(yFull.final).toBeGreaterThan(xFull.final);
      expect(maxSeen).toBeLessThanOrEqual(250);
      expect(allSamples.includes(256)).toBe(false);

      // Reset to the base 4x4 for the remaining scenarios.
      await setSets(page, BASE_COLS, 1500);
      expect(await cellCount(page)).toBe(16);
    });

    await softStep('Scenario 2 — hover a cell: tooltip identity and expand icon reveal', async () => {
      const result = await page.evaluate(async () => {
        const cells = document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer');
        const cell = cells[1] as HTMLElement; // off-diagonal X=HEIGHT x Y=AGE
        const parent = cell.parentElement!;
        const iconBefore = parent.querySelector('[name="icon-expand-arrows"]') as HTMLElement | null;
        const visBefore = iconBefore ? getComputedStyle(iconBefore).visibility : null;
        const r = cell.getBoundingClientRect();
        const cx = r.x + r.width / 2, cy = r.y + r.height / 2;
        cell.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, clientX: cx, clientY: cy}));
        cell.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: cx, clientY: cy}));
        await new Promise((rr) => setTimeout(rr, 600));
        const tt = document.querySelector('.d4-tooltip') as HTMLElement | null;
        const ttWidth = tt ? tt.getBoundingClientRect().width : 0;
        const ttText = (tt && ttWidth > 0) ? tt.innerText : null;
        const iconAfter = parent.querySelector('[name="icon-expand-arrows"]') as HTMLElement | null;
        const visAfter = iconAfter ? getComputedStyle(iconAfter).visibility : null;
        return {visBefore, visAfter, ttWidth, ttText};
      });
      // The tooltip is read only while its rectangle has a non-zero width
      // (it retains stale text when hidden).
      expect(result.ttWidth).toBeGreaterThan(0);
      expect(result.ttText).toMatch(/HEIGHT/);
      expect(result.ttText).toMatch(/AGE/);
      // The expand icon is always in the DOM, hidden until hover.
      expect(result.visBefore).toBe('hidden');
      expect(result.visAfter).toBe('visible');
    });

    await softStep('Scenario 3 — open a cell as a matching standalone viewer', async () => {
      const before = await page.evaluate(() => grok.shell.tv.viewers.map((vw: any) => vw.type));

      const expandCell = (idx: number) => page.evaluate(async (i: number) => {
        const cells = document.querySelectorAll('[name="viewer-Matrix-plot"] canvas.d4-matrix-plot-inner-viewer');
        const cell = cells[i] as HTMLElement;
        const r = cell.getBoundingClientRect();
        const cx = r.x + r.width / 2, cy = r.y + r.height / 2;
        cell.dispatchEvent(new MouseEvent('mouseenter', {bubbles: true, clientX: cx, clientY: cy}));
        cell.dispatchEvent(new MouseEvent('mousemove', {bubbles: true, clientX: cx, clientY: cy}));
        await new Promise((rr) => setTimeout(rr, 500));
        const icon = cell.parentElement!.querySelector('[name="icon-expand-arrows"]') as HTMLElement | null;
        if (icon) icon.click();
        await new Promise((rr) => setTimeout(rr, 1200));
        return grok.shell.tv.viewers.map((vw: any) => vw.type);
      }, idx);
      const closeAdded = () => page.evaluate(async () => {
        for (const vw of grok.shell.tv.viewers.slice())
          if (vw.type !== 'Grid' && vw.type !== 'Matrix plot') vw.close();
        await new Promise((r) => setTimeout(r, 600));
      });

      // Off-diagonal cell (Cell Plot Type Density plot) -> standalone Density plot.
      const afterOff = await expandCell(1);
      expect(afterOff).toContain('Density plot');
      await closeAdded();

      // Diagonal cell -> standalone Histogram.
      const afterDiag = await expandCell(0);
      expect(afterDiag).toContain('Histogram');
      await closeAdded();

      const finalSet = await page.evaluate(() => grok.shell.tv.viewers.map((vw: any) => vw.type));
      expect(finalSet.sort()).toEqual([...before].sort());
    });

    await softStep('Scenario 4 — wheel zoom is per-cell', async () => {
      const errBefore = errCount();
      const targetBase = await settledCellInk(page, 1);
      const neighbourBase = await settledCellInk(page, 2);
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
      await page.waitForTimeout(400);
      const targetZoom = await settledCellInk(page, 1);
      const neighbourZoom = await settledCellInk(page, 2);
      console.log(`MatrixPlot wheel zoom: target ${targetBase}->${targetZoom} neighbour ${neighbourBase}->${neighbourZoom}`);
      // The zoom repaints ONLY the hovered cell; the neighbour is the control.
      expect(Math.abs(targetZoom - targetBase)).toBeGreaterThan(300);
      expect(Math.abs(neighbourZoom - neighbourBase)).toBeLessThan(100);

      await wheel(1, 300); // zoom back out (approximate)
      await page.waitForTimeout(400);
      await settledCellInk(page, 1);
      expect(errCount()).toBe(errBefore);
    });
  } finally {
    await page.evaluate(async (names: string[]) => {
      const df = grok.shell.tv?.dataFrame;
      if (df) for (const nm of names) if (df.columns.names().includes(nm)) df.columns.remove(nm);
      grok.shell.closeAll();
    }, fixtureCols);
  }

  v.finishSpec();
});
