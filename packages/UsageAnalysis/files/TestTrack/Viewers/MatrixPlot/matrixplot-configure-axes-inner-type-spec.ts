/* ---
realizes: [matrixplot.cp.configure-axes-inner-type, matrixplot.int.axes-drive-inner-grid, viewers.matrix-plot]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep, isLocalBootNoise} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;

// Scenarios 1-4 of the configure scenario. The layout and project round-trips (Scenario 5)
// live in matrixplot-configure-axes-inner-type-server-spec.ts on the server lane.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isBenignError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text);

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

const columnLabels = (page: Page) => page.evaluate(() => {
  const root = document.querySelector('[name="viewer-Matrix-plot"]')!;
  return [...root.querySelectorAll('div')]
    .filter((e) => e.children.length === 0 && e.textContent!.trim().length > 0
      && getComputedStyle(e).writingMode === 'horizontal-tb')
    .map((e) => e.textContent!.trim());
});

const readSets = (page: Page) => page.evaluate(() => {
  const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
  return {x: mp.props.xColumnNames, y: mp.props.yColumnNames, cellPlotType: mp.props.cellPlotType};
});

async function setSets(page: Page, x: string[] | null, y: string[] | null, settleMs = 900) {
  const set: Record<string, string[]> = {};
  if (x) set.xColumnNames = x;
  if (y) set.yColumnNames = y;
  await v.setViewerProps(page, 'Matrix plot', [{set}], 200);
  await v.waitForViewerRendered(page, 'Matrix plot', settleMs);
}

async function setCellPlotType(page: Page, value: string) {
  await page.evaluate(() => {
    const label = document.querySelector('[name="prop-view-cell-plot-type"]') as HTMLElement | null;
    if (!label) throw new Error('cell-plot-type label not found');
    label.scrollIntoView({block: 'center'});
    label.click();
  });

  await v.pollValue(() => page.locator('select.property-grid-item-editor-spinner').count(),
    (n) => n > 0, 300, 100);
  await page.evaluate((v: string) => {
    const sel = document.querySelector('select.property-grid-item-editor-spinner') as HTMLSelectElement | null;
    if (!sel) throw new Error('cell-plot-type select editor not found');
    sel.value = v;
    sel.dispatchEvent(new Event('change', {bubbles: true}));
  }, value);

  await v.waitForViewerRendered(page, 'Matrix plot', 800);
}

// The gear's click handler is AppEvents.setCurrentObject(viewer), which is what grok.shell.o
// does; the gear itself is not the subject here and a click on it, synthetic or real, opens
// nothing on a shared local page in ~1 run in 3 (measured 2026-09-03). The panel is proven
// bound to THIS viewer, since the previous spec's grid outlives its closed viewer.
async function openGear(page: Page) {
  const bound = () => page.evaluate(() => {
    const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot');
    return !!document.querySelector('.property-grid [name="prop-view-x"]') && !!mp && grok.shell.o?.dart === mp.dart;
  });
  await page.evaluate(() => { grok.shell.o = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot'); });
  if (await v.pollValue(bound, (b) => b, 3000, 100)) return;
  await v.clickViewerTitlebarIcon(page, 'Matrix-plot', 'icon-font-icon-settings');
  expect(await v.pollValue(bound, (b) => b, 3000, 100)).toBe(true);
}

test('Matrix Plot — Column Sets, Cell Plot Type', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e)) && !isLocalBootNoise(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text()) && !isLocalBootNoise(m.text())) consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await openDatagrok(page);
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'matrix-plot', 'Matrix-plot');

  await softStep('Scenario 1 — default state: 16 cells, demog auto-pick, Density plot', async () => {
    const cells = await cellCount(page);
    const sets = await readSets(page);
    expect(cells).toBe(16);
    expect(sets.x).toEqual(['AGE', 'HEIGHT', 'WEIGHT', 'STARTED']);
    expect(sets.y).toEqual(['AGE', 'HEIGHT', 'WEIGHT', 'STARTED']);
    expect(sets.cellPlotType).toBe('Density plot');
  });

  await openGear(page);

  await softStep('Scenario 2 — change X via the Select columns dialog (real clicks; GROK-20438 labels)', async () => {
    await page.locator('[name="prop-view-x"] button').click();
    await page.locator('[name^="dialog-Select-columns"]').waitFor({timeout: 8000});
    const rect = await page.evaluate(() => {
      const g = document.querySelector('[name^="dialog-Select-columns"] [name="viewer-Grid"]')!;
      const r = g.getBoundingClientRect();
      return {top: r.top, right: r.right, height: r.height};
    });

    const headerH = 24;
    const rowH = (rect.height - headerH) / 4;
    const clickX = rect.right - 38;
    const rowY = (i: number) => rect.top + headerH + rowH * i + rowH / 2;

    await page.mouse.click(clickX, rowY(2));
    await v.pollValue(() => readSets(page), (s) => s.x.length === 3, 500, 100);
    await page.mouse.click(clickX, rowY(3));
    await v.pollValue(() => readSets(page), (s) => s.x.length === 2, 500, 100);

    const live = await readSets(page);
    expect(live.x.length).toBeGreaterThan(0);
    await page.locator('[name^="dialog-Select-columns"] [name="button-OK"]').click();
    await v.waitForViewerRendered(page, 'Matrix plot', 900);

    const sets = await v.pollValue(() => readSets(page), (s) => s.x.length === 2, 3000, 150);
    expect(sets.x).toEqual(['AGE', 'HEIGHT']);
    const cells = await v.pollValue(() => cellCount(page), (n) => n === 8, 3000, 150);
    expect(cells).toBe(8);

    const labels = new Set(await v.pollValue(() => columnLabels(page),
      (l) => l.includes('AGE') && l.includes('HEIGHT') && !l.includes('WEIGHT') && !l.includes('STARTED'),
      3000, 150));
    expect(labels.has('AGE')).toBe(true);
    expect(labels.has('HEIGHT')).toBe(true);
    expect(labels.has('WEIGHT')).toBe(false);
    expect(labels.has('STARTED')).toBe(false);
  });

  await softStep('Scenario 3 — cycle the column sets: no new console/page error (GROK-16473)', async () => {
    const errBefore = errCount();
    await setSets(page, ['AGE', 'HEIGHT', 'WEIGHT'], null);
    await setSets(page, null, ['AGE', 'HEIGHT']);
    await setSets(page, ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED'], ['AGE', 'HEIGHT', 'WEIGHT', 'STARTED']);

    const cells = await v.pollValue(() => cellCount(page), (n) => n === 16, 3000, 150);
    expect(cells).toBe(16);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 4 — switch Cell Plot Type: off-diagonal cell repaints', async () => {
    const densInk = await settledCellInk(page, 1);
    await setCellPlotType(page, 'Scatter plot');
    const scatInk = await settledCellInk(page, 1);
    await setCellPlotType(page, 'Density plot');
    const backInk = await settledCellInk(page, 1);
    console.log(`MatrixPlot cellPlotType ink: dens=${densInk} scat=${scatInk} back=${backInk}`);
    expect(densInk).toBeGreaterThan(0);
    expect(scatInk).toBeGreaterThan(0);

    expect(Math.abs(scatInk - densInk)).toBeGreaterThan(500);
    expect(Math.abs(backInk - scatInk)).toBeGreaterThan(500);
  });

  await v.closeAllAndWait(page);
  v.finishSpec();
});
