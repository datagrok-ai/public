/* ---
realizes: [matrixplot.cp.configure-axes-inner-type, matrixplot-axes-drive-inner-grid, viewers.matrix-plot]
--- */
import {test, expect, Page} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';

const isAmbientError = (text: string) =>
  /WebSocket/.test(text) || /Failed to load resource/.test(text) || /404 \(\)/.test(text) ||
  /favicon/.test(text);

// NullError / "Stack trace" / cloned-iframe messages are benign ONLY during the
// ribbon project save; elsewhere that same Dart class is what the error guards catch.
let inProjectSaveWindow = false;
const isBenignError = (text: string) => {
  if (isAmbientError(text)) return true;
  if (inProjectSaveWindow)
    return /Unable to find element in cloned iframe/.test(text) ||
      /Stack trace [A-Za-z]+/.test(text) ||
      /NullError: method not found: '\w+' on null/.test(text);
  return false;
};

/** Settle-gated non-white ink over one inner cell canvas (row-major index). */
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

/** Re-read one cell's ink until two consecutive readings agree (settle the render tail). */
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

/** Leaf DOM X-axis label divs above the matrix (the top strip). The X labels
 * render horizontally (writing-mode horizontal-tb); the Y labels on the left
 * edge are vertical (vertical-rl) and are excluded — GROK-20438 is about the
 * X-axis labels matching the chosen X set. */
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
  await page.evaluate(async (args: {x: string[] | null; y: string[] | null}) => {
    const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
    if (args.x) mp.props.xColumnNames = args.x;
    if (args.y) mp.props.yColumnNames = args.y;
    await new Promise((r) => setTimeout(r, 200));
  }, {x, y});
  await page.waitForTimeout(settleMs);
}

/** Drive the Cell Plot Type through the settings editor (not the JS-API prop):
 * clicking the label swaps in a select editor; set `.value` + dispatch `change`. */
async function setCellPlotType(page: Page, value: string) {
  await page.evaluate(() => {
    const label = document.querySelector('[name="prop-view-cell-plot-type"]') as HTMLElement | null;
    if (!label) throw new Error('cell-plot-type label not found');
    label.scrollIntoView({block: 'center'});
    label.click();
  });
  await page.waitForTimeout(300);
  await page.evaluate((v: string) => {
    const sel = document.querySelector('select.property-grid-item-editor-spinner') as HTMLSelectElement | null;
    if (!sel) throw new Error('cell-plot-type select editor not found');
    sel.value = v;
    sel.dispatchEvent(new Event('change', {bubbles: true}));
  }, value);
  await page.waitForTimeout(800);
}

test('Matrix Plot — Column Sets, Cell Plot Type, Persistence', async ({page}: {page: Page}) => {
  test.setTimeout(600_000);

  const pageErrors: string[] = [];
  page.on('pageerror', (e) => { if (!isBenignError(String(e))) pageErrors.push(String(e)); });
  const consoleErrors: string[] = [];
  page.on('console', (m) => {
    if (m.type() === 'error' && !isBenignError(m.text())) consoleErrors.push(m.text());
  });
  const errCount = () => pageErrors.length + consoleErrors.length;

  await loginToDatagrok(page);
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

  // Open the settings panel for the dialog + cell-plot-type rows.
  await v.openViewerGear(page, 'Matrix plot');

  await softStep('Scenario 2 — change X via the Select columns dialog (real clicks; GROK-20438 labels)', async () => {
    // The "..." multi-column editor opens ONLY on a real (trusted) click.
    await page.locator('[name="prop-view-x"] button').click();
    await page.locator('[name^="dialog-Select-columns"]').waitFor({timeout: 8000});
    const rect = await page.evaluate(() => {
      const g = document.querySelector('[name^="dialog-Select-columns"] [name="viewer-Grid"]')!;
      const r = g.getBoundingClientRect();
      return {top: r.top, right: r.right, height: r.height};
    });
    // Dialog grid geometry: a single column-header band above four equal data
    // rows (AGE, HEIGHT, WEIGHT, STARTED in dataframe order); the header height
    // and row pitch are derived from the grid rect below. The checkbox is the
    // LAST grid column near the right edge and toggles only on a REAL coordinate
    // click.
    const headerH = 24;
    const rowH = (rect.height - headerH) / 4;
    const clickX = rect.right - 38;
    const rowY = (i: number) => rect.top + headerH + rowH * i + rowH / 2;
    // Uncheck WEIGHT (row 2) and STARTED (row 3) — never empty the set.
    await page.mouse.click(clickX, rowY(2));
    await page.waitForTimeout(500);
    await page.mouse.click(clickX, rowY(3));
    await page.waitForTimeout(500);
    // Toggles live-apply before OK; guard against ever committing an empty set.
    const live = await readSets(page);
    expect(live.x.length).toBeGreaterThan(0);
    await page.locator('[name^="dialog-Select-columns"] [name="button-OK"]').click();
    await page.waitForTimeout(900);

    const sets = await readSets(page);
    expect(sets.x).toEqual(['AGE', 'HEIGHT']);
    const cells = await cellCount(page);
    expect(cells).toBe(8); // 2 X columns x 4 Y columns
    // GROK-20438: the rendered DOM labels must match the chosen X set.
    const labels = new Set(await columnLabels(page));
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
    const cells = await cellCount(page);
    expect(cells).toBe(16);
    expect(errCount()).toBe(errBefore);
  });

  await softStep('Scenario 4 — switch Cell Plot Type: off-diagonal cell repaints', async () => {
    // The matrix is a full 4x4 here; cell index 1 is the off-diagonal
    // X=HEIGHT x Y=AGE cell. Density and Scatter paint this cell with clearly
    // different ink coverage, so the repaint is asserted against a coverage floor.
    const densInk = await settledCellInk(page, 1);
    await setCellPlotType(page, 'Scatter plot');
    const scatInk = await settledCellInk(page, 1);
    await setCellPlotType(page, 'Density plot');
    const backInk = await settledCellInk(page, 1);
    console.log(`MatrixPlot cellPlotType ink: dens=${densInk} scat=${scatInk} back=${backInk}`);
    expect(densInk).toBeGreaterThan(0);
    expect(scatInk).toBeGreaterThan(0);
    // Density<->Scatter is a strong off-diagonal repaint; the 500-non-white-px
    // floor sits far above the render-tail settle jitter (<40 px) while staying
    // well under the observed swing, giving a wide margin in each direction.
    expect(Math.abs(scatInk - densInk)).toBeGreaterThan(500);
    expect(Math.abs(backInk - scatInk)).toBeGreaterThan(500);
  });

  await softStep('Scenario 5a — layout round-trip restores the saved viewer set and config', async () => {
    const layoutId = await page.evaluate(async () => {
      const mp = grok.shell.tv.viewers.find((vw: any) => vw.type === 'Matrix plot') as any;
      mp.props.xColumnNames = ['AGE', 'HEIGHT'];
      mp.props.yColumnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      mp.props.cellPlotType = 'Scatter plot';
      await new Promise((r) => setTimeout(r, 900));
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      return layout.id as string;
    });
    try {
      const result = await page.evaluate(async (id) => {
        const tv = grok.shell.tv;
        await new Promise((r) => setTimeout(r, 1000));
        tv.addViewer('Scatter plot');
        await new Promise((r) => setTimeout(r, 800));
        const saved = await grok.dapi.layouts.find(id);
        tv.loadLayout(saved);
        await new Promise((r) => setTimeout(r, 3000));
        const hasScatter = tv.viewers.some((vw: any) => vw.type === 'Scatter plot');
        const hasMatrix = tv.viewers.some((vw: any) => vw.type === 'Matrix plot');
        const mp = tv.viewers.find((vw: any) => vw.type === 'Matrix plot');
        return {
          hasScatter, hasMatrix,
          x: mp?.props.xColumnNames, y: mp?.props.yColumnNames, cellPlotType: mp?.props.cellPlotType,
        };
      }, layoutId);
      // The restored set equals the SAVED set: Matrix present, later Scatter absent.
      expect(result.hasMatrix).toBe(true);
      expect(result.hasScatter).toBe(false);
      expect(result.x).toEqual(['AGE', 'HEIGHT']);
      expect(result.y).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
      expect(result.cellPlotType).toBe('Scatter plot');
    } finally {
      await page.evaluate(async (id) => {
        try {
          const saved = await grok.dapi.layouts.find(id);
          if (saved) await grok.dapi.layouts.delete(saved);
        } catch (_) {}
      }, layoutId);
    }
  });

  await softStep('Scenario 5b — project save / Close All / reopen restores the Matrix plot (GROK-10925)', async () => {
    const projName = 'zz-matrixplot-persistence-probe-' + Date.now();
    let projectId: string | null = null;
    inProjectSaveWindow = true;
    try {
      const saved = await saveProjectViaUI(page, projName);
      projectId = saved.projectId;
      expect(projectId).toBeTruthy();

      const result = await page.evaluate(async (id) => {
        grok.shell.closeAll();
        await new Promise((r) => setTimeout(r, 1500));
        const full = await grok.dapi.projects.find(id);
        await full.open();
        let types: string[] = [];
        for (let t = 0; t < 20; t++) {
          await new Promise((r) => setTimeout(r, 1000));
          types = [];
          for (const view of grok.shell.tableViews)
            for (const vw of view.viewers) types.push(vw.type);
          if (types.includes('Matrix plot')) break;
        }
        let mp: any = null;
        for (const view of grok.shell.tableViews)
          for (const vw of view.viewers)
            if (vw.type === 'Matrix plot') mp = vw;
        return {
          types,
          x: mp?.props.xColumnNames, y: mp?.props.yColumnNames, cellPlotType: mp?.props.cellPlotType,
        };
      }, projectId);

      expect(result.types).toContain('Matrix plot');
      expect(result.x).toEqual(['AGE', 'HEIGHT']);
      expect(result.y).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
      expect(result.cellPlotType).toBe('Scatter plot');
    } finally {
      inProjectSaveWindow = false;
      if (projectId) await deleteProjectWithCleanup(page, {projectId});
    }
  });

  await page.evaluate(() => grok.shell.closeAll());
  v.finishSpec();
});
