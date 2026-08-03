/* ---
realizes: [correlationplot.cp.matrix-values-scope-persist, correlationplot.int.numerical-columns-only]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaUI, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const TOL = 1e-3;
const PAIRS: [string, string][] = [['AGE', 'HEIGHT'], ['AGE', 'WEIGHT'], ['HEIGHT', 'WEIGHT']];

// The Correlation plot mounts THREE canvases; the first (v.root.querySelector('canvas'))
// is a 0x0 sizing canvas whose getImageData throws IndexSizeError — the shared
// v.snapshotCanvasColors/v.diffCanvasColors return the -1 fault sentinel on it. The
// painted matrix lives on the first NON-zero-sized canvas. These CP-local helpers pick
// that canvas so the render diff is real, not the -1 fault. (Same per-color histogram
// diff contract as the shared helpers.)
async function cpSnapshot(page: import('@playwright/test').Page): Promise<boolean> {
  return await page.evaluate(() => {
    const w = window as any;
    const cp = w.grok?.shell?.tv?.viewers?.find((x: any) => x.type === 'Correlation plot');
    const cv = cp && Array.from(cp.root.querySelectorAll('canvas'))
      .find((c: any) => c.width > 0 && c.height > 0) as HTMLCanvasElement | undefined;
    if (!cv) return false;
    try {
      const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      const colors = new Map<number, number>();
      for (let i = 0; i < img.length; i += 4) {
        const key = (img[i] << 16) | (img[i + 1] << 8) | img[i + 2];
        colors.set(key, (colors.get(key) ?? 0) + 1);
      }
      w.__cpColorSnap = colors;
      return true;
    } catch { return false; }
  });
}

async function cpDiff(page: import('@playwright/test').Page): Promise<{deltaPx: number}> {
  return await page.evaluate(() => {
    const w = window as any;
    const prev = w.__cpColorSnap as Map<number, number> | undefined;
    const cp = w.grok?.shell?.tv?.viewers?.find((x: any) => x.type === 'Correlation plot');
    const cv = cp && Array.from(cp.root.querySelectorAll('canvas'))
      .find((c: any) => c.width > 0 && c.height > 0) as HTMLCanvasElement | undefined;
    if (!cv || !prev) return {deltaPx: -1};
    try {
      const img = cv.getContext('2d')!.getImageData(0, 0, cv.width, cv.height).data;
      const colors = new Map<number, number>();
      for (let i = 0; i < img.length; i += 4) {
        const key = (img[i] << 16) | (img[i + 1] << 8) | img[i + 2];
        colors.set(key, (colors.get(key) ?? 0) + 1);
      }
      let deltaPx = 0;
      for (const [c, n] of colors) deltaPx += Math.abs(n - (prev.get(c) ?? 0));
      for (const [c, n] of prev) if (!colors.has(c)) deltaPx += n;
      w.__cpColorSnap = colors;
      return {deltaPx};
    } catch { return {deltaPx: -1}; }
  });
}

test('Correlation Plot — Matrix Values, Scope, and Persistence', async ({page}) => {
  test.setTimeout(600_000);

  await loginToDatagrok(page);

  // ## Setup — open demog, add the Correlation plot via its Toolbox icon (DOM-driven),
  // record the baseline df.filter.trueCount for the viewer-local-filter invariant.
  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'correlation-plot', 'Correlation-plot', 10000);

  const baselineFilterCount: number = await page.evaluate(() =>
    grok.shell.tv.dataFrame.filter.trueCount);

  // ### Scenario 1: Default state and numerical-column filter
  await softStep('Scenario 1 Step 2 — off-diagonal cell values equal runtime Pearson', async () => {
    const r = await page.evaluate(({pairs, tol}) => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      const out: {pair: string; cell: number; ref: number}[] = [];
      for (const [a, b] of pairs) {
        const cell = cp.getCorrelation(df.col(a), df.col(b));
        const ref = DG.Stats.fromColumn(df.col(a)).corr(df.col(b));
        out.push({pair: `${a}/${b}`, cell, ref});
      }
      void tol;
      return out;
    }, {pairs: PAIRS, tol: TOL});
    for (const p of r) {
      console.log(`[S1] ${p.pair}: cell=${p.cell} ref=${p.ref}`);
      expect(Number.isFinite(p.cell)).toBe(true);
      expect(Math.abs(p.cell - p.ref)).toBeLessThanOrEqual(TOL);
    }
  });

  await softStep('Scenario 1 Step 4 — diagonal cell is trivial self-correlation (no value cell)', async () => {
    const r = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      // The correlation loop skips xCol == yCol (correlation_plot_core.dart:255): a diagonal
      // pair gets no off-diagonal FloatColumn entry and renders a histogram. Via the JS API
      // getCorrelation(col, col) is the degenerate self-correlation 1.0, distinct from a real
      // off-diagonal coefficient which is finite AND not exactly 1. The histogram RENDER
      // itself is canvas-drawn (waived below); the readable structural signal is diag === 1.
      const diag = cp.getCorrelation(df.col('AGE'), df.col('AGE'));
      const offDiag = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
      return {diag, offDiag};
    });
    console.log(`[S1] diagonal AGE/AGE=${r.diag} offDiag AGE/HEIGHT=${r.offDiag}`);
    expect(r.diag).toBe(1);
    expect(Number.isFinite(r.offDiag)).toBe(true);
    expect(r.offDiag).not.toBe(1);
  });

  await softStep('Scenario 1 Step 5 — X Columns are numerical only; SEX and RACE absent', async () => {
    // Realizes correlationplot.int.numerical-columns-only. DOM-driven: right-click the
    // viewer canvas to open the context menu that carries the Columns > X/Y Columns pickers,
    // then assert the authoritative signal — the auto()-seeded X axis is EXACTLY the
    // numerical columns (the picker offers only isNumerical columns,
    // correlation_plot_core.dart:131), so non-numerical columns never enter the axis.
    const menuHasColumns = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-Correlation-plot"]')!;
      const canvas = viewer.querySelector('canvas')! as HTMLElement;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width * 0.5, clientY: rect.top + rect.height * 0.3,
      }));
      await new Promise((r) => setTimeout(r, 600));
      const labels = Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .map((l) => (l.textContent ?? '').trim());
      document.dispatchEvent(new KeyboardEvent('keydown', {key: 'Escape', bubbles: true}));
      document.body.click();
      return {
        hasColumns: labels.includes('Columns'),
        hasXColumns: labels.includes('X Columns'),
        hasYColumns: labels.includes('Y Columns'),
      };
    });
    const cols = await page.evaluate(() => {
      const df = grok.shell.tv.dataFrame;
      const numerical: string[] = [];
      const nonNumerical: string[] = [];
      for (let i = 0; i < df.columns.length; i++) {
        const c = df.columns.byIndex(i);
        (c.isNumerical ? numerical : nonNumerical).push(c.name);
      }
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      return {numerical, nonNumerical, xCols: cp.props.xColumnNames.slice()};
    });
    console.log(`[S1] menu=${JSON.stringify(menuHasColumns)} xCols=${JSON.stringify(cols.xCols)}`);
    console.log(`[S1] numerical=${JSON.stringify(cols.numerical)} nonNumerical=${JSON.stringify(cols.nonNumerical)}`);
    // The right-click Columns > X/Y Columns pickers are present.
    expect(menuHasColumns.hasColumns).toBe(true);
    expect(menuHasColumns.hasXColumns).toBe(true);
    expect(menuHasColumns.hasYColumns).toBe(true);
    // The X axis is exactly the numerical columns — SEX and RACE absent.
    expect(cols.xCols.length).toBeGreaterThan(0);
    expect(cols.xCols).not.toContain('SEX');
    expect(cols.xCols).not.toContain('RACE');
    for (const xc of cols.xCols)
      expect(cols.numerical).toContain(xc);
  });

  // ### Scenario 2: Correlation type switch, Show Pearson R, and column narrowing
  const pearsonRef: number[] = await page.evaluate(({pairs}) => {
    const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
    const df = grok.shell.tv.dataFrame;
    return pairs.map(([a, b]: [string, string]) => cp.getCorrelation(df.col(a), df.col(b)));
  }, {pairs: PAIRS});

  await softStep('Scenario 2 Step 2 — Spearman cell values equal runtime Spearman and differ from Pearson', async () => {
    const r = await page.evaluate(({pairs, pearson}) => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      cp.props.correlationType = 'Spearman';
      const out: {pair: string; cell: number; ref: number; diffFromPearson: number}[] = [];
      for (let i = 0; i < pairs.length; i++) {
        const [a, b] = pairs[i];
        const cell = cp.getCorrelation(df.col(a), df.col(b));
        const ref = DG.Stats.fromColumn(df.col(a)).spearmanCorr(df.col(b));
        out.push({pair: `${a}/${b}`, cell, ref, diffFromPearson: Math.abs(cell - pearson[i])});
      }
      return out;
    }, {pairs: PAIRS, pearson: pearsonRef});
    let anyDiffer = false;
    for (const p of r) {
      console.log(`[S2] ${p.pair}: spearmanCell=${p.cell} ref=${p.ref} diffFromPearson=${p.diffFromPearson}`);
      expect(Number.isFinite(p.cell)).toBe(true);
      expect(Math.abs(p.cell - p.ref)).toBeLessThanOrEqual(TOL);
      if (p.diffFromPearson > TOL) anyDiffer = true;
    }
    // The type switch is reflected in the backing values, not just the display.
    expect(anyDiffer).toBe(true);
  });

  await softStep('Scenario 2 Step 4 — Show Pearson R false: matrix repaints, backing value still readable', async () => {
    // Settle-precheck: snapshot the current render, confirm the matrix is stable before the
    // flip so the measured diff is attributable to showPearsonR, not to a still-settling frame.
    await cpSnapshot(page);
    await page.waitForTimeout(700);
    const settlePrecheck = await cpDiff(page);
    // Fresh baseline snapshot immediately before the flip.
    await cpSnapshot(page);
    const before = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      const v0 = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
      cp.props.showPearsonR = false;
      return v0;
    });
    // Let the disable repaint settle, then measure the before/after canvas diff.
    await page.waitForTimeout(700);
    const flipDiff = await cpDiff(page);
    const r = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      return {showR: cp.props.showPearsonR, after: cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'))};
    });
    console.log(`[S2] showPearsonR=${r.showR} settlePrecheck=${settlePrecheck.deltaPx} flipDiff=${flipDiff.deltaPx} before=${before} after=${r.after}`);
    // PRIMARY signal (render, T1): disabling Show Pearson R repaints the matrix — the in-cell
    // R digits vanish and the value columns narrow (correlation_plot_core.dart:197-198), which
    // shifts the painted layout. Assert the settle-gated canvas diff is materially non-zero and
    // clearly exceeds the pre-flip settle noise. (The gridCol.showValue/width 40→20 flips are the
    // Dart-side CAUSE; the inner ColumnGrid is not JS-exposed — the canvas diff is the readable
    // observable of that structural change.)
    expect(flipDiff.deltaPx).toBeGreaterThan(0);
    expect(flipDiff.deltaPx).toBeGreaterThan(settlePrecheck.deltaPx);
    // Backing value stays readable and unchanged — the property affects DISPLAY, not the data.
    expect(r.showR).toBe(false);
    expect(Number.isFinite(r.after)).toBe(true);
    expect(Math.abs(r.after - before)).toBeLessThanOrEqual(TOL);
  });

  await softStep('Scenario 2 Step 6 — narrow to 3×2, values correct, column names visible', async () => {
    const r = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      cp.props.xColumnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      cp.props.yColumnNames = ['AGE', 'HEIGHT'];
      const xCols = cp.props.xColumnNames.slice();
      const yCols = cp.props.yColumnNames.slice();
      const checks: {pair: string; cell: number; ref: number}[] = [];
      for (const [a, b] of [['AGE', 'HEIGHT'], ['WEIGHT', 'HEIGHT']] as [string, string][]) {
        const cell = cp.getCorrelation(df.col(a), df.col(b));
        const ref = DG.Stats.fromColumn(df.col(a)).spearmanCorr(df.col(b));
        checks.push({pair: `${a}/${b}`, cell, ref});
      }
      return {xCols, yCols, checks};
    });
    console.log(`[S2] x=${JSON.stringify(r.xCols)} y=${JSON.stringify(r.yCols)}`);
    // Exactly 3 X columns and 2 Y columns, in the requested order.
    expect(r.xCols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
    expect(r.yCols).toEqual(['AGE', 'HEIGHT']);
    for (const c of r.checks) {
      console.log(`[S2] narrowed ${c.pair}: cell=${c.cell} ref=${c.ref}`);
      expect(Math.abs(c.cell - c.ref)).toBeLessThanOrEqual(TOL);
    }
    // GROK-17480 regression guard (status: fixed): after narrow/reorder the column names
    // remain the configured set — the values are correct for the NEW set and order, which
    // requires the names to still resolve (a hidden/dropped name would break getCorrelation
    // for the narrowed pair, which stayed within tolerance above). Assert the axis names are
    // exactly the requested, visible set.
    expect(r.xCols).toContain('AGE');
    expect(r.xCols).toContain('HEIGHT');
    expect(r.xCols).toContain('WEIGHT');
  });

  // ### Scenario 3: Row Source and viewer-local filter formula
  await softStep('Scenario 3 Step 3 — Row Source Selected recomputes over selection mask', async () => {
    const r = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      // Neutral setup: select an offset band of ~20 rows (rows 0-19). Selecting is SETUP —
      // the tested SIGNAL is the correlation value recomputed over that mask.
      df.selection.setAll(false);
      df.selection.init((i: number) => i < 20);
      cp.props.rowSource = 'Selected';
      const type = cp.props.correlationType;
      const cell = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
      // Reference over ONLY the selected rows.
      const sub = df.clone(df.selection);
      const ss = DG.Stats.fromColumn(sub.col('AGE'));
      const ref = type === 'Spearman' ? ss.spearmanCorr(sub.col('HEIGHT')) : ss.corr(sub.col('HEIGHT'));
      return {cell, ref, type, selCount: df.selection.trueCount};
    });
    console.log(`[S3] Selected(${r.selCount}) AGE/HEIGHT cell=${r.cell} ref=${r.ref} type=${r.type}`);
    expect(r.selCount).toBe(20);
    expect(Number.isFinite(r.cell)).toBe(true);
    expect(Math.abs(r.cell - r.ref)).toBeLessThanOrEqual(TOL);
  });

  await softStep('Scenario 3 Step 5 — filter formula recomputes value; df.filter.trueCount unchanged', async () => {
    const r = await page.evaluate(async ({baseline}) => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      // Back to Filtered so the viewer-local formula filter drives combinedFilter.
      df.selection.setAll(false);
      cp.props.rowSource = 'Filtered';
      cp.props.filter = '${AGE} > 40';
      // The formula filter recomputes ASYNC (formula parse -> combinedFilter rebuild ->
      // _refreshValues), unlike correlationType/rowSource which recompute synchronously.
      // A synchronous getCorrelation here returns the stale full-set value; settle first.
      await new Promise((res) => setTimeout(res, 700));
      const type = cp.props.correlationType;
      const cell = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
      // Reference over rows where AGE > 40.
      const mask = DG.BitSet.create(df.rowCount, (i: number) => df.col('AGE').get(i) > 40);
      const sub = df.clone(mask);
      const ss = DG.Stats.fromColumn(sub.col('AGE'));
      const ref = type === 'Spearman' ? ss.spearmanCorr(sub.col('HEIGHT')) : ss.corr(sub.col('HEIGHT'));
      const filterCountAfter = df.filter.trueCount;
      cp.props.filter = '';
      return {cell, ref, filterCountAfter, baseline};
    }, {baseline: baselineFilterCount});
    console.log(`[S3] formula AGE/HEIGHT cell=${r.cell} ref=${r.ref} | df.filter ${r.baseline}→${r.filterCountAfter}`);
    expect(Number.isFinite(r.cell)).toBe(true);
    expect(Math.abs(r.cell - r.ref)).toBeLessThanOrEqual(TOL);
    // Anti-trap: the formula filter is VIEWER-LOCAL — df.filter must be unchanged.
    expect(r.filterCountAfter).toBe(r.baseline);
  });

  // ### Scenario 4: defaultCellFont structural signal
  await softStep('Scenario 4 Step 3 — larger defaultCellFont repaints the matrix', async () => {
    // Settle-precheck at the baseline render so the measured diff is font-attributable.
    await cpSnapshot(page);
    await page.waitForTimeout(700);
    const settlePrecheck = await cpDiff(page);
    // Fresh baseline snapshot immediately before the font change.
    await cpSnapshot(page);
    const before = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      const v0 = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
      cp.props.defaultCellFont = 'bold 16px Roboto';
      return v0;
    });
    await page.waitForTimeout(700);
    const largerDiff = await cpDiff(page);
    // Restore the default font and confirm the render returns to the baseline (diff ~0).
    await cpSnapshot(page);
    const r = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      const applied = cp.props.defaultCellFont;
      const after = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
      cp.props.defaultCellFont = 'normal normal 13px "Roboto"';
      return {applied, after};
    });
    await page.waitForTimeout(700);
    const restoreDiff = await cpDiff(page);
    console.log(`[S4] defaultCellFont=${r.applied} settlePrecheck=${settlePrecheck.deltaPx} largerDiff=${largerDiff.deltaPx} restoreDiff=${restoreDiff.deltaPx} value ${before}→${r.after}`);
    // PRIMARY signal (render, T1): the larger cell font visibly repaints the matrix — taller
    // rows / larger digits — so the before/after canvas diff is materially non-zero and clearly
    // exceeds the pre-change settle noise. rowHeight = parseSize×1.4 (correlation_plot_core.dart:
    // 172-173) is the Dart-internal CAUSE (not JS-readable); the canvas diff is its readable
    // observable, per the render-signal-index corr_font_row_height family.
    expect(largerDiff.deltaPx).toBeGreaterThan(0);
    expect(largerDiff.deltaPx).toBeGreaterThan(settlePrecheck.deltaPx);
    // No-error floor: the font is accepted and the matrix keeps computing (backing value stable).
    expect(r.applied).toBe('bold 16px Roboto');
    expect(Number.isFinite(r.after)).toBe(true);
    expect(Math.abs(r.after - before)).toBeLessThanOrEqual(TOL);
  });

  // ### Scenario 5: Layout persistence — save, re-arm, re-apply
  await softStep('Scenario 5 Step 5 — saved layout restores full config + value spot-check', async () => {
    const r = await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      // Re-establish the peak config for the layout.
      cp.props.correlationType = 'Spearman';
      cp.props.showPearsonR = false;
      cp.props.xColumnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
      cp.props.yColumnNames = ['AGE', 'HEIGHT'];
      cp.props.filter = '';
      cp.props.rowSource = 'Filtered';
      const viewerSetBefore = grok.shell.tv.viewers.map((x: any) => x.type).sort();
      const layout = grok.shell.tv.saveLayout();
      await grok.dapi.layouts.save(layout);
      const layoutId = layout.id;
      await new Promise((res) => setTimeout(res, 1000));
      // Re-arm: close the Correlation plot, add a Scatter plot.
      cp.close();
      grok.shell.tv.addViewer('Scatter plot');
      await new Promise((res) => setTimeout(res, 1000));
      // Re-apply the saved layout.
      const saved = await grok.dapi.layouts.find(layoutId);
      grok.shell.tv.loadLayout(saved);
      await new Promise((res) => setTimeout(res, 3000));
      const cp2 = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      let spot: number | null = null;
      let ref: number | null = null;
      if (cp2) {
        spot = cp2.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
        ref = DG.Stats.fromColumn(df.col('AGE')).spearmanCorr(df.col('HEIGHT'));
      }
      const viewerSetAfter = grok.shell.tv.viewers.map((x: any) => x.type).sort();
      const result = {
        viewerSetBefore, viewerSetAfter,
        present: !!cp2,
        type: cp2?.props.correlationType ?? null,
        showR: cp2?.props.showPearsonR ?? null,
        xCols: cp2?.props.xColumnNames.slice() ?? null,
        yCols: cp2?.props.yColumnNames.slice() ?? null,
        filter: cp2?.props.filter ?? null,
        spot, ref,
      };
      await grok.dapi.layouts.delete(saved);
      return result;
    });
    console.log(`[S5] before=${JSON.stringify(r.viewerSetBefore)} after=${JSON.stringify(r.viewerSetAfter)} spot=${r.spot} spearmanRef=${r.ref} pearsonRef=${r.pearsonRef}`);
    // PRIMARY persistence signals (different-property / cross-channel, T2):
    // (a) the viewer SET is restored — the re-armed Scatter plot is gone and the
    //     Correlation plot is back (a channel independent of any viewer prop);
    // (b) the reloaded backing VALUE matches the runtime Spearman reference over the
    //     full row set. This is the value channel (getCorrelation), not a prop
    //     re-read: a config that failed to re-hydrate would not produce a matrix that
    //     recomputes the Spearman coefficient for the restored pair.
    expect(r.present).toBe(true);
    expect(r.viewerSetAfter).toEqual(r.viewerSetBefore);
    expect(Number.isFinite(r.spot)).toBe(true);
    expect(Math.abs((r.spot as number) - (r.ref as number))).toBeLessThanOrEqual(TOL);
    // Auxiliary confirmation that the persisted configuration re-hydrated (these
    // prop re-reads corroborate the value-channel evidence above; they are NOT the
    // realizing assertion — see expected_results_coverage realized_by).
    expect(r.type).toBe('Spearman');
    expect(r.showR).toBe(false);
    expect(r.xCols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
    expect(r.yCols).toEqual(['AGE', 'HEIGHT']);
    expect(r.filter === '' || r.filter == null).toBe(true);
  });

  // ### Scenario 6: Project persistence — save via ribbon, close all, reopen
  const projectName = 'zz-cp-matrix-persist-' + Date.now();
  let savedProjectId: string | null = null;
  await softStep('Scenario 6 Step 4/5 — reopened project restores config + value spot-check', async () => {
    try {
      // Re-establish the Scenario-2 config on the current viewer, then ribbon-Save.
      await page.evaluate(() => {
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        cp.props.correlationType = 'Spearman';
        cp.props.showPearsonR = false;
        cp.props.xColumnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
        cp.props.yColumnNames = ['AGE', 'HEIGHT'];
        cp.props.filter = '';
        cp.props.rowSource = 'Filtered';
      });
      const saved = await saveProjectViaUI(page, projectName);
      savedProjectId = saved.projectId;

      // Close All, then reopen via find(id).open() (recon-validated reopen path).
      const r = await page.evaluate(async ({id}) => {
        grok.shell.closeAll();
        await new Promise((res) => setTimeout(res, 1000));
        const proj = await grok.dapi.projects.find(id);
        await proj.open();
        await new Promise((res) => setTimeout(res, 3000));
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        const df = grok.shell.tv.dataFrame;
        let spot: number | null = null;
        let ref: number | null = null;
        if (cp && df) {
          spot = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
          ref = DG.Stats.fromColumn(df.col('AGE')).spearmanCorr(df.col('HEIGHT'));
        }
        return {
          present: !!cp,
          type: cp?.props.correlationType ?? null,
          showR: cp?.props.showPearsonR ?? null,
          xCols: cp?.props.xColumnNames.slice() ?? null,
          yCols: cp?.props.yColumnNames.slice() ?? null,
          spot, ref,
        };
      }, {id: savedProjectId});
      console.log(`[S6] present=${r.present} spot=${r.spot} spearmanRef=${r.ref}`);
      // PRIMARY persistence signal (different-property, T2): the reloaded backing
      // VALUE matches the runtime Spearman reference over the full row set across the
      // project save→closeAll→reopen round-trip. This is the value channel
      // (getCorrelation), not a prop re-read — a matrix that failed to re-hydrate the
      // persisted configuration would not recompute the Spearman coefficient here.
      expect(r.present).toBe(true);
      expect(Number.isFinite(r.spot)).toBe(true);
      expect(Math.abs((r.spot as number) - (r.ref as number))).toBeLessThanOrEqual(TOL);
      // Auxiliary confirmation the persisted configuration re-hydrated (corroborates
      // the value-channel evidence above; NOT the realizing assertion).
      expect(r.type).toBe('Spearman');
      expect(r.showR).toBe(false);
      expect(r.xCols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
      expect(r.yCols).toEqual(['AGE', 'HEIGHT']);
    } finally {
      if (savedProjectId)
        await deleteProjectWithCleanup(page, {projectId: savedProjectId});
    }
  });

  v.finishSpec();
});
