/* ---
realizes: [correlationplot.cp.matrix-values-scope-persist, correlationplot.int.numerical-columns-only]
--- */
import {test, expect} from '@playwright/test';
import {loginToDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';
import {saveProjectViaApi, deleteProjectWithCleanup} from '../../helpers/projects';

declare const grok: any;
declare const DG: any;

test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const TOL = 1e-3;
const PAIRS: [string, string][] = [['AGE', 'HEIGHT'], ['AGE', 'WEIGHT'], ['HEIGHT', 'WEIGHT']];

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

  await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
  await v.addViewerByIcon(page, 'correlation-plot', 'Correlation-plot', 10000);

  const baselineFilterCount: number = await page.evaluate(() =>
    grok.shell.tv.dataFrame.filter.trueCount);

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

    const menuHasColumns = await page.evaluate(async () => {
      const viewer = document.querySelector('[name="viewer-Correlation-plot"]')!;
      const canvas = viewer.querySelector('canvas')! as HTMLElement;
      const rect = canvas.getBoundingClientRect();
      canvas.dispatchEvent(new MouseEvent('contextmenu', {
        bubbles: true, cancelable: true, button: 2,
        clientX: rect.left + rect.width * 0.5, clientY: rect.top + rect.height * 0.3,
      }));

      const readLabels = () => Array.from(document.querySelectorAll('.d4-menu-item-label'))
        .map((l) => (l.textContent ?? '').trim());
      const menuDeadline = Date.now() + 600;
      let labels = readLabels();
      while (!(labels.includes('Columns') && labels.includes('X Columns') &&
               labels.includes('Y Columns')) && Date.now() < menuDeadline) {
        await new Promise((r) => setTimeout(r, 25));
        labels = readLabels();
      }
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

    expect(menuHasColumns.hasColumns).toBe(true);
    expect(menuHasColumns.hasXColumns).toBe(true);
    expect(menuHasColumns.hasYColumns).toBe(true);

    expect(cols.xCols.length).toBeGreaterThan(0);
    expect(cols.xCols).not.toContain('SEX');
    expect(cols.xCols).not.toContain('RACE');
    for (const xc of cols.xCols)
      expect(cols.numerical).toContain(xc);
  });

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

    expect(anyDiffer).toBe(true);
  });

  await softStep('Scenario 2 Step 4 — Show Pearson R false: matrix repaints, backing value still readable', async () => {

    await cpSnapshot(page);

    await page.waitForTimeout(700);
    const settlePrecheck = await cpDiff(page);

    await cpSnapshot(page);
    const before = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      const v0 = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
      cp.props.showPearsonR = false;
      return v0;
    });

    await v.waitForViewerRendered(page, 'Correlation plot', 700);
    const flipDiff = await cpDiff(page);
    const r = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      return {showR: cp.props.showPearsonR, after: cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'))};
    });
    console.log(`[S2] showPearsonR=${r.showR} settlePrecheck=${settlePrecheck.deltaPx} flipDiff=${flipDiff.deltaPx} before=${before} after=${r.after}`);

    expect(flipDiff.deltaPx).toBeGreaterThan(0);
    expect(flipDiff.deltaPx).toBeGreaterThan(settlePrecheck.deltaPx);

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

    expect(r.xCols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
    expect(r.yCols).toEqual(['AGE', 'HEIGHT']);
    for (const c of r.checks) {
      console.log(`[S2] narrowed ${c.pair}: cell=${c.cell} ref=${c.ref}`);
      expect(Math.abs(c.cell - c.ref)).toBeLessThanOrEqual(TOL);
    }

    expect(r.xCols).toContain('AGE');
    expect(r.xCols).toContain('HEIGHT');
    expect(r.xCols).toContain('WEIGHT');
  });

  await softStep('Scenario 3 Step 3 — Row Source Selected recomputes over selection mask', async () => {
    const r = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;

      df.selection.setAll(false);
      df.selection.init((i: number) => i < 20);
      cp.props.rowSource = 'Selected';
      const type = cp.props.correlationType;
      const cell = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));

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

      df.selection.setAll(false);
      cp.props.rowSource = 'Filtered';
      cp.props.filter = '${AGE} > 40';

      await new Promise((res) => setTimeout(res, 700));
      const type = cp.props.correlationType;
      const cell = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));

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

    expect(r.filterCountAfter).toBe(r.baseline);
  });

  await softStep('Scenario 4 Step 3 — larger defaultCellFont repaints the matrix', async () => {

    await cpSnapshot(page);

    await page.waitForTimeout(700);
    const settlePrecheck = await cpDiff(page);

    await cpSnapshot(page);
    const before = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      const v0 = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
      cp.props.defaultCellFont = 'bold 16px Roboto';
      return v0;
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 700);
    const largerDiff = await cpDiff(page);

    await cpSnapshot(page);
    const r = await page.evaluate(() => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
      const df = grok.shell.tv.dataFrame;
      const applied = cp.props.defaultCellFont;
      const after = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
      cp.props.defaultCellFont = 'normal normal 13px "Roboto"';
      return {applied, after};
    });
    await v.waitForViewerRendered(page, 'Correlation plot', 700);
    const restoreDiff = await cpDiff(page);
    console.log(`[S4] defaultCellFont=${r.applied} settlePrecheck=${settlePrecheck.deltaPx} largerDiff=${largerDiff.deltaPx} restoreDiff=${restoreDiff.deltaPx} value ${before}→${r.after}`);

    expect(largerDiff.deltaPx).toBeGreaterThan(0);
    expect(largerDiff.deltaPx).toBeGreaterThan(settlePrecheck.deltaPx);

    expect(r.applied).toBe('bold 16px Roboto');
    expect(Number.isFinite(r.after)).toBe(true);
    expect(Math.abs(r.after - before)).toBeLessThanOrEqual(TOL);
  });

  await softStep('Scenario 5 Step 5 — saved layout restores full config + value spot-check', async () => {
    const r = await page.evaluate(async () => {
      const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');

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

      cp.close();
      grok.shell.tv.addViewer('Scatter plot');
      const rearmDeadline = Date.now() + 1000;
      while (Date.now() < rearmDeadline) {
        const types = grok.shell.tv.viewers.map((x: any) => x.type);
        if (types.includes('Scatter plot') && !types.includes('Correlation plot')) break;
        await new Promise((res) => setTimeout(res, 50));
      }

      const saved = await grok.dapi.layouts.find(layoutId);
      const applyDeadline = Date.now() + 3000;
      grok.shell.tv.loadLayout(saved);
      while (Date.now() < applyDeadline) {
        const back = grok.shell.tv?.viewers?.find((x: any) => x.type === 'Correlation plot');
        if (back && back.props) break;
        await new Promise((res) => setTimeout(res, 50));
      }
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

    expect(r.present).toBe(true);
    expect(r.viewerSetAfter).toEqual(r.viewerSetBefore);
    expect(Number.isFinite(r.spot)).toBe(true);
    expect(Math.abs((r.spot as number) - (r.ref as number))).toBeLessThanOrEqual(TOL);

    expect(r.type).toBe('Spearman');
    expect(r.showR).toBe(false);
    expect(r.xCols).toEqual(['AGE', 'HEIGHT', 'WEIGHT']);
    expect(r.yCols).toEqual(['AGE', 'HEIGHT']);
    expect(r.filter === '' || r.filter == null).toBe(true);
  });

  const projectName = 'zz-cp-matrix-persist-' + Date.now();
  let savedProjectId: string | null = null;
  await softStep('Scenario 6 Step 4/5 — reopened project restores config + value spot-check', async () => {
    try {

      await page.evaluate(() => {
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        cp.props.correlationType = 'Spearman';
        cp.props.showPearsonR = false;
        cp.props.xColumnNames = ['AGE', 'HEIGHT', 'WEIGHT'];
        cp.props.yColumnNames = ['AGE', 'HEIGHT'];
        cp.props.filter = '';
        cp.props.rowSource = 'Filtered';
      });
      const saved = await saveProjectViaApi(page, projectName);
      savedProjectId = saved.projectId;

      await v.closeAllAndWait(page);
      const r = await page.evaluate(async ({id}) => {
        const proj = await grok.dapi.projects.find(id);
        await proj.open();

        const openDeadline = Date.now() + 3000;
        while (Date.now() < openDeadline) {
          const back = grok.shell.tv?.viewers?.find((x: any) => x.type === 'Correlation plot');
          if (back && back.props && grok.shell.tv?.dataFrame?.rowCount > 0) break;
          await new Promise((res) => setTimeout(res, 50));
        }
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

      expect(r.present).toBe(true);
      expect(Number.isFinite(r.spot)).toBe(true);
      expect(Math.abs((r.spot as number) - (r.ref as number))).toBeLessThanOrEqual(TOL);

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
