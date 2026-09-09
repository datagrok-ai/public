/* ---
realizes: [correlationplot.cp.matrix-values-scope-persist, correlationplot.int.numerical-columns-only]
--- */
import {expect, Page} from '@playwright/test';
import {localTest as test} from '../../shared-page';
import {openDatagrok, specTestOptions, softStep} from '../../spec-login';
import * as v from '../../helpers/viewers';

declare const grok: any;
declare const DG: any;

// The client-side half of the matrix-values scenario. The layout and project round-trips
// (Scenarios 5 and 6) live in correlation-plot-matrix-values-server-spec.ts.
test.use(specTestOptions);

const datasetPath = 'System:DemoFiles/demog.csv';
const TOL = 1e-3;
const PAIRS: [string, string][] = [['AGE', 'HEIGHT'], ['AGE', 'WEIGHT'], ['HEIGHT', 'WEIGHT']];
const CP = 'Correlation plot';
const CANVAS = 'canvas[name="canvas"]';

// The idle noise floor of a canvas that has stopped painting; the assertion that follows demands
// a repaint larger than it.
async function settledSnap(page: Page): Promise<number> {
  await v.waitForCanvasQuiet(page, CP, {canvasSelector: CANVAS, timeoutMs: 1500, optional: true});
  await v.snapshotCanvasColors(page, CP, CANVAS);
  return (await v.diffCanvasColors(page, CP, CANVAS)).deltaPx;
}

async function repaintOver(page: Page, floor: number, capMs: number): Promise<number> {
  return v.waitForCanvasChange(page, CP, {canvasSelector: CANVAS, minDelta: floor + 1, timeoutMs: capMs})
    .catch(() => -1);
}

test('Correlation Plot — Matrix Values, Scope, and Persistence', async ({page}) => {
  test.setTimeout(600_000);

  await openDatagrok(page);

  try {
    await v.openTable(page, {path: datasetPath, semTypeTimeoutMs: 3000});
    await v.addViewerByIcon(page, 'correlation-plot', 'Correlation-plot', 10000);

    const baselineFilterCount: number = await page.evaluate(() =>
      grok.shell.tv.dataFrame.filter.trueCount);

    await softStep('Scenario 1 Step 2 — off-diagonal cell values equal runtime Pearson', async () => {
      const r = await page.evaluate(({pairs}) => {
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        const df = grok.shell.tv.dataFrame;
        const out: {pair: string; cell: number; ref: number}[] = [];
        for (const [a, b] of pairs) {
          const cell = cp.getCorrelation(df.col(a), df.col(b));
          const ref = DG.Stats.fromColumn(df.col(a)).corr(df.col(b));
          out.push({pair: `${a}/${b}`, cell, ref});
        }
        return out;
      }, {pairs: PAIRS});
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
        const w = window as any;
        const viewer = document.querySelector('[name="viewer-Correlation-plot"]')!;
        const canvas = viewer.querySelector('canvas')! as HTMLElement;
        const rect = canvas.getBoundingClientRect();
        canvas.dispatchEvent(new MouseEvent('contextmenu', {
          bubbles: true, cancelable: true, button: 2,
          clientX: rect.left + rect.width * 0.5, clientY: rect.top + rect.height * 0.3,
        }));

        const readLabels = () => Array.from(document.querySelectorAll('.d4-menu-item-label'))
          .map((l) => (l.textContent ?? '').trim());
        const labels: string[] = await w.__poll(readLabels,
          (ls: string[]) => ls.includes('Columns') && ls.includes('X Columns') && ls.includes('Y Columns'), 600, 25);
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
      const settlePrecheck = await settledSnap(page);
      const before = await page.evaluate(() => {
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        const df = grok.shell.tv.dataFrame;
        const v0 = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
        cp.props.showPearsonR = false;
        return v0;
      });
      const flipDiff = await repaintOver(page, settlePrecheck, 2000);
      const r = await page.evaluate(() => {
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        const df = grok.shell.tv.dataFrame;
        return {showR: cp.props.showPearsonR, after: cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'))};
      });
      console.log(`[S2] showPearsonR=${r.showR} settlePrecheck=${settlePrecheck} flipDiff=${flipDiff} before=${before} after=${r.after}`);

      expect(flipDiff).toBeGreaterThan(0);
      expect(flipDiff).toBeGreaterThan(settlePrecheck);

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
        const w = window as any;
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        const df = grok.shell.tv.dataFrame;

        df.selection.setAll(false);
        cp.props.rowSource = 'Filtered';
        const unfiltered = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
        cp.props.filter = '${AGE} > 40';
        const cell = await w.__poll(() => cp.getCorrelation(df.col('AGE'), df.col('HEIGHT')),
          (c: number) => c !== unfiltered, 700, 25);
        const type = cp.props.correlationType;

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
      const settlePrecheck = await settledSnap(page);
      const before = await page.evaluate(() => {
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        const df = grok.shell.tv.dataFrame;
        const v0 = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
        cp.props.defaultCellFont = 'bold 16px Roboto';
        return v0;
      });
      const largerDiff = await repaintOver(page, settlePrecheck, 2000);

      await settledSnap(page);
      const r = await page.evaluate(() => {
        const cp = grok.shell.tv.viewers.find((x: any) => x.type === 'Correlation plot');
        const df = grok.shell.tv.dataFrame;
        const applied = cp.props.defaultCellFont;
        const after = cp.getCorrelation(df.col('AGE'), df.col('HEIGHT'));
        cp.props.defaultCellFont = 'normal normal 13px "Roboto"';
        return {applied, after};
      });
      const restoreDiff = await repaintOver(page, 0, 2000);
      console.log(`[S4] defaultCellFont=${r.applied} settlePrecheck=${settlePrecheck} largerDiff=${largerDiff} restoreDiff=${restoreDiff} value ${before}→${r.after}`);

      expect(largerDiff).toBeGreaterThan(0);
      expect(largerDiff).toBeGreaterThan(settlePrecheck);

      expect(r.applied).toBe('bold 16px Roboto');
      expect(Number.isFinite(r.after)).toBe(true);
      expect(Math.abs(r.after - before)).toBeLessThanOrEqual(TOL);
    });
  } finally {
    await page.evaluate(() => { delete (window as any).__canvasColorSnap; }).catch(() => {});
    await v.closeAllAndWait(page);
  }

  v.finishSpec();
});
