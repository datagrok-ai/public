// Tests for control comparisons (Dunnett + Holm-Welch).
//
// Dunnett fixtures are the SciPy-validated cases from @datagrok-libraries/sci-comp
// (stats/__tests__/fixtures/dunnett.json, scipy 1.17 scipy.stats.dunnett) — the port must
// reproduce them after the gammaln/normalCdf swap to jStat. Holm fixtures are cross-checked
// against R `p.adjust(p, method = "holm")`. The Holm-Welch wiring is pinned to the package's
// own `twoSampleTTest` via a canary on the k=2 reduction.

import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {
  controlComparisons, holmCorrect, ControlComparison, ControlComparisonsMethod,
} from '../control-comparisons/control-comparisons-tools';
import {twoSampleTTest} from '../ttest/ttest-tools';
import {
  conclusionColumnPerRow, CONCLUSION_LABEL, CONCLUSION_COL_NAME,
} from '../group-comparison/conclusion-column';
import {FIXTURES, HolmWelchExpected} from './control-comparisons-fixtures';

const ALPHA = 0.05;
const CONTROL = 'ctrl';

/** Build a double-precision value column. */
function valueCol(name: string, data: readonly number[]): DG.Column<DG.COLUMN_TYPE.FLOAT> {
  return DG.Column.fromFloat64Array(name, Float64Array.from(data)) as unknown as DG.Column<DG.COLUMN_TYPE.FLOAT>;
}

/** Approximate-equality check with relative + absolute tolerances. */
function close(actual: number, expected: number, label: string, rtol = 1e-4, atol = 1e-6): void {
  const diff = Math.abs(actual - expected);
  const threshold = atol + rtol * Math.abs(expected);
  expect(diff <= threshold, true,
    `${label}: |${actual} - ${expected}| = ${diff} > ${threshold} (rtol=${rtol}, atol=${atol})`);
}

/**
 * Run control comparisons on a control array + a list of treated arrays.
 * Control is labelled `ctrl`; treated groups `t0, t1, ...`. Returns the per-treated comparison
 * keyed by its original treated index, plus the control mean/n.
 */
function runCC(control: readonly number[], treated: readonly (readonly number[])[],
  method: ControlComparisonsMethod, alpha = ALPHA) {
  const labels: string[] = [];
  const vals: number[] = [];
  for (const v of control) {labels.push(CONTROL); vals.push(v);}
  treated.forEach((arr, i) => {
    for (const v of arr) {labels.push(`t${i}`); vals.push(v);}
  });

  const cats = DG.Column.fromStrings('group', labels);
  const values = valueCol('feature', vals);
  const controlCode = cats.categories.indexOf(CONTROL);

  const report = controlComparisons(cats, values, controlCode, cats.categories.length, {method, alpha});

  // Map each treated index to its comparison via the category code of `t{i}`.
  const byTreated = treated.map((_, i) => {
    const code = cats.categories.indexOf(`t${i}`);
    return report.comparisons.find((c) => c.groupCode === code) as ControlComparison;
  });
  return {report, byTreated};
}

// ── SciPy-validated Dunnett fixtures (sci-comp dunnett.json) ──────

interface DunnettCase {
  name: string;
  control: number[];
  treated: number[][];
  expected: {statistic: number, pValueAdj: number}[];
}

const DUNNETT_CASES: DunnettCase[] = [
  {
    name: '[Dunnett 1955] blood counts',
    control: [7.4, 8.5, 7.2, 8.24, 9.84, 8.32],
    treated: [[9.76, 8.8, 7.68, 9.36], [12.8, 9.68, 12.16, 9.2, 10.55]],
    expected: [
      {statistic: 0.8570321259343803, pValueAdj: 0.6200769777606123},
      {statistic: 3.693752311235726, pValueAdj: 0.0058540018515012004},
    ],
  },
  {
    name: '3 treated, balanced',
    control: [1.0, 2.0, 3.0, 2.5, 1.5],
    treated: [[3.0, 4.0, 5.0, 4.5], [5.0, 6.0, 7.0, 6.5], [7.0, 8.0, 9.0, 8.5]],
    expected: [
      {statistic: 3.7940263593345884, pValueAdj: 0.006095222685964408},
      {statistic: 7.364874697531849, pValueAdj: 1.0626325918750013e-05},
      {statistic: 10.93572303572911, pValueAdj: 5.261925783983656e-09},
    ],
  },
  {
    name: '2 treated, unbalanced',
    control: [10.0, 11.0, 9.5, 10.5, 10.2, 11.3, 9.8],
    treated: [[12.0, 11.5, 13.0], [14.0, 15.0, 14.5, 13.5]],
    expected: [
      {statistic: 3.9784013687041444, pValueAdj: 0.004172949665478098},
      {statistic: 9.34453984503784, pValueAdj: 7.816759652889971e-07},
    ],
  },
];

category('Control comparisons: Dunnett (SciPy fixtures)', () => {
  for (const c of DUNNETT_CASES) {
    test(c.name, async () => {
      const {byTreated} = runCC(c.control, c.treated, 'Dunnett');
      c.expected.forEach((exp, i) => {
        // sci-comp fixture tolerances: statistic 0.01, p 0.001 (QMC vs Simpson integration floor).
        close(byTreated[i].statistic, exp.statistic, `${c.name}[${i}] statistic`, 0, 0.01);
        close(byTreated[i].pValueAdj, exp.pValueAdj, `${c.name}[${i}] pValueAdj`, 0, 1e-3);
      });
    }, {timeout: 30000});
  }
});

// ── Holm correction (R p.adjust(..., "holm")) ────────────────────

category('Control comparisons: Holm correction', () => {
  test('matches R p.adjust holm', async () => {
    // R: p.adjust(c(0.01, 0.04, 0.03, 0.005), method = "holm") -> 0.03 0.06 0.06 0.02
    const adj = holmCorrect([0.01, 0.04, 0.03, 0.005]);
    close(adj[0], 0.03, 'holm[0]', 0, 1e-12);
    close(adj[1], 0.06, 'holm[1]', 0, 1e-12);
    close(adj[2], 0.06, 'holm[2]', 0, 1e-12);
    close(adj[3], 0.02, 'holm[3]', 0, 1e-12);
  });

  test('clamps at 1 and stays monotone', async () => {
    // p.adjust(c(0.5, 0.6, 0.7), "holm") -> 1 1 1 (3*0.5=1.5 clamps; monotone after).
    const adj = holmCorrect([0.5, 0.6, 0.7]);
    expect(adj.every((p) => p === 1), true, `expected all 1, got ${adj}`);
  });

  test('single comparison is unchanged', async () => {
    const adj = holmCorrect([0.023]);
    close(adj[0], 0.023, 'holm single', 0, 1e-12);
  });
});

// ── Holm-Welch wiring pinned to the package t-test ───────────────

category('Control comparisons: Holm-Welch', () => {
  test('k=2 reduces to the two-sample Welch t-test', async () => {
    const control = [5.1, 4.8, 5.5, 5.0, 4.9, 5.2];
    const treated = [6.0, 6.3, 5.8, 6.1, 6.4, 5.9];

    // Control comparisons, one treated group.
    const {byTreated} = runCC(control, [treated], 'Holm-Welch');
    const cc = byTreated[0];

    // Direct two-sample Welch t-test on the same data (control = code 0, treated = code 1).
    const labels = [...control.map(() => 'ctrl'), ...treated.map(() => 'trt')];
    const cats = DG.Column.fromStrings('g', labels);
    const vals = valueCol('f', [...control, ...treated]);
    const tt = twoSampleTTest(cats, vals, cats.categories.length, {method: 'Welch', alpha: ALPHA});

    // meanDiff/t/df/CI/hedgesG must match; with one comparison Holm leaves the raw p untouched.
    close(cc.meanDiff, tt.meanDiff, 'meanDiff', 0, 1e-12);
    close(cc.statistic, tt.t, 't', 0, 1e-12);
    close(cc.df, tt.df, 'df', 0, 1e-12);
    close(cc.ciLow, tt.ciLow, 'ciLow', 0, 1e-12);
    close(cc.ciHigh, tt.ciHigh, 'ciHigh', 0, 1e-12);
    close(cc.hedgesG, tt.hedgesG, 'hedgesG', 0, 1e-12);
    close(cc.pValueRaw, tt.pValue, 'pValueRaw', 0, 1e-12);
    close(cc.pValueAdj, tt.pValue, 'pValueAdj (m=1)', 0, 1e-12);
  });

  test('meanDiff is signed group − control; Holm ≥ raw', async () => {
    const control = [1, 2, 3, 2, 1, 2];
    const treated = [[10, 11, 12, 11], [1.1, 2.1, 1.9, 2.0]];
    const {byTreated} = runCC(control, treated, 'Holm-Welch');

    expect(byTreated[0].meanDiff > 0, true, 'group above control -> positive meanDiff');
    for (const c of byTreated)
      expect(c.pValueAdj >= c.pValueRaw - 1e-12, true, 'adjusted p must be ≥ raw p');
  });
});

// ── Report shape / metadata ──────────────────────────────────────

category('Control comparisons: report', () => {
  test('Dunnett exposes pooled MSE/df; control metadata correct', async () => {
    const control = [1.0, 2.0, 3.0, 2.5, 1.5];
    const treated = [[3.0, 4.0, 5.0, 4.5], [5.0, 6.0, 7.0, 6.5]];
    const {report} = runCC(control, treated, 'Dunnett');

    expect(report.method, 'Dunnett');
    expect(report.controlN, 5);
    close(report.controlMean, 2.0, 'controlMean', 0, 1e-12);
    expect(report.comparisons.length, 2);
    // Pooled df = N − k = (5 + 4 + 4) − 3 = 10.
    expect(report.pooledDF, 10);
    expect(report.pooledMSE !== undefined && report.pooledMSE > 0, true, 'pooledMSE present');
  });

  test('Holm-Welch omits pooled metadata', async () => {
    const {report} = runCC([1, 2, 3, 2, 1], [[3, 4, 5, 4]], 'Holm-Welch');
    expect(report.method, 'Holm-Welch');
    expect(report.pooledDF === undefined, true, 'no pooled df for Holm-Welch');
    expect(report.pooledMSE === undefined, true, 'no pooled MSE for Holm-Welch');
  });
});

// ── Edge cases ───────────────────────────────────────────────────

category('Control comparisons: edge cases', () => {
  // k=2 (one treated group): pooled-variance Dunnett reduces to Student's two-sample t-test, and the
  // single-contrast max distribution is exactly Student-t — so the adjusted p equals the two-sided
  // Student p. Cross-checks the ported dunnettMaxCdf against the package's own (R-validated) t-test.
  test('Dunnett k=2 equals Student two-sample t-test', async () => {
    const control = [5.1, 4.8, 5.5, 5.0, 4.9, 5.2, 5.3];
    const treated = [6.0, 6.3, 5.8, 6.1, 6.4, 5.9];
    const {byTreated} = runCC(control, [treated], 'Dunnett');
    const cc = byTreated[0];

    const labels = [...control.map(() => 'ctrl'), ...treated.map(() => 'trt')];
    const cats = DG.Column.fromStrings('g', labels);
    const vals = valueCol('f', [...control, ...treated]);
    const student = twoSampleTTest(cats, vals, cats.categories.length, {method: 'Student', alpha: ALPHA});

    close(cc.statistic, student.t, 'Dunnett t == Student t', 0, 1e-9);
    close(cc.df, student.df, 'Dunnett df == pooled df', 0, 1e-12);
    // Single contrast: adjusted p == two-sided Student p (within the Simpson integration floor).
    close(cc.pValueAdj, student.pValue, 'Dunnett adj p == Student two-sided p', 0, 2e-3);
  }, {timeout: 30000});

  // Near-zero difference: group mean ≈ control mean → meanDiff ≈ 0, not significant, p ≈ 1.
  test('near-zero difference is not significant (both methods)', async () => {
    const control = [10.0, 10.5, 9.5, 10.2, 9.8, 10.1];
    const group = [10.1, 9.9, 10.3, 9.7, 10.2, 10.0]; // essentially the same distribution
    for (const method of ['Dunnett', 'Holm-Welch'] as ControlComparisonsMethod[]) {
      const {byTreated} = runCC(control, [group], method);
      const c = byTreated[0];
      expect(Math.abs(c.meanDiff) < 0.3, true, `${method}: meanDiff ≈ 0`);
      expect(c.significant, false, `${method}: not significant`);
      expect(c.pValueAdj > 0.3, true, `${method}: large p`);
    }
  }, {timeout: 30000});

  // Extreme variance ratio: Welch–Satterthwaite df is pushed toward its minimum. The per-pair stats
  // come from the R-validated twoSampleTTest kernel, so we pin Holm-Welch to it exactly (m=1 here).
  test('Holm-Welch extreme variance ratio matches Welch t-test (df → small)', async () => {
    const control = [50.0, 50.1, 49.9, 50.05, 49.95, 50.02]; // tiny variance
    const group = [40.0, 60.0, 30.0, 70.0, 45.0]; // huge variance
    const {byTreated} = runCC(control, [group], 'Holm-Welch');
    const cc = byTreated[0];

    const labels = [...control.map(() => 'ctrl'), ...group.map(() => 'trt')];
    const cats = DG.Column.fromStrings('g', labels);
    const vals = valueCol('f', [...control, ...group]);
    const welch = twoSampleTTest(cats, vals, cats.categories.length, {method: 'Welch', alpha: ALPHA});

    close(cc.df, welch.df, 'Welch df', 0, 1e-10);
    close(cc.statistic, welch.t, 'Welch t', 0, 1e-10);
    close(cc.pValueAdj, welch.pValue, 'adj p == raw (m=1)', 0, 1e-12);
    // Unequal-variance Welch df is far below the pooled n0+n1−2 = 9.
    expect(cc.df < 6, true, `df pushed toward small, got ${cc.df}`);
  }, {timeout: 30000});

  // For balanced equal-variance data both methods agree on the verdict; Dunnett is slightly tighter
  // (smaller adjusted p) because it exploits the shared-control correlation.
  test('Dunnett vs Holm-Welch: agree on verdict, Dunnett tighter', async () => {
    const control = [10.0, 11.0, 9.0, 10.5, 9.5, 10.2];
    const treated = [[13.0, 14.0, 12.0, 13.5, 12.5, 13.2], [16.0, 17.0, 15.0, 16.5, 15.5, 16.2]];
    const dun = runCC(control, treated, 'Dunnett').byTreated;
    const hw = runCC(control, treated, 'Holm-Welch').byTreated;

    for (let i = 0; i < treated.length; ++i) {
      expect(dun[i].significant, hw[i].significant, `verdict agrees for group ${i}`);
      expect(dun[i].pValueAdj <= hw[i].pValueAdj + 1e-9, true,
        `Dunnett adj p (${dun[i].pValueAdj}) ≤ Holm-Welch (${hw[i].pValueAdj}) for group ${i}`);
    }
  }, {timeout: 30000});

  // Mixed verdicts: one group clearly differs, another does not — conclusionColumnPerRow must render
  // both labels in row order.
  test('mixed verdicts populate conclusionColumnPerRow', async () => {
    const control = [10.0, 10.2, 9.8, 10.1, 9.9, 10.0];
    const treated = [[20.0, 21.0, 19.0, 20.5, 19.5, 20.2], [10.1, 9.9, 10.2, 9.8, 10.0, 10.1]];
    const {byTreated} = runCC(control, treated, 'Dunnett');

    expect(byTreated[0].significant, true, 'far group significant');
    expect(byTreated[1].significant, false, 'near group not significant');

    const col = conclusionColumnPerRow(byTreated.map((c) => c.significant));
    expect(col.name, CONCLUSION_COL_NAME);
    expect(col.get(0), CONCLUSION_LABEL.significant);
    expect(col.get(1), CONCLUSION_LABEL.notSignificant);
  }, {timeout: 30000});
});

// ── Holm family size: m = k−1, not k(k−1)/2 ──────────────────────

/** Holm step-down applied with an explicit (possibly wrong) family size `m`, ascending input order. */
function holmWithFamily(pValues: readonly number[], m: number): number[] {
  const order = pValues.map((_, i) => i).sort((a, b) => pValues[a] - pValues[b]);
  const adj = new Array<number>(pValues.length);
  let running = 0;
  for (let rank = 0; rank < pValues.length; ++rank) {
    running = Math.max(running, (m - rank) * pValues[order[rank]]);
    adj[order[rank]] = Math.min(running, 1);
  }
  return adj;
}

// ── External references (scipy.stats / statsmodels) ──────────────
// Generated by tools/generate-control-comparisons-fixtures.py (scipy 1.17, statsmodels 0.14).
// Dunnett adjusted p uses scipy's QMC, so its tolerance is loose (~1e-3) against our Simpson port;
// everything else (Welch t/df/p/CI, Holm) matches scipy/statsmodels near machine precision.

/** Per-field tolerances for the Holm-Welch comparison against scipy/statsmodels. */
const HW_TOL: Record<keyof HolmWelchExpected, [number, number]> = {
  meanDiff: [1e-9, 1e-11],
  statistic: [1e-9, 1e-11],
  df: [1e-9, 1e-11],
  pValueRaw: [1e-7, 1e-12],
  pValueAdj: [1e-7, 1e-12],
  ciLow: [1e-6, 1e-7],
  ciHigh: [1e-6, 1e-7],
  hedgesG: [1e-7, 1e-9],
};

category('Control comparisons: external references', () => {
  for (const fx of FIXTURES) {
    test(`${fx.name} — Dunnett (scipy.stats.dunnett)`, async () => {
      const {byTreated} = runCC(fx.control, fx.treated, 'Dunnett');
      fx.dunnett.forEach((exp, i) => {
        close(byTreated[i].statistic, exp.statistic, `${fx.name}[${i}] statistic`, 1e-9, 1e-9);
        // scipy QMC vs our Simpson integration: match within the shared precision floor.
        close(byTreated[i].pValueAdj, exp.pValueAdj, `${fx.name}[${i}] pValueAdj`, 0, 1.5e-3);
      });
    }, {timeout: 30000});

    test(`${fx.name} — Holm-Welch (ttest_ind + statsmodels holm)`, async () => {
      const {byTreated} = runCC(fx.control, fx.treated, 'Holm-Welch');
      fx.holmWelch.forEach((exp, i) => {
        const c = byTreated[i];
        const actual: HolmWelchExpected = {
          meanDiff: c.meanDiff, statistic: c.statistic, df: c.df,
          pValueRaw: c.pValueRaw, pValueAdj: c.pValueAdj,
          ciLow: c.ciLow, ciHigh: c.ciHigh, hedgesG: c.hedgesG,
        };
        for (const k of Object.keys(exp) as (keyof HolmWelchExpected)[]) {
          const [rtol, atol] = HW_TOL[k];
          close(actual[k], exp[k], `${fx.name}[${i}] ${k}`, rtol, atol);
        }
      });
    }, {timeout: 30000});
  }
});

category('Control comparisons: Holm family size', () => {
  // k=5 (control + 4 groups): Holm must use m = k−1 = 4, NOT k(k−1)/2 = 10. A regression here would
  // over-correct the p-values. We re-derive the adjusted p's from the report's own raw p's under both
  // family sizes and assert the report matches m=4 and differs from m=10.
  test('k=5: family size is k−1, not k(k−1)/2', async () => {
    const control = [5.0, 5.2, 4.8, 5.1, 4.9, 5.05];
    const treated = [
      [9.0, 9.5, 8.5, 9.2, 8.8, 9.1], // strongly different → tiny raw p
      [6.0, 6.3, 5.7, 6.1, 5.9, 6.05],
      [5.6, 5.9, 5.3, 5.7, 5.5, 5.65],
      [5.3, 5.5, 5.1, 5.4, 5.2, 5.35],
    ];
    const {byTreated} = runCC(control, treated, 'Holm-Welch');
    expect(byTreated.length, 4, 'four comparisons');

    const rawPs = byTreated.map((c) => c.pValueRaw);
    const adjPs = byTreated.map((c) => c.pValueAdj);

    const expectedCorrect = holmWithFamily(rawPs, 4); // m = k−1

    for (let i = 0; i < 4; ++i)
      close(adjPs[i], expectedCorrect[i], `adj p[${i}] uses m=k−1`, 0, 1e-12);

    // The smallest raw p pins the family size: its adjusted value is exactly m × raw, so the
    // ratio is k−1 = 4 — and unmistakably not k(k−1)/2 = 10. (Absolute differences are useless
    // here: both products are minuscule.)
    const minIdx = rawPs.indexOf(Math.min(...rawPs));
    const ratio = adjPs[minIdx] / rawPs[minIdx];
    close(ratio, 4, 'smallest adj/raw ratio == k−1', 0, 1e-9);
    expect(Math.abs(ratio - 10) > 1, true, `ratio ${ratio} must not be the k(k−1)/2 family size`);
  }, {timeout: 30000});
});
