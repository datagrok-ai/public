// Tests for the two-sample t-test.
//
// Numeric fixtures are cross-checked against R 4.6 `stats::t.test` (and `oneway.test`
// for the t² = F identity). Hedges' g / g_s* expected values use the exact-Γ correction
// J(df) = exp(lgamma(df/2) − ½·ln(df/2) − lgamma((df−1)/2)) evaluated at the test's own df
// (n0 + n1 − 2 for Student, the Welch–Satterthwaite ν for Welch) — i.e. Delacre et al. 2021.
// Regenerate with `Rscript` using `t.test(x1, x0, var.equal=<Student>, conf.level=0.95)`.

import * as DG from 'datagrok-api/dg';

import {category, expect, test} from '@datagrok-libraries/test/src/test';

import {twoSampleTTest, TwoSampleTTest, TTestMethod} from '../ttest/ttest-tools';
import {oneWayAnova} from '../anova/anova-tools';

const ALPHA = 0.05;

/** Build a string category column from integer-encoded labels.
 *  The 'g' prefix is required: DG.Column.fromStrings auto-promotes purely numeric
 *  strings to an INT column, which breaks the factorization indexing. */
function factorCol(name: string, data: readonly number[]): DG.Column {
  return DG.Column.fromStrings(name, data.map((n) => `g${n}`));
}

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

/** Run a two-sample t-test on (g0, g1): g0 values get category code 0, g1 values code 1.
 *  `meanDiff` and `t` are therefore signed as mean(g1) − mean(g0). */
function runTTest(g0: readonly number[], g1: readonly number[], method: TTestMethod): TwoSampleTTest {
  const codes = [...g0.map(() => 0), ...g1.map(() => 1)];
  const cats = factorCol('group', codes);
  const vals = valueCol('feature', [...g0, ...g1]);
  return twoSampleTTest(cats, vals, cats.categories.length, {method, alpha: ALPHA});
}

/** Per-field [rtol, atol]. Defaults reflect the jStat accuracy observed against R:
 *  arithmetic quantities ~1e-13, p-value ~1e-10, CI ~1e-9 (limited by `studentt.inv`). */
type FieldTol = Partial<Record<keyof TwoSampleTTest, [number, number]>>;
const DEFAULT_TOL: FieldTol = {
  mean0: [1e-10, 1e-12], mean1: [1e-10, 1e-12],
  var0: [1e-10, 1e-12], var1: [1e-10, 1e-12],
  meanDiff: [1e-10, 1e-12], se: [1e-10, 1e-12],
  t: [1e-9, 1e-11], df: [1e-9, 1e-11],
  pValue: [1e-7, 1e-12],
  ciLow: [1e-7, 1e-8], ciHigh: [1e-7, 1e-8],
  cohenD: [1e-8, 1e-10], hedgesG: [1e-8, 1e-10],
};

/** Assert every field present in `exp`, using `tol` overrides on top of `DEFAULT_TOL`. */
function checkAll(r: TwoSampleTTest, exp: Partial<TwoSampleTTest>, tol: FieldTol = {}): void {
  for (const k of Object.keys(exp) as (keyof TwoSampleTTest)[]) {
    const expected = exp[k];
    if (typeof expected !== 'number')
      continue;
    const [rtol, atol] = tol[k] ?? DEFAULT_TOL[k] ?? [1e-8, 1e-10];
    close(r[k] as number, expected, k, rtol, atol);
  }
}

interface Fixture {
  name: string;
  g0: number[];
  g1: number[];
  method: TTestMethod;
  exp: Partial<TwoSampleTTest>;
  tol?: FieldTol;
}

// Iris Petal.Length (R `datasets::iris`): versicolor (group 0) vs virginica (group 1).
const IRIS_VERSICOLOR = [4.7, 4.5, 4.9, 4.0, 4.6, 4.5, 4.7, 3.3, 4.6, 3.9, 3.5, 4.2, 4.0, 4.7, 3.6,
  4.4, 4.5, 4.1, 4.5, 3.9, 4.8, 4.0, 4.9, 4.7, 4.3, 4.4, 4.8, 5.0, 4.5, 3.5,
  3.8, 3.7, 3.9, 5.1, 4.5, 4.5, 4.7, 4.4, 4.1, 4.0, 4.4, 4.6, 4.0, 3.3, 4.2, 4.2, 4.2, 4.3, 3.0, 4.1];
const IRIS_VIRGINICA = [6.0, 5.1, 5.9, 5.6, 5.8, 6.6, 4.5, 6.3, 5.8, 6.1, 5.1, 5.3, 5.5, 5.0, 5.1,
  5.3, 5.5, 6.7, 6.9, 5.0, 5.7, 4.9, 6.7, 4.9, 5.7, 6.0, 4.8, 4.9, 5.6, 5.8,
  6.1, 6.4, 5.6, 5.1, 5.6, 6.1, 5.6, 5.5, 4.8, 5.4, 5.6, 5.1, 5.1, 5.9, 5.7, 5.2, 5.0, 5.2, 5.4, 5.1];

// 30 deterministic values (n0=3 vs n1=30); R `var` = 1.5720574712643678.
const N30 = [0.5, 1.2, -0.3, 2.1, 0.0, 1.7, -1.1, 0.8, 3.2, -0.5, 1.0, 2.4, 0.3, -0.7, 1.9,
  0.6, 2.8, -1.4, 1.3, 0.1, 2.0, -0.2, 1.6, 0.9, 3.0, -0.9, 1.1, 0.4, 2.6, 0.7];

const FIXTURES: Fixture[] = [
  // --- 1. Canonical -------------------------------------------------------------------
  {
    name: '1.1 student_5v5 — smallest exact example',
    g0: [1, 2, 3, 4, 5], g1: [6, 7, 8, 9, 10], method: 'Student',
    // R: t.test(6:10, 1:5, var.equal=TRUE)
    exp: {mean0: 3, mean1: 8, var0: 2.5, var1: 2.5, meanDiff: 5, se: 1, t: 5, df: 8,
      pValue: 0.001052825793366539, ciLow: 2.6939958647958338, ciHigh: 7.3060041352041658,
      cohenD: 3.1622776601683791, hedgesG: 2.8545985858444332},
  },
  {
    name: '1.2 welch_5v5_equalvar — Welch ν collapses to 8',
    g0: [1, 2, 3, 4, 5], g1: [6, 7, 8, 9, 10], method: 'Welch',
    // R: t.test(6:10, 1:5) → df = 8 (exact, balanced + equal var)
    exp: {meanDiff: 5, se: 1, t: 5, df: 8, pValue: 0.0010528257933665405,
      cohenD: 3.1622776601683791, hedgesG: 2.8545985858444332},
    tol: {df: [0, 1e-12]}, // must be machine-exact
  },
  {
    name: '1.3 iris_versicolor_vs_virginica — canonical dataset',
    g0: IRIS_VERSICOLOR, g1: IRIS_VIRGINICA, method: 'Welch',
    // R: t.test(virginica, versicolor)
    exp: {mean0: 4.26, mean1: 5.552, var0: 0.22081632653061228, var1: 0.30458775510204084,
      meanDiff: 1.2919999999999998, t: 12.60377944138498, df: 95.570435006454076,
      pValue: 4.9002875273982702e-22, ciLow: 1.0885097008057105, ciHigh: 1.4954902991942889,
      cohenD: 2.5207558882769958, hedgesG: 2.5009133923884228},
    tol: {pValue: [1e-6, 0]}, // deep tail: jStat ibeta vs R pt
  },
  // --- 2. Welch, unequal variances ----------------------------------------------------
  {
    name: '2.1 var_ratio_100_to_1',
    g0: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
    g1: [10, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9], method: 'Welch',
    // R: t.test(seq(10,10.9,0.1), 0:9)
    exp: {var0: 9.1666666666666661, var1: 0.091666666666666743, meanDiff: 5.9499999999999993,
      t: 6.183730577279805, df: 9.1799820017998179, pValue: 0.00014905402565134119,
      ciLow: 3.7798349899719264, ciHigh: 8.1201650100280727,
      cohenD: 2.765448385068332, hedgesG: 2.5320963157394512},
  },
  // --- 3. Unbalanced n ----------------------------------------------------------------
  {
    name: '3.1 n_3_vs_30',
    g0: [1, 2, 3], g1: N30, method: 'Welch',
    exp: {mean0: 2, var0: 1, var1: 1.5720574712643678, meanDiff: -1.0966666666666667,
      t: -1.7657531092747702, df: 2.6736932677843557, pValue: 0.18669571774168153,
      ciLow: -3.2169344643439155, ciHigh: 1.0236011310105824,
      cohenD: -0.9670508407520696, hedgesG: -0.66342644147419794},
  },
  // --- 4. Null / identical means ------------------------------------------------------
  {
    name: '4.1 identical_groups — t = 0, p = 1',
    g0: [1, 2, 3, 4, 5], g1: [1, 2, 3, 4, 5], method: 'Student',
    exp: {meanDiff: 0, t: 0, df: 8, pValue: 1, ciLow: -2.3060041352041662,
      ciHigh: 2.3060041352041662, cohenD: 0, hedgesG: 0},
    tol: {t: [0, 1e-15], pValue: [0, 1e-12], meanDiff: [0, 1e-15]},
  },
  {
    name: '4.2 identical_means_unequal_var — t ≈ 0, Welch ν ≠ 8',
    g0: [4.99, 5.0, 5.01, 4.98, 5.02], g1: [3, 4, 5, 6, 7], method: 'Welch',
    // Means are exactly 5.0; variances differ ~10000×, so ν ≈ 4.0008 (not 8).
    exp: {meanDiff: 0, t: 0, df: 4.000799999992001, pValue: 1,
      ciLow: -1.9631865012897638, ciHigh: 1.9631865012897638, cohenD: 0, hedgesG: 0},
    tol: {t: [0, 1e-12], pValue: [0, 1e-12], meanDiff: [0, 1e-12]},
  },
  // --- 5. Extreme-tail p-value --------------------------------------------------------
  {
    name: '5.1 large_t_tail — |t| ~ 7385, p ~ 1e-60',
    g0: [0, 0.0001, -0.0001, 0.0002, -0.0002, 0.0003, -0.0003, 0.0004, -0.0004, 0.0005],
    g1: [1, 1.0001, 0.9999, 1.0002, 0.9998, 1.0003, 0.9997, 1.0004, 0.9996, 1.0005], method: 'Student',
    // R: t.test(g0+1, g0, var.equal=TRUE) → p = 8.6079406922972295e-60
    exp: {meanDiff: 0.99999999999999989, t: 7385.4894587603703, df: 18, pValue: 8.6079406922972295e-60},
    tol: {pValue: [1e-6, 0]}, // different algorithm paths in the deep tail
  },
  // --- 6. Welch df at the pooled boundary ---------------------------------------------
  {
    name: '6.2 welch_df_at_pooled — ν = 18 exactly',
    g0: [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], g1: [11, 12, 13, 14, 15, 16, 17, 18, 19, 20], method: 'Welch',
    exp: {meanDiff: 10, t: 7.385489458759964, df: 18, pValue: 7.5031381721089191e-07,
      ciLow: 7.1553381099892785, ciHigh: 12.84466189001072,
      cohenD: 3.3028912953790819, hedgesG: 3.1630020478675243},
    tol: {df: [1e-12, 1e-11]},
  },
  // --- 7. Sign convention (signed mean(g1) − mean(g0)) --------------------------------
  {
    name: '7.1 group0_smaller — meanDiff = +9',
    g0: [1, 2, 3], g1: [10, 11, 12], method: 'Welch',
    exp: {meanDiff: 9, t: 11.022703842524301, df: 4.0000000000000009, pValue: 0.00038506771136654082,
      cohenD: 9, hedgesG: 7.1809610472257877},
  },
  // --- 8. Mixed sign / large magnitude ------------------------------------------------
  {
    name: '8.1 negative_values — mirror of 1.1, meanDiff = −5',
    g0: [-5, -4, -3, -2, -1], g1: [-10, -9, -8, -7, -6], method: 'Student',
    exp: {meanDiff: -5, t: -5, df: 8, pValue: 0.001052825793366539, cohenD: -3.1622776601683791,
      hedgesG: -2.8545985858444332},
  },
  {
    name: '8.2 large_magnitudes — Welford under near-cancellation',
    g0: [1e9, 1e9 + 1, 1e9 + 2], g1: [1e9 + 10, 1e9 + 11, 1e9 + 12], method: 'Student',
    exp: {var0: 1, var1: 1, meanDiff: 10, t: 12.24744871391589, df: 4,
      pValue: 0.00025521674944192687, cohenD: 10, hedgesG: 7.9788456080286529},
    tol: {var0: [1e-7, 1e-7], var1: [1e-7, 1e-7]}, // 1e9 scale
  },
  // --- 9. Welford robustness corner ---------------------------------------------------
  {
    name: '9.1 welford_stability — shifted+scaled copy of 1.1; t, df must match',
    g0: [1 + 1e-10, 1 + 2e-10, 1 + 3e-10, 1 + 4e-10, 1 + 5e-10],
    g1: [1 + 6e-10, 1 + 7e-10, 1 + 8e-10, 1 + 9e-10, 1 + 10e-10], method: 'Student',
    // t and df are shift/scale-invariant → identical to 1.1. Naïve Σx²−n·x̄² loses ~10 digits here.
    exp: {t: 5, df: 8, pValue: 0.001052825793366539, cohenD: 3.1622776601683795,
      hedgesG: 2.8545985858444332},
    tol: {t: [1e-9, 1e-11], cohenD: [1e-7, 1e-9]},
  },
  // --- 10. One constant group is NOT degenerate — the test stays well-defined ----------
  {
    name: '10.4 one_constant_group — finite result (only the OTHER group supplies variance)',
    g0: [5, 5, 5, 5, 5], g1: [1, 2, 3, 4, 5], method: 'Welch',
    // se = √(0 + var1/n1) > 0, so NO_VARIATION does NOT fire (matches R/SciPy, not a throw).
    // R: t.test(1:5, rep(5,5))
    exp: {var0: 0, var1: 2.5, meanDiff: -2, se: 0.70710678118654757, t: -2.8284271247461898,
      df: 4, pValue: 0.04742065558431962, ciLow: -3.963243161477557, ciHigh: -0.03675683852244279,
      cohenD: -1.7888543819998317, hedgesG: -1.4272992929222168},
  },
];

category('t-test: fixtures', () => {
  for (const fx of FIXTURES) {
    test(fx.name, async () => {
      const r = runTTest(fx.g0, fx.g1, fx.method);
      expect(r.method, fx.method);
      checkAll(r, fx.exp, fx.tol);
    }, {timeout: 10000});
  }
});

category('t-test: t² = F identity', () => {
  test('11.1 Student t² == Fisher ANOVA F (unbalanced 3 vs 4)', async () => {
    const g0 = [1, 2, 3]; const g1 = [4, 5, 6, 7];
    const r = runTTest(g0, g1, 'Student');
    // R: t.test(4:7, 1:3, var.equal=TRUE) and oneway.test(..., var.equal=TRUE)
    close(r.t, 3.872983346207417, 't', 1e-9);
    close(r.df, 5, 'df', 0, 1e-12);
    close(r.pValue, 0.011724811003954635, 'pValue', 1e-7);

    const codes = [0, 0, 0, 1, 1, 1, 1];
    const cats = factorCol('group', codes);
    const vals = valueCol('feature', [...g0, ...g1]);
    const anova = oneWayAnova(cats, vals, ALPHA, {method: 'Fisher', toValidate: false});
    if (anova.method !== 'Fisher')
      throw new Error('Expected Fisher report');
    close(anova.anovaTable.fStat, r.t * r.t, 'F == t²', 1e-12);
    expect(anova.anovaTable.dfBn, 1, 'Fisher dfBn must be 1 for k=2');
    expect(anova.anovaTable.dfWn, 5, 'Fisher dfWn');
    close(anova.anovaTable.pValue, r.pValue, 'p-value (Fisher == Student t)', 1e-9);
  }, {timeout: 10000});

  test('11.2 Welch t² == Welch ANOVA W (unbalanced 3 vs 4)', async () => {
    const g0 = [1, 2, 3]; const g1 = [4, 5, 6, 7];
    const r = runTTest(g0, g1, 'Welch');

    const codes = [0, 0, 0, 1, 1, 1, 1];
    const cats = factorCol('group', codes);
    const vals = valueCol('feature', [...g0, ...g1]);
    const anova = oneWayAnova(cats, vals, ALPHA, {method: 'Welch', toValidate: false});
    if (anova.method !== 'Welch')
      throw new Error('Expected Welch report');
    // R: oneway.test(..., var.equal=FALSE) → W = 16.3333…, denom df = 4.95918…, p = 0.010076943…
    close(anova.anovaTable.fStat, r.t * r.t, 'W == t²', 1e-12);
    close(anova.anovaTable.dfWn, r.df, 'Welch dfWn == t-test df', 1e-9);
    close(anova.anovaTable.pValue, r.pValue, 'p-value (Welch ANOVA == Welch t)', 1e-9);
    close(r.df, 4.9591836734693864, 't-test ν (R oneway.test)', 1e-9);
  }, {timeout: 10000});
});

category('t-test: guards', () => {
  test('10.2 more than two groups throws', async () => {
    const cats = factorCol('group', [0, 0, 1, 1, 2, 2]);
    const vals = valueCol('feature', [1, 2, 3, 4, 5, 6]);
    let threw = false;
    try {
      twoSampleTTest(cats, vals, cats.categories.length, {method: 'Welch', alpha: ALPHA});
    } catch (e) {
      threw = true;
      expect((e as Error).message.includes('exactly 2'), true, 'message should mention exactly 2 groups');
    }
    expect(threw, true, 'three groups must throw');
  }, {timeout: 10000});

  test('10.3 degenerate group (n < 2) throws', async () => {
    const cats = factorCol('group', [0, 1, 1, 1]);
    const vals = valueCol('feature', [1, 2, 3, 4]);
    let threw = false;
    try {
      twoSampleTTest(cats, vals, cats.categories.length, {method: 'Welch', alpha: ALPHA});
    } catch (e) {
      threw = true;
      expect((e as Error).message.includes('2 observations'), true, 'message should mention 2 observations');
    }
    expect(threw, true, 'a singleton group must throw');
  }, {timeout: 10000});

  test('10.5 both zero-variance groups throw NO_VARIATION', async () => {
    // NO_VARIATION fires only when *both* groups are constant (se = 0). A single constant
    // group still yields a finite test — see fixture 10.4.
    const cats = factorCol('group', [0, 0, 0, 1, 1, 1]);
    const vals = valueCol('feature', [5, 5, 5, 7, 7, 7]);
    let threw = false;
    try {
      twoSampleTTest(cats, vals, cats.categories.length, {method: 'Student', alpha: ALPHA});
    } catch (e) {
      threw = true;
      expect((e as Error).message.includes('variation'), true, 'message should mention no within-group variation');
    }
    expect(threw, true, 'both-constant groups must throw');
  }, {timeout: 10000});
});
