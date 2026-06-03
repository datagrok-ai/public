// Two-sample t-test - computations.
//
// Reuses the ANOVA *factorization* (the Welford split of a feature column by a
// categorical column) but NOT the ANOVA engine — the two-sample report is computed
// independently from the two per-group `GroupStats`.
//
// The p-value uses the identity t² = F(1, df). The two-sided t p-value equals the
// upper-tail probability of an F(1, df) variate at t², so `fSurvival` is reused
// (numerically stable in the tail).

/* REFERENCES

   [1] Student's t-test, https://en.wikipedia.org/wiki/Student%27s_t-test
   [2] Welch's t-test, https://en.wikipedia.org/wiki/Welch%27s_t-test
   [3] Effect size (Cohen's d, Hedges' g), https://en.wikipedia.org/wiki/Effect_size
   [4] M. Delacre, D. Lakens, et al. Why Hedges' g_s* based on the non-pooled
       standard deviation should be reported with Welch's t-test. 2021.
*/

//@ts-ignore: no types
import * as jStat from 'jstat';

import {
  CatCol, NumCol, GroupStats,
  factorize, groupStats, getVariance, fSurvival, checkSignificanceLevel,
} from '../anova/anova-tools';

/** t-test method. `Welch` (default) is robust to unequal variances. */
export type TTestMethod = 'Welch' | 'Student';

enum ERROR_MSG {
  NOT_TWO_GROUPS = 'two-sample t-test requires exactly 2 non-empty groups',
  DEGENERATE_GROUP = 'each group needs at least 2 observations to estimate variance',
  NO_VARIATION = 'no within-group variation — the t-test is undefined',
}

/** Two-sample t-test result. Group order follows the category codes (group 0, group 1);
 *  `meanDiff` and `t` are signed as `group 1 − group 0`. */
export interface TwoSampleTTest {
  /** Method actually used. */
  method: TTestMethod,
  /** Category code of group 0 (maps to a label in the UI). */
  code0: number,
  /** Category code of group 1. */
  code1: number,
  /** Group 0 sample mean. */
  mean0: number,
  /** Group 1 sample mean. */
  mean1: number,
  /** Group 0 sample size. */
  n0: number,
  /** Group 1 sample size. */
  n1: number,
  /** Group 0 unbiased sample variance. */
  var0: number,
  /** Group 1 unbiased sample variance. */
  var1: number,
  /** Mean difference, `mean1 − mean0` (signed). */
  meanDiff: number,
  /** Standard error of the mean difference. */
  se: number,
  /** t-statistic (signed). */
  t: number,
  /** Degrees of freedom (integer for Student, fractional for Welch–Satterthwaite). */
  df: number,
  /** Two-sided p-value, `P(|T_df| > |t|)`. */
  pValue: number,
  /** Significance level used for the confidence interval. */
  alpha: number,
  /** Lower bound of the `(1 − alpha)` CI for the mean difference. */
  ciLow: number,
  /** Upper bound of the `(1 − alpha)` CI for the mean difference. */
  ciHigh: number,
  /** Effect size: Cohen's d (Student, pooled SD) / d_s (Welch, non-pooled SD √((v0+v1)/2)). */
  cohenD: number,
  /** Bias-corrected effect size: Hedges' g (Student) / g_s* (Welch). */
  hedgesG: number,
}

/** Options for `twoSampleTTest`. */
export interface TwoSampleTTestOptions {
  /** Method: 'Welch' (default; robust to unequal variances) or 'Student' (pooled variance). */
  method?: TTestMethod,
  /** Significance level for the confidence interval. Default 0.05. */
  alpha?: number,
}

/**
 * Hedges' small-sample bias correction factor, exact Γ form:
 *   J = Γ(df/2) / (√(df/2) · Γ((df−1)/2)) = exp(lgamma(df/2) − ½·ln(df/2) − lgamma((df−1)/2)).
 *
 * Evaluated at the degrees of freedom of the test itself: n0 + n1 − 2 for Student, and the
 * fractional Welch–Satterthwaite ν for Welch. Using ν (not n0 + n1 − 2) in the Welch branch is
 * what makes the corrected effect size the g_s* of Delacre et al. 2021 [4]; it matches
 * `effectsize::hedges_g(..., pooled_sd = FALSE)` and TOSTER. The exact Γ form (vs. the classic
 * `1 − 3/(4·df − 1)` approximation) matters most at the small fractional df Welch can produce.
 *
 * Exported so Control comparisons reuses the same single effect-size correction for both methods.
 */
export function hedgesCorrection(df: number): number {
  return Math.exp(jStat.gammaln(df / 2) - 0.5 * Math.log(df / 2) - jStat.gammaln((df - 1) / 2));
}

/**
 * Two-sample t-test (Welch by default).
 *
 * Reuses the ANOVA factorization; computes t, df, two-sided p-value, the
 * `(1 − alpha)` confidence interval for the mean difference, and the effect size
 * independently from the two per-group `GroupStats`.
 */
export function twoSampleTTest(categories: CatCol, values: NumCol, uniqueCount: number,
  opts: TwoSampleTTestOptions = {}): TwoSampleTTest {
  const f = factorize(categories, values, uniqueCount);

  // Collect the non-empty groups together with their category codes.
  const present: {code: number, stats: GroupStats}[] = [];
  for (let i = 0; i < f.catCount; ++i) {
    if (f.sizes[i] > 0)
      present.push({code: i, stats: groupStats(f, i)});
  }

  if (present.length !== 2)
    throw new Error(`${ERROR_MSG.NOT_TWO_GROUPS} (got ${present.length}).`);

  return twoSampleTTestFromStats(present[0].stats, present[1].stats,
    present[0].code, present[1].code, opts);
} // twoSampleTTest

/**
 * Two-sample t-test computed directly from two per-group `GroupStats` (Welford mean/M2/size).
 *
 * The column-based `twoSampleTTest` is a thin factorization wrapper over this. Exposed so
 * Control comparisons can reuse the exact same Welch kernel per (group vs control) pair —
 * `g0` = control, `g1` = group, so `meanDiff` and `t` are signed as group − control.
 */
export function twoSampleTTestFromStats(g0: GroupStats, g1: GroupStats, code0: number, code1: number,
  opts: TwoSampleTTestOptions = {}): TwoSampleTTest {
  const method = opts.method ?? 'Welch';
  const alpha = opts.alpha ?? 0.05;
  checkSignificanceLevel(alpha);

  const n0 = g0.size;
  const n1 = g1.size;

  if (n0 < 2 || n1 < 2)
    throw new Error(ERROR_MSG.DEGENERATE_GROUP);

  const var0 = getVariance(g0);
  const var1 = getVariance(g1);
  const meanDiff = g1.mean - g0.mean;

  let se: number;
  let df: number;
  let stdForEffect: number; // standardizer for the effect size

  if (method === 'Student') {
    df = n0 + n1 - 2;
    const sp2 = ((n0 - 1) * var0 + (n1 - 1) * var1) / df;
    se = Math.sqrt(sp2 * (1 / n0 + 1 / n1));
    stdForEffect = Math.sqrt(sp2); // pooled SD
  } else {
    const a = var0 / n0;
    const b = var1 / n1;
    se = Math.sqrt(a + b);
    // Welch–Satterthwaite degrees of freedom (fractional).
    df = (a + b) * (a + b) / (a * a / (n0 - 1) + b * b / (n1 - 1));
    stdForEffect = Math.sqrt((var0 + var1) / 2); // non-pooled averaged SD (g_s*, see [4])
  }

  if (!(se > 0) || !(stdForEffect > 0))
    throw new Error(ERROR_MSG.NO_VARIATION);

  const t = meanDiff / se;
  // Two-sided t p-value via t² = F(1, df); fSurvival is stable in the tail.
  const pValue = fSurvival(t * t, 1, df);

  const tCrit = jStat.studentt.inv(1 - alpha / 2, df);
  const ciLow = meanDiff - tCrit * se;
  const ciHigh = meanDiff + tCrit * se;

  const cohenD = meanDiff / stdForEffect;
  const hedgesG = cohenD * hedgesCorrection(df);

  return {
    method,
    code0,
    code1,
    mean0: g0.mean,
    mean1: g1.mean,
    n0, n1,
    var0, var1,
    meanDiff,
    se,
    t,
    df,
    pValue,
    alpha,
    ciLow,
    ciHigh,
    cohenD,
    hedgesG,
  };
} // twoSampleTTestFromStats
