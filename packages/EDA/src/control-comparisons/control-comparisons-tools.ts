// Control comparisons (many-to-one) - computations.
//
// Each non-control group is compared against a single control. Two methods:
//
//  * Dunnett — classical pooled-variance procedure. The per-comparison t-statistics share the
//    control, so their adjusted p-values come from the joint distribution of the contrasts under
//    H0, computed via 2D numerical integration (`dunnettMaxCdf`).
//
//  * Holm-Welch — pairwise Welch's t-tests of each group vs the control, with Holm step-down
//    correction applied to the resulting k−1 raw p-values.
//
// The Welch kernel is reused as-is from the package's t-test (`twoSampleTTestFromStats`) — not
// reimplemented. The Dunnett joint-distribution CDF is ported from `@datagrok-libraries/sci-comp`
// (`stats/tests/dunnett.ts`), with sci-comp's `gammaln`/`normalCdf` swapped for the package's
// `jStat`. Effect size uses the package's single Γ-exact `hedgesCorrection` for both methods.

/* REFERENCES

   [1] C.W. Dunnett. A multiple comparison procedure for comparing several treatments with a
       control. JASA 50 (272), 1955, pp. 1096-1121.
   [2] S. Holm. A simple sequentially rejective multiple test procedure. Scandinavian Journal of
       Statistics 6 (2), 1979, pp. 65-70.
   [3] Welch's t-test, https://en.wikipedia.org/wiki/Welch%27s_t-test
*/

//@ts-ignore: no types
import * as jStat from 'jstat';

import {
  CatCol, NumCol, GroupStats,
  factorize, groupStats, fSurvival, checkSignificanceLevel,
} from '../anova/anova-tools';
import {twoSampleTTestFromStats, hedgesCorrection} from '../ttest/ttest-tools';

/** Control-comparisons method. `Dunnett` (default) pools variance; `Holm-Welch` does not. */
export type ControlComparisonsMethod = 'Dunnett' | 'Holm-Welch';

enum ERROR_MSG {
  TOO_FEW_GROUPS = 'control comparisons require a control plus at least one other group',
  CONTROL_ABSENT = 'the selected control group is not present in the data',
  CONTROL_DEGENERATE = 'the control group needs at least 2 observations to estimate variance',
  GROUP_DEGENERATE = 'each group needs at least 2 observations to estimate variance',
  POOLED_DF = 'not enough residual degrees of freedom for the pooled-variance (Dunnett) estimate',
}

/** Result of one (group vs control) comparison. Numbers are signed as `group − control`. */
export interface ControlComparison {
  /** Category code of the group (maps to a label in the UI). */
  groupCode: number;
  /** Group sample size. */
  n: number;
  /** Group sample mean. */
  mean: number;
  /** Mean difference, `mean(group) − mean(control)`. */
  meanDiff: number;
  /** Lower bound of the per-comparison `(1 − alpha)` CI for the mean difference. */
  ciLow: number;
  /** Upper bound of the per-comparison `(1 − alpha)` CI for the mean difference. */
  ciHigh: number;
  /** Test statistic: Dunnett's t (pooled SE) or Welch's t. Signed. */
  statistic: number;
  /** Degrees of freedom: pooled (Dunnett) or Welch–Satterthwaite (Holm-Welch). */
  df: number;
  /** Uncorrected two-sided p-value for this comparison. */
  pValueRaw: number;
  /** Adjusted p-value: Dunnett joint distribution / Holm step-down. */
  pValueAdj: number;
  /** Per-pair Hedges' g (Γ-exact bias correction). */
  hedgesG: number;
  /** `pValueAdj < alpha`. */
  significant: boolean;
}

/** Control-comparisons report. `comparisons` is in category-code order, control excluded. */
export interface ControlComparisonsReport {
  method: ControlComparisonsMethod;
  /** Category code of the control group. */
  controlCode: number;
  /** Control sample size. */
  controlN: number;
  /** Control sample mean. */
  controlMean: number;
  comparisons: ControlComparison[];
  alpha: number;
  /** Pooled within-group mean squared error (Dunnett only). */
  pooledMSE?: number;
  /** Pooled degrees of freedom, `N − k` (Dunnett only). */
  pooledDF?: number;
}

/** Options for `controlComparisons`. */
export interface ControlComparisonsOptions {
  /** Method: 'Dunnett' (default) or 'Holm-Welch'. */
  method?: ControlComparisonsMethod;
  /** Significance level. Default 0.05. */
  alpha?: number;
}

/**
 * Holm's step-down correction (Holm 1979 [2]).
 *
 * Sort p ascending, scale the i-th smallest by `(m − i)`, enforce monotonicity (each adjusted
 * value ≥ the previous), clamp to 1, and return in the original index order. Applied to exactly
 * the k−1 control-vs-each p-values — the family size `m` is their count, not k(k−1)/2.
 */
export function holmCorrect(pValues: readonly number[]): number[] {
  const m = pValues.length;
  const order = Array.from({length: m}, (_, i) => i).sort((a, b) => pValues[a] - pValues[b]);

  const adj = new Array<number>(m);
  let running = 0;
  for (let rank = 0; rank < m; ++rank) {
    const idx = order[rank];
    const scaled = (m - rank) * pValues[idx];
    running = Math.max(running, scaled); // step-down monotonicity
    adj[idx] = Math.min(running, 1);
  }
  return adj;
}

/**
 * Run control comparisons of `values` grouped by `categories` against the group `controlCode`.
 *
 * `uniqueCount` is the factorization bin count (category codes upper bound). Reuses the ANOVA
 * factorization; both methods are computed from per-group Welford stats (mean, M2, size).
 */
export function controlComparisons(categories: CatCol, values: NumCol, controlCode: number,
  uniqueCount: number, opts: ControlComparisonsOptions = {}): ControlComparisonsReport {
  const method = opts.method ?? 'Dunnett';
  const alpha = opts.alpha ?? 0.05;
  checkSignificanceLevel(alpha);

  const f = factorize(categories, values, uniqueCount);

  // Non-empty groups in category-code order.
  const present: {code: number, stats: GroupStats}[] = [];
  for (let i = 0; i < f.catCount; ++i) {
    if (f.sizes[i] > 0)
      present.push({code: i, stats: groupStats(f, i)});
  }

  const control = present.find((p) => p.code === controlCode);
  if (control === undefined)
    throw new Error(ERROR_MSG.CONTROL_ABSENT);
  if (control.stats.size < 2)
    throw new Error(ERROR_MSG.CONTROL_DEGENERATE);

  const groups = present.filter((p) => p.code !== controlCode);
  if (groups.length < 1)
    throw new Error(ERROR_MSG.TOO_FEW_GROUPS);
  if (groups.some((g) => g.stats.size < 2))
    throw new Error(ERROR_MSG.GROUP_DEGENERATE);

  return method === 'Dunnett' ?
    dunnett(control.stats, controlCode, groups, alpha) :
    holmWelch(control.stats, controlCode, groups, alpha);
} // controlComparisons

// ── Dunnett ──────────────────────────────────────────────────────

/** Dunnett's many-to-one comparison: pooled variance + joint reference distribution. */
function dunnett(control: GroupStats, controlCode: number,
  groups: {code: number, stats: GroupStats}[], alpha: number): ControlComparisonsReport {
  const n0 = control.size;
  const meanCtrl = control.mean;

  // Pooled within-group variance across control + all groups.
  const k = groups.length + 1;
  let n = n0;
  let ss = control.m2;
  for (const g of groups) {n += g.stats.size; ss += g.stats.m2;}
  const pooledDF = n - k;
  if (pooledDF < 1)
    throw new Error(ERROR_MSG.POOLED_DF);

  const mse = ss / pooledDF;
  const s = Math.sqrt(mse);
  const tCrit = jStat.studentt.inv(1 - alpha / 2, pooledDF);

  // Per-contrast t-statistics and the (alpha, beta) weights of the Dunnett structure.
  const ts: number[] = [];
  const alphas: number[] = [];
  const betas: number[] = [];
  for (const g of groups) {
    const ni = g.stats.size;
    const se = s * Math.sqrt(1 / ni + 1 / n0);
    ts.push((g.stats.mean - meanCtrl) / se);
    alphas.push(Math.sqrt(n0 / (n0 + ni)));
    betas.push(Math.sqrt(ni / (n0 + ni)));
  }

  const comparisons: ControlComparison[] = groups.map((g, i) => {
    const ni = g.stats.size;
    const se = s * Math.sqrt(1 / ni + 1 / n0);
    const meanDiff = g.stats.mean - meanCtrl;
    const t = ts[i];
    // Adjusted p: 1 − P(|T_max| < |t_i|) under the joint distribution.
    const pValueAdj = Math.min(1, 1 - dunnettMaxCdf(Math.abs(t), pooledDF, alphas, betas));
    // Raw p: two-sided pooled-t p-value via the t² = F(1, df) identity (stable in the tail).
    const pValueRaw = fSurvival(t * t, 1, pooledDF);
    const hedgesG = (meanDiff / s) * hedgesCorrection(pooledDF);
    const pValueAdjClamped = Math.max(pValueAdj, 0);
    return {
      groupCode: g.code,
      n: ni,
      mean: g.stats.mean,
      meanDiff,
      ciLow: meanDiff - tCrit * se,
      ciHigh: meanDiff + tCrit * se,
      statistic: t,
      df: pooledDF,
      pValueRaw,
      pValueAdj: pValueAdjClamped,
      hedgesG,
      significant: pValueAdjClamped < alpha,
    };
  });

  return {
    method: 'Dunnett',
    controlCode,
    controlN: n0,
    controlMean: meanCtrl,
    comparisons,
    alpha,
    pooledMSE: mse,
    pooledDF,
  };
} // dunnett

// ── Holm-Welch ───────────────────────────────────────────────────

/** Pairwise Welch's t-tests vs control + Holm step-down correction. */
function holmWelch(control: GroupStats, controlCode: number,
  groups: {code: number, stats: GroupStats}[], alpha: number): ControlComparisonsReport {
  // Welch kernel reused from the t-test: control = g0, group = g1, so meanDiff = group − control.
  const welch = groups.map((g) =>
    twoSampleTTestFromStats(control, g.stats, controlCode, g.code, {method: 'Welch', alpha}));

  const rawPs = welch.map((w) => w.pValue);
  const adjPs = holmCorrect(rawPs);

  const comparisons: ControlComparison[] = groups.map((g, i) => {
    const w = welch[i];
    return {
      groupCode: g.code,
      n: g.stats.size,
      mean: g.stats.mean,
      meanDiff: w.meanDiff,
      ciLow: w.ciLow,
      ciHigh: w.ciHigh,
      statistic: w.t,
      df: w.df,
      pValueRaw: rawPs[i],
      pValueAdj: adjPs[i],
      hedgesG: w.hedgesG,
      significant: adjPs[i] < alpha,
    };
  });

  return {
    method: 'Holm-Welch',
    controlCode,
    controlN: control.size,
    controlMean: control.mean,
    comparisons,
    alpha,
  };
} // holmWelch

// ── Ported Dunnett joint-distribution CDF ────────────────────────
// Ported from @datagrok-libraries/sci-comp stats/tests/dunnett.ts. gammaln/normalCdf -> jStat.
// The outer scale window is computed adaptively (see below) rather than the sci-comp fixed
// [1e-3, 6] range, which silently breaks for large df.

/**
 * P(|T_1| < c, ..., |T_k| < c) under H0 with the Dunnett structure.
 *
 *   T_i = (alpha_i Z_i − beta_i Z_0) / s,  (Z_i, Z_0) i.i.d. N(0,1),  s = S/σ ~ chi(df)/√df.
 *
 * Outer Simpson over s, inner Simpson over Z_0 ∈ [−zMax, zMax]; the integrand is the product of
 * the conditional rectangle probabilities.
 *
 * The outer variable s = S/σ ~ chi(df)/√df concentrates around 1 as df grows, with width ~1/√(2·df).
 * A fixed [1e-3, 6] grid (as in sci-comp) then puts only ~2 nodes inside the peak at df ≈ 5·10³ and
 * the integral collapses to garbage (adjusted p-values inconsistent with the raw ones). Instead the
 * window is bracketed by a Wilson–Hilferty normal approximation of the χ² tails, so the panel step
 * tracks the peak width at any df; 200 panels then resolve it across the full df range.
 */
function dunnettMaxCdf(c: number, df: number, alphas: number[], betas: number[]): number {
  const zMax = 8.0;
  const nOuter = 200; // Simpson panels for s
  const nInner = 160; // Simpson panels for z

  // Wilson–Hilferty: χ²_df ≈ df·(1 − 2/(9df) ± Z·√(2/(9df)))³, so s = √(χ²/df) = w^{3/2}.
  // Z = 8 brackets ~1e-15 of the tail mass — far below the ~1e-4 integration precision floor.
  const h = 2 / (9 * df);
  const sd = Math.sqrt(h);
  const loW = 1 - h - 8 * sd;
  const hiW = 1 - h + 8 * sd;
  const sMin = Math.max(loW > 0 ? Math.pow(loW, 1.5) : 0, 1e-4);
  const sMax = Math.pow(hiW, 1.5);

  const outerH = (sMax - sMin) / nOuter;
  const innerH = (2 * zMax) / nInner;

  let total = 0;
  for (let i = 0; i <= nOuter; i++) {
    const s = sMin + i * outerH;
    const fs = chiDensity(s, df);
    if (fs === 0) continue;
    const wOuter = (i === 0 || i === nOuter) ? 1 : (i % 2 === 1 ? 4 : 2);

    let inner = 0;
    for (let j = 0; j <= nInner; j++) {
      const z = -zMax + j * innerH;
      const phi = Math.exp(-0.5 * z * z) * 0.39894228040143268; // 1/√(2π)
      let prod = 1;
      for (let m = 0; m < alphas.length; m++) {
        const upper = (c * s + betas[m] * z) / alphas[m];
        const lower = (-c * s + betas[m] * z) / alphas[m];
        prod *= jStat.normal.cdf(upper, 0, 1) - jStat.normal.cdf(lower, 0, 1);
        if (prod === 0) break;
      }
      const wInner = (j === 0 || j === nInner) ? 1 : (j % 2 === 1 ? 4 : 2);
      inner += wInner * phi * prod;
    }
    inner *= innerH / 3;

    total += wOuter * fs * inner;
  }
  return total * outerH / 3;
} // dunnettMaxCdf

/** Density of s = S/σ where S² · df / σ² ~ χ²_df. */
function chiDensity(s: number, df: number): number {
  if (s <= 0) return 0;
  const logF =
    (df / 2) * Math.log(df) -
    (df / 2 - 1) * Math.LN2 -
    jStat.gammaln(df / 2) +
    (df - 1) * Math.log(s) -
    s * s * df / 2;
  return Math.exp(logF);
} // chiDensity
