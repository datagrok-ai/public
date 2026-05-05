/**
 * Boschloo's exact test for a 2×2 contingency table.
 *
 * Boschloo's test is an unconditional exact test that uses Fisher's exact
 * one-sided p-value as the test statistic and maximises the rejection
 * probability over the nuisance parameter π = P(success). It is uniformly
 * more powerful than Fisher's conditional exact test on the same data, at
 * the cost of higher computational expense.
 *
 * **Convention.** Following `scipy.stats.boschloo_exact`, columns of the
 * input table are interpreted as the two binomial experiments (groups).
 * For `[[a, b], [c, d]]`:
 *
 * ```
 *               group 1   group 2 │
 *  success         a         b    │
 *  failure         c         d    │
 *  ─────────────────────────────────
 *               n1 = a+c   n2 = b+d
 * ```
 *
 * For two-sided tests the result is invariant under transposition, so
 * SEND-style tables (rows-as-groups) yield the same `pValue` either way.
 * One-sided alternatives `'less'` / `'greater'` reference the column
 * orientation: `'less'` ⇔ `p_group1 < p_group2`.
 *
 * References:
 *   - Boschloo R. D. (1970). "Raised conditional level of significance for
 *     the 2×2-table when testing the equality of two probabilities."
 *     Statistica Neerlandica 24(1), 1–9.
 *   - Wikipedia, "Boschloo's test", https://en.wikipedia.org/wiki/Boschloo%27s_test
 *   - SciPy documentation for `scipy.stats.boschloo_exact`.
 */

import {combinationln, hypgeomCdf} from '../distributions';
import {fisherExact2x2} from './fisher-exact';

/** Direction of the alternative hypothesis. */
export type BoschlooAlternative = 'two-sided' | 'less' | 'greater';

/** Optional settings for `boschlooExact`. */
export interface BoschlooSettings {
  /** Direction of the alternative hypothesis. Default `'two-sided'`. */
  alternative?: BoschlooAlternative;
  /**
   * Number of grid points over `π ∈ (0, 1)` used in the supremum search,
   * before a golden-section refinement around the best grid point. Default
   * `32` (matches `scipy.stats.boschloo_exact`'s default Sobol budget). Must
   * be at least 4. Larger values cost CPU but rarely change the p-value
   * beyond ~1e-6 — the rejection-probability surface is smooth and
   * approximately unimodal in π.
   */
  nGrid?: number;
}

/** Result of `boschlooExact`. */
export interface BoschlooResult {
  /**
   * Test statistic — Fisher's exact one-sided p-value for the observed
   * table under the chosen `alternative`. For `'two-sided'`, the statistic
   * is taken from the side that produced the smaller one-sided Boschloo
   * p-value. `NaN` when a column total is zero.
   */
  statistic: number;
  /**
   * Boschloo p-value: `sup_π P(reject | π)` under H₀ : `p_group1 = p_group2`,
   * where the rejection region is `{T ≤ T_obs}` and `T` is Fisher's exact
   * one-sided p-value. For `'two-sided'`, it is `min(1, 2·p_min_one_sided)`.
   * `NaN` when a column total is zero.
   */
  pValue: number;
}

/**
 * Boschloo's exact test on a 2×2 contingency table.
 *
 * @param table 2×2 array of non-negative integer counts. Columns are groups.
 * @param settings Optional `{alternative, nGrid}`.
 * @returns `{statistic, pValue}`. Both are `NaN` when either group total is
 *          zero (matches `scipy.stats.boschloo_exact`).
 *
 * @throws `Error` for malformed (non-2×2) input or negative cells, or for
 *         an unknown `alternative`.
 */
export function boschlooExact(
  table: ReadonlyArray<ReadonlyArray<number>>,
  settings: BoschlooSettings = {},
): BoschlooResult {
  if (table.length !== 2 || table[0].length !== 2 || table[1].length !== 2)
    throw new Error('boschlooExact: expected 2×2 table.');
  for (const row of table) {
    for (const v of row) {
      if (v < 0 || !Number.isFinite(v))
        throw new Error('boschlooExact: cells must be non-negative finite numbers.');
    }
  }
  const alternative: BoschlooAlternative = settings.alternative ?? 'two-sided';
  const nGrid = settings.nGrid ?? 32;
  if (nGrid < 4)
    throw new Error('boschlooExact: nGrid must be >= 4.');

  if (alternative === 'two-sided') {
    const less = oneSided(table, 'less', nGrid);
    const greater = oneSided(table, 'greater', nGrid);
    if (Number.isNaN(less.pValue)) return less;
    const minSide = less.pValue < greater.pValue ? less : greater;
    return {
      statistic: minSide.statistic,
      pValue: clamp01(2 * minSide.pValue),
    };
  }
  if (alternative !== 'less' && alternative !== 'greater')
    throw new Error(`boschlooExact: unknown alternative ${JSON.stringify(alternative)}.`);
  return oneSided(table, alternative, nGrid);
}

/**
 * Convenience wrapper: compute Boschloo (primary) and Fisher (alternative)
 * two-sided p-values together with the sample odds ratio.
 *
 * Matches the SEND `incidence_exact_both` shape: callers can switch
 * between the two p-values without recomputing the sample odds ratio
 * or rebuilding the table.
 */
export function incidenceExactBoth(
  table: ReadonlyArray<ReadonlyArray<number>>,
): {oddsRatio: number; pValue: number; pValueFisher: number} {
  const boschloo = boschlooExact(table, {alternative: 'two-sided'});
  const fisher = fisherExact2x2(table);
  // Per SEND policy: degenerate Boschloo (NaN on zero column) → 1.0
  const pValue = Number.isNaN(boschloo.pValue) ? 1.0 : boschloo.pValue;
  return {
    oddsRatio: fisher.oddsRatio,
    pValue,
    pValueFisher: fisher.pValue,
  };
}

// ─── Internal ───────────────────────────────────────────────────────────

function oneSided(
  table: ReadonlyArray<ReadonlyArray<number>>,
  side: 'less' | 'greater',
  nGrid: number,
): BoschlooResult {
  const a = table[0][0];
  const b = table[0][1];
  const c = table[1][0];
  const d = table[1][1];

  // Column totals are the binomial sample sizes (groups).
  const n1 = a + c;
  const n2 = b + d;
  if (n1 === 0 || n2 === 0)
    return {statistic: NaN, pValue: NaN};

  const total = n1 + n2;

  // Build the (n1+1) × (n2+1) table of Fisher one-sided p-values.
  // For 'less':    pvalues[x1][x2] = P(X ≤ x1 | total, x1+x2, n1)
  // For 'greater': pvalues[x1][x2] = P(X ≤ x2 | total, x1+x2, n2)
  // (Each conditions on the sum x1+x2 and uses one of the two groups.)
  const nRows = n1 + 1;
  const nCols = n2 + 1;
  const pvalues = new Float64Array(nRows * nCols);
  for (let x1 = 0; x1 < nRows; x1++) {
    for (let x2 = 0; x2 < nCols; x2++) {
      const k = x1 + x2;
      const p = side === 'less' ?
        hypgeomCdf(x1, total, k, n1) :
        hypgeomCdf(x2, total, k, n2);
      pvalues[x1 * nCols + x2] = p;
    }
  }

  // The Fisher p-value at the observed cell is the test statistic.
  const fisherStat = pvalues[a * nCols + b];

  // Rejection region: {(x1, x2) : pvalues[x1, x2] ≤ fisherStat · (1 + 1e-11)}.
  // The slack absorbs floating-point noise so that mathematically equivalent
  // hypergeometric CDFs computed via different sums (e.g. cdf(1, 24, 7, 17)
  // vs cdf(3, 24, 10, 17), which both equal 120/C(24,17)) end up in the same
  // rejection region. scipy uses 1e-13 because its hypergeom.cdf is
  // numerically stable to ~1 ulp; jstat's CDF accumulates ~ 1e-12 relative
  // error on tables with several PMF terms, so we use 1e-11 to leave a 10×
  // margin while still rejecting genuinely larger p-values.
  const threshold = fisherStat * (1 + 1e-11);
  const inR = new Uint8Array(nRows * nCols);
  for (let i = 0; i < pvalues.length; i++) inR[i] = pvalues[i] <= threshold ? 1 : 0;

  // Precompute log-binomial coefficients for the two binomial PMFs.
  const logC1 = new Float64Array(nRows);
  for (let x1 = 0; x1 < nRows; x1++) logC1[x1] = combinationln(n1, x1);
  const logC2 = new Float64Array(nCols);
  for (let x2 = 0; x2 < nCols; x2++) logC2[x2] = combinationln(n2, x2);

  // log P(reject | π) via log-sum-exp over the rejection region:
  //   log Bin(x1; n1, π) · Bin(x2; n2, π)
  //     = lnC1[x1] + lnC2[x2] + (x1+x2)·ln π + (total − x1 − x2)·ln(1−π)
  const logPRejection = (pi: number): number => {
    if (pi <= 0 || pi >= 1) return -Infinity;
    const lpi = Math.log(pi);
    const l1mpi = Math.log(1 - pi);
    let logMax = -Infinity;
    // First pass: find max log-term for log-sum-exp stability
    for (let x1 = 0; x1 < nRows; x1++) {
      const base1 = logC1[x1];
      const off = x1 * nCols;
      for (let x2 = 0; x2 < nCols; x2++) {
        if (!inR[off + x2]) continue;
        const k = x1 + x2;
        const t = base1 + logC2[x2] + k * lpi + (total - k) * l1mpi;
        if (t > logMax) logMax = t;
      }
    }
    if (logMax === -Infinity) return -Infinity;
    // Second pass: accumulate exp(t − logMax)
    let sum = 0;
    for (let x1 = 0; x1 < nRows; x1++) {
      const base1 = logC1[x1];
      const off = x1 * nCols;
      for (let x2 = 0; x2 < nCols; x2++) {
        if (!inR[off + x2]) continue;
        const k = x1 + x2;
        const t = base1 + logC2[x2] + k * lpi + (total - k) * l1mpi;
        sum += Math.exp(t - logMax);
      }
    }
    return logMax + Math.log(sum);
  };

  // Grid search over π ∈ (0, 1). The endpoints contribute 0 mass to the
  // interior of R, so we skip i = 0 and i = nGrid.
  let bestPi = 0.5;
  let bestLogP = -Infinity;
  for (let i = 1; i < nGrid; i++) {
    const pi = i / nGrid;
    const lp = logPRejection(pi);
    if (lp > bestLogP) {
      bestLogP = lp;
      bestPi = pi;
    }
  }

  // Golden-section refinement in a window of ±2 grid steps around bestPi.
  const step = 1 / nGrid;
  const lo = Math.max(1e-10, bestPi - 2 * step);
  const hi = Math.min(1 - 1e-10, bestPi + 2 * step);
  const refined = goldenSectionMax(logPRejection, lo, hi, 1e-10, 100);
  if (refined.value > bestLogP) {
    bestLogP = refined.value;
    bestPi = refined.x;
  }

  return {
    statistic: fisherStat,
    pValue: clamp01(Math.exp(bestLogP)),
  };
}

function clamp01(x: number): number {
  if (Number.isNaN(x)) return x;
  return Math.min(1, Math.max(0, x));
}

/**
 * Golden-section search for the maximum of a unimodal function on `[a, b]`.
 *
 * The Boschloo rejection-probability surface is smooth and approximately
 * unimodal in π over a small interval around the grid maximum, so 100
 * iterations (≈ 50 function evaluations after the initial 4) of the golden
 * ratio cut suffice to refine the maximum to 1e-10 in π.
 */
function goldenSectionMax(
  f: (x: number) => number,
  a: number,
  b: number,
  tolerance: number,
  maxIter: number,
): {x: number; value: number} {
  const phi = (Math.sqrt(5) - 1) / 2; // 0.61803...
  let lo = a;
  let hi = b;
  let x1 = hi - phi * (hi - lo);
  let x2 = lo + phi * (hi - lo);
  let f1 = f(x1);
  let f2 = f(x2);
  for (let iter = 0; iter < maxIter && hi - lo > tolerance; iter++) {
    if (f1 > f2) {
      hi = x2;
      x2 = x1;
      f2 = f1;
      x1 = hi - phi * (hi - lo);
      f1 = f(x1);
    } else {
      lo = x1;
      x1 = x2;
      f1 = f2;
      x2 = lo + phi * (hi - lo);
      f2 = f(x2);
    }
  }
  if (f1 > f2) return {x: x1, value: f1};
  return {x: x2, value: f2};
}
