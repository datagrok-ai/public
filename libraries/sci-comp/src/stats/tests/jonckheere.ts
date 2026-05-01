/**
 * Jonckheere-Terpstra trend test for ordered independent groups.
 *
 * Three p-value methods are supported:
 *   - `approximate`: asymptotic normal approximation with optional
 *     continuity correction and tie-aware variance (default).
 *   - `permutation`: Monte Carlo p-value over `nperm` random permutations
 *     of the pooled data.
 *   - `exact`: complete enumeration of distinct permutations of the pooled
 *     data. Cost is `N!` (or `N!/Π(t_h!)` with ties), so the caller is
 *     responsible for keeping `N` small — for `N > ~10` the runtime grows
 *     to seconds or minutes, and beyond `N = 13–14` it becomes infeasible.
 *     Use `approximate` or `permutation` outside that range.
 *
 * References:
 *   - Jonckheere AR (1954). "A distribution-free k-sample test against
 *     ordered alternatives." Biometrika 41, 133-145.
 *   - Terpstra TJ (1952). "The asymptotic normality and consistency of
 *     Kendall's test against trend, when ties are present in one ranking."
 *     Indagationes Mathematicae 14, 327-333.
 *   - Lehmann EL (1975). "Nonparametrics: Statistical Methods Based on
 *     Ranks." Holden-Day, §6.2 (tie-corrected variance).
 */

import {normalCdf} from '../distributions';
import {stripNaN, toFloat64} from '../internal/normalize';
import {shuffleInPlace} from '../internal/random';
import {Alternative, NumericInput} from '../types';

/**
 * Method used to compute the p-value.
 *
 * - `'approximate'`: asymptotic normal approximation; cheap, recommended
 *   default.
 * - `'permutation'`: Monte Carlo over `nperm` shuffles; recommended when
 *   ties are present and `N` is too large for `'exact'`.
 * - `'exact'`: full enumeration; cost ~`N!`, only practical for `N ≲ 10`.
 *   The caller is responsible for not selecting it on large samples.
 */
export type JonckheereMethod = 'approximate' | 'permutation' | 'exact';

/** Optional settings for the Jonckheere-Terpstra test. */
export interface JonckheereSettings {
  /** Direction of the alternative hypothesis. Default `'two-sided'`. */
  alternative?: Alternative;
  /**
   * P-value method. Default `'approximate'`.
   *
   * Note: `'exact'` enumerates `N!` permutations and is only practical for
   * `N ≲ 10` — it does not self-limit, so do not select it on large
   * samples. See {@link JonckheereMethod} for guidance.
   */
  method?: JonckheereMethod;
  /** Number of permutations. Required when `method === 'permutation'`. */
  nperm?: number;
  /**
   * Apply continuity correction in the asymptotic approximation. Default
   * `true`. Has no effect on `'permutation'` / `'exact'` methods.
   */
  continuity?: boolean;
  /**
   * PRNG returning values in `[0, 1)` for the permutation method. Defaults
   * to `Math.random`. Pass `mulberry32(seed)` for reproducible output.
   */
  rng?: () => number;
}

/** Result of the Jonckheere-Terpstra test. */
export interface JonckheereResult {
  /**
   * Z-statistic from the normal approximation (`approximate`,
   * `permutation`) or `null` for the `exact` method.
   */
  statistic: number | null;
  /** Two-sided / one-sided p-value, depending on `alternative`. */
  pValue: number | null;
  /**
   * J statistic under the midrank convention: pairs `(a, b)` with
   * `a ∈ groupᵢ`, `b ∈ groupⱼ`, `i < j`, contribute 1 when `b > a`,
   * ½ when `b == a`, 0 when `b < a`. May be a half-integer when
   * cross-group ties are present.
   */
  jStatistic: number | null;
}

const NULL_RESULT: JonckheereResult = {statistic: null, pValue: null, jStatistic: null};

const MIN_GROUPS = 3;

/**
 * Jonckheere-Terpstra trend test against an ordered alternative.
 *
 * @param groups Ordered list of independent samples (each group is the data
 *               at one level of the ordered factor). Requires at least
 *               `MIN_GROUPS = 3` groups.
 * @param settings Optional configuration: `alternative`, `method`, `nperm`,
 *               `continuity`, `rng`.
 * @returns `{statistic, pValue, jStatistic}` — `statistic` is `null` for
 *          the exact method; `pValue`/`statistic` are `null` when the
 *          variance is non-positive or the pooled sample is too small.
 *
 * @throws `Error` when fewer than `MIN_GROUPS` groups are supplied, when
 *         `method === 'permutation'` and `nperm` is missing, or when an
 *         unknown `method` / `alternative` is provided.
 */
export function jonckheere(
  groups: readonly NumericInput[],
  settings: JonckheereSettings = {},
): JonckheereResult {
  const alternative: Alternative = settings.alternative ?? 'two-sided';
  const method: JonckheereMethod = settings.method ?? 'approximate';
  const continuity = settings.continuity ?? true;

  if (groups.length < MIN_GROUPS) {
    throw new Error(
      `jonckheere: requires at least ${MIN_GROUPS} ordered groups (got ${groups.length}).`,
    );
  }

  const cleaned: Float64Array[] = groups.map((g) => stripNaN(toFloat64(g)));
  const k = cleaned.length;
  const ns = new Int32Array(k);
  let N = 0;
  for (let i = 0; i < k; i++) {
    ns[i] = cleaned[i].length;
    N += ns[i];
  }
  if (N < MIN_GROUPS) return NULL_RESULT;

  const J = jStatistic(cleaned);
  const {mean, variance} = meanVarianceTieAware(cleaned, ns, N);

  switch (method) {
  case 'approximate': {
    if (variance <= 0) return {statistic: null, pValue: null, jStatistic: J};
    const {z, p} = pvalueApproximate(J, mean, variance, alternative, continuity);
    return {statistic: z, pValue: p, jStatistic: J};
  }
  case 'permutation': {
    if (settings.nperm === undefined || settings.nperm <= 0)
      throw new Error('jonckheere: a positive `nperm` is required when method = "permutation".');

    const {z, p} = pvaluePermutation(
      cleaned, ns, N, J, alternative, settings.nperm, settings.rng ?? Math.random,
    );
    return {statistic: z, pValue: p, jStatistic: J};
  }
  case 'exact': {
    const p = pvalueExact(cleaned, ns, N, J, alternative);
    return {statistic: null, pValue: p, jStatistic: J};
  }
  default:
    throw new Error(
      `jonckheere: invalid method ${JSON.stringify(method)}. ` +
        `Must be one of 'approximate', 'permutation', 'exact'.`,
    );
  }
}

/**
 * Compute the J statistic under the midrank convention.
 *
 * For each pair `(a, b)` with `a ∈ groupᵢ`, `b ∈ groupⱼ`, `i < j`:
 *   `b > a` → 1, `b == a` → ½, `b < a` → 0.
 *
 * This matches the convention used by `clinfun::jonckheere.test`,
 * `DescTools::JonckheereTerpstraTest`, `PMCMRplus::jonckheereTest`, and
 * SAS PROC FREQ. It is also the convention under which the closed-form
 * E[J] = (N² − Σnᵢ²)/4 and the tie-corrected Var[J] (Lehmann §6.2) are
 * derived; previously this function used strict `<` (ties → 0), which
 * was inconsistent with the variance formula and biased z under ties.
 */
function jStatistic(groups: readonly Float64Array[]): number {
  let J = 0;
  for (let i = 0; i < groups.length - 1; i++) {
    const a = groups[i];
    if (a.length === 0) continue;
    for (let j = i + 1; j < groups.length; j++) {
      const b = groups[j];
      if (b.length === 0) continue;
      for (let p = 0; p < b.length; p++) {
        const bv = b[p];
        for (let q = 0; q < a.length; q++) {
          if (bv > a[q]) J += 1;
          else if (bv === a[q]) J += 0.5;
        }
      }
    }
  }
  return J;
}

/** J statistic (midrank) computed from a pooled buffer with explicit group offsets. */
function jStatisticFromPooled(pooled: Float64Array, offsets: Int32Array): number {
  let J = 0;
  const k = offsets.length - 1;
  for (let i = 0; i < k - 1; i++) {
    const aStart = offsets[i];
    const aEnd = offsets[i + 1];
    if (aStart === aEnd) continue;
    for (let j = i + 1; j < k; j++) {
      const bStart = offsets[j];
      const bEnd = offsets[j + 1];
      if (bStart === bEnd) continue;
      for (let p = bStart; p < bEnd; p++) {
        const bv = pooled[p];
        for (let q = aStart; q < aEnd; q++) {
          if (bv > pooled[q]) J += 1;
          else if (bv === pooled[q]) J += 0.5;
        }
      }
    }
  }
  return J;
}

/**
 * Mean and variance of the J statistic under H₀ (Lehmann 1975 §6.2).
 *
 * Without ties:
 *   E[J]   = (N² − Σ nᵢ²) / 4
 *   Var[J] = [N(N−1)(2N+5) − Σ nᵢ(nᵢ−1)(2nᵢ+5)] / 72
 *
 * With ties (let `t_h` be the multiplicity of the h-th tied value across the
 * pooled sample), Var[J] gets an additional `−Σ t_h(t_h−1)(2t_h+5)/72` term
 * plus two cross-product corrections (the A and B terms below) that vanish
 * in the no-ties case.
 */
function meanVarianceTieAware(
  groups: readonly Float64Array[], ns: Int32Array, N: number,
): {mean: number; variance: number} {
  let sumNi2 = 0;
  for (let i = 0; i < ns.length; i++) sumNi2 += ns[i] * ns[i];
  const mean = (N * N - sumNi2) / 4;

  const ts = tiedMultiplicities(groups, N);

  const SR = N * (N - 1) * (2 * N + 5);
  let SC = 0;
  for (let i = 0; i < ns.length; i++) {
    const n = ns[i];
    SC += n * (n - 1) * (2 * n + 5);
  }
  let ST = 0;
  for (const t of ts) ST += t * (t - 1) * (2 * t + 5);
  let variance = (SR - SC - ST) / 72;

  if (ts.length > 0 && N >= 2) {
    let sumNi2Pair = 0;
    let sumNi3Triple = 0;
    for (let i = 0; i < ns.length; i++) {
      const n = ns[i];
      sumNi2Pair += n * (n - 1);
      sumNi3Triple += n * (n - 1) * (n - 2);
    }
    let sumT2Pair = 0;
    let sumT3Triple = 0;
    for (const t of ts) {
      sumT2Pair += t * (t - 1);
      sumT3Triple += t * (t - 1) * (t - 2);
    }
    if (N >= 3) {
      const denomA = 36 * N * (N - 1) * (N - 2);
      if (denomA > 0) variance += (sumNi3Triple * sumT3Triple) / denomA;
    }
    const denomB = 8 * N * (N - 1);
    if (denomB > 0) variance += (sumNi2Pair * sumT2Pair) / denomB;
  }
  return {mean, variance};
}

/** Multiplicities `t_h` of values that appear more than once in the pooled sample. */
function tiedMultiplicities(groups: readonly Float64Array[], N: number): number[] {
  const all = new Float64Array(N);
  let off = 0;
  for (const g of groups) {
    all.set(g, off);
    off += g.length;
  }
  const sorted = Array.from(all).sort((a, b) => a - b);
  const ts: number[] = [];
  let i = 0;
  while (i < sorted.length) {
    let j = i + 1;
    while (j < sorted.length && sorted[j] === sorted[i]) j++;
    if (j - i > 1) ts.push(j - i);
    i = j;
  }
  return ts;
}

/** Asymptotic normal approximation, with optional continuity correction. */
function pvalueApproximate(
  J: number, mean: number, variance: number,
  alternative: Alternative, continuity: boolean,
): {z: number; p: number} {
  const delta = continuity ? 0.5 * Math.sign(J - mean) : 0;
  const z = (J - mean - delta) / Math.sqrt(variance);
  return {z, p: tailProbability(z, alternative)};
}

/** Monte Carlo permutation p-value with `nperm` shuffles of the pooled data. */
function pvaluePermutation(
  groups: readonly Float64Array[], ns: Int32Array, N: number,
  observedJ: number, alternative: Alternative, nperm: number, rng: () => number,
): {z: number; p: number} {
  const pooled = new Float64Array(N);
  let off = 0;
  for (const g of groups) {
    pooled.set(g, off);
    off += g.length;
  }
  const offsets = new Int32Array(ns.length + 1);
  for (let i = 0; i < ns.length; i++) offsets[i + 1] = offsets[i] + ns[i];

  const work = new Float64Array(pooled);
  const dist = new Float64Array(nperm);
  for (let i = 0; i < nperm; i++) {
    shuffleInPlace(work, rng);
    dist[i] = jStatisticFromPooled(work, offsets);
  }

  let mean = 0;
  for (let i = 0; i < nperm; i++) mean += dist[i];
  mean /= nperm;
  let ss = 0;
  for (let i = 0; i < nperm; i++) {
    const d = dist[i] - mean;
    ss += d * d;
  }
  const std = nperm > 1 ? Math.sqrt(ss / (nperm - 1)) : 0;
  const z = std === 0 ? 0 : (observedJ - mean) / std;

  // Include the observed statistic in the reference distribution
  let geCount = 1;
  let leCount = 1;
  for (let i = 0; i < nperm; i++) {
    if (dist[i] >= observedJ) geCount++;
    if (dist[i] <= observedJ) leCount++;
  }
  const total = nperm + 1;
  const pInc = geCount / total;
  const pDec = leCount / total;

  let p: number;
  switch (alternative) {
  case 'two-sided':
    p = Math.min(1, 2 * Math.min(pInc, pDec, 0.5));
    break;
  case 'increasing':
    p = pInc;
    break;
  case 'decreasing':
    p = pDec;
    break;
  default:
    throw new Error(`jonckheere: invalid alternative ${JSON.stringify(alternative)}.`);
  }
  return {z, p};
}

/**
 * Exact p-value via complete enumeration of distinct permutations of the
 * pooled sample using Knuth's lexicographic next-permutation algorithm.
 */
function pvalueExact(
  groups: readonly Float64Array[], ns: Int32Array, N: number,
  observedJ: number, alternative: Alternative,
): number {
  const offsets = new Int32Array(ns.length + 1);
  for (let i = 0; i < ns.length; i++) offsets[i + 1] = offsets[i] + ns[i];

  const work = new Float64Array(N);
  let off = 0;
  for (const g of groups) {
    work.set(g, off);
    off += g.length;
  }
  // Sort ascending — the lexicographically smallest permutation
  const tmp = Array.from(work).sort((a, b) => a - b);
  for (let i = 0; i < N; i++) work[i] = tmp[i];

  let total = 0;
  let geCount = 0;
  let leCount = 0;
  do {
    const J = jStatisticFromPooled(work, offsets);
    total++;
    if (J >= observedJ) geCount++;
    if (J <= observedJ) leCount++;
  } while (nextPermutation(work));

  const pInc = geCount / total;
  const pDec = leCount / total;
  switch (alternative) {
  case 'two-sided':
    return Math.min(1, 2 * Math.min(pInc, pDec));
  case 'increasing':
    return pInc;
  case 'decreasing':
    return pDec;
  default:
    throw new Error(`jonckheere: invalid alternative ${JSON.stringify(alternative)}.`);
  }
}

/** Tail probability of the standard normal for the given alternative. */
function tailProbability(z: number, alternative: Alternative): number {
  switch (alternative) {
  case 'two-sided':
    return 2 * (1 - normalCdf(Math.abs(z)));
  case 'increasing':
    return 1 - normalCdf(z);
  case 'decreasing':
    return normalCdf(z);
  default:
    throw new Error(`jonckheere: invalid alternative ${JSON.stringify(alternative)}.`);
  }
}

/**
 * Knuth Algorithm L: in-place next permutation in lexicographic order.
 * Returns `false` when the input is the last (descending) permutation.
 * Naturally skips equivalent permutations when ties are present.
 */
function nextPermutation(a: Float64Array): boolean {
  const n = a.length;
  let k = n - 2;
  while (k >= 0 && a[k] >= a[k + 1]) k--;
  if (k < 0) return false;
  let l = n - 1;
  while (a[l] <= a[k]) l--;
  const tmp = a[k];
  a[k] = a[l];
  a[l] = tmp;
  for (let i = k + 1, j = n - 1; i < j; i++, j--) {
    const t = a[i];
    a[i] = a[j];
    a[j] = t;
  }
  return true;
}
