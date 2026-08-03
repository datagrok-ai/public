/**
 * Mann-Whitney U test (two-sample, non-parametric).
 *
 * Matches scipy's `mannwhitneyu(method='auto')`: exact distribution for
 * `min(n1, n2) ≤ 8` without ties, asymptotic with continuity correction and
 * tie adjustment otherwise.
 */

import {normalSf} from '../distributions';
import {stripNaN, toFloat64} from '../internal/normalize';
import {rankdata, tieCorrection} from '../internal/rank';
import {Alternative, NumericInput, TestResult} from '../types';

/**
 * Mann-Whitney U test.
 *
 * The returned `statistic` is U for the first sample (matching scipy).
 * Returns `null` fields when either sample has no finite values.
 */
export function mannWhitneyU(
  x: NumericInput,
  y: NumericInput,
  alternative: Alternative = 'two-sided',
): TestResult {
  const a = stripNaN(toFloat64(x));
  const b = stripNaN(toFloat64(y));
  const n1 = a.length;
  const n2 = b.length;
  if (n1 < 1 || n2 < 1)
    return {statistic: null, pValue: null};

  const combined = new Float64Array(n1 + n2);
  for (let i = 0; i < n1; i++) combined[i] = a[i];
  for (let i = 0; i < n2; i++) combined[n1 + i] = b[i];
  const ranks = rankdata(combined);

  let R1 = 0;
  for (let i = 0; i < n1; i++) R1 += ranks[i];
  const U1 = R1 - n1 * (n1 + 1) / 2;

  const T = tieCorrection(combined);
  const hasTies = T > 0;

  const p = (!hasTies && Math.min(n1, n2) <= 8) ?
    exactP(U1, n1, n2, alternative) :
    asymptoticP(U1, n1, n2, T, alternative);

  return {statistic: U1, pValue: p};
}

// ── Exact distribution (DP over the discrete null) ────────────────

function exactP(u: number, n1: number, n2: number, alt: Alternative): number {
  const pmf = exactPmf(n1, n2);
  const max = n1 * n2;

  if (alt === 'increasing') {
    let p = 0;
    for (let k = Math.ceil(u); k <= max; k++) p += pmf[k];
    return p;
  }
  if (alt === 'decreasing') {
    let p = 0;
    for (let k = 0; k <= Math.floor(u); k++) p += pmf[k];
    return p;
  }
  // two-sided
  const mid = max / 2;
  if (u < mid) {
    let plow = 0;
    for (let k = 0; k <= Math.floor(u); k++) plow += pmf[k];
    return Math.min(2 * plow, 1.0);
  }
  if (u > mid) {
    let phigh = 0;
    for (let k = Math.ceil(u); k <= max; k++) phigh += pmf[k];
    return Math.min(2 * phigh, 1.0);
  }
  return 1.0;
}

/**
 * Exact PMF of U₁ under H₀ for sample sizes (n1, n2) with no ties.
 *
 * DP: place ranks one by one. When the next rank goes to group 1, it
 * contributes `j` (current count of group-2 ranks already placed) to U₁;
 * when it goes to group 2, it contributes 0.
 */
function exactPmf(n1: number, n2: number): Float64Array {
  const max = n1 * n2;
  const cnt: Float64Array[][] = [];
  for (let i = 0; i <= n1; i++) {
    const row: Float64Array[] = [];
    for (let j = 0; j <= n2; j++)
      row.push(new Float64Array(max + 1));
    cnt.push(row);
  }
  cnt[0][0][0] = 1;
  for (let i = 0; i <= n1; i++) {
    for (let j = 0; j <= n2; j++) {
      for (let u = 0; u <= max; u++) {
        const c = cnt[i][j][u];
        if (c === 0) continue;
        if (i < n1 && u + j <= max) cnt[i + 1][j][u + j] += c;
        if (j < n2) cnt[i][j + 1][u] += c;
      }
    }
  }
  let total = 0;
  for (let u = 0; u <= max; u++) total += cnt[n1][n2][u];
  const pmf = new Float64Array(max + 1);
  for (let u = 0; u <= max; u++) pmf[u] = cnt[n1][n2][u] / total;
  return pmf;
}

// ── Asymptotic with continuity correction and tie adjustment ──────

function asymptoticP(
  u: number, n1: number, n2: number, T: number, alt: Alternative,
): number {
  const n = n1 + n2;
  const mu = n1 * n2 / 2;
  const varU = n1 * n2 / 12 * ((n + 1) - T / (n * (n - 1)));
  if (varU <= 0) return 1.0;
  const sigma = Math.sqrt(varU);

  if (alt === 'two-sided') {
    const cc = u > mu ? -0.5 : (u < mu ? 0.5 : 0);
    const z = (u - mu + cc) / sigma;
    return 2 * normalSf(Math.abs(z));
  }
  if (alt === 'increasing') {
    const z = (u - mu - 0.5) / sigma;
    return normalSf(z);
  }
  // decreasing
  const z = (u - mu + 0.5) / sigma;
  return 1 - normalSf(z);
}
