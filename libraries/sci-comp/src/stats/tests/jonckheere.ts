/**
 * Jonckheere-Terpstra trend test for ordered independent groups.
 */

import {normalCdf} from '../distributions';
import {stripNaN, toFloat64} from '../internal/normalize';
import {NumericInput, TestResult} from '../types';

/**
 * Jonckheere-Terpstra two-sided trend test using the normal approximation.
 *
 * J = Σ_{i<j} U_ij where U_ij counts pairs (a, b) with a ∈ groupᵢ, b ∈ groupⱼ
 * such that b > a (ties contribute 0.5).
 *
 * Returns `null` when k < 2, total N < 4, or the variance is non-positive.
 */
export function jonckheere(groups: readonly NumericInput[]): TestResult {
  const cleaned: Float64Array[] = groups.map((g) => stripNaN(toFloat64(g)));
  const k = cleaned.length;
  const ns = cleaned.map((g) => g.length);
  let N = 0;
  for (const n of ns) N += n;
  if (k < 2 || N < 4) return {statistic: null, pValue: null};

  let J = 0;
  for (let i = 0; i < k; i++) {
    if (ns[i] === 0) continue;
    for (let j = i + 1; j < k; j++) {
      if (ns[j] === 0) continue;
      const a = cleaned[i];
      const b = cleaned[j];
      // Outer comparison: count b > a (strictly) + 0.5 × ties
      let strict = 0;
      let ties = 0;
      for (let p = 0; p < b.length; p++) {
        const bv = b[p];
        for (let q = 0; q < a.length; q++) {
          const d = bv - a[q];
          if (d > 0) strict++;
          else if (d === 0) ties++;
        }
      }
      J += strict + 0.5 * ties;
    }
  }

  let sumSq = 0;
  let sumWeighted = 0;
  for (const n of ns) {
    sumSq += n * n;
    sumWeighted += n * n * (2 * n + 3);
  }
  const E = (N * N - sumSq) / 4;
  const Var = (N * N * (2 * N + 3) - sumWeighted) / 72;
  if (Var <= 0) return {statistic: null, pValue: null};

  const Z = (J - E) / Math.sqrt(Var);
  const p = 2 * (1 - normalCdf(Math.abs(Z)));
  return {statistic: Z, pValue: p};
}
