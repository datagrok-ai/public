/**
 * One-way ANOVA — fixed-effects F-test for a difference in group means.
 *
 * `oneWayAnova(values, groups)` partitions the total sum of squares into
 * between-group and within-group components and reports the F-test. It is the
 * no-covariate complement to `runAncova` (which requires a covariate), used by
 * the NCA dose-proportionality adapter for the *secondary* dose-normalized-AUC
 * comparison across dose groups (FR-412). NaN responses are dropped pairwise
 * with their group label.
 *
 * Note (caller's responsibility): a non-significant one-way ANOVA is *absence of
 * evidence*, not evidence of proportionality — the adapter frames it strictly
 * secondary to the equivalence-based power-model CI.
 */

import {fCdf} from '../distributions';
import {toFloat64} from '../internal/normalize';
import {NumericInput} from '../types';

export interface OneWayAnovaResult {
  fStatistic: number;
  /** Between-group degrees of freedom, `k − 1`. */
  dfBetween: number;
  /** Within-group degrees of freedom, `N − k`. */
  dfWithin: number;
  pValue: number;
  /** Number of distinct groups after NaN removal. */
  groups: number;
  /** Total number of observations used after NaN removal. */
  n: number;
}

/**
 * One-way ANOVA F-test of `values ~ C(groups)`.
 *
 * Returns `null` when there are insufficient data for the test: fewer than two
 * groups, `N ≤ k` (no within-group df), or zero within-group variation.
 *
 * @throws if `values` and `groups` differ in length.
 */
export function oneWayAnova(
  values: NumericInput, groups: NumericInput,
): OneWayAnovaResult | null {
  const vRaw = toFloat64(values);
  const gRaw = toFloat64(groups);
  if (vRaw.length !== gRaw.length)
    throw new Error('oneWayAnova: values and groups length mismatch.');

  // Pairwise NaN removal (a NaN group label is also dropped).
  const byGroup = new Map<number, number[]>();
  let n = 0;
  let grandSum = 0;
  for (let i = 0; i < vRaw.length; i++) {
    const v = vRaw[i];
    const g = gRaw[i];
    if (Number.isNaN(v) || Number.isNaN(g)) continue;
    let arr = byGroup.get(g);
    if (arr === undefined) {arr = []; byGroup.set(g, arr);}
    arr.push(v);
    grandSum += v;
    n++;
  }

  const k = byGroup.size;
  if (k < 2 || n <= k) return null;

  const grandMean = grandSum / n;
  let ssBetween = 0;
  let ssWithin = 0;
  for (const arr of byGroup.values()) {
    let sum = 0;
    for (const v of arr) sum += v;
    const gMean = sum / arr.length;
    ssBetween += arr.length * (gMean - grandMean) ** 2;
    for (const v of arr) ssWithin += (v - gMean) ** 2;
  }

  const dfBetween = k - 1;
  const dfWithin = n - k;
  if (!(ssWithin > 0)) return null; // zero within-group variation → F undefined

  const msBetween = ssBetween / dfBetween;
  const msWithin = ssWithin / dfWithin;
  const fStatistic = msBetween / msWithin;
  const pValue = 1 - fCdf(fStatistic, dfBetween, dfWithin);

  return {fStatistic, dfBetween, dfWithin, pValue, groups: k, n};
}
