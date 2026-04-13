/**
 * Single-pass OLS linear regression against t = 0, 1, ..., n-1.
 *
 * Matches `scipy.stats.linregress` behaviour.
 */

import {NumericArray} from './types';
import {twoTailPvalue} from './stats-utils';

/** Result of ordinary least-squares linear regression. */
export interface LinearTrendResult {
  /** Slope of the fitted line. */
  slope: number;
  /** Y-intercept of the fitted line. */
  intercept: number;
  /** Pearson correlation coefficient. */
  rvalue: number;
  /** Two-tailed p-value for a hypothesis test with H₀: slope = 0. */
  pvalue: number;
  /** Standard error of the slope estimate. */
  stderr: number;
}

/**
 * Compute OLS linear trend for a time-series slice.
 *
 * Regresses `x[0..n-1]` against `t = 0, 1, ..., n-1`. Uses analytic
 * sums for the time index (no loop needed for sumT, sumT²).
 *
 * Edge cases:
 * - `n < 2`: returns all zeros with `pvalue = 1.0`.
 * - `n == 2`: `stderr = 0`, `pvalue = 1.0` (zero residual degrees of freedom).
 * - Constant series (`SS_xx == 0`): `rvalue = 0`, `slope = 0`, `pvalue = 1.0`.
 */
export function linearTrend(x: NumericArray, n: number): LinearTrendResult {
  if (n < 2)
    return {slope: 0, intercept: 0, rvalue: 0, pvalue: 1.0, stderr: 0};

  // Analytic sums for t = 0..n-1
  const sumT = n * (n - 1) / 2;
  const sumT2 = n * (n - 1) * (2 * n - 1) / 6;

  let sumX = 0;
  let sumX2 = 0;
  let sumTX = 0;
  for (let i = 0; i < n; i++) {
    const v = x[i];
    sumX += v;
    sumX2 += v * v;
    sumTX += i * v;
  }

  const meanT = sumT / n;
  const meanX = sumX / n;

  const ssTT = sumT2 - n * meanT * meanT;
  const ssXX = sumX2 - n * meanX * meanX;
  const ssTX = sumTX - n * meanT * meanX;

  const slope = ssTX / ssTT;
  const intercept = meanX - slope * meanT;

  let rvalue: number;
  if (ssXX <= 0)
    rvalue = 0;
  else
    rvalue = ssTX / Math.sqrt(ssTT * ssXX);

  const rss = Math.max(0, ssXX - slope * ssTX);

  let stderr: number;
  if (n > 2 && ssTT > 0)
    stderr = Math.sqrt(rss / (n - 2)) / Math.sqrt(ssTT);
  else
    stderr = 0;

  let pvalue: number;
  if (stderr > 0) {
    const tStat = slope / stderr;
    pvalue = twoTailPvalue(tStat, n - 2);
  } else
    pvalue = 1.0;

  return {slope, intercept, rvalue, pvalue, stderr};
}
