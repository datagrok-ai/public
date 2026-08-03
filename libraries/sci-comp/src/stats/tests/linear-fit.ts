/**
 * Simple linear regression with a slope confidence interval.
 *
 * `linearFit(x, y, {ciLevel})` fits `y = intercept + slope·x` by OLS (the shared
 * `fitOls` engine) and returns the slope CI from the Student-t quantile. Pairs
 * where either coordinate is NaN are dropped.
 *
 * This is the general primitive the NCA dose-proportionality adapter consumes
 * for the Smith power model (`ln(AUC) ~ ln(Dose)`): the regression math lives
 * here in sci-comp, never in the adapter. It is equally usable for any
 * single-covariate linear fit (e.g. log-linear λz).
 *
 * Degenerate inputs are returned, not thrown:
 *  - `n < 2` or zero x-spread (single distinct x / collinear) → all-NaN slope.
 *  - `n = 2` (df = 0) → slope/intercept defined, but `slopeSe`/`slopeCI` NaN
 *    (no residual df for a standard error).
 */

import {studentTInv} from '../distributions';
import {Matrix} from '../internal/matrix';
import {mean, pairwiseNonNaN, toFloat64} from '../internal/normalize';
import {fitOls} from '../internal/ols';
import {NumericInput} from '../types';

export interface LinearFitOptions {
  /** Two-sided confidence level for the slope CI. Default `0.90` (Smith 1−2α). */
  ciLevel?: number;
}

export interface LinearFitResult {
  slope: number;
  intercept: number;
  /** Standard error of the slope. NaN when `df ≤ 0`. */
  slopeSe: number;
  /** `[lo, hi]` slope CI at `ciLevel`. `[NaN, NaN]` when `df ≤ 0` or singular. */
  slopeCI: [number, number];
  /** Coefficient of determination. NaN when the response has zero variance. */
  rSquared: number;
  /** Residual degrees of freedom, `n − 2`. */
  df: number;
  /** Number of `(x, y)` pairs used after NaN removal. */
  n: number;
}

function degenerate(n: number, df: number): LinearFitResult {
  return {slope: NaN, intercept: NaN, slopeSe: NaN, slopeCI: [NaN, NaN], rSquared: NaN, df, n};
}

/**
 * Fit `y = intercept + slope·x` and return the slope CI at `options.ciLevel`
 * (default 0.90). NaN pairs are dropped before fitting.
 *
 * @throws if `ciLevel` is not strictly inside `(0, 1)`.
 */
export function linearFit(
  x: NumericInput, y: NumericInput, options: LinearFitOptions = {},
): LinearFitResult {
  const ciLevel = options.ciLevel ?? 0.90;
  if (!(ciLevel > 0 && ciLevel < 1))
    throw new Error(`linearFit: ciLevel must be in (0, 1), got ${ciLevel}.`);

  const {x: xf, y: yf} = pairwiseNonNaN(toFloat64(x), toFloat64(y));
  const n = xf.length;
  const df = n - 2;

  // A defined slope needs ≥2 points and non-degenerate x-spread.
  const xbar = mean(xf);
  let sxx = 0;
  for (let i = 0; i < n; i++) {const d = xf[i] - xbar; sxx += d * d;}
  if (n < 2 || !(sxx > 0)) return degenerate(n, df);

  const X: Matrix = [];
  const yArr: number[] = [];
  for (let i = 0; i < n; i++) {X.push([1, xf[i]]); yArr.push(yf[i]);}
  const fit = fitOls(X, yArr);
  const intercept = fit.beta[0];
  const slope = fit.beta[1];

  const ybar = mean(yf);
  let tss = 0;
  for (let i = 0; i < n; i++) {const d = yf[i] - ybar; tss += d * d;}
  const rSquared = tss > 0 ? 1 - fit.rss / tss : NaN;

  // df ≤ 0 (n = 2): slope is defined, but there is no residual df for an SE/CI.
  if (df <= 0)
    return {slope, intercept, slopeSe: NaN, slopeCI: [NaN, NaN], rSquared, df, n};

  const slopeSe = Math.sqrt(Math.max(0, fit.vcov[1][1]));
  const tCrit = studentTInv(1 - (1 - ciLevel) / 2, df);
  const half = tCrit * slopeSe;
  return {slope, intercept, slopeSe, slopeCI: [slope - half, slope + half], rSquared, df, n};
}
