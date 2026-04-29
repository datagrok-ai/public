/**
 * Typed wrappers around jstat's distributions and special functions.
 *
 * jstat ships without TypeScript types, so this module is the only place that
 * uses `@ts-ignore` / `any` against it. The rest of the stats module imports
 * from here.
 */

// @ts-ignore: no types
import * as jStat from 'jstat';

const J: any = jStat;

// ── Standard normal ──────────────────────────────────────────────

/** Standard normal CDF: Φ(x). */
export function normalCdf(x: number): number {
  return J.normal.cdf(x, 0, 1);
}

/** Standard normal survival function: 1 − Φ(x). */
export function normalSf(x: number): number {
  return 1 - J.normal.cdf(x, 0, 1);
}

/** Standard normal inverse CDF (quantile): Φ⁻¹(p). */
export function normalInv(p: number): number {
  return J.normal.inv(p, 0, 1);
}

// ── Student's t ──────────────────────────────────────────────────

/** Student's t CDF at `x` with `df` degrees of freedom. */
export function studentTCdf(x: number, df: number): number {
  return J.studentt.cdf(x, df);
}

/** Student's t inverse CDF (quantile): t-value at probability `p`. */
export function studentTInv(p: number, df: number): number {
  return J.studentt.inv(p, df);
}

/** Two-tailed p-value for a t-statistic. */
export function studentTTwoTail(t: number, df: number): number {
  return 2 * (1 - J.studentt.cdf(Math.abs(t), df));
}

// ── F (central) ──────────────────────────────────────────────────

/** Central F CDF at `x` with `df1`, `df2` degrees of freedom. */
export function fCdf(x: number, df1: number, df2: number): number {
  return J.centralF.cdf(x, df1, df2);
}

// ── Chi-square ───────────────────────────────────────────────────

/** Chi-square CDF at `x` with `df` degrees of freedom. */
export function chi2Cdf(x: number, df: number): number {
  return J.chisquare.cdf(x, df);
}

/** Chi-square survival function: 1 − F(x). */
export function chi2Sf(x: number, df: number): number {
  return 1 - J.chisquare.cdf(x, df);
}

// ── Hypergeometric ───────────────────────────────────────────────

/**
 * Hypergeometric PMF: P(X = k) where X is the count of successes in a draw of
 * `size` items from a population of `N` containing `K` successes.
 *
 * jstat's signature: pdf(k, N, K, n).
 */
export function hypgeomPmf(k: number, N: number, K: number, n: number): number {
  return J.hypgeom.pdf(k, N, K, n);
}

// ── Special functions ────────────────────────────────────────────

/** Error function. */
export function erf(x: number): number {
  return J.erf(x);
}

/** Log-gamma: ln Γ(x). */
export function gammaln(x: number): number {
  return J.gammaln(x);
}

/** Gamma function: Γ(x). */
export function gammafn(x: number): number {
  return J.gammafn(x);
}

/** Regularised incomplete beta: I_x(a, b). */
export function ibeta(x: number, a: number, b: number): number {
  return J.ibeta(x, a, b);
}

/** Inverse regularised incomplete beta: I⁻¹_p(a, b). */
export function ibetainv(p: number, a: number, b: number): number {
  return J.ibetainv(p, a, b);
}

/** Log of binomial coefficient: ln C(n, k). */
export function combinationln(n: number, k: number): number {
  return J.combinationln(n, k);
}

/** Log-factorial: ln n!. */
export function factorialln(n: number): number {
  return J.factorialln(n);
}
