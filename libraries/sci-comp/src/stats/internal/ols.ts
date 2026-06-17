/**
 * Ordinary least squares via the normal equations — shared internal engine.
 *
 * Lifted verbatim from `tests/ancova.ts` so that ANCOVA and the public
 * `linearFit` consume a single implementation rather than duplicating the
 * normal-equations math. `X` must already include the intercept column.
 */

import {inverse, matmul, Matrix, matvec, transpose} from './matrix';

export interface OlsFit {
  /** Coefficient vector, aligned to the columns of `X`. */
  beta: number[];
  /** Residual sum of squares. */
  rss: number;
  /** Residual degrees of freedom, `n − p`. */
  df: number;
  /** Coefficient variance-covariance matrix (`(XᵀX)⁻¹ · mse`). */
  vcov: Matrix;
}

/** OLS fit of `y ~ X` where `X` includes the intercept column. */
export function fitOls(X: Matrix, y: number[]): OlsFit {
  const n = X.length;
  const p = X[0].length;
  const Xt = transpose(X);
  const XtX = matmul(Xt, X);
  const Xty = matvec(Xt, y);
  const XtXinv = inverse(XtX);
  const beta = matvec(XtXinv, Xty);
  const yPred = matvec(X, beta);
  let rss = 0;
  for (let i = 0; i < n; i++) {
    const r = y[i] - yPred[i];
    rss += r * r;
  }
  const df = n - p;
  const mse = df > 0 ? rss / df : Infinity;
  const vcov: Matrix = [];
  for (let i = 0; i < p; i++) {
    const row = new Array<number>(p);
    for (let j = 0; j < p; j++) row[j] = XtXinv[i][j] * mse;
    vcov.push(row);
  }
  return {beta, rss, df, vcov};
}
