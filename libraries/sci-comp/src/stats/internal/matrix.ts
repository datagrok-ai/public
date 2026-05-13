/**
 * Minimal matrix utilities for OLS in ANCOVA.
 *
 * Matrix inverse uses `jStat.inv`; everything else is hand-written for
 * the (typically small, p ≤ 6) dimensions encountered.
 */

// @ts-ignore: no types
import * as jStat from 'jstat';

const J: any = jStat;

export type Matrix = number[][];

/** Matrix transpose. */
export function transpose(A: Matrix): Matrix {
  const m = A.length;
  const n = A[0].length;
  const T: Matrix = [];
  for (let i = 0; i < n; i++) {
    const row = new Array(m);
    for (let j = 0; j < m; j++) row[j] = A[j][i];
    T.push(row);
  }
  return T;
}

/** Matrix-matrix product. */
export function matmul(A: Matrix, B: Matrix): Matrix {
  const m = A.length;
  const k = A[0].length;
  const n = B[0].length;
  const C: Matrix = [];
  for (let i = 0; i < m; i++) {
    const row = new Array(n).fill(0);
    for (let l = 0; l < k; l++) {
      const a = A[i][l];
      if (a === 0) continue;
      for (let j = 0; j < n; j++) row[j] += a * B[l][j];
    }
    C.push(row);
  }
  return C;
}

/** Matrix-vector product. */
export function matvec(A: Matrix, x: number[]): number[] {
  const m = A.length;
  const n = A[0].length;
  const y = new Array<number>(m);
  for (let i = 0; i < m; i++) {
    let s = 0;
    for (let j = 0; j < n; j++) s += A[i][j] * x[j];
    y[i] = s;
  }
  return y;
}

/** Matrix inverse via jstat. */
export function inverse(A: Matrix): Matrix {
  return J.inv(A);
}

/** Inner product xᵀ A x. */
export function quadForm(x: number[], A: Matrix): number {
  let s = 0;
  for (let i = 0; i < x.length; i++) {
    let inner = 0;
    for (let j = 0; j < x.length; j++) inner += A[i][j] * x[j];
    s += x[i] * inner;
  }
  return s;
}
