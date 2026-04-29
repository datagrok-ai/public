/**
 * Internal helpers for input normalisation and NaN handling.
 */

import {NumericInput} from '../types';

/** Coerce any `NumericInput` to a `Float64Array` (copying when necessary). */
export function toFloat64(x: NumericInput): Float64Array {
  if (x instanceof Float64Array)
    return x;
  return Float64Array.from(x as ArrayLike<number>);
}

/** Strip NaN values from a typed array, returning a new Float64Array. */
export function stripNaN(x: Float64Array): Float64Array {
  let n = 0;
  for (let i = 0; i < x.length; i++)
    if (!Number.isNaN(x[i])) n++;
  if (n === x.length) return x;
  const out = new Float64Array(n);
  let j = 0;
  for (let i = 0; i < x.length; i++)
    if (!Number.isNaN(x[i])) out[j++] = x[i];
  return out;
}

/** Pairwise NaN removal: keeps positions where both arrays have finite values. */
export function pairwiseNonNaN(
  x: Float64Array, y: Float64Array,
): {x: Float64Array; y: Float64Array} {
  if (x.length !== y.length)
    throw new Error(`pairwiseNonNaN: length mismatch (${x.length} vs ${y.length})`);
  let n = 0;
  for (let i = 0; i < x.length; i++)
    if (!Number.isNaN(x[i]) && !Number.isNaN(y[i])) n++;
  const xo = new Float64Array(n);
  const yo = new Float64Array(n);
  let j = 0;
  for (let i = 0; i < x.length; i++) {
    if (!Number.isNaN(x[i]) && !Number.isNaN(y[i])) {
      xo[j] = x[i];
      yo[j] = y[i];
      j++;
    }
  }
  return {x: xo, y: yo};
}

/** Sample mean. Returns NaN for empty input. */
export function mean(x: Float64Array): number {
  if (x.length === 0) return NaN;
  let s = 0;
  for (let i = 0; i < x.length; i++) s += x[i];
  return s / x.length;
}

/** Sample variance with `ddof=1` (Bessel-corrected). Returns NaN for n<2. */
export function variance(x: Float64Array): number {
  const n = x.length;
  if (n < 2) return NaN;
  const m = mean(x);
  let s = 0;
  for (let i = 0; i < n; i++) {
    const d = x[i] - m;
    s += d * d;
  }
  return s / (n - 1);
}

/** Sample standard deviation with `ddof=1`. */
export function stddev(x: Float64Array): number {
  return Math.sqrt(variance(x));
}

/** Sum of an array. */
export function sum(x: Float64Array | readonly number[]): number {
  let s = 0;
  for (let i = 0; i < x.length; i++) s += x[i];
  return s;
}
