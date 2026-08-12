/**
 * Area-under-the-curve (AUC) — trapezoidal integration kernels.
 *
 * Three integration schemes matching the PKNCA `auc.method` options:
 * - `linear`              — straight trapezoid;
 * - `log-linear`          — exponential decay assumption between points;
 * - `linear-up/log-down`  — linear on ascending intervals, log-linear on
 *                            descending intervals (PKNCA default).
 *
 * Each scheme comes in a naive and a Neumaier-compensated summation variant.
 *
 * Shared signature: `time` is sorted ascending, `conc` has the same length,
 * and integration runs over the inclusive index range `[startIdx, endIdx]`,
 * summing the contribution of each pair `(i, i+1)` — so the result is `0`
 * when `endIdx <= startIdx`. Callers handle BLQ pre-processing themselves
 * (see {@link applyBlqStrategy}).
 */

/**
 * Linear trapezoidal AUC.
 *
 * Formula per interval: `(t_{i+1} − t_i) · (c_i + c_{i+1}) / 2`.
 */
export function aucLinearNaive(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let auc = 0;
  for (let i = startIdx; i < endIdx; i++) {
    const dt = time[i + 1] - time[i];
    auc += dt * (conc[i] + conc[i + 1]) / 2;
  }
  return auc;
}

/**
 * Log-linear trapezoidal AUC.
 *
 * Formula per interval: `(t_{i+1} − t_i) · (c_i − c_{i+1}) / ln(c_i / c_{i+1})`,
 * which integrates an exponential interpolant `c(t) = c_i · exp(−k·(t−t_i))`
 * exactly. Falls back to the linear formula when the log-linear formula is
 * undefined or numerically unstable:
 * - either concentration ≤ 0 (log of non-positive),
 * - or `c_i == c_{i+1}` (division by `ln(1) = 0`).
 */
export function aucLogLinearNaive(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let auc = 0;
  for (let i = startIdx; i < endIdx; i++) {
    const dt = time[i + 1] - time[i];
    const c1 = conc[i];
    const c2 = conc[i + 1];
    if (c1 > 0 && c2 > 0 && c1 !== c2)
      auc += dt * (c1 - c2) / Math.log(c1 / c2);
    else
      auc += dt * (c1 + c2) / 2;
  }
  return auc;
}

/**
 * Linear-up / log-down trapezoidal AUC (PKNCA default `lin up/log down`).
 *
 * Per interval:
 * - linear if `c_{i+1} ≥ c_i` (ascending or flat);
 * - log-linear if `c_{i+1} < c_i` AND both concentrations are positive;
 * - linear otherwise (descending with a zero or negative).
 *
 * Rationale: drug absorption phases are well-approximated by linear
 * interpolation, while elimination decays follow exponential kinetics and
 * are better captured by the log-linear formula.
 */
export function aucLinearUpLogDownNaive(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let auc = 0;
  for (let i = startIdx; i < endIdx; i++) {
    const dt = time[i + 1] - time[i];
    const c1 = conc[i];
    const c2 = conc[i + 1];
    if (c2 < c1 && c1 > 0 && c2 > 0)
      auc += dt * (c1 - c2) / Math.log(c1 / c2);
    else
      auc += dt * (c1 + c2) / 2;
  }
  return auc;
}

/**
 * Extrapolate AUC from the last observed time to infinity, assuming
 * exponential decay at rate `lambdaZ` from `cLast` onwards.
 *
 * Formula: `cLast / lambdaZ`.
 *
 * Returns `+Infinity` when `lambdaZ == 0` and `cLast > 0`, `NaN` when both
 * are zero (`0/0`), and `0` when `cLast == 0`. Callers should guard
 * upstream — `lambdaZ <= 0` is a sign that the terminal slope was not
 * estimable and AUCinf should not be reported.
 */
export function aucExtrapolateToInfinity(cLast: number, lambdaZ: number): number {
  return cLast / lambdaZ;
}

// Compensated (Neumaier) variants. Identical to the naive ones on
// well-conditioned inputs (< 1e-15 relative); they earn their keep when
// contributions span many orders of magnitude — a Cmax-dominated spike
// alongside a long tail of small terminal terms — where naive Float64
// summation accumulates catastrophic cancellation error.
//
// Neumaier, A. (1974). "Rundungsfehleranalyse einiger Verfahren zur
// Summation endlicher Summen." ZAMM 54: 39–51.

/**
 * Neumaier-compensated summation of an iterable of numbers.
 *
 * Drop-in replacement for `[...values].reduce((s, v) => s + v, 0)` that
 * tracks a running compensation term so that round-off errors of opposite
 * sign cancel out instead of accumulating. Single-pass, O(n) time, O(1)
 * extra memory. `values` is iterated once.
 */
export function neumaierSum(values: Iterable<number>): number {
  let sum = 0;
  let c = 0;
  for (const v of values) {
    const t = sum + v;
    if (Math.abs(sum) >= Math.abs(v))
      c += (sum - t) + v;
    else
      c += (v - t) + sum;
    sum = t;
  }
  return sum + c;
}

/**
 * Linear trapezoidal AUC with Neumaier-compensated summation.
 *
 * Same semantics as {@link aucLinearNaive}; differs only in the summation
 * algorithm. Use when the profile contains contributions spanning many
 * orders of magnitude.
 */
export function aucLinearCompensated(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let sum = 0;
  let c = 0;
  for (let i = startIdx; i < endIdx; i++) {
    const v = (time[i + 1] - time[i]) * (conc[i] + conc[i + 1]) / 2;
    const t = sum + v;
    if (Math.abs(sum) >= Math.abs(v))
      c += (sum - t) + v;
    else
      c += (v - t) + sum;
    sum = t;
  }
  return sum + c;
}

/**
 * Log-linear trapezoidal AUC with Neumaier-compensated summation.
 *
 * Same fallback rules as {@link aucLogLinearNaive} (linear when either
 * concentration is non-positive or `c1 == c2`).
 */
export function aucLogLinearCompensated(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let sum = 0;
  let c = 0;
  for (let i = startIdx; i < endIdx; i++) {
    const dt = time[i + 1] - time[i];
    const c1 = conc[i];
    const c2 = conc[i + 1];
    const v = (c1 > 0 && c2 > 0 && c1 !== c2) ?
      dt * (c1 - c2) / Math.log(c1 / c2) :
      dt * (c1 + c2) / 2;
    const t = sum + v;
    if (Math.abs(sum) >= Math.abs(v))
      c += (sum - t) + v;
    else
      c += (v - t) + sum;
    sum = t;
  }
  return sum + c;
}

/**
 * Linear-up / log-down trapezoidal AUC with Neumaier-compensated summation.
 *
 * Same semantics as {@link aucLinearUpLogDownNaive}; differs only in the
 * summation algorithm.
 */
export function aucLinearUpLogDownCompensated(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let sum = 0;
  let c = 0;
  for (let i = startIdx; i < endIdx; i++) {
    const dt = time[i + 1] - time[i];
    const c1 = conc[i];
    const c2 = conc[i + 1];
    const v = (c2 < c1 && c1 > 0 && c2 > 0) ?
      dt * (c1 - c2) / Math.log(c1 / c2) :
      dt * (c1 + c2) / 2;
    const t = sum + v;
    if (Math.abs(sum) >= Math.abs(v))
      c += (sum - t) + v;
    else
      c += (v - t) + sum;
    sum = t;
  }
  return sum + c;
}
