/**
 * Area-under-the-first-moment-curve (AUMC) — trapezoidal moment kernels.
 *
 * AUMC integrates `t·C(t)` (the first moment of the concentration curve) and
 * is the basis for the mean residence time (MRT) and the steady-state volume
 * (Vss). It is **not** the AUC of a `t·C` array: the log-linear interval has a
 * distinct closed form derived from the exponential interpolant. The three
 * schemes mirror the {@link module:auc} options and are tied to the same
 * `rules.aucMethod` choice:
 * - `linear`              — trapezoid on the moment curve;
 * - `log-linear`          — exact moment of an exponential interpolant;
 * - `linear-up/log-down`  — linear on ascending intervals, log-linear on
 *                            descending ones (PKNCA default).
 *
 * Closed forms over an interval `[t₁, t₂]` (Gabrielsson & Weiss §2.8.1;
 * Yamaoka et al. 1978):
 * - linear:     `(t₂ − t₁)·(t₁·c₁ + t₂·c₂) / 2`
 * - log-linear: `(t₁·c₁ − t₂·c₂)/k + (c₁ − c₂)/k²`,
 *               where `k = ln(c₁/c₂)/(t₂ − t₁)` is the per-interval elimination
 *               rate (`k > 0` when `c₁ > c₂ > 0`) — the **same orientation**
 *               `auc.ts` uses inline (`Math.log(c1/c2)`), so the second term
 *               carries its stated sign. A flipped `ln(c₂/c₁)` would silently
 *               sign-flip the `1/k²` term and fail PKNCA parity on descending
 *               segments.
 *
 * The log-linear moment falls back to the linear formula when the exponential
 * interpolant is undefined or numerically unstable (`c₁ ≤ 0`, `c₂ ≤ 0`, or
 * `c₁ = c₂`), matching the `auc.ts` fallback rules.
 *
 * Compensated (Neumaier) variants live alongside the naive ones: the
 * `(c₁−c₂)/k²` term is sensitive to catastrophic cancellation as `k → 0`
 * (near-flat segments) under naive Float64 summation — the same regime that
 * motivated the compensated AUC kernels.
 */

/**
 * Trapezoidal moment contribution over one interval (linear interpolant).
 *
 * `(t₂ − t₁)·(t₁·c₁ + t₂·c₂)/2` — the trapezoidal rule applied to the moment
 * curve `t·C`, used both for the `linear` method and as the fallback for the
 * log-linear method.
 */
function linMoment(t1: number, c1: number, t2: number, c2: number): number {
  return (t2 - t1) * (t1 * c1 + t2 * c2) / 2;
}

/**
 * Moment contribution over one interval assuming exponential decay
 * `C(t) = c₁·exp(−k·(t−t₁))` between the points.
 *
 * `(t₁·c₁ − t₂·c₂)/k + (c₁ − c₂)/k²`, `k = ln(c₁/c₂)/(t₂ − t₁)`.
 *
 * Caller guarantees `c₁ > 0`, `c₂ > 0`, `c₁ ≠ c₂` (otherwise use
 * {@link linMoment}); under those conditions `k` is finite and non-zero.
 */
function logMoment(t1: number, c1: number, t2: number, c2: number): number {
  const k = Math.log(c1 / c2) / (t2 - t1);
  return (t1 * c1 - t2 * c2) / k + (c1 - c2) / (k * k);
}

/**
 * Linear trapezoidal AUMC over the inclusive index range `[startIdx, endIdx]`.
 *
 * @param time - Time vector, sorted ascending.
 * @param conc - Concentration vector, same length as `time`.
 * @param startIdx - Inclusive start index.
 * @param endIdx - Inclusive end index. Returns `0` when `endIdx <= startIdx`.
 * @returns Sum of the per-interval moment contributions.
 */
export function aumcLinearNaive(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let aumc = 0;
  for (let i = startIdx; i < endIdx; i++)
    aumc += linMoment(time[i], conc[i], time[i + 1], conc[i + 1]);
  return aumc;
}

/**
 * Log-linear trapezoidal AUMC. Falls back to the linear moment formula when
 * either concentration is non-positive or `c₁ == c₂` (mirrors
 * {@link aumcLogLinearNaive}'s AUC counterpart).
 */
export function aumcLogLinearNaive(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let aumc = 0;
  for (let i = startIdx; i < endIdx; i++) {
    const t1 = time[i]; const c1 = conc[i];
    const t2 = time[i + 1]; const c2 = conc[i + 1];
    aumc += (c1 > 0 && c2 > 0 && c1 !== c2) ?
      logMoment(t1, c1, t2, c2) : linMoment(t1, c1, t2, c2);
  }
  return aumc;
}

/**
 * Linear-up / log-down trapezoidal AUMC (PKNCA default `lin up/log down`).
 * Linear moment on ascending/flat intervals, log-linear moment on descending
 * intervals where both concentrations are positive.
 */
export function aumcLinearUpLogDownNaive(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let aumc = 0;
  for (let i = startIdx; i < endIdx; i++) {
    const t1 = time[i]; const c1 = conc[i];
    const t2 = time[i + 1]; const c2 = conc[i + 1];
    aumc += (c2 < c1 && c1 > 0 && c2 > 0) ?
      logMoment(t1, c1, t2, c2) : linMoment(t1, c1, t2, c2);
  }
  return aumc;
}

// ──────────────────────────────────────────────────────────────────────────
// Compensated (Neumaier) variants — see auc.ts for the summation rationale.
// ──────────────────────────────────────────────────────────────────────────

/** Linear trapezoidal AUMC with Neumaier-compensated summation. */
export function aumcLinearCompensated(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let sum = 0;
  let c = 0;
  for (let i = startIdx; i < endIdx; i++) {
    const v = linMoment(time[i], conc[i], time[i + 1], conc[i + 1]);
    const t = sum + v;
    if (Math.abs(sum) >= Math.abs(v))
      c += (sum - t) + v;
    else
      c += (v - t) + sum;
    sum = t;
  }
  return sum + c;
}

/** Log-linear trapezoidal AUMC with Neumaier-compensated summation. */
export function aumcLogLinearCompensated(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let sum = 0;
  let c = 0;
  for (let i = startIdx; i < endIdx; i++) {
    const t1 = time[i]; const c1 = conc[i];
    const t2 = time[i + 1]; const c2 = conc[i + 1];
    const v = (c1 > 0 && c2 > 0 && c1 !== c2) ?
      logMoment(t1, c1, t2, c2) : linMoment(t1, c1, t2, c2);
    const t = sum + v;
    if (Math.abs(sum) >= Math.abs(v))
      c += (sum - t) + v;
    else
      c += (v - t) + sum;
    sum = t;
  }
  return sum + c;
}

/** Linear-up / log-down trapezoidal AUMC with Neumaier-compensated summation. */
export function aumcLinearUpLogDownCompensated(
  time: Float64Array, conc: Float64Array,
  startIdx: number, endIdx: number,
): number {
  let sum = 0;
  let c = 0;
  for (let i = startIdx; i < endIdx; i++) {
    const t1 = time[i]; const c1 = conc[i];
    const t2 = time[i + 1]; const c2 = conc[i + 1];
    const v = (c2 < c1 && c1 > 0 && c2 > 0) ?
      logMoment(t1, c1, t2, c2) : linMoment(t1, c1, t2, c2);
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
 * Extrapolate the first moment from the last observed time to infinity,
 * assuming exponential decay at rate `lambdaZ` from `(tLast, cLast)` onwards.
 *
 * **Two-term tail** (the canonical AUMC guard): the moment of the tail is
 * `(tLast·cLast)/lambdaZ + cLast/lambdaZ²`, **not** the one-term `cLast/λz²`.
 * The one-term mistake silently under-reports AUMC/MRT/Vss.
 *
 * Mirrors {@link aucExtrapolateToInfinity}'s degenerate-case semantics:
 * returns `0` when `cLast == 0`, `+Infinity` when `lambdaZ == 0` with
 * `cLast > 0`, and `NaN` for the `0/0` case. Callers guard upstream —
 * `lambdaZ <= 0` means the terminal slope was not estimable.
 *
 * @param tLast - Time of the last observed point.
 * @param cLast - Concentration at the last observed point.
 * @param lambdaZ - Terminal-phase rate constant (1/time-unit).
 * @returns The extrapolated tail moment.
 */
export function aumcExtrapolateToInfinity(
  tLast: number, cLast: number, lambdaZ: number,
): number {
  return (tLast * cLast) / lambdaZ + cLast / (lambdaZ * lambdaZ);
}
