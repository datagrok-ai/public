import type {LambdaZResult, LambdaZStrategy} from './types';

/**
 * Auto best-fit terminal slope (lambda_z).
 *
 * The WINDOW-SELECTION tie-break matches PKNCA / WinNonlin's documented best-fit rule:
 * iterate over subsets of the last `k` eligible post-Cmax points, fit a log-linear OLS
 * regression on each, and pick the subset with the MOST points whose adjusted R² is
 * within `adjRSquaredFactor` of the maximum adjusted R² across candidates (see the loop
 * below and {@link LambdaZStrategy.adjRSquaredFactor}). A point is **eligible** when it
 * is post-Cmax (or Cmax itself, if `excludeCmax` is `false`), measurable
 * (`blqMask[i] === 0`), positive, and finite.
 *
 * A subset is **valid** when:
 * - it has at least `options.minPoints` points;
 * - the fitted slope is negative (so `lambda_z = -slope > 0`);
 * - `adjRSquared >= options.minRSquared`.
 *
 * Returns the best valid candidate, or `null` if none qualifies.
 *
 * **Scope of the PKNCA claim (deliberate divergence, NOT full parity).** Only the
 * tie-break above matches PKNCA. Two guards are sci-comp-specific and STRICTER than
 * PKNCA, and — unlike PKNCA — they drop candidates BEFORE the maximum adjusted R² is
 * computed (peer review R1, VAL-01-LZ-R019):
 * - `minRSquared` has **no** PKNCA equivalent — `pk.calc.half.life` carries no adj-R²
 *   floor parameter; it is an intentional sci-comp guardrail.
 * - PKNCA computes its adj-R² maximum over the sign-**unfiltered** candidate set and
 *   rejects `lambda.z <= 0` only at the final selection mask, so a spurious high-R²,
 *   wrong-signed short window can raise PKNCA's ceiling enough that no window survives →
 *   PKNCA declines to auto-select (NA / manual review). sci-comp drops such windows
 *   before the max and so always auto-selects the best VALID window instead.
 * These stricter guards are by design; aligning the max-over-unfiltered ordering (and a
 * PKNCA-style degenerate-2-point-fit warning) is tracked as a follow-up, not fixed here.
 *
 * Reference: <https://billdenney.github.io/pknca/articles/Selection-of-Calculation-Intervals.html>
 *
 * @param time - Time vector, sorted ascending.
 * @param conc - Concentration vector, same length as `time`.
 * @param blqMask - 1 = BLQ, 0 = measurable. Same length as `time`.
 * @param cmaxIdx - Index of the observed Cmax in `time` / `conc`.
 * @param options - Strategy parameters (`minPoints`, `minRSquared`, `excludeCmax`).
 * @returns The fit result, or `null` if no eligible subset qualifies.
 */
export function lambdaZBestFit(
  time: Float64Array, conc: Float64Array, blqMask: Uint8Array,
  cmaxIdx: number, options: LambdaZStrategy,
): LambdaZResult | null {
  const startIdx = options.excludeCmax ? cmaxIdx + 1 : cmaxIdx;
  const eligible: number[] = [];
  for (let i = startIdx; i < time.length; i++) {
    if (blqMask[i] !== 0) continue;
    const c = conc[i];
    if (!Number.isFinite(c) || c <= 0) continue;
    eligible.push(i);
  }
  if (eligible.length < options.minPoints) return null;

  const factor = options.adjRSquaredFactor ?? 0;
  // PKNCA / WinNonlin "best fit" window selection: fit every eligible terminal
  // window, then pick the one with the MOST points whose adjusted R² is within
  // `factor` of the MAXIMUM adjusted R² across all candidates. This is a FLAT
  // tolerance measured off the single global maximum — NOT an additive per-point
  // bonus.
  //
  // The previous implementation scored each window as `adjRSquared + factor·n`
  // and took the max. That accumulates the bonus with window length, so a longer
  // window only has to sit within `factor·Δn` of the best adj-R² to win — it
  // diverges from PKNCA on any profile with a candidate between (max − factor)
  // and (max − factor·Δn). Verified against PKNCA 0.12.1 on rat-IV R019
  // (VAL-01-LZ-R019): adj-R² peaks at n=4 (0.998164); every longer window is
  // >1e-4 below it, so PKNCA picks n=4 (λz=0.24047). The additive rule picked
  // n=8 (adjR²=0.997855, score 0.998655 > n=4's 0.998564) → λz=0.23494, off 2.3%.
  const candidates: {fit: LambdaZResult; n: number}[] = [];
  let maxAdjRSquared = -Infinity;
  for (let k = options.minPoints; k <= eligible.length; k++) {
    const subset = eligible.slice(eligible.length - k);
    const fit = fitLogLinear(time, conc, subset);
    if (fit === null) continue;
    if (fit.lambdaZ <= 0) continue;
    if (fit.adjRSquared < options.minRSquared) continue;
    candidates.push({fit, n: k});
    if (fit.adjRSquared > maxAdjRSquared) maxAdjRSquared = fit.adjRSquared;
  }
  // Among windows within `factor` of the global-max adj-R², keep the most points.
  // `factor = 0` collapses to "pick the window with the maximum adj-R², ties →
  // more points" — PKNCA's behaviour with tie-breaking disabled.
  let best: LambdaZResult | null = null;
  let bestN = -1;
  for (const {fit, n} of candidates) {
    if (fit.adjRSquared >= maxAdjRSquared - factor && n > bestN) {
      best = fit;
      bestN = n;
    }
  }
  return best;
}

/**
 * Manual lambda_z fit on a caller-supplied set of point indices.
 *
 * No best-fit search, no validity checks beyond `n >= 2` and
 * positivity / finiteness — the caller is presumed to know what they
 * picked. The returned slope can have any sign (caller decides whether to
 * reject `lambdaZ <= 0`).
 *
 * @param time - Time vector, sorted ascending.
 * @param conc - Concentration vector, same length as `time`.
 * @param pointIndices - Indices into `time` / `conc` of the chosen points.
 *                       Out-of-range or non-positive / non-finite entries
 *                       are silently dropped.
 * @returns The fit result, or `null` if fewer than 2 valid points remain.
 */
export function lambdaZManual(
  time: Float64Array, conc: Float64Array, pointIndices: Int32Array,
): LambdaZResult | null {
  const valid: number[] = [];
  for (let k = 0; k < pointIndices.length; k++) {
    const i = pointIndices[k];
    if (i < 0 || i >= time.length) continue;
    const c = conc[i];
    if (!Number.isFinite(c) || c <= 0) continue;
    valid.push(i);
  }
  if (valid.length < 2) return null;
  return fitLogLinear(time, conc, valid);
}

/**
 * Numerically stable closed-form OLS of `ln(conc) = intercept + slope·time`
 * on the supplied set of original-array indices.
 *
 * Uses centered sums to avoid catastrophic cancellation when `time` values
 * are large. Returns `null` when the design is degenerate (`n < 2` or all
 * `time` values equal in the subset).
 */
function fitLogLinear(
  time: Float64Array, conc: Float64Array, indices: number[],
): LambdaZResult | null {
  const n = indices.length;
  if (n < 2) return null;

  let sumX = 0;
  let sumY = 0;
  for (const i of indices) {
    sumX += time[i];
    sumY += Math.log(conc[i]);
  }
  const meanX = sumX / n;
  const meanY = sumY / n;

  let sxx = 0;
  let sxy = 0;
  let syy = 0;
  for (const i of indices) {
    const dx = time[i] - meanX;
    const dy = Math.log(conc[i]) - meanY;
    sxx += dx * dx;
    sxy += dx * dy;
    syy += dy * dy;
  }
  if (sxx === 0) return null;

  const slope = sxy / sxx;
  const intercept = meanY - slope * meanX;

  const ssRes = syy - slope * sxy; // = Σ(y - ŷ)²
  const rSquared = (syy === 0) ? 1 : 1 - ssRes / syy;
  const adjRSquared = (n > 2) ?
    1 - (1 - rSquared) * (n - 1) / (n - 2) :
    rSquared;

  return {
    lambdaZ: -slope,
    intercept,
    rSquared,
    adjRSquared,
    pointsUsed: Int32Array.from(indices),
    tStart: time[indices[0]],
    tEnd: time[indices[n - 1]],
  };
}
