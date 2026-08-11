import type {LambdaZResult, LambdaZStrategy} from './types';
import {halfLifeFromLambdaZ} from './derived';

/**
 * Auto best-fit terminal slope (lambda_z).
 *
 * Iterates over subsets of the last `k` eligible post-Cmax points, fits a log-linear
 * OLS regression on each, and picks the subset with the MOST points whose adjusted R²
 * is within `adjRSquaredFactor` of the maximum adjusted R² across candidates — the
 * PKNCA / WinNonlin best-fit tie-break. A point is **eligible** when it is post-Cmax
 * (or Cmax itself, if `excludeCmax` is `false`), measurable (`blqMask[i] === 0`),
 * positive, and finite.
 *
 * A subset is **valid** when it has at least `options.minPoints` points, the fitted
 * slope is negative (so `lambda_z = -slope > 0`), and `adjRSquared >= minRSquared`.
 * Returns the best valid candidate, or `null` if none qualifies.
 *
 * **Only the tie-break matches PKNCA.** Two guards are stricter, and both drop
 * candidates BEFORE the maximum adjusted R² is taken:
 * - `minRSquared` has no PKNCA equivalent — `pk.calc.half.life` carries no adj-R²
 *   floor; it is a sci-comp guardrail.
 * - PKNCA takes its maximum over the sign-unfiltered candidate set and rejects
 *   `lambda.z <= 0` only at the final selection mask, so a spurious high-R²,
 *   wrong-signed short window can raise its ceiling until nothing survives and PKNCA
 *   declines to auto-select. Dropping those windows first means sci-comp always
 *   auto-selects the best valid window instead.
 *
 * Reference: <https://billdenney.github.io/pknca/articles/Selection-of-Calculation-Intervals.html>
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
  // `factor` is a FLAT tolerance off the single global maximum adj-R², not an
  // additive per-point bonus: an additive `adjRSquared + factor·n` score
  // accumulates with window length and over-selects long windows.
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
  // `factor = 0` collapses to strict max adj-R², ties → more points.
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
 * positivity / finiteness — the caller is presumed to know what they picked.
 * Out-of-range or non-positive / non-finite `pointIndices` are silently
 * dropped; `null` is returned if fewer than 2 valid points remain. The
 * returned slope can have any sign (caller decides whether to reject
 * `lambdaZ <= 0`).
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

  const lambdaZ = -slope;
  const tStart = time[indices[0]];
  const tEnd = time[indices[n - 1]];

  // Diagnostic only — see LambdaZResult.spanRatio. The NaN branch is reachable
  // only via `lambdaZManual`, which permits a non-decaying fit; such a fit has
  // no half-life for a window to span.
  const spanRatio = lambdaZ > 0 ? (tEnd - tStart) / halfLifeFromLambdaZ(lambdaZ) : NaN;

  return {
    lambdaZ,
    intercept,
    rSquared,
    adjRSquared,
    pointsUsed: Int32Array.from(indices),
    tStart,
    tEnd,
    spanRatio,
  };
}
