import type {LambdaZResult, LambdaZStrategy} from './types';

/**
 * Auto best-fit terminal slope (lambda_z).
 *
 * PKNCA-compatible algorithm: iterate over subsets of the last `k` eligible
 * post-Cmax points, fit a log-linear OLS regression on each subset, and
 * pick the candidate that maximises adjusted R². A point is **eligible**
 * when it is post-Cmax (or Cmax itself, if `excludeCmax` is `false`),
 * measurable (`blqMask[i] === 0`), positive, and finite.
 *
 * A subset is **valid** when:
 * - it has at least `options.minPoints` points;
 * - the fitted slope is negative (so `lambda_z = -slope > 0`);
 * - `adjRSquared >= options.minRSquared`.
 *
 * Returns the best valid candidate, or `null` if none qualifies.
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
  let best: LambdaZResult | null = null;
  let bestScore = -Infinity;
  for (let k = options.minPoints; k <= eligible.length; k++) {
    const subset = eligible.slice(eligible.length - k);
    const fit = fitLogLinear(time, conc, subset);
    if (fit === null) continue;
    if (fit.lambdaZ <= 0) continue;
    if (fit.adjRSquared < options.minRSquared) continue;
    // Score favours larger subsets when adj-R² values are close
    // (PKNCA convention; factor = 1e-4 in their defaults).
    const score = fit.adjRSquared + factor * subset.length;
    if (score > bestScore) {
      best = fit;
      bestScore = score;
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
