import {lambdaZBestFit, lambdaZManual} from '../lambda-z';
import type {LambdaZStrategy} from '../types';

const DEFAULT_STRATEGY: LambdaZStrategy = {
  mode: 'auto-best-fit',
  minPoints: 3,
  minRSquared: 0.85,
  excludeCmax: true,
};

/** Build a clean exponential profile c(t) = c0·exp(-k·t). */
function exponential(N: number, k: number, dt: number, c0 = 1) {
  const time = new Float64Array(N);
  const conc = new Float64Array(N);
  for (let i = 0; i < N; i++) {
    time[i] = i * dt;
    conc[i] = c0 * Math.exp(-k * time[i]);
  }
  return {time, conc};
}

describe('lambdaZBestFit', () => {
  it('recovers the rate constant exactly on a clean exponential', () => {
    const k = 0.25;
    const {time, conc} = exponential(10, k, 1);
    const blqMask = new Uint8Array(10);
    // cmaxIdx = 0 (decay starts at t=0); excludeCmax=true → use indices 1..9.
    const r = lambdaZBestFit(time, conc, blqMask, 0, DEFAULT_STRATEGY)!;
    expect(r).not.toBeNull();
    expect(Math.abs(r.lambdaZ - k)).toBeLessThan(1e-10);
    expect(r.rSquared).toBeGreaterThan(0.9999999);
  });

  it('on a noisy exponential lambda_z stays close to truth and adjR² stays high', () => {
    const k = 0.1;
    const {time, conc} = exponential(20, k, 1, 100);
    // Add 1% multiplicative noise (deterministic — reproducible via mulberry-style step).
    const noisy = new Float64Array(conc.length);
    let seed = 0;
    for (let i = 0; i < conc.length; i++) {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      const u = (seed / 0x7fffffff) - 0.5; // u ∈ [-0.5, 0.5]
      noisy[i] = conc[i] * (1 + 0.01 * u);
    }
    const blqMask = new Uint8Array(20);
    const r = lambdaZBestFit(time, noisy, blqMask, 0, DEFAULT_STRATEGY)!;
    expect(r).not.toBeNull();
    expect(Math.abs(r.lambdaZ - k) / k).toBeLessThan(0.01); // within 1%
    expect(r.adjRSquared).toBeGreaterThan(0.99);
  });

  it('returns null on a strictly ascending profile (slope >= 0)', () => {
    const time = new Float64Array([0, 1, 2, 3, 4]);
    const conc = new Float64Array([1, 2, 4, 8, 16]);
    const blqMask = new Uint8Array(5);
    // cmaxIdx = 4 (last). With excludeCmax=true, eligible starts at 5 → empty.
    // But that's a separate failure mode. To exercise the slope rejection,
    // pretend Cmax is at the start.
    const r = lambdaZBestFit(time, conc, blqMask, 0, DEFAULT_STRATEGY);
    expect(r).toBeNull();
  });

  it('excludeCmax=true keeps Cmax out of the eligible window', () => {
    // 4-point profile, cmaxIdx=1. minPoints=3.
    //   excludeCmax=true  → eligible = {2, 3} → 2 < 3 → null.
    //   excludeCmax=false → eligible = {1, 2, 3} → fit succeeds, includes Cmax.
    const time = new Float64Array([0, 1, 2, 3]);
    const conc = new Float64Array([1, 100, 50, 25]);
    const blqMask = new Uint8Array(4);
    const cmaxIdx = 1;

    expect(lambdaZBestFit(time, conc, blqMask, cmaxIdx,
      {...DEFAULT_STRATEGY, excludeCmax: true})).toBeNull();

    const withInclude = lambdaZBestFit(time, conc, blqMask, cmaxIdx,
      {...DEFAULT_STRATEGY, excludeCmax: false})!;
    expect(withInclude).not.toBeNull();
    expect(Array.from(withInclude.pointsUsed)).toContain(cmaxIdx);
    expect(Math.abs(withInclude.lambdaZ - Math.log(2))).toBeLessThan(1e-10);
  });

  it('honours minPoints — returns null when fewer eligible points than required', () => {
    // Only 2 post-Cmax measurable points, but minPoints = 3.
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([10, 5, 2.5]);
    const blqMask = new Uint8Array(3);
    const r = lambdaZBestFit(time, conc, blqMask, 0,
      {...DEFAULT_STRATEGY, minPoints: 3});
    expect(r).toBeNull();
  });

  it('honours minRSquared — returns null when no subset reaches the threshold', () => {
    // Random-walk-ish data: no strong log-linear structure.
    const time = new Float64Array([0, 1, 2, 3, 4, 5]);
    const conc = new Float64Array([10, 7, 9, 5, 8, 4]);
    const blqMask = new Uint8Array(6);
    const r = lambdaZBestFit(time, conc, blqMask, 0,
      {...DEFAULT_STRATEGY, minRSquared: 0.999});
    expect(r).toBeNull();
  });

  it('skips BLQ points and NaN concentrations', () => {
    // 6 points; index 1 is BLQ, index 3 is NaN. Both must be ignored.
    const time = new Float64Array([0, 1, 2, 3, 4, 5]);
    const conc = new Float64Array([1, 999, 0.5, NaN, 0.125, 0.0625]);
    const blqMask = new Uint8Array([0, 1, 0, 0, 0, 0]);
    const r = lambdaZBestFit(time, conc, blqMask, 0, DEFAULT_STRATEGY)!;
    expect(r).not.toBeNull();
    // Subset must not contain index 1 or 3
    const used = Array.from(r.pointsUsed);
    expect(used).not.toContain(1);
    expect(used).not.toContain(3);
  });

  it('picks the subset with the highest adjusted R² when several are valid', () => {
    // Construct a profile where the last-3 fit is cleaner than longer ones.
    // First 3 points have small drift; last 3 are perfectly log-linear.
    const time = new Float64Array([0, 1, 2, 3, 4, 5]);
    const k = 0.5;
    const c0 = 100;
    const conc = new Float64Array([
      c0, // index 0
      c0 * Math.exp(-k * 1) * 1.05, // index 1: 5% off
      c0 * Math.exp(-k * 2) * 0.95, // index 2: -5% off
      c0 * Math.exp(-k * 3),
      c0 * Math.exp(-k * 4),
      c0 * Math.exp(-k * 5),
    ]);
    const blqMask = new Uint8Array(6);
    const r = lambdaZBestFit(time, conc, blqMask, 0, DEFAULT_STRATEGY)!;
    expect(r).not.toBeNull();
    // Best subset should include the last 3 perfectly-aligned points.
    const used = Array.from(r.pointsUsed);
    expect(used).toContain(3);
    expect(used).toContain(4);
    expect(used).toContain(5);
    expect(Math.abs(r.lambdaZ - k)).toBeLessThan(0.05);
  });

  it('PKNCA flat-tolerance window selection: real rat-IV R019 profile → n=4 (VAL-01-LZ-R019)', () => {
    // Regression for the additive-vs-flat adjRSquaredFactor bug. Ground truth from
    // PKNCA 0.12.1 (auc.method="lin up/log down", min.hl.points=3, min.hl.r.squared=0.85,
    // allow.tmax.in.half.life=FALSE, min.span.ratio=2) on rat dataset 03 subject R019:
    //   lambda.z = 0.2404699, half.life = 2.88247, lambda.z.n.points = 4, window t=[6,24].
    // Per-window adj-R² peaks at n=4 (0.998164); every longer window is > 1e-4 below it,
    // so PKNCA's FLAT tolerance keeps n=4. The prior additive score `adjR²+1e-4·n` instead
    // peaked at n=8 (adjR²=0.997855, score 0.998655 > n=4's 0.998564) → λz=0.23494, 2.3% off
    // (exceeds NFR-05's 0.5%). This test fails on the additive rule, passes on the flat rule.
    const time = Float64Array.of(0.0, 0.083, 0.25, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0);
    const conc = Float64Array.of(5656.92, 4937.10, 5405.76, 4160.80, 3382.91, 2474.69,
      1591.91, 1177.43, 646.56, 289.48, 14.85);
    const blqMask = new Uint8Array(time.length);
    const strategy: LambdaZStrategy = {
      mode: 'auto-best-fit', minPoints: 3, minRSquared: 0.85, excludeCmax: true,
      adjRSquaredFactor: 1e-4, // PKNCA default — where additive and flat diverge
    };
    // Cmax is the t=0 IV-bolus point (index 0); excludeCmax drops it.
    const r = lambdaZBestFit(time, conc, blqMask, 0, strategy)!;
    expect(r).not.toBeNull();
    expect(Array.from(r.pointsUsed)).toEqual([7, 8, 9, 10]); // t = 6, 8, 12, 24
    expect(r.tStart).toBe(6);
    expect(r.tEnd).toBe(24);
    expect(Math.abs(r.lambdaZ - 0.2404699) / 0.2404699).toBeLessThan(5e-4); // within PKNCA 0.5%
    expect(Math.abs(Math.LN2 / r.lambdaZ - 2.88247)).toBeLessThan(0.01); // t½ within 0.01 h
  });

  it('factor=0 tie-break: on exact adjR² ties, the window with MORE points wins (PKNCA convention)', () => {
    // A perfect exponential makes every candidate window fit with adjR² = 1 (an exact,
    // bit-level tie). With adjRSquaredFactor unset (factor=0), the flat-tolerance rule
    // keeps the MOST points among the tied windows — PKNCA's documented "ties → more
    // points" direction (half.life.R: n.points == max(n.points[mask_best])). The prior
    // additive score kept the FIRST/smallest window on exact ties; this locks the corrected
    // direction so it can't silently regress.
    const k = 0.25;
    const {time, conc} = exponential(10, k, 1); // indices 0..9, clean decay from t=0
    const blqMask = new Uint8Array(10);
    // cmaxIdx=0, excludeCmax=true → eligible = indices 1..9 (9 points); every window adjR²=1.
    const r = lambdaZBestFit(time, conc, blqMask, 0, DEFAULT_STRATEGY)!;
    expect(r).not.toBeNull();
    expect(r.pointsUsed.length).toBe(9); // most points wins the exact tie, not the shortest
    expect(Array.from(r.pointsUsed)).toEqual([1, 2, 3, 4, 5, 6, 7, 8, 9]);
    expect(Math.abs(r.lambdaZ - k)).toBeLessThan(1e-10);
  });

  it('reports tStart, tEnd, pointsUsed, intercept consistently', () => {
    const k = 0.2;
    const c0 = 50;
    const {time, conc} = exponential(8, k, 1, c0);
    const blqMask = new Uint8Array(8);
    const r = lambdaZBestFit(time, conc, blqMask, 0, DEFAULT_STRATEGY)!;
    expect(r.pointsUsed[0]).toBe(used0(r.pointsUsed));
    expect(r.tStart).toBe(time[r.pointsUsed[0]]);
    expect(r.tEnd).toBe(time[r.pointsUsed[r.pointsUsed.length - 1]]);
    // Intercept of OLS on log space ≈ ln(c0)
    expect(Math.abs(r.intercept - Math.log(c0))).toBeLessThan(1e-10);
  });

  function used0(arr: Int32Array): number {
    return arr[0];
  }
});

describe('lambdaZManual', () => {
  it('fits the exact provided indices', () => {
    const k = 0.3;
    const c0 = 25;
    const {time, conc} = exponential(10, k, 1, c0);
    const r = lambdaZManual(time, conc, Int32Array.of(5, 6, 7, 8, 9))!;
    expect(r).not.toBeNull();
    expect(Math.abs(r.lambdaZ - k)).toBeLessThan(1e-10);
    expect(Array.from(r.pointsUsed)).toEqual([5, 6, 7, 8, 9]);
  });

  it('returns null when fewer than 2 valid indices remain', () => {
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([1, 0.5, 0.25]);
    expect(lambdaZManual(time, conc, Int32Array.of(0))).toBeNull();
    expect(lambdaZManual(time, conc, new Int32Array(0))).toBeNull();
  });

  it('drops out-of-range and non-positive indices silently', () => {
    const time = new Float64Array([0, 1, 2, 3]);
    const conc = new Float64Array([1, 0.5, 0, 0.125]); // index 2 is zero
    // Provide [99, 0, 2, 1, 3] — 99 out of range, 2 is non-positive.
    const r = lambdaZManual(time, conc, Int32Array.of(99, 0, 2, 1, 3))!;
    expect(r).not.toBeNull();
    expect(Array.from(r.pointsUsed)).toEqual([0, 1, 3]);
  });

  it('does not auto-reject ascending fits — caller decides on lambdaZ <= 0', () => {
    const time = new Float64Array([0, 1, 2]);
    const conc = new Float64Array([1, 2, 4]);
    const r = lambdaZManual(time, conc, Int32Array.of(0, 1, 2))!;
    expect(r).not.toBeNull();
    expect(r.lambdaZ).toBeLessThan(0); // slope is positive → lambda_z < 0
  });
});

describe('numerical stability', () => {
  it('OLS is stable when time values are large', () => {
    // Time values around 1e6, exponential decay over a small range.
    const k = 0.1;
    const N = 10;
    const time = new Float64Array(N);
    const conc = new Float64Array(N);
    for (let i = 0; i < N; i++) {
      time[i] = 1e6 + i;
      conc[i] = Math.exp(-k * (time[i] - 1e6));
    }
    const blqMask = new Uint8Array(N);
    const r = lambdaZBestFit(time, conc, blqMask, 0, DEFAULT_STRATEGY)!;
    expect(r).not.toBeNull();
    expect(Math.abs(r.lambdaZ - k)).toBeLessThan(1e-9);
  });

  it('all-equal time → returns null (degenerate design)', () => {
    const time = new Float64Array([5, 5, 5]);
    const conc = new Float64Array([1, 0.5, 0.25]);
    expect(lambdaZManual(time, conc, Int32Array.of(0, 1, 2))).toBeNull();
  });

  it('all-equal conc → lambdaZ=0 (rejected by best-fit, returned by manual)', () => {
    // syy=0 branch in OLS: all log(conc) equal → R² = 1 by convention.
    const time = new Float64Array([0, 1, 2, 3]);
    const conc = new Float64Array([5, 5, 5, 5]);
    const blqMask = new Uint8Array(4);
    expect(lambdaZBestFit(time, conc, blqMask, 0, DEFAULT_STRATEGY)).toBeNull();

    const r = lambdaZManual(time, conc, Int32Array.of(0, 1, 2, 3))!;
    expect(r).not.toBeNull();
    expect(Math.abs(r.lambdaZ)).toBe(0); // ±0 both acceptable
    expect(r.rSquared).toBe(1);
  });

  it('n=2 manual fit reports adjR² = rSquared (no DoF adjustment)', () => {
    const time = new Float64Array([0, 1]);
    const conc = new Float64Array([1, 0.5]);
    const r = lambdaZManual(time, conc, Int32Array.of(0, 1))!;
    expect(r.rSquared).toBe(r.adjRSquared);
    expect(r.rSquared).toBe(1); // exactly two points → perfect fit
  });
});
