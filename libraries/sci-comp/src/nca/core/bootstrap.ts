/**
 * Stratified-by-timepoint bootstrap for sparse NCA (UC-04 / FR-304).
 *
 * The closed-form Holder variance ({@link sparseAuc}) is the regulatory AUC
 * interval, but the **nonlinear** parameters (Cmax, t½, CL) have no closed-form
 * sparse variance. For those, this module resamples animals **within each
 * nominal timepoint** (the analogue of the destructive design), re-aggregates,
 * recomputes the statistic, and reports a **BCa** interval (the sparse-AUC
 * bootstrap is right-skewed and bounded at 0, so the bias-corrected accelerated
 * interval is the right one — AD-6).
 *
 * ## Honest-failure gate (thesis #4 / Bonate 1998)
 *
 * The stratified bootstrap **degenerates combinatorially** at small n: the number
 * of distinct stratified resamples is `h = Π_i C(2n_i − 1, n_i)`. A 5-timepoint
 * × 2-animal destructive design has `h = 3⁵ = 243 ≤ 360` (Bonate 1998 eq 8) — the
 * "bootstrap distribution" is then a handful of repeated points, not a
 * distribution. {@link summarizeBootstrap} **suppresses** the interval in that
 * regime (`suppressed: true`) and the caller defers to the closed-form CI. Two
 * gates apply:
 *  - **Unconditional hard floor** (AC-U4 / R1-plan NB-2): `h ≤ 360` is suppressed
 *    *regardless* of the calibrated gate. Calibration can only make the gate
 *    stricter, never admit this degenerate case.
 *  - **Calibrated min-n** (`minNPerTimepoint`, default 3): any timepoint with
 *    fewer animals suppresses. Calibrated against the G-SP-08 fixture at build
 *    (Rule 16) — not a hardcoded scientific constant.
 *
 * Reproducibility: at a fixed `masterSeed` the resampling is deterministic
 * ({@link mulberry32}). The closed-form CI is **never** gated (Holder holds to n=2).
 *
 * @see Bonate PL (1998) "Coverage and precision of confidence intervals for area
 *   under the curve using parametric and non-parametric methods in a sparse
 *   sampling design." Pharm Res 15(3):405-410.
 * @see Efron B, Tibshirani RJ (1993) — BCa intervals.
 */

import {mulberry32} from './prng';
import {normalInv, normalCdf} from '../../stats/distributions';
import type {SparseInput} from './sparse';

/** Unconditional hard floor on distinct stratified resamples (Bonate eq 8). */
export const BOOTSTRAP_DISTINCT_RESAMPLE_FLOOR = 360;
/** Default calibrated minimum animals per timepoint (calibrate at build, Rule 16). */
export const DEFAULT_MIN_N_PER_TIMEPOINT = 3;
/** Default bootstrap iterations. */
export const DEFAULT_BOOTSTRAP_ITERATIONS = 1000;

export interface BootstrapOptions {
  /** Resamples to draw. Default 1000. */
  readonly iterations?: number;
  /** PRNG master seed (reproducibility). Default 12345. */
  readonly masterSeed?: number;
  /** Two-sided CI confidence level. Default 0.95. */
  readonly ciLevel?: number;
  /** Calibrated minimum animals per nominal timepoint. Default 3. */
  readonly minNPerTimepoint?: number;
}

export interface BootstrapSummary {
  /** The statistic on the original (un-resampled) data. */
  readonly estimate: number;
  /** Median of the bootstrap replicates (`NaN` when suppressed). */
  readonly median: number;
  /** Confidence interval (BCa, or percentile fallback). `[NaN, NaN]` when suppressed. */
  readonly ci: readonly [number, number];
  /** Which interval method was used. */
  readonly ciMethod: 'BCa' | 'percentile' | 'none';
  /** Replicates actually computed (0 when suppressed). */
  readonly iterations: number;
  /** True when the min-n gate (or hard floor) suppressed the bootstrap. */
  readonly suppressed: boolean;
  /** Human-readable reason when suppressed, else null. */
  readonly suppressReason: string | null;
}

/** One nominal timepoint's row indices into the flat SparseInput. */
interface TimeStratum {
  readonly time: number;
  readonly rows: number[];
}

function stratify(input: SparseInput): TimeStratum[] {
  const map = new Map<number, number[]>();
  for (let i = 0; i < input.nominalTime.length; i++) {
    const t = input.nominalTime[i];
    if (!Number.isFinite(t)) continue;
    let rows = map.get(t);
    if (rows === undefined) {rows = []; map.set(t, rows);}
    rows.push(i);
  }
  return [...map.entries()].sort((a, b) => a[0] - b[0]).map(([time, rows]) => ({time, rows}));
}

/** Number of distinct stratified resamples `h = Π_i C(2n_i − 1, n_i)` (Bonate eq 8). */
function distinctResamples(strata: TimeStratum[]): number {
  let h = 1;
  for (const s of strata) {
    const n = s.rows.length;
    h *= binomial(2 * n - 1, n);
    if (!Number.isFinite(h) || h > 1e15) return Infinity; // large enough — not degenerate
  }
  return h;
}

function binomial(n: number, k: number): number {
  if (k < 0 || k > n) return 0;
  k = Math.min(k, n - k);
  let c = 1;
  for (let i = 0; i < k; i++) c = c * (n - i) / (i + 1);
  return Math.round(c);
}

/** Build a resampled SparseInput by drawing rows with replacement within strata. */
function resample(input: SparseInput, strata: TimeStratum[], rand: () => number): SparseInput {
  const total = input.nominalTime.length;
  const pick: number[] = [];
  for (const s of strata) {
    const n = s.rows.length;
    for (let k = 0; k < n; k++) pick.push(s.rows[Math.floor(rand() * n)]);
  }
  const m = pick.length;
  const nominalTime = new Float64Array(m);
  const conc = new Float64Array(m);
  const blqMask = new Uint8Array(m);
  const animalId = input.animalId !== null ? new Int32Array(m) : null;
  const scalarLloq = typeof input.lloq === 'number';
  const lloq = scalarLloq ? input.lloq : new Float64Array(m);
  for (let i = 0; i < m; i++) {
    const src = pick[i];
    nominalTime[i] = input.nominalTime[src];
    conc[i] = input.conc[src];
    blqMask[i] = input.blqMask[src];
    if (animalId !== null) animalId[i] = (input.animalId as Int32Array)[src];
    if (!scalarLloq) (lloq as Float64Array)[i] = (input.lloq as Float64Array)[src];
  }
  void total;
  return {nominalTime, conc, blqMask, lloq, animalId};
}

function percentile(sorted: number[], p: number): number {
  if (sorted.length === 0) return NaN;
  if (p <= 0) return sorted[0];
  if (p >= 1) return sorted[sorted.length - 1];
  const idx = p * (sorted.length - 1);
  const lo = Math.floor(idx);
  const hi = Math.ceil(idx);
  if (lo === hi) return sorted[lo];
  return sorted[lo] + (idx - lo) * (sorted[hi] - sorted[lo]);
}

/**
 * Bootstrap a statistic on a sparse design with stratified-by-timepoint
 * resampling and a BCa interval, min-n gated.
 *
 * @param input - The sparse observations.
 * @param statistic - Computes the scalar of interest from a (resampled) input.
 *   Return `NaN` for a degenerate resample; NaN replicates are dropped.
 * @param options - Iterations, seed, CI level, calibrated min-n.
 */
export function summarizeBootstrap(
  input: SparseInput,
  statistic: (resampled: SparseInput) => number,
  options: BootstrapOptions = {},
): BootstrapSummary {
  const iterations = options.iterations ?? DEFAULT_BOOTSTRAP_ITERATIONS;
  const masterSeed = options.masterSeed ?? 12345;
  const ciLevel = options.ciLevel ?? 0.95;
  const minN = options.minNPerTimepoint ?? DEFAULT_MIN_N_PER_TIMEPOINT;

  const strata = stratify(input);
  const estimate = statistic(input);

  // --- Gates.
  const minObserved = strata.reduce((acc, s) => Math.min(acc, s.rows.length), Infinity);
  const h = distinctResamples(strata);
  const suppressed = ((): string | null => {
    if (h <= BOOTSTRAP_DISTINCT_RESAMPLE_FLOOR) {
      return `only ${h} distinct stratified resamples (≤ ${BOOTSTRAP_DISTINCT_RESAMPLE_FLOOR}); ` +
        'the bootstrap is combinatorially degenerate (Bonate 1998 eq 8)';
    }
    if (minObserved < minN) {
      return `a nominal timepoint has ${minObserved} animals (< ${minN}); ` +
        'the stratified bootstrap is unreliable at this n';
    }
    return null;
  })();
  if (suppressed !== null) {
    return {
      estimate, median: NaN, ci: [NaN, NaN], ciMethod: 'none',
      iterations: 0, suppressed: true, suppressReason: suppressed,
    };
  }

  // --- Resample.
  const rand = mulberry32(masterSeed);
  const reps: number[] = [];
  for (let b = 0; b < iterations; b++) {
    const v = statistic(resample(input, strata, rand));
    if (Number.isFinite(v)) reps.push(v);
  }
  if (reps.length < 2) {
    return {
      estimate, median: NaN, ci: [NaN, NaN], ciMethod: 'none',
      iterations: reps.length, suppressed: true,
      suppressReason: 'too few finite bootstrap replicates to form an interval',
    };
  }
  reps.sort((a, b) => a - b);
  const median = percentile(reps, 0.5);

  // --- BCa.
  const alpha = (1 - ciLevel) / 2;
  const bca = bcaInterval(reps, estimate, input, statistic, alpha);
  if (bca !== null) {
    return {
      estimate, median, ci: bca, ciMethod: 'BCa',
      iterations: reps.length, suppressed: false, suppressReason: null,
    };
  }
  // --- Percentile fallback (BCa undefined — e.g. zero jackknife variance).
  return {
    estimate, median,
    ci: [percentile(reps, alpha), percentile(reps, 1 - alpha)],
    ciMethod: 'percentile', iterations: reps.length, suppressed: false, suppressReason: null,
  };
}

/**
 * BCa interval (Efron). Bias-correction `z0` from the replicate rank of the
 * estimate; acceleration `a` from a leave-one-observation-out jackknife. Returns
 * null when `z0`/`a` are not finite (caller falls back to percentile).
 */
function bcaInterval(
  reps: number[], estimate: number, input: SparseInput,
  statistic: (resampled: SparseInput) => number, alpha: number,
): [number, number] | null {
  const B = reps.length;
  // Bias correction.
  let less = 0;
  for (const v of reps) if (v < estimate) less++;
  const prop = less / B;
  if (prop <= 0 || prop >= 1) return null;
  const z0 = normalInv(prop);
  if (!Number.isFinite(z0)) return null;

  // Acceleration via jackknife (leave one observation out).
  const n = input.nominalTime.length;
  const jack: number[] = [];
  for (let i = 0; i < n; i++) {
    const v = statistic(omitRow(input, i));
    if (Number.isFinite(v)) jack.push(v);
  }
  if (jack.length < 2) return null;
  let jbar = 0;
  for (const v of jack) jbar += v;
  jbar /= jack.length;
  let num = 0; let den = 0;
  for (const v of jack) {const d = jbar - v; num += d * d * d; den += d * d;}
  if (den === 0) return null;
  const a = num / (6 * Math.pow(den, 1.5));
  if (!Number.isFinite(a)) return null;

  const zLo = normalInv(alpha);
  const zHi = normalInv(1 - alpha);
  const adj = (z: number): number => {
    const denom = 1 - a * (z0 + z);
    if (denom === 0) return NaN;
    return normalCdf(z0 + (z0 + z) / denom);
  };
  const a1 = adj(zLo);
  const a2 = adj(zHi);
  if (!Number.isFinite(a1) || !Number.isFinite(a2)) return null;
  return [percentile(reps, a1), percentile(reps, a2)];
}

/** SparseInput with row `omit` removed (jackknife). */
function omitRow(input: SparseInput, omit: number): SparseInput {
  const n = input.nominalTime.length;
  const m = n - 1;
  const nominalTime = new Float64Array(m);
  const conc = new Float64Array(m);
  const blqMask = new Uint8Array(m);
  const animalId = input.animalId !== null ? new Int32Array(m) : null;
  const scalarLloq = typeof input.lloq === 'number';
  const lloq = scalarLloq ? input.lloq : new Float64Array(m);
  let j = 0;
  for (let i = 0; i < n; i++) {
    if (i === omit) continue;
    nominalTime[j] = input.nominalTime[i];
    conc[j] = input.conc[i];
    blqMask[j] = input.blqMask[i];
    if (animalId !== null) animalId[j] = (input.animalId as Int32Array)[i];
    if (!scalarLloq) (lloq as Float64Array)[j] = (input.lloq as Float64Array)[i];
    j++;
  }
  return {nominalTime, conc, blqMask, lloq, animalId};
}
