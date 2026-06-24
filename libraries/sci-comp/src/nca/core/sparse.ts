/**
 * Sparse / destructive-sampling NCA — design-aware closed-form AUC with an
 * honest standard error and degrees of freedom (UC-04 / FR-301..306).
 *
 * When each animal contributes only one (destructive) or a few (batch) samples,
 * per-subject NCA is undefined. The standard approach builds a **composite**
 * mean concentration-time profile per nominal timepoint and integrates that.
 * The statistically correct AUC variance for such a design is **Bailer's method**
 * (destructive) generalised to **Holder's covariance estimator** (batch), with a
 * **Nedelman-Jia correlated Satterthwaite df** for the t-based CI.
 *
 * ## What this module computes (validated against PKNCA 0.12.1, AC-U1..U3)
 *
 * - {@link buildCompositeProfile}: arithmetic mean / SD / %CV / n / %BLQ per
 *   nominal timepoint, with the PKNCA `arithmetic mean, <=50% BLQ` rule (a
 *   timepoint whose BLQ fraction exceeds 50 % contributes a zero mean to the
 *   AUClast point estimate). BLQ is imputed **before** averaging (AD-5).
 * - {@link sparseAuc}: the composite **AUClast** (linear-trapezoidal, anchored
 *   at t=0, C=0), its Holder variance, the Nedelman-Jia df, and the Student-t CI.
 *
 * ## The math (reproduced verbatim from the primary sources)
 *
 * Lever-arm linear-trapezoidal weights over the design timepoints (with the
 * t=0 anchor): `w = c(0, diff(t)/2) + c(diff(t)/2, 0)` — PKNCA
 * `sparse_auc_weight_linear`. The AUC estimator is `AÛC = Σ_i w_i · C̄_i`.
 *
 * Holder 2001 eq (A1) variance (`ψ`):
 * ```
 * ψ = Σ_i w_i² σ̂_i²/r_i  +  2 Σ_{i<j} w_i w_j r_ij σ̂_ij /(r_i r_j)
 * ```
 * with the **unbiased** covariance estimator, Holder 2001 eq (A3) (the biased
 * eq A2 overestimates ψ by 2–6 % at ρ=0.6–0.9 — Holder Table 2):
 * ```
 * σ̂_ij = Σ_k (x_ik − x̄_i)(x_jk − x̄_j) / [ (r_ij − 1) + (1 − r_ij/r_i)(1 − r_ij/r_j) ]
 * ```
 * where r_i = animals at t_i, r_ij = animals sampled at **both** t_i and t_j (over
 * the r_ij animals common to both). When all r_ij = 0 this **reduces to Bailer**.
 *
 * Nedelman-Jia 1998 correlated Satterthwaite df (J Biopharm Stat 8(2):320-321),
 * matrix form (reduces to the scalar independence formula when r_ij = 0):
 * ```
 * df = ψ² / [ Σ_i (w_i² Σ_ii)²/(r_i−1) + 2 Σ_{i<j} (w_i w_j Σ_ij)²/(r_ij−1) ]
 * ```
 * where Σ is the covariance matrix of the timepoint means (Σ_ii = σ̂_i²/r_i,
 * Σ_ij = r_ij σ̂_ij/(r_i r_j)).
 *
 * ## Scientific guards (the honesty core — thesis #4)
 *
 * - **Linear-trapezoidal only** (AD-3): the Holder variance is valid only where
 *   AUC is linear in the timepoint means. A steep, wide-gap terminal phase makes
 *   the linear rule over-estimate a convex decline; {@link sparseAuc} raises
 *   `SPARSE_TERMINAL_OVEREST` (Jia-Nedelman 1996).
 * - **Topology from the data, not the label** (AD-2): destructive / batch / serial
 *   is derived from the r_ij overlap matrix and cross-checked against the declared
 *   label (`SPARSE_TOPOLOGY_MISMATCH`).
 * - **Variance-borrowing for n=1** (AD-7 / FR-306): a singleton timepoint has no
 *   sample variance; rather than silently zeroing it (→ false-narrow CI), a
 *   Nedelman variance-vs-mean power law fit across the n≥2 timepoints models it
 *   (`SPARSE_VARIANCE_MODELED`).
 *
 * @see Holder DJ (2001) "Comments on Nedelman and Jia's Extension of Satterthwaite's
 *   Approximation Applied to Pharmacokinetics." J Biopharm Stat 11(1-2):75-79.
 * @see Nedelman JR, Jia X (1998) J Biopharm Stat 8(2):317-328.
 * @see Bailer AJ (1988) J Pharmacokinet Biopharm 16(3):303-309.
 * @see Jia X, Nedelman J (1996) — log-trap nonlinearity / terminal-gap AUC bias.
 */

import type {BlqRule} from './types';
import {studentTInv} from '../../stats/distributions';

// ──────────────────────────────────────────────────────────────────────────
// Public types
// ──────────────────────────────────────────────────────────────────────────

/** Data-derived sampling topology (from the r_ij overlap matrix). */
export type SamplingTopology = 'destructive' | 'batch' | 'serial';

/**
 * Flat, columnar sparse input — one entry per (animal, nominal-time) observation.
 * The adapter copies these out of a DataFrame; {@link sparseAuc} and
 * {@link buildCompositeProfile} group internally by nominal time. The bootstrap
 * resamples this same shape.
 */
export interface SparseInput {
  /** Nominal time of each observation (the binning key). */
  readonly nominalTime: Float64Array;
  /** Observed concentration of each observation (raw; BLQ rows flagged by blqMask). */
  readonly conc: Float64Array;
  /** 1 = below LOQ at this observation, 0 = measurable. Same length as conc. */
  readonly blqMask: Uint8Array;
  /** Per-row LLOQ, or a scalar applied to all rows (used by `set-half-lloq`). */
  readonly lloq: Float64Array | number;
  /**
   * Integer animal id per observation; `null` → the design is assumed
   * **destructive** (each observation is its own animal — r_ij = 0). When
   * present, r_ij is derived from the animal × nominal-time table.
   */
  readonly animalId: Int32Array | null;
}

export interface SparseAucOptions {
  /** Two-sided CI confidence level. Default 0.95. */
  readonly ciLevel?: number;
  /** Per-record BLQ rule applied before averaging (AD-5). Default `'set-zero'`. */
  readonly blqRule?: BlqRule;
  /**
   * Estimated terminal half-life, for the `SPARSE_TERMINAL_OVEREST` heuristic
   * (AD-3). When absent, the heuristic falls back to the concentration-only
   * criterion (terminal mean < 20 % Cmax) so the flag never silently disappears.
   */
  readonly estimatedHalfLife?: number;
  /** Enable Nedelman variance-borrowing for n=1 timepoints (AD-7). Default true. */
  readonly varianceBorrow?: boolean;
  /** Declared topology label, cross-checked against the data (AD-2). */
  readonly declaredTopology?: SamplingTopology;
}

export type SparseWarningCode =
  | 'SPARSE_TOPOLOGY_MISMATCH'
  | 'SPARSE_TERMINAL_OVEREST'
  | 'SPARSE_VARIANCE_MODELED'
  | 'SPARSE_ANIMAL_ID_ABSENT_ASSUMED_DESTRUCTIVE'
  | 'SPARSE_DF_UNAVAILABLE';

export interface SparseWarning {
  readonly code: SparseWarningCode;
  readonly message: string;
}

/** One nominal timepoint of the composite profile. */
export interface CompositeTimepoint {
  readonly time: number;
  /** Arithmetic mean of the imputed concentrations (for display). */
  readonly mean: number;
  /** Sample SD (n−1 denominator); 0 when n < 2. */
  readonly sd: number;
  /** Coefficient of variation (%), `NaN` when mean ≤ 0. */
  readonly cvPct: number;
  /** Animals sampled at this nominal time. */
  readonly n: number;
  /** BLQ observations at this nominal time. */
  readonly nBlq: number;
}

export interface CompositeProfileResult {
  readonly timepoints: ReadonlyArray<CompositeTimepoint>;
}

export interface SparseAucResult {
  /** Composite AUClast (linear trapezoidal, anchored at t=0, C=0). */
  readonly auc: number;
  /** Holder standard error of the AUC. */
  readonly se: number;
  /** Nedelman-Jia correlated Satterthwaite degrees of freedom. */
  readonly df: number;
  /** Student-t two-sided CI at `ciLevel`. */
  readonly ci: readonly [number, number];
  /** Topology derived from the r_ij matrix. */
  readonly topology: SamplingTopology;
  /** The composite profile (mean/SD/%CV/n/%BLQ per nominal time). */
  readonly composite: CompositeProfileResult;
  readonly warnings: ReadonlyArray<SparseWarning>;
}

// ──────────────────────────────────────────────────────────────────────────
// Internal: grouping by nominal time
// ──────────────────────────────────────────────────────────────────────────

/** One nominal timepoint's imputed samples + animal ids, post BLQ rule. */
interface TimeGroup {
  readonly time: number;
  /** Imputed concentrations (BLQ rule applied). */
  readonly conc: number[];
  /** Animal ids aligned with conc (synthetic unique ids when animalId is null). */
  readonly animals: number[];
  readonly nBlq: number;
}

/**
 * Per-record BLQ imputation for the composite pool (AD-5, impute-before-average).
 *
 * Deliberately NOT `applyBlqStrategy` from `blq.ts`: that function's rules are
 * keyed on a single profile's *phase* (preFirstMeasurable / embedded / afterLast),
 * which is a per-profile-time-ordering concept. The sparse composite pools many
 * animals at ONE nominal timepoint — there is no within-record phase, so the
 * phase model does not apply. We reuse `blq.ts`'s `BlqRule` **type** and its exact
 * set-zero / set-half-lloq / exclude / missing **semantics** at the pooling step;
 * sharing the phase machinery would force an ill-fitting model. (Reviewer note,
 * UC-04 build: this is a justified deviation from the §1a "reuse blq.ts" anchor.)
 */
function imputeBlq(conc: number, blq: number, lloqVal: number, rule: BlqRule): number | null {
  if (blq === 0) return conc;
  switch (rule) {
  case 'set-zero': return 0;
  case 'set-half-lloq': return lloqVal / 2;
  case 'missing': return NaN;
  case 'exclude': return null; // drop the row from the pool
  }
}

function groupByTime(input: SparseInput, rule: BlqRule): TimeGroup[] {
  const n = input.nominalTime.length;
  const map = new Map<number, {conc: number[]; animals: number[]; nBlq: number}>();
  let syntheticAnimal = -1;
  for (let i = 0; i < n; i++) {
    const t = input.nominalTime[i];
    if (!Number.isFinite(t)) continue;
    const lloqVal = typeof input.lloq === 'number' ? input.lloq : input.lloq[i];
    const imputed = imputeBlq(input.conc[i], input.blqMask[i], lloqVal, rule);
    let g = map.get(t);
    if (g === undefined) {g = {conc: [], animals: [], nBlq: 0}; map.set(t, g);}
    if (input.blqMask[i] === 1) g.nBlq++;
    if (imputed === null) continue; // excluded record
    g.conc.push(imputed);
    // A null animalId column means destructive — each row is a distinct animal.
    g.animals.push(input.animalId !== null ? input.animalId[i] : syntheticAnimal--);
  }
  return [...map.entries()]
    .sort((a, b) => a[0] - b[0])
    .map(([time, g]) => ({time, conc: g.conc, animals: g.animals, nBlq: g.nBlq}));
}

// ──────────────────────────────────────────────────────────────────────────
// Internal: per-timepoint statistics
// ──────────────────────────────────────────────────────────────────────────

function mean(xs: ReadonlyArray<number>): number {
  if (xs.length === 0) return NaN;
  let s = 0;
  for (const x of xs) s += x;
  return s / xs.length;
}

/** Sample variance (n−1). Returns 0 when n < 2 (caller may model it — AD-7). */
function sampleVar(xs: ReadonlyArray<number>, m: number): number {
  const n = xs.length;
  if (n < 2) return 0;
  let ss = 0;
  for (const x of xs) {const d = x - m; ss += d * d;}
  return ss / (n - 1);
}

// ──────────────────────────────────────────────────────────────────────────
// Internal: lever-arm weights (PKNCA sparse_auc_weight_linear)
// ──────────────────────────────────────────────────────────────────────────

/**
 * Lever-arm linear-trapezoidal weights over `times` (which must start with the
 * t=0 anchor): `w = c(0, diff/2) + c(diff/2, 0)`. The returned weight at index i
 * is the AUC coefficient on the mean at `times[i]`.
 */
function leverArmWeights(times: ReadonlyArray<number>): number[] {
  const n = times.length;
  const w = new Array<number>(n).fill(0);
  for (let i = 0; i < n; i++) {
    const lo = i > 0 ? (times[i] - times[i - 1]) / 2 : 0;
    const hi = i < n - 1 ? (times[i + 1] - times[i]) / 2 : 0;
    w[i] = lo + hi;
  }
  return w;
}

// ──────────────────────────────────────────────────────────────────────────
// Public: composite profile builder
// ──────────────────────────────────────────────────────────────────────────

/**
 * Build the composite mean concentration-time profile from sparse observations:
 * arithmetic mean / sample-SD / %CV / n / %BLQ per nominal time (AD-4 arithmetic
 * mean is load-bearing for the closed-form variance). BLQ is imputed per
 * `blqRule` **before** pooling (AD-5).
 */
export function buildCompositeProfile(
  input: SparseInput, blqRule: BlqRule = 'set-zero',
): CompositeProfileResult {
  const groups = groupByTime(input, blqRule);
  const timepoints: CompositeTimepoint[] = groups.map((g) => {
    const m = mean(g.conc);
    const v = sampleVar(g.conc, m);
    const sd = Math.sqrt(v);
    return {
      time: g.time,
      mean: m,
      sd,
      cvPct: m > 0 ? 100 * sd / m : NaN,
      n: g.conc.length,
      nBlq: g.nBlq,
    };
  });
  return {timepoints};
}

// ──────────────────────────────────────────────────────────────────────────
// Public: design-aware closed-form sparse AUC + SE + df + CI
// ──────────────────────────────────────────────────────────────────────────

/**
 * Compute the composite sparse AUClast with its Holder standard error,
 * Nedelman-Jia degrees of freedom, and Student-t confidence interval.
 *
 * Validated against PKNCA 0.12.1 `pk.calc.sparse_auclast` (destructive) and a
 * hand-derived Holder oracle (batch) — see `__tests__/fixtures/05_mouse_sparse.json`.
 */
export function sparseAuc(input: SparseInput, options: SparseAucOptions = {}): SparseAucResult {
  const ciLevel = options.ciLevel ?? 0.95;
  const blqRule = options.blqRule ?? 'set-zero';
  const varianceBorrow = options.varianceBorrow ?? true;
  const warnings: SparseWarning[] = [];

  const groups = groupByTime(input, blqRule);
  const composite = buildCompositeProfile(input, blqRule);
  const K = groups.length;

  // --- Per-timepoint mean, n, sample variance.
  const means = groups.map((g) => mean(g.conc));
  const r = groups.map((g) => g.conc.length);
  const variances = groups.map((g, i) => sampleVar(g.conc, means[i]));

  // --- AD-7: variance-borrowing for n=1 timepoints (Nedelman var-vs-mean fit).
  if (varianceBorrow) modelSingletonVariances(means, r, variances, warnings);

  // --- r_ij overlap matrix + topology (AD-2 / AD-11).
  const rij = overlapMatrix(groups);
  const topology = deriveTopology(r, rij, K);
  if (input.animalId === null) {
    warnings.push({
      code: 'SPARSE_ANIMAL_ID_ABSENT_ASSUMED_DESTRUCTIVE',
      message: 'No animal-ID column: the design is assumed destructive (Bailer, ' +
        'r_ij = 0). If samples are actually paired across timepoints (batch/serial), ' +
        'the closed-form SE is too narrow — supply animal IDs.',
    });
  } else if (options.declaredTopology !== undefined && options.declaredTopology !== topology) {
    warnings.push({
      code: 'SPARSE_TOPOLOGY_MISMATCH',
      message: `Declared sampling type '${options.declaredTopology}' but the animal × ` +
        `nominal-time table is '${topology}'. Using the data-derived '${topology}'.`,
    });
  }

  // --- Mean profile for AUClast with the >50%-BLQ-zeroing rule.
  const meansForAuc = means.map((m, i) => (r[i] > 0 && groups[i].nBlq / r[i] > 0.5 ? 0 : m));

  // --- Anchored timepoint vector (t=0, C=0) + lever-arm weights.
  const times = [0, ...groups.map((g) => g.time)];
  const wAll = leverArmWeights(times);
  const w = wAll.slice(1); // weight on each measured timepoint mean

  // --- AUClast = linear trapezoidal on the anchored mean profile to the last
  // measurable (non-zero) timepoint.
  const auc = aucLast([0, ...meansForAuc], times);

  // --- Holder A1 covariance matrix of the timepoint means (Σ).
  const Sigma = covarianceMatrix(groups, means, r, variances, rij);

  // --- ψ = wᵀ Σ w  (A1 variance).
  let psi = 0;
  for (let i = 0; i < K; i++)
    for (let j = 0; j < K; j++) psi += w[i] * w[j] * Sigma[i][j];

  const se = Math.sqrt(Math.max(psi, 0));

  // --- Nedelman-Jia correlated Satterthwaite df.
  const df = satterthwaiteDf(w, Sigma, r, rij, psi);
  if (!Number.isFinite(df) || df <= 0) {
    warnings.push({
      code: 'SPARSE_DF_UNAVAILABLE',
      message: 'Degrees of freedom could not be computed (all timepoint variances ' +
        'are zero); the closed-form CI is unavailable.',
    });
  }

  // --- AD-3: terminal-overestimation heuristic.
  flagTerminalOverest(groups, meansForAuc, options.estimatedHalfLife, warnings);

  // --- Student-t CI.
  let ci: [number, number] = [NaN, NaN];
  if (Number.isFinite(df) && df > 0 && se > 0) {
    const t = studentTInv(1 - (1 - ciLevel) / 2, df);
    ci = [auc - t * se, auc + t * se];
  }

  return {auc, se, df, ci, topology, composite, warnings};
}

// ──────────────────────────────────────────────────────────────────────────
// Internal: AUClast, overlap matrix, covariance, df, topology, guards
// ──────────────────────────────────────────────────────────────────────────

/** Linear-trapezoidal AUC to the last measurable (non-zero) concentration. */
function aucLast(conc: ReadonlyArray<number>, time: ReadonlyArray<number>): number {
  let last = -1;
  for (let i = conc.length - 1; i >= 0; i--)
    if (conc[i] > 0) {last = i; break;}

  if (last <= 0) return 0;
  let auc = 0;
  for (let i = 0; i < last; i++) auc += (time[i + 1] - time[i]) * (conc[i] + conc[i + 1]) / 2;
  return auc;
}

/** r_ij[i][j] = animals sampled at **both** timepoint i and timepoint j. */
function overlapMatrix(groups: TimeGroup[]): number[][] {
  const K = groups.length;
  const sets = groups.map((g) => new Set(g.animals));
  const m = Array.from({length: K}, () => new Array<number>(K).fill(0));
  for (let i = 0; i < K; i++) {
    for (let j = i + 1; j < K; j++) {
      let c = 0;
      const si = sets[i];
      for (const a of sets[j]) if (si.has(a)) c++;
      m[i][j] = c; m[j][i] = c;
    }
  }
  return m;
}

function deriveTopology(r: number[], rij: number[][], K: number): SamplingTopology {
  let anyOverlap = false;
  let allFull = true;
  for (let i = 0; i < K; i++) {
    for (let j = i + 1; j < K; j++) {
      if (rij[i][j] > 0) anyOverlap = true;
      // "serial" ⇔ every pair fully overlaps (the same animals at every time).
      if (rij[i][j] !== Math.min(r[i], r[j]) || r[i] !== r[j]) allFull = false;
    }
  }
  if (!anyOverlap) return 'destructive';
  if (allFull) return 'serial';
  return 'batch';
}

/**
 * Holder A1 + A3 covariance matrix of the timepoint means:
 *   Σ_ii = σ̂_i²/r_i ; Σ_ij = r_ij σ̂_ij/(r_i r_j) with σ̂_ij the unbiased A3.
 */
function covarianceMatrix(
  groups: TimeGroup[], means: number[], r: number[], variances: number[], rij: number[][],
): number[][] {
  const K = groups.length;
  const Sigma = Array.from({length: K}, () => new Array<number>(K).fill(0));
  for (let i = 0; i < K; i++) Sigma[i][i] = r[i] > 0 ? variances[i] / r[i] : 0;

  for (let i = 0; i < K; i++) {
    for (let j = i + 1; j < K; j++) {
      if (rij[i][j] < 1) continue;
      // Animals common to both timepoints, with their paired deviations.
      const idxJ = new Map<number, number>();
      groups[j].animals.forEach((a, k) => idxJ.set(a, k));
      let sumXiXj = 0;
      groups[i].animals.forEach((a, k) => {
        const kj = idxJ.get(a);
        if (kj !== undefined) sumXiXj += (groups[i].conc[k] - means[i]) * (groups[j].conc[kj] - means[j]);
      });
      const rij_ = rij[i][j];
      // A3 unbiased denominator.
      const denom = (rij_ - 1) + (1 - rij_ / r[i]) * (1 - rij_ / r[j]);
      const sij = denom !== 0 ? sumXiXj / denom : 0;
      const covMeans = rij_ * sij / (r[i] * r[j]);
      Sigma[i][j] = covMeans; Sigma[j][i] = covMeans;
    }
  }
  return Sigma;
}

/**
 * Nedelman-Jia correlated Satterthwaite df:
 *   df = ψ² / [ Σ_i (w_i² Σ_ii)²/(r_i−1) + 2 Σ_{i<j,Σ_ij≠0} (w_i w_j Σ_ij)²/(r_ij−1) ]
 * (ν_ii = r_i−1, ν_ij = r_ij−1, floored at 1). Reduces to the scalar
 * Satterthwaite df in the destructive (r_ij = 0) case.
 */
function satterthwaiteDf(
  w: number[], Sigma: number[][], r: number[], rij: number[][], psi: number,
): number {
  const K = w.length;
  let den = 0;
  for (let i = 0; i < K; i++) {
    const term = w[i] * w[i] * Sigma[i][i];
    den += (term * term) / Math.max(r[i] - 1, 1);
  }
  for (let i = 0; i < K; i++) {
    for (let j = i + 1; j < K; j++) {
      if (Sigma[i][j] === 0) continue;
      const term = w[i] * w[j] * Sigma[i][j];
      den += 2 * (term * term) / Math.max(rij[i][j] - 1, 1);
    }
  }
  return den > 0 ? (psi * psi) / den : NaN;
}

/**
 * AD-7 — model the variance of n=1 timepoints from a Nedelman variance-vs-mean
 * power law fitted across the n≥2 timepoints (`ln σ̂² = a + b·ln μ`). Mutates
 * `variances` in place and raises `SPARSE_VARIANCE_MODELED`. No-op when there
 * are no singletons, or fewer than two n≥2 timepoints with a positive mean.
 */
function modelSingletonVariances(
  means: number[], r: number[], variances: number[], warnings: SparseWarning[],
): void {
  const singletons: number[] = [];
  for (let i = 0; i < r.length; i++) if (r[i] === 1) singletons.push(i);
  if (singletons.length === 0) return;

  const lx: number[] = [];
  const ly: number[] = [];
  for (let i = 0; i < r.length; i++) {
    if (r[i] >= 2 && means[i] > 0 && variances[i] > 0) {
      lx.push(Math.log(means[i]));
      ly.push(Math.log(variances[i]));
    }
  }
  if (lx.length < 2) return; // cannot fit — leave variance at 0 (no false modeling)

  // OLS ln σ̂² = a + b·ln μ.
  const nFit = lx.length;
  let sx = 0; let sy = 0; let sxx = 0; let sxy = 0;
  for (let i = 0; i < nFit; i++) {sx += lx[i]; sy += ly[i]; sxx += lx[i] * lx[i]; sxy += lx[i] * ly[i];}
  const denom = nFit * sxx - sx * sx;
  if (denom === 0) return;
  const b = (nFit * sxy - sx * sy) / denom;
  const a = (sy - b * sx) / nFit;

  const modeled: number[] = [];
  for (const i of singletons) {
    if (means[i] > 0) {
      variances[i] = Math.exp(a + b * Math.log(means[i]));
      modeled.push(i);
    }
  }
  if (modeled.length > 0) {
    warnings.push({
      code: 'SPARSE_VARIANCE_MODELED',
      message: `Variance modeled (Nedelman variance-vs-mean power law) at ${modeled.length} ` +
        `singleton timepoint(s) — the CI there reflects a modeled, not observed, variance.`,
    });
  }
}

/**
 * AD-3 — flag a likely AUC over-estimation on a steep, wide-gap terminal phase
 * (the linear rule over-estimates a convex decline). Trigger: terminal gap
 * `(t_K − t_{K−1}) > 2 × t½` AND terminal mean `< 20 % Cmax`. With no t½, fall
 * back to the concentration-only criterion so the flag never disappears silently.
 */
function flagTerminalOverest(
  groups: TimeGroup[], meansForAuc: number[], halfLife: number | undefined,
  warnings: SparseWarning[],
): void {
  // Last measurable index.
  let last = -1;
  for (let i = meansForAuc.length - 1; i >= 0; i--) if (meansForAuc[i] > 0) {last = i; break;}
  if (last < 1) return;
  const cmax = Math.max(...meansForAuc.slice(0, last + 1));
  if (!(cmax > 0)) return;
  const terminalLow = meansForAuc[last] < 0.2 * cmax;
  if (!terminalLow) return;
  const gap = groups[last].time - groups[last - 1].time;
  const wideGap = halfLife !== undefined && halfLife > 0 ? gap > 2 * halfLife : true;
  if (wideGap) {
    warnings.push({
      code: 'SPARSE_TERMINAL_OVEREST',
      message: 'The terminal phase is steep with a wide last-interval gap; the ' +
        'linear-trapezoidal rule may over-estimate AUC over a convex decline ' +
        (halfLife !== undefined ? '(gap > 2×t½, ' : '(no t½ supplied; ') +
        'terminal mean < 20% Cmax). Interpret AUClast with care.',
    });
  }
}
