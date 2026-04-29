/**
 * Box-bounds helpers for L-BFGS-B.
 *
 * Internal representation mirrors Nocedal's Fortran convention:
 *   - `lower`, `upper`: `Float64Array(n)` with `±Infinity` marking the
 *     respective side as unbounded;
 *   - `nbd`: `Uint8Array(n)` encoding the bound type per variable
 *     (0 = Free, 1 = Lower only, 2 = Both, 3 = Upper only).
 *
 * The user-facing `LBFGSBBounds` object (scalar or per-dim input for
 * either side, both optional) is translated into this triple exactly
 * once per call by {@link normalizeBounds}.
 *
 * Reference: Byrd, Lu, Nocedal, Zhu (1995) §8; L-BFGS-B v3.0 `cmprlb`.
 */
import type {LBFGSBBounds, BoundType} from './types';
import {BOUND_FREE, BOUND_LOWER, BOUND_BOTH, BOUND_UPPER} from './types';

/* ================================================================== */
/*  Normalisation from user input                                      */
/* ================================================================== */

export interface NormalizedBounds {
  /** Per-dimension lower bound; `-Infinity` means unbounded below. */
  lower: Float64Array;
  /** Per-dimension upper bound; `+Infinity` means unbounded above. */
  upper: Float64Array;
  /** Per-dimension {@link BoundType} encoding. */
  nbd: Uint8Array;
  /** True iff at least one variable has any finite bound. */
  anyFinite: boolean;
}

/**
 * Convert the public `{lower?, upper?}` input (scalar or per-dimension)
 * into the internal `NormalizedBounds` triple. Throws on length mismatch
 * or when `lower[i] > upper[i]`.
 *
 * A missing side defaults to `∓Infinity`. When the entire input is
 * absent, all variables are Free.
 */
export function normalizeBounds(
  input: LBFGSBBounds | undefined,
  n: number,
): NormalizedBounds {
  const lower = new Float64Array(n);
  const upper = new Float64Array(n);
  lower.fill(-Infinity);
  upper.fill(+Infinity);

  if (input) {
    fillSide(lower, input.lower, n, 'lower');
    fillSide(upper, input.upper, n, 'upper');
  }

  // Validate & detect at least one finite bound.
  let anyFinite = false;
  for (let i = 0; i < n; i++) {
    const lo = lower[i];
    const hi = upper[i];
    if (lo > hi)
      throw new Error(`L-BFGS-B: bounds[${i}]: lower (${lo}) > upper (${hi})`);
    if (Number.isFinite(lo) || Number.isFinite(hi)) anyFinite = true;
  }

  const nbd = classifyBounds(lower, upper);
  return {lower, upper, nbd, anyFinite};
}

function fillSide(
  out: Float64Array,
  input: number | ArrayLike<number> | undefined,
  n: number,
  side: 'lower' | 'upper',
): void {
  if (input === undefined) return;
  if (typeof input === 'number') {
    if (Number.isNaN(input))
      throw new Error(`L-BFGS-B: bounds.${side} must not be NaN`);
    out.fill(input);
    return;
  }
  if (input.length !== n) {
    throw new Error(
      `L-BFGS-B: bounds.${side}.length (${input.length}) !== x0.length (${n})`,
    );
  }
  for (let i = 0; i < n; i++) {
    const v = input[i];
    if (Number.isNaN(v))
      throw new Error(`L-BFGS-B: bounds.${side}[${i}] is NaN`);
    out[i] = v;
  }
}

/* ================================================================== */
/*  classifyBounds                                                     */
/* ================================================================== */

/**
 * Build the `nbd` encoding from `(lower, upper)` arrays.
 *
 * Per-entry rules (matching Nocedal's `nbd[i]`):
 *   - `lower = -∞`, `upper = +∞` → 0 (Free)
 *   - `lower` finite, `upper = +∞` → 1 (Lower only)
 *   - `lower` finite, `upper` finite → 2 (Both)
 *   - `lower = -∞`, `upper` finite → 3 (Upper only)
 */
export function classifyBounds(lower: Float64Array, upper: Float64Array): Uint8Array {
  const n = lower.length;
  const out = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    const hasLo = Number.isFinite(lower[i]);
    const hasHi = Number.isFinite(upper[i]);
    let code: BoundType;
    if (hasLo && hasHi) code = BOUND_BOTH;
    else if (hasLo) code = BOUND_LOWER;
    else if (hasHi) code = BOUND_UPPER;
    else code = BOUND_FREE;
    out[i] = code;
  }
  return out;
}

/* ================================================================== */
/*  project                                                            */
/* ================================================================== */

/**
 * Project `x` onto the feasible box `[lower, upper]`. Writes into `out`;
 * `x` may alias `out` (in-place projection).
 */
export function project(
  x: Float64Array,
  lower: Float64Array,
  upper: Float64Array,
  nbd: Uint8Array,
  out: Float64Array,
): void {
  const n = x.length;
  for (let i = 0; i < n; i++) {
    let v = x[i];
    const code = nbd[i];
    if (code === BOUND_LOWER || code === BOUND_BOTH) {
      const lo = lower[i];
      if (v < lo) v = lo;
    }
    if (code === BOUND_UPPER || code === BOUND_BOTH) {
      const hi = upper[i];
      if (v > hi) v = hi;
    }
    out[i] = v;
  }
}

/* ================================================================== */
/*  projectedGradient                                                  */
/* ================================================================== */

/**
 * Compute the projected gradient `P(x − g, l, u) − x` and return its
 * ∞-norm.
 *
 * Writes the component-wise projected-gradient vector into `out` (`x`
 * and `out` must not alias). The return value is the primary L-BFGS-B
 * convergence metric (`gradTolerance` in settings).
 */
export function projectedGradient(
  x: Float64Array,
  g: Float64Array,
  lower: Float64Array,
  upper: Float64Array,
  nbd: Uint8Array,
  out: Float64Array,
): number {
  const n = x.length;
  let infNorm = 0;
  for (let i = 0; i < n; i++) {
    let projected = x[i] - g[i];
    const code = nbd[i];
    if (code === BOUND_LOWER || code === BOUND_BOTH) {
      const lo = lower[i];
      if (projected < lo) projected = lo;
    }
    if (code === BOUND_UPPER || code === BOUND_BOTH) {
      const hi = upper[i];
      if (projected > hi) projected = hi;
    }
    const pg = projected - x[i];
    out[i] = pg;
    const abs = pg < 0 ? -pg : pg;
    if (abs > infNorm) infNorm = abs;
  }
  return infNorm;
}

/* ================================================================== */
/*  maxFeasibleStep                                                    */
/* ================================================================== */

/**
 * Compute the largest `α > 0` for which `x + α·d` stays feasible.
 *
 * Per coordinate:
 *   - `d[i] > 0`, upper finite → `(u[i] − x[i]) / d[i]`
 *   - `d[i] < 0`, lower finite → `(l[i] − x[i]) / d[i]`
 *   - otherwise → `+∞`.
 *
 * Returns `+Infinity` if no coordinate caps the step.
 *
 * Numerical guard: any negative ratio (can happen if `x` has drifted
 * slightly infeasible by rounding) is clamped to zero so the caller
 * does not take a backward step.
 */
export function maxFeasibleStep(
  x: Float64Array,
  d: Float64Array,
  lower: Float64Array,
  upper: Float64Array,
  nbd: Uint8Array,
): number {
  const n = x.length;
  let stpmax = Infinity;
  for (let i = 0; i < n; i++) {
    const di = d[i];
    if (di === 0) continue;
    const code = nbd[i];
    let cap: number;
    if (di > 0) {
      if (code !== BOUND_UPPER && code !== BOUND_BOTH) continue;
      cap = (upper[i] - x[i]) / di;
    } else {
      if (code !== BOUND_LOWER && code !== BOUND_BOTH) continue;
      cap = (lower[i] - x[i]) / di;
    }
    if (cap < 0) cap = 0; // x drifted past the bound
    if (cap < stpmax) stpmax = cap;
  }
  return stpmax;
}
