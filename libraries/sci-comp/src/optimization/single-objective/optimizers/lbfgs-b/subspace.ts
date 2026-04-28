/**
 * Subspace minimisation for L-BFGS-B with the Morales–Nocedal (2011)
 * project-and-angle-test refinement.
 *
 * Given the Cauchy point `xᶜ` and free set `F = {i : lᵢ < xᶜᵢ < uᵢ}`,
 * we approximately solve
 *
 *    min_{d̄ ∈ ℝᵗ}  r̄ᶜᵀd̄ + ½ d̄ᵀ B̄ d̄,
 *    B̄ = θ Iₜ − Nᵀ M N,       N = Wᵀ Z ∈ ℝ²ᵐˣᵗ,
 *
 * where `Z` is the selector of free coordinates (`t = |F|`). The
 * reduced gradient at `xᶜ` (BLNZ'95 eq. 5.4) is
 *
 *    r̄ᶜ = Zᵀ [ g + θ(xᶜ − xₖ) − W M c ].
 *
 * The unconstrained subspace step is obtained in closed form via
 * Sherman–Morrison–Woodbury,
 *
 *    d̄ᵘ = −(1/θ) r̄ᶜ − (1/θ²) Nᵀ (I₂ₘ − (1/θ) M N Nᵀ)⁻¹ M N r̄ᶜ.
 *
 * The 2m×2m middle matrix is factored by LU with partial pivoting
 * (it is not generally SPD).
 *
 * Endpoint selection (Morales–Nocedal 2011, `c-jlm-jn` block in
 * `subsm.f`):
 *   1. Form the unconstrained Newton-like point `x̂ = xᶜ + Z d̄ᵘ`.
 *   2. Project once onto the box: `x̄ = P_{[l, u]}(x̂)`.
 *   3. Form the direction back to the iterate: `d_k = x̄ − x_k`.
 *   4. **Angle test** on the true gradient:
 *         `gₖᵀ d_k ≤ −η · ‖d_k‖ · ‖gₖ‖`.
 *      If passing, return `xHat = x̄`; the outer-loop line search
 *      then operates on `d = xHat − x_k`.
 *   5. **Truncation fallback** (1997 scheme): otherwise find the
 *      maximum `α ∈ (0, 1]` such that `x_k + α(x̂ − x_k)` stays
 *      feasible, and return that point. If `α = 0`, fall back to xᶜ.
 *
 * The 1997 code truncated unconditionally and the 2011 paper showed
 * that path can produce a search direction `d` that is nearly
 * orthogonal to `−g`, stalling the outer line search. The angle
 * test prefers the projected point whenever it is meaningfully
 * descending and only falls back to truncation when it is not.
 *
 * References:
 *  - Byrd, Lu, Nocedal, Zhu (1995) §5.
 *  - Zhu, Byrd, Lu, Nocedal (1997) — Algorithm 778.
 *  - Morales, Nocedal (2011) — Remark on Algorithm 778 (the
 *    `c-jlm-jn` fix in `subsm.f`).
 */
import {BOUND_LOWER, BOUND_BOTH, BOUND_UPPER} from './types';
import type {BFGSMat} from './bfgs-mat';

/**
 * Strong-descent / angle-test threshold from `subsm.f`. The
 * projected endpoint is preferred when its direction back to xₖ
 * makes an angle with `−gₖ` whose cosine is at least η. Calibrated
 * against bounded-benchmarks; see specs-n-plans/lbfgs-b-baseline.
 */
const ETA = 1e-2;

/* ================================================================== */
/*  Workspace                                                          */
/* ================================================================== */

export interface SubspaceWorkspace {
  /** Reduced gradient; `r[0..freeCount)` is the dense r̄ᶜ. */
  r: Float64Array; // length n (oversized)
  /** Subspace direction, extended to full `n` (0 on non-free coords). */
  du: Float64Array; // length n
  /** Output. */
  xHat: Float64Array; // length n

  // 2m-sized scratches
  mc: Float64Array;
  nr: Float64Array;
  mnr: Float64Array;
  v: Float64Array;
  nntCol: Float64Array;

  /** `2m × 2m` middle matrix (column-major). */
  kmat: Float64Array;
  piv: Int32Array;

  /** `n`-sized scratches used inside `modelValue` / `applyB`. */
  wMc: Float64Array;
  delta: Float64Array;
  bDelta: Float64Array;
  bz: Float64Array;
  h: Float64Array;
  scratchN: Float64Array;
  scratch2m: Float64Array;
}

export function makeSubspaceWorkspace(n: number, m: number): SubspaceWorkspace {
  const twoM = 2 * m;
  return {
    r: new Float64Array(n),
    du: new Float64Array(n),
    xHat: new Float64Array(n),
    mc: new Float64Array(twoM),
    nr: new Float64Array(twoM),
    mnr: new Float64Array(twoM),
    v: new Float64Array(twoM),
    nntCol: new Float64Array(twoM),
    kmat: new Float64Array(twoM * twoM),
    piv: new Int32Array(twoM),
    wMc: new Float64Array(n),
    delta: new Float64Array(n),
    bDelta: new Float64Array(n),
    bz: new Float64Array(n),
    h: new Float64Array(n),
    scratchN: new Float64Array(n),
    scratch2m: new Float64Array(twoM),
  };
}

export interface SubspaceResult {
  /** `true` if subspace min produced a distinct `xHat` (caller uses it
   *  as the outer-loop search direction endpoint). `false` when we
   *  fell back to `xᶜ` (caller should still take the Cauchy step). */
  improved: boolean;
}

/* ================================================================== */
/*  LU with partial pivoting (private)                                 */
/* ================================================================== */

/**
 * Solve `A x = b` with `A` k×k column-major; mutates `A` (stores LU
 * factors in place). Returns `false` on singular pivot. `x` may alias
 * `b`; on entry `x` must contain `b`.
 */
function luSolveInPlace(
  A: Float64Array,
  x: Float64Array,
  k: number,
  stride: number,
  piv: Int32Array,
): boolean {
  for (let i = 0; i < k; i++) piv[i] = i;

  for (let c = 0; c < k; c++) {
    // Find pivot (max |A[i, c]| in column c, i ≥ c).
    let maxRow = c;
    let maxVal = Math.abs(A[c + c * stride]);
    for (let i = c + 1; i < k; i++) {
      const v = Math.abs(A[i + c * stride]);
      if (v > maxVal) {
        maxVal = v;
        maxRow = i;
      }
    }
    if (maxVal === 0) return false;

    if (maxRow !== c) {
      for (let j = 0; j < k; j++) {
        const tmp = A[c + j * stride];
        A[c + j * stride] = A[maxRow + j * stride];
        A[maxRow + j * stride] = tmp;
      }
      const tx = x[c]; x[c] = x[maxRow]; x[maxRow] = tx;
      const tp = piv[c]; piv[c] = piv[maxRow]; piv[maxRow] = tp;
    }

    const pivot = A[c + c * stride];
    for (let i = c + 1; i < k; i++) {
      const factor = A[i + c * stride] / pivot;
      A[i + c * stride] = factor;
      for (let j = c + 1; j < k; j++)
        A[i + j * stride] -= factor * A[c + j * stride];
      x[i] -= factor * x[c];
    }
  }

  for (let i = k - 1; i >= 0; i--) {
    let s = x[i];
    for (let j = i + 1; j < k; j++) s -= A[i + j * stride] * x[j];
    x[i] = s / A[i + i * stride];
  }
  return true;
}

/* ================================================================== */
/*  Subspace minimisation                                              */
/* ================================================================== */

/**
 * Run subspace minimisation from `(xᶜ, c)` produced by the Cauchy sweep.
 *
 * On success writes the refined candidate into `ws.xHat`. The
 * endpoint selection follows Morales–Nocedal 2011: project the
 * unconstrained Newton-like point onto the box and prefer it as
 * long as the angle test on the true gradient passes; otherwise
 * fall back to the 1997 line-truncation scheme along `xₖ → x̂`. If
 * truncation gives `α = 0`, `xHat` is set to `xᶜ` and `improved = false`.
 */
export function subspaceMin(
  x: Float64Array,
  g: Float64Array,
  lower: Float64Array,
  upper: Float64Array,
  nbd: Uint8Array,
  xc: Float64Array,
  c: Float64Array,
  freeSet: Int32Array,
  freeCount: number,
  mat: BFGSMat,
  ws: SubspaceWorkspace,
): SubspaceResult {
  const n = x.length;
  const col = mat.col;
  const theta = mat.theta;
  const t = freeCount;
  const twoCol = 2 * col;

  // Default: xHat = xc.
  for (let i = 0; i < n; i++) ws.xHat[i] = xc[i];

  // Nothing to refine: no free variables or no stored pairs.
  if (t === 0 || col === 0) return {improved: false};

  const {r, du, xHat, mc, nr, mnr, v, nntCol, kmat, piv, wMc} = ws;
  // Reuse existing n-sized scratches (delta, bDelta) for the new
  // endpoint-selection logic; the old Morales–Nocedal backtrack
  // helpers are no longer used.
  const xHatFull = ws.delta;
  const xBar = ws.bDelta;

  // ---- (1) Reduced gradient r̄ᶜ on free coords ---------------------
  for (let k = 0; k < twoCol; k++) mc[k] = c[k];
  mat.solveM(mc);
  mat.applyW(mc, wMc);
  for (let k = 0; k < t; k++) {
    const i = freeSet[k];
    r[k] = g[i] + theta * (xc[i] - x[i]) - wMc[i];
  }

  // ---- (2) N · r̄ᶜ  (length 2m) ------------------------------------
  for (let j = 0; j < twoCol; j++) nr[j] = 0;
  for (let k = 0; k < t; k++) {
    const i = freeSet[k];
    const rk = r[k];
    for (let j = 0; j < col; j++) {
      const jP = mat.physOf(j);
      nr[j] += mat.yStore[i + jP * n] * rk;
      nr[col + j] += theta * mat.sStore[i + jP * n] * rk;
    }
  }

  // ---- (3) M · (N · r̄ᶜ) --------------------------------------------
  for (let k = 0; k < twoCol; k++) mnr[k] = nr[k];
  mat.solveM(mnr);

  // ---- (4) Build K = I − (1/θ) M N Nᵀ, column by column ------------
  for (let a = 0; a < twoCol; a++) {
    // Column a of N Nᵀ: (N Nᵀ)[b, a] = Σ_k W[i_k, b] · W[i_k, a].
    for (let b = 0; b < twoCol; b++) nntCol[b] = 0;
    for (let k = 0; k < t; k++) {
      const i = freeSet[k];
      const Wia = a < col ?
        mat.yStore[i + mat.physOf(a) * n] :
        theta * mat.sStore[i + mat.physOf(a - col) * n];
      for (let b = 0; b < twoCol; b++) {
        const Wib = b < col ?
          mat.yStore[i + mat.physOf(b) * n] :
          theta * mat.sStore[i + mat.physOf(b - col) * n];
        nntCol[b] += Wib * Wia;
      }
    }
    // M · (column a of N Nᵀ).
    mat.solveM(nntCol);
    // K[:, a] = −(1/θ) · (M N Nᵀ)[:, a]; then K[a, a] += 1.
    for (let b = 0; b < twoCol; b++)
      kmat[b + a * twoCol] = -nntCol[b] / theta;
    kmat[a + a * twoCol] += 1;
  }

  // ---- (5) Solve K v = M N r̄ᶜ (in v; mnr copied in) ---------------
  for (let k = 0; k < twoCol; k++) v[k] = mnr[k];
  const okLU = luSolveInPlace(kmat, v, twoCol, twoCol, piv);
  if (!okLU) return {improved: false};

  // ---- (6) d̄ᵘ = −(1/θ) r̄ᶜ − (1/θ²) Nᵀ v  ----------------------
  for (let i = 0; i < n; i++) du[i] = 0;
  for (let k = 0; k < t; k++) {
    const i = freeSet[k];
    let nTv = 0;
    for (let j = 0; j < col; j++) {
      const jP = mat.physOf(j);
      nTv += mat.yStore[i + jP * n] * v[j];
      nTv += theta * mat.sStore[i + jP * n] * v[col + j];
    }
    du[i] = -r[k] / theta - nTv / (theta * theta);
  }

  // ---- (7) Form the unconstrained Newton-like point xHatFull = xc + du.
  //          On non-free coords du is zero, so xHatFull[i] = xc[i] there.
  for (let i = 0; i < n; i++) xHatFull[i] = xc[i] + du[i];

  // ---- (8) Project once onto [l, u].
  for (let i = 0; i < n; i++) {
    let xi = xHatFull[i];
    const code = nbd[i];
    if ((code === BOUND_LOWER || code === BOUND_BOTH) && xi < lower[i]) xi = lower[i];
    if ((code === BOUND_UPPER || code === BOUND_BOTH) && xi > upper[i]) xi = upper[i];
    xBar[i] = xi;
  }

  // ---- (9) Angle test on the true gradient: is x̄ − xₖ a direction
  //          of strong descent on the original problem?
  let gd = 0;
  let dd = 0;
  let gg = 0;
  for (let i = 0; i < n; i++) {
    const di = xBar[i] - x[i];
    gd += g[i] * di;
    dd += di * di;
    gg += g[i] * g[i];
  }
  const dNorm = Math.sqrt(dd);
  const gNorm = Math.sqrt(gg);

  if (dNorm > 0 && gNorm > 0 && gd <= -ETA * dNorm * gNorm) {
    // Strong descent on the original f → endpoint = x̄.
    for (let i = 0; i < n; i++) xHat[i] = xBar[i];
    return {improved: true};
  }

  // ---- (10) Truncation fallback (1997): largest α ∈ (0, 1] keeping
  //          xₖ + α · (x̂ − xₖ) feasible. The model-improving direction
  //          is x̂ − xₖ; we just stop at the first hitting bound.
  let alpha = 1;
  for (let i = 0; i < n; i++) {
    const stepDir = xHatFull[i] - x[i];
    if (stepDir > 0 && (nbd[i] === BOUND_UPPER || nbd[i] === BOUND_BOTH)) {
      const cap = (upper[i] - x[i]) / stepDir;
      if (cap < alpha) alpha = Math.max(0, cap);
    } else if (stepDir < 0 && (nbd[i] === BOUND_LOWER || nbd[i] === BOUND_BOTH)) {
      const cap = (lower[i] - x[i]) / stepDir;
      if (cap < alpha) alpha = Math.max(0, cap);
    }
  }

  if (alpha > 0) {
    for (let i = 0; i < n; i++)
      xHat[i] = x[i] + alpha * (xHatFull[i] - x[i]);
    return {improved: true};
  }

  // ---- (11) α = 0 → fall back to xᶜ.
  for (let i = 0; i < n; i++) xHat[i] = xc[i];
  return {improved: false};
}
