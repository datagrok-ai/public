/**
 * Subspace minimisation for L-BFGS-B with the Morales–Nocedal (2011)
 * endpoint-selection refinement.
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
 * `subsm.f` lines 51–80):
 *   1. Form the unconstrained Newton-like point `x̂ = xᶜ + Z d̄ᵘ`.
 *   2. Project once onto the box: `x̄ = P_{[l, u]}(x̂)`.
 *   3. **Directional-derivative test** on the original gradient at
 *      `xₖ`: accept `x̄` iff `gₖ · (x̄ − xₖ) < 0` and `‖x̄ − xₖ‖² > 0`.
 *      No angle threshold or normalisation; the outer-loop line
 *      search handles step sizing. The `dd > 0` clause is a TS-only
 *      guard against an FP-cancelled `gd < 0` with `x̄ ≡ xₖ`; it is a
 *      strict subset of the canonical Fortran condition (Fortran's
 *      `dcsrch` absorbs the same case at higher cost).
 *   4. **Truncation fallback** (1997 scheme): back-track from `xᶜ`
 *      along `d̄ᵘ`, only on free coordinates. The maximum `α ∈ [0, 1]`
 *      that keeps `xᶜ + α · d̄ᵘ` feasible is the truncation factor.
 *      On `α < 1` the binding coordinate is snapped to its bound to
 *      absorb ε-drift from `room / dk`.
 *   5. If `α = 0`, fall back to `xᶜ`.
 *
 * The 1997 code truncated unconditionally; the 2011 paper showed
 * that path can produce a search direction `d` that is nearly
 * orthogonal to `−g`, stalling the outer line search. The
 * directional-derivative test prefers the projected point whenever
 * it is descending and only falls back to truncation when it is
 * not. Crucially the truncation reference is `xᶜ`, not `xₖ` — when
 * the Cauchy phase has moved several variables to bounds the two
 * differ and only the former matches Fortran v3.0.
 *
 * References:
 *  - Byrd, Lu, Nocedal, Zhu (1995) §5.
 *  - Zhu, Byrd, Lu, Nocedal (1997) — Algorithm 778.
 *  - Morales, Nocedal (2011) — Remark on Algorithm 778 (the
 *    `c-jlm-jn` fix in `subsm.f`).
 */
import {BOUND_FREE, BOUND_LOWER, BOUND_BOTH, BOUND_UPPER} from './types';
import type {BFGSMat} from './bfgs-mat';

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
 * unconstrained Newton-like point onto the box and prefer it when
 * the directional-derivative test `gₖ · (x̄ − xₖ) < 0` holds;
 * otherwise fall back to the 1997 line-truncation scheme along
 * `xᶜ → xᶜ + d̄ᵘ` on free coordinates. If truncation gives `α = 0`,
 * `xHat` is set to `xᶜ` and `improved = false`.
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

  // ---- (9) Directional-derivative test on the original gradient
  //          (subsm.f c-jlm-jn lines 51–62): accept x̄ iff
  //          gₖ · (x̄ − xₖ) < 0. No angle/normalisation; the outer
  //          line search handles step sizing. The `dd > 0` clause
  //          guards against an FP-cancelled gd < 0 with x̄ ≡ xₖ.
  let gd = 0;
  let dd = 0;
  for (let i = 0; i < n; i++) {
    const di = xBar[i] - x[i];
    gd += g[i] * di;
    dd += di * di;
  }

  if (gd < 0 && dd > 0) {
    for (let i = 0; i < n; i++) xHat[i] = xBar[i];
    return {improved: true};
  }

  // ---- (10) 1997 truncation fallback (subsm.f c-jlm-jn lines 56–80).
  //          Back-track from xᶜ along d̄ᵘ, only on free coordinates. The
  //          maximum α ∈ [0, 1] that keeps xᶜ + α·d̄ᵘ feasible is the
  //          truncation factor. On α < 1 the first-hitting coordinate
  //          is snapped to its bound to absorb the ε-rounding from
  //          `room / dk`.
  let alpha = 1;
  let iexit = -1;
  let iexitDir = 0;
  for (let kk = 0; kk < t; kk++) {
    const i = freeSet[kk];
    const code = nbd[i];
    if (code === BOUND_FREE) continue;
    const dk = du[i];
    if (dk === 0) continue;

    let cap = alpha;
    if (dk < 0 && (code === BOUND_LOWER || code === BOUND_BOTH)) {
      const room = lower[i] - xc[i]; // ≤ 0 on feasible xᶜ
      if (room >= 0) cap = 0;
      else if (dk * alpha < room) cap = room / dk;
    } else if (dk > 0 && (code === BOUND_UPPER || code === BOUND_BOTH)) {
      const room = upper[i] - xc[i]; // ≥ 0 on feasible xᶜ
      if (room <= 0) cap = 0;
      else if (dk * alpha > room) cap = room / dk;
    }
    if (cap < alpha) {
      alpha = cap;
      iexit = i;
      iexitDir = dk;
    }
  }

  if (alpha === 0) {
    // ---- (11) Truncation degenerate → fall back to xᶜ.
    for (let i = 0; i < n; i++) xHat[i] = xc[i];
    return {improved: false};
  }

  // Build xHat: non-free coords stay at xᶜ; free coords advance by
  // α·d̄ᵘ, except the binding coord `iexit` which is snapped to its
  // bound (mirrors subsm.f's `if (alpha .lt. one)` block: x(iexit)
  // is set to the bound, d(iexit) is zeroed; we achieve the same by
  // writing the bound into xHat directly and skipping the standard
  // `xc + α·du` formula on that coord).
  for (let i = 0; i < n; i++) xHat[i] = xc[i];
  for (let kk = 0; kk < t; kk++) {
    const i = freeSet[kk];
    if (i === iexit && alpha < 1)
      xHat[i] = iexitDir > 0 ? upper[i] : lower[i];
    else
      xHat[i] = xc[i] + alpha * du[i];
  }
  return {improved: true};
}
