/**
 * Subspace minimisation for L-BFGS-B, with the Morales–Nocedal (2011)
 * project-and-backtrack refinement.
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
 * Morales–Nocedal 2011 fix: if `xᶜ + Z d̄ᵘ` violates a bound, instead
 * of truncating along the straight line (the 1997 code) we project
 * and backtrack in λ ∈ {1, ½, ¼, …} until the quadratic model
 * strictly decreases. If no backtrack succeeds we fall back to `xᶜ`.
 *
 * References:
 *  - Byrd, Lu, Nocedal, Zhu (1995) §5.
 *  - Morales, Nocedal (2011) — the `c-jlm-jn` fix in `subsm.f`.
 */
import {BOUND_LOWER, BOUND_BOTH, BOUND_UPPER} from './types';
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
/*  Model evaluation helper                                            */
/* ================================================================== */

/**
 * Compute `Δm = m_k(point) − m_k(ref)` using
 *   Δm = δᵀ(g + B z) + ½ δᵀ B δ,
 *   δ = point − ref,   z = ref − xₖ.
 *
 * `bz = B z` and `h = g + B z` must be precomputed by the caller.
 */
function modelDelta(
  point: Float64Array,
  ref: Float64Array,
  h: Float64Array,
  mat: BFGSMat,
  ws: SubspaceWorkspace,
): number {
  const n = point.length;
  const {delta, bDelta, scratch2m, scratchN} = ws;
  for (let i = 0; i < n; i++) delta[i] = point[i] - ref[i];
  mat.applyB(delta, bDelta, scratch2m, scratchN);
  let dot1 = 0;
  let dot2 = 0;
  for (let i = 0; i < n; i++) {
    dot1 += delta[i] * h[i];
    dot2 += delta[i] * bDelta[i];
  }
  return dot1 + 0.5 * dot2;
}

/* ================================================================== */
/*  Subspace minimisation                                              */
/* ================================================================== */

/**
 * Run subspace minimisation from `(xᶜ, c)` produced by the Cauchy sweep.
 *
 * On success writes the refined candidate into `ws.xHat`. If the
 * unconstrained subspace minimiser is infeasible and no backtrack
 * produces descent, `xHat` is set to `xᶜ` and `improved = false`.
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
  maxBacktrack = 20,
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

  const {r, du, xHat, mc, nr, mnr, v, nntCol, kmat, piv, wMc, h, bz,
    scratch2m: sm2, scratchN: smN} = ws;

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

  // ---- (7) Try the unconstrained subspace step ---------------------
  let feasible = true;
  for (let i = 0; i < n; i++) {
    xHat[i] = xc[i] + du[i];
    const code = nbd[i];
    if (code === BOUND_LOWER || code === BOUND_BOTH)
      if (xHat[i] < lower[i]) {feasible = false; break;}

    if (code === BOUND_UPPER || code === BOUND_BOTH)
      if (xHat[i] > upper[i]) {feasible = false; break;}
  }
  if (feasible) return {improved: true};

  // ---- (8) Morales–Nocedal backtrack -------------------------------
  // Precompute h = g + B·z where z = xc − x_k.
  for (let i = 0; i < n; i++) smN[i] = xc[i] - x[i];
  mat.applyB(smN, bz, sm2, h); // h is a scratch for applyB's inner use
  // h gets used as applyB's scratchN above — overwrite it with g + bz.
  for (let i = 0; i < n; i++) h[i] = g[i] + bz[i];

  let lambda = 1;
  for (let step = 0; step < maxBacktrack; step++) {
    // x̂ = project(xc + λ · du, l, u). Non-free coords have du=0, so they
    // stay at xc (which is already feasible).
    for (let i = 0; i < n; i++) {
      let xi = xc[i] + lambda * du[i];
      const code = nbd[i];
      if (code === BOUND_LOWER || code === BOUND_BOTH) {
        const lo = lower[i];
        if (xi < lo) xi = lo;
      }
      if (code === BOUND_UPPER || code === BOUND_BOTH) {
        const hi = upper[i];
        if (xi > hi) xi = hi;
      }
      xHat[i] = xi;
    }
    const dm = modelDelta(xHat, xc, h, mat, ws);
    if (dm < 0) return {improved: true};
    lambda *= 0.5;
  }

  // Exhausted backtracks → fall back to xᶜ.
  for (let i = 0; i < n; i++) xHat[i] = xc[i];
  return {improved: false};
}
