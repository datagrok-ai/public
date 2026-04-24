/**
 * Generalized Cauchy point for L-BFGS-B.
 *
 * The Cauchy point `xᶜ` is the first local minimiser of the quadratic
 * model `mₖ(x) = fₖ + gᵀ(x−xₖ) + ½(x−xₖ)ᵀBₖ(x−xₖ)` along the
 * piecewise-linear projected-gradient path
 *
 *    x(t) = P(xₖ − t gₖ, l, u),     t ≥ 0.
 *
 * Along each linear segment we maintain the 2m-vectors
 *    p = Wᵀ d,     c = Wᵀ z,      z = x(t) − xₖ,
 * and the segment-local directional derivatives (BLNZ'95 eqs. 4.9–4.10)
 *
 *    f'  = gᵀd + θ dᵀz − pᵀM c,
 *    f'' = θ dᵀd − pᵀM p.
 *
 * When variable `b` reaches its bound at breakpoint `t⁽ʲ⁾`, `c, z`
 * are advanced by `Δt · (p, d)`, `xᶜ_b` is snapped, and the
 * rank-1 updates of eqs. 4.11–4.12 are applied to `p, f', f''` in
 * O(m²). Breakpoints are processed in ascending order via a binary
 * min-heap.
 *
 * References:
 *  - R. H. Byrd, P. Lu, J. Nocedal, C. Zhu (1995) §4, Algorithm CP.
 */
import {BOUND_FREE, BOUND_LOWER, BOUND_BOTH, BOUND_UPPER} from './types';
import type {BFGSMat} from './bfgs-mat';

/* ================================================================== */
/*  Binary min-heap keyed by breakpoint value                          */
/* ================================================================== */

/**
 * Allocation-free binary min-heap storing `(t, i)` pairs keyed on `t`.
 * Used by the Cauchy sweep; reused across outer iterations via
 * `clear()`.
 */
export class MinHeap {
  readonly capacity: number;
  readonly ts: Float64Array;
  readonly idx: Int32Array;
  size: number;

  constructor(capacity: number) {
    this.capacity = capacity;
    this.ts = new Float64Array(capacity);
    this.idx = new Int32Array(capacity);
    this.size = 0;
  }

  clear(): void {
    this.size = 0;
  }

  isEmpty(): boolean {
    return this.size === 0;
  }

  push(t: number, i: number): void {
    let k = this.size++;
    this.ts[k] = t;
    this.idx[k] = i;
    while (k > 0) {
      const parent = (k - 1) >> 1;
      if (this.ts[parent] <= this.ts[k]) break;
      this.swap(k, parent);
      k = parent;
    }
  }

  /** Remove and return the min element. Assumes `!isEmpty()`. */
  pop(): {t: number; i: number} {
    const t = this.ts[0];
    const i = this.idx[0];
    this.size--;
    if (this.size > 0) {
      this.ts[0] = this.ts[this.size];
      this.idx[0] = this.idx[this.size];
      let k = 0;
      for (; ;) {
        const left = 2 * k + 1;
        const right = left + 1;
        let smallest = k;
        if (left < this.size && this.ts[left] < this.ts[smallest]) smallest = left;
        if (right < this.size && this.ts[right] < this.ts[smallest]) smallest = right;
        if (smallest === k) break;
        this.swap(k, smallest);
        k = smallest;
      }
    }
    return {t, i};
  }

  private swap(a: number, b: number): void {
    const tt = this.ts[a]; this.ts[a] = this.ts[b]; this.ts[b] = tt;
    const ii = this.idx[a]; this.idx[a] = this.idx[b]; this.idx[b] = ii;
  }
}

/* ================================================================== */
/*  Cauchy workspace                                                   */
/* ================================================================== */

/**
 * Pre-allocated scratch for `cauchyPoint`. All buffers are sized once
 * in `makeCauchyWorkspace(n, m)` and reused across outer iterations.
 */
export interface CauchyWorkspace {
  /** Cauchy point (length `n`). Written by `cauchyPoint`. */
  xc: Float64Array;
  /** Direction along current segment (length `n`). */
  d: Float64Array;
  /** `c = Wᵀ (xᶜ − xₖ)` (length `2m`). Written by `cauchyPoint`. */
  c: Float64Array;
  /** `p = Wᵀ d` maintained along the sweep (length `2m`). */
  p: Float64Array;
  /** `M · p` (length `2m`). */
  mp: Float64Array;
  /** Scratch for `M · c` (length `2m`). */
  mc: Float64Array;
  /** Row of `W` for the variable hitting a bound (length `2m`). */
  wb: Float64Array;
  /** Scratch for `M · w_b` (length `2m`). */
  mwb: Float64Array;
  /** Free-set indices (length `n`). */
  freeSet: Int32Array;
  /** Per-variable breakpoint values (length `n`). */
  t: Float64Array;
  /** Breakpoint min-heap (capacity `n`). */
  heap: MinHeap;
}

export function makeCauchyWorkspace(n: number, m: number): CauchyWorkspace {
  const twoM = 2 * m;
  return {
    xc: new Float64Array(n),
    d: new Float64Array(n),
    c: new Float64Array(twoM),
    p: new Float64Array(twoM),
    mp: new Float64Array(twoM),
    mc: new Float64Array(twoM),
    wb: new Float64Array(twoM),
    mwb: new Float64Array(twoM),
    freeSet: new Int32Array(n),
    t: new Float64Array(n),
    heap: new MinHeap(n),
  };
}

export interface CauchyResult {
  /** Number of free variables (`xᶜ` strictly inside bounds). */
  freeCount: number;
  /**
   * `false` when `f''` became non-positive or non-finite and the sweep
   * fell back to `xᶜ = x` (Cauchy is undefined for that iterate —
   * caller should restart with memory reset).
   */
  ok: boolean;
}

/* ================================================================== */
/*  Cauchy point                                                       */
/* ================================================================== */

const FDP_EPS = Number.EPSILON;

/**
 * Run Algorithm CP to compute `xᶜ` and `c = Wᵀ(xᶜ − xₖ)`. Writes into
 * `ws.xc`, `ws.c`, `ws.freeSet`; returns `freeCount` and an `ok` flag.
 *
 * The workspace's `d`, `p`, `t`, `heap` buffers are used internally
 * and left in an indeterminate state on return.
 */
export function cauchyPoint(
  x: Float64Array,
  g: Float64Array,
  lower: Float64Array,
  upper: Float64Array,
  nbd: Uint8Array,
  mat: BFGSMat,
  ws: CauchyWorkspace,
): CauchyResult {
  const n = x.length;
  const col = mat.col;
  const twoCol = 2 * col;
  const theta = mat.theta;
  const {xc, d, c, p, mp, mc, wb, mwb, t, heap, freeSet} = ws;

  heap.clear();
  for (let i = 0; i < twoCol; i++) c[i] = 0;

  // ---- (1) Breakpoints, initial d, xc ------------------------------
  let fPrime = 0;
  let dTd = 0;
  for (let i = 0; i < n; i++) {
    xc[i] = x[i];
    const code = nbd[i];
    const gi = g[i];
    let ti: number;
    if (code === BOUND_FREE || gi === 0)
      ti = Infinity;
    else if (gi < 0 && (code === BOUND_UPPER || code === BOUND_BOTH))
      ti = (x[i] - upper[i]) / gi;
    else if (gi > 0 && (code === BOUND_LOWER || code === BOUND_BOTH))
      ti = (x[i] - lower[i]) / gi;
    else
      ti = Infinity;

    t[i] = ti;
    if (ti > 0) {
      d[i] = -gi;
      fPrime -= gi * gi;
      dTd += gi * gi;
      if (ti !== Infinity) heap.push(ti, i);
    } else {
      // Variable is on its bound with g pointing into infeasible region —
      // cannot move. xc stays at x[i] (which is the bound).
      d[i] = 0;
    }
  }

  // ---- (2) p = Wᵀ d, Mp, initial f'' -------------------------------
  let fDoublePrime: number;
  if (col > 0) {
    mat.applyWt(d, p);
    for (let k = 0; k < twoCol; k++) mp[k] = p[k];
    mat.solveM(mp);
    let pMp = 0;
    for (let k = 0; k < twoCol; k++) pMp += p[k] * mp[k];
    fDoublePrime = theta * dTd - pMp;
  } else
    fDoublePrime = theta * dTd;


  if (!Number.isFinite(fDoublePrime) || !(fDoublePrime > 0)) {
    const freeCount = buildFreeSet(xc, lower, upper, freeSet, n);
    return {freeCount, ok: false};
  }
  {
    const floor = FDP_EPS * Math.max(1, fDoublePrime);
    if (fDoublePrime < floor) fDoublePrime = floor;
  }

  // ---- (3) Main sweep ----------------------------------------------
  let dtMin = -fPrime / fDoublePrime;
  let tOld = 0;

  if (!heap.isEmpty()) {
    let first = heap.pop();
    let tCurrent = first.t;
    let b = first.i;
    let deltaT = tCurrent - tOld;

    while (dtMin >= deltaT) {
      // Snap variable b to the bound its direction points to.
      xc[b] = d[b] > 0 ? upper[b] : lower[b];
      const zb = xc[b] - x[b];
      tOld = tCurrent;

      // Advance c by Δt · p using the *pre-update* p.
      for (let k = 0; k < twoCol; k++) c[k] += deltaT * p[k];

      const gb = g[b];

      if (col > 0) {
        // Build w_b = row b of W = [Y | θS].
        for (let j = 0; j < col; j++) {
          const jP = mat.physOf(j);
          wb[j] = mat.yStore[b + jP * n];
          wb[col + j] = theta * mat.sStore[b + jP * n];
        }
        // M c and M w_b (mp is still current for the pre-update p).
        for (let k = 0; k < twoCol; k++) mc[k] = c[k];
        mat.solveM(mc);
        for (let k = 0; k < twoCol; k++) mwb[k] = wb[k];
        mat.solveM(mwb);

        let wbMc = 0;
        let wbMp = 0;
        let wbMwb = 0;
        for (let k = 0; k < twoCol; k++) {
          wbMc += wb[k] * mc[k];
          wbMp += wb[k] * mp[k];
          wbMwb += wb[k] * mwb[k];
        }
        fPrime += deltaT * fDoublePrime + gb * gb + theta * gb * zb - gb * wbMc;
        fDoublePrime += -theta * gb * gb - 2 * gb * wbMp - gb * gb * wbMwb;
      } else {
        fPrime += deltaT * fDoublePrime + gb * gb + theta * gb * zb;
        fDoublePrime += -theta * gb * gb;
      }

      {
        const floor = FDP_EPS * Math.max(1, Math.abs(fDoublePrime));
        if (fDoublePrime < floor) fDoublePrime = floor;
      }

      // Update p ← p + g_b · w_b, then refresh Mp.
      if (col > 0) {
        for (let k = 0; k < twoCol; k++) p[k] += gb * wb[k];
        for (let k = 0; k < twoCol; k++) mp[k] = p[k];
        mat.solveM(mp);
      }

      d[b] = 0;
      dtMin = -fPrime / fDoublePrime;

      if (heap.isEmpty()) break;
      if (fPrime >= 0) {
        dtMin = 0;
        break;
      }

      first = heap.pop();
      tCurrent = first.t;
      b = first.i;
      deltaT = tCurrent - tOld;
    }
  }

  // ---- (4) Finalise -------------------------------------------------
  if (dtMin < 0) dtMin = 0;
  const finalT = tOld + dtMin;

  // All variables with d[i] ≠ 0 are still "active" (not snapped) and
  // advance by the remaining partial step.
  for (let i = 0; i < n; i++) {
    if (d[i] !== 0) {
      let xi = x[i] + finalT * d[i];
      // Absorb eps-scale drift so xc remains feasible.
      const code = nbd[i];
      if (code === BOUND_LOWER || code === BOUND_BOTH) {
        const lo = lower[i];
        if (xi < lo) xi = lo;
      }
      if (code === BOUND_UPPER || code === BOUND_BOTH) {
        const hi = upper[i];
        if (xi > hi) xi = hi;
      }
      xc[i] = xi;
    }
  }

  // c ← c + dtMin · p  (final partial segment; p already reflects all prior snaps).
  for (let k = 0; k < twoCol; k++) c[k] += dtMin * p[k];

  const freeCount = buildFreeSet(xc, lower, upper, freeSet, n);
  return {freeCount, ok: true};
}

function buildFreeSet(
  xc: Float64Array,
  lower: Float64Array,
  upper: Float64Array,
  out: Int32Array,
  n: number,
): number {
  let count = 0;
  for (let i = 0; i < n; i++)
    if (xc[i] > lower[i] && xc[i] < upper[i]) out[count++] = i;

  return count;
}
