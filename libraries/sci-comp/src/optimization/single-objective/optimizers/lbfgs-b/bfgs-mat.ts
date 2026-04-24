/**
 * Compact limited-memory BFGS representation.
 *
 * Stores the last `m` curvature pairs `(sᵢ, yᵢ)` in ring-buffered
 * column-major layouts and applies the BFGS matrix
 *
 *    Bₖ = θₖ I − Wₖ Mₖ Wₖᵀ,      Wₖ = [Yₖ | θₖ Sₖ] ∈ ℝⁿˣ²ᵐ,
 *
 * in O(mn + m²) flops per product. The inverse middle matrix Mₖ = Kₖ⁻¹
 * is applied via the block factorisation
 *
 *    Kₖ = [ −D        Lᵀ  ]  =  [  D^½         0 ] · [ −D^½   D^{−½}Lᵀ ]
 *         [  L     θ SᵀS ]     [ −L D^{−½}    J ]   [  0        Jᵀ    ]
 *
 * with `L` the strict lower triangle of `G = SᵀY`, `D = diag(G)`, and
 * `J` the lower Cholesky factor of `T = θ SᵀS + L D⁻¹ Lᵀ`.
 *
 * Only the physical (ring-buffer) indexing of `G` and `SᵀS` is
 * persisted between iterations; the logical-order Cholesky `J` is
 * rebuilt after every accepted update. For `m ≤ 20` the rebuild is
 * O(m³) ≈ a few thousand flops and is negligible next to the O(mn)
 * curvature-pair bookkeeping.
 *
 * References:
 *  - R. H. Byrd, J. Nocedal, R. B. Schnabel (1994), *Representations
 *    of Quasi-Newton Matrices and their Use in Limited Memory Methods*,
 *    Math. Programming 63:129–156.
 *  - R. H. Byrd, P. Lu, J. Nocedal, C. Zhu (1995) §3.
 */

/* ================================================================== */
/*  BFGSMat                                                            */
/* ================================================================== */

/**
 * Compact L-BFGS memory with the operations needed by the Cauchy-point
 * routine, subspace minimisation, and the outer driver.
 */
export class BFGSMat {
  readonly n: number;
  readonly m: number;

  /** `s` columns in column-major `n × m`, physical indexing. */
  readonly sStore: Float64Array;
  /** `y` columns in column-major `n × m`, physical indexing. */
  readonly yStore: Float64Array;
  /** `gPhys[iPhys + jPhys·m] = sᵢ · yⱼ` (non-symmetric). */
  readonly gPhys: Float64Array;
  /** `ssPhys[iPhys + jPhys·m] = sᵢ · sⱼ` (symmetric). */
  readonly ssPhys: Float64Array;
  /** Lower Cholesky of `T` in column-major, logical indexing. */
  readonly lChol: Float64Array;

  /** Shanno–Phua scalar `θ = yᵀy / sᵀy` of the newest pair. */
  theta: number;
  /** Number of pairs currently stored (0..m). */
  col: number;
  /** Physical index of the oldest pair (ring-buffer head). */
  head: number;
  /** `true` once Cholesky is up-to-date for the current state. */
  private choleskyValid: boolean;

  // Scratch
  private readonly tMat: Float64Array; // m×m, logical-indexed, symmetric
  private readonly dSqrt: Float64Array; // m: √(G_log[i,i])

  constructor(n: number, m: number) {
    if (!Number.isInteger(n) || n < 1) throw new Error('BFGSMat: n must be ≥ 1');
    if (!Number.isInteger(m) || m < 1) throw new Error('BFGSMat: m must be ≥ 1');
    this.n = n;
    this.m = m;
    this.sStore = new Float64Array(n * m);
    this.yStore = new Float64Array(n * m);
    this.gPhys = new Float64Array(m * m);
    this.ssPhys = new Float64Array(m * m);
    this.lChol = new Float64Array(m * m);
    this.tMat = new Float64Array(m * m);
    this.dSqrt = new Float64Array(m);
    this.theta = 1;
    this.col = 0;
    this.head = 0;
    this.choleskyValid = true; // vacuously (col=0)
  }

  /** Wipe memory: `col = 0`, `theta = 1`. Used on line-search reset. */
  reset(): void {
    this.theta = 1;
    this.col = 0;
    this.head = 0;
    this.choleskyValid = true;
  }

  /** Physical-slot index for logical position `k ∈ [0, col)`. */
  physOf(k: number): number {
    return (this.head + k) % this.m;
  }

  /* ----------------------------------------------------------------- */
  /*  Pair update                                                      */
  /* ----------------------------------------------------------------- */

  /**
   * Try to push a new curvature pair `(sNew, yNew)`.
   * Curvature gate: accept iff `sᵀy > curvatureEps · max(1, yᵀy)`.
   *
   * On accept: updates `θ`, `G`, `SᵀS`, rebuilds Cholesky. Returns true.
   * On reject: memory unchanged. Returns false.
   */
  update(sNew: Float64Array, yNew: Float64Array, curvatureEps: number): boolean {
    const n = this.n;
    const m = this.m;

    if (sNew.length !== n || yNew.length !== n)
      throw new Error('BFGSMat.update: sNew/yNew length mismatch');

    let ys = 0;
    let yy = 0;
    let ssNew = 0;
    for (let i = 0; i < n; i++) {
      const si = sNew[i];
      const yi = yNew[i];
      ys += si * yi;
      yy += yi * yi;
      ssNew += si * si;
    }
    if (!(yy > 0) || !(ys > curvatureEps * Math.max(1, yy)))
      return false;

    const newP = this.col < m ? (this.head + this.col) % m : this.head;
    if (this.col < m) this.col++;
    else this.head = (this.head + 1) % m;

    const off = newP * n;
    for (let i = 0; i < n; i++) {
      this.sStore[off + i] = sNew[i];
      this.yStore[off + i] = yNew[i];
    }

    // Fill G[*, newP], G[newP, *], SS[*, newP] (symmetric).
    // This iterates over all current pairs including the new one.
    for (let kl = 0; kl < this.col; kl++) {
      const kP = this.physOf(kl);
      if (kP === newP) continue;
      const kOff = kP * n;
      let sNewYk = 0; // sNew · y_k  → gPhys[newP, kP]
      let skYNew = 0; // s_k · yNew  → gPhys[kP, newP]
      let sNewSk = 0; // sNew · s_k  → ssPhys[newP, kP] = ssPhys[kP, newP]
      for (let i = 0; i < n; i++) {
        sNewYk += sNew[i] * this.yStore[kOff + i];
        skYNew += this.sStore[kOff + i] * yNew[i];
        sNewSk += sNew[i] * this.sStore[kOff + i];
      }
      this.gPhys[newP + kP * m] = sNewYk;
      this.gPhys[kP + newP * m] = skYNew;
      this.ssPhys[newP + kP * m] = sNewSk;
      this.ssPhys[kP + newP * m] = sNewSk;
    }
    this.gPhys[newP + newP * m] = ys;
    this.ssPhys[newP + newP * m] = ssNew;

    this.theta = yy / ys;

    this.rebuildCholesky();
    return true;
  }

  /* ----------------------------------------------------------------- */
  /*  Cholesky of T = θ SᵀS + L D⁻¹ Lᵀ                                  */
  /* ----------------------------------------------------------------- */

  private rebuildCholesky(): void {
    const m = this.m;
    const col = this.col;
    const theta = this.theta;

    for (let i = 0; i < col; i++) {
      const iP = this.physOf(i);
      const di = this.gPhys[iP + iP * m];
      // Curvature gate ensured d > 0; guard anyway.
      this.dSqrt[i] = di > 0 ? Math.sqrt(di) : 0;
    }

    // Build T (logical, column-major, upper triangle used; symmetric).
    for (let j = 0; j < col; j++) {
      const jP = this.physOf(j);
      for (let i = 0; i <= j; i++) {
        const iP = this.physOf(i);
        let sum = theta * this.ssPhys[iP + jP * m];
        const kMax = i < j ? i : j;
        for (let k = 0; k < kMax; k++) {
          const kP = this.physOf(k);
          const dk = this.gPhys[kP + kP * m];
          if (dk > 0)
            sum += this.gPhys[iP + kP * m] * this.gPhys[jP + kP * m] / dk;
        }
        this.tMat[i + j * m] = sum;
        if (i !== j) this.tMat[j + i * m] = sum;
      }
    }

    this.choleskyValid = choleskyLower(this.tMat, this.lChol, col, m);
    if (!this.choleskyValid)
      throw new Error('BFGSMat: Cholesky of T failed — memory likely inconsistent');
  }

  /* ----------------------------------------------------------------- */
  /*  Compact-form products                                            */
  /* ----------------------------------------------------------------- */

  /**
   * `result ← Wᵀ v` where `result` has length `2 · col`. Logical order:
   * `result[0..col) = Yᵀv`, `result[col..2col) = θ Sᵀv`.
   */
  applyWt(v: Float64Array, result: Float64Array): void {
    const n = this.n;
    const col = this.col;
    const theta = this.theta;
    for (let j = 0; j < col; j++) {
      const jP = this.physOf(j);
      const off = jP * n;
      let yDot = 0;
      let sDot = 0;
      for (let i = 0; i < n; i++) {
        yDot += this.yStore[off + i] * v[i];
        sDot += this.sStore[off + i] * v[i];
      }
      result[j] = yDot;
      result[col + j] = theta * sDot;
    }
  }

  /**
   * `result ← W u` where `u` has length `2 · col`:
   * `result = Y · u[0..col) + θ S · u[col..2col)`.
   */
  applyW(u: Float64Array, result: Float64Array): void {
    const n = this.n;
    const col = this.col;
    const theta = this.theta;
    for (let i = 0; i < n; i++) result[i] = 0;
    for (let j = 0; j < col; j++) {
      const jP = this.physOf(j);
      const off = jP * n;
      const uTop = u[j];
      const uBot = theta * u[col + j];
      for (let i = 0; i < n; i++)
        result[i] += this.yStore[off + i] * uTop + this.sStore[off + i] * uBot;
    }
  }

  /**
   * In-place: `z ← Mₖ z` where `z` has length `2 · col`.
   * Uses the block factorisation `K = A · B` described in the file header.
   *
   * Layout: `z[0..col) = z₁`, `z[col..2col) = z₂`.
   */
  solveM(z: Float64Array): void {
    const col = this.col;
    if (col === 0) return;
    if (!this.choleskyValid)
      throw new Error('BFGSMat.solveM: Cholesky not valid');

    const m = this.m;
    const dSqrt = this.dSqrt;
    const lChol = this.lChol;

    // ---- Solve A y = z ---------------------------------------------
    // y1[i] = z1[i] / dSqrt[i]
    for (let i = 0; i < col; i++) z[i] /= dSqrt[i];

    // rhs2[i] = z2[i] + sum_{j<i} G_log[i,j] * y1[j] / dSqrt[j]
    for (let i = 0; i < col; i++) {
      const iP = this.physOf(i);
      let s = z[col + i];
      for (let j = 0; j < i; j++) {
        const jP = this.physOf(j);
        s += this.gPhys[iP + jP * m] * z[j] / dSqrt[j];
      }
      z[col + i] = s;
    }

    // Forward solve L_chol y2 = rhs2.
    for (let i = 0; i < col; i++) {
      let s = z[col + i];
      for (let j = 0; j < i; j++) s -= lChol[i + j * m] * z[col + j];
      z[col + i] = s / lChol[i + i * m];
    }

    // ---- Solve B x = y ---------------------------------------------
    // Back solve L_chol^T x2 = y2.
    for (let i = col - 1; i >= 0; i--) {
      let s = z[col + i];
      for (let j = i + 1; j < col; j++) s -= lChol[j + i * m] * z[col + j];
      z[col + i] = s / lChol[i + i * m];
    }

    // x1[i] = (sum_{j>i} G_log[j,i] * x2[j]) / G_log[i,i] − y1[i] / dSqrt[i]
    for (let i = 0; i < col; i++) {
      const iP = this.physOf(i);
      let ltx = 0;
      for (let j = i + 1; j < col; j++) {
        const jP = this.physOf(j);
        ltx += this.gPhys[jP + iP * m] * z[col + j];
      }
      const di = this.gPhys[iP + iP * m];
      z[i] = ltx / di - z[i] / dSqrt[i];
    }
  }

  /**
   * `result ← Bₖ v`.  Uses `scratch2m` (length `≥ 2m`) and `scratchN`
   * (length `≥ n`) as workspace; caller owns both.
   */
  applyB(v: Float64Array, result: Float64Array, scratch2m: Float64Array, scratchN: Float64Array): void {
    const n = this.n;
    const theta = this.theta;
    for (let i = 0; i < n; i++) result[i] = theta * v[i];
    if (this.col === 0) return;
    this.applyWt(v, scratch2m);
    this.solveM(scratch2m);
    this.applyW(scratch2m, scratchN);
    for (let i = 0; i < n; i++) result[i] -= scratchN[i];
  }
}

/* ================================================================== */
/*  Private linalg helpers                                             */
/* ================================================================== */

/**
 * In-place lower Cholesky: given SPD `A` (`n × n`, column-major, stride
 * `stride`), writes `L` such that `A = L Lᵀ` into `out`. Upper triangle
 * of `out` is zeroed. `A` may alias `out`. Returns `false` when `A` is
 * not positive definite (some diagonal pivot ≤ 0).
 *
 * Exported for unit testing; not part of the public `BFGSMat` surface.
 */
export function choleskyLower(
  A: Float64Array,
  out: Float64Array,
  n: number,
  stride: number,
): boolean {
  for (let j = 0; j < n; j++) {
    let diag = A[j + j * stride];
    for (let k = 0; k < j; k++) {
      const ljk = out[j + k * stride];
      diag -= ljk * ljk;
    }
    if (!(diag > 0)) return false;
    const ljj = Math.sqrt(diag);
    out[j + j * stride] = ljj;
    for (let i = j + 1; i < n; i++) {
      let sum = A[i + j * stride];
      for (let k = 0; k < j; k++) sum -= out[i + k * stride] * out[j + k * stride];
      out[i + j * stride] = sum / ljj;
    }
    for (let i = 0; i < j; i++) out[i + j * stride] = 0;
  }
  return true;
}
