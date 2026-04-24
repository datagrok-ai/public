/**
 * BFGSMat compact-representation tests.
 *
 * Covers:
 *   - lower Cholesky correctness against a brute-force reference;
 *   - single-pair BFGSMat agrees with the closed-form BFGS rank-2
 *     update from `B₀ = θI`;
 *   - secant equation `B · s = y` holds for the most recent stored pair;
 *   - `solveM` inverts `K` (block factorisation round-trip);
 *   - `applyW` / `applyWt` agree with direct column-wise formulas;
 *   - ring-buffer eviction at capacity;
 *   - curvature gate rejects bad pairs;
 *   - `reset()` wipes state cleanly.
 */
import {BFGSMat, choleskyLower} from '../optimizers/lbfgs-b/bfgs-mat';

/* ================================================================== */
/*  Helpers                                                            */
/* ================================================================== */

function mulMatMat(A: Float64Array, At: Float64Array, out: Float64Array, n: number): void {
  // out = A · A^T, A is n×n lower tri (column-major).
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      let s = 0;
      for (let k = 0; k < n; k++) s += A[i + k * n] * At[j + k * n];
      out[i + j * n] = s;
    }
  }
}

/** Compute B · v using the compact state directly (forward K application). */
function applyKForward(m: BFGSMat, u: Float64Array, out: Float64Array): void {
  // out[0..col) = −D·u[0..col) + L^T·u[col..2col)
  // out[col..2col) = L·u[0..col) + θ·SᵀS·u[col..2col)
  const col = m.col;
  const M = m.m;
  const theta = m.theta;
  for (let i = 0; i < col; i++) {
    const iP = m.physOf(i);
    // Top block: −D·u1 + L^T·u2
    // (L^T)[i,j] = G_log[j,i] for j > i, else 0.
    const di = m.gPhys[iP + iP * M];
    let top = -di * u[i];
    for (let j = i + 1; j < col; j++) {
      const jP = m.physOf(j);
      top += m.gPhys[jP + iP * M] * u[col + j];
    }
    out[i] = top;
    // Bottom block: L·u1 + θ·SᵀS·u2
    // L[i,j] = G_log[i,j] for j < i, else 0.
    let bot = 0;
    for (let j = 0; j < i; j++) {
      const jP = m.physOf(j);
      bot += m.gPhys[iP + jP * M] * u[j];
    }
    for (let j = 0; j < col; j++) {
      const jP = m.physOf(j);
      bot += theta * m.ssPhys[iP + jP * M] * u[col + j];
    }
    out[col + i] = bot;
  }
}

function explicitBFGSRank2(
  theta: number,
  s: Float64Array,
  y: Float64Array,
  v: Float64Array,
  out: Float64Array,
): void {
  // B₁·v = θv + (yᵀv / sᵀy)·y − (θ sᵀv / sᵀs)·s
  const n = s.length;
  let sy = 0;
  let ss = 0;
  let yv = 0;
  let sv = 0;
  for (let i = 0; i < n; i++) {
    sy += s[i] * y[i];
    ss += s[i] * s[i];
    yv += y[i] * v[i];
    sv += s[i] * v[i];
  }
  const coefY = yv / sy;
  const coefS = (theta * sv) / ss;
  for (let i = 0; i < n; i++)
    out[i] = theta * v[i] + coefY * y[i] - coefS * s[i];
}

/* ================================================================== */
/*  Lower Cholesky                                                     */
/* ================================================================== */

describe('choleskyLower', () => {
  it('factors a 1×1 SPD matrix', () => {
    const A = new Float64Array([4]);
    const L = new Float64Array(1);
    expect(choleskyLower(A, L, 1, 1)).toBe(true);
    expect(L[0]).toBeCloseTo(2, 12);
  });

  it('factors a 3×3 SPD matrix and L·Lᵀ recovers A', () => {
    // A = [[4, 2, 2]; [2, 10, 5]; [2, 5, 11]]
    const A = new Float64Array([
      4, 2, 2,
      2, 10, 5,
      2, 5, 11,
    ]);
    const L = new Float64Array(9);
    expect(choleskyLower(A, L, 3, 3)).toBe(true);
    // Upper triangle zeroed
    expect(L[0 + 1 * 3]).toBe(0);
    expect(L[0 + 2 * 3]).toBe(0);
    expect(L[1 + 2 * 3]).toBe(0);
    // L · Lᵀ ≈ A
    const recon = new Float64Array(9);
    mulMatMat(L, L, recon, 3);
    for (let i = 0; i < 9; i++) expect(recon[i]).toBeCloseTo(A[i], 10);
  });

  it('returns false on non-SPD (negative pivot)', () => {
    const A = new Float64Array([-1]);
    const L = new Float64Array(1);
    expect(choleskyLower(A, L, 1, 1)).toBe(false);
  });

  it('returns false on indefinite matrix', () => {
    // [1 2; 2 1]: eigenvalues 3 and -1 → indefinite.
    const A = new Float64Array([1, 2, 2, 1]);
    const L = new Float64Array(4);
    expect(choleskyLower(A, L, 2, 2)).toBe(false);
  });
});

/* ================================================================== */
/*  Construction / trivial state                                       */
/* ================================================================== */

describe('BFGSMat construction and reset', () => {
  it('throws on n < 1 or m < 1', () => {
    expect(() => new BFGSMat(0, 5)).toThrow('n must be ≥ 1');
    expect(() => new BFGSMat(5, 0)).toThrow('m must be ≥ 1');
  });

  it('applyB with col=0 gives B = θI (identity × θ)', () => {
    const mat = new BFGSMat(3, 5);
    const v = new Float64Array([1, -2, 3]);
    const out = new Float64Array(3);
    const s2m = new Float64Array(10);
    const sN = new Float64Array(3);
    mat.applyB(v, out, s2m, sN);
    expect(out[0]).toBeCloseTo(1, 12);
    expect(out[1]).toBeCloseTo(-2, 12);
    expect(out[2]).toBeCloseTo(3, 12);
  });

  it('reset() restores trivial state', () => {
    const mat = new BFGSMat(3, 5);
    mat.update(
      new Float64Array([1, 0, 0]),
      new Float64Array([2, 0, 0]),
      Number.EPSILON,
    );
    expect(mat.col).toBe(1);
    mat.reset();
    expect(mat.col).toBe(0);
    expect(mat.theta).toBe(1);
    expect(mat.head).toBe(0);
  });
});

/* ================================================================== */
/*  Single-pair compact form vs explicit formula                       */
/* ================================================================== */

describe('BFGSMat single-pair agrees with explicit BFGS rank-2 update', () => {
  it('matches closed form for random s, y, v in n=5', () => {
    const n = 5;
    // s·y = 0.4·1.0 + 0.7·0.5 + 1.1·0.6 + 0.2·0.9 + 0.3·0.2 = 1.77 > 0 ✓
    const s = new Float64Array([0.4, 0.7, 1.1, 0.2, 0.3]);
    const y = new Float64Array([1.0, 0.5, 0.6, 0.9, 0.2]);
    const v = new Float64Array([-0.3, 0.8, 0.1, -0.9, 0.5]);

    const mat = new BFGSMat(n, 5);
    const ok = mat.update(s, y, 0);
    expect(ok).toBe(true);
    expect(mat.col).toBe(1);

    const out = new Float64Array(n);
    const s2m = new Float64Array(10);
    const sN = new Float64Array(n);
    mat.applyB(v, out, s2m, sN);

    const ref = new Float64Array(n);
    explicitBFGSRank2(mat.theta, s, y, v, ref);

    for (let i = 0; i < n; i++) expect(out[i]).toBeCloseTo(ref[i], 10);
  });
});

/* ================================================================== */
/*  Secant equation                                                    */
/* ================================================================== */

describe('BFGSMat secant equation', () => {
  it('B · s = y for the most recent pair (single pair)', () => {
    const n = 4;
    const s = new Float64Array([1, -0.5, 0.3, 0.8]);
    const y = new Float64Array([2, 0.7, -0.4, 1.1]);
    const mat = new BFGSMat(n, 5);
    mat.update(s, y, 0);

    const out = new Float64Array(n);
    const s2m = new Float64Array(10);
    const sN = new Float64Array(n);
    mat.applyB(s, out, s2m, sN);
    for (let i = 0; i < n; i++) expect(out[i]).toBeCloseTo(y[i], 10);
  });

  it('B · sₖ = yₖ for the most recent pair after 3 updates', () => {
    const n = 4;
    const pairs = [
      [new Float64Array([1, 0, 0, 0]), new Float64Array([1.5, 0.2, 0, 0.1])],
      [new Float64Array([0, 1, 0, 0]), new Float64Array([0.3, 2.0, -0.1, 0])],
      [new Float64Array([0.5, 0.5, 0.5, 0.5]), new Float64Array([1.0, 1.5, 0.8, 0.9])],
    ];
    const mat = new BFGSMat(n, 5);
    for (const [s, y] of pairs)
      expect(mat.update(s, y, 0)).toBe(true);

    const last = pairs[pairs.length - 1];
    const out = new Float64Array(n);
    const s2m = new Float64Array(10);
    const sN = new Float64Array(n);
    mat.applyB(last[0], out, s2m, sN);
    for (let i = 0; i < n; i++) expect(out[i]).toBeCloseTo(last[1][i], 9);
  });
});

/* ================================================================== */
/*  solveM inverse check                                               */
/* ================================================================== */

describe('BFGSMat.solveM inverts K', () => {
  it('K · (M · z) ≈ z for random z (m=3)', () => {
    const n = 5;
    const m = 3;
    const mat = new BFGSMat(n, m);
    // Generate 3 (s,y) pairs with s random, y = H·s for some SPD H.
    const Hdiag = [1.5, 2.0, 0.8, 1.2, 1.7];
    const ss = [
      new Float64Array([0.3, -0.5, 0.2, 0.7, 0.1]),
      new Float64Array([1.0, 0.4, -0.6, 0.2, 0.9]),
      new Float64Array([-0.2, 0.7, 0.3, -0.5, 0.8]),
    ];
    for (const s of ss) {
      const y = new Float64Array(n);
      for (let i = 0; i < n; i++) y[i] = Hdiag[i] * s[i];
      expect(mat.update(s, y, 0)).toBe(true);
    }
    const col = mat.col;
    expect(col).toBe(m);

    const z = new Float64Array(2 * col);
    for (let i = 0; i < 2 * col; i++) z[i] = Math.sin(i * 1.3) + 0.1 * i;
    const zOrig = new Float64Array(z);

    mat.solveM(z); // z ← M · z
    const zApplied = new Float64Array(2 * col);
    applyKForward(mat, z, zApplied); // K · (M · z)
    for (let i = 0; i < 2 * col; i++)
      expect(zApplied[i]).toBeCloseTo(zOrig[i], 8);
  });
});

/* ================================================================== */
/*  applyW / applyWt correctness                                       */
/* ================================================================== */

describe('BFGSMat applyW / applyWt', () => {
  it('Wᵀ eⱼ recovers logical column j of W', () => {
    const n = 4;
    const m = 3;
    const mat = new BFGSMat(n, m);
    const ss = [
      new Float64Array([1, 0, 0, 0]),
      new Float64Array([0, 1, 0, 0]),
      new Float64Array([0, 0, 1, 0]),
    ];
    for (const s of ss) {
      const y = new Float64Array(n);
      for (let i = 0; i < n; i++) y[i] = 2 * s[i] + 0.1;
      mat.update(s, y, 0);
    }

    // W = [Y | θS].  W e_0 = Y column 0 (logical).
    const col = mat.col;
    const e = new Float64Array(n);
    const out = new Float64Array(2 * col);
    for (let i = 0; i < n; i++) {
      e[i] = 1;
      mat.applyWt(e, out);
      // out[j] = Y[i, phys(j)], out[col+j] = θ S[i, phys(j)]
      for (let j = 0; j < col; j++) {
        const jP = mat.physOf(j);
        expect(out[j]).toBeCloseTo(mat.yStore[i + jP * n], 12);
        expect(out[col + j]).toBeCloseTo(mat.theta * mat.sStore[i + jP * n], 12);
      }
      e[i] = 0;
    }
  });

  it('applyW ∘ applyWt matches direct W · (Wᵀ v)', () => {
    const n = 5;
    const m = 4;
    const mat = new BFGSMat(n, m);
    for (let k = 0; k < m; k++) {
      const s = new Float64Array(n);
      const y = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        s[i] = Math.sin(k + i);
        y[i] = Math.cos(k + i) + 1.5 * s[i];
      }
      mat.update(s, y, 0);
    }
    const v = new Float64Array([0.5, -0.3, 1.1, 0.2, -0.7]);

    const wtv = new Float64Array(2 * mat.col);
    const wwtv = new Float64Array(n);
    mat.applyWt(v, wtv);
    mat.applyW(wtv, wwtv);

    // Brute force: (W · Wᵀ · v)_i = Σ_j W[i,j] · (Σ_k W[k,j] v[k])
    // where W[:, j] = Y[:, phys(j)] for j < col, θ S[:, phys(j-col)] for j ≥ col.
    const ref = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < mat.col; j++) {
        const jP = mat.physOf(j);
        ref[i] += mat.yStore[i + jP * n] * wtv[j];
        ref[i] += mat.theta * mat.sStore[i + jP * n] * wtv[mat.col + j];
      }
    }
    for (let i = 0; i < n; i++) expect(wwtv[i]).toBeCloseTo(ref[i], 10);
  });
});

/* ================================================================== */
/*  Ring-buffer eviction                                               */
/* ================================================================== */

describe('BFGSMat ring-buffer semantics', () => {
  it('caps col at m and advances head on overflow', () => {
    const n = 3;
    const m = 2;
    const mat = new BFGSMat(n, m);

    for (let k = 0; k < m + 2; k++) {
      const s = new Float64Array([1 + k, 0, 0]);
      const y = new Float64Array([2 + k, 0, 0]);
      mat.update(s, y, 0);
    }
    expect(mat.col).toBe(m);
    expect(mat.head).toBe(2 % m);

    // Newest pair: last one pushed.
    const newestPhys = mat.physOf(mat.col - 1);
    expect(mat.sStore[newestPhys * n]).toBeCloseTo(1 + (m + 1), 12);
    expect(mat.yStore[newestPhys * n]).toBeCloseTo(2 + (m + 1), 12);

    // Secant equation still holds for most recent.
    const out = new Float64Array(n);
    const s2m = new Float64Array(2 * m);
    const sN = new Float64Array(n);
    const lastS = new Float64Array([1 + (m + 1), 0, 0]);
    const lastY = new Float64Array([2 + (m + 1), 0, 0]);
    mat.applyB(lastS, out, s2m, sN);
    for (let i = 0; i < n; i++) expect(out[i]).toBeCloseTo(lastY[i], 10);
  });
});

/* ================================================================== */
/*  Curvature gate                                                     */
/* ================================================================== */

describe('BFGSMat curvature gate', () => {
  it('rejects pair with sᵀy ≤ 0', () => {
    const mat = new BFGSMat(3, 5);
    const s = new Float64Array([1, 0, 0]);
    const y = new Float64Array([-1, 0, 0]); // sᵀy = -1
    const ok = mat.update(s, y, 1e-12);
    expect(ok).toBe(false);
    expect(mat.col).toBe(0);
  });

  it('rejects pair with sᵀy too small relative to yᵀy', () => {
    const mat = new BFGSMat(3, 5);
    const s = new Float64Array([1e-10, 0, 0]);
    const y = new Float64Array([1, 0, 0]);
    // sᵀy = 1e-10, yᵀy = 1. curvatureEps = 1e-8 → threshold = 1e-8.
    // sᵀy = 1e-10 ≤ 1e-8 · 1 ✓ rejected.
    const ok = mat.update(s, y, 1e-8);
    expect(ok).toBe(false);
  });

  it('accepts pair at the threshold +ε', () => {
    const mat = new BFGSMat(3, 5);
    const s = new Float64Array([1, 0, 0]);
    const y = new Float64Array([1, 0, 0]); // sᵀy = 1, yᵀy = 1
    expect(mat.update(s, y, 1e-12)).toBe(true);
  });
});
