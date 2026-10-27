/* Correlation-aware MPO weights (caMPO): Garellick et al., ChemRxiv 2026, doi:10.26434/chemrxiv.15008121/v1.
   Pure functions over typed arrays, no platform dependencies. Matrices are row-major Float64Arrays. */

/** Desirability columns reduced to the rows where every property has a value, without constant columns. */
export interface PreparedDesirabilities {
  /** Complete rows of the columns in `active`, in the same order. */
  cols: Float32Array[];
  /** Indices (in property order) of the columns that enter the correlation matrix. */
  active: Int32Array;
  rows: number;
}

/** `D[k]` is null for a property with no column; NaN in `D[k][i]` means row i does not contribute. */
export function prepareDesirabilities(D: (Float32Array | null)[], n: number): PreparedDesirabilities {
  const d = D.length;
  const mask = new Uint8Array(n).fill(1);
  for (const a of D) {
    if (!a)
      continue;
    for (let i = 0; i < n; i++) {
      if (Number.isNaN(a[i]))
        mask[i] = 0;
    }
  }
  let rows = 0;
  for (let i = 0; i < n; i++)
    rows += mask[i];
  if (rows <= d)
    return {cols: [], active: new Int32Array(0), rows};

  const active: number[] = [];
  const cols: Float32Array[] = [];
  for (let k = 0; k < d; k++) {
    const a = D[k];
    if (!a)
      continue;
    const out = new Float32Array(rows);
    let sum = 0;
    for (let i = 0, r = 0; i < n; i++) {
      if (mask[i]) {
        out[r++] = a[i];
        sum += a[i];
      }
    }
    const mean = sum / rows;
    let v = 0;
    for (let r = 0; r < rows; r++)
      v += (out[r] - mean) * (out[r] - mean);
    if (v > 1e-12) {
      active.push(k);
      cols.push(out);
    }
  }
  return {cols, active: Int32Array.from(active), rows};
}

/** Pearson correlation matrix of equal-length columns with no missing values and non-zero variance. */
export function pearsonMatrix(cols: Float32Array[], rows: number): Float64Array {
  const m = cols.length;
  const centered = cols.map((a) => {
    let sum = 0;
    for (let r = 0; r < rows; r++)
      sum += a[r];
    const mean = sum / rows;
    const c = new Float64Array(rows);
    for (let r = 0; r < rows; r++)
      c[r] = a[r] - mean;
    return c;
  });
  const norm = centered.map((c) => {
    let s = 0;
    for (let r = 0; r < rows; r++)
      s += c[r] * c[r];
    return Math.sqrt(s);
  });

  const R = new Float64Array(m * m);
  for (let p = 0; p < m; p++) {
    R[p * m + p] = 1;
    for (let q = p + 1; q < m; q++) {
      const a = centered[p];
      const b = centered[q];
      let s = 0;
      for (let r = 0; r < rows; r++)
        s += a[r] * b[r];
      R[p * m + q] = R[q * m + p] = Math.max(-1, Math.min(1, s / (norm[p] * norm[q])));
    }
  }
  return R;
}

/** Lower-triangular L with A = L·Lᵀ, or null if A is not positive definite. */
export function cholesky(A: Float64Array, m: number): Float64Array | null {
  const L = new Float64Array(m * m);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j <= i; j++) {
      let s = A[i * m + j];
      for (let k = 0; k < j; k++)
        s -= L[i * m + k] * L[j * m + k];
      if (i === j) {
        if (s <= 1e-12)
          return null;
        L[i * m + i] = Math.sqrt(s);
      } else
        L[i * m + j] = s / L[j * m + j];
    }
  }
  return L;
}

/** Solves L·Lᵀ·x = b. */
export function choleskySolve(L: Float64Array, m: number, b: Float64Array): Float64Array {
  const y = new Float64Array(m);
  for (let i = 0; i < m; i++) {
    let s = b[i];
    for (let k = 0; k < i; k++)
      s -= L[i * m + k] * y[k];
    y[i] = s / L[i * m + i];
  }
  const x = new Float64Array(m);
  for (let i = m - 1; i >= 0; i--) {
    let s = y[i];
    for (let k = i + 1; k < m; k++)
      s -= L[k * m + i] * x[k];
    x[i] = s / L[i * m + i];
  }
  return x;
}

export function choleskyInverse(L: Float64Array, m: number): Float64Array {
  const inv = new Float64Array(m * m);
  const e = new Float64Array(m);
  for (let c = 0; c < m; c++) {
    e.fill(0);
    e[c] = 1;
    const col = choleskySolve(L, m, e);
    for (let r = 0; r < m; r++)
      inv[r * m + c] = col[r];
  }
  return inv;
}

/** Shrinkage towards the identity: (1 − λ)·R + λ·I. */
export function regularise(R: Float64Array, m: number, lambda: number): Float64Array {
  const out = new Float64Array(m * m);
  for (let i = 0; i < m * m; i++)
    out[i] = (1 - lambda) * R[i];
  for (let k = 0; k < m; k++)
    out[k * m + k] += lambda;
  return out;
}

function norm1(A: Float64Array, m: number): number {
  let max = 0;
  for (let j = 0; j < m; j++) {
    let s = 0;
    for (let i = 0; i < m; i++)
      s += Math.abs(A[i * m + j]);
    max = Math.max(max, s);
  }
  return max;
}

export const MAX_CONDITION = 1e6;
export const LAMBDAS = [0, 1e-3, 1e-2, 5e-2, 0.1];

export interface InverseResult {
  Rinv: Float64Array;
  /** Shrinkage that was needed; 0 when R was inverted as it is. */
  lambda: number;
}

/** Inverse of a correlation matrix, regularised with the smallest λ that makes it well-conditioned. */
export function robustInverse(R: Float64Array, m: number): InverseResult | null {
  for (const lambda of LAMBDAS) {
    const A = lambda === 0 ? R : regularise(R, m, lambda);
    const L = cholesky(A, m);
    if (!L)
      continue;
    const Rinv = choleskyInverse(L, m);
    if (norm1(A, m) * norm1(Rinv, m) <= MAX_CONDITION)
      return {Rinv, lambda};
  }
  return null;
}

export interface BaseWeights {
  /** s_i / S over the properties in R; may be negative; NaN for properties outside R. */
  raw: Float64Array;
  /** Negative raw weights set to 0 and the rest renormalised to sum to 1; 0 outside R. */
  base: Float64Array;
  /** True if at least one raw weight was negative. */
  clipped: boolean;
}

/** w = R⁻¹·1 / (1ᵀ·R⁻¹·1): row sums of the inverse over the sum of all its elements. */
export function campoWeights(Rinv: Float64Array, m: number, active: Int32Array, d: number): BaseWeights {
  const s = new Float64Array(m);
  let S = 0;
  for (let i = 0; i < m; i++) {
    for (let j = 0; j < m; j++)
      s[i] += Rinv[i * m + j];
    S += s[i];
  }

  const raw = new Float64Array(d).fill(NaN);
  const base = new Float64Array(d);
  let clipped = false;
  let pos = 0;
  for (let i = 0; i < m; i++) {
    const r = s[i] / S;
    raw[active[i]] = r;
    if (r < 0)
      clipped = true;
    else
      pos += r;
  }
  for (let i = 0; i < m; i++) {
    const r = raw[active[i]];
    base[active[i]] = r > 0 ? r / pos : 0;
  }
  return {raw, base, clipped};
}

export type FinalWeights = {final: Float64Array, error: null} | {final: null, error: string};

/** Applies priorities, then pins: pinned weights are fixed, the rest scaled to fill 1 − Σ pinned.
 *  `pinned[i]` is NaN for a property that is not pinned. */
export function finalWeights(base: Float64Array, prio: Float64Array, pinned: Float64Array): FinalWeights {
  const d = base.length;
  const b = new Float64Array(d);
  let B = 0;
  let P = 0;
  for (let i = 0; i < d; i++) {
    if (Number.isNaN(pinned[i])) {
      b[i] = base[i] * prio[i];
      B += b[i];
    } else
      P += pinned[i];
  }
  if (P > 1 + 1e-9)
    return {final: null, error: `Pinned weights sum to ${P.toFixed(3)} (> 1)`};

  const rest = 1 - P;
  if (rest > 1e-9 && B <= 1e-12) {
    return {final: null, error: P > 0 ? 'Nothing to scale: all unpinned weights are 0' :
      'All weighted priorities are 0'};
  }

  const k = B > 0 ? rest / B : 0;
  const final = new Float64Array(d);
  for (let i = 0; i < d; i++)
    final[i] = Number.isNaN(pinned[i]) ? b[i] * k : pinned[i];
  return {final, error: null};
}

/** Share of property `pos`'s variance explained by the others: R²ᵢ = 1 − 1 / R⁻¹ᵢᵢ. */
export function explainedShare(Rinv: Float64Array, m: number, pos: number): number {
  return 1 - 1 / Rinv[pos * m + pos];
}
