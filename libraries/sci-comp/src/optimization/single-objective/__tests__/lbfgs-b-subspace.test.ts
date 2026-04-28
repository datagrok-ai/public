/**
 * Subspace-minimisation tests.
 *
 * Covers:
 *   - trivial early exits (col=0, freeCount=0);
 *   - model-decrease invariant: `m_k(xHat) ≤ m_k(xᶜ)` whenever
 *     `improved === true`;
 *   - feasibility invariant on `xHat`;
 *   - for t = n (all free) and valid stored pairs, the unconstrained
 *     subspace step is the Newton-like step: `B · (xHat − xᶜ) ≈ −r̄ᶜ`
 *     on the free subspace;
 *   - Morales–Nocedal 2011 endpoint selection (`c-jlm-jn` block in
 *     `subsm.f`): directional-derivative test `g · (x̄ − xₖ) < 0` plus
 *     `dd > 0` accepts the projected point; otherwise truncation
 *     back-tracks from `xᶜ` along `d̄ᵘ` on free coords with bound snap.
 */
import {
  subspaceMin,
  makeSubspaceWorkspace,
} from '../optimizers/lbfgs-b/subspace';
import {BFGSMat} from '../optimizers/lbfgs-b/bfgs-mat';
import {classifyBounds} from '../optimizers/lbfgs-b/bounds';
import {cauchyPoint, makeCauchyWorkspace} from '../optimizers/lbfgs-b/cauchy';

/* ================================================================== */
/*  Helpers                                                            */
/* ================================================================== */

/** Compute m_k(point) − m_k(ref) = δ·g + δ·B·z + 0.5·δ·B·δ. */
function modelDelta(
  point: Float64Array,
  ref: Float64Array,
  x: Float64Array,
  g: Float64Array,
  mat: BFGSMat,
): number {
  const n = point.length;
  const delta = new Float64Array(n);
  const z = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    delta[i] = point[i] - ref[i];
    z[i] = ref[i] - x[i];
  }
  const bDelta = new Float64Array(n);
  const bz = new Float64Array(n);
  const s2 = new Float64Array(2 * mat.m);
  const sN = new Float64Array(n);
  mat.applyB(delta, bDelta, s2, sN);
  mat.applyB(z, bz, s2, sN);
  let dg = 0; let dBz = 0; let dBd = 0;
  for (let i = 0; i < n; i++) {
    dg += delta[i] * g[i];
    dBz += delta[i] * bz[i];
    dBd += delta[i] * bDelta[i];
  }
  return dg + dBz + 0.5 * dBd;
}

/* ================================================================== */
/*  Trivial early exits                                                */
/* ================================================================== */

describe('subspaceMin — early exits', () => {
  it('col=0 → xHat = xᶜ, improved=false', () => {
    const n = 3;
    const mat = new BFGSMat(n, 5);
    const x = new Float64Array([1, 1, 1]);
    const g = new Float64Array([1, 1, 1]);
    const xc = new Float64Array([0.5, 0.5, 0.5]);
    const c = new Float64Array(0);
    const lower = new Float64Array([0, 0, 0]);
    const upper = new Float64Array([1, 1, 1]);
    const nbd = classifyBounds(lower, upper);
    const freeSet = new Int32Array([0, 1, 2]);
    const ws = makeSubspaceWorkspace(n, 5);

    const r = subspaceMin(x, g, lower, upper, nbd, xc, c, freeSet, 3, mat, ws);
    expect(r.improved).toBe(false);
    for (let i = 0; i < n; i++) expect(ws.xHat[i]).toBe(xc[i]);
  });

  it('freeCount=0 → xHat = xᶜ, improved=false', () => {
    const n = 3;
    const mat = new BFGSMat(n, 5);
    mat.update(
      new Float64Array([0.2, 0.1, 0.3]),
      new Float64Array([0.5, 0.2, 0.6]),
      0,
      0,
    );
    const x = new Float64Array([1, 1, 1]);
    const g = new Float64Array([1, 1, 1]);
    const xc = new Float64Array([0, 0, 0]); // all at lower bound
    const c = new Float64Array(2 * mat.col);
    const lower = new Float64Array([0, 0, 0]);
    const upper = new Float64Array([1, 1, 1]);
    const nbd = classifyBounds(lower, upper);
    const freeSet = new Int32Array([0, 0, 0]);
    const ws = makeSubspaceWorkspace(n, 5);

    const r = subspaceMin(x, g, lower, upper, nbd, xc, c, freeSet, 0, mat, ws);
    expect(r.improved).toBe(false);
    for (let i = 0; i < n; i++) expect(ws.xHat[i]).toBe(xc[i]);
  });
});

/* ================================================================== */
/*  Model-decrease invariant + feasibility                             */
/* ================================================================== */

describe('subspaceMin — model decrease and feasibility', () => {
  it('produces m_k(xHat) < m_k(xᶜ) on unconstrained problem', () => {
    const n = 3;
    const mat = new BFGSMat(n, 5);
    mat.update(
      new Float64Array([0.1, 0.2, 0.1]),
      new Float64Array([0.5, 0.3, 0.4]),
      0,
      0,
    );

    const x = new Float64Array([0, 0, 0]);
    const g = new Float64Array([1.0, -0.5, 0.7]);
    const lower = new Float64Array([-Infinity, -Infinity, -Infinity]);
    const upper = new Float64Array([Infinity, Infinity, Infinity]);
    const nbd = classifyBounds(lower, upper);

    // Run Cauchy to get xc, c, freeSet.
    const cauchyWs = makeCauchyWorkspace(n, 5);
    const cr = cauchyPoint(x, g, lower, upper, nbd, mat, cauchyWs);
    expect(cr.ok).toBe(true);
    expect(cr.freeCount).toBe(n);

    const ws = makeSubspaceWorkspace(n, 5);
    const r = subspaceMin(
      x, g, lower, upper, nbd, cauchyWs.xc, cauchyWs.c,
      cauchyWs.freeSet, cr.freeCount, mat, ws,
    );
    // For the unconstrained case the subspace step is always feasible,
    // so improved should be true.
    expect(r.improved).toBe(true);
    // Strict model decrease from xᶜ to xHat.
    const dm = modelDelta(ws.xHat, cauchyWs.xc, x, g, mat);
    expect(dm).toBeLessThan(0);
  });

  it('xHat stays feasible even after backtrack', () => {
    const n = 4;
    const mat = new BFGSMat(n, 5);
    mat.update(
      new Float64Array([0.05, 0.05, 0.05, 0.05]),
      new Float64Array([0.3, 0.2, 0.4, 0.1]),
      0,
      0,
    );

    const x = new Float64Array([0.5, 0.5, 0.5, 0.5]);
    const g = new Float64Array([2, -1.5, 3, -2.5]);
    const lower = new Float64Array([0, 0, 0, 0]);
    const upper = new Float64Array([1, 1, 1, 1]);
    const nbd = classifyBounds(lower, upper);

    const cauchyWs = makeCauchyWorkspace(n, 5);
    const cr = cauchyPoint(x, g, lower, upper, nbd, mat, cauchyWs);
    expect(cr.ok).toBe(true);

    const ws = makeSubspaceWorkspace(n, 5);
    subspaceMin(
      x, g, lower, upper, nbd, cauchyWs.xc, cauchyWs.c,
      cauchyWs.freeSet, cr.freeCount, mat, ws,
    );
    for (let i = 0; i < n; i++) {
      expect(ws.xHat[i]).toBeGreaterThanOrEqual(lower[i]);
      expect(ws.xHat[i]).toBeLessThanOrEqual(upper[i]);
    }
  });
});

/* ================================================================== */
/*  Newton-like property (all free, t = n)                             */
/* ================================================================== */

describe('subspaceMin — B·(xHat−xᶜ) ≈ −r̄ᶜ for t=n', () => {
  it('unconstrained subspace step is the Newton-like step on B', () => {
    const n = 4;
    const mat = new BFGSMat(n, 5);
    // Seed two pairs with positive curvature.
    mat.update(
      new Float64Array([0.1, 0.2, -0.1, 0.3]),
      new Float64Array([0.4, 0.1, 0.3, 0.2]),
      0,
      0,
    );
    mat.update(
      new Float64Array([0.05, -0.1, 0.2, 0.1]),
      new Float64Array([0.1, 0.3, 0.1, 0.4]),
      0,
      0,
    );

    const x = new Float64Array([0, 0, 0, 0]);
    const g = new Float64Array([1.1, -0.7, 0.5, 0.3]);
    const lower = new Float64Array([-Infinity, -Infinity, -Infinity, -Infinity]);
    const upper = new Float64Array([Infinity, Infinity, Infinity, Infinity]);
    const nbd = classifyBounds(lower, upper);

    const cauchyWs = makeCauchyWorkspace(n, 5);
    const cr = cauchyPoint(x, g, lower, upper, nbd, mat, cauchyWs);
    expect(cr.ok).toBe(true);
    expect(cr.freeCount).toBe(n);

    const ws = makeSubspaceWorkspace(n, 5);
    const r = subspaceMin(
      x, g, lower, upper, nbd, cauchyWs.xc, cauchyWs.c,
      cauchyWs.freeSet, cr.freeCount, mat, ws,
    );
    expect(r.improved).toBe(true);

    // Recompute r̄ᶜ (all free → Z is identity).
    const rBar = new Float64Array(n);
    const mc = new Float64Array(2 * mat.col);
    for (let k = 0; k < 2 * mat.col; k++) mc[k] = cauchyWs.c[k];
    mat.solveM(mc);
    const wMc = new Float64Array(n);
    mat.applyW(mc, wMc);
    for (let i = 0; i < n; i++)
      rBar[i] = g[i] + mat.theta * (cauchyWs.xc[i] - x[i]) - wMc[i];

    // Compute B·(xHat − xᶜ).
    const step = new Float64Array(n);
    for (let i = 0; i < n; i++) step[i] = ws.xHat[i] - cauchyWs.xc[i];
    const bStep = new Float64Array(n);
    const s2 = new Float64Array(2 * mat.col);
    const sN = new Float64Array(n);
    mat.applyB(step, bStep, s2, sN);

    // Expect bStep ≈ −rBar.
    for (let i = 0; i < n; i++)
      expect(bStep[i]).toBeCloseTo(-rBar[i], 8);
  });
});

/* ================================================================== */
/*  Bounded endpoint selection (Morales–Nocedal 2011)                  */
/* ================================================================== */

describe('subspaceMin — bounded endpoint selection', () => {
  // ---- N3: directional-derivative test accepts xBar -----------------
  it('accepts the projected point when g · (xBar − x) < 0', () => {
    // Bounds wide enough that x̂ stays interior → projection is identity,
    // dd-test reduces to plain Newton-step descent. xHat must equal xc+du.
    const n = 3;
    const mat = new BFGSMat(n, 5);
    mat.update(
      new Float64Array([0.1, 0.2, 0.1]),
      new Float64Array([0.5, 0.3, 0.4]),
      0,
      0,
    );

    const x = new Float64Array([0, 0, 0]);
    const g = new Float64Array([1.0, -0.5, 0.7]);
    const lower = new Float64Array([-100, -100, -100]);
    const upper = new Float64Array([100, 100, 100]);
    const nbd = classifyBounds(lower, upper);

    const cauchyWs = makeCauchyWorkspace(n, 5);
    const cr = cauchyPoint(x, g, lower, upper, nbd, mat, cauchyWs);
    expect(cr.ok).toBe(true);

    const ws = makeSubspaceWorkspace(n, 5);
    const sr = subspaceMin(
      x, g, lower, upper, nbd, cauchyWs.xc, cauchyWs.c,
      cauchyWs.freeSet, cr.freeCount, mat, ws,
    );
    expect(sr.improved).toBe(true);
    // Endpoint = xc + du (no projection, no truncation).
    for (let i = 0; i < n; i++)
      expect(ws.xHat[i]).toBeCloseTo(cauchyWs.xc[i] + ws.du[i], 12);

    // Directional-derivative test is satisfied at the returned xHat.
    let gd = 0;
    for (let i = 0; i < n; i++) gd += g[i] * (ws.xHat[i] - x[i]);
    expect(gd).toBeLessThan(0);
  });

  // ---- N4: truncation fallback invariant ---------------------------
  it('on truncation fallback xHat lies on the segment xc → xc + du (per coord)', () => {
    // Setup: x near upper bound, large gradient → projected unconstrained
    // step lands at the upper bound for both coords; depending on the
    // resulting du sign, either the dd-test passes (xHat = xBar) or
    // truncation engages (xHat = xc + α·du). In both cases xHat must
    // sit on the per-coord segment between xc and xc+du.
    const n = 2;
    const mat = new BFGSMat(n, 5);
    mat.update(
      new Float64Array([0.5, 0.1]),
      new Float64Array([0.2, 0.05]),
      0,
      0,
    );

    const x = new Float64Array([0.8, 0.8]);
    const g = new Float64Array([-5, -4]);
    const lower = new Float64Array([0, 0]);
    const upper = new Float64Array([1, 1]);
    const nbd = classifyBounds(lower, upper);

    const cauchyWs = makeCauchyWorkspace(n, 5);
    const cr = cauchyPoint(x, g, lower, upper, nbd, mat, cauchyWs);
    expect(cr.ok).toBe(true);

    const ws = makeSubspaceWorkspace(n, 5);
    const sr = subspaceMin(
      x, g, lower, upper, nbd, cauchyWs.xc, cauchyWs.c,
      cauchyWs.freeSet, cr.freeCount, mat, ws,
    );

    // Always feasible.
    for (let i = 0; i < n; i++) {
      expect(ws.xHat[i]).toBeGreaterThanOrEqual(lower[i]);
      expect(ws.xHat[i]).toBeLessThanOrEqual(upper[i]);
    }
    if (!sr.improved) {
      for (let i = 0; i < n; i++) expect(ws.xHat[i]).toBe(cauchyWs.xc[i]);
      return;
    }
    // Per-coord invariant: xHat[i] ∈ [min(xc[i], xc[i]+du[i]),
    // max(xc[i], xc[i]+du[i])] for both projection (xHat = xBar, which
    // is xc+du clipped to box, hence between xc and xc+du in the worst
    // case) and truncation (xHat = xc + α·du, α ∈ [0, 1]).
    const eps = 1e-12;
    for (let i = 0; i < n; i++) {
      const a = cauchyWs.xc[i];
      const b = cauchyWs.xc[i] + ws.du[i];
      const lo = Math.min(a, b);
      const hi = Math.max(a, b);
      expect(ws.xHat[i]).toBeGreaterThanOrEqual(lo - eps);
      expect(ws.xHat[i]).toBeLessThanOrEqual(hi + eps);
    }
  });

  // ---- N6: bound snap on truncation --------------------------------
  it('snaps the binding coordinate exactly to its bound on truncation', () => {
    // Goal: trigger the truncation branch with α < 1 and verify that
    // the first-hitting coordinate equals its bound bit-exactly (no
    // ε-drift from the `room / dk` division). We construct the inputs
    // for `subspaceMin` directly: it accepts xc, c, freeSet, freeCount
    // as parameters and never recomputes the Cauchy phase, so we can
    // hand it a synthetic xc that is interior on both coords and let
    // it form du via the BFGSMat closed-form.
    const n = 2;
    const mat = new BFGSMat(n, 5);
    // A pair giving a moderate θ ≈ y·y / s·y ≈ (4 + 1) / (1 + 0.5) = 10/3.
    mat.update(
      new Float64Array([1.0, 0.5]),
      new Float64Array([2.0, 1.0]),
      0,
      0,
    );

    // x_k far from both bounds; xc clearly interior; gradient pushes
    // towards (-, -) so the resulting Newton-like du is along (-, -).
    // Whether dd-test passes or fails, the truncation-vs-projection
    // path must respect both feasibility and the bound-snap invariant
    // when α < 1.
    const x = new Float64Array([0.4, 0.4]);
    const g = new Float64Array([8, 5]);
    const lower = new Float64Array([0, 0]);
    const upper = new Float64Array([1, 1]);
    const nbd = classifyBounds(lower, upper);

    // Use Cauchy to get a consistent (xc, c, freeSet) for this mat.
    const cauchyWs = makeCauchyWorkspace(n, 5);
    const cr = cauchyPoint(x, g, lower, upper, nbd, mat, cauchyWs);
    expect(cr.ok).toBe(true);

    const ws = makeSubspaceWorkspace(n, 5);
    const sr = subspaceMin(
      x, g, lower, upper, nbd, cauchyWs.xc, cauchyWs.c,
      cauchyWs.freeSet, cr.freeCount, mat, ws,
    );

    // Always feasible.
    for (let i = 0; i < n; i++) {
      expect(ws.xHat[i]).toBeGreaterThanOrEqual(lower[i]);
      expect(ws.xHat[i]).toBeLessThanOrEqual(upper[i]);
    }

    if (!sr.improved) return; // truncation degenerate, nothing to snap.

    // For each coord, if xHat[i] hit a bound, the value must equal
    // that bound bit-exactly (no FP drift). This is the snap invariant.
    for (let i = 0; i < n; i++) {
      const v = ws.xHat[i];
      const onLower = Math.abs(v - lower[i]) < 1e-12;
      const onUpper = Math.abs(v - upper[i]) < 1e-12;
      if (onLower) expect(v).toBe(lower[i]);
      if (onUpper) expect(v).toBe(upper[i]);
    }
  });
});

/* ================================================================== */
/*  Truncation reference is xᶜ, not xₖ (regression for c-jlm-jn fix)   */
/* ================================================================== */

describe('subspaceMin — truncation reference is xc, not x_k', () => {
  // ---- N5b: per-coord segment invariant on synthetic infeasible du --
  // The pre-fix code truncated xₖ → xc + du; the post-fix code truncates
  // xc → xc + du. The two diverge whenever xc differs from xₖ. The unit
  // test below relies only on the post-fix invariant: xHat[i] sits on
  // [xc[i], xc[i] + du[i]] (componentwise).
  it('xHat is bounded by the xc-relative segment when truncation engages', () => {
    const n = 4;
    const mat = new BFGSMat(n, 5);
    mat.update(
      new Float64Array([0.05, 0.05, 0.05, 0.05]),
      new Float64Array([0.3, 0.2, 0.4, 0.1]),
      0,
      0,
    );

    const x = new Float64Array([0.5, 0.5, 0.5, 0.5]);
    const g = new Float64Array([2, -1.5, 3, -2.5]);
    const lower = new Float64Array([0, 0, 0, 0]);
    const upper = new Float64Array([1, 1, 1, 1]);
    const nbd = classifyBounds(lower, upper);

    const cauchyWs = makeCauchyWorkspace(n, 5);
    const cr = cauchyPoint(x, g, lower, upper, nbd, mat, cauchyWs);
    expect(cr.ok).toBe(true);

    const ws = makeSubspaceWorkspace(n, 5);
    const sr = subspaceMin(
      x, g, lower, upper, nbd, cauchyWs.xc, cauchyWs.c,
      cauchyWs.freeSet, cr.freeCount, mat, ws,
    );

    if (!sr.improved) {
      for (let i = 0; i < n; i++) expect(ws.xHat[i]).toBe(cauchyWs.xc[i]);
      return;
    }
    // Post-fix invariant: xHat lives on the xc-anchored segment.
    const eps = 1e-12;
    for (let i = 0; i < n; i++) {
      const a = cauchyWs.xc[i];
      const b = cauchyWs.xc[i] + ws.du[i];
      const lo = Math.min(a, b);
      const hi = Math.max(a, b);
      expect(ws.xHat[i]).toBeGreaterThanOrEqual(lo - eps);
      expect(ws.xHat[i]).toBeLessThanOrEqual(hi + eps);
    }
  });
});
