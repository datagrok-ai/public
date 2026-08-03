/**
 * Cauchy-point tests.
 *
 * Covers:
 *   - MinHeap basics (empty, single, sorted extraction);
 *   - Unconstrained case (all Free, col=0 / col=1): `xᶜ = x + (1/f'') d`;
 *   - 1-D hitting a bound (closed-form comparison);
 *   - Fixed variables (l = u) stay at the bound;
 *   - Mixed 2-D with known symbolic result;
 *   - Free-set correctness (strict interior only);
 *   - Post-condition `c = Wᵀ(xᶜ − x)` for col ≥ 1;
 *   - Feasibility invariant: `l ≤ xᶜ ≤ u` always.
 */
import {
  MinHeap,
  cauchyPoint,
  makeCauchyWorkspace,
} from '../optimizers/lbfgs-b/cauchy';
import {BFGSMat} from '../optimizers/lbfgs-b/bfgs-mat';
import {classifyBounds} from '../optimizers/lbfgs-b/bounds';

/* ================================================================== */
/*  MinHeap                                                            */
/* ================================================================== */

describe('MinHeap', () => {
  it('starts empty', () => {
    const h = new MinHeap(4);
    expect(h.isEmpty()).toBe(true);
    expect(h.size).toBe(0);
  });

  it('push and pop in sorted order', () => {
    const h = new MinHeap(8);
    const vals = [5, 2, 7, 1, 6, 3, 4];
    for (let i = 0; i < vals.length; i++) h.push(vals[i], i);
    expect(h.size).toBe(vals.length);
    const popped: number[] = [];
    while (!h.isEmpty()) popped.push(h.pop().t);
    expect(popped).toEqual([1, 2, 3, 4, 5, 6, 7]);
  });

  it('preserves index with key', () => {
    const h = new MinHeap(4);
    h.push(10, 100);
    h.push(5, 50);
    expect(h.pop()).toEqual({t: 5, i: 50});
    expect(h.pop()).toEqual({t: 10, i: 100});
  });

  it('clear() resets', () => {
    const h = new MinHeap(4);
    h.push(3, 0); h.push(1, 1);
    h.clear();
    expect(h.isEmpty()).toBe(true);
  });
});

/* ================================================================== */
/*  Cauchy — unconstrained                                             */
/* ================================================================== */

describe('cauchyPoint — unconstrained (all Free)', () => {
  it('col=0: xᶜ = x − g/θ', () => {
    const n = 3;
    const mat = new BFGSMat(n, 5);
    const x = new Float64Array([1, -2, 3]);
    const g = new Float64Array([0.5, -1, 2]);
    const lower = new Float64Array([-Infinity, -Infinity, -Infinity]);
    const upper = new Float64Array([Infinity, Infinity, Infinity]);
    const nbd = classifyBounds(lower, upper);
    const ws = makeCauchyWorkspace(n, 5);

    const r = cauchyPoint(x, g, lower, upper, nbd, mat, ws);
    expect(r.ok).toBe(true);
    expect(r.freeCount).toBe(n);

    // With col=0, θ=1: f' = −‖g‖², f'' = θ‖g‖² = ‖g‖², dtMin = 1/θ = 1.
    for (let i = 0; i < n; i++)
      expect(ws.xc[i]).toBeCloseTo(x[i] - g[i] / mat.theta, 12);
  });

  it('col=1: xᶜ = x + dtMin · (−g), c = Wᵀ(xᶜ − x)', () => {
    const n = 3;
    const mat = new BFGSMat(n, 5);
    // Push a pair with positive curvature.
    const s0 = new Float64Array([0.3, 0.4, 0.1]);
    const y0 = new Float64Array([0.6, 0.5, 0.2]);
    mat.update(s0, y0, 0, 0);

    const x = new Float64Array([0, 0, 0]);
    const g = new Float64Array([0.5, -0.3, 0.2]);
    const lower = new Float64Array([-Infinity, -Infinity, -Infinity]);
    const upper = new Float64Array([Infinity, Infinity, Infinity]);
    const nbd = classifyBounds(lower, upper);
    const ws = makeCauchyWorkspace(n, 5);

    const r = cauchyPoint(x, g, lower, upper, nbd, mat, ws);
    expect(r.ok).toBe(true);

    // xᶜ − x proportional to −g (direction preserved in unconstrained case).
    const alpha = -(ws.xc[0] - x[0]) / g[0];
    for (let i = 1; i < n; i++)
      expect(-(ws.xc[i] - x[i]) / g[i]).toBeCloseTo(alpha, 10);

    expect(alpha).toBeGreaterThan(0);

    // Verify c = Wᵀ(xᶜ − x).
    const delta = new Float64Array(n);
    for (let i = 0; i < n; i++) delta[i] = ws.xc[i] - x[i];
    const cRef = new Float64Array(2 * mat.col);
    mat.applyWt(delta, cRef);
    for (let k = 0; k < 2 * mat.col; k++)
      expect(ws.c[k]).toBeCloseTo(cRef[k], 10);
  });
});

/* ================================================================== */
/*  Cauchy — 1-D bounded                                               */
/* ================================================================== */

describe('cauchyPoint — 1-D bounded', () => {
  it('snaps to lower when gradient is positive and step overshoots', () => {
    const n = 1;
    const mat = new BFGSMat(n, 5);
    const x = new Float64Array([0.5]);
    const g = new Float64Array([1]);
    const lower = new Float64Array([0]);
    const upper = new Float64Array([Infinity]);
    const nbd = classifyBounds(lower, upper);
    const ws = makeCauchyWorkspace(n, 5);

    const r = cauchyPoint(x, g, lower, upper, nbd, mat, ws);
    expect(r.ok).toBe(true);
    // With θ=1: dtMin = 1 > breakpoint 0.5 → snaps to lower bound.
    expect(ws.xc[0]).toBe(0);
    expect(r.freeCount).toBe(0);
  });

  it('stops at interior minimum when dtMin < first breakpoint', () => {
    const n = 1;
    const mat = new BFGSMat(n, 5);
    // Force θ > 1 by seeding a pair: θ = yᵀy / sᵀy.
    // With s=[1], y=[4]: θ = 16/4 = 4.
    mat.update(new Float64Array([1]), new Float64Array([4]), 0, 0);

    const x = new Float64Array([0.5]);
    const g = new Float64Array([1]);
    const lower = new Float64Array([-10]); // far from interior minimum
    const upper = new Float64Array([Infinity]);
    const nbd = classifyBounds(lower, upper);
    const ws = makeCauchyWorkspace(n, 5);

    const r = cauchyPoint(x, g, lower, upper, nbd, mat, ws);
    expect(r.ok).toBe(true);
    // Interior minimum should be within (lower, x).
    expect(ws.xc[0]).toBeGreaterThan(lower[0]);
    expect(ws.xc[0]).toBeLessThan(x[0]);
    expect(r.freeCount).toBe(1);
  });
});

/* ================================================================== */
/*  Cauchy — fixed variables                                           */
/* ================================================================== */

describe('cauchyPoint — fixed variables', () => {
  it('l == u variables stay at the bound', () => {
    const n = 3;
    const mat = new BFGSMat(n, 5);
    const x = new Float64Array([0.5, 0.5, 0.5]);
    const g = new Float64Array([1, 1, 1]);
    // Coord 1 fixed at 0.5.
    const lower = new Float64Array([-Infinity, 0.5, -Infinity]);
    const upper = new Float64Array([Infinity, 0.5, Infinity]);
    const nbd = classifyBounds(lower, upper);
    const ws = makeCauchyWorkspace(n, 5);

    // Project x first — coord 1 is already at the bound.
    const r = cauchyPoint(x, g, lower, upper, nbd, mat, ws);
    expect(r.ok).toBe(true);
    // Coord 1 must stay at 0.5 (never in free set).
    expect(ws.xc[1]).toBe(0.5);
    // Free coords 0 and 2 move along −g (both unbounded).
    expect(ws.xc[0]).toBeLessThan(x[0]);
    expect(ws.xc[2]).toBeLessThan(x[2]);
    // Free set excludes coord 1.
    const freeSetArr = Array.from(ws.freeSet.subarray(0, r.freeCount));
    expect(freeSetArr).toEqual([0, 2]);
  });
});

/* ================================================================== */
/*  Feasibility invariant                                              */
/* ================================================================== */

describe('cauchyPoint — feasibility invariant', () => {
  it('l ≤ xᶜ ≤ u under random breakpoint ordering', () => {
    const n = 6;
    const mat = new BFGSMat(n, 5);
    // Seed with two pairs to produce non-trivial W.
    mat.update(
      new Float64Array([0.1, 0.2, -0.1, 0.3, 0.05, -0.2]),
      new Float64Array([0.4, 0.1, 0.3, 0.2, 0.1, 0.1]),
      0,
      0,
    );
    mat.update(
      new Float64Array([0.05, -0.1, 0.2, 0.1, -0.05, 0.3]),
      new Float64Array([0.1, 0.3, 0.1, 0.4, 0.2, 0.15]),
      0,
      0,
    );

    const x = new Float64Array([0.3, 0.5, 0.7, 0.4, 0.6, 0.2]);
    const g = new Float64Array([0.8, -0.6, 0.4, -0.9, 0.7, 0.3]);
    const lower = new Float64Array([0, 0, 0, 0, 0, 0]);
    const upper = new Float64Array([1, 1, 1, 1, 1, 1]);
    const nbd = classifyBounds(lower, upper);
    const ws = makeCauchyWorkspace(n, 5);

    const r = cauchyPoint(x, g, lower, upper, nbd, mat, ws);
    expect(r.ok).toBe(true);
    for (let i = 0; i < n; i++) {
      expect(ws.xc[i]).toBeGreaterThanOrEqual(lower[i]);
      expect(ws.xc[i]).toBeLessThanOrEqual(upper[i]);
    }
  });
});

/* ================================================================== */
/*  Post-condition: c = Wᵀ (xᶜ − x)                                    */
/* ================================================================== */

describe('cauchyPoint — c = Wᵀ(xᶜ − x) invariant', () => {
  it('holds after multi-breakpoint sweep', () => {
    const n = 4;
    const mat = new BFGSMat(n, 5);
    mat.update(
      new Float64Array([0.2, 0.3, -0.1, 0.4]),
      new Float64Array([0.5, 0.2, 0.3, 0.1]),
      0,
      0,
    );

    const x = new Float64Array([0.5, 0.5, 0.5, 0.5]);
    const g = new Float64Array([1, -0.5, 0.8, 1.2]);
    const lower = new Float64Array([0, 0, 0, 0]);
    const upper = new Float64Array([1, 1, 1, 1]);
    const nbd = classifyBounds(lower, upper);
    const ws = makeCauchyWorkspace(n, 5);

    const r = cauchyPoint(x, g, lower, upper, nbd, mat, ws);
    expect(r.ok).toBe(true);

    const delta = new Float64Array(n);
    for (let i = 0; i < n; i++) delta[i] = ws.xc[i] - x[i];
    const cRef = new Float64Array(2 * mat.col);
    mat.applyWt(delta, cRef);
    for (let k = 0; k < 2 * mat.col; k++)
      expect(ws.c[k]).toBeCloseTo(cRef[k], 10);
  });
});

/* ================================================================== */
/*  Free-set correctness                                               */
/* ================================================================== */

describe('cauchyPoint — free-set correctness', () => {
  it('free set excludes snapped and fixed variables', () => {
    const n = 4;
    const mat = new BFGSMat(n, 5);
    const x = new Float64Array([0.1, 0.5, 0.9, 0.5]);
    const g = new Float64Array([10, 0, -10, 0]); // coord 0 snaps lower, coord 2 snaps upper
    const lower = new Float64Array([0, 0, 0, 0.5]); // coord 3 fixed at 0.5
    const upper = new Float64Array([1, 1, 1, 0.5]);
    const nbd = classifyBounds(lower, upper);
    const ws = makeCauchyWorkspace(n, 5);

    const r = cauchyPoint(x, g, lower, upper, nbd, mat, ws);
    expect(r.ok).toBe(true);
    // Coord 0: huge gradient, snaps to 0.
    expect(ws.xc[0]).toBe(0);
    // Coord 2: huge negative gradient, snaps to 1.
    expect(ws.xc[2]).toBe(1);
    // Coord 3 fixed at 0.5.
    expect(ws.xc[3]).toBe(0.5);
    // Coord 1: g=0 → no movement.
    expect(ws.xc[1]).toBe(0.5);
    // Free set: only coord 1 is strictly inside.
    const freeArr = Array.from(ws.freeSet.subarray(0, r.freeCount));
    expect(freeArr).toEqual([1]);
  });
});
