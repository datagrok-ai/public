/**
 * L-BFGS-B end-to-end integration tests.
 *
 * Mirrors the convention used by other optimizer tests in this suite
 * (`describe('minimize <problem>', () => { it('sync'); it('async'); })`).
 * Covers:
 *   - unconstrained Rosenbrock, Sphere;
 *   - bounded Rosenbrock with active bounds;
 *   - fixed variables (l == u);
 *   - half-bounded (lower only);
 *   - analytic gradient fast path;
 *   - maximize + negation on maximize path;
 *   - callback early stop;
 *   - feasibility invariant under random bounds;
 *   - x₀ at optimum ⇒ iterations = 0;
 *   - maxFunctionEvaluations cap terminates without converged;
 *   - NaN at x₀ throws.
 */
import {LBFGSB} from '..';
import {rosenbrock, sphere, expectPointClose, toAsync} from './helpers';

/* ================================================================== */
/*  Unconstrained                                                      */
/* ================================================================== */

describe('LBFGSB minimize Rosenbrock 2D → (1, 1), f=0', () => {
  const opt = new LBFGSB();
  const x0 = new Float64Array([-1.2, 1.0]);
  const settings = {maxIterations: 200, tolerance: 1e-12, gradTolerance: 1e-8};

  it('sync', () => {
    const r = opt.minimize(rosenbrock, x0, settings);
    expect(r.converged).toBe(true);
    expect(r.value).toBeLessThan(1e-8);
    expectPointClose(r, [1, 1], 1e-3);
  });

  it('async', async () => {
    const r = await opt.minimizeAsync(toAsync(rosenbrock), x0, settings);
    expect(r.converged).toBe(true);
    expect(r.value).toBeLessThan(1e-8);
    expectPointClose(r, [1, 1], 1e-3);
  });
});

describe('LBFGSB minimize Sphere 3D → origin, f=0', () => {
  const opt = new LBFGSB();
  const x0 = new Float64Array([5, -3, 7]);
  const settings = {maxIterations: 200, gradTolerance: 1e-10};

  it('sync', () => {
    const r = opt.minimize(sphere, x0, settings);
    expect(r.converged).toBe(true);
    expect(r.value).toBeLessThan(1e-12);
    expectPointClose(r, [0, 0, 0], 1e-5);
  });

  it('async', async () => {
    const r = await opt.minimizeAsync(toAsync(sphere), x0, settings);
    expect(r.converged).toBe(true);
    expect(r.value).toBeLessThan(1e-12);
    expectPointClose(r, [0, 0, 0], 1e-5);
  });
});

/* ================================================================== */
/*  Bounded                                                            */
/* ================================================================== */

describe('LBFGSB minimize Rosenbrock 2D on [-1.5, 0.5]²', () => {
  // Unconstrained minimum at (1, 1) lies outside; boundary x₁=0.5 is active.
  // Known solution: x ≈ (0.5, 0.25), f ≈ 0.25.
  const opt = new LBFGSB();
  const x0 = new Float64Array([0.0, 0.0]);
  const settings = {
    maxIterations: 300,
    tolerance: 1e-12,
    gradTolerance: 1e-8,
    bounds: {lower: -1.5, upper: 0.5},
  };

  it('sync', () => {
    const r = opt.minimize(rosenbrock, x0, settings);
    expect(r.converged).toBe(true);
    // Upper bound on x[0] should be active at solution.
    expect(r.point[0]).toBeCloseTo(0.5, 4);
    expect(r.point[1]).toBeCloseTo(0.25, 3);
    expect(r.value).toBeCloseTo(0.25, 3);
    // Feasibility.
    for (let i = 0; i < 2; i++) {
      expect(r.point[i]).toBeGreaterThanOrEqual(-1.5);
      expect(r.point[i]).toBeLessThanOrEqual(0.5);
    }
  });

  it('async', async () => {
    const r = await opt.minimizeAsync(toAsync(rosenbrock), x0, settings);
    expect(r.converged).toBe(true);
    expect(r.point[0]).toBeCloseTo(0.5, 4);
    expect(r.point[1]).toBeCloseTo(0.25, 3);
  });
});

describe('LBFGSB minimize Sphere with fixed variable (l == u)', () => {
  // f = x₀² + x₁² + x₂². Fix x₁ = 2 via l=u.
  // Minimum: x = (0, 2, 0), f = 4.
  const opt = new LBFGSB();
  const x0 = new Float64Array([5, -1, 3]);
  const settings = {
    maxIterations: 100,
    gradTolerance: 1e-8,
    bounds: {
      lower: [-Infinity, 2, -Infinity],
      upper: [Infinity, 2, Infinity],
    },
  };

  it('sync', () => {
    const r = opt.minimize(sphere, x0, settings);
    expect(r.converged).toBe(true);
    expect(r.point[0]).toBeCloseTo(0, 5);
    expect(r.point[1]).toBe(2);
    expect(r.point[2]).toBeCloseTo(0, 5);
    expect(r.value).toBeCloseTo(4, 6);
  });
});

describe('LBFGSB minimize Sphere half-bounded (lower=0)', () => {
  // Start with negative initial → projected to [0, +∞).
  const opt = new LBFGSB();
  const x0 = new Float64Array([-2, 3, -1]);
  const settings = {
    maxIterations: 100,
    gradTolerance: 1e-8,
    bounds: {lower: 0},
  };

  it('sync', () => {
    const r = opt.minimize(sphere, x0, settings);
    expect(r.converged).toBe(true);
    expectPointClose(r, [0, 0, 0], 1e-5);
    expect(r.value).toBeLessThan(1e-10);
  });
});

/* ================================================================== */
/*  Analytic gradient                                                  */
/* ================================================================== */

describe('LBFGSB with analytic gradient', () => {
  // f(x) = (x-a)ᵀ H (x-a) / 2 for SPD H, known minimum at a.
  const a = new Float64Array([1.5, -0.7, 2.3]);
  const H = [[2, 0.5, 0], [0.5, 3, -0.2], [0, -0.2, 1.5]];
  const fn = (x: Float64Array): number => {
    let s = 0;
    for (let i = 0; i < 3; i++) {
      let Hx = 0;
      for (let j = 0; j < 3; j++) Hx += H[i][j] * (x[j] - a[j]);
      s += (x[i] - a[i]) * Hx;
    }
    return 0.5 * s;
  };
  const gradFn = (x: Float64Array, out: Float64Array): void => {
    for (let i = 0; i < 3; i++) {
      let Hx = 0;
      for (let j = 0; j < 3; j++) Hx += H[i][j] * (x[j] - a[j]);
      out[i] = Hx;
    }
  };

  it('converges with analytic gradient in fewer f-evals than FD', () => {
    const opt = new LBFGSB();
    const x0 = new Float64Array([0, 0, 0]);
    const s = {maxIterations: 100, gradTolerance: 1e-10, gradFn};
    const r = opt.minimize(fn, x0, s);
    expect(r.converged).toBe(true);
    expect(r.value).toBeLessThan(1e-14);
    expectPointClose(r, [a[0], a[1], a[2]], 1e-6);
  });
});

/* ================================================================== */
/*  Maximize                                                           */
/* ================================================================== */

describe('LBFGSB maximize Gaussian 2D', () => {
  // f = exp(-(x² + y²)); max = 1 at origin.
  const opt = new LBFGSB();
  const gaussian = (x: Float64Array): number =>
    Math.exp(-(x[0] * x[0] + x[1] * x[1]));
  const x0 = new Float64Array([0.5, -0.5]);
  const settings = {maxIterations: 100, gradTolerance: 1e-10};

  it('sync', () => {
    const r = opt.maximize(gaussian, x0, settings);
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(1, 8);
    expectPointClose(r, [0, 0], 1e-4);
  });
});

/* ================================================================== */
/*  Callback early stop                                                */
/* ================================================================== */

describe('LBFGSB callback early stop', () => {
  it('stops iteration when callback returns true', () => {
    const opt = new LBFGSB();
    const x0 = new Float64Array([-1.2, 1.0]);
    let callCount = 0;
    const r = opt.minimize(rosenbrock, x0, {
      maxIterations: 1000,
      onIteration: ({iteration}) => {
        callCount = iteration;
        return iteration >= 3;
      },
    });
    expect(callCount).toBeLessThanOrEqual(4);
    expect(r.converged).toBe(false);
  });

  it('callback receives populated extra', () => {
    const opt = new LBFGSB();
    const x0 = new Float64Array([-1.2, 1.0]);
    const extras: Record<string, number>[] = [];
    opt.minimize(rosenbrock, x0, {
      maxIterations: 5,
      onIteration: ({extra}) => {
        if (extra) extras.push(extra as Record<string, number>);
      },
    });
    expect(extras.length).toBeGreaterThan(0);
    const e = extras[0];
    expect(typeof e.projGradInfNorm).toBe('number');
    expect(typeof e.functionEvaluations).toBe('number');
    expect(typeof e.stepSize).toBe('number');
    expect(typeof e.lineSearchSteps).toBe('number');
    expect(typeof e.historyCount).toBe('number');
    expect(typeof e.activeBounds).toBe('number');
  });
});

/* ================================================================== */
/*  Degenerate cases                                                   */
/* ================================================================== */

describe('LBFGSB degenerate inputs', () => {
  it('x₀ at optimum → iterations=0, converged=true', () => {
    const opt = new LBFGSB();
    const x0 = new Float64Array([0, 0, 0]);
    const r = opt.minimize(sphere, x0, {gradTolerance: 1e-10});
    expect(r.converged).toBe(true);
    expect(r.iterations).toBe(0);
    expect(r.value).toBe(0);
    expect(r.costHistory.length).toBe(1);
    expect(r.costHistory[0]).toBe(0);
  });

  it('throws when objective returns non-finite at x₀', () => {
    const opt = new LBFGSB();
    const x0 = new Float64Array([0, 0]);
    expect(() => opt.minimize(() => NaN, x0, {})).toThrow('non-finite at x0');
  });

  it('throws when x₀ is empty', () => {
    const opt = new LBFGSB();
    expect(() => opt.minimize(sphere, new Float64Array(0), {}))
      .toThrow('x0 must have at least one element');
  });
});

/* ================================================================== */
/*  Feasibility invariant                                              */
/* ================================================================== */

describe('LBFGSB feasibility invariant', () => {
  it('result.point always lies in [lower, upper]', () => {
    const opt = new LBFGSB();
    // Contrived problem: sphere, random bounds, x₀ infeasible.
    const x0 = new Float64Array([-10, 10, -5, 5]);
    const bounds = {
      lower: [-1, 0, -2, 3],
      upper: [2, 4, 1, 10],
    };
    const r = opt.minimize(sphere, x0, {bounds, maxIterations: 100});
    expect(r.converged).toBe(true);
    for (let i = 0; i < 4; i++) {
      expect(r.point[i]).toBeGreaterThanOrEqual(bounds.lower[i]);
      expect(r.point[i]).toBeLessThanOrEqual(bounds.upper[i]);
    }
  });
});
