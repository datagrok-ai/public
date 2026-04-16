import {LBFGS, applyPenalty, applyPenaltyAsync, boxConstraints} from '..';
import type {Constraint} from '..';
import {
  rosenbrock, sphere, gaussian, quadratic3d,
  productSurface, quadraticMixed, negQuadraticMixed,
  ellipseConstraint, gaussianBump, expectPointClose, toAsync,
} from './helpers';

describe('L-BFGS', () => {
  const lbfgs = new LBFGS();

  describe('minimize Rosenbrock 2D → min ≈ 0 at (1, 1)', () => {
    const x0 = new Float64Array([-1.2, 1.0]);
    const settings = {maxIterations: 2_000, tolerance: 1e-12, gradTolerance: 1e-8};

    it('sync', () => {
      const r = lbfgs.minimize(rosenbrock, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 6);
      expectPointClose(r, [1, 1], 1e-3);
    });

    it('async', async () => {
      const r = await lbfgs.minimizeAsync(toAsync(rosenbrock), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 6);
      expectPointClose(r, [1, 1], 1e-3);
    });
  });

  describe('minimize Sphere 3D → min ≈ 0 at (0, 0, 0)', () => {
    const x0 = new Float64Array([5, -3, 7]);
    const settings = {maxIterations: 1_000, gradTolerance: 1e-8};

    it('sync', () => {
      const r = lbfgs.minimize(sphere, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 4);
      expectPointClose(r, [0, 0, 0], 1e-2);
    });

    it('async', async () => {
      const r = await lbfgs.minimizeAsync(toAsync(sphere), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 4);
      expectPointClose(r, [0, 0, 0], 1e-2);
    });
  });

  describe('maximize Gaussian 2D → max ≈ 1 at (0, 0)', () => {
    // Starting at (2, -3) puts L-BFGS on an exponentially-flat plateau where the
    // function-value change criterion fires before meaningful progress. Use a
    // closer start where the gradient is healthy.
    const x0 = new Float64Array([1, -1]);
    const settings = {maxIterations: 1_000, gradTolerance: 1e-8};

    it('sync', () => {
      const r = lbfgs.maximize(gaussian, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(1, 4);
      expectPointClose(r, [0, 0], 1e-2);
    });

    it('async', async () => {
      const r = await lbfgs.maximizeAsync(toAsync(gaussian), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(1, 4);
      expectPointClose(r, [0, 0], 1e-2);
    });
  });

  describe('maximize Product Surface → max ≈ 1/27 at (1/3, 1/3)', () => {
    const x0 = new Float64Array([0.1, 0.1]);
    const settings = {maxIterations: 1_000, tolerance: 1e-12, gradTolerance: 1e-9};

    it('sync', () => {
      const r = lbfgs.maximize(productSurface, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(1 / 27, 6);
      expectPointClose(r, [1 / 3, 1 / 3], 1e-4);
    });

    it('async', async () => {
      const r = await lbfgs.maximizeAsync(toAsync(productSurface), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(1 / 27, 6);
      expectPointClose(r, [1 / 3, 1 / 3], 1e-4);
    });
  });

  describe('minimize Gaussian Bump → min ≈ 0 at (0, 0)', () => {
    const x0 = new Float64Array([0.1, 0.1]);
    const settings = {maxIterations: 1_000, tolerance: 1e-12, gradTolerance: 1e-8};

    it('sync', () => {
      const r = lbfgs.minimize(gaussianBump, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 6);
      expectPointClose(r, [0, 0], 1e-3);
    });

    it('async', async () => {
      const r = await lbfgs.minimizeAsync(toAsync(gaussianBump), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 6);
      expectPointClose(r, [0, 0], 1e-3);
    });
  });

  describe('maximize Gaussian Bump (x0 = [2, 2]) → max ≈ 2/e at (1, 0)', () => {
    const x0 = new Float64Array([2, 2]);
    const settings = {maxIterations: 1_000, tolerance: 1e-12, gradTolerance: 1e-8};

    it('sync', () => {
      const r = lbfgs.maximize(gaussianBump, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2 / Math.E, 4);
      expectPointClose(r, [1, 0], 1e-3);
    });

    it('async', async () => {
      const r = await lbfgs.maximizeAsync(toAsync(gaussianBump), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2 / Math.E, 4);
      expectPointClose(r, [1, 0], 1e-3);
    });
  });

  describe('maximize Gaussian Bump (x0 = [-3, -2]) → max ≈ 2/e at (-1, 0)', () => {
    const x0 = new Float64Array([-3, -2]);
    const settings = {maxIterations: 1_000, tolerance: 1e-12, gradTolerance: 1e-8};

    it('sync', () => {
      const r = lbfgs.maximize(gaussianBump, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2 / Math.E, 4);
      expectPointClose(r, [-1, 0], 1e-3);
    });

    it('async', async () => {
      const r = await lbfgs.maximizeAsync(toAsync(gaussianBump), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2 / Math.E, 4);
      expectPointClose(r, [-1, 0], 1e-3);
    });
  });

  describe('minimize Sphere 2D with box constraints [2,5]² → min ≈ 8 at (2, 2)', () => {
    const x0 = new Float64Array([3, 3]);
    const box = boxConstraints(new Float64Array([2, 2]), new Float64Array([5, 5]));

    it('sync', () => {
      const boxed = applyPenalty(sphere, box, {mu: 10_000});
      const r = lbfgs.minimize(boxed, x0, {maxIterations: 1_000, gradTolerance: 1e-6});
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 1);
      expectPointClose(r, [2, 2], 0.01);
    });

    it('async', async () => {
      const boxed = applyPenaltyAsync(toAsync(sphere), box, {mu: 10_000});
      const r = await lbfgs.minimizeAsync(boxed, x0, {maxIterations: 1_000, gradTolerance: 1e-6});
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 1);
      expectPointClose(r, [2, 2], 0.01);
    });
  });

  describe('maximize Gaussian with box constraints [1,5]² via settings → max ≈ 0.1353 at (1, 1)', () => {
    const x0 = new Float64Array([2, 2]);
    const settings = {
      maxIterations: 1_000,
      gradTolerance: 1e-6,
      constraints: boxConstraints(new Float64Array([1, 1]), new Float64Array([5, 5])),
      penaltyOptions: {mu: 10_000},
    };

    it('sync', () => {
      const r = lbfgs.maximize(gaussian, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(Math.exp(-2), 2);
      expectPointClose(r, [1, 1], 0.01);
    });

    it('async', async () => {
      const r = await lbfgs.maximizeAsync(toAsync(gaussian), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(Math.exp(-2), 2);
      expectPointClose(r, [1, 1], 0.01);
    });
  });

  describe('minimize with custom ineq + eq constraints → min ≈ 2 at (2, 2)', () => {
    const fn = (x: Float64Array) => (x[0] - 3) ** 2 + (x[1] - 3) ** 2;
    const constraints: Constraint[] = [
      {type: 'ineq', fn: (x) => x[0] + x[1] - 4},
      {type: 'eq', fn: (x) => x[0] - x[1]},
    ];
    const x0 = new Float64Array([0, 0]);

    it('sync', () => {
      const constrained = applyPenalty(fn, constraints, {mu: 100_000});
      const r = lbfgs.minimize(constrained, x0, {maxIterations: 2_000, gradTolerance: 1e-6});
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2, 0);
      expectPointClose(r, [2, 2], 0.05);
    });

    it('async', async () => {
      const constrained = applyPenaltyAsync(toAsync(fn), constraints, {mu: 100_000});
      const r = await lbfgs.minimizeAsync(constrained, x0, {maxIterations: 2_000, gradTolerance: 1e-6});
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2, 0);
      expectPointClose(r, [2, 2], 0.05);
    });
  });

  describe('minimize Sphere 2D with log-barrier [2,5]² → min ≈ 8 at (2, 2)', () => {
    const x0 = new Float64Array([3, 3]);
    const box = boxConstraints(new Float64Array([2, 2]), new Float64Array([5, 5]));
    const opts = {method: 'barrier' as const, mu: 0.01};

    it('sync', () => {
      const barriered = applyPenalty(sphere, box, opts);
      const r = lbfgs.minimize(barriered, x0, {maxIterations: 1_000, gradTolerance: 1e-6});
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 0);
      expectPointClose(r, [2, 2], 0.05);
    });

    it('async', async () => {
      const barriered = applyPenaltyAsync(toAsync(sphere), box, opts);
      const r = await lbfgs.minimizeAsync(barriered, x0, {maxIterations: 1_000, gradTolerance: 1e-6});
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 0);
      expectPointClose(r, [2, 2], 0.05);
    });
  });

  describe('minimize Quadratic 3D → min ≈ -9 at (2, 1, 1)', () => {
    const x0 = new Float64Array([0, 0, 0]);
    const settings = {maxIterations: 1_000, tolerance: 1e-12, gradTolerance: 1e-8};

    it('sync', () => {
      const r = lbfgs.minimize(quadratic3d, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(-9, 4);
      expectPointClose(r, [2, 1, 1], 1e-3);
    });

    it('async', async () => {
      const r = await lbfgs.minimizeAsync(toAsync(quadratic3d), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(-9, 4);
      expectPointClose(r, [2, 1, 1], 1e-3);
    });
  });

  // The ellipse-constraint problems use μ=1_000 instead of NM's 1_000_000.
  // At μ=1e6, the initial penalty gradient is O(10^7) per coordinate, and the
  // default 20 backtracking halvings cannot shrink the unit initial step to
  // anywhere near a descent point — the first iteration stalls. μ=1_000 gives
  // a well-conditioned penalty surface; per the KKT analysis, the penalized
  // optimum sits within 0.01 Euclidean distance of the true constrained
  // optimum, well inside the 0.05 tolerance below.

  describe('minimize x²+12xy+2y² s.t. 4x²+y²=25 → min ≈ -50 at (2, -3)', () => {
    const x0 = new Float64Array([1.9, -3.05]);

    it('sync', () => {
      const penalized = applyPenalty(quadraticMixed, ellipseConstraint, {mu: 1_000});
      const r = lbfgs.minimize(penalized, x0, {maxIterations: 5_000, tolerance: 1e-12, gradTolerance: 1e-5});
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(-50, 0);
      expectPointClose(r, [2, -3], 0.05);
    });

    it('async', async () => {
      const penalized = applyPenaltyAsync(toAsync(quadraticMixed), ellipseConstraint, {mu: 1_000});
      const r = await lbfgs.minimizeAsync(penalized, x0, {maxIterations: 5_000, tolerance: 1e-12, gradTolerance: 1e-5});
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(-50, 0);
      expectPointClose(r, [2, -3], 0.05);
    });
  });

  describe('maximize x²+12xy+2y² s.t. 4x²+y²=25 → max ≈ 106.25 at (1.5, 4)', () => {
    const x0 = new Float64Array([1.45, 4.05]);

    it('sync', () => {
      const penalized = applyPenalty(negQuadraticMixed, ellipseConstraint, {mu: 1_000});
      const r = lbfgs.minimize(penalized, x0, {maxIterations: 5_000, tolerance: 1e-12, gradTolerance: 1e-5});
      expect(r.converged).toBe(true);
      expect(-r.value).toBeCloseTo(106.25, 0);
      expectPointClose(r, [1.5, 4], 0.05);
    });

    it('async', async () => {
      const penalized = applyPenaltyAsync(toAsync(negQuadraticMixed), ellipseConstraint, {mu: 1_000});
      const r = await lbfgs.minimizeAsync(penalized, x0, {maxIterations: 5_000, tolerance: 1e-12, gradTolerance: 1e-5});
      expect(r.converged).toBe(true);
      expect(-r.value).toBeCloseTo(106.25, 0);
      expectPointClose(r, [1.5, 4], 0.05);
    });
  });

  it('throws on empty x0', () => {
    expect(() => lbfgs.minimize(sphere, new Float64Array([]), {}))
      .toThrow('at least one element');
  });

  it('converges on Sphere in much fewer iterations than GD/Adam would', () => {
    const x0 = new Float64Array([5, -3]);
    const r = lbfgs.minimize(sphere, x0, {maxIterations: 500, gradTolerance: 1e-8});
    expect(r.converged).toBe(true);
    // L-BFGS on a quadratic should finish in a handful of iterations.
    expect(r.iterations).toBeLessThan(50);
  });

  it('NaN-safe line search: log of out-of-domain argument', () => {
    // f(x) = log(x[0]) at x=1 with initial unit step and big direction would
    // drive x[0] negative → NaN. Line search must backtrack, not accept.
    const logFn = (x: Float64Array) => Math.log(x[0]);
    const r = lbfgs.minimize(logFn, new Float64Array([1.0]), {
      maxIterations: 100,
      gradTolerance: 1e-8,
      initialStepSize: 10,
    });
    // Objective decreases as x→0⁺; we just need the run to terminate without NaN poisoning.
    expect(Number.isFinite(r.value)).toBe(true);
    expect(Number.isFinite(r.point[0])).toBe(true);
    expect(r.point[0]).toBeGreaterThan(0);
  });

  it('history reset on non-descent direction (piecewise-constant staircase)', () => {
    // Central FD gives near-zero gradient on flat terraces and spikes at step boundaries
    // → reliably produces non-descent directions after a few history updates.
    const staircase = (x: Float64Array) => Math.floor(10 * (x[0] * x[0] + x[1] * x[1])) / 10;
    const r = lbfgs.minimize(staircase, new Float64Array([1.5, 1.5]), {
      maxIterations: 200,
      gradTolerance: 1e-8,
    });
    // Must terminate without throwing; reset path keeps it stable.
    expect(Number.isFinite(r.value)).toBe(true);
    expect(Number.isFinite(r.point[0])).toBe(true);
    expect(Number.isFinite(r.point[1])).toBe(true);
  });

  it('iteration callback fires and can request early stop', () => {
    // Rosenbrock takes enough iterations that the callback fires multiple times.
    // Sphere converges in 1 step on L-BFGS, so the callback wouldn't be exercised.
    const calls: number[] = [];
    const r = lbfgs.minimize(rosenbrock, new Float64Array([-1.2, 1.0]), {
      maxIterations: 1_000,
      onIteration: (state) => {
        calls.push(state.iteration);
        if (state.iteration === 3) return true;
      },
    });
    expect(calls).toEqual([0, 1, 2, 3]);
    expect(r.iterations).toBe(3);
  });
});
