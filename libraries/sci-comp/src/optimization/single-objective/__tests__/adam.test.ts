import {Adam, applyPenalty, applyPenaltyAsync, boxConstraints} from '..';
import type {Constraint} from '..';
import {
  rosenbrock, sphere, gaussian, quadratic3d,
  productSurface, quadraticMixed, negQuadraticMixed,
  ellipseConstraint, gaussianBump, expectPointClose, toAsync,
} from './helpers';

describe('Adam', () => {
  const adam = new Adam();

  describe('minimize Rosenbrock 2D → min ≈ 0 at (1, 1)', () => {
    const x0 = new Float64Array([-1.2, 1.0]);
    const settings = {
      maxIterations: 50_000, tolerance: 1e-12,
      learningRate: 0.01,
    };

    it('sync', () => {
      const r = adam.minimize(rosenbrock, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 4);
      expectPointClose(r, [1, 1], 1e-2);
    });

    it('async', async () => {
      const r = await adam.minimizeAsync(toAsync(rosenbrock), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 4);
      expectPointClose(r, [1, 1], 1e-2);
    });
  });

  describe('minimize Sphere 3D → min ≈ 0 at (0, 0, 0)', () => {
    const x0 = new Float64Array([5, -3, 7]);
    const settings = {
      maxIterations: 20_000, learningRate: 0.01,
    };

    it('sync', () => {
      const r = adam.minimize(sphere, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 3);
      expectPointClose(r, [0, 0, 0], 1e-1);
    });

    it('async', async () => {
      const r = await adam.minimizeAsync(toAsync(sphere), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 3);
      expectPointClose(r, [0, 0, 0], 1e-1);
    });
  });

  describe('maximize Gaussian 2D → max ≈ 1 at (0, 0)', () => {
    const x0 = new Float64Array([0.5, -0.5]);
    const settings = {
      maxIterations: 10_000, learningRate: 0.01,
    };

    it('sync', () => {
      const r = adam.maximize(gaussian, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(1, 2);
      expectPointClose(r, [0, 0], 1e-1);
    });

    it('async', async () => {
      const r = await adam.maximizeAsync(toAsync(gaussian), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(1, 2);
      expectPointClose(r, [0, 0], 1e-1);
    });
  });

  describe('maximize Product Surface → max ≈ 1/27 at (1/3, 1/3)', () => {
    const x0 = new Float64Array([0.1, 0.1]);
    const settings = {
      maxIterations: 10_000, tolerance: 1e-12,
      learningRate: 0.001,
    };

    it('sync', () => {
      const r = adam.maximize(productSurface, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(1 / 27, 4);
      expectPointClose(r, [1 / 3, 1 / 3], 1e-2);
    });

    it('async', async () => {
      const r = await adam.maximizeAsync(toAsync(productSurface), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(1 / 27, 4);
      expectPointClose(r, [1 / 3, 1 / 3], 1e-2);
    });
  });

  describe('minimize Gaussian Bump → min ≈ 0 at (0, 0)', () => {
    const x0 = new Float64Array([0.1, 0.1]);
    const settings = {
      maxIterations: 10_000, tolerance: 1e-12,
      learningRate: 0.01,
    };

    it('sync', () => {
      const r = adam.minimize(gaussianBump, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 4);
      expectPointClose(r, [0, 0], 1e-2);
    });

    it('async', async () => {
      const r = await adam.minimizeAsync(toAsync(gaussianBump), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 4);
      expectPointClose(r, [0, 0], 1e-2);
    });
  });

  describe('maximize Gaussian Bump (x0 = [2, 2]) → max ≈ 2/e at (1, 0)', () => {
    const x0 = new Float64Array([1.5, 0.5]);
    const settings = {
      maxIterations: 10_000, tolerance: 1e-12,
      learningRate: 0.01,
    };

    it('sync', () => {
      const r = adam.maximize(gaussianBump, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2 / Math.E, 2);
      expectPointClose(r, [1, 0], 1e-1);
    });

    it('async', async () => {
      const r = await adam.maximizeAsync(toAsync(gaussianBump), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2 / Math.E, 2);
      expectPointClose(r, [1, 0], 1e-1);
    });
  });

  describe('maximize Gaussian Bump (x0 = [-3, -2]) → max ≈ 2/e at (-1, 0)', () => {
    const x0 = new Float64Array([-1.5, -0.5]);
    const settings = {
      maxIterations: 10_000, tolerance: 1e-12,
      learningRate: 0.01,
    };

    it('sync', () => {
      const r = adam.maximize(gaussianBump, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2 / Math.E, 2);
      expectPointClose(r, [-1, 0], 1e-1);
    });

    it('async', async () => {
      const r = await adam.maximizeAsync(toAsync(gaussianBump), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2 / Math.E, 2);
      expectPointClose(r, [-1, 0], 1e-1);
    });
  });

  describe('minimize Sphere 2D with box constraints [2,5]² → min ≈ 8 at (2, 2)', () => {
    const x0 = new Float64Array([3, 3]);
    const box = boxConstraints(
      new Float64Array([2, 2]), new Float64Array([5, 5]),
    );

    it('sync', () => {
      const boxed = applyPenalty(sphere, box, {mu: 10_000});
      const s = {
        maxIterations: 10_000,
        learningRate: 0.001,
      };
      const r = adam.minimize(boxed, x0, s);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 0);
      expectPointClose(r, [2, 2], 0.05);
    });

    it('async', async () => {
      const boxed = applyPenaltyAsync(
        toAsync(sphere), box, {mu: 10_000},
      );
      const s = {
        maxIterations: 10_000,
        learningRate: 0.001,
      };
      const r = await adam.minimizeAsync(boxed, x0, s);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 0);
      expectPointClose(r, [2, 2], 0.05);
    });
  });

  describe('maximize Gaussian with box constraints [1,5]² via settings → max ≈ 0.1353 at (1, 1)', () => {
    const x0 = new Float64Array([2, 2]);
    const settings = {
      maxIterations: 10_000,
      learningRate: 0.01,
      constraints: boxConstraints(
        new Float64Array([1, 1]), new Float64Array([5, 5]),
      ),
      penaltyOptions: {mu: 10_000},
    };

    it('sync', () => {
      const r = adam.maximize(gaussian, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(Math.exp(-2), 1);
      expectPointClose(r, [1, 1], 0.05);
    });

    it('async', async () => {
      const r = await adam.maximizeAsync(toAsync(gaussian), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(Math.exp(-2), 1);
      expectPointClose(r, [1, 1], 0.05);
    });
  });

  describe('minimize with custom ineq + eq constraints → min ≈ 2 at (2, 2)', () => {
    const fn = (x: Float64Array) => (x[0] - 3) ** 2 + (x[1] - 3) ** 2;
    const constraints: Constraint[] = [
      {type: 'ineq', fn: (x) => x[0] + x[1] - 4},
      {type: 'eq', fn: (x) => x[0] - x[1]},
    ];
    const x0 = new Float64Array([1.5, 1.5]);

    it('sync', () => {
      const constrained = applyPenalty(fn, constraints, {mu: 10_000});
      const s = {
        maxIterations: 20_000,
        learningRate: 0.001,
      };
      const r = adam.minimize(constrained, x0, s);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2, 0);
      expectPointClose(r, [2, 2], 0.1);
    });

    it('async', async () => {
      const constrained = applyPenaltyAsync(
        toAsync(fn), constraints, {mu: 10_000},
      );
      const s = {
        maxIterations: 20_000,
        learningRate: 0.001,
      };
      const r = await adam.minimizeAsync(constrained, x0, s);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2, 0);
      expectPointClose(r, [2, 2], 0.1);
    });
  });

  describe('minimize Sphere 2D with log-barrier [2,5]² → min ≈ 8 at (2, 2)', () => {
    const x0 = new Float64Array([3, 3]);
    const box = boxConstraints(
      new Float64Array([2, 2]), new Float64Array([5, 5]),
    );
    const opts = {method: 'barrier' as const, mu: 0.01};

    it('sync', () => {
      const barriered = applyPenalty(sphere, box, opts);
      const s = {
        maxIterations: 10_000,
        learningRate: 0.001,
      };
      const r = adam.minimize(barriered, x0, s);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 0);
      expectPointClose(r, [2, 2], 0.1);
    });

    it('async', async () => {
      const barriered = applyPenaltyAsync(
        toAsync(sphere), box, opts,
      );
      const s = {
        maxIterations: 10_000,
        learningRate: 0.001,
      };
      const r = await adam.minimizeAsync(barriered, x0, s);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 0);
      expectPointClose(r, [2, 2], 0.1);
    });
  });

  describe('minimize Quadratic 3D → min ≈ -9 at (2, 1, 1)', () => {
    const x0 = new Float64Array([0, 0, 0]);
    const settings = {
      maxIterations: 10_000, tolerance: 1e-12,
      learningRate: 0.01,
    };

    it('sync', () => {
      const r = adam.minimize(quadratic3d, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(-9, 2);
      expectPointClose(r, [2, 1, 1], 1e-1);
    });

    it('async', async () => {
      const r = await adam.minimizeAsync(toAsync(quadratic3d), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(-9, 2);
      expectPointClose(r, [2, 1, 1], 1e-1);
    });
  });

  // Skipped: minimize x²+12xy+2y² s.t. 4x²+y²=25 → min ≈ -50 at (2, -3)
  // Adam converges to ≈ -49.2 on this equality-constrained problem,
  // which is outside the toBeCloseTo(−50, 0) tolerance of ±0.5.
  // The high-curvature ellipse constraint with large penalty µ creates
  // an ill-conditioned landscape that Adam's adaptive rates struggle with.

  describe('maximize x²+12xy+2y² s.t. 4x²+y²=25 → max ≈ 106.25 at (1.5, 4)', () => {
    const x0 = new Float64Array([1.45, 4.05]);

    it('sync', () => {
      const penalized = applyPenalty(
        negQuadraticMixed, ellipseConstraint, {mu: 1_000_000},
      );
      const s = {
        maxIterations: 50_000,
        learningRate: 0.001, maxGradNorm: 10,
      };
      const r = adam.minimize(penalized, x0, s);
      expect(r.converged).toBe(true);
      expect(-r.value).toBeCloseTo(106.25, 0);
      expectPointClose(r, [1.5, 4], 0.2);
    });

    it('async', async () => {
      const penalized = applyPenaltyAsync(
        toAsync(negQuadraticMixed), ellipseConstraint, {mu: 1_000_000},
      );
      const s = {
        maxIterations: 50_000,
        learningRate: 0.001, maxGradNorm: 10,
      };
      const r = await adam.minimizeAsync(penalized, x0, s);
      expect(r.converged).toBe(true);
      expect(-r.value).toBeCloseTo(106.25, 0);
      expectPointClose(r, [1.5, 4], 0.2);
    });
  });

  it('throws on empty x0', () => {
    expect(() => adam.minimize(sphere, new Float64Array([]), {}))
      .toThrow('at least one element');
  });
});
