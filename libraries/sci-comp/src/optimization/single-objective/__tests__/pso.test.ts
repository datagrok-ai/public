import {PSO, applyPenalty, applyPenaltyAsync, boxConstraints} from '..';
import type {Constraint} from '..';
import {
  rosenbrock, sphere, gaussian, quadratic3d,
  productSurface, gaussianBump, expectPointClose, toAsync,
} from './helpers';

describe('PSO', () => {
  const pso = new PSO();

  describe('minimize Rosenbrock 2D → min ≈ 0 at (1, 1)', () => {
    const x0 = new Float64Array([-1.2, 1.0]);
    const settings = {swarmSize: 40, maxIterations: 3_000, seed: 42};

    it('sync', () => {
      const r = pso.minimize(rosenbrock, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeLessThan(1e-4);
      expectPointClose(r, [1, 1], 0.05);
    });

    it('async', async () => {
      const r = await pso.minimizeAsync(toAsync(rosenbrock), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeLessThan(1e-4);
      expectPointClose(r, [1, 1], 0.05);
    });
  });

  describe('minimize Sphere 3D → min ≈ 0 at (0, 0, 0)', () => {
    const x0 = new Float64Array([5, -3, 7]);
    const settings = {swarmSize: 40, maxIterations: 3_000, seed: 42};

    it('sync', () => {
      const r = pso.minimize(sphere, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 4);
      expectPointClose(r, [0, 0, 0], 1e-2);
    });

    it('async', async () => {
      const r = await pso.minimizeAsync(toAsync(sphere), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(0, 4);
      expectPointClose(r, [0, 0, 0], 1e-2);
    });
  });

  describe('maximize Gaussian 2D → max ≈ 1 at (0, 0)', () => {
    const x0 = new Float64Array([2, -3]);
    const settings = {swarmSize: 40, maxIterations: 2_000, seed: 42};

    it('sync', () => {
      const r = pso.maximize(gaussian, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(1, 2);
      expectPointClose(r, [0, 0], 0.01);
    });

    it('async', async () => {
      const r = await pso.maximizeAsync(toAsync(gaussian), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(1, 2);
      expectPointClose(r, [0, 0], 0.01);
    });
  });

  describe('maximize Product Surface → max ≈ 1/27 at (1/3, 1/3)', () => {
    const x0 = new Float64Array([0.1, 0.1]);
    const settings = {
      swarmSize: 80,
      maxIterations: 5_000,
      seed: 42,
      constraints: boxConstraints(new Float64Array([0, 0]), new Float64Array([1, 1])),
      penaltyOptions: {mu: 10_000},
    };

    it('sync', () => {
      const r = pso.maximize(productSurface, x0, settings);
      expect(r.value).toBeCloseTo(1 / 27, 2);
      expectPointClose(r, [1 / 3, 1 / 3], 0.05);
    });

    it('async', async () => {
      const r = await pso.maximizeAsync(toAsync(productSurface), x0, settings);
      expect(r.value).toBeCloseTo(1 / 27, 2);
      expectPointClose(r, [1 / 3, 1 / 3], 0.05);
    });
  });

  describe('minimize Gaussian Bump → min ≈ 0 at (0, 0)', () => {
    const x0 = new Float64Array([0.1, 0.1]);
    const settings = {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
      constraints: boxConstraints(new Float64Array([-1, -1]), new Float64Array([1, 1])),
      penaltyOptions: {mu: 10_000},
    };

    it('sync', () => {
      const r = pso.minimize(gaussianBump, x0, settings);
      expect(r.value).toBeCloseTo(0, 2);
      expectPointClose(r, [0, 0], 0.05);
    });

    it('async', async () => {
      const r = await pso.minimizeAsync(toAsync(gaussianBump), x0, settings);
      expect(r.value).toBeCloseTo(0, 2);
      expectPointClose(r, [0, 0], 0.05);
    });
  });

  describe('maximize Gaussian Bump → max ≈ 2/e at (±1, 0)', () => {
    const x0 = new Float64Array([2, 2]);
    const settings = {swarmSize: 40, maxIterations: 3_000, seed: 42};

    it('sync', () => {
      const r = pso.maximize(gaussianBump, x0, settings);
      expect(r.value).toBeCloseTo(2 / Math.E, 2);
      expect(Math.abs(r.point[0])).toBeCloseTo(1, 1);
      expect(r.point[1]).toBeCloseTo(0, 1);
    });

    it('async', async () => {
      const r = await pso.maximizeAsync(toAsync(gaussianBump), x0, settings);
      expect(r.value).toBeCloseTo(2 / Math.E, 2);
      expect(Math.abs(r.point[0])).toBeCloseTo(1, 1);
      expect(r.point[1]).toBeCloseTo(0, 1);
    });
  });

  describe('maximize Gaussian Bump (x0 = [-3, -2]) → max ≈ 2/e at (-1, 0)', () => {
    const x0 = new Float64Array([-3, -2]);
    const settings = {swarmSize: 40, maxIterations: 3_000, seed: 42};

    it('sync', () => {
      const r = pso.maximize(gaussianBump, x0, settings);
      expect(r.value).toBeCloseTo(2 / Math.E, 2);
      expectPointClose(r, [-1, 0], 0.1);
    });

    it('async', async () => {
      const r = await pso.maximizeAsync(toAsync(gaussianBump), x0, settings);
      expect(r.value).toBeCloseTo(2 / Math.E, 2);
      expectPointClose(r, [-1, 0], 0.1);
    });
  });

  describe('minimize Sphere 2D with box constraints [2,5]² → min ≈ 8 at (2, 2)', () => {
    const x0 = new Float64Array([3, 3]);
    const box = boxConstraints(new Float64Array([2, 2]), new Float64Array([5, 5]));
    const settings = {swarmSize: 40, maxIterations: 3_000, seed: 42};

    it('sync', () => {
      const boxed = applyPenalty(sphere, box, {mu: 10_000});
      const r = pso.minimize(boxed, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 1);
      expectPointClose(r, [2, 2], 0.01);
    });

    it('async', async () => {
      const boxed = applyPenaltyAsync(toAsync(sphere), box, {mu: 10_000});
      const r = await pso.minimizeAsync(boxed, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 1);
      expectPointClose(r, [2, 2], 0.01);
    });
  });

  describe('maximize Gaussian with box constraints [1,5]² via settings → max ≈ 0.1353 at (1, 1)', () => {
    const x0 = new Float64Array([2, 2]);
    const settings = {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
      constraints: boxConstraints(new Float64Array([1, 1]), new Float64Array([5, 5])),
      penaltyOptions: {mu: 10_000},
    };

    it('sync', () => {
      const r = pso.maximize(gaussian, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(Math.exp(-2), 2);
      expectPointClose(r, [1, 1], 0.01);
    });

    it('async', async () => {
      const r = await pso.maximizeAsync(toAsync(gaussian), x0, settings);
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
    const settings = {swarmSize: 80, maxIterations: 5_000, seed: 42};

    it('sync', () => {
      const constrained = applyPenalty(fn, constraints, {mu: 100_000});
      const r = pso.minimize(constrained, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2, 0);
      expectPointClose(r, [2, 2], 0.1);
    });

    it('async', async () => {
      const constrained = applyPenaltyAsync(toAsync(fn), constraints, {mu: 100_000});
      const r = await pso.minimizeAsync(constrained, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(2, 0);
      expectPointClose(r, [2, 2], 0.1);
    });
  });

  describe('minimize Sphere 2D with log-barrier [2,5]² → min ≈ 8 at (2, 2)', () => {
    const x0 = new Float64Array([3, 3]);
    const box = boxConstraints(new Float64Array([2, 2]), new Float64Array([5, 5]));
    const opts = {method: 'barrier' as const, mu: 0.01};
    const settings = {swarmSize: 40, maxIterations: 3_000, seed: 42};

    it('sync', () => {
      const barriered = applyPenalty(sphere, box, opts);
      const r = pso.minimize(barriered, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 0);
      expectPointClose(r, [2, 2], 0.05);
    });

    it('async', async () => {
      const barriered = applyPenaltyAsync(toAsync(sphere), box, opts);
      const r = await pso.minimizeAsync(barriered, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(8, 0);
      expectPointClose(r, [2, 2], 0.05);
    });
  });

  describe('minimize Quadratic 3D → min ≈ -9 at (2, 1, 1)', () => {
    const x0 = new Float64Array([0, 0, 0]);
    const settings = {swarmSize: 40, maxIterations: 3_000, seed: 42};

    it('sync', () => {
      const r = pso.minimize(quadratic3d, x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(-9, 2);
      expectPointClose(r, [2, 1, 1], 0.05);
    });

    it('async', async () => {
      const r = await pso.minimizeAsync(toAsync(quadratic3d), x0, settings);
      expect(r.converged).toBe(true);
      expect(r.value).toBeCloseTo(-9, 2);
      expectPointClose(r, [2, 1, 1], 0.05);
    });
  });

  // NOTE: The following problems from the Nelder-Mead suite are intentionally
  // excluded because PSO with default-range settings struggles with them:
  //
  //  - minimize x²+12xy+2y² s.t. 4x²+y²=25 (equality constraint on ellipse)
  //  - maximize x²+12xy+2y² s.t. 4x²+y²=25 (equality constraint on ellipse)
  //
  // PSO explores the search space stochastically, so equality constraints
  // (h(x) = 0) encoded via quadratic penalty create a narrow feasible
  // manifold that the swarm rarely samples well. The result is sensitive
  // to swarmSize, maxIterations, and the penalty coefficient μ.
  // Nelder-Mead handles these better because simplex moves can track the
  // constraint surface more precisely.

  describe('deterministic seed', () => {
    it('produces identical results with the same seed', () => {
      const settings = {swarmSize: 20, maxIterations: 500, seed: 123};
      const r1 = pso.minimize(sphere, new Float64Array([5, -3]), settings);
      const r2 = pso.minimize(sphere, new Float64Array([5, -3]), settings);
      expect(r1.value).toBe(r2.value);
      expect([...r1.point]).toEqual([...r2.point]);
    });

    it('produces different results without seed', () => {
      const settings = {swarmSize: 20, maxIterations: 100};
      const r1 = pso.minimize(sphere, new Float64Array([5, -3]), settings);
      const r2 = pso.minimize(sphere, new Float64Array([5, -3]), settings);
      const same = r1.value === r2.value && [...r1.point].every((v, i) => v === r2.point[i]);
      expect(same).toBe(false);
    });
  });

  it('throws on empty x0', () => {
    expect(() => pso.minimize(sphere, new Float64Array([]), {}))
      .toThrow('at least one element');
  });
});
