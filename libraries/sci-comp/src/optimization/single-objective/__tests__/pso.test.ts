import {PSO, applyPenalty, boxConstraints} from '..';
import type {Constraint} from '..';
import {
  rosenbrock, sphere, gaussian, quadratic3d,
  productSurface, gaussianBump, expectPointClose,
} from './helpers';

describe('PSO', () => {
  const pso = new PSO();

  it('minimize Rosenbrock 2D → min ≈ 0 at (1, 1)', () => {
    const r = pso.minimize(rosenbrock, new Float64Array([-1.2, 1.0]), {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeLessThan(1e-4);
    expectPointClose(r, [1, 1], 0.05);
  });

  it('minimize Sphere 3D → min = 0 at (0, 0, 0)', () => {
    const r = pso.minimize(sphere, new Float64Array([5, -3, 7]), {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(0, 4);
    expectPointClose(r, [0, 0, 0], 1e-2);
  });

  it('maximize Gaussian 2D → max = 1 at (0, 0)', () => {
    const r = pso.maximize(gaussian, new Float64Array([2, -3]), {
      swarmSize: 40,
      maxIterations: 2_000,
      seed: 42,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(1, 2);
    expectPointClose(r, [0, 0], 0.01);
  });

  it('maximize Product Surface → max = 1/27 at (1/3, 1/3)', () => {
    const r = pso.maximize(productSurface, new Float64Array([0.1, 0.1]), {
      swarmSize: 80,
      maxIterations: 5_000,
      seed: 42,
      constraints: boxConstraints(new Float64Array([0, 0]), new Float64Array([1, 1])),
      penaltyOptions: {mu: 10_000},
    });
    expect(r.value).toBeCloseTo(1 / 27, 2);
    expectPointClose(r, [1 / 3, 1 / 3], 0.05);
  });

  it('minimize Gaussian Bump → min = 0 at (0, 0)', () => {
    const r = pso.minimize(gaussianBump, new Float64Array([0.1, 0.1]), {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
      constraints: boxConstraints(new Float64Array([-1, -1]), new Float64Array([1, 1])),
      penaltyOptions: {mu: 10_000},
    });
    expect(r.value).toBeCloseTo(0, 2);
    expectPointClose(r, [0, 0], 0.05);
  });

  it('maximize Gaussian Bump → max = 2/e at (±1, 0)', () => {
    const r = pso.maximize(gaussianBump, new Float64Array([2, 2]), {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
    });
    expect(r.value).toBeCloseTo(2 / Math.E, 2);
    expect(Math.abs(r.point[0])).toBeCloseTo(1, 1);
    expect(r.point[1]).toBeCloseTo(0, 1);
  });

  it('minimize Sphere 2D with box constraints [2,5]² → min = 8 at (2, 2)', () => {
    const boxed = applyPenalty(
      sphere,
      boxConstraints(new Float64Array([2, 2]), new Float64Array([5, 5])),
      {mu: 10_000},
    );
    const r = pso.minimize(boxed, new Float64Array([3, 3]), {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(8, 1);
    expectPointClose(r, [2, 2], 0.01);
  });

  it('maximize Gaussian with box constraints [1,5]² via settings → max ≈ 0.1353 at (1, 1)', () => {
    const r = pso.maximize(gaussian, new Float64Array([2, 2]), {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
      constraints: boxConstraints(new Float64Array([1, 1]), new Float64Array([5, 5])),
      penaltyOptions: {mu: 10_000},
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(Math.exp(-2), 2);
    expectPointClose(r, [1, 1], 0.01);
  });

  it('minimize with custom ineq + eq constraints → min = 2 at (2, 2)', () => {
    const fn = (x: Float64Array) => (x[0] - 3) ** 2 + (x[1] - 3) ** 2;
    const constraints: Constraint[] = [
      {type: 'ineq', fn: (x) => x[0] + x[1] - 4},
      {type: 'eq', fn: (x) => x[0] - x[1]},
    ];
    const constrained = applyPenalty(fn, constraints, {mu: 100_000});
    const r = pso.minimize(constrained, new Float64Array([0, 0]), {
      swarmSize: 80,
      maxIterations: 5_000,
      seed: 42,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(2, 0);
    expectPointClose(r, [2, 2], 0.1);
  });

  it('minimize Sphere 2D with log-barrier [2,5]² → min ≈ 8 at (2, 2)', () => {
    const barriered = applyPenalty(
      sphere,
      boxConstraints(new Float64Array([2, 2]), new Float64Array([5, 5])),
      {method: 'barrier', mu: 0.01},
    );
    const r = pso.minimize(barriered, new Float64Array([3, 3]), {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(8, 0);
    expectPointClose(r, [2, 2], 0.05);
  });

  it('minimize Quadratic 3D → min = -9 at (2, 1, 1)', () => {
    const r = pso.minimize(quadratic3d, new Float64Array([0, 0, 0]), {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(-9, 2);
    expectPointClose(r, [2, 1, 1], 0.05);
  });

  it('produces deterministic results with the same seed', () => {
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
    // With Math.random two runs are extremely unlikely to match exactly
    const same = r1.value === r2.value && [...r1.point].every((v, i) => v === r2.point[i]);
    expect(same).toBe(false);
  });

  it('throws on empty x0', () => {
    expect(() => pso.minimize(sphere, new Float64Array([]), {}))
      .toThrow('at least one element');
  });
});
