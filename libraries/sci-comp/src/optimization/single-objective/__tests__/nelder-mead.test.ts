import {NelderMead, applyPenalty, boxConstraints} from '..';
import type {Constraint} from '..';
import {
  rosenbrock, sphere, gaussian, quadratic3d,
  productSurface, quadraticMixed, negQuadraticMixed,
  ellipseConstraint, gaussianBump, expectPointClose,
} from './helpers';

describe('Nelder-Mead', () => {
  const nm = new NelderMead();

  it('minimize Rosenbrock 2D → min = 0 at (1, 1)', () => {
    const r = nm.minimize(rosenbrock, new Float64Array([-1.2, 1.0]), {
      maxIterations: 5_000,
      tolerance: 1e-12,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(0, 6);
    expectPointClose(r, [1, 1], 1e-4);
  });

  it('minimize Sphere 3D → min = 0 at (0, 0, 0)', () => {
    const r = nm.minimize(sphere, new Float64Array([5, -3, 7]), {
      maxIterations: 5_000,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(0, 4);
    expectPointClose(r, [0, 0, 0], 1e-2);
  });

  it('maximize Gaussian 2D → max = 1 at (0, 0)', () => {
    const r = nm.maximize(gaussian, new Float64Array([2, -3]), {
      maxIterations: 5_000,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(1, 4);
    expectPointClose(r, [0, 0], 1e-2);
  });

  it('maximize Product Surface → max = 1/27 at (1/3, 1/3)', () => {
    const r = nm.maximize(productSurface, new Float64Array([0.1, 0.1]), {
      maxIterations: 5_000,
      tolerance: 1e-12,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(1 / 27, 6);
    expectPointClose(r, [1 / 3, 1 / 3], 1e-4);
  });

  it('minimize Gaussian Bump → min = 0 at (0, 0)', () => {
    const r = nm.minimize(gaussianBump, new Float64Array([0.1, 0.1]), {
      maxIterations: 5_000,
      tolerance: 1e-12,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(0, 6);
    expectPointClose(r, [0, 0], 1e-4);
  });

  it('maximize Gaussian Bump (x0 = [2, 2]) → max = 2/e at (1, 0)', () => {
    const r = nm.maximize(gaussianBump, new Float64Array([2, 2]), {
      maxIterations: 5_000,
      tolerance: 1e-12,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(2 / Math.E, 4);
    expectPointClose(r, [1, 0], 1e-3);
  });

  it('maximize Gaussian Bump (x0 = [-3, -2]) → max = 2/e at (-1, 0)', () => {
    const r = nm.maximize(gaussianBump, new Float64Array([-3, -2]), {
      maxIterations: 5_000,
      tolerance: 1e-12,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(2 / Math.E, 4);
    expectPointClose(r, [-1, 0], 1e-3);
  });

  it('minimize Sphere 2D with box constraints [2,5]² → min = 8 at (2, 2)', () => {
    const boxed = applyPenalty(
      sphere,
      boxConstraints(new Float64Array([2, 2]), new Float64Array([5, 5])),
      {mu: 10_000},
    );
    const r = nm.minimize(boxed, new Float64Array([3, 3]), {
      maxIterations: 5_000,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(8, 1);
    expectPointClose(r, [2, 2], 0.01);
  });

  it('maximize Gaussian with box constraints [1,5]² via settings → max ≈ 0.1353 at (1, 1)', () => {
    const r = nm.maximize(gaussian, new Float64Array([2, 2]), {
      maxIterations: 5_000,
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
    const r = nm.minimize(constrained, new Float64Array([0, 0]), {
      maxIterations: 10_000,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(2, 0);
    expectPointClose(r, [2, 2], 0.05);
  });

  it('minimize Sphere 2D with log-barrier [2,5]² → min ≈ 8 at (2, 2)', () => {
    const barriered = applyPenalty(
      sphere,
      boxConstraints(new Float64Array([2, 2]), new Float64Array([5, 5])),
      {method: 'barrier', mu: 0.01},
    );
    const r = nm.minimize(barriered, new Float64Array([3, 3]), {
      maxIterations: 5_000,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(8, 0);
    expectPointClose(r, [2, 2], 0.05);
  });

  it('minimize Quadratic 3D → min = -9 at (2, 1, 1)', () => {
    const r = nm.minimize(quadratic3d, new Float64Array([0, 0, 0]), {
      maxIterations: 5_000,
      tolerance: 1e-12,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(-9, 4);
    expectPointClose(r, [2, 1, 1], 1e-3);
  });

  it('minimize x²+12xy+2y² s.t. 4x²+y²=25 → min = -50 at (2, -3)', () => {
    const penalized = applyPenalty(quadraticMixed, ellipseConstraint, {mu: 1_000_000});
    const r = nm.minimize(penalized, new Float64Array([1, -2]), {
      maxIterations: 10_000,
      tolerance: 1e-12,
    });
    expect(r.converged).toBe(true);
    expect(r.value).toBeCloseTo(-50, 0);
    expectPointClose(r, [2, -3], 0.05);
  });

  it('maximize x²+12xy+2y² s.t. 4x²+y²=25 → max = 106.25 at (1.5, 4)', () => {
    const penalized = applyPenalty(negQuadraticMixed, ellipseConstraint, {mu: 1_000_000});
    const r = nm.minimize(penalized, new Float64Array([1, 2]), {
      maxIterations: 10_000,
      tolerance: 1e-12,
    });
    expect(r.converged).toBe(true);
    expect(-r.value).toBeCloseTo(106.25, 0);
    expectPointClose(r, [1.5, 4], 0.05);
  });

  it('throws on empty x0', () => {
    expect(() => nm.minimize(sphere, new Float64Array([]), {}))
      .toThrow('at least one element');
  });
});
