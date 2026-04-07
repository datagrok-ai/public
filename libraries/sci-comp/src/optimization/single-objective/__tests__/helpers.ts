import type {Constraint, OptimizationResult, AsyncObjectiveFunction, ObjectiveFunction} from '..';

/* ================================================================== */
/*  Test functions                                                     */
/* ================================================================== */

export const rosenbrock = (x: Float64Array): number => {
  let sum = 0;
  for (let i = 0; i < x.length - 1; i++)
    sum += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2;
  return sum;
};

export const sphere = (x: Float64Array): number => {
  let sum = 0;
  for (let i = 0; i < x.length; i++) sum += x[i] ** 2;
  return sum;
};

export const gaussian = (x: Float64Array): number =>
  Math.exp(-(x[0] ** 2 + x[1] ** 2));

/** w = 3x^2 + 5y^2 + 2z^2 - 6xy + 2yz - 6x - 6z: min = -9 at (2, 1, 1) */
export const quadratic3d = (x: Float64Array): number =>
  3 * x[0] ** 2 + 5 * x[1] ** 2 + 2 * x[2] ** 2 -
  6 * x[0] * x[1] + 2 * x[1] * x[2] -
  6 * x[0] - 6 * x[2];

/** z = xy(1 - x - y): max = 1/27 at (1/3, 1/3) */
export const productSurface = (x: Float64Array): number =>
  x[0] * x[1] * (1 - x[0] - x[1]);

/** z = x^2 + 12xy + 2y^2: constrained by 4x^2 + y^2 = 25 */
export const quadraticMixed = (x: Float64Array): number =>
  x[0] ** 2 + 12 * x[0] * x[1] + 2 * x[1] ** 2;

export const negQuadraticMixed = (x: Float64Array): number =>
  -(x[0] ** 2 + 12 * x[0] * x[1] + 2 * x[1] ** 2);

export const ellipseConstraint: Constraint[] = [
  {type: 'eq', fn: (x) => 4 * x[0] ** 2 + x[1] ** 2 - 25},
];

/** z = exp(-x^2 - y^2) * (2x^2 + y^2): min = 0 at (0,0), max = 2/e at (±1, 0) */
export const gaussianBump = (x: Float64Array): number =>
  Math.exp(-(x[0] ** 2 + x[1] ** 2)) * (2 * x[0] ** 2 + x[1] ** 2);

/* ================================================================== */
/*  Helpers                                                            */
/* ================================================================== */

/** Wrap a sync objective function as async. */
export const toAsync = (fn: ObjectiveFunction): AsyncObjectiveFunction =>
  async (x: Float64Array) => fn(x);

/** Check that every component of result.point is close to the expected value. */
export function expectPointClose(r: OptimizationResult, expected: number[], tol: number) {
  expect(r.point.length).toBe(expected.length);
  for (let i = 0; i < expected.length; i++)
    expect(r.point[i]).toBeCloseTo(expected[i], -Math.log10(tol));
}
