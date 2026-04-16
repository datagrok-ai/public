import type {ObjectiveFunction} from '../types';

/**
 * Shared suite of classical test functions for unconstrained-optimization
 * benchmarks. Each function exposes its global minimum in TSDoc; the
 * accompanying `TestProblem` record bundles a function with its known
 * optimum for reporting.
 */

/** Sphere: f(x) = Σ xᵢ²  — global min f(0,…,0) = 0 (unimodal, convex). */
export const sphere: ObjectiveFunction = (x) => {
  let s = 0;
  for (let i = 0; i < x.length; i++) s += x[i] * x[i];
  return s;
};

/** Rosenbrock: f(x) = Σᵢ [100·(xᵢ₊₁ − xᵢ²)² + (1 − xᵢ)²] — min f(1,…,1) = 0. */
export const rosenbrock: ObjectiveFunction = (x) => {
  let s = 0;
  for (let i = 0; i < x.length - 1; i++)
    s += 100 * (x[i + 1] - x[i] * x[i]) ** 2 + (1 - x[i]) ** 2;
  return s;
};

/** Beale: min f(3, 0.5) = 0. */
export const beale: ObjectiveFunction = (x) =>
  (1.5 - x[0] + x[0] * x[1]) ** 2 +
  (2.25 - x[0] + x[0] * x[1] * x[1]) ** 2 +
  (2.625 - x[0] + x[0] * x[1] * x[1] * x[1]) ** 2;

/** Booth: min f(1, 3) = 0. */
export const booth: ObjectiveFunction = (x) =>
  (x[0] + 2 * x[1] - 7) ** 2 + (2 * x[0] + x[1] - 5) ** 2;

/** Matyas: min f(0, 0) = 0. */
export const matyas: ObjectiveFunction = (x) =>
  0.26 * (x[0] * x[0] + x[1] * x[1]) - 0.48 * x[0] * x[1];

/** Himmelblau: 4 equivalent minima, all f = 0. */
export const himmelblau: ObjectiveFunction = (x) =>
  (x[0] * x[0] + x[1] - 11) ** 2 + (x[0] + x[1] * x[1] - 7) ** 2;

/** Three-Hump Camel: min f(0, 0) = 0. */
export const threeHumpCamel: ObjectiveFunction = (x) =>
  2 * x[0] ** 2 - 1.05 * x[0] ** 4 + x[0] ** 6 / 6 + x[0] * x[1] + x[1] ** 2;

/** Rastrigin: min f(0,…,0) = 0 (highly multimodal). */
export const rastrigin: ObjectiveFunction = (x) => {
  const A = 10;
  let s = A * x.length;
  for (let i = 0; i < x.length; i++)
    s += x[i] * x[i] - A * Math.cos(2 * Math.PI * x[i]);
  return s;
};

/** Ackley: min f(0, 0) = 0 (multimodal). */
export const ackley: ObjectiveFunction = (x) =>
  -20 * Math.exp(-0.2 * Math.sqrt(0.5 * (x[0] ** 2 + x[1] ** 2))) -
  Math.exp(0.5 * (Math.cos(2 * Math.PI * x[0]) + Math.cos(2 * Math.PI * x[1]))) +
  Math.E + 20;

/** Lévi N.13: min f(1, 1) = 0 (multimodal). */
export const levi13: ObjectiveFunction = (x) =>
  Math.sin(3 * Math.PI * x[0]) ** 2 +
  (x[0] - 1) ** 2 * (1 + Math.sin(3 * Math.PI * x[1]) ** 2) +
  (x[1] - 1) ** 2 * (1 + Math.sin(2 * Math.PI * x[1]) ** 2);

/** Griewank: min f(0,…,0) = 0 (multimodal, many shallow ripples). */
export const griewank: ObjectiveFunction = (x) => {
  let sum = 0;
  let prod = 1;
  for (let i = 0; i < x.length; i++) {
    sum += x[i] * x[i];
    prod *= Math.cos(x[i] / Math.sqrt(i + 1));
  }
  return 1 + sum / 4000 - prod;
};

/** Styblinski-Tang: min f(−2.9035,…) ≈ −39.16617·n (multimodal). */
export const styblinskiTang: ObjectiveFunction = (x) => {
  let s = 0;
  for (let i = 0; i < x.length; i++)
    s += x[i] ** 4 - 16 * x[i] ** 2 + 5 * x[i];
  return s / 2;
};

/** Easom: min f(π, π) = −1 (nearly flat outside narrow peak). */
export const easom: ObjectiveFunction = (x) =>
  -Math.cos(x[0]) * Math.cos(x[1]) *
  Math.exp(-((x[0] - Math.PI) ** 2 + (x[1] - Math.PI) ** 2));

/** Goldstein-Price: min f(0, −1) = 3 (multimodal). */
export const goldsteinPrice: ObjectiveFunction = (x) => {
  const a = 1 + (x[0] + x[1] + 1) ** 2 *
    (19 - 14 * x[0] + 3 * x[0] ** 2 - 14 * x[1] + 6 * x[0] * x[1] + 3 * x[1] ** 2);
  const b = 30 + (2 * x[0] - 3 * x[1]) ** 2 *
    (18 - 32 * x[0] + 12 * x[0] ** 2 + 48 * x[1] - 36 * x[0] * x[1] + 27 * x[1] ** 2);
  return a * b;
};

/** McCormick: min f(−0.54719, −1.54719) ≈ −1.9133. */
export const mccormick: ObjectiveFunction = (x) =>
  Math.sin(x[0] + x[1]) + (x[0] - x[1]) ** 2 - 1.5 * x[0] + 2.5 * x[1] + 1;

/** The four Himmelblau minima (all f = 0), for nearest-optimum distance reporting. */
export const HIMMELBLAU_MINIMA: ReadonlyArray<Float64Array> = [
  new Float64Array([3, 2]),
  new Float64Array([-2.805118, 3.131312]),
  new Float64Array([-3.779310, -3.283186]),
  new Float64Array([3.584428, -1.848126]),
];
