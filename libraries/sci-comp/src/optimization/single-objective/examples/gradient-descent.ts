import {GradientDescent, applyPenalty, boxConstraints} from '..';
import type {OptimizationResult} from '..';

/* ================================================================== */
/*  Test functions                                                     */
/* ================================================================== */

/** Rosenbrock: min = 0 at (1, 1, ..., 1) */
const rosenbrock = (x: Float64Array): number => {
  let sum = 0;
  for (let i = 0; i < x.length - 1; i++)
    sum += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2;
  return sum;
};

/** Sphere: min = 0 at (0, 0, ..., 0) */
const sphere = (x: Float64Array): number => {
  let sum = 0;
  for (let i = 0; i < x.length; i++) sum += x[i] ** 2;
  return sum;
};

/** Gaussian: max = 1 at (0, 0) */
const gaussian = (x: Float64Array): number =>
  Math.exp(-(x[0] ** 2 + x[1] ** 2));

/* ================================================================== */
/*  Helper                                                             */
/* ================================================================== */

function printProblem(lines: string[]) {
  for (const line of lines)
    console.log(`  ${line}`);
}

function printResult(label: string, r: OptimizationResult) {
  console.log(`\n--- ${label} ---`);
  console.log(`  value:      ${r.value.toFixed(6)}`);
  console.log(`  point:      [${[...r.point].map((v) => v.toFixed(6)).join(', ')}]`);
  console.log(`  iterations: ${r.iterations}`);
  console.log(`  converged:  ${r.converged}`);
  if (!r.converged) console.warn('  WARNING: optimizer did not converge');
}

/* ================================================================== */
/*  1. minimize                                                        */
/* ================================================================== */

const gd = new GradientDescent();

console.log('========== GRADIENT DESCENT: MINIMIZE ==========');

printProblem([
  'f(x) = sum[ 100*(x[i+1] - x[i]^2)^2 + (1 - x[i])^2 ]  (Rosenbrock)',
  'goal:     minimize',
  'x0:       [-1.2, 1.0]',
  'expected: min = 0 at (1, 1)',
]);

printResult(
  'GradientDescent → Rosenbrock 2D',
  gd.minimize(rosenbrock, new Float64Array([-1.2, 1.0]), {
    maxIterations: 50_000,
    tolerance: 1e-12,
    learningRate: 0.001,
    momentum: 0.9,
  }),
);

printProblem([
  '',
  'f(x, y, z) = x^2 + y^2 + z^2  (Sphere 3D)',
  'goal:     minimize',
  'x0:       [5, -3, 7]',
  'expected: min = 0 at (0, 0, 0)',
]);

printResult(
  'GradientDescent → Sphere 3D',
  gd.minimize(sphere, new Float64Array([5, -3, 7]), {
    maxIterations: 10_000,
    learningRate: 0.01,
    momentum: 0.9,
  }),
);

/* ================================================================== */
/*  2. maximize                                                        */
/* ================================================================== */

console.log('\n========== GRADIENT DESCENT: MAXIMIZE ==========');

printProblem([
  'f(x, y) = exp(-(x^2 + y^2))  (Gaussian)',
  'goal:     maximize',
  'x0:       [2, -3]',
  'expected: max = 1 at (0, 0)',
]);

printResult(
  'GradientDescent → maximize Gaussian 2D',
  gd.maximize(gaussian, new Float64Array([2, -3]), {
    maxIterations: 10_000,
    learningRate: 0.1,
    momentum: 0.9,
  }),
);

/* ================================================================== */
/*  3. constrained                                                     */
/* ================================================================== */

console.log('\n========== GRADIENT DESCENT: CONSTRAINED ==========');

printProblem([
  'f(x, y) = x^2 + y^2  (Sphere)',
  'subject to: 2 <= x <= 5, 2 <= y <= 5',
  'goal:     minimize',
  'x0:       [3, 3]',
  'expected: min = 8 at (2, 2)',
]);

const box = boxConstraints(new Float64Array([2, 2]), new Float64Array([5, 5]));
const boxed = applyPenalty(sphere, box, {mu: 10_000});

printResult(
  'GradientDescent → Sphere 2D (box [2,5]²)',
  gd.minimize(boxed, new Float64Array([3, 3]), {
    maxIterations: 10_000,
    learningRate: 0.0001,
    momentum: 0.9,
  }),
);
