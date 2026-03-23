import {getOptimizer, listOptimizers} from '..';
import type {OptimizationResult} from '..';

/* ================================================================== */
/*  Test functions                                                     */
/* ================================================================== */

const sphere = (x: Float64Array): number => {
  let sum = 0;
  for (let i = 0; i < x.length; i++) sum += x[i] ** 2;
  return sum;
};

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
/*  Registry usage                                                     */
/* ================================================================== */

console.log('========== REGISTRY ==========');
console.log('Available optimizers:', listOptimizers().join(', '));

const optimizer = getOptimizer('nelder-mead');

printProblem([
  'f(x, y, z) = x^2 + y^2 + z^2  (Sphere 3D)',
  'goal:     minimize',
  'x0:       [5, -3, 7]',
  'expected: min = 0 at (0, 0, 0)',
]);

printResult(
  'minimize via registry → Sphere 3D',
  optimizer.minimize(sphere, new Float64Array([5, -3, 7]), {maxIterations: 5_000}),
);

printProblem([
  '',
  'f(x, y) = exp(-(x^2 + y^2))  (Gaussian)',
  'goal:     maximize',
  'x0:       [4, -2]',
  'expected: max = 1 at (0, 0)',
]);

printResult(
  'maximize via registry → Gaussian 2D',
  optimizer.maximize(gaussian, new Float64Array([4, -2]), {maxIterations: 5_000}),
);
