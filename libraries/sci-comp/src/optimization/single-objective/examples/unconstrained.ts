import {NelderMead, PSO} from '..';
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

const nm = new NelderMead();
const pso = new PSO();

console.log('========== MINIMIZE: UNCONSTRAINED ==========');

printProblem([
  'f(x) = sum[ 100*(x[i+1] - x[i]^2)^2 + (1 - x[i])^2 ]  (Rosenbrock)',
  'goal:     minimize',
  'x0:       [-1.2, 1.0]',
  'expected: min = 0 at (1, 1)',
]);

printResult(
  'Nelder-Mead → Rosenbrock 2D',
  nm.minimize(rosenbrock, new Float64Array([-1.2, 1.0]), {
    maxIterations: 5_000,
    tolerance: 1e-12,
  }),
);

printResult(
  'PSO → Rosenbrock 2D',
  pso.minimize(rosenbrock, new Float64Array([-1.2, 1.0]), {
    swarmSize: 40,
    maxIterations: 3_000,
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
  'Nelder-Mead → Sphere 3D',
  nm.minimize(sphere, new Float64Array([5, -3, 7]), {maxIterations: 5_000}),
);

/* ================================================================== */
/*  2. maximize                                                        */
/* ================================================================== */

console.log('\n========== MAXIMIZE: UNCONSTRAINED ==========');

printProblem([
  'f(x, y) = exp(-(x^2 + y^2))  (Gaussian)',
  'goal:     maximize',
  'x0:       [2, -3]',
  'expected: max = 1 at (0, 0)',
]);

printResult(
  'Nelder-Mead → maximize Gaussian 2D',
  nm.maximize(gaussian, new Float64Array([2, -3]), {
    maxIterations: 5_000,
  }),
);

printResult(
  'PSO → maximize Gaussian 2D',
  pso.maximize(gaussian, new Float64Array([2, -3]), {
    swarmSize: 40,
    maxIterations: 2_000,
  }),
);
