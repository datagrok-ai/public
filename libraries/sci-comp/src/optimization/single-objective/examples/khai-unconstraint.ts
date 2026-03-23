import {NelderMead, PSO} from '..';
import type {OptimizationResult} from '..';

/* ================================================================== */
/*  Test functions                                                     */
/* ================================================================== */

/** w = 3x^2 + 5y^2 + 2z^2 - 6xy + 2yz - 6x - 6z: min = -9 at (2, 1, 1) */
const quadratic3d = (x: Float64Array): number =>
  3 * x[0] ** 2 + 5 * x[1] ** 2 + 2 * x[2] ** 2 -
  6 * x[0] * x[1] + 2 * x[1] * x[2] -
  6 * x[0] - 6 * x[2];

/** z = xy(1 - x - y): max = 1/27 at (1/3, 1/3) */
const productSurface = (x: Float64Array): number =>
  x[0] * x[1] * (1 - x[0] - x[1]);

/** z = exp(-x^2 - y^2) * (2x^2 + y^2): min = 0 at (0,0), max = 2/e at (±1, 0) */
const gaussianBump = (x: Float64Array): number =>
  Math.exp(-(x[0] ** 2 + x[1] ** 2)) * (2 * x[0] ** 2 + x[1] ** 2);

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
/*  1. minimize quadratic 3D                                           */
/* ================================================================== */

const nm = new NelderMead();
const pso = new PSO();

console.log('========== MINIMIZE: QUADRATIC 3D ==========');

printProblem([
  'w(x,y,z) = 3x^2 + 5y^2 + 2z^2 - 6xy + 2yz - 6x - 6z',
  'goal:     minimize',
  'x0:       [0, 0, 0]',
  'expected: min = -9 at (2, 1, 1)',
]);

printResult(
  'Nelder-Mead → Quadratic 3D',
  nm.minimize(quadratic3d, new Float64Array([0, 0, 0]), {
    maxIterations: 5_000,
    tolerance: 1e-12,
  }),
);

printResult(
  'PSO → Quadratic 3D',
  pso.minimize(quadratic3d, new Float64Array([0, 0, 0]), {
    swarmSize: 40,
    maxIterations: 3_000,
  }),
);

/* ================================================================== */
/*  2. maximize product surface                                        */
/* ================================================================== */

console.log('\n========== MAXIMIZE: PRODUCT SURFACE ==========');

printProblem([
  'z(x,y) = xy(1 - x - y)',
  'goal:     maximize',
  'x0:       [0.1, 0.1]',
  'expected: max = 1/27 ≈ 0.037037 at (1/3, 1/3)',
]);

printResult(
  'Nelder-Mead → maximize Product Surface',
  nm.maximize(productSurface, new Float64Array([0.1, 0.1]), {
    maxIterations: 5_000,
    tolerance: 1e-12,
  }),
);

printResult(
  'PSO → maximize Product Surface',
  pso.maximize(productSurface, new Float64Array([0.1, 0.1]), {
    swarmSize: 40,
    maxIterations: 3_000,
  }),
);

/* ================================================================== */
/*  3. gaussian bump: exp(-x^2 - y^2) * (2x^2 + y^2)                  */
/* ================================================================== */

console.log('\n========== MINIMIZE: GAUSSIAN BUMP ==========');

printProblem([
  'z(x,y) = exp(-x^2 - y^2) * (2x^2 + y^2)',
  'goal:     minimize',
  'x0:       [0.1, 0.1]',
  'expected: min = 0 at (0, 0)',
]);

printResult(
  'Nelder-Mead → minimize Gaussian Bump',
  nm.minimize(gaussianBump, new Float64Array([0.1, 0.1]), {
    maxIterations: 5_000,
    tolerance: 1e-12,
  }),
);

printResult(
  'PSO → minimize Gaussian Bump',
  pso.minimize(gaussianBump, new Float64Array([0.1, 0.1]), {
    swarmSize: 40,
    maxIterations: 3_000,
  }),
);

console.log('\n========== MAXIMIZE: GAUSSIAN BUMP (x0 = [2, 2]) ==========');

printProblem([
  'z(x,y) = exp(-x^2 - y^2) * (2x^2 + y^2)',
  'goal:     maximize',
  'x0:       [2, 2]',
  'expected: max = 2/e ≈ 0.7358 at (1, 0)',
]);

printResult(
  'Nelder-Mead → maximize Gaussian Bump (x0 = [2, 2])',
  nm.maximize(gaussianBump, new Float64Array([2, 2]), {
    maxIterations: 5_000,
    tolerance: 1e-12,
  }),
);

console.log('\n========== MAXIMIZE: GAUSSIAN BUMP (x0 = [-3, -2]) ==========');

printProblem([
  'z(x,y) = exp(-x^2 - y^2) * (2x^2 + y^2)',
  'goal:     maximize',
  'x0:       [-3, -2]',
  'expected: max = 2/e ≈ 0.7358 at (-1, 0)',
]);

printResult(
  'Nelder-Mead → maximize Gaussian Bump (x0 = [-3, -2])',
  nm.maximize(gaussianBump, new Float64Array([-3, -2]), {
    maxIterations: 5_000,
    tolerance: 1e-12,
  }),
);

printResult(
  'PSO → maximize Gaussian Bump (x0 = [-3, -2])',
  pso.maximize(gaussianBump, new Float64Array([-3, -2]), {
    swarmSize: 40,
    maxIterations: 3_000,
  }),
);
