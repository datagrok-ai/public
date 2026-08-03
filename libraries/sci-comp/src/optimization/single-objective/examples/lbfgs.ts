import {LBFGS, GradientDescent, Adam, boxConstraints, applyPenalty} from '..';
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

const lbfgs = new LBFGS();

console.log('========== L-BFGS: MINIMIZE UNCONSTRAINED ==========');

printProblem([
  'f(x) = sum[ 100*(x[i+1] - x[i]^2)^2 + (1 - x[i])^2 ]  (Rosenbrock)',
  'goal:     minimize',
  'x0:       [-1.2, 1.0]',
  'expected: min = 0 at (1, 1)',
]);

printResult(
  'L-BFGS → Rosenbrock 2D',
  lbfgs.minimize(rosenbrock, new Float64Array([-1.2, 1.0]), {
    maxIterations: 1_000,
    tolerance: 1e-12,
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
  'L-BFGS → Sphere 3D',
  lbfgs.minimize(sphere, new Float64Array([5, -3, 7]), {maxIterations: 500}),
);

/* ================================================================== */
/*  2. maximize                                                        */
/* ================================================================== */

console.log('\n========== L-BFGS: MAXIMIZE UNCONSTRAINED ==========');

printProblem([
  'f(x, y) = exp(-(x^2 + y^2))  (Gaussian)',
  'goal:     maximize',
  'x0:       [2, -3]',
  'expected: max = 1 at (0, 0)',
]);

printResult(
  'L-BFGS → maximize Gaussian 2D',
  lbfgs.maximize(gaussian, new Float64Array([2, -3]), {maxIterations: 500}),
);

/* ================================================================== */
/*  3. constrained minimize                                            */
/* ================================================================== */

console.log('\n========== L-BFGS: MINIMIZE CONSTRAINED ==========');

printProblem([
  'f(x, y) = x^2 + y^2  (Sphere 2D)',
  'subject to: 2 ≤ x, y ≤ 5',
  'goal:     minimize',
  'x0:       [3, 3]',
  'expected: min = 8 at (2, 2)',
]);

const box = boxConstraints(new Float64Array([2, 2]), new Float64Array([5, 5]));
const boxedSphere = applyPenalty(sphere, box, {mu: 10_000});

printResult(
  'L-BFGS → Sphere 2D with box [2,5]²',
  lbfgs.minimize(boxedSphere, new Float64Array([3, 3]), {
    maxIterations: 1_000,
    gradTolerance: 1e-6,
  }),
);

/* ================================================================== */
/*  4. historySize tuning                                              */
/* ================================================================== */

console.log('\n========== L-BFGS: historySize TUNING ==========');

printProblem([
  'Same Rosenbrock problem, varying history size m ∈ {3, 10, 20}.',
  'Larger m = better Hessian approximation, more memory + per-iter cost.',
]);

for (const m of [3, 10, 20]) {
  printResult(
    `L-BFGS (historySize=${m}) → Rosenbrock 2D`,
    lbfgs.minimize(rosenbrock, new Float64Array([-1.2, 1.0]), {
      maxIterations: 1_000,
      historySize: m,
      tolerance: 1e-12,
    }),
  );
}

/* ================================================================== */
/*  5. iteration count vs GD / Adam                                    */
/* ================================================================== */

console.log('\n========== L-BFGS vs GradientDescent vs Adam ==========');
console.log('  Same problem, same termination tolerance — compare iteration counts.');

const gd = new GradientDescent();
const adam = new Adam();

printProblem([
  '',
  'f(x, y, z) = x^2 + y^2 + z^2  (Sphere 3D)',
  'x0:       [5, -3, 7]',
]);

printResult(
  'L-BFGS',
  lbfgs.minimize(sphere, new Float64Array([5, -3, 7]), {
    maxIterations: 10_000,
    gradTolerance: 1e-6,
  }),
);

printResult(
  'GradientDescent',
  gd.minimize(sphere, new Float64Array([5, -3, 7]), {
    maxIterations: 10_000,
    learningRate: 0.1,
  }),
);

printResult(
  'Adam',
  adam.minimize(sphere, new Float64Array([5, -3, 7]), {
    maxIterations: 10_000,
    learningRate: 0.1,
  }),
);
