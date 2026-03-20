import {NelderMead, applyPenalty, boxConstraints} from '..';
import type {OptimizationResult, Constraint} from '..';

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
  console.log(`  value:      ${r.value.toExponential(6)}`);
  console.log(`  point:      [${[...r.point].map((v) => v.toFixed(6)).join(', ')}]`);
  console.log(`  iterations: ${r.iterations}`);
  console.log(`  converged:  ${r.converged}`);
}

/* ================================================================== */
/*  1. minimize — box constraints (quadratic penalty)                  */
/* ================================================================== */

const nm = new NelderMead();

console.log('========== MINIMIZE: BOX CONSTRAINTS (quadratic penalty) ==========');

printProblem([
  'f(x, y) = x^2 + y^2  (Sphere)',
  'goal:        minimize',
  'subject to:  2 <= x <= 5, 2 <= y <= 5',
  'method:      quadratic penalty, mu = 10000',
  'x0:          [3, 3]',
  'expected:    min = 8 at (2, 2)',
]);

const boxed = applyPenalty(
  sphere,
  boxConstraints(
    new Float64Array([2, 2]),
    new Float64Array([5, 5]),
  ),
  {mu: 10_000},
);

printResult(
  'Nelder-Mead → Sphere 2D, box [2,5]^2',
  nm.minimize(boxed, new Float64Array([3, 3]), {maxIterations: 5_000}),
);

/* ================================================================== */
/*  2. maximize — box constraints (via settings)                       */
/* ================================================================== */

console.log('\n========== MAXIMIZE: BOX CONSTRAINTS (via settings) ==========');

printProblem([
  'f(x, y) = exp(-(x^2 + y^2))  (Gaussian)',
  'goal:        maximize',
  'subject to:  1 <= x <= 5, 1 <= y <= 5',
  'method:      quadratic penalty (via settings), mu = 10000',
  'x0:          [2, 2]',
  'expected:    max = exp(-2) ~ 0.1353 at (1, 1)',
]);

printResult(
  'Nelder-Mead → maximize Gaussian, box [1,5]^2',
  nm.maximize(gaussian, new Float64Array([2, 2]), {
    maxIterations: 5_000,
    constraints: boxConstraints(
      new Float64Array([1, 1]),
      new Float64Array([5, 5]),
    ),
    penaltyOptions: {mu: 10_000},
  }),
);

/* ================================================================== */
/*  3. minimize — custom inequality + equality constraints             */
/* ================================================================== */

console.log('\n========== MINIMIZE: CUSTOM CONSTRAINTS ==========');

printProblem([
  'f(x, y) = (x - 3)^2 + (y - 3)^2',
  'goal:        minimize',
  'subject to:  x + y <= 4  (inequality)',
  '             x = y       (equality)',
  'method:      quadratic penalty, mu = 100000',
  'x0:          [0, 0]',
  'expected:    min = 2 at (2, 2)',
]);

const quadratic = (x: Float64Array) => (x[0] - 3) ** 2 + (x[1] - 3) ** 2;

const constraints: Constraint[] = [
  {type: 'ineq', fn: (x) => x[0] + x[1] - 4},
  {type: 'eq', fn: (x) => x[0] - x[1]},
];

const constrained = applyPenalty(quadratic, constraints, {mu: 100_000});

printResult(
  'Nelder-Mead → custom ineq + eq',
  nm.minimize(constrained, new Float64Array([0, 0]), {maxIterations: 10_000}),
);

/* ================================================================== */
/*  4. minimize — log-barrier                                          */
/* ================================================================== */

console.log('\n========== MINIMIZE: LOG-BARRIER ==========');

printProblem([
  'f(x, y) = x^2 + y^2  (Sphere)',
  'goal:        minimize',
  'subject to:  2 <= x <= 5, 2 <= y <= 5',
  'method:      log-barrier, mu = 0.01',
  'x0:          [3, 3]  (must be strictly feasible)',
  'expected:    min ~ 8 at (2, 2)',
]);

const barriered = applyPenalty(
  sphere,
  boxConstraints(
    new Float64Array([2, 2]),
    new Float64Array([5, 5]),
  ),
  {method: 'barrier', mu: 0.01},
);

printResult(
  'Nelder-Mead → Sphere 2D, barrier, box [2,5]^2',
  nm.minimize(barriered, new Float64Array([3, 3]), {maxIterations: 5_000}),
);
