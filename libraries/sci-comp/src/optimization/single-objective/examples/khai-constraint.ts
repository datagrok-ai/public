import {NelderMead, PSO, applyPenalty} from '..';
import type {OptimizationResult, Constraint} from '..';

/* ================================================================== */
/*  Test functions                                                     */
/* ================================================================== */

/** z = x^2 + 12xy + 2y^2 */
const quadraticMixed = (x: Float64Array): number =>
  x[0] ** 2 + 12 * x[0] * x[1] + 2 * x[1] ** 2;

/** Constraint: 4x^2 + y^2 = 25 */
const ellipseConstraint: Constraint[] = [
  {type: 'eq', fn: (x) => 4 * x[0] ** 2 + x[1] ** 2 - 25},
];

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
/*  1. minimize z = x^2 + 12xy + 2y^2 subject to 4x^2 + y^2 = 25     */
/* ================================================================== */

const nm = new NelderMead();
const pso = new PSO();

console.log('========== MINIMIZE: z = x^2 + 12xy + 2y^2, s.t. 4x^2 + y^2 = 25 ==========');

printProblem([
  'z(x,y) = x^2 + 12xy + 2y^2',
  'subject to: 4x^2 + y^2 = 25',
  'goal:     minimize',
  'expected: min = -50 at (2, -3) or (-2, 3)',
]);

const minPenalized = applyPenalty(quadraticMixed, ellipseConstraint, {mu: 1_000_000});

// near the ellipse, in the basin of (2, -3)
printResult(
  'Nelder-Mead → minimize (x0 = [1, -2])',
  nm.minimize(minPenalized, new Float64Array([1, -2]), {
    maxIterations: 10_000,
    tolerance: 1e-12,
  }),
);

// near the ellipse, in the basin of (-2, 3)
printResult(
  'Nelder-Mead → minimize (x0 = [-1, 2])',
  nm.minimize(minPenalized, new Float64Array([-1, 2]), {
    maxIterations: 10_000,
    tolerance: 1e-12,
  }),
);

// (0, 5) is on the ellipse: 0 + 25 = 25
printResult(
  'PSO → minimize (x0 = [0, 5])',
  pso.minimize(minPenalized, new Float64Array([0, 5]), {
    swarmSize: 80,
    maxIterations: 5_000,
  }),
);

/* ================================================================== */
/*  2. maximize z = x^2 + 12xy + 2y^2 subject to 4x^2 + y^2 = 25     */
/* ================================================================== */

console.log('\n========== MAXIMIZE: z = x^2 + 12xy + 2y^2, s.t. 4x^2 + y^2 = 25 ==========');

printProblem([
  'z(x,y) = x^2 + 12xy + 2y^2',
  'subject to: 4x^2 + y^2 = 25',
  'goal:     maximize',
  'expected: max = 425/4 = 106.25 at (3/2, 4) or (-3/2, -4)',
]);

const maxSettings = {
  maxIterations: 10_000,
  tolerance: 1e-12,
  constraints: ellipseConstraint,
  penaltyOptions: {mu: 1_000_000},
};

// near the ellipse, in the basin of (1.5, 4)
printResult(
  'Nelder-Mead → maximize (x0 = [1, 2])',
  nm.maximize(quadraticMixed, new Float64Array([1, 2]), maxSettings),
);

// near the ellipse, in the basin of (-1.5, -4)
printResult(
  'Nelder-Mead → maximize (x0 = [-1, -2])',
  nm.maximize(quadraticMixed, new Float64Array([-1, -2]), maxSettings),
);

// (0, -5) is on the ellipse: 0 + 25 = 25
printResult(
  'PSO → maximize (x0 = [0, -5])',
  pso.maximize(quadraticMixed, new Float64Array([0, -5]), {
    ...maxSettings,
    swarmSize: 80,
    maxIterations: 5_000,
  }),
);
