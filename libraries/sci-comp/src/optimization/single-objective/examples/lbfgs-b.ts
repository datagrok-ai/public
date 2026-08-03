/**
 * L-BFGS-B — native box constraints, analytic gradient, fixed variables.
 *
 * Run with:
 *   npx tsx src/optimization/single-objective/examples/lbfgs-b.ts
 */
import {LBFGSB} from '..';
import type {OptimizationResult} from '..';

/* ================================================================== */
/*  Test functions                                                     */
/* ================================================================== */

const rosenbrock = (x: Float64Array): number => {
  let sum = 0;
  for (let i = 0; i < x.length - 1; i++)
    sum += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2;
  return sum;
};

const rosenbrockGrad = (x: Float64Array, gOut: Float64Array): void => {
  const n = x.length;
  for (let i = 0; i < n; i++) gOut[i] = 0;
  for (let i = 0; i < n - 1; i++) {
    const t = x[i + 1] - x[i] * x[i];
    gOut[i] += -400 * x[i] * t - 2 * (1 - x[i]);
    gOut[i + 1] += 200 * t;
  }
};

/* ================================================================== */
/*  Reporter                                                           */
/* ================================================================== */

function printResult(label: string, r: OptimizationResult, fEvals?: number): void {
  console.log(`\n--- ${label} ---`);
  console.log(`  value:      ${r.value.toFixed(10)}`);
  console.log(`  point:      [${[...r.point].map((v) => v.toFixed(6)).join(', ')}]`);
  console.log(`  iterations: ${r.iterations}`);
  if (fEvals !== undefined) console.log(`  f evals:    ${fEvals}`);
  console.log(`  converged:  ${r.converged}`);
}

function withCounter(fn: (x: Float64Array) => number) {
  let n = 0;
  const wrapped = (x: Float64Array): number => {n++; return fn(x);};
  return {fn: wrapped, count: () => n};
}

/* ================================================================== */
/*  Example 1 — unconstrained Rosenbrock with analytic gradient        */
/* ================================================================== */

function example1(): void {
  console.log('\n=============================================================');
  console.log(' Example 1: Unconstrained Rosenbrock 10D, analytic gradient');
  console.log('=============================================================');
  const n = 10;
  const x0 = new Float64Array(n).fill(-1.2);
  const opt = new LBFGSB();
  const c = withCounter(rosenbrock);
  const r = opt.minimize(c.fn, x0, {
    gradFn: rosenbrockGrad,
    maxIterations: 500,
    gradTolerance: 1e-8,
  });
  printResult('Analytic gradient', r, c.count());
  console.log('  Expected: x ≈ (1, 1, …, 1), f = 0.');
}

/* ================================================================== */
/*  Example 2 — analytic vs finite-difference gradient                 */
/* ================================================================== */

function example2(): void {
  console.log('\n=============================================================');
  console.log(' Example 2: Analytic gradient vs FD fallback (Rosenbrock 2D)');
  console.log('=============================================================');
  const x0 = new Float64Array([-1.2, 1.0]);
  const opt = new LBFGSB();
  const settings = {maxIterations: 200, gradTolerance: 1e-10};

  const c1 = withCounter(rosenbrock);
  const r1 = opt.minimize(c1.fn, x0, {...settings, gradFn: rosenbrockGrad});
  printResult('Analytic gradient', r1, c1.count());

  const c2 = withCounter(rosenbrock);
  const r2 = opt.minimize(c2.fn, x0, settings);
  printResult('Finite-difference fallback', r2, c2.count());
  console.log('  FD costs 2n extra fn calls per gradient → many more evals.');
}

/* ================================================================== */
/*  Example 3 — bounded Rosenbrock                                     */
/* ================================================================== */

function example3(): void {
  console.log('\n=============================================================');
  console.log(' Example 3: Rosenbrock 2D on [-1.5, 0.5]² (upper bound active)');
  console.log('=============================================================');
  const x0 = new Float64Array([0.0, 0.0]);
  const opt = new LBFGSB();
  const r = opt.minimize(rosenbrock, x0, {
    bounds: {lower: -1.5, upper: 0.5},
    gradFn: rosenbrockGrad,
    maxIterations: 300,
    gradTolerance: 1e-8,
  });
  printResult('Bounded Rosenbrock', r);
  console.log('  Unconstrained min at (1, 1) is infeasible.');
  console.log('  Bounded min: x ≈ (0.5, 0.25), f ≈ 0.25 with x[0] pinned at upper.');
}

/* ================================================================== */
/*  Example 4 — fixed variables (l == u) and half-bounded              */
/* ================================================================== */

function example4(): void {
  console.log('\n=============================================================');
  console.log(' Example 4: Mixed bounds — fixed, half-bounded, unbounded');
  console.log('=============================================================');
  // f(x) = Σ xᵢ². Coord 1 fixed at 2 (l == u). Coord 2 has lower bound 0.
  // Coord 0 unbounded. Coord 3 bounded both sides. Coord 4 unbounded.
  const n = 5;
  const x0 = new Float64Array([5, 99, -1, 3, -2]);
  const sphereFn = (x: Float64Array): number => {
    let s = 0;
    for (let i = 0; i < n; i++) s += x[i] * x[i];
    return s;
  };
  const sphereGrad = (x: Float64Array, gOut: Float64Array): void => {
    for (let i = 0; i < n; i++) gOut[i] = 2 * x[i];
  };
  const opt = new LBFGSB();
  const r = opt.minimize(sphereFn, x0, {
    bounds: {
      lower: [-Infinity, 2, 0, -1, -Infinity],
      upper: [Infinity, 2, Infinity, 1, Infinity],
    },
    gradFn: sphereGrad,
    maxIterations: 100,
    gradTolerance: 1e-10,
  });
  printResult('Mixed bounds', r);
  console.log('  Coord 0: unbounded → 0.');
  console.log('  Coord 1: fixed at 2 → stays 2 regardless of x0.');
  console.log('  Coord 2: lower=0 active → pinned to 0.');
  console.log('  Coord 3: interior → 0 (within [-1, 1]).');
  console.log('  Coord 4: unbounded → 0.');
  console.log('  Expected f = 4 (from the fixed coord 1).');
}

/* ================================================================== */

console.log('L-BFGS-B examples');
console.log('──────────────────');
example1();
example2();
example3();
example4();
console.log('\nDone.\n');
