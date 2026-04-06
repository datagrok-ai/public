import {NelderMead, PSO} from '..';
import type {OptimizationResult, IterationState} from '..';

/* ================================================================== */
/*  Test functions                                                     */
/* ================================================================== */

/** Rosenbrock: min = 0 at (1, 1) */
const rosenbrock = (x: Float64Array): number => {
  let sum = 0;
  for (let i = 0; i < x.length - 1; i++)
    sum += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2;
  return sum;
};

/** Sphere: min = 0 at origin */
const sphere = (x: Float64Array): number => {
  let sum = 0;
  for (let i = 0; i < x.length; i++) sum += x[i] ** 2;
  return sum;
};

/* ================================================================== */
/*  Helper                                                             */
/* ================================================================== */

function printResult(label: string, r: OptimizationResult) {
  console.log(`\n--- ${label} ---`);
  console.log(`  value:      ${r.value.toFixed(6)}`);
  console.log(`  point:      [${[...r.point].map((v) => v.toFixed(6)).join(', ')}]`);
  console.log(`  iterations: ${r.iterations}`);
  console.log(`  converged:  ${r.converged}`);
  if (!r.converged) console.warn('  WARNING: optimizer did not converge');
}

/* ================================================================== */
/*  Main (wrapped in async IIFE for await support)                     */
/* ================================================================== */

(async () => {
  /* ================================================================ */
  /*  1. Async objective function                                      */
  /* ================================================================ */

  console.log('========== ASYNC: minimizeAsync ==========');

  console.log('  Simulating an async objective (e.g. remote evaluation)');

  /** Wraps a sync function as async — simulates an API call or simulation. */
  const asyncRosenbrock = async (x: Float64Array): Promise<number> => rosenbrock(x);

  const nm = new NelderMead();

  printResult(
    'Nelder-Mead → async Rosenbrock 2D',
    await nm.minimizeAsync(asyncRosenbrock, new Float64Array([-1.2, 1.0]), {
      maxIterations: 5_000,
      tolerance: 1e-12,
    }),
  );

  const pso = new PSO();

  printResult(
    'PSO → async Rosenbrock 2D',
    await pso.minimizeAsync(asyncRosenbrock, new Float64Array([-1.2, 1.0]), {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
    }),
  );

  /* ================================================================ */
  /*  2. onIteration callback — progress monitoring                    */
  /* ================================================================ */

  console.log('\n========== CALLBACK: onIteration (progress) ==========');

  console.log('  Logging every 100th iteration:');

  printResult(
    'Nelder-Mead → Sphere 3D with progress',
    nm.minimize(sphere, new Float64Array([5, -3, 7]), {
      maxIterations: 5_000,
      onIteration: (state: IterationState) => {
        if (state.iteration % 100 === 0)
          console.log(`    iter ${state.iteration}: best = ${state.bestValue.toFixed(8)}`);
      },
    }),
  );

  /* ================================================================ */
  /*  3. onIteration callback — early stopping                         */
  /* ================================================================ */

  console.log('\n========== CALLBACK: onIteration (early stopping) ==========');

  console.log('  Stop as soon as value < 0.001:');

  printResult(
    'PSO → Sphere 3D with early stop',
    pso.minimize(sphere, new Float64Array([5, -3, 7]), {
      swarmSize: 40,
      maxIterations: 3_000,
      seed: 42,
      onIteration: (state: IterationState) => {
        if (state.bestValue < 0.001) {
          console.log(`    early stop at iter ${state.iteration}, value = ${state.bestValue.toFixed(8)}`);
          return true; // stop
        }
      },
    }),
  );
})();
