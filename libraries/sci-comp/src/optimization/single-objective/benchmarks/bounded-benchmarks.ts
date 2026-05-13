/**
 * Bounded-optimization benchmarks.
 *
 * Runs every optimizer in the library on a suite of box-constrained
 * problems, comparing:
 *   - **L-BFGS-B** — native box constraints via `settings.bounds`
 *     (geometric handling through Cauchy point + subspace minimisation),
 *   - **others** (L-BFGS, GD, Adam, PSO, Nelder-Mead) — bounds translated
 *     to inequality constraints via `boxConstraints()` and handled by
 *     the quadratic-penalty layer.
 *
 * A solution counts as "correct" only when it is (i) within tolerance of
 * the known minimum AND (ii) exactly feasible. A penalty-based solver
 * that violates a bound, even slightly, fails the feasibility check —
 * this is the main motivation for native bounds.
 *
 * The problem list (`BOUNDED_PROBLEMS`) is sorted by
 * `activeBoundsAtOptimum` ascending so that "bounds inactive" sanity
 * checks surface first; any regression there indicates L-BFGS-B is
 * mishandling a scenario that is effectively unconstrained.
 */
import {NelderMead} from '../optimizers/nelder-mead';
import {PSO} from '../optimizers/pso';
import {GradientDescent} from '../optimizers/gradient-descent';
import {Adam} from '../optimizers/adam';
import {LBFGS} from '../optimizers/lbfgs';
import {LBFGSB} from '../optimizers/lbfgs-b';
import {boxConstraints} from '../penalty';
import type {ObjectiveFunction, OptimizationResult} from '../types';
import type {Optimizer} from '../optimizer';
import {BOUNDED_PROBLEMS} from './test-functions';
import type {BoundedTestProblem} from './test-functions';

/* ================================================================== */
/*  Counter                                                            */
/* ================================================================== */

interface Counted {
  fn: ObjectiveFunction;
  count: () => number;
}

function withCounter(fn: ObjectiveFunction): Counted {
  let calls = 0;
  return {
    fn: (x: Float64Array) => {calls++; return fn(x);},
    count: () => calls,
  };
}

/* ================================================================== */
/*  Optimizer configurations                                           */
/* ================================================================== */

type BoundsMode = 'native' | 'penalty';

interface OptimizerConfig {
  name: string;
  optimizer: Optimizer<any>;
  boundsMode: BoundsMode;
  settings: Record<string, unknown>;
}

const optimizers: OptimizerConfig[] = [
  {
    name: 'L-BFGS-B',
    optimizer: new LBFGSB(),
    boundsMode: 'native',
    settings: {maxIterations: 1_000, historySize: 10},
  },
  {
    name: 'L-BFGS',
    optimizer: new LBFGS(),
    boundsMode: 'penalty',
    settings: {maxIterations: 1_000, historySize: 10},
  },
  {
    name: 'Adam',
    optimizer: new Adam(),
    boundsMode: 'penalty',
    settings: {maxIterations: 10_000, learningRate: 0.1},
  },
  {
    name: 'GradientDescent',
    optimizer: new GradientDescent(),
    boundsMode: 'penalty',
    settings: {maxIterations: 10_000, learningRate: 0.001, momentum: 0.9},
  },
  {
    name: 'PSO',
    optimizer: new PSO(),
    boundsMode: 'penalty',
    settings: {maxIterations: 10_000, swarmSize: 50, seed: 42},
  },
  {
    name: 'Nelder-Mead',
    optimizer: new NelderMead(),
    boundsMode: 'penalty',
    settings: {maxIterations: 10_000},
  },
];

/* ================================================================== */
/*  Settings construction per problem                                  */
/* ================================================================== */

function buildSettings(
  cfg: OptimizerConfig,
  problem: BoundedTestProblem,
): Record<string, unknown> {
  if (cfg.boundsMode === 'native') {
    return {
      ...cfg.settings,
      bounds: {
        lower: Array.from(problem.lower),
        upper: Array.from(problem.upper),
      },
    };
  }
  return {
    ...cfg.settings,
    constraints: boxConstraints(problem.lower, problem.upper),
  };
}

/* ================================================================== */
/*  Feasibility check                                                  */
/* ================================================================== */

function feasibilityViolation(
  point: Float64Array,
  lower: Float64Array,
  upper: Float64Array,
): number {
  let worst = 0;
  for (let i = 0; i < point.length; i++) {
    if (point[i] < lower[i]) {
      const v = lower[i] - point[i];
      if (v > worst) worst = v;
    }
    if (point[i] > upper[i]) {
      const v = point[i] - upper[i];
      if (v > worst) worst = v;
    }
  }
  return worst;
}

/* ================================================================== */
/*  Formatting                                                         */
/* ================================================================== */

function pad(s: string, len: number): string {
  return s.length >= len ? s : s + ' '.repeat(len - s.length);
}

function padLeft(s: string, len: number): string {
  return s.length >= len ? s : ' '.repeat(len - s.length) + s;
}

function fmtNum(v: number, decimals = 4): string {
  if (!Number.isFinite(v)) return String(v);
  if (v === 0) return '0';
  if (Math.abs(v) < 1e-3) return v.toExponential(2);
  if (Math.abs(v) >= 1e4) return v.toExponential(2);
  return v.toFixed(decimals);
}

function distance(a: Float64Array, b: Float64Array): number {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += (a[i] - b[i]) ** 2;
  return Math.sqrt(s);
}

/* ================================================================== */
/*  Runner                                                             */
/* ================================================================== */

interface RowData {
  icon: string;
  method: string;
  mode: string;
  foundValue: string;
  error: string;
  distToOpt: string;
  feasVio: string;
  converged: string;
  iterations: string;
  fnEvals: string;
  wallTimeMs: string;
}

const COL_HEADERS: Record<keyof RowData, string> = {
  icon: '',
  method: 'Method',
  mode: 'Mode',
  foundValue: 'Found Value',
  error: '|Error|',
  distToOpt: 'Dist to Opt',
  feasVio: 'Feas Vio',
  converged: 'Conv',
  iterations: 'Iters',
  fnEvals: 'Fn Evals',
  wallTimeMs: 'Time (ms)',
};

const RIGHT_ALIGNED: Set<keyof RowData> = new Set([
  'foundValue', 'error', 'distToOpt', 'feasVio', 'iterations', 'fnEvals', 'wallTimeMs',
]);

function printTable(rows: RowData[]): void {
  const keys = Object.keys(COL_HEADERS) as (keyof RowData)[];
  const widths: Record<string, number> = {};

  for (const k of keys) {
    widths[k] = COL_HEADERS[k].length;
    for (const r of rows) {
      const len = k === 'icon' ? 2 : r[k].length;
      if (len > widths[k]) widths[k] = len;
    }
  }

  const header = keys.map((k) => pad(COL_HEADERS[k], widths[k])).join(' | ');
  const separator = keys.map((k) => '-'.repeat(widths[k])).join('-+-');
  console.log('  ' + header);
  console.log('  ' + separator);

  for (const r of rows) {
    const line = keys.map((k) => {
      if (k === 'icon') return r[k] + ' ';
      return RIGHT_ALIGNED.has(k) ? padLeft(r[k], widths[k]) : pad(r[k], widths[k]);
    }).join(' | ');
    console.log('  ' + line);
  }
}

function runBenchmarks(): void {
  console.log('');
  console.log('╔══════════════════════════════════════════════════════════════════╗');
  console.log('║              BOUNDED OPTIMIZATION BENCHMARKS                     ║');
  console.log('║              L-BFGS-B native vs others via penalty               ║');
  console.log('╚══════════════════════════════════════════════════════════════════╝');
  console.log('');

  const successCounts: Record<string, number> = {};
  for (const cfg of optimizers) successCounts[cfg.name] = 0;

  // Sort by activeBoundsAtOptimum (ascending): inactive-bound regressions first.
  const sortedProblems = [...BOUNDED_PROBLEMS].sort(
    (a, b) => a.activeBoundsAtOptimum - b.activeBoundsAtOptimum,
  );

  for (const problem of sortedProblems) {
    console.log('━'.repeat(90));
    console.log(`  ${problem.name}`);
    console.log(`  Dim: ${problem.dimension} | Type: ${problem.type}`);
    console.log(`  Active bounds at optimum: ${problem.activeBoundsAtOptimum}/${problem.dimension}`);
    console.log(`  Known min: f = ${fmtNum(problem.knownMin)} at ` +
      '(' + Array.from(problem.knownPoint).map((v) => fmtNum(v, 3)).join(', ') + ')');
    console.log('━'.repeat(90));

    const rows: RowData[] = [];

    for (const cfg of optimizers) {
      const counted = withCounter(problem.fn);
      const settings = buildSettings(cfg, problem);

      const t0 = performance.now();
      let result: OptimizationResult;
      try {
        result = cfg.optimizer.minimize(
          counted.fn,
          Float64Array.from(problem.x0),
          settings as any,
        );
      } catch (e: any) {
        rows.push({
          icon: '❌',
          method: cfg.name,
          mode: cfg.boundsMode,
          foundValue: 'ERROR',
          error: '-',
          distToOpt: '-',
          feasVio: '-',
          converged: '-',
          iterations: '-',
          fnEvals: '-',
          wallTimeMs: '-',
        });
        continue;
      }
      const elapsed = performance.now() - t0;

      const absError = Math.abs(result.value - problem.knownMin);
      const distOpt = distance(result.point, problem.knownPoint);
      const feasVio = feasibilityViolation(result.point, problem.lower, problem.upper);

      const success =
        result.converged &&
        feasVio < 1e-6 &&
        (absError < 1e-3 || distOpt < 1e-3);

      if (success) successCounts[cfg.name]++;

      let icon: string;
      if (feasVio >= 1e-6) icon = '❌'; // infeasible — fail outright
      else if (success) icon = '✅';
      else if (result.converged) icon = '⚠';
      else icon = '❌';

      rows.push({
        icon,
        method: cfg.name,
        mode: cfg.boundsMode,
        foundValue: fmtNum(result.value),
        error: fmtNum(absError),
        distToOpt: fmtNum(distOpt),
        feasVio: feasVio === 0 ? '0' : fmtNum(feasVio),
        converged: result.converged ? 'yes' : 'no',
        iterations: String(result.iterations),
        fnEvals: String(counted.count()),
        wallTimeMs: elapsed.toFixed(1),
      });
    }

    printTable(rows);
    console.log('');
  }

  console.log('━'.repeat(90));
  console.log('  Success summary (feasible + accurate):');
  console.log('━'.repeat(90));
  const total = sortedProblems.length;
  const rankings = Object.entries(successCounts)
    .sort((a, b) => b[1] - a[1]);
  for (const [name, count] of rankings)
    console.log(`  ${pad(name, 20)}  ${count}/${total}`);
  console.log('');
  console.log('  Notes:');
  console.log('  · Success = (converged) ∧ (feas_vio < 1e-6) ∧ (|Error| < 1e-3 ∨ Dist < 1e-3)');
  console.log('  · "Mode": native = L-BFGS-B via settings.bounds (geometric);');
  console.log('            penalty = other optimizers via boxConstraints (quadratic).');
  console.log('  · Feas Vio = max overshoot past the bound in any coordinate.');
  console.log('  · ❌ iff infeasible OR did not converge. ⚠ = converged but wrong answer.');
}

runBenchmarks();
