import {NelderMead} from '../optimizers/nelder-mead';
import {PSO} from '../optimizers/pso';
import {GradientDescent} from '../optimizers/gradient-descent';
import {Adam} from '../optimizers/adam';
import {LBFGS} from '../optimizers/lbfgs';
import {LBFGSB} from '../optimizers/lbfgs-b';
import type {ObjectiveFunction, OptimizationResult} from '../types';
import type {Optimizer} from '../optimizer';
import {
  sphere, rosenbrock, beale, booth, matyas, himmelblau, threeHumpCamel,
  rastrigin, ackley, levi13, griewank, styblinskiTang, easom, goldsteinPrice,
  mccormick, HIMMELBLAU_MINIMA,
} from './test-functions';

/* ================================================================== */
/*  Function evaluation counter                                        */
/* ================================================================== */

interface Counted {
  fn: ObjectiveFunction;
  count: () => number;
  reset: () => void;
}

function withCounter(fn: ObjectiveFunction): Counted {
  let calls = 0;
  return {
    fn: (x: Float64Array) => {calls++; return fn(x);},
    count: () => calls,
    reset: () => {calls = 0;},
  };
}

/* ================================================================== */
/*  Benchmark problem definition                                       */
/* ================================================================== */

interface BenchmarkProblem {
  name: string;
  dimension: number;
  type: string;
  fn: ObjectiveFunction;
  x0: Float64Array;
  knownMin: number;
  knownPoint: Float64Array;
}

/* ================================================================== */
/*  Problem suite                                                      */
/* ================================================================== */

const problems: BenchmarkProblem[] = [
  {
    name: 'Sphere',
    dimension: 2,
    type: 'Unimodal, convex',
    fn: sphere,
    x0: new Float64Array([5, -3]),
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Rosenbrock',
    dimension: 2,
    type: 'Unimodal, narrow valley',
    fn: rosenbrock,
    x0: new Float64Array([-1.2, 1.0]),
    knownMin: 0,
    knownPoint: new Float64Array([1, 1]),
  },
  {
    name: 'Beale',
    dimension: 2,
    type: 'Unimodal, non-convex',
    fn: beale,
    x0: new Float64Array([0, 0]),
    knownMin: 0,
    knownPoint: new Float64Array([3, 0.5]),
  },
  {
    name: 'Booth',
    dimension: 2,
    type: 'Unimodal, convex',
    fn: booth,
    x0: new Float64Array([0, 0]),
    knownMin: 0,
    knownPoint: new Float64Array([1, 3]),
  },
  {
    name: 'Matyas',
    dimension: 2,
    type: 'Unimodal, convex',
    fn: matyas,
    x0: new Float64Array([5, -5]),
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Himmelblau',
    dimension: 2,
    type: 'Multimodal (4 minima)',
    fn: himmelblau,
    x0: new Float64Array([0, 0]),
    knownMin: 0,
    knownPoint: new Float64Array([3, 2]),
  },
  {
    name: 'Three-Hump Camel',
    dimension: 2,
    type: 'Multimodal (3 local minima)',
    fn: threeHumpCamel,
    x0: new Float64Array([2, -1]),
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Rastrigin',
    dimension: 2,
    type: 'Multimodal (highly)',
    fn: rastrigin,
    x0: new Float64Array([2.5, -3.5]),
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Ackley',
    dimension: 2,
    type: 'Multimodal',
    fn: ackley,
    x0: new Float64Array([2, -2]),
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Lévi N.13',
    dimension: 2,
    type: 'Multimodal',
    fn: levi13,
    x0: new Float64Array([-4, 5]),
    knownMin: 0,
    knownPoint: new Float64Array([1, 1]),
  },
  {
    name: 'Griewank',
    dimension: 2,
    type: 'Multimodal',
    fn: griewank,
    x0: new Float64Array([100, -200]),
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Styblinski-Tang',
    dimension: 2,
    type: 'Multimodal',
    fn: styblinskiTang,
    x0: new Float64Array([0, 0]),
    knownMin: -39.16617 * 2,
    knownPoint: new Float64Array([-2.903534, -2.903534]),
  },
  {
    name: 'Easom',
    dimension: 2,
    type: 'Nearly flat, narrow peak',
    fn: easom,
    x0: new Float64Array([1, 1]),
    knownMin: -1,
    knownPoint: new Float64Array([Math.PI, Math.PI]),
  },
  {
    name: 'Goldstein-Price',
    dimension: 2,
    type: 'Multimodal',
    fn: goldsteinPrice,
    x0: new Float64Array([0, 0]),
    knownMin: 3,
    knownPoint: new Float64Array([0, -1]),
  },
  {
    name: 'McCormick',
    dimension: 2,
    type: 'Multimodal',
    fn: mccormick,
    x0: new Float64Array([0, 0]),
    knownMin: -1.9133,
    knownPoint: new Float64Array([-0.54719, -1.54719]),
  },
];

/* ================================================================== */
/*  Optimizer configurations (default hyperparameters)                 */
/* ================================================================== */

interface OptimizerConfig {
  name: string;
  optimizer: Optimizer<any>;
  settings: Record<string, unknown>;
}

const optimizers: OptimizerConfig[] = [
  {
    name: 'Nelder-Mead',
    optimizer: new NelderMead(),
    settings: {maxIterations: 10_000},
  },
  {
    name: 'PSO',
    optimizer: new PSO(),
    settings: {maxIterations: 10_000, swarmSize: 50, seed: 42},
  },
  {
    name: 'GradientDescent',
    optimizer: new GradientDescent(),
    settings: {maxIterations: 10_000, learningRate: 0.001, momentum: 0.9},
  },
  {
    name: 'Adam',
    optimizer: new Adam(),
    settings: {maxIterations: 10_000, learningRate: 0.1},
  },
  {
    name: 'L-BFGS',
    optimizer: new LBFGS(),
    settings: {maxIterations: 1_000, historySize: 10},
  },
  {
    name: 'L-BFGS-B',
    optimizer: new LBFGSB(),
    settings: {maxIterations: 1_000, historySize: 10},
  },
];

/* ================================================================== */
/*  Formatting helpers                                                 */
/* ================================================================== */

function pad(s: string, len: number): string {
  return s.length >= len ? s : s + ' '.repeat(len - s.length);
}

function padLeft(s: string, len: number): string {
  return s.length >= len ? s : ' '.repeat(len - s.length) + s;
}

function fmtNum(v: number, decimals: number = 6): string {
  if (Math.abs(v) < 1e-3 && v !== 0) return v.toExponential(decimals);
  return v.toFixed(decimals);
}

function fmtPoint(p: Float64Array): string {
  return '(' + Array.from(p).map((v) => fmtNum(v, 4)).join(', ') + ')';
}

function fmtSettings(s: Record<string, unknown>): string {
  return Object.entries(s)
    .filter(([k]) => k !== 'onIteration')
    .map(([k, v]) => `${k}=${v}`)
    .join(', ');
}

function distance(a: Float64Array, b: Float64Array): number {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += (a[i] - b[i]) ** 2;
  return Math.sqrt(s);
}

/* ================================================================== */
/*  Himmelblau: distance to nearest of 4 known minima                  */
/* ================================================================== */

function distToOptimum(problem: BenchmarkProblem, point: Float64Array): number {
  if (problem.name === 'Himmelblau')
    return Math.min(...HIMMELBLAU_MINIMA.map((m) => distance(point, m)));

  return distance(point, problem.knownPoint);
}

/* ================================================================== */
/*  Table printer                                                      */
/* ================================================================== */

interface RowData {
  icon: string;
  method: string;
  settings: string;
  foundValue: string;
  foundPoint: string;
  error: string;
  distToOpt: string;
  converged: string;
  iterations: string;
  fnEvals: string;
  wallTimeMs: string;
}

const COL_HEADERS: Record<keyof RowData, string> = {
  icon: '',
  method: 'Method',
  settings: 'Settings',
  foundValue: 'Found Value',
  foundPoint: 'Found Point',
  error: '|Error|',
  distToOpt: 'Dist to Opt',
  converged: 'Conv',
  iterations: 'Iters',
  fnEvals: 'Fn Evals',
  wallTimeMs: 'Time (ms)',
};

const RIGHT_ALIGNED: Set<keyof RowData> = new Set([
  'foundValue', 'error', 'distToOpt', 'iterations', 'fnEvals', 'wallTimeMs',
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

/* ================================================================== */
/*  Run benchmarks                                                     */
/* ================================================================== */

function runBenchmarks(): void {
  console.log('');
  console.log('╔══════════════════════════════════════════════════════════════════╗');
  console.log('║          UNCONSTRAINED OPTIMIZATION BENCHMARKS                   ║');
  console.log('║          Default hyperparameters · 5 optimizers · 15 problems    ║');
  console.log('╚══════════════════════════════════════════════════════════════════╝');
  console.log('');

  for (const problem of problems) {
    console.log('━'.repeat(78));
    console.log(`  ${problem.name}`);
    console.log(`  Dim: ${problem.dimension} | Type: ${problem.type}`);
    console.log(`  Known min: f${fmtPoint(problem.knownPoint)} = ${fmtNum(problem.knownMin)}`);
    console.log(`  x0 = ${fmtPoint(problem.x0)}`);
    console.log('━'.repeat(78));

    const rows: RowData[] = [];

    for (const cfg of optimizers) {
      const counted = withCounter(problem.fn);

      const t0 = performance.now();
      let result: OptimizationResult;

      try {
        result = cfg.optimizer.minimize(
          counted.fn,
          Float64Array.from(problem.x0),
          cfg.settings as any,
        );
      } catch (e: any) {
        rows.push({
          icon: '\u274C',
          method: cfg.name,
          settings: fmtSettings(cfg.settings),
          foundValue: 'ERROR',
          foundPoint: '-',
          error: (e.message ?? 'unknown').slice(0, 40),
          distToOpt: '-',
          converged: '-',
          iterations: '-',
          fnEvals: '-',
          wallTimeMs: '-',
        });
        continue;
      }

      const elapsed = performance.now() - t0;
      const distOpt = distToOptimum(problem, result.point);

      const absError = Math.abs(result.value - problem.knownMin);
      const quality = Math.max(absError, distOpt);
      let icon: string;
      if (!result.converged)
        icon = '\u274C';
      else if (quality < 1e-3)
        icon = '\u2705';
      else
        icon = '\u26A0';

      rows.push({
        icon,
        method: cfg.name,
        settings: fmtSettings(cfg.settings),
        foundValue: fmtNum(result.value),
        foundPoint: fmtPoint(result.point),
        error: fmtNum(absError),
        distToOpt: fmtNum(distOpt),
        converged: result.converged ? 'yes' : 'no',
        iterations: String(result.iterations),
        fnEvals: String(counted.count()),
        wallTimeMs: elapsed.toFixed(1),
      });
    }

    printTable(rows);
    console.log('');
  }

  // Summary
  console.log('━'.repeat(78));
  console.log('  Notes:');
  console.log('  · All optimizers use default hyperparameters (see Settings column)');
  console.log('  · Fn Evals = total calls to objective function (via counter wrapper)');
  console.log('  · |Error| = |found_value − known_minimum|');
  console.log('  · Dist to Opt = ‖found_point − known_point‖₂');
  console.log('  · PSO uses seed=42 for reproducibility');
  console.log('  · Himmelblau: Dist to Opt uses nearest of 4 known minima');
  console.log('━'.repeat(78));
}

runBenchmarks();
