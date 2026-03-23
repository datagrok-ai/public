import {NelderMead} from '../optimizers/nelder-mead';
import {PSO} from '../optimizers/pso';
import {GradientDescent} from '../optimizers/gradient-descent';
import {Adam} from '../optimizers/adam';
import type {ObjectiveFunction, OptimizationResult} from '../types';
import type {Optimizer} from '../optimizer';

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
/*  Test functions                                                     */
/* ================================================================== */

// Sphere: f(x) = Σ xᵢ²
// Dimension: n (tested: 2) | Type: unimodal, convex
// Search domain: −∞ ≤ xᵢ ≤ ∞
// Global minimum: f(0, …, 0) = 0
const sphere: ObjectiveFunction = (x) => {
  let s = 0;
  for (let i = 0; i < x.length; i++) s += x[i] * x[i];
  return s;
};

// Rosenbrock: f(x) = Σᵢ₌₁ⁿ⁻¹ [100·(xᵢ₊₁ − xᵢ²)² + (1 − xᵢ)²]
// Dimension: n (tested: 2) | Type: unimodal, non-convex
// Search domain: −∞ ≤ xᵢ ≤ ∞
// Global minimum: f(1, …, 1) = 0
const rosenbrock: ObjectiveFunction = (x) => {
  let s = 0;
  for (let i = 0; i < x.length - 1; i++)
    s += 100 * (x[i + 1] - x[i] * x[i]) ** 2 + (1 - x[i]) ** 2;
  return s;
};

// Beale: f(x,y) = (1.5 − x + xy)² + (2.25 − x + xy²)² + (2.625 − x + xy³)²
// Dimension: 2 | Type: unimodal, non-convex
// Search domain: −4.5 ≤ x, y ≤ 4.5
// Global minimum: f(3, 0.5) = 0
const beale: ObjectiveFunction = (x) =>
  (1.5 - x[0] + x[0] * x[1]) ** 2 +
  (2.25 - x[0] + x[0] * x[1] * x[1]) ** 2 +
  (2.625 - x[0] + x[0] * x[1] * x[1] * x[1]) ** 2;

// Booth: f(x,y) = (x + 2y − 7)² + (2x + y − 5)²
// Dimension: 2 | Type: unimodal, convex
// Search domain: −10 ≤ x, y ≤ 10
// Global minimum: f(1, 3) = 0
const booth: ObjectiveFunction = (x) =>
  (x[0] + 2 * x[1] - 7) ** 2 + (2 * x[0] + x[1] - 5) ** 2;

// Matyas: f(x,y) = 0.26·(x² + y²) − 0.48·x·y
// Dimension: 2 | Type: unimodal, convex
// Search domain: −10 ≤ x, y ≤ 10
// Global minimum: f(0, 0) = 0
const matyas: ObjectiveFunction = (x) =>
  0.26 * (x[0] * x[0] + x[1] * x[1]) - 0.48 * x[0] * x[1];

// Himmelblau: f(x,y) = (x² + y − 11)² + (x + y² − 7)²
// Dimension: 2 | Type: multimodal (4 identical minima)
// Search domain: −5 ≤ x, y ≤ 5
// Global minima: f(3,2) = f(−2.805118,3.131312) = f(−3.779310,−3.283186) = f(3.584428,−1.848126) = 0
const himmelblau: ObjectiveFunction = (x) =>
  (x[0] * x[0] + x[1] - 11) ** 2 + (x[0] + x[1] * x[1] - 7) ** 2;

// Three-Hump Camel: f(x,y) = 2x² − 1.05x⁴ + x⁶/6 + xy + y²
// Dimension: 2 | Type: multimodal (3 local minima)
// Search domain: −5 ≤ x, y ≤ 5
// Global minimum: f(0, 0) = 0
const threeHumpCamel: ObjectiveFunction = (x) =>
  2 * x[0] ** 2 - 1.05 * x[0] ** 4 + x[0] ** 6 / 6 + x[0] * x[1] + x[1] ** 2;

// Rastrigin: f(x) = 10n + Σ [xᵢ² − 10·cos(2πxᵢ)]
// Dimension: n (tested: 2) | Type: multimodal (highly)
// Search domain: −5.12 ≤ xᵢ ≤ 5.12
// Global minimum: f(0, …, 0) = 0
const rastrigin: ObjectiveFunction = (x) => {
  const A = 10;
  let s = A * x.length;
  for (let i = 0; i < x.length; i++)
    s += x[i] * x[i] - A * Math.cos(2 * Math.PI * x[i]);
  return s;
};

// Ackley: f(x,y) = −20·exp[−0.2·√(0.5·(x²+y²))] − exp[0.5·(cos(2πx)+cos(2πy))] + e + 20
// Dimension: 2 | Type: multimodal
// Search domain: −5 ≤ x, y ≤ 5
// Global minimum: f(0, 0) = 0
const ackley: ObjectiveFunction = (x) =>
  -20 * Math.exp(-0.2 * Math.sqrt(0.5 * (x[0] ** 2 + x[1] ** 2))) -
  Math.exp(0.5 * (Math.cos(2 * Math.PI * x[0]) + Math.cos(2 * Math.PI * x[1]))) +
  Math.E + 20;

// Lévi N.13: f(x,y) = sin²(3πx) + (x−1)²·(1+sin²(3πy)) + (y−1)²·(1+sin²(2πy))
// Dimension: 2 | Type: multimodal
// Search domain: −10 ≤ x, y ≤ 10
// Global minimum: f(1, 1) = 0
const levi13: ObjectiveFunction = (x) =>
  Math.sin(3 * Math.PI * x[0]) ** 2 +
  (x[0] - 1) ** 2 * (1 + Math.sin(3 * Math.PI * x[1]) ** 2) +
  (x[1] - 1) ** 2 * (1 + Math.sin(2 * Math.PI * x[1]) ** 2);

// Griewank: f(x) = 1 + (1/4000)·Σ xᵢ² − Π cos(xᵢ/√i)
// Dimension: n (tested: 2) | Type: multimodal
// Search domain: −600 ≤ xᵢ ≤ 600
// Global minimum: f(0, …, 0) = 0
const griewank: ObjectiveFunction = (x) => {
  let sum = 0;
  let prod = 1;
  for (let i = 0; i < x.length; i++) {
    sum += x[i] * x[i];
    prod *= Math.cos(x[i] / Math.sqrt(i + 1));
  }
  return 1 + sum / 4000 - prod;
};

// Styblinski-Tang: f(x) = (1/2)·Σ (xᵢ⁴ − 16xᵢ² + 5xᵢ)
// Dimension: n (tested: 2) | Type: multimodal
// Search domain: −5 ≤ xᵢ ≤ 5
// Global minimum: f(−2.903534, …, −2.903534) ≈ −39.16617·n
const styblinskiTang: ObjectiveFunction = (x) => {
  let s = 0;
  for (let i = 0; i < x.length; i++)
    s += x[i] ** 4 - 16 * x[i] ** 2 + 5 * x[i];
  return s / 2;
};

// Easom: f(x,y) = −cos(x)·cos(y)·exp(−((x−π)²+(y−π)²))
// Dimension: 2 | Type: unimodal (nearly flat everywhere except near (π,π))
// Search domain: −100 ≤ x, y ≤ 100
// Global minimum: f(π, π) = −1
const easom: ObjectiveFunction = (x) =>
  -Math.cos(x[0]) * Math.cos(x[1]) *
  Math.exp(-((x[0] - Math.PI) ** 2 + (x[1] - Math.PI) ** 2));

// Goldstein-Price: f(x,y) = [1+(x+y+1)²·(19−14x+3x²−14y+6xy+3y²)]·[30+(2x−3y)²·(18−32x+12x²+48y−36xy+27y²)]
// Dimension: 2 | Type: multimodal
// Search domain: −2 ≤ x, y ≤ 2
// Global minimum: f(0, −1) = 3
const goldsteinPrice: ObjectiveFunction = (x) => {
  const a = 1 + (x[0] + x[1] + 1) ** 2 *
    (19 - 14 * x[0] + 3 * x[0] ** 2 - 14 * x[1] + 6 * x[0] * x[1] + 3 * x[1] ** 2);
  const b = 30 + (2 * x[0] - 3 * x[1]) ** 2 *
    (18 - 32 * x[0] + 12 * x[0] ** 2 + 48 * x[1] - 36 * x[0] * x[1] + 27 * x[1] ** 2);
  return a * b;
};

// McCormick: f(x,y) = sin(x+y) + (x−y)² − 1.5x + 2.5y + 1
// Dimension: 2 | Type: multimodal
// Search domain: −1.5 ≤ x ≤ 4, −3 ≤ y ≤ 4
// Global minimum: f(−0.54719, −1.54719) ≈ −1.9133
const mccormick: ObjectiveFunction = (x) =>
  Math.sin(x[0] + x[1]) + (x[0] - x[1]) ** 2 - 1.5 * x[0] + 2.5 * x[1] + 1;

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
    settings: {maxIterations: 5_000},
  },
  {
    name: 'PSO',
    optimizer: new PSO(),
    settings: {maxIterations: 3_000, swarmSize: 30, seed: 42},
  },
  {
    name: 'GradientDescent',
    optimizer: new GradientDescent(),
    settings: {maxIterations: 10_000, learningRate: 0.01},
  },
  {
    name: 'Adam',
    optimizer: new Adam(),
    settings: {maxIterations: 10_000, learningRate: 0.001},
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

const HIMMELBLAU_MINIMA = [
  new Float64Array([3, 2]),
  new Float64Array([-2.805118, 3.131312]),
  new Float64Array([-3.779310, -3.283186]),
  new Float64Array([3.584428, -1.848126]),
];

function distToOptimum(problem: BenchmarkProblem, point: Float64Array): number {
  if (problem.name === 'Himmelblau')
    return Math.min(...HIMMELBLAU_MINIMA.map((m) => distance(point, m)));

  return distance(point, problem.knownPoint);
}

/* ================================================================== */
/*  Table printer                                                      */
/* ================================================================== */

interface RowData {
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
    for (const r of rows)
      if (r[k].length > widths[k]) widths[k] = r[k].length;
  }

  const header = keys.map((k) => pad(COL_HEADERS[k], widths[k])).join(' | ');
  const separator = keys.map((k) => '-'.repeat(widths[k])).join('-+-');
  console.log('  ' + header);
  console.log('  ' + separator);

  for (const r of rows) {
    const line = keys.map((k) =>
      RIGHT_ALIGNED.has(k) ? padLeft(r[k], widths[k]) : pad(r[k], widths[k]),
    ).join(' | ');
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
  console.log('║          Default hyperparameters · 4 optimizers · 15 problems    ║');
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

      rows.push({
        method: cfg.name,
        settings: fmtSettings(cfg.settings),
        foundValue: fmtNum(result.value),
        foundPoint: fmtPoint(result.point),
        error: fmtNum(Math.abs(result.value - problem.knownMin)),
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
