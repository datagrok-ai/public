/**
 * Multi-start unconstrained benchmarks.
 *
 * Complements `unconstrained-benchmarks.ts` (one x₀ per problem) by running
 * every optimizer from **three** carefully-chosen starting points per problem,
 * specifically picked to expose failure modes:
 *
 *   - A (baseline) — the x₀ used in the single-start benchmark
 *   - B (perturbed / adversarial) — a small non-integer shift of A that
 *     breaks algorithm-specific luck (e.g. zero sin() gradient contributions
 *     on Rastrigin / Lévi with integer- or half-integer x₀)
 *   - C (near-optimum) — a probe close to the known global minimum, so we
 *     can distinguish "no chance from this x₀" vs "algorithm can't finish"
 *
 * The output highlights x₀-sensitivity rather than absolute performance.
 */
import {NelderMead} from '../optimizers/nelder-mead';
import {PSO} from '../optimizers/pso';
import {GradientDescent} from '../optimizers/gradient-descent';
import {Adam} from '../optimizers/adam';
import {LBFGS} from '../optimizers/lbfgs';
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
}

function withCounter(fn: ObjectiveFunction): Counted {
  let calls = 0;
  return {
    fn: (x: Float64Array) => {calls++; return fn(x);},
    count: () => calls,
  };
}

/* ================================================================== */
/*  Problem + x₀ triple                                                */
/* ================================================================== */

interface MultistartProblem {
  name: string;
  type: string;
  fn: ObjectiveFunction;
  x0s: [Float64Array, Float64Array, Float64Array];
  /** Short per-x₀ commentary shown in the markdown. */
  notes: [string, string, string];
  knownMin: number;
  knownPoint: Float64Array;
}

const problems: MultistartProblem[] = [
  {
    name: 'Sphere',
    type: 'Unimodal, convex',
    fn: sphere,
    x0s: [
      new Float64Array([5, -3]),
      new Float64Array([-10, 10]),
      new Float64Array([100, 100]),
    ],
    notes: ['baseline', 'mirrored far', 'very far (numerical stress)'],
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Rosenbrock',
    type: 'Unimodal, narrow valley',
    fn: rosenbrock,
    x0s: [
      new Float64Array([-1.2, 1.0]),
      new Float64Array([0, 0]),
      new Float64Array([2, -2]),
    ],
    notes: ['classic', 'origin', 'past optimum'],
    knownMin: 0,
    knownPoint: new Float64Array([1, 1]),
  },
  {
    name: 'Beale',
    type: 'Unimodal, non-convex',
    fn: beale,
    x0s: [
      new Float64Array([0, 0]),
      new Float64Array([-2, -2]),
      new Float64Array([4, 0.4]),
    ],
    notes: ['baseline', 'opposite quadrant', 'near-optimum'],
    knownMin: 0,
    knownPoint: new Float64Array([3, 0.5]),
  },
  {
    name: 'Booth',
    type: 'Unimodal, convex',
    fn: booth,
    x0s: [
      new Float64Array([0, 0]),
      new Float64Array([-5, -5]),
      new Float64Array([10, 10]),
    ],
    notes: ['baseline', 'opposite quadrant', 'far'],
    knownMin: 0,
    knownPoint: new Float64Array([1, 3]),
  },
  {
    name: 'Matyas',
    type: 'Unimodal, convex',
    fn: matyas,
    x0s: [
      new Float64Array([5, -5]),
      new Float64Array([-8, 8]),
      new Float64Array([0.1, 0.1]),
    ],
    notes: ['baseline', 'mirrored', 'near-optimum'],
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Himmelblau',
    type: 'Multimodal (4 minima)',
    fn: himmelblau,
    x0s: [
      new Float64Array([0, 0]),
      new Float64Array([-4, -4]),
      new Float64Array([3, -3]),
    ],
    notes: [
      'toward basin of (3, 2)',
      'toward basin of (−3.78, −3.28)',
      'toward basin of (3.58, −1.85)',
    ],
    knownMin: 0,
    knownPoint: new Float64Array([3, 2]),
  },
  {
    name: 'Three-Hump Camel',
    type: 'Multimodal (3 local minima)',
    fn: threeHumpCamel,
    x0s: [
      new Float64Array([2, -1]),
      new Float64Array([-2, 1]),
      new Float64Array([0.5, -0.5]),
    ],
    notes: ['baseline (local basin)', 'mirror (local basin)', 'near global'],
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Rastrigin',
    type: 'Multimodal (highly)',
    fn: rastrigin,
    x0s: [
      new Float64Array([2.5, -3.5]),
      new Float64Array([2.6, -3.4]),
      new Float64Array([0.3, 0.4]),
    ],
    notes: [
      'half-integer aligned — sin(2πxᵢ)=0 zeroes out the multimodal term',
      '+0.1 perturbation — breaks the alignment',
      'near global optimum',
    ],
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Ackley',
    type: 'Multimodal',
    fn: ackley,
    x0s: [
      new Float64Array([2, -2]),
      new Float64Array([3.5, 3.5]),
      new Float64Array([0.5, 0.5]),
    ],
    notes: ['baseline (integer)', 'mid-range non-integer', 'near global'],
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Lévi N.13',
    type: 'Multimodal',
    fn: levi13,
    x0s: [
      new Float64Array([-4, 5]),
      new Float64Array([-3.9, 5.1]),
      new Float64Array([1.5, 1.5]),
    ],
    notes: [
      'integer aligned — all sin(kπxᵢ) terms zero',
      '+0.1 perturbation — breaks the alignment',
      'near global optimum',
    ],
    knownMin: 0,
    knownPoint: new Float64Array([1, 1]),
  },
  {
    name: 'Griewank',
    type: 'Multimodal (shallow ripples)',
    fn: griewank,
    x0s: [
      new Float64Array([100, -200]),
      new Float64Array([5, 5]),
      new Float64Array([0.5, -0.5]),
    ],
    notes: ['baseline (far)', 'moderate', 'near global'],
    knownMin: 0,
    knownPoint: new Float64Array([0, 0]),
  },
  {
    name: 'Styblinski-Tang',
    type: 'Multimodal',
    fn: styblinskiTang,
    x0s: [
      new Float64Array([0, 0]),
      new Float64Array([3, 3]),
      new Float64Array([-2, -2]),
    ],
    notes: ['baseline', 'opposite-sign basin', 'near global'],
    knownMin: -39.16617 * 2,
    knownPoint: new Float64Array([-2.903534, -2.903534]),
  },
  {
    name: 'Easom',
    type: 'Nearly flat, narrow peak',
    fn: easom,
    x0s: [
      new Float64Array([1, 1]),
      new Float64Array([2.5, 2.5]),
      new Float64Array([3, 3]),
    ],
    notes: ['baseline (flat region)', 'closer to peak', 'near peak'],
    knownMin: -1,
    knownPoint: new Float64Array([Math.PI, Math.PI]),
  },
  {
    name: 'Goldstein-Price',
    type: 'Multimodal',
    fn: goldsteinPrice,
    x0s: [
      new Float64Array([0, 0]),
      new Float64Array([1, 1]),
      new Float64Array([-0.5, -1]),
    ],
    notes: ['baseline', 'away from global', 'near global'],
    knownMin: 3,
    knownPoint: new Float64Array([0, -1]),
  },
  {
    name: 'McCormick',
    type: 'Multimodal',
    fn: mccormick,
    x0s: [
      new Float64Array([0, 0]),
      new Float64Array([-1, -2]),
      new Float64Array([2, 2]),
    ],
    notes: [
      'baseline — PSO diverges here without bounds',
      'near global',
      'opposite side',
    ],
    knownMin: -1.9133,
    knownPoint: new Float64Array([-0.54719, -1.54719]),
  },
];

/* ================================================================== */
/*  Optimizer configurations                                           */
/* ================================================================== */

interface OptimizerConfig {
  name: string;
  optimizer: Optimizer<any>;
  settings: Record<string, unknown>;
}

const optimizers: OptimizerConfig[] = [
  {name: 'Nelder-Mead', optimizer: new NelderMead(), settings: {maxIterations: 10_000}},
  {name: 'PSO', optimizer: new PSO(), settings: {maxIterations: 10_000, swarmSize: 50, seed: 42}},
  {
    name: 'GradientDescent', optimizer: new GradientDescent(),
    settings: {maxIterations: 10_000, learningRate: 0.001, momentum: 0.9},
  },
  {name: 'Adam', optimizer: new Adam(), settings: {maxIterations: 10_000, learningRate: 0.1}},
  {name: 'L-BFGS', optimizer: new LBFGS(), settings: {maxIterations: 1_000, historySize: 10}},
];

/* ================================================================== */
/*  Formatting                                                         */
/* ================================================================== */

function pad(s: string, len: number): string {
  return s.length >= len ? s : s + ' '.repeat(len - s.length);
}

function padLeft(s: string, len: number): string {
  return s.length >= len ? s : ' '.repeat(len - s.length) + s;
}

function fmtNum(v: number, decimals: number = 4): string {
  if (!isFinite(v)) return String(v);
  if (v === 0) return '0';
  if (Math.abs(v) < 1e-3) return v.toExponential(2);
  if (Math.abs(v) >= 1e4) return v.toExponential(2);
  return v.toFixed(decimals);
}

function fmtPoint(p: Float64Array): string {
  return '(' + Array.from(p).map((v) => fmtNum(v, 3)).join(', ') + ')';
}

function distance(a: Float64Array, b: Float64Array): number {
  let s = 0;
  for (let i = 0; i < a.length; i++) s += (a[i] - b[i]) ** 2;
  return Math.sqrt(s);
}

function distToOptimum(p: MultistartProblem, point: Float64Array): number {
  if (p.name === 'Himmelblau')
    return Math.min(...HIMMELBLAU_MINIMA.map((m) => distance(point, m)));
  return distance(point, p.knownPoint);
}

function statusIcon(r: OptimizationResult, p: MultistartProblem): '\u2705' | '\u26A0' | '\u274C' {
  if (!r.converged) return '\u274C';
  const absErr = Math.abs(r.value - p.knownMin);
  const distOpt = distToOptimum(p, r.point);
  return Math.max(absErr, distOpt) < 1e-3 ? '\u2705' : '\u26A0';
}

/* ================================================================== */
/*  Per-problem run                                                    */
/* ================================================================== */

interface RunOutcome {
  icon: string;
  value: string;
  iters: string;
  fnEvals: string;
  timeMs: string;
  foundGlobal: boolean;
}

function runOne(
  cfg: OptimizerConfig, fn: ObjectiveFunction, x0: Float64Array, p: MultistartProblem,
): RunOutcome {
  const counted = withCounter(fn);
  const t0 = performance.now();
  let r: OptimizationResult;
  try {
    r = cfg.optimizer.minimize(counted.fn, Float64Array.from(x0), cfg.settings as any);
  } catch (e: any) {
    const msg = (e.message ?? '').slice(0, 20);
    return {
      icon: '\u274C', value: `ERR: ${msg}`,
      iters: '-', fnEvals: '-', timeMs: '-', foundGlobal: false,
    };
  }
  const elapsed = performance.now() - t0;
  const icon = statusIcon(r, p);
  return {
    icon,
    value: fmtNum(r.value),
    iters: String(r.iterations),
    fnEvals: String(counted.count()),
    timeMs: elapsed.toFixed(1),
    foundGlobal: icon === '\u2705',
  };
}

/* ================================================================== */
/*  Table printer                                                      */
/* ================================================================== */

function printProblemTable(p: MultistartProblem, results: RunOutcome[][]): void {
  console.log('\u2501'.repeat(78));
  console.log(`  ${p.name}`);
  console.log(`  Type: ${p.type}`);
  console.log(`  Known min: f${fmtPoint(p.knownPoint)} = ${fmtNum(p.knownMin)}`);
  for (let k = 0; k < 3; k++)
    console.log(`  x₀${String.fromCharCode(65 + k)} = ${fmtPoint(p.x0s[k])}   // ${p.notes[k]}`);
  console.log('\u2501'.repeat(78));

  const header = [
    'Method',
    'x₀A status · value (its, fn)',
    'x₀B status · value (its, fn)',
    'x₀C status · value (its, fn)',
  ];
  const rows: string[][] = [header];
  for (let i = 0; i < optimizers.length; i++) {
    const row = [optimizers[i].name];
    for (let k = 0; k < 3; k++) {
      const o = results[i][k];
      row.push(`${o.icon} ${o.value} (${o.iters} it, ${o.fnEvals} fn)`);
    }
    rows.push(row);
  }

  const widths = header.map((_, col) =>
    rows.reduce((w, r) => Math.max(w, r[col].length), 0));

  for (let r = 0; r < rows.length; r++) {
    const line = rows[r].map((cell, col) =>
      col === 0 ? pad(cell, widths[col]) : padLeft(cell, widths[col])).join(' | ');
    console.log('  ' + line);
    if (r === 0)
      console.log('  ' + widths.map((w) => '-'.repeat(w)).join('-+-'));
  }
  console.log('');
}

/* ================================================================== */
/*  Aggregate summary                                                  */
/* ================================================================== */

function printSummary(successByOpt: Map<string, number>, totalRuns: number): void {
  console.log('\u2501'.repeat(78));
  console.log('  SUMMARY — global optimum hit count across all (problem, x₀) pairs');
  console.log(`  Total runs per optimizer: ${totalRuns}`);
  console.log('\u2501'.repeat(78));
  const header = ['Method', 'Global opts found', 'Success rate'];
  const rows: string[][] = [header];
  for (const cfg of optimizers) {
    const hits = successByOpt.get(cfg.name) ?? 0;
    const pct = (100 * hits / totalRuns).toFixed(1);
    rows.push([cfg.name, `${hits} / ${totalRuns}`, `${pct}%`]);
  }
  const widths = header.map((_, col) =>
    rows.reduce((w, r) => Math.max(w, r[col].length), 0));
  for (let r = 0; r < rows.length; r++) {
    const line = rows[r].map((cell, col) =>
      col === 0 ? pad(cell, widths[col]) : padLeft(cell, widths[col])).join(' | ');
    console.log('  ' + line);
    if (r === 0)
      console.log('  ' + widths.map((w) => '-'.repeat(w)).join('-+-'));
  }
  console.log('');
}

/* ================================================================== */
/*  Main                                                               */
/* ================================================================== */

function runBenchmarks(): void {
  console.log('');
  console.log('\u2554' + '\u2550'.repeat(66) + '\u2557');
  console.log('\u2551          MULTI-START UNCONSTRAINED BENCHMARKS                    \u2551');
  console.log('\u2551          5 optimizers · 15 problems · 3 starting points = 225 runs\u2551');
  console.log('\u255A' + '\u2550'.repeat(66) + '\u255D');
  console.log('');

  const successByOpt = new Map<string, number>();
  for (const cfg of optimizers) successByOpt.set(cfg.name, 0);

  let totalRuns = 0;
  for (const p of problems) {
    const results: RunOutcome[][] = [];
    for (const cfg of optimizers) {
      const row: RunOutcome[] = [];
      for (let k = 0; k < 3; k++) {
        const o = runOne(cfg, p.fn, p.x0s[k], p);
        row.push(o);
        if (o.foundGlobal) successByOpt.set(cfg.name, (successByOpt.get(cfg.name) ?? 0) + 1);
        totalRuns++;
      }
      results.push(row);
    }
    printProblemTable(p, results);
  }

  printSummary(successByOpt, totalRuns / optimizers.length);

  console.log('\u2501'.repeat(78));
  console.log('  Legend:  \u2705 found global (converged & quality < 1e-3)');
  console.log('           \u26A0  converged but landed far from global (local minimum)');
  console.log('           \u274C did not converge');
  console.log('\u2501'.repeat(78));
}

runBenchmarks();
