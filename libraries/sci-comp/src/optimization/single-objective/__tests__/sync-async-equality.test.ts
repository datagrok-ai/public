// Bit-for-bit sync↔async equality test.
//
// The codegen produces a sync mirror that is textually identical to the
// async source modulo `await` removal and `Promise<T>` unwrap. Both paths
// therefore execute the same arithmetic in the same order, and under
// IEEE-754 must produce identical results — not just close, *equal*.
//
// This test is the load-bearing safety net for the rollout: if a codegen
// transform ever introduces semantic drift (or if a hand-written sync sibling
// has silently diverged from the async source), this test fails.
//
// Stochastic optimizers (PSO) consume `Math.random()`. Both runs are wrapped
// in `withSeed` so they consume the same deterministic PRNG sequence.

import {
  Adam, GradientDescent, PSO, NelderMead, LBFGS, LBFGSB,
  Optimizer, CommonSettings,
} from '..';
import {rosenbrock, sphere, gaussian, quadratic3d, toAsync} from './helpers';

/* ----- seeded Math.random (Mulberry32) ----- */

function seededRandom(seed: number): () => number {
  let s = seed | 0;
  return () => {
    s = (s + 0x6D2B79F5) | 0;
    let t = Math.imul(s ^ (s >>> 15), 1 | s);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function withSeed<T>(seed: number, fn: () => T): T {
  const original = Math.random;
  Math.random = seededRandom(seed);
  try {return fn();} finally {Math.random = original;}
}

async function withSeedAsync<T>(seed: number, fn: () => Promise<T>): Promise<T> {
  const original = Math.random;
  Math.random = seededRandom(seed);
  try {return await fn();} finally {Math.random = original;}
}

/* ----- fixtures ----- */

const settings: CommonSettings = {maxIterations: 200, tolerance: 1e-10};

const problems = [
  {name: 'rosenbrock-2d', fn: rosenbrock, x0: new Float64Array([-1.2, 1.0])},
  {name: 'sphere-3d', fn: sphere, x0: new Float64Array([5, -3, 7])},
  {name: 'gaussian-2d', fn: gaussian, x0: new Float64Array([0.5, -0.5])},
  {name: 'quadratic-3d', fn: quadratic3d, x0: new Float64Array([0, 0, 0])},
];

const optimizers: Array<{name: string; make: () => Optimizer<any>}> = [
  {name: 'NelderMead', make: () => new NelderMead()},
  {name: 'PSO', make: () => new PSO()},
  {name: 'GradientDescent', make: () => new GradientDescent()},
  {name: 'Adam', make: () => new Adam()},
  {name: 'LBFGS', make: () => new LBFGS()},
  {name: 'LBFGSB', make: () => new LBFGSB()},
];

/* ----- the test ----- */

describe.each(optimizers)('$name: sync ≡ async (bit-equal)', ({make}) => {
  test.each(problems)('$name', async ({fn, x0}) => {
    const opt = make();
    const rs = withSeed(42, () => opt.minimize(fn, x0, settings));
    const ra = await withSeedAsync(42, () => opt.minimizeAsync(toAsync(fn), x0, settings));
    expect(rs.value).toBe(ra.value);
    expect([...rs.point]).toEqual([...ra.point]);
    expect(rs.iterations).toBe(ra.iterations);
    expect(rs.converged).toBe(ra.converged);
    expect([...rs.costHistory]).toEqual([...ra.costHistory]);
  });
});
