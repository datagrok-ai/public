// Levins Metapopulation Model — Math tests

import {category, test, expect, expectFloat} from '@datagrok-libraries/utils/src/test';
import {mrt, ODEs} from 'diff-grok';

import {createLevinsODE, LevinsParams} from '../levins/model';
import {DEFAULTS, solve, getEquilibrium} from '../levins/core';

// ── Helpers ──

/** Max absolute error between numerical and exact solutions across all grid points */
function getMaxError(odes: ODEs, exact: (t: number) => number): number {
  const solution = mrt(odes);
  const tArr = solution[0];
  const yArr = solution[1];
  let error = 0;

  for (let i = 0; i < tArr.length; i++)
    error = Math.max(error, Math.abs(exact(tArr[i]) - yArr[i]));

  return error;
}

/** Evaluates func at given p and returns dp/dt */
function evalFunc(params: LevinsParams, p: number): number {
  const ode = createLevinsODE(params);
  const y = new Float64Array([p]);
  const out = new Float64Array(1);
  ode.func(0, y, out);
  return out[0];
}

// ── Correctness: MRT solver ──

const TINY = 0.1;

category('Math: MRT solver', () => {
  test('Non-stiff 1D: dy/dt = 4·exp(0.8t) − 0.5y', async () => {
    // Reference: Chapra & Canale, p. 736
    const odes: ODEs = {
      name: 'Non-stiff 1D',
      arg: {name: 't', start: 0, finish: 4, step: 0.01},
      initial: [2],
      func: (_t: number, y: Float64Array, out: Float64Array) => {
        out[0] = 4 * Math.exp(0.8 * _t) - 0.5 * y[0];
      },
      tolerance: 1e-6,
      solutionColNames: ['y'],
    };

    const exact = (t: number) =>
      (4 / 1.3) * (Math.exp(0.8 * t) - Math.exp(-0.5 * t)) + 2 * Math.exp(-0.5 * t);

    const error = getMaxError(odes, exact);
    expectFloat(error, 0, TINY);
  });

  test('Stiff 1D: dy/dt = −1000y + 3000 − 2000·exp(−t)', async () => {
    // Reference: Chapra & Canale, p. 767
    const odes: ODEs = {
      name: 'Stiff 1D',
      arg: {name: 't', start: 0, finish: 4, step: 0.01},
      initial: [0],
      func: (_t: number, y: Float64Array, out: Float64Array) => {
        out[0] = -1000 * y[0] + 3000 - 2000 * Math.exp(-_t);
      },
      tolerance: 5e-7,
      solutionColNames: ['y'],
    };

    const exact = (t: number) =>
      3 - 0.998 * Math.exp(-1000 * t) - 2.002 * Math.exp(-t);

    const error = getMaxError(odes, exact);
    expectFloat(error, 0, TINY);
  });

  test('Stiff 2D: VDPOL (van der Pol, µ=1000)', async () => {
    // Reference: https://archimede.uniba.it/~testset/report/vdpol.pdf
    const vdpol: ODEs = {
      name: 'van der Pol',
      arg: {name: 't', start: 0, finish: 2000, step: 0.1},
      initial: [-1, 1],
      func: (_t: number, y: Float64Array, out: Float64Array) => {
        out[0] = y[1];
        out[1] = -y[0] + 1000 * (1 - y[0] * y[0]) * y[1];
      },
      tolerance: 1e-12,
      solutionColNames: ['x1', 'x2'],
    };

    mrt(vdpol);
  }, {benchmark: true, timeout: 2000});
});

// ── Correctness: Levins func ──

const BASE: LevinsParams = {
  p0: 0.5, m: 0.5, e0: 0.2, rescueEffect: false,
  t_start: 0, t_end: 50, t_step: 0.1, tolerance: 1e-7,
};

category('Math: Levins func', () => {
  // dp/dt = m·p·(1−p) − e₀·p = 0.5·0.5·0.5 − 0.2·0.5 = 0.125 − 0.1 = 0.025
  test('func_01: base model, p=0.5', async () => {
    expectFloat(evalFunc({...BASE, m: 0.5, e0: 0.2, rescueEffect: false}, 0.5), 0.025, 1e-12);
  });

  // dp/dt = 1.0·0.1·0.9 − 0.3·0.1 = 0.09 − 0.03 = 0.06
  test('func_02: base model, low p=0.1', async () => {
    expectFloat(evalFunc({...BASE, m: 1.0, e0: 0.3, rescueEffect: false}, 0.1), 0.06, 1e-12);
  });

  // At equilibrium p*=1−e₀/m=0.6: dp/dt = 0.5·0.6·0.4 − 0.2·0.6 = 0.12 − 0.12 = 0.0
  test('func_03: equilibrium p*=0.6, dp/dt=0', async () => {
    expectFloat(evalFunc({...BASE, m: 0.5, e0: 0.2, rescueEffect: false}, 0.6), 0.0, 1e-12);
  });

  // Rescue: e=e₀·(1−p)=0.2·0.5=0.1; dp/dt = 0.5·0.5·0.5 − 0.1·0.5 = 0.125 − 0.05 = 0.075
  test('func_04: rescue effect, p=0.5', async () => {
    expectFloat(evalFunc({...BASE, m: 0.5, e0: 0.2, rescueEffect: true}, 0.5), 0.075, 1e-12);
  });

  // Rescue: e=0.5·(1−0.8)=0.1; dp/dt = 0.3·0.8·0.2 − 0.1·0.8 = 0.048 − 0.08 = −0.032
  test('func_05: rescue + decline, p=0.8', async () => {
    expectFloat(evalFunc({...BASE, m: 0.3, e0: 0.5, rescueEffect: true}, 0.8), -0.032, 1e-12);
  });
});

// ── Correctness: Levins equilibrium ──

category('Math: Equilibrium', () => {
  test('p* = 1 - e0/m (base model)', async () => {
    expectFloat(getEquilibrium(0.5, 0.2, false), 0.6, 1e-10);
  });

  test('p* = 0 when m <= e0', async () => {
    expectFloat(getEquilibrium(0.2, 0.5, false), 0, 1e-10);
  });

  test('p* = 0 when m = e0', async () => {
    expectFloat(getEquilibrium(0.5, 0.5, false), 0, 1e-10);
  });

  test('p* = NaN with rescue effect', async () => {
    expect(isNaN(getEquilibrium(0.5, 0.2, true)), true, 'Should be NaN with rescue');
  });
});

// ── Output property verification: solve ──

category('Math: Solve output properties', () => {
  test('solve_01: default parameters produce non-empty arrays of equal length', async () => {
    const result = solve(DEFAULTS);
    expect(result.t.length > 0, true, 't should be non-empty');
    expect(result.p.length > 0, true, 'p should be non-empty');
    expect(result.t.length, result.p.length, 't and p should have equal length');
  });

  test('solve_02: p values in [0, 1]', async () => {
    const result = solve(DEFAULTS);
    for (let i = 0; i < result.p.length; i++)
      expect(result.p[i] >= 0 && result.p[i] <= 1, true, `p[${i}] = ${result.p[i]} out of [0, 1]`);
  });

  test('solve_03: p(0) = p0', async () => {
    const result = solve(DEFAULTS);
    expectFloat(result.p[0], DEFAULTS.p0, 1e-6);
  });

  test('solve_04: t starts at t_start', async () => {
    const result = solve(DEFAULTS);
    expectFloat(result.t[0], DEFAULTS.t_start, 1e-12);
  });

  test('solve_05: convergence to p*', async () => {
    const params = {...DEFAULTS, m: 0.5, e0: 0.2, rescueEffect: false, t_end: 200};
    const result = solve(params);
    const pStar = getEquilibrium(params.m, params.e0, params.rescueEffect);
    expectFloat(result.p[result.p.length - 1], pStar, 0.01);
  });

  test('solve_06: rescue effect — p in [0, 1]', async () => {
    const params = {...DEFAULTS, m: 0.3, e0: 0.5, rescueEffect: true, t_end: 100};
    const result = solve(params);
    for (let i = 0; i < result.p.length; i++)
      expect(result.p[i] >= 0 && result.p[i] <= 1, true, `p[${i}] = ${result.p[i]} out of [0, 1]`);
  });

  test('solve_07: higher m → higher p(t_end)', async () => {
    const r1 = solve({...DEFAULTS, m: 0.5, e0: 0.2});
    const r2 = solve({...DEFAULTS, m: 1.0, e0: 0.2});
    expect(r2.p[r2.p.length - 1] > r1.p[r1.p.length - 1], true,
      'p(t_end) with m=1.0 should exceed p(t_end) with m=0.5');
  });

  test('solve_08: custom p0', async () => {
    const result = solve({...DEFAULTS, p0: 0.9});
    expectFloat(result.p[0], 0.9, 1e-6);
  });
});
