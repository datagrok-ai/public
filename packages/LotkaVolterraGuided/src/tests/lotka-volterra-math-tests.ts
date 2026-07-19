// Lotka-Volterra — Math tests

import {category, test, expect, expectFloat} from '@datagrok-libraries/utils/src/test';
import {mrt, ODEs} from 'diff-grok';

import {createLotkaVolterraODE, LotkaVolterraParams} from '../lotka-volterra/model';
import {DEFAULTS, solve, getEquilibrium} from '../lotka-volterra/core';

// -- Helpers --

/** Evaluates the ODE right-hand side at given (x, y) and returns [dx/dt, dy/dt] */
function evalFunc(params: LotkaVolterraParams, x: number, y: number): [number, number] {
  const ode = createLotkaVolterraODE(params);
  const state = new Float64Array([x, y]);
  const out = new Float64Array(2);
  ode.func(0, state, out);
  return [out[0], out[1]];
}

// -- ODE right-hand side verification --

const BASE: LotkaVolterraParams = {
  alpha: 1.0, beta: 0.1, delta: 0.075, gamma: 1.5,
  x0: 10, y0: 5, T: 100,
};

category('Math: LV func', () => {
  // dx/dt = 1.0·10 − 0.1·10·5 = 10 − 5 = 5.0
  // dy/dt = 0.075·10·5 − 1.5·5 = 3.75 − 7.5 = −3.75
  test('func_01: default params, (x=10, y=5)', async () => {
    const [dx, dy] = evalFunc(BASE, 10, 5);
    expectFloat(dx, 5.0, 1e-12);
    expectFloat(dy, -3.75, 1e-12);
  });

  // At equilibrium x*=γ/δ=20, y*=α/β=10: dx/dt=0, dy/dt=0
  test('func_02: equilibrium (x*=20, y*=10)', async () => {
    const [dx, dy] = evalFunc(BASE, 20, 10);
    expectFloat(dx, 0.0, 1e-10);
    expectFloat(dy, 0.0, 1e-10);
  });

  // dx/dt = 1.0·30 − 0.1·30·4 = 30 − 12 = 18
  // dy/dt = 0.075·30·4 − 1.5·4 = 9 − 6 = 3.0
  test('func_03: (x=30, y=4)', async () => {
    const [dx, dy] = evalFunc(BASE, 30, 4);
    expectFloat(dx, 18.0, 1e-12);
    expectFloat(dy, 3.0, 1e-12);
  });

  // α=2.0, β=0.5, δ=0.3, γ=1.0
  // dx/dt = 2.0·5 − 0.5·5·2 = 10 − 5 = 5.0
  // dy/dt = 0.3·5·2 − 1.0·2 = 3 − 2 = 1.0
  test('func_04: different params', async () => {
    const params: LotkaVolterraParams = {...BASE, alpha: 2.0, beta: 0.5, delta: 0.3, gamma: 1.0};
    const [dx, dy] = evalFunc(params, 5, 2);
    expectFloat(dx, 5.0, 1e-12);
    expectFloat(dy, 1.0, 1e-12);
  });

  // α=0.5, β=0.02, δ=0.01, γ=0.3
  // dx/dt = 0.5·50 − 0.02·50·10 = 25 − 10 = 15.0
  // dy/dt = 0.01·50·10 − 0.3·10 = 5 − 3 = 2.0
  test('func_05: another param set', async () => {
    const params: LotkaVolterraParams = {...BASE, alpha: 0.5, beta: 0.02, delta: 0.01, gamma: 0.3};
    const [dx, dy] = evalFunc(params, 50, 10);
    expectFloat(dx, 15.0, 1e-12);
    expectFloat(dy, 2.0, 1e-12);
  });
});

// -- Equilibrium verification --

category('Math: LV Equilibrium', () => {
  test('eq_01: default params x*=20, y*=10', async () => {
    const eq = getEquilibrium(1.0, 0.1, 0.075, 1.5);
    expectFloat(eq.xStar, 20.0, 1e-10);
    expectFloat(eq.yStar, 10.0, 1e-10);
  });

  test('eq_02: α=2.0, β=0.5, δ=0.3, γ=1.0 → x*=10/3, y*=4', async () => {
    const eq = getEquilibrium(2.0, 0.5, 0.3, 1.0);
    expectFloat(eq.xStar, 10 / 3, 1e-10);
    expectFloat(eq.yStar, 4.0, 1e-10);
  });
});

// -- Output property verification --

category('Math: LV Solve properties', () => {
  test('solve_01: default parameters produce non-empty arrays of equal length', async () => {
    const result = solve(DEFAULTS);
    expect(result.t.length > 0, true, 't should be non-empty');
    expect(result.x.length > 0, true, 'x should be non-empty');
    expect(result.y.length > 0, true, 'y should be non-empty');
    expect(result.t.length, result.x.length, 't and x should have equal length');
    expect(result.t.length, result.y.length, 't and y should have equal length');
  });

  test('solve_02: x,y values are non-negative', async () => {
    const result = solve(DEFAULTS);
    for (let i = 0; i < result.x.length; i++) {
      expect(result.x[i] >= -0.01, true, `x[${i}] = ${result.x[i]} is negative`);
      expect(result.y[i] >= -0.01, true, `y[${i}] = ${result.y[i]} is negative`);
    }
  });

  test('solve_03: initial conditions preserved', async () => {
    const result = solve(DEFAULTS);
    expectFloat(result.x[0], DEFAULTS.x0, 1e-6);
    expectFloat(result.y[0], DEFAULTS.y0, 1e-6);
  });

  test('solve_04: t starts at 0', async () => {
    const result = solve(DEFAULTS);
    expectFloat(result.t[0], 0, 1e-12);
  });

  test('solve_05: equilibrium start → stays near equilibrium', async () => {
    const eq = getEquilibrium(DEFAULTS.alpha, DEFAULTS.beta, DEFAULTS.delta, DEFAULTS.gamma);
    const params = {...DEFAULTS, x0: eq.xStar, y0: eq.yStar, T: 50};
    const result = solve(params);
    const lastX = result.x[result.x.length - 1];
    const lastY = result.y[result.y.length - 1];
    expectFloat(lastX, eq.xStar, 1.0);
    expectFloat(lastY, eq.yStar, 1.0);
  });

  test('solve_06: summary stats computed correctly', async () => {
    const result = solve(DEFAULTS);
    expect(result.maxPrey > 0, true, 'maxPrey should be positive');
    expect(result.maxPredators > 0, true, 'maxPredators should be positive');
    expect(result.stepCount, result.t.length, 'stepCount should match t.length');
    expect(result.maxPrey >= DEFAULTS.x0, true, 'maxPrey should be at least x0');
  });

  test('solve_07: custom initial conditions', async () => {
    const result = solve({...DEFAULTS, x0: 40, y0: 15});
    expectFloat(result.x[0], 40, 1e-6);
    expectFloat(result.y[0], 15, 1e-6);
  });
});

// -- MRT solver verification --

category('Math: MRT solver', () => {
  test('Non-stiff 1D: dy/dt = 4\u00B7exp(0.8t) \u2212 0.5y', async () => {
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

    const solution = mrt(odes);
    const tArr = solution[0];
    const yArr = solution[1];
    let maxError = 0;
    for (let i = 0; i < tArr.length; i++)
      maxError = Math.max(maxError, Math.abs(exact(tArr[i]) - yArr[i]));

    expect(maxError < 0.1, true, `Max error ${maxError} exceeds 0.1`);
  });

  test('Stiff 1D: dy/dt = \u22121000y + 3000 \u2212 2000\u00B7exp(\u2212t)', async () => {
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

    const solution = mrt(odes);
    const tArr = solution[0];
    const yArr = solution[1];
    let maxError = 0;
    for (let i = 0; i < tArr.length; i++)
      maxError = Math.max(maxError, Math.abs(exact(tArr[i]) - yArr[i]));

    expect(maxError < 0.1, true, `Max error ${maxError} exceeds 0.1`);
  });
});
