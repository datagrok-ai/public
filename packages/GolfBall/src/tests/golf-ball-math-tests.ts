// Golf Ball Flight — Math tests

import {category, test, expect, expectFloat} from '@datagrok-libraries/utils/src/test';
import {mrt, ODEs} from 'diff-grok';

import {createGolfBallODE, GolfBallParams} from '../golf-ball/model';
import {DEFAULTS, solve} from '../golf-ball/core';

// -- Helpers --

/** Evaluates the ODE right-hand side at given state and returns [dx/dt, dy/dt, dvx/dt, dvy/dt] */
function evalFunc(params: GolfBallParams, x: number, y: number,
  vx: number, vy: number): [number, number, number, number] {
  const ode = createGolfBallODE(params);
  const state = new Float64Array([x, y, vx, vy]);
  const out = new Float64Array(4);
  ode.func(0, state, out);
  return [out[0], out[1], out[2], out[3]];
}

// -- ODE right-hand side verification --

category('Math: GB func', () => {
  // v0=70, θ=45: vx0 = 70·cos(45°) ≈ 49.497, vy0 = 70·sin(45°) ≈ 49.497
  // A = π·(0.0427/2)² ≈ 0.001432
  // v = √(49.497² + 49.497²) ≈ 69.997
  // drag = 0.25·1.225·0.001432·69.997 / (2·0.0459) ≈ 0.3364
  // dvx/dt = -0.3364·49.497 ≈ -16.65
  // dvy/dt = -9.81 - 0.3364·49.497 ≈ -26.46
  test('func_01: default params at launch', async () => {
    const thetaRad = 45 * Math.PI / 180;
    const vx0 = 70 * Math.cos(thetaRad);
    const vy0 = 70 * Math.sin(thetaRad);
    const [dx, dy, dvx, dvy] = evalFunc(DEFAULTS, 0, 0, vx0, vy0);
    expectFloat(dx, vx0, 1e-6);
    expectFloat(dy, vy0, 1e-6);

    // Verify drag calculation manually
    const A = Math.PI * (0.0427 / 2) ** 2;
    const speed = Math.sqrt(vx0 * vx0 + vy0 * vy0);
    const drag = 0.25 * 1.225 * A * speed / (2 * 0.0459);
    expectFloat(dvx, -drag * vx0, 0.01);
    expectFloat(dvy, -9.81 - drag * vy0, 0.01);
  });

  // At zero velocity → drag = 0, only gravity
  test('func_02: zero velocity → only gravity', async () => {
    const [dx, dy, dvx, dvy] = evalFunc(DEFAULTS, 100, 50, 0, 0);
    expectFloat(dx, 0, 1e-12);
    expectFloat(dy, 0, 1e-12);
    expectFloat(dvx, 0, 1e-12);
    expectFloat(dvy, -9.81, 1e-6);
  });

  // v0=50, θ=30: vx0 = 50·cos(30°) ≈ 43.301, vy0 = 50·sin(30°) = 25.0
  test('func_03: v0=50, θ=30', async () => {
    const params: GolfBallParams = {...DEFAULTS, v0: 50, theta: 30};
    const thetaRad = 30 * Math.PI / 180;
    const vx0 = 50 * Math.cos(thetaRad);
    const vy0 = 50 * Math.sin(thetaRad);
    const [dx, dy, dvx, dvy] = evalFunc(params, 0, 0, vx0, vy0);
    expectFloat(dx, vx0, 1e-6);
    expectFloat(dy, vy0, 1e-6);

    const A = Math.PI * (0.0427 / 2) ** 2;
    const speed = Math.sqrt(vx0 * vx0 + vy0 * vy0);
    const drag = 0.25 * 1.225 * A * speed / (2 * 0.0459);
    expectFloat(dvx, -drag * vx0, 0.01);
    expectFloat(dvy, -DEFAULTS.g - drag * vy0, 0.01);
  });

  // Horizontal-only velocity: vy=0 → dvy = -g (no vertical drag component)
  test('func_04: horizontal velocity only', async () => {
    const [_dx, _dy, dvx, dvy] = evalFunc(DEFAULTS, 0, 10, 50, 0);
    expect(dvx < 0, true, 'dvx should be negative (drag opposes motion)');
    expectFloat(dvy, -DEFAULTS.g, 1e-6);
  });
});

// -- Output property verification --

category('Math: GB Solve properties', () => {
  test('solve_01: default parameters produce non-empty arrays of equal length', async () => {
    const result = solve(DEFAULTS);
    expect(result.t.length > 0, true, 't should be non-empty');
    expect(result.x.length > 0, true, 'x should be non-empty');
    expect(result.y.length > 0, true, 'y should be non-empty');
    expect(result.v.length > 0, true, 'v should be non-empty');
    expect(result.t.length, result.x.length, 't and x should have equal length');
    expect(result.t.length, result.y.length, 't and y should have equal length');
    expect(result.t.length, result.v.length, 't and v should have equal length');
  });

  test('solve_02: y values are non-negative during flight', async () => {
    const result = solve(DEFAULTS);
    for (let i = 0; i < result.y.length; i++)
      expect(result.y[i] >= -0.01, true, `y[${i}] = ${result.y[i]} is negative`);
  });

  test('solve_03: initial conditions preserved', async () => {
    const result = solve(DEFAULTS);
    expectFloat(result.x[0], 0, 1e-6);
    expectFloat(result.y[0], 0, 1e-6);
  });

  test('solve_04: t starts at 0', async () => {
    const result = solve(DEFAULTS);
    expectFloat(result.t[0], 0, 1e-12);
  });

  test('solve_05: v(t) = sqrt(vx^2 + vy^2) consistency', async () => {
    const result = solve(DEFAULTS);
    // v[0] should be v0
    expectFloat(result.v[0], DEFAULTS.v0, 0.5);
  });

  test('solve_06: summary stats are reasonable', async () => {
    const result = solve(DEFAULTS);
    expect(result.maxHeight > 0, true, 'maxHeight should be positive');
    expect(result.maxHeightTime > 0, true, 'maxHeightTime should be positive');
    expect(result.flightDistance > 0, true, 'flightDistance should be positive');
    expect(result.flightTime > 0, true, 'flightTime should be positive');
    expect(result.landingSpeed > 0, true, 'landingSpeed should be positive');
    expect(result.landingAngle > 0, true, 'landingAngle should be positive');
    expect(result.landingAngle < 90, true, 'landingAngle should be < 90');
  });

  test('solve_07: landing y is approximately 0', async () => {
    const result = solve(DEFAULTS);
    const lastY = result.y[result.y.length - 1];
    expectFloat(lastY, 0, 0.1);
  });

  test('solve_08: custom params produce different results', async () => {
    const result1 = solve(DEFAULTS);
    const result2 = solve({...DEFAULTS, theta: 30});
    expect(Math.abs(result1.flightDistance - result2.flightDistance) > 1, true,
      'Different theta should produce different distances');
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
