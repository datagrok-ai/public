// CR3BP — Math tests

import {category, test, expect} from '@datagrok-libraries/utils/src/test';
import {mrt, ODEs} from 'diff-grok';

import {evalRHS, computeLagrangePoints, jacobiConstant, effectivePotential} from '../cr3bp/model';
import {DEFAULTS, solve} from '../cr3bp/core';

// -- Helpers --

/** Asserts that actual is within tolerance of expected */
function expectFloat(actual: number, expected: number, tol: number): void {
  const diff = Math.abs(actual - expected);
  expect(diff < tol, true, `Expected ${expected}, got ${actual}, diff ${diff} exceeds tolerance ${tol}`);
}

// -- ODE right-hand side verification --

const EARTH_MOON_MU = 0.01215;

category('Math: CR3BP RHS', () => {
  // At (x=0.8, y=0, vx=0, vy=0), mu=0.01215:
  // r1 = |0.8 + 0.01215| = 0.81215
  // r2 = |0.8 - 1 + 0.01215| = 0.18785
  // dx/dt = 0, dy/dt = 0
  // dvx/dt = 2*0 + 0.8 - 0.98785*0.81215/0.81215^3 - 0.01215*(-0.18785)/0.18785^3
  // dvy/dt = -2*0 + 0 - 0 = 0
  test('rhs_01: mu=0.01215, (0.8, 0, 0, 0)', async () => {
    const [dx, dy, dvx, dvy] = evalRHS(EARTH_MOON_MU, 0.8, 0, 0, 0);
    expectFloat(dx, 0, 1e-12);
    expectFloat(dy, 0, 1e-12);
    // dvx = 0.8 - 0.98785/0.81215^2 + 0.01215/0.18785^2
    const r1 = 0.81215;
    const r2 = 0.18785;
    const expectedDvx = 0.8 - 0.98785 / (r1 * r1) + 0.01215 / (r2 * r2);
    expectFloat(dvx, expectedDvx, 1e-4);
    expectFloat(dvy, 0, 1e-12);
  });

  // At L4 with mu=0.01215: (0.48785, sqrt(3)/2, 0, 0)
  // At a Lagrange point with zero velocity, acceleration should be zero
  test('rhs_02: at L4 equilibrium (mu=0.01215)', async () => {
    const xL4 = 0.5 - EARTH_MOON_MU;
    const yL4 = Math.sqrt(3) / 2;
    const [dx, dy, dvx, dvy] = evalRHS(EARTH_MOON_MU, xL4, yL4, 0, 0);
    expectFloat(dx, 0, 1e-12);
    expectFloat(dy, 0, 1e-12);
    // At L4, dvx and dvy should be near zero (equilibrium)
    expectFloat(dvx, 0, 1e-6);
    expectFloat(dvy, 0, 1e-6);
  });

  // At L5 with mu=0.01215: (0.48785, -sqrt(3)/2, 0, 0)
  test('rhs_03: at L5 equilibrium (mu=0.01215)', async () => {
    const xL5 = 0.5 - EARTH_MOON_MU;
    const yL5 = -Math.sqrt(3) / 2;
    const [dx, dy, dvx, dvy] = evalRHS(EARTH_MOON_MU, xL5, yL5, 0, 0);
    expectFloat(dx, 0, 1e-12);
    expectFloat(dy, 0, 1e-12);
    expectFloat(dvx, 0, 1e-6);
    expectFloat(dvy, 0, 1e-6);
  });

  // Equal mass case: mu=0.5, at L4 (0, sqrt(3)/2, 0, 0)
  test('rhs_04: equal mass mu=0.5, at L4', async () => {
    const [dx, dy, dvx, dvy] = evalRHS(0.5, 0.0, Math.sqrt(3) / 2, 0, 0);
    expectFloat(dx, 0, 1e-12);
    expectFloat(dy, 0, 1e-12);
    expectFloat(dvx, 0, 1e-6);
    expectFloat(dvy, 0, 1e-6);
  });

  // With nonzero velocity
  test('rhs_05: nonzero velocity adds Coriolis terms', async () => {
    const [dx, dy, dvx, dvy] = evalRHS(EARTH_MOON_MU, 0.5, 0, 1, 0);
    expectFloat(dx, 1.0, 1e-12); // dx/dt = vx
    expectFloat(dy, 0.0, 1e-12); // dy/dt = vy
    // dvx/dt = 2*vy + ... → Coriolis term is 2*0 = 0
    // dvy/dt = -2*vx + ... → Coriolis term is -2*1 = -2
    // Just check that dvy has the -2 Coriolis contribution
    expect(dvy < -1.5, true, 'dvy should include -2*vx Coriolis term');
  });
});

// -- Lagrange point verification --

category('Math: Lagrange Points', () => {
  test('lp_01: L4 analytical (mu=0.01215)', async () => {
    const lps = computeLagrangePoints(EARTH_MOON_MU);
    const l4 = lps.find((lp) => lp.name === 'L4')!;
    expectFloat(l4.x, 0.5 - EARTH_MOON_MU, 1e-10);
    expectFloat(l4.y, Math.sqrt(3) / 2, 1e-10);
  });

  test('lp_02: L5 analytical (mu=0.01215)', async () => {
    const lps = computeLagrangePoints(EARTH_MOON_MU);
    const l5 = lps.find((lp) => lp.name === 'L5')!;
    expectFloat(l5.x, 0.5 - EARTH_MOON_MU, 1e-10);
    expectFloat(l5.y, -Math.sqrt(3) / 2, 1e-10);
  });

  test('lp_03: L1 between bodies (mu=0.01215)', async () => {
    const lps = computeLagrangePoints(EARTH_MOON_MU);
    const l1 = lps.find((lp) => lp.name === 'L1')!;
    // L1 must be between the two bodies: -mu < x < 1-mu
    expect(l1.x > -EARTH_MOON_MU, true, 'L1 should be > -mu');
    expect(l1.x < 1 - EARTH_MOON_MU, true, 'L1 should be < 1-mu');
    expectFloat(l1.y, 0, 1e-12);
    // L1 for Earth-Moon is approximately 0.8369 (near literature value)
    expectFloat(l1.x, 0.8369, 0.01);
  });

  test('lp_04: L2 beyond smaller body (mu=0.01215)', async () => {
    const lps = computeLagrangePoints(EARTH_MOON_MU);
    const l2 = lps.find((lp) => lp.name === 'L2')!;
    // L2 must be beyond the smaller body: x > 1-mu
    expect(l2.x > 1 - EARTH_MOON_MU, true, 'L2 should be > 1-mu');
    expectFloat(l2.y, 0, 1e-12);
  });
});

// -- Lagrange equilibrium verification --

category('Math: Lagrange Equilibrium', () => {
  test('eq_01: L1 equilibrium condition (mu=0.01215)', async () => {
    const lps = computeLagrangePoints(EARTH_MOON_MU);
    const l1 = lps.find((lp) => lp.name === 'L1')!;
    // At equilibrium, f(x) = x - (1-mu)(x+mu)/|x+mu|^3 - mu(x-1+mu)/|x-1+mu|^3 should be ~0
    const r1 = Math.abs(l1.x + EARTH_MOON_MU);
    const r2 = Math.abs(l1.x - 1 + EARTH_MOON_MU);
    const f = l1.x - (1 - EARTH_MOON_MU) * (l1.x + EARTH_MOON_MU) / (r1 ** 3)
      - EARTH_MOON_MU * (l1.x - 1 + EARTH_MOON_MU) / (r2 ** 3);
    expectFloat(Math.abs(f), 0, 1e-10);
  });

  test('eq_02: L2 equilibrium condition (mu=0.01215)', async () => {
    const lps = computeLagrangePoints(EARTH_MOON_MU);
    const l2 = lps.find((lp) => lp.name === 'L2')!;
    const r1 = Math.abs(l2.x + EARTH_MOON_MU);
    const r2 = Math.abs(l2.x - 1 + EARTH_MOON_MU);
    const f = l2.x - (1 - EARTH_MOON_MU) * (l2.x + EARTH_MOON_MU) / (r1 ** 3)
      - EARTH_MOON_MU * (l2.x - 1 + EARTH_MOON_MU) / (r2 ** 3);
    expectFloat(Math.abs(f), 0, 1e-10);
  });

  test('eq_03: L3 equilibrium condition (mu=0.01215)', async () => {
    const lps = computeLagrangePoints(EARTH_MOON_MU);
    const l3 = lps.find((lp) => lp.name === 'L3')!;
    const r1 = Math.abs(l3.x + EARTH_MOON_MU);
    const r2 = Math.abs(l3.x - 1 + EARTH_MOON_MU);
    const f = l3.x - (1 - EARTH_MOON_MU) * (l3.x + EARTH_MOON_MU) / (r1 ** 3)
      - EARTH_MOON_MU * (l3.x - 1 + EARTH_MOON_MU) / (r2 ** 3);
    expectFloat(Math.abs(f), 0, 1e-10);
  });
});

// -- Jacobi constant verification --

category('Math: Jacobi Constant', () => {
  // Jacobi constant at x=0.8369, y=0, vx=0, vy=0, mu=0.01215 should be ~3.188
  test('cj_01: near L1 state (mu=0.01215)', async () => {
    const cj = jacobiConstant(0.8369, 0, 0, 0, EARTH_MOON_MU);
    expectFloat(cj, 3.188, 0.05);
  });

  test('cj_02: Jacobi constant conservation during integration', async () => {
    const result = solve(DEFAULTS);
    const cj0 = result.cj[0];
    let maxDrift = 0;
    for (let i = 0; i < result.cj.length; i++)
      maxDrift = Math.max(maxDrift, Math.abs(result.cj[i] - cj0));
    expect(maxDrift < 1e-4, true, `Jacobi constant drift ${maxDrift} exceeds 1e-4`);
  });
});

// -- ZVC grid verification --

category('Math: ZVC Grid', () => {
  test('zvc_01: forbidden region consistency', async () => {
    const result = solve(DEFAULTS);
    const grid = result.zvcGrid;
    const cj0 = grid.cj0;

    // All forbidden points should have 2*U > cj0
    for (let i = 0; i < grid.xArr.length; i++)
      expect(grid.U[i] > cj0, true, `Point (${grid.xArr[i]}, ${grid.yArr[i]}): 2*U=${grid.U[i]} should be > C_J=${cj0}`);
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

// -- Solve output properties --

category('Math: CR3BP Solve properties', () => {
  test('solve_01: default params produce non-empty arrays of equal length', async () => {
    const result = solve(DEFAULTS);
    expect(result.t.length > 0, true, 't should be non-empty');
    expect(result.x.length, result.t.length, 'x and t should have equal length');
    expect(result.y.length, result.t.length, 'y and t should have equal length');
    expect(result.vx.length, result.t.length, 'vx and t should have equal length');
    expect(result.vy.length, result.t.length, 'vy and t should have equal length');
    expect(result.vmag.length, result.t.length, 'vmag and t should have equal length');
    expect(result.cj.length, result.t.length, 'cj and t should have equal length');
  });

  test('solve_02: initial conditions preserved', async () => {
    const result = solve(DEFAULTS);
    expectFloat(result.x[0], DEFAULTS.x0, 1e-6);
    expectFloat(result.y[0], DEFAULTS.y0, 1e-6);
    expectFloat(result.vx[0], DEFAULTS.vx0, 1e-6);
    expectFloat(result.vy[0], DEFAULTS.vy0, 1e-6);
  });

  test('solve_03: t starts at 0', async () => {
    const result = solve(DEFAULTS);
    expectFloat(result.t[0], 0, 1e-12);
  });

  test('solve_04: vmag computed correctly', async () => {
    const result = solve(DEFAULTS);
    for (let i = 0; i < Math.min(10, result.t.length); i++) {
      const expectedVmag = Math.sqrt(result.vx[i] ** 2 + result.vy[i] ** 2);
      expectFloat(result.vmag[i], expectedVmag, 1e-12);
    }
  });

  test('solve_05: 5 Lagrange points computed', async () => {
    const result = solve(DEFAULTS);
    expect(result.lagrangePoints.length, 5, 'Should have 5 Lagrange points');
    const names = result.lagrangePoints.map((lp) => lp.name).sort();
    expect(names.join(','), 'L1,L2,L3,L4,L5', 'Should have L1-L5');
  });

  test('solve_06: L4 stability — start near L4 with small mu, stay near L4', async () => {
    const mu = 0.01;
    const xL4 = 0.5 - mu;
    const yL4 = Math.sqrt(3) / 2;
    const result = solve({mu, x0: xL4, y0: yL4, vx0: 0, vy0: 0, T: 100});
    const lastX = result.x[result.x.length - 1];
    const lastY = result.y[result.y.length - 1];
    const dist = Math.sqrt((lastX - xL4) ** 2 + (lastY - yL4) ** 2);
    expect(dist < 0.5, true, `Body drifted ${dist} from L4, should stay near L4 for small mu`);
  });
});
