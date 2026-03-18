// SIR/SEIR — Math tests

import {category, test, expect, expectFloat} from '@datagrok-libraries/utils/src/test';
import {mrt, ODEs} from 'diff-grok';

import {
  EpidemicParams, createSIR_ODE, createSEIR_ODE, POPULATION,
} from '../sir-seir/model';
import {DEFAULTS, solve} from '../sir-seir/core';

// -- Helpers --

/** Evaluates the SIR ODE RHS at given state and returns [dS/dt, dI/dt, dR/dt] */
function evalSIR(params: EpidemicParams, S: number, I: number, R: number): [number, number, number] {
  const ode = createSIR_ODE(params);
  const state = new Float64Array([S, I, R]);
  const out = new Float64Array(3);
  ode.func(0, state, out);
  return [out[0], out[1], out[2]];
}

/** Evaluates the SEIR ODE RHS at given state and returns [dS/dt, dE/dt, dI/dt, dR/dt] */
function evalSEIR(params: EpidemicParams, S: number, E: number, I: number, R: number): [number, number, number, number] {
  const ode = createSEIR_ODE(params);
  const state = new Float64Array([S, E, I, R]);
  const out = new Float64Array(4);
  ode.func(0, state, out);
  return [out[0], out[1], out[2], out[3]];
}

// -- SIR ODE RHS verification --

category('Math: SIR func', () => {
  // Reference example 1: SIR, R₀=3, γ=0.1, S=9999, I=1, R=0
  // β = 3·0.1 = 0.3
  // dS/dt = -0.3·9999·1/10000 ≈ -0.02999
  // dI/dt = 0.3·9999·1/10000 - 0.1·1 = 0.02999 - 0.1 ≈ -0.07001 ... wait
  // Actually from spec ref example 1: dI/dt = 0.3·9999·1/10000 - 0.1·1 ≈ 0.200
  // Hmm, spec says dI/dt ≈ 0.200 but β·S·I/N = 0.3·9999/10000 ≈ 0.29997, minus γ·I = 0.1
  // So dI/dt = 0.29997 - 0.1 = 0.19997 ≈ 0.200 ✓
  // But dS/dt = -0.29997 ≈ -0.300. Spec table row 1 says dS/dt ≈ -0.300 ✓
  test('func_01: SIR default, S=9999, I=1, R=0', async () => {
    const params: EpidemicParams = {modelType: 'SIR', r0: 3.0, gamma: 0.1, sigma: 0.2, vaccination: 0};
    const [dS, dI, dR] = evalSIR(params, 9999, 1, 0);
    expectFloat(dS, -0.3 * 9999 / 10000, 1e-6);
    expectFloat(dI, 0.3 * 9999 / 10000 - 0.1, 1e-6);
    expectFloat(dR, 0.1, 1e-6);
  });

  // Reference example 3 from ODE table: SIR, R₀=3, γ=0.1, S=5000, I=2000, R=3000
  // β = 0.3
  // dS/dt = -0.3·5000·2000/10000 = -300/10000·10000000/10000... let me compute:
  // -0.3 * 5000 * 2000 / 10000 = -0.3 * 10000000 / 10000 = -0.3 * 1000 = -300
  // Wait, -0.3 * 5000 * 2000 / 10000 = -0.3 * 10_000_000 / 10000 = -300
  // Spec row 3 says dS/dt = -30.0. Let me recheck:
  // -β·S·I/N = -0.3 · 5000 · 2000 / 10000 = -0.3 · 1000 = -300
  // But spec says -30.0. Hmm.
  // Wait, the ODE spec line 219 says: β=0.3; -0.3·5000·2000/10000
  // 5000*2000 = 10,000,000. /10000 = 1000. *0.3 = 300. So dS/dt = -300.
  // Spec says -30.0 for expected dS/dt... that might be a typo in the spec.
  // Actually looking more carefully: the spec row 3: S=5000, I=2000, R=3000
  // β=0.3; -0.3·5000·2000/10000 = -300. Spec says -30.0 which is incorrect.
  // Let me use the formula directly.
  test('func_02: SIR mid-epidemic, S=5000, I=2000, R=3000', async () => {
    const params: EpidemicParams = {modelType: 'SIR', r0: 3.0, gamma: 0.1, sigma: 0.2, vaccination: 0};
    const [dS, dI, dR] = evalSIR(params, 5000, 2000, 3000);
    // β = 0.3, dS/dt = -0.3 * 5000 * 2000 / 10000 = -300
    // dI/dt = 0.3 * 5000 * 2000 / 10000 - 0.1 * 2000 = 300 - 200 = 100
    // dR/dt = 0.1 * 2000 = 200
    expectFloat(dS, -300.0, 1e-6);
    expectFloat(dI, 100.0, 1e-6);
    expectFloat(dR, 200.0, 1e-6);
  });

  // Reference example 5: SIR, R₀=1.5, γ=0.1, S=3300, I=1, R=6699
  // β = 0.15
  // dS/dt = -0.15 * 3300 * 1 / 10000 = -0.0495
  // dI/dt = 0.15 * 3300 / 10000 - 0.1 = 0.0495 - 0.1 = -0.0505
  // R_eff = 1.5 * 3300/10000 = 0.495 < 1 → I decreases
  test('func_03: SIR herd immunity, S=3300, I=1, R=6699', async () => {
    const params: EpidemicParams = {modelType: 'SIR', r0: 1.5, gamma: 0.1, sigma: 0.2, vaccination: 0};
    const [dS, dI, dR] = evalSIR(params, 3300, 1, 6699);
    expectFloat(dS, -0.15 * 3300 / 10000, 1e-6);
    expectFloat(dI, 0.15 * 3300 / 10000 - 0.1, 1e-6);
    expectFloat(dR, 0.1, 1e-6);
    expect(dI < 0, true, 'I should decrease when R_eff < 1');
  });
});

// -- SEIR ODE RHS verification --

category('Math: SEIR func', () => {
  // Reference example 2: SEIR, R₀=3, γ=0.1, σ=0.2, S=9999, E=0, I=1, R=0
  // β = 0.3
  // dS/dt = -0.3·9999·1/10000 ≈ -0.02999
  // dE/dt = 0.3·9999·1/10000 - 0.2·0 ≈ 0.02999
  // dI/dt = 0.2·0 - 0.1·1 = -0.1
  // dR/dt = 0.1·1 = 0.1
  test('func_01: SEIR default, S=9999, E=0, I=1, R=0', async () => {
    const params: EpidemicParams = {modelType: 'SEIR', r0: 3.0, gamma: 0.1, sigma: 0.2, vaccination: 0};
    const [dS, dE, dI, dR] = evalSEIR(params, 9999, 0, 1, 0);
    expectFloat(dS, -0.3 * 9999 / 10000, 1e-6);
    expectFloat(dE, 0.3 * 9999 / 10000, 1e-6);
    expectFloat(dI, -0.1, 1e-6);
    expectFloat(dR, 0.1, 1e-6);
  });

  // Reference example 4: SEIR, R₀=2, γ=0.2, σ=0.5, S=8000, E=500, I=1000, R=500
  // β = 2·0.2 = 0.4
  // dS/dt = -0.4·8000·1000/10000 = -320
  // dE/dt = 0.4·8000·1000/10000 - 0.5·500 = 320 - 250 = 70
  // dI/dt = 0.5·500 - 0.2·1000 = 250 - 200 = 50
  // dR/dt = 0.2·1000 = 200
  // Note: spec says dS=-32, dE=-218, dI=150, dR=200 which doesn't match the formula.
  // Let me use formula directly.
  test('func_02: SEIR mid-epidemic, S=8000, E=500, I=1000, R=500', async () => {
    const params: EpidemicParams = {modelType: 'SEIR', r0: 2.0, gamma: 0.2, sigma: 0.5, vaccination: 0};
    const [dS, dE, dI, dR] = evalSEIR(params, 8000, 500, 1000, 500);
    // β = 0.4
    expectFloat(dS, -0.4 * 8000 * 1000 / 10000, 1e-6);  // -320
    expectFloat(dE, 0.4 * 8000 * 1000 / 10000 - 0.5 * 500, 1e-6);  // 320 - 250 = 70
    expectFloat(dI, 0.5 * 500 - 0.2 * 1000, 1e-6);  // 250 - 200 = 50
    expectFloat(dR, 0.2 * 1000, 1e-6);  // 200
  });

  // Simple case with all zeros except S and I
  test('func_03: SEIR no exposed, S=10000, E=0, I=0, R=0', async () => {
    const params: EpidemicParams = {modelType: 'SEIR', r0: 3.0, gamma: 0.1, sigma: 0.2, vaccination: 0};
    const [dS, dE, dI, dR] = evalSEIR(params, 10000, 0, 0, 0);
    expectFloat(dS, 0, 1e-6);
    expectFloat(dE, 0, 1e-6);
    expectFloat(dI, 0, 1e-6);
    expectFloat(dR, 0, 1e-6);
  });
});

// -- R_eff verification --

category('Math: R_eff', () => {
  test('R_eff at start: R₀=3, S=9999 → R_eff ≈ 3.0', async () => {
    const rEff = 3.0 * 9999 / 10000;
    expectFloat(rEff, 2.9997, 1e-4);
  });

  test('R_eff near herd immunity: R₀=3, S=3333 → R_eff ≈ 1.0', async () => {
    const rEff = 3.0 * 3333 / 10000;
    expectFloat(rEff, 0.9999, 1e-3);
  });
});

// -- Summary stats verification --

category('Math: Summary stats', () => {
  test('HIT = 1 - 1/R₀ for R₀=1.5 → 33.3%', async () => {
    const hit = (1 - 1 / 1.5) * 100;
    expectFloat(hit, 33.333, 0.01);
  });

  test('HIT = 1 - 1/R₀ for R₀=3.0 → 66.7%', async () => {
    const hit = (1 - 1 / 3.0) * 100;
    expectFloat(hit, 66.667, 0.01);
  });
});

// -- SIR solve property verification --

category('Math: SIR solve properties', () => {
  const sirDefaults: EpidemicParams = {...DEFAULTS, modelType: 'SIR'};

  test('solve_01: non-empty arrays of equal length', async () => {
    const result = solve(sirDefaults);
    expect(result.t.length > 0, true, 't should be non-empty');
    expect(result.S.length, result.t.length, 'S and t should have equal length');
    expect(result.I.length, result.t.length, 'I and t should have equal length');
    expect(result.R.length, result.t.length, 'R and t should have equal length');
    expect(result.E, null, 'E should be null for SIR');
  });

  test('solve_02: non-negativity of all compartments', async () => {
    const result = solve(sirDefaults);
    for (let i = 0; i < result.t.length; i++) {
      expect(result.S[i] >= -0.01, true, `S[${i}] = ${result.S[i]} is negative`);
      expect(result.I[i] >= -0.01, true, `I[${i}] = ${result.I[i]} is negative`);
      expect(result.R[i] >= -0.01, true, `R[${i}] = ${result.R[i]} is negative`);
    }
  });

  test('solve_03: initial conditions preserved', async () => {
    const result = solve(sirDefaults);
    expectFloat(result.S[0], POPULATION * (1 - 0) - 1, 1e-3);
    expectFloat(result.I[0], 1, 1e-3);
    expectFloat(result.R[0], 0, 1e-3);
  });

  test('solve_04: population conservation S+I+R=N', async () => {
    const result = solve(sirDefaults);
    for (let i = 0; i < result.t.length; i++) {
      const total = result.S[i] + result.I[i] + result.R[i];
      expectFloat(total, POPULATION, 1.0);
    }
  });
});

// -- SEIR solve property verification --

category('Math: SEIR solve properties', () => {
  const seirDefaults: EpidemicParams = {...DEFAULTS, modelType: 'SEIR'};

  test('solve_01: non-empty arrays with E present', async () => {
    const result = solve(seirDefaults);
    expect(result.t.length > 0, true, 't should be non-empty');
    expect(result.E !== null, true, 'E should not be null for SEIR');
    expect(result.E!.length, result.t.length, 'E and t should have equal length');
  });

  test('solve_02: non-negativity of all compartments', async () => {
    const result = solve(seirDefaults);
    for (let i = 0; i < result.t.length; i++) {
      expect(result.S[i] >= -0.01, true, `S[${i}] = ${result.S[i]} is negative`);
      expect(result.E![i] >= -0.01, true, `E[${i}] = ${result.E![i]} is negative`);
      expect(result.I[i] >= -0.01, true, `I[${i}] = ${result.I[i]} is negative`);
      expect(result.R[i] >= -0.01, true, `R[${i}] = ${result.R[i]} is negative`);
    }
  });

  test('solve_03: initial conditions preserved', async () => {
    const result = solve(seirDefaults);
    expectFloat(result.S[0], POPULATION - 1, 1e-3);
    expectFloat(result.E![0], 0, 1e-3);
    expectFloat(result.I[0], 1, 1e-3);
    expectFloat(result.R[0], 0, 1e-3);
  });

  test('solve_04: population conservation S+E+I+R=N', async () => {
    const result = solve(seirDefaults);
    for (let i = 0; i < result.t.length; i++) {
      const total = result.S[i] + result.E![i] + result.I[i] + result.R[i];
      expectFloat(total, POPULATION, 1.0);
    }
  });
});

// -- R₀ < 1 behavior --

category('Math: R₀ < 1 behavior', () => {
  test('I(t) decreases monotonically when R₀ < 1', async () => {
    const params: EpidemicParams = {
      modelType: 'SIR', r0: 0.8, gamma: 0.1, sigma: 0.2, vaccination: 0,
    };
    const result = solve(params);
    // I should be monotonically non-increasing (with small tolerance for solver)
    for (let i = 1; i < result.I.length; i++)
      expect(result.I[i] <= result.I[i - 1] + 0.01, true, `I[${i}] = ${result.I[i]} > I[${i - 1}] = ${result.I[i - 1]}`);
  });
});

// -- Herd immunity --

category('Math: Herd immunity', () => {
  test('Near-zero epidemic when vaccination >= HIT', async () => {
    // R₀=3, HIT = 1 - 1/3 ≈ 66.7%. Set vaccination to 70%.
    const params: EpidemicParams = {
      modelType: 'SIR', r0: 3.0, gamma: 0.1, sigma: 0.2, vaccination: 70,
    };
    const result = solve(params);
    // Peak I should be very small (close to initial I₀=1)
    expect(result.peakCount <= 10, true, `Peak count ${result.peakCount} too high for herd immunity`);
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
