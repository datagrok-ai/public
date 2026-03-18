// Multi-Compartment PK — Math tests

import {category, test, expect, expectFloat} from '@datagrok-libraries/utils/src/test';
import {mrt, ODEs} from 'diff-grok';

import {
  PkParams, DEFAULTS, solvePk, computeWssr,
  computeMicroConstants, createPkODE,
} from '../compartmental-pk/core';

// -- Helpers --

/** Evaluate ODE RHS at given state and return derivatives */
function evalRHS(
  params: PkParams, state: number[], rIn: number,
): number[] {
  const micro = computeMicroConstants(params);
  const is3comp = params.modelType === '3-Compartment';
  const isOral = params.inputMode === 'Oral';
  const nVars = (is3comp ? 3 : 2) + (isOral ? 1 : 0);

  const ode = createPkODE(params, micro, state, 0, 1, rIn);
  const y = new Float64Array(state);
  const out = new Float64Array(nVars);
  ode.func(0, y, out);
  return Array.from(out);
}

// -- Micro-constants verification --

category('Math: Micro-constants', () => {
  test('micro_01: k10 = CL/Vc', async () => {
    const params = {...DEFAULTS, cl: 5, vc: 20};
    const micro = computeMicroConstants(params);
    expectFloat(micro.k10, 0.25, 1e-12);
  });

  test('micro_02: k12 = Q/Vc, k21 = Q/Vp', async () => {
    const params = {...DEFAULTS, q: 10, vc: 20, vp: 40};
    const micro = computeMicroConstants(params);
    expectFloat(micro.k12, 0.5, 1e-12);
    expectFloat(micro.k21, 0.25, 1e-12);
  });

  test('micro_03: k13, k31 for 3-comp', async () => {
    const params: PkParams = {...DEFAULTS, modelType: '3-Compartment', q2: 2, vc: 20, vd: 60};
    const micro = computeMicroConstants(params);
    expectFloat(micro.k13, 0.1, 1e-12);
    expectFloat(micro.k31, 2 / 60, 1e-12);
  });
});

// -- ODE RHS verification --

category('Math: PK ODE RHS (2-comp, IV Bolus)', () => {
  // Ref #1: CL=5, Vc=20, Vp=40, Q=10, Dose=100, t=0
  // k10=0.25, k12=0.5, k21=0.25
  // dAc/dt = -0.25*100 -0.5*100 +0.25*0 +0 = -75
  // dAp/dt = 0.5*100 -0.25*0 = 50
  test('rhs_01: 2-comp IV bolus at t=0', async () => {
    const params: PkParams = {...DEFAULTS, modelType: '2-Compartment', inputMode: 'IV Bolus',
      cl: 5, vc: 20, vp: 40, q: 10, dose: 100};
    const derivatives = evalRHS(params, [100, 0], 0);
    expectFloat(derivatives[0], -75, 1e-8); // dAc/dt
    expectFloat(derivatives[1], 50, 1e-8);  // dAp/dt
  });

  test('rhs_02: 2-comp IV bolus at later time', async () => {
    // Ac=50, Ap=30: dAc/dt = -0.25*50 -0.5*50 +0.25*30 = -12.5 -25 +7.5 = -30
    // dAp/dt = 0.5*50 -0.25*30 = 25 -7.5 = 17.5
    const params: PkParams = {...DEFAULTS, modelType: '2-Compartment', inputMode: 'IV Bolus',
      cl: 5, vc: 20, vp: 40, q: 10};
    const derivatives = evalRHS(params, [50, 30], 0);
    expectFloat(derivatives[0], -30, 1e-8);
    expectFloat(derivatives[1], 17.5, 1e-8);
  });
});

category('Math: PK ODE RHS (3-comp, IV Bolus)', () => {
  // Ref #2: CL=5, Vc=20, Vp=40, Q=10, Vd=60, Q2=2, Dose=100, t=0
  // k10=0.25, k12=0.5, k21=0.25, k13=0.1, k31=2/60≈0.0333
  // dAc/dt = -0.25*100 -0.5*100 -0.1*100 +0.25*0 +0.0333*0 = -85
  // dAp/dt = 0.5*100 -0.25*0 = 50
  // dAd/dt = 0.1*100 -0.0333*0 = 10
  test('rhs_03: 3-comp IV bolus at t=0', async () => {
    const params: PkParams = {...DEFAULTS, modelType: '3-Compartment', inputMode: 'IV Bolus',
      cl: 5, vc: 20, vp: 40, q: 10, vd: 60, q2: 2, dose: 100};
    const derivatives = evalRHS(params, [100, 0, 0], 0);
    expectFloat(derivatives[0], -85, 1e-8);
    expectFloat(derivatives[1], 50, 1e-8);
    expectFloat(derivatives[2], 10, 1e-8);
  });
});

category('Math: PK ODE RHS (2-comp, Oral)', () => {
  // Ref #3: CL=5, Vc=20, Vp=40, Q=10, ka=1.0, F=0.8, Dose=100, t=0
  // Agut(0)=F*Dose=80, Ac(0)=0, Ap(0)=0
  // dAgut/dt = -1.0*80 = -80
  // R_in = ka*Agut = 1.0*80 = 80
  // dAc/dt = -0.25*0 -0.5*0 +0.25*0 +80 = 80
  // dAp/dt = 0.5*0 -0.25*0 = 0
  test('rhs_04: 2-comp oral at t=0', async () => {
    const params: PkParams = {...DEFAULTS, modelType: '2-Compartment', inputMode: 'Oral',
      cl: 5, vc: 20, vp: 40, q: 10, ka: 1.0, f: 0.8, dose: 100};
    // State: [Ac, Ap, Agut]
    const derivatives = evalRHS(params, [0, 0, 80], 0);
    expectFloat(derivatives[0], 80, 1e-8);   // dAc/dt
    expectFloat(derivatives[1], 0, 1e-8);    // dAp/dt
    expectFloat(derivatives[2], -80, 1e-8);  // dAgut/dt
  });
});

category('Math: PK ODE RHS (2-comp, IV Infusion)', () => {
  // Ref #4: CL=5, Vc=20, Vp=40, Q=10, Dose=100, Tinf=1, t=0.5
  // R_in = 100/1 = 100
  // Ac(0)=0, dAc/dt|0 = -0.25*0 -0.5*0 +0.25*0 +100 = 100
  test('rhs_05: 2-comp infusion at t=0', async () => {
    const params: PkParams = {...DEFAULTS, modelType: '2-Compartment', inputMode: 'IV Infusion',
      cl: 5, vc: 20, vp: 40, q: 10, dose: 100, tInf: 1};
    const derivatives = evalRHS(params, [0, 0], 100);
    expectFloat(derivatives[0], 100, 1e-8);  // dAc/dt
    expectFloat(derivatives[1], 0, 1e-8);    // dAp/dt
  });
});

// -- WSSR computation --

category('Math: WSSR', () => {
  test('wssr_01: uniform weights', async () => {
    const cObs = [5, 3, 1];
    const cPred = [4, 3, 2];
    // WSSR = 1*(5-4)^2 + 1*(3-3)^2 + 1*(1-2)^2 = 1 + 0 + 1 = 2
    const wssr = computeWssr(cObs, cPred, 'Uniform');
    expectFloat(wssr, 2.0, 1e-10);
  });

  test('wssr_02: 1/Cobs² weights', async () => {
    const cObs = [4, 2];
    const cPred = [5, 1];
    // w1 = 1/16, w2 = 1/4
    // WSSR = (1/16)*(4-5)^2 + (1/4)*(2-1)^2 = 1/16 + 1/4 = 5/16 = 0.3125
    const wssr = computeWssr(cObs, cPred, '1/C_obs\u00B2');
    expectFloat(wssr, 0.3125, 1e-10);
  });
});

// -- Solve output properties --

category('Math: PK Solve properties', () => {
  test('solve_01: non-empty arrays of equal length', async () => {
    const result = solvePk(DEFAULTS);
    expect(result.t.length > 0, true, 't should be non-empty');
    expect(result.conc.length > 0, true, 'conc should be non-empty');
    expect(result.t.length, result.ac.length, 't and ac should match');
    expect(result.t.length, result.ap.length, 't and ap should match');
    expect(result.t.length, result.conc.length, 't and conc should match');
  });

  test('solve_02: non-negative amounts', async () => {
    const result = solvePk(DEFAULTS);
    for (let i = 0; i < result.ac.length; i++) {
      expect(result.ac[i] >= -0.01, true, `ac[${i}] = ${result.ac[i]} is negative`);
      expect(result.ap[i] >= -0.01, true, `ap[${i}] = ${result.ap[i]} is negative`);
    }
  });

  test('solve_03: IV bolus initial concentration = Dose/Vc', async () => {
    const params = {...DEFAULTS, inputMode: 'IV Bolus' as const, dose: 100, vc: 20};
    const result = solvePk(params);
    expectFloat(result.conc[0], 5.0, 0.1); // C(0) = 100/20 = 5
  });

  test('solve_04: infusion and oral start at C(0) = 0', async () => {
    const infResult = solvePk({...DEFAULTS, inputMode: 'IV Infusion'});
    expectFloat(infResult.conc[0], 0, 0.01);

    const oralResult = solvePk({...DEFAULTS, inputMode: 'Oral'});
    expectFloat(oralResult.conc[0], 0, 0.01);
  });

  test('solve_05: C_max > 0 and t_max > 0', async () => {
    const result = solvePk(DEFAULTS);
    expect(result.cMax > 0, true, 'C_max should be positive');
    expect(result.tMax >= 0, true, 't_max should be non-negative');
  });

  test('solve_06: AUC > 0', async () => {
    const result = solvePk(DEFAULTS);
    expect(result.auc > 0, true, 'AUC should be positive');
  });

  test('solve_07: concentration decays toward zero', async () => {
    const params = {...DEFAULTS, tEnd: 200};
    const result = solvePk(params);
    const lastConc = result.conc[result.conc.length - 1];
    expect(lastConc < 0.01, true, `C(T_end) = ${lastConc} should be near zero`);
  });

  test('solve_08: repeated dosing: Css_max > Css_min', async () => {
    const params: PkParams = {...DEFAULTS, repeatedDosing: true, nDoses: 5, tau: 12, tEnd: 72};
    const result = solvePk(params);
    expect(result.cssMax !== null, true);
    expect(result.cssMin !== null, true);
    if (result.cssMax !== null && result.cssMin !== null)
      expect(result.cssMax > result.cssMin, true, 'Css_max should exceed Css_min');
  });
});

// -- MRT solver verification --

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

    const solution = mrt(odes);
    const tArr = solution[0];
    const yArr = solution[1];
    let maxError = 0;
    for (let i = 0; i < tArr.length; i++)
      maxError = Math.max(maxError, Math.abs(exact(tArr[i]) - yArr[i]));

    expect(maxError < 0.1, true, `Max error ${maxError} exceeds 0.1`);
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

    const solution = mrt(odes);
    const tArr = solution[0];
    const yArr = solution[1];
    let maxError = 0;
    for (let i = 0; i < tArr.length; i++)
      maxError = Math.max(maxError, Math.abs(exact(tArr[i]) - yArr[i]));

    expect(maxError < 0.1, true, `Max error ${maxError} exceeds 0.1`);
  });
});
