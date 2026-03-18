// SIR/SEIR Epidemic Model — ODE specification

import {ODEs} from 'diff-grok';

/** Parameters for the SIR/SEIR ODE system */
export interface EpidemicParams {
  modelType: 'SIR' | 'SEIR';
  r0: number;
  gamma: number;
  sigma: number;
  vaccination: number;
}

/** Solution of the SIR/SEIR ODE system */
export interface EpidemicSolution {
  t: Float64Array;
  S: Float64Array;
  E: Float64Array | null;
  I: Float64Array;
  R: Float64Array;
  rEff: Float64Array;
  peakDay: number;
  peakCount: number;
  finalRecoveredPct: number;
  herdImmunityThreshold: number;
}

const N = 10000;
const I0 = 1;

/** Creates the ODEs specification for the SIR model */
export function createSIR_ODE(params: EpidemicParams): ODEs {
  const beta = params.r0 * params.gamma;
  const v = params.vaccination / 100;
  const S0 = N * (1 - v) - I0;
  const R0_init = N * v;

  return {
    name: 'SIR',
    arg: {name: 't', start: 0, finish: 300, step: 0.1},
    initial: [S0, I0, R0_init],
    func: (_t: number, y: Float64Array, out: Float64Array) => {
      out[0] = -beta * y[0] * y[1] / N;                  // dS/dt
      out[1] = beta * y[0] * y[1] / N - params.gamma * y[1]; // dI/dt
      out[2] = params.gamma * y[1];                        // dR/dt
    },
    tolerance: 1e-6,
    solutionColNames: ['S', 'I', 'R'],
  };
}

/** Creates the ODEs specification for the SEIR model */
export function createSEIR_ODE(params: EpidemicParams): ODEs {
  const beta = params.r0 * params.gamma;
  const v = params.vaccination / 100;
  const S0 = N * (1 - v) - I0;
  const E0 = 0;
  const R0_init = N * v;

  return {
    name: 'SEIR',
    arg: {name: 't', start: 0, finish: 300, step: 0.1},
    initial: [S0, E0, I0, R0_init],
    func: (_t: number, y: Float64Array, out: Float64Array) => {
      out[0] = -beta * y[0] * y[2] / N;                        // dS/dt
      out[1] = beta * y[0] * y[2] / N - params.sigma * y[1];   // dE/dt
      out[2] = params.sigma * y[1] - params.gamma * y[2];       // dI/dt
      out[3] = params.gamma * y[2];                              // dR/dt
    },
    tolerance: 1e-6,
    solutionColNames: ['S', 'E', 'I', 'R'],
  };
}

/** Total population constant */
export const POPULATION = N;

/** Initial infected count */
export const INITIAL_INFECTED = I0;
