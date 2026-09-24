// Lotka-Volterra Predator-Prey Model — ODE specification

import {ODEs} from 'diff-grok';

/** Parameters for the Lotka-Volterra ODE system */
export interface LotkaVolterraParams {
  alpha: number;
  beta: number;
  delta: number;
  gamma: number;
  x0: number;
  y0: number;
  T: number;
}

/** Creates the ODEs specification for the Lotka-Volterra model */
export function createLotkaVolterraODE(params: LotkaVolterraParams): ODEs {
  const {alpha, beta, delta, gamma, x0, y0, T} = params;

  return {
    name: 'LotkaVolterra',
    arg: {name: 't', start: 0, finish: T, step: 0.05},
    initial: [x0, y0],
    func: (_t: number, y: Float64Array, out: Float64Array) => {
      out[0] = alpha * y[0] - beta * y[0] * y[1];
      out[1] = delta * y[0] * y[1] - gamma * y[1];
    },
    tolerance: 1e-6,
    solutionColNames: ['x(t)', 'y(t)'],
  };
}

/** Computes the non-trivial equilibrium point (x*, y*) */
export function getEquilibrium(
  alpha: number, beta: number, delta: number, gamma: number,
): {xStar: number; yStar: number} {
  return {
    xStar: gamma / delta,
    yStar: alpha / beta,
  };
}
