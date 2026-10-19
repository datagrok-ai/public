// Levins Metapopulation Model — ODE specification

import {ODEs} from 'diff-grok';

/** Parameters for the Levins metapopulation ODE */
export interface LevinsParams {
  p0: number;
  m: number;
  e0: number;
  rescueEffect: boolean;
  t_start: number;
  t_end: number;
  t_step: number;
  tolerance: number;
}

/** Creates the ODEs specification for the Levins model, usable in both main thread and workers */
export function createLevinsODE(params: LevinsParams): ODEs {
  const {p0, m, e0, rescueEffect, t_start, t_end, t_step, tolerance} = params;

  return {
    name: 'Levins',
    arg: {name: 't', start: t_start, finish: t_end, step: t_step},
    initial: [p0],
    func: (_t: number, y: Float64Array, out: Float64Array) => {
      const e = rescueEffect ? e0 * (1 - y[0]) : e0;
      out[0] = m * y[0] * (1 - y[0]) - e * y[0];
    },
    tolerance: tolerance,
    solutionColNames: ['p(t)'],
  };
}

/** Computes the analytical equilibrium p* for the base Levins model */
export function getEquilibrium(m: number, e0: number, rescueEffect: boolean): number {
  return rescueEffect ? NaN : Math.max(0, 1 - e0 / m);
}
