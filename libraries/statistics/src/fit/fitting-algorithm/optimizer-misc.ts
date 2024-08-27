// Optimizer routine
import * as DG from 'datagrok-api/dg';


export type Extremum = {
  point: Float32Array,
  cost: number,
  iterCosts: number[],
  iterCount: number,
}

export type OptimizationResult = {
  extremums: Extremum[],
  fails: DG.DataFrame | null,
};

export interface IOptimizer {
  (objectiveFunc: (x: Float32Array) => {likelihood: number, residuals: number[]},
  paramsInitial: Float32Array,
  settings: any,
  restrictionsBottom: Float32Array,
  restrictionsTop: Float32Array): Extremum;
}

/** Inconsistent tables error */
export class InconsistentTables extends Error {
  constructor(msg: string) {
    super(msg);
  }
}
