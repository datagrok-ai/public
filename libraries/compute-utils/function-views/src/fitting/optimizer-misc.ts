/* eslint-disable valid-jsdoc */
// Oprimizer routine

import * as DG from 'datagrok-api/dg';

export type Extremum = {
  point: Float32Array,
  cost: number,
  iterCosts: number[],
  iterCount: number,
}

export type Setting = {
  default: number,
  min: number,
  max: number,
}

export type OptimizationResult = {
  extremums: Extremum[],
  fails: DG.DataFrame | null,
};

export type OptimizationTask = {
  costFunc: (x: Float32Array) => Promise<number>,
  minVals: Float32Array,
  maxVals: Float32Array,
  samplesCount: number
};

export interface IOptimizer {
  (objectiveFunc: (x: Float32Array) => Promise<number>,
    paramsInitial: Float32Array,
    settings: Map<string, number>,
    restrictionsBottom: Float32Array,
    restrictionsTop: Float32Array): Promise<Extremum>;
};

/** Inconsistent tables error */
export class InconsistentTables extends Error {
  constructor(msg: string) {
    super(msg);
  }
}

/** Sleep function */
export function sleep(ms: number) {
  return new Promise((resolve, reject) => setTimeout(resolve, ms));
};

/** Target dataframe */
export type TargetTableOutput = {
  name: string,
  target: DG.DataFrame,
  argColName: string,
};
