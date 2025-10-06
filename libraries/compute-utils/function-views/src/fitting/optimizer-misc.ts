/* eslint-disable valid-jsdoc */
// Oprimizer routine

import * as DG from 'datagrok-api/dg';

export type Extremum = {
  point: Float64Array,
  cost: number,
  iterCosts: number[],
  iterCount: number,
}

export type Setting = {
  default: number,
  min: number,
  max: number,
  caption: string,
  tooltipText: string,
  inputType: string,
}

export type OptimizationResult = {
  extremums: Extremum[],
  fails: DG.DataFrame | null,
};

export type OptimizationTask = {
  costFunc: (x: Float64Array) => Promise<number>,
  minVals: Float64Array,
  maxVals: Float64Array,
  samplesCount: number
};

export interface IOptimizer {
  (objectiveFunc: (x: Float64Array) => Promise<number>,
    paramsInitial: Float64Array,
    settings: Map<string, number>,
    restrictionsBottom: Float64Array,
    restrictionsTop: Float64Array,
    threshold?: number,
  ): Promise<Extremum>;
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
