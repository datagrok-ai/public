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
  (objectiveFunc: (x: Float64Array) => Promise<number|undefined>,
    paramsInitial: Float64Array,
    settings: Map<string, number>,
    threshold?: number,
  ): Promise<Extremum>;
};

/** JS package level API (should expose some in the next compute api) */

interface BoundCommon {
  name: string,
}

export interface BoundValue extends BoundCommon {
  type: 'value',
  value: number,
}

export interface BoundFormula extends BoundCommon {
  type: 'formula',
  formula: string,
}

export type BoundData = BoundValue | BoundFormula;

export interface ConstValue {
  type: 'const',
  value: any,
}

export interface ChangingValue {
  type: 'changing',
  top: BoundData,
  bottom: BoundData
}

export type ValueBoundsData = ConstValue | ChangingValue;
export type OptimizerInputsConfig = Record<string, ValueBoundsData>;

interface OutputTargetCommon {
  propName: string
}

export interface OutputTargetScalar extends OutputTargetCommon {
  type: DG.TYPE.INT | DG.TYPE.BIG_INT | DG.TYPE.FLOAT,
  target: number,
}

export interface OutputTargetDataFrame extends OutputTargetCommon {
  type: DG.TYPE.DATA_FRAME,
  target: DG.DataFrame,
  argName: string,
  cols: DG.Column[],
}

export type OutputTargetItem = OutputTargetScalar | OutputTargetDataFrame;

export type OptimizerOutputsConfig = OutputTargetItem[];


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


export async function throttle(
  desiredWorkIntervalMs: number, sleepIntervalMs: number, lastWorkStartTs: number = performance.now()
): Promise<number> {
  const now = performance.now();
  if (now - lastWorkStartTs >= desiredWorkIntervalMs) {
    await sleep(sleepIntervalMs);
    return performance.now();
  } else {
    return lastWorkStartTs;
  }
}

/** Target dataframe */
export type TargetTableOutput = {
  name: string,
  target: DG.DataFrame,
  argColName: string,
};
