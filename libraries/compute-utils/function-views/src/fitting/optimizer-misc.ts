// Oprimizer routine

export type Extremum = {
  point: Float32Array,
  cost: number,
  iterCosts: number[],
  iterCount: number,
}

export type OptimizationResult = {
  extremums: Extremum[],
  warnings: string[],
};

export type OptimizationTask = {
  costFunc: (x: Float32Array) => Promise<number>,
  minVals: Float32Array,
  maxVals: Float32Array,
  samplesCount: number
};

export interface IOptimizer {
  (objectiveFunc: (x: Float32Array) => Promise<number>, paramsInitial: Float32Array, settings: any):
  Promise<Extremum>;
};
