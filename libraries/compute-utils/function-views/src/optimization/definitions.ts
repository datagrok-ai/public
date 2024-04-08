// Optimization structures definitions

export type Extremum = {
  arg: Float32Array,
  cost: number,
}

export type OptimizationTask = {
  costFunc: (x: Float32Array) => Promise<number>,
  minVals: Float32Array,
  maxVals: Float32Array,
  samplesCount: number,
};
