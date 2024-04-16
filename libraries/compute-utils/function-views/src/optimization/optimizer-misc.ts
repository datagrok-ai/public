export type Extremum = {
    arg: Float32Array,
    cost: number,
    iterationProfile: number[]
}

export type OptimizationTask = {
    costFunc: (x: Float32Array) => Promise<number>,
    minVals: Float32Array,
    maxVals: Float32Array,
    samplesCount: number
};

export interface IOptimizer {
    (objectiveFunc: (x: Float32Array) => Promise<number>, paramsInitial: Float32Array): Promise<Extremum>;
};
