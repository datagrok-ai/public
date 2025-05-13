// The MOEA/D method for multi-objective optimization: https://ieeexplore.ieee.org/document/4358754
// The method specific definitions

export type Func = (x: Float32Array) => Float32Array;

export type InputOptions = {
  dim: number,
  mins: Float32Array,
  maxs: Float32Array,
};

//n_weights=100, generations=100, neighbors=10, mutation_rate=0.2
export type MoeadOptions = {
  nWeights: number,
  generations: number,
  neighbors: number,
  mutationRate: number,
};

export enum MOEAD_DEFAULTS {
  N_WEIGHTS = 5, //100,
  GENERATIONS = 100,
  NEIGHBORS = 10,
  MUTATION_RATE = 0.2,
};

export const DEFAULT_SETTINGS: MoeadOptions = {
  nWeights: MOEAD_DEFAULTS.N_WEIGHTS,
  generations: MOEAD_DEFAULTS.GENERATIONS,
  neighbors: MOEAD_DEFAULTS.NEIGHBORS,
  mutationRate: MOEAD_DEFAULTS.MUTATION_RATE,
};
