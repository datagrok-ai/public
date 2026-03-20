// The MOEA/D method for multi-objective optimization: https://ieeexplore.ieee.org/document/4358754
// The method specific definitions

export type Func = (x: Float32Array) => Promise<Float32Array>;

export type InputOptions = {
  dim: number,
  mins: Float32Array,
  maxs: Float32Array,
};

export type MoeadOptions = {
  nWeights: number,
  generations: number,
  neighbors: number,
  mutationRate: number,
};

export type OptResult = Promise<Float32Array[]>;

export enum MOEAD_DEFAULTS {
  N_WEIGHTS = 99,
  GENERATIONS = 10,
  NEIGHBORS = 10,
  MUTATION_RATE = 0.2,
};

export const DEFAULT_SETTINGS: MoeadOptions = {
  nWeights: MOEAD_DEFAULTS.N_WEIGHTS,
  generations: MOEAD_DEFAULTS.GENERATIONS,
  neighbors: MOEAD_DEFAULTS.NEIGHBORS,
  mutationRate: MOEAD_DEFAULTS.MUTATION_RATE,
};

export type Validation = {
  res: boolean,
  msg: string,
};

export const SLEEP_TIME = 5;

export enum OPT_TYPE {
  MIN = 'Minimize',
  MAX = 'Maximize',
};

export type TooltipInfo = {
  msg: string,
  isInput: boolean,
};
