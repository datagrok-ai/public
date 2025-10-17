// Pareto front constants & definitions

export enum OPT_TYPE {
  MIN = 'minimize',
  MAX = 'maximize',
};

export type NumericFeature = {
  toOptimize: boolean,
  optType: OPT_TYPE,
};

export type NumericArray = Float32Array | Float64Array | Int32Array | Uint32Array;
export type ParetoLabel = 'optimal' | 'non-optimal';

export const DIFFERENCE = 2;
export enum DOCK_RATIO {
  FORM = 0.2,
  VIEWER = 0.5,
};

export const OPTIMALITY_COL_NAME = 'Pareto optimality';
