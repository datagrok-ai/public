// Pareto front constants & definitions

export enum OPT_TYPE {
  MIN = 'minimize',
  MAX = 'maximize',
};

export enum LABEL {
  OPTIMAL = 'optimal',
  NON_OPT = 'non-optimal',
};

export type NumericFeature = {
  toOptimize: boolean,
  optType: OPT_TYPE,
};

export type NumericArray = Float32Array | Float64Array | Int32Array | Uint32Array;

export const DIFFERENCE = 3;
export enum RATIO {
  FORM = 0.2,
  VIEWER = 0.5,
};

export const OPTIMALITY_COL_NAME = 'Pareto optimality';
export const PC_MAX_COLS = 10;

export const AXIS_NAMES = ['xColumnName', 'yColumnName'];
export const AXIS_NAMES_3D = ['xColumnName', 'yColumnName', 'zColumnName'];
export type ColorOpt = Record<string, string | undefined>;
