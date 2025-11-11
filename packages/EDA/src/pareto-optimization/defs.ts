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

export const DIFFERENCE = 2;
export enum RATIO {
  FORM = 0.15,
  VIEWER = 0.5,
};

export enum COL_NAME {
  OPT = 'Pareto optimality',
  SIZE = '_size',
};

export const PC_MAX_COLS = 10;
export enum SIZE {
  OPTIMAL = 8,
  NON_OPT = 4,
};

export const AXIS_NAMES = ['xColumnName', 'yColumnName'];
export const AXIS_NAMES_3D = ['xColumnName', 'yColumnName', 'zColumnName'];
export type ColorOpt = Record<string, string | undefined>;

export const SCATTER_ROW_LIM = 5000;
export const SCATTER3D_ROW_LIM = 1000;

export const USE_PARETO_AXES = true;
export const USE_AUTO_SELECTION = true;
export const PARETO_COLOR_CODING = true;
