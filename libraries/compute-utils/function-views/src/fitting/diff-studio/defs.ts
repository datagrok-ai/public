import {IVP2WebWorker} from '@datagrok/diff-studio-tools';

export const ARG_IDX = 0;
export const NO_ERRORS = 'Success';

export const DEFAULT_SET_VAL = 0;
export const MIN_TARGET_COLS_COUNT = 2;
export const ARG_INP_COUNT = 3;
export const ARG_COL_IDX = 0;

/** */
export type NelderMeadInput = {
  settingNames: string[],
  settingVals: number[],
  loss: string,
  ivp2ww: IVP2WebWorker,
  nonParamNames: string[],
  fixedInputsNames: string[],
  fixedInputsVals: number[],
  variedInputNames: string[],
  variedInpMin: Float32Array,
  variedInpMax: Float32Array,
  targetNames: string[],
  targetVals: Array<Float64Array>,
  scaleVals: Float64Array,
  samplesCount: number,
};
