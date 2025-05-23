// Sensitivity Analysis definitions for Diff Studio models

import {IVP2WebWorker, Pipeline} from '@datagrok/diff-grok';

/** Task for sensitivity analysis of Diff Studio model */
export type AnalysisTask = {
  ivp2ww: IVP2WebWorker,
  pipeline: Pipeline,
  inputVec: Float64Array,
  rowIdx: number,
  argColIdx: number,
  argVal: number,
};

export const DEFAULT_NUM = -10;
export const MIN_WORKERS_COUNT = 1;
export const WORKERS_COUNT_DOWNSHIFT = 2;

/** Task to pass to the workers */
export type WorkerTask = {
  ivp2ww: IVP2WebWorker,
  pipelines: Pipeline[],
  inputVecs: Float64Array[],
  indices: number[],
  rowIdx: number,
  argColIdx: number,
  argVal: number,
};

export enum RESULT_CODE {
  SUCCEED = 0,
  FAILED = 1,
};

export const NO_ERRORS = 'Success';
