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
