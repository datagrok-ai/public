// Sensitivity analysis routine for Diff Studio models

import * as DGL from '@datagrok/diff-grok';
import {AnalysisTask, DEFAULT_NUM} from './defs';
import {LAST_ROW_IDX, FIRST_ROW_IDX} from '../constants';

function getSolutionRow(solution: Float64Array[], rowIdx: number): Float64Array {
  const outputsCount = solution.length;
  const result = new Float64Array(outputsCount);

  for (let i = 0; i < outputsCount; ++i)
    result[i] = solution[i][rowIdx];

  return result;
}

function getIndexOfClosestValue(args: Float64Array, targetVal: number): number {
  let targetIdx = 0;
  let deviation = Math.abs(args[targetIdx] - targetVal);
  let current: number;

  for (let curIdx = 1; curIdx < args.length; ++curIdx) {
    current = Math.abs(args[curIdx] - targetVal);

    if (current < deviation) {
      deviation = current;
      targetIdx = curIdx;
    }
  }

  return targetIdx;
}

export function getAnalysis(task: AnalysisTask): Float64Array {
  // evaluate model
  const solution = DGL.applyPipeline(
    task.pipeline,
    task.ivp2ww,
    task.inputVec,
  );

  const solutionTableRowCount = solution[0].length;
  const rowIdx = task.rowIdx;

  // extract target row
  if (rowIdx !== DEFAULT_NUM) {
    switch (rowIdx) {
    case FIRST_ROW_IDX:
      return getSolutionRow(solution, 0);

    case LAST_ROW_IDX:
      return getSolutionRow(solution, solutionTableRowCount - 1);

    default:
      throw new Error(`Incorrect row idx: ${rowIdx}`);
    }
  }

  // extract row by the target argument
  const index = getIndexOfClosestValue(solution[task.argColIdx], task.argVal);

  return getSolutionRow(solution, index);
}
