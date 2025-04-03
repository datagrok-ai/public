/* eslint-disable valid-jsdoc */
// Utilities for similarity computation

import * as DG from 'datagrok-api/dg';
import {Extremum, TargetTableOutput} from './optimizer-misc';

/** Similarity constants */
export enum SIMILARITY {
    MAX = 1,
    MIN = 0,
    SCALE = 100, // %
  };

/** Return the Maximum Relative Deviation (%) between vectors */
function getMRD(source: Float32Array, updated: Float32Array): number {
  let mrd = 0;

  const len = source.length;

  if (len !== updated.length)
    throw new Error(`Inconsistent length of found points: source - ${source.length}, updated - ${updated.length}`);

  let divisor: number;

  for (let i = 0; i < len; ++i) {
    divisor = source[i] !== 0 ? source[i] : 1;

    mrd = Math.max(mrd, Math.abs((source[i] - updated[i]) / divisor) * SIMILARITY.SCALE);
  }

  return mrd;
}

/** Return the Intersection over Union (IoU) of two ranges */
function getIoU(min1: number, max1: number, min2: number, max2: number): number {
  if ((min1 > max2) || (max1 < min2))
    return SIMILARITY.MIN;

  const intersectionLength = Math.min(max1, max2) - Math.max(min1, min2);
  const unionLength = Math.max(max1, max2) - Math.min(min1, min2);

  return unionLength > 0 ? intersectionLength / unionLength : SIMILARITY.MAX;
}

/** Return the similarity between argument columns */
function getArgIoU(argCol1: DG.Column, argCol2: DG.Column): number {
  if (!argCol1.isNumerical)
    throw new Error(`Non-numerical argument column: ${argCol1.name}, type - ${argCol1.type}`);

  if (!argCol2.isNumerical)
    throw new Error(`Non-numerical argument column: ${argCol2.name}, type - ${argCol2.type}`);

  return getIoU(argCol1.stats.min, argCol1.stats.max, argCol2.stats.min, argCol2.stats.max);
}

/** Return similarity of 2 functions */
function getFuncSimilarity(argCol1: DG.Column, funcCol1: DG.Column, argCol2: DG.Column, funcCol2: DG.Column,
  argMin: number, argMax: number): number {
  return 0;
}

/** Check whether df-s are similar */
function areSimilar(df1: DG.DataFrame, df2: DG.DataFrame,
  argName: string, argSimTreshold: number, funcSimTreshold: number): boolean {
  const argCol1 = df1.col(argName);
  const argCol2 = df2.col(argName);

  if (argCol1 === null)
    throw new Error(`Inconsistent dataframe: ${df1.name}, no argument column "${argName}"`);

  if (argCol2 === null)
    throw new Error(`Inconsistent dataframe: ${df2.name}, no argument column "${argName}"`);

  console.log('Arg sim: ', getArgIoU(argCol1, argCol2));

  // Check similarity of args: no reason to analyze funcs if args are different
  if (getArgIoU(argCol1, argCol2) < argSimTreshold)
    return false;

  const argMin = Math.max(argCol1.stats.min, argCol2.stats.min);
  const argMax = Math.max(argCol1.stats.max, argCol2.stats.max);

  let maxDeviation = 0;

  for (const funcCol1 of df1.columns) {
    if (funcCol1.name === argName)
      continue;

    if (!funcCol1.isNumerical)
      continue;

    const funcCol2 = df2.col(funcCol1.name);

    if (funcCol2 === null)
      throw new Error(`Inconsistent dataframe: ${df2.name}, no function column "${funcCol1.name}"`);

    if (!funcCol2.isNumerical)
      throw new Error(`Inconsistent dataframe: ${df2.name}, non-numerical function column "${funcCol2.name}"`);

    maxDeviation = Math.max(maxDeviation, getFuncSimilarity(argCol1, funcCol1, argCol2, funcCol2, argMin, argMax));
  }

  console.log('MAD:', maxDeviation);

  return maxDeviation < funcSimTreshold;
} // areSimilar

/** Return non-similar extrema */
export async function getNonSimilar(extrema: Extremum[],
  similarity: number,
  getCalledFuncCall: (x: Float32Array) => Promise<DG.FuncCall>,
  targetDfs: TargetTableOutput[],
): Promise<Extremum[]> {
  // Get computations at extrema points
  const pointsCount = extrema.length;
  const calledFuncCalls = new Array<DG.FuncCall>(pointsCount);

  for (let i = 0; i < pointsCount; ++i)
    calledFuncCalls[i] = await getCalledFuncCall(extrema[i].point);

  const nonSimilar: Extremum[] = [];
  let toPush: boolean;
  const targetDfsCount = targetDfs.length;
  const scaledSimilarity = similarity / SIMILARITY.SCALE;
  console.log('Scaled sim: ', scaledSimilarity);

  for (let extrIdx = 0; extrIdx < pointsCount; ++extrIdx) {
    const extr = extrema[extrIdx];
    const funcCall = calledFuncCalls[extrIdx];

    toPush = true;

    for (let prevIdx = 0; prevIdx < nonSimilar.length; ++prevIdx) {
      toPush &&= getMRD(nonSimilar[prevIdx].point, extr.point) > similarity;

      if (!toPush)
        break;

      for (let idx = 0; idx < targetDfsCount; ++idx) {
        toPush &&= !areSimilar(
          calledFuncCalls[prevIdx].getParamValue(targetDfs[idx].name) as DG.DataFrame,
          funcCall.getParamValue(targetDfs[idx].name) as DG.DataFrame,
          targetDfs[idx].argColName,
          SIMILARITY.MAX - scaledSimilarity,
          scaledSimilarity,
        );

        if (!toPush)
          break;
      }

      if (!toPush)
        break;
    }

    if (toPush)
      nonSimilar.push(extr);
  }

  console.log(targetDfs);

  return nonSimilar;
} // getNonSimilar
