/* eslint-disable valid-jsdoc */
import * as DG from 'datagrok-api/dg';

import {IVP, IVP2WebWorker} from '@datagrok/diff-studio-tools';
import {LOSS} from '../constants';
import {ARG_IDX} from './constants';
import {sampleParams} from '../optimizer-sampler';

const DEFAULT_SET_VAL = 0;
const MIN_TARGET_COLS_COUNT = 2;

/** */
type NelderMeadInput = {
  settingNames: string[],
  settingVals: number[],
  loss: string,
  ivp2ww: IVP2WebWorker,
  nonParamNames: string[],
  fixedInputsNames: string[],
  fixedInputsVals: number[],
  variedInputNames: string[],
  variedStart: Float32Array,
  variedInpMin: Float32Array,
  variedInpMax: Float32Array,
  targetRowCount: number,
  targetNames: string[],
  targetVals: Array<Float64Array>,
  scaleVals: Float64Array,
  samplesCount: number,
};

/** */
async function fit(task: NelderMeadInput): Promise<Float32Array> {
  const res = new Float32Array(1);

  console.log(task);

  return res;
}

/** Return fitted params of Diff Studio model using the Nelder-Mead method */
export async function getFittedParams(
  loss: LOSS,
  ivp: IVP,
  ivp2ww: IVP2WebWorker,
  settings: Map<string, number>,
  variedInputNames: string[],
  minVals: Float32Array,
  maxVals: Float32Array,
  fixedInputs: Record<string, number>,
  argColName: string,
  target: DG.DataFrame,
  samplesCount: number): Promise<Float32Array[]> {
  // Extract settings names & values
  const settingNames: string[] = [...settings.keys()];
  const settingVals = settingNames.map((name) => settings.get(name) ?? DEFAULT_SET_VAL);

  // Extract fixed inputs names & vals
  const fixedInputsNames: string[] = Object.keys(fixedInputs);
  const fixedInputsVals: number[] = Object.values(fixedInputs);

  // Extract target data
  const cols = target.columns;
  const colsCount = cols.length;
  const targetRowCount = target.rowCount;

  if (colsCount < MIN_TARGET_COLS_COUNT)
    throw new Error(`Not enough of target columns: ${colsCount}. Minimum: ${MIN_TARGET_COLS_COUNT}`);

  const targetNames = new Array<string>(cols.length);
  targetNames[ARG_IDX] = argColName;

  let idx = ARG_IDX + 1;
  for (const name of cols.names()) {
    if (name !== argColName) {
      targetNames[idx] = name;
      ++idx;
    }
  }

  const scaleVals = new Float64Array(colsCount);
  const targetVals = targetNames.map((name, idx) => {
    const col = cols.byName(name);
    const absMax = Math.max(Math.abs(col.stats.max), Math.abs(col.stats.min));
    scaleVals[idx] = absMax > 0 ? absMax : 1;

    const res = new Float64Array(targetRowCount);
    const raw = col.getRawData();

    for (let i = 0; i < targetRowCount; ++i)
      res[i] = raw[i];

    return res;
  });

  // Extract non-parameters names
  const t0 = `_${ivp.arg.name}0`;
  const t1 = `_${ivp.arg.name}1`;
  const h = `_h`;
  const nonParamNames = [t0, t1, h].concat(ivp.deqs.solutionNames);

  // Generate starting points
  const startingPoints = sampleParams(samplesCount, minVals, maxVals);

  const res: Float32Array[] = [];

  // Run fitting
  for (let i = 0; i < 1; ++i) {//TODO: replace 1 with samplesCount
    // Inputs for Nelder-Mead method
    const task: NelderMeadInput = {
      settingNames: settingNames,
      settingVals: settingVals,
      loss: loss,
      ivp2ww: ivp2ww,
      nonParamNames: nonParamNames,
      fixedInputsNames: fixedInputsNames,
      fixedInputsVals: fixedInputsVals,
      variedInputNames: variedInputNames,
      variedStart: startingPoints[i],
      variedInpMin: minVals,
      variedInpMax: maxVals,
      targetRowCount: targetRowCount,
      targetNames: targetNames,
      targetVals: targetVals,
      scaleVals: scaleVals,
      samplesCount: samplesCount,
    };

    res.push(await fit(task));
  }

  return res;
} // getFittedParams
