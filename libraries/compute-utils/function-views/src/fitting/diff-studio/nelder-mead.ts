/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {IVP, IVP2WebWorker, solveIvp} from '@datagrok/diff-studio-tools';
import {LOSS} from '../constants';
import {ARG_IDX} from './constants';
import {sampleParams} from '../optimizer-sampler';
import {mad, rmse} from './utils';
import {optimizeNM} from '../optimizer-nelder-mead';
import {OptimizationResult, Extremum} from '../optimizer-misc';

const DEFAULT_SET_VAL = 0;
const MIN_TARGET_COLS_COUNT = 2;
const ARG_INP_COUNT = 3;
const ARG_COL_IDX = 0;

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
  variedStart: Float32Array,
  variedInpMin: Float32Array,
  variedInpMax: Float32Array,
  targetNames: string[],
  targetVals: Array<Float64Array>,
  scaleVals: Float64Array,
  samplesCount: number,
};

/** */
export async function fit(task: NelderMeadInput): Promise<Extremum> {
  //console.log(task);

  const ivp = task.ivp2ww;

  const inputSize = ARG_INP_COUNT + ivp.deqsCount + ivp.paramNames.length;
  const ivpInputVals = new Float64Array(inputSize);
  const ivpInputNames = task.nonParamNames.concat(ivp.paramNames);
  const funcNames = task.nonParamNames.slice(ARG_INP_COUNT - 1);

  let idx = 0;
  let cur = 0;

  for (const name of ivpInputNames) {
    cur = task.fixedInputsNames.indexOf(name);
    if (cur > -1)
      ivpInputVals[idx] = task.fixedInputsVals[cur];
    else {
      cur = task.variedInputNames.indexOf(name);
      ivpInputVals[idx] = task.variedStart[cur];
    }
    ++idx;
  }

  const inpIndex = task.variedInputNames.map((name) => ivpInputNames.indexOf(name));

  const outIndex = task.targetNames.slice(1).map((name) => {
    const idx = funcNames.indexOf(name);

    if (idx < 0)
      throw new Error(`Inconsistent target dataframe: no ${name} column in the model's output.`);

    return idx;
  });
  //console.log('Out idx-s: ', outIndex);

  const dim = task.variedStart.length;
  const targetArgVals = task.targetVals[0];
  const targetFuncVals = task.targetVals.slice(1);
  const scaleVals = task.scaleVals.slice(1);

  const metric = (task.loss == LOSS.MAD) ? mad : rmse;

  /*for (let k = 0; k < 1; ++k) {
    const vec = new Float64Array(dim);

    for (let i = 0; i < dim; ++i)
      vec[i] = task.variedInpMin[i] + Math.random() * (task.variedInpMax[i] - task.variedInpMin[i]);

    console.log(vec);

    for (let i = 0; i < dim; ++i)
      ivpInputVals[inpIndex[i]] = vec[i];

    const solution = solveIvp(ivp, ivpInputVals);

    const modelArgVals = solution[ARG_COL_IDX];
    const modelFuncVals = outIndex.map((idx) => solution[idx]);

    grok.shell.addTableView(DG.DataFrame.fromColumns(solution.map((arr, idx) => {
      return DG.Column.fromFloat64Array(`${idx}`, arr);
    })));

    console.log(modelArgVals);
    console.log(modelFuncVals);

    metric(targetArgVals, targetFuncVals, scaleVals, modelArgVals, modelFuncVals);
  }*/

  let solution: Float64Array[];

  const costFunc = async (x: Float32Array): Promise<number> => {
    for (let i = 0; i < dim; ++i)
      ivpInputVals[inpIndex[i]] = x[i];

    solution = solveIvp(ivp, ivpInputVals);

    return metric(
      targetArgVals,
      targetFuncVals,
      scaleVals,
      solution[ARG_COL_IDX],
      outIndex.map((idx) => solution[idx]),
    );
  };

  const vec = new Float32Array(dim);

  for (let i = 0; i < dim; ++i)
    vec[i] = task.variedInpMin[i] + Math.random() * (task.variedInpMax[i] - task.variedInpMin[i]);

  costFunc(vec);

  const settings = new Map<string, number>(task.settingNames.map((name, idx) => [name, task.settingVals[idx]]));
  //console.log(settings);

  const res = await optimizeNM(costFunc, task.variedStart, settings, task.variedInpMin, task.variedInpMax);

  return res;
} // fit

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
  samplesCount: number): Promise<OptimizationResult> {
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

  // Run fitting
  /* for (let i = 0; i < samplesCount; ++i) {//TODO: replace 1 with samplesCount
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
      targetNames: targetNames,
      targetVals: targetVals,
      scaleVals: scaleVals,
      samplesCount: samplesCount,
    };

    res.push(await fit(task));
  }*/

  const nThreads = Math.min(Math.max(1, navigator.hardwareConcurrency - 2), samplesCount);
  const workers = new Array(nThreads).fill(null).map((_) => new Worker(new URL('workers/basic.ts', import.meta.url)));

  const resultsArray: Extremum[] = [];

  const promises = workers.map((w, idx) => {
    return new Promise<void>((resolve, reject) => {
      w.postMessage({
        task: {
          settingNames: settingNames,
          settingVals: settingVals,
          loss: loss,
          ivp2ww: ivp2ww,
          nonParamNames: nonParamNames,
          fixedInputsNames: fixedInputsNames,
          fixedInputsVals: fixedInputsVals,
          variedInputNames: variedInputNames,
          variedStart: startingPoints[idx],
          variedInpMin: minVals,
          variedInpMax: maxVals,
          targetNames: targetNames,
          targetVals: targetVals,
          scaleVals: scaleVals,
          samplesCount: samplesCount,
        },
      });

      w.onmessage = (e: any) => {
        w.terminate();
        if (e.data.callResult === 0)
          resultsArray.push(e.data.res);
        else {
          reject(e.data.msg ?? 'error in calculation');
          return;
        }
        resolve();
      };

      w.onerror = (e) => {
        w.terminate();
        console.error(e);
        reject(e);
      };
    });
  });

  await Promise.all(promises);

  console.log(resultsArray);

  return {
    extremums: resultsArray,
    fails: null,
  };
} // getFittedParams
