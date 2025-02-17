/* eslint-disable valid-jsdoc */
import {IVP} from '@datagrok/diff-studio-tools';

export function isWorkerApplicable(ivp: IVP): boolean {
  return (ivp.loop === null) && (ivp.updates === null) && (ivp.outputs === null);
}

/** Returns indices corresponding to the closest items */
function getIndices(targetArg: Float64Array, modelArg: Float64Array): Uint32Array {
  const modelCount = modelArg.length;
  const elemsCount = Math.min(targetArg.length, modelCount);
  const indices = new Uint32Array(elemsCount);
  let idxModel = 0;
  let difPrev = 0;
  let difCur = 0;

  for (let idxExp = 0; idxExp < elemsCount; ++idxExp) {
    while (true) {
      difPrev = Math.abs(modelArg[idxModel] - targetArg[idxExp]);
      ++idxModel;

      if (idxModel < modelCount) {
        difCur = Math.abs(modelArg[idxModel] - targetArg[idxExp]);

        if (difCur > difPrev) {
          --idxModel;
          break;
        } else
          difPrev = difCur;
      } else {
        --idxModel;
        break;
      }
    }

    indices[idxExp] = idxModel;
  }

  return indices;
}; // getIndices

export function mad(targetArg: Float64Array, targetFuncs: Float64Array[], scaleVals: Float64Array,
  modelArg: Float64Array, modelFuncs: Float64Array[]): number {
  let res = 0;
  const dim = scaleVals.length;
  let localMax = 0;

  const indices = getIndices(targetArg, modelArg);
  const argsCount = indices.length;

  for (let i = 0; i < dim; ++i) {
    localMax = 0;

    for (let k = 0; k < argsCount; ++k)
      localMax = Math.max(localMax, Math.abs(targetFuncs[i][k] - modelFuncs[i][indices[k]]));

    //console.log(`${targetArg[k]} <-> ${targetFuncs[i][k]}; ${modelArg[indices[k]]} <-> ${modelFuncs[i][indices[k]]}`);
    //console.log(`local = ${localMax}`);


    res = Math.max(res, localMax / scaleVals[i]);
  }

  //console.log(`MAD = ${res}`);

  return res;
} // mad

export function rmse(targetArg: Float64Array, targetFuncs: Float64Array[], scaleVals: Float64Array,
  modelArg: Float64Array, modelFuncs: Float64Array[]): number {
  let res = 0;
  const dim = scaleVals.length;
  let scale = 0;
  let sum = 0;

  const indices = getIndices(targetArg, modelArg);
  const argsCount = indices.length;

  for (let i = 0; i < dim; ++i) {
    scale = scaleVals[i];
    sum = 0;

    for (let k = 0; k < argsCount; ++k)
      sum += ((targetFuncs[i][k] - modelFuncs[i][indices[k]]) / scale)**2;

    //console.log(`${targetArg[k]} <-> ${targetFuncs[i][k]}; ${modelArg[indices[k]]} <-> ${modelFuncs[i][indices[k]]}`);
    //console.log(`sum = ${sum}`);


    res += sum;
  }

  //console.log(`RMSE = ${Math.sqrt(res)}`);

  return Math.sqrt(res);
} // rmse
