/* eslint-disable valid-jsdoc */
// Fitting utilities

import * as DG from 'datagrok-api/dg';

import {InconsistentTables} from './optimizer-misc';

/** Returns indeces corresponding to the closest items */
export function getIndeces(expArg: DG.Column, simArg: DG.Column): Uint32Array {
  const expArgRaw = expArg.getRawData();
  const simArgRaw = simArg.getRawData();
  const simCount = simArg.length;
  const elemsCount = Math.min(expArg.length, simCount);
  const indeces = new Uint32Array(elemsCount);
  let idxSim = 0;
  let difPrev = 0;
  let difCur = 0;

  for (let idxExp = 0; idxExp < elemsCount; ++idxExp) {
    while (true) {
      difPrev = Math.abs(simArgRaw[idxSim] - expArgRaw[idxExp]);
      ++idxSim;

      if (idxSim < simCount) {
        difCur = Math.abs(simArgRaw[idxSim] - expArgRaw[idxExp]);

        if (difCur > difPrev) {
          --idxSim;
          break;
        } else
          difPrev = difCur;
      } else {
        --idxSim;
        break;
      }
    }

    indeces[idxExp] = idxSim;
  }

  return indeces;
};

export function getErrors(arg: string, expDf: DG.DataFrame, simDf: DG.DataFrame, toScale: boolean): Float32Array {
  const expArg = expDf.col(arg);
  const simArg = simDf.col(arg);

  if (expArg === null)
    throw new InconsistentTables(`no "${arg}" column in the target dataframe "${expDf.name}"`);

  if (simArg === null)
    throw new InconsistentTables(`no "${arg}" column in the output dataframe "${simDf.name}"`);

  const indeces = getIndeces(expArg, simArg);

  const expColumns = expDf.columns;
  const expColsCount = expColumns.length;
  const errors = new Float32Array((expColsCount - 1) * indeces.length);
  let errIdx = 0;

  for (const expCol of expColumns) {
    if (expCol.name !== arg) {
      const simCol = simDf.col(expCol.name);

      if (simCol === null)
        throw new InconsistentTables(`no "${expCol.name}" column in the output dataframe "${simDf.name}"`);

      const simRaw = simCol.getRawData();
      const expRaw = expCol.getRawData();

      if (toScale) {
        const expScale = expCol.stats.max - expCol.stats.min;
        const coef = (expScale > 0) ? expScale : 1;

        indeces.forEach((simIdx, expIdx) => {
          errors[errIdx] = (simRaw[simIdx] - expRaw[expIdx]) / coef;
          ++errIdx;
        });
      } else {
        indeces.forEach((simIdx, expIdx) => {
          errors[errIdx] = simRaw[simIdx] - expRaw[expIdx];
          ++errIdx;
        });
      }
    }
  }

  return errors;
} // getErrors
