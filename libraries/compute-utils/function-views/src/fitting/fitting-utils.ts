/* eslint-disable valid-jsdoc */
// Fitting utilities

import * as DG from 'datagrok-api/dg';

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

export function getErrors(arg: string, expDf: DG.DataFrame, simDf: DG.DataFrame): Float32Array {
  const expArg = expDf.col(arg);
  const simArg = simDf.col(arg);

  if (expArg === null)
    throw new Error(`No "${arg}" column in the "${expDf.name} table"`);

  if (simArg === null)
    throw new Error(`No "${arg}" column in the "${simDf.name} table"`);

  const indeces = getIndeces(expArg, simArg);

  const expColumns = expDf.columns;
  const expColsCount = expColumns.length;
  const errors = new Float32Array((expColsCount - 1) * indeces.length);
  let errIdx = 0;

  for (const expCol of expColumns) {
    if (expCol.name !== arg) {
      const simCol = simDf.col(expCol.name);

      if (simCol === null) {
        throw new Error(`Inconsistent dataframes "${expDf.name}" & "${simDf.name}":
          no "${expCol.name}" column in "${expDf.name}".`);
      }

      const simRaw = simCol.getRawData();
      const expRaw = expCol.getRawData();

      indeces.forEach((simIdx, expIdx) => {
        errors[errIdx] = simRaw[simIdx] - expRaw[expIdx];
        ++errIdx;
      });
    }
  }

  return errors;
} // getErrors
