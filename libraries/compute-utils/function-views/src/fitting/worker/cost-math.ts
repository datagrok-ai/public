// Pure cost-function math, runnable on either arm.
//
// Extracted so worker/fitting.worker.ts can compute MAD/RMSE without pulling
// in fitting-utils.ts (which imports `datagrok-api/dg`, `grok`, `ui`, CSS).
// Typed against a structural ColLike/DfLike — DG.Column / DG.DataFrame and
// LiteColumn / LiteDataFrame both satisfy it, so a single implementation is
// shared between main thread and worker.

export interface ColLike {
  readonly name: string;
  readonly length: number;
  getRawData(): ArrayLike<number>;
  readonly stats: {readonly min: number; readonly max: number};
}

export interface DfLike {
  col(name: string): ColLike | null;
  readonly name?: string;
}

// Mirror of optimizer-misc.InconsistentTables. Worker code uses this local
// class; the main-arm executor inspects err.name === 'InconsistentTables' to
// decide whether to re-throw or fold into the fails DF, matching today's
// instanceof check in optimizer.ts.
export class InconsistentTablesError extends Error {
  constructor(msg: string) {
    super(msg);
    this.name = 'InconsistentTables';
  }
}

export function getIndices(expArg: ColLike, simArg: ColLike): Uint32Array {
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
}

export function getErrors(
  expArg: ColLike | null, expFuncs: ColLike[], simDf: DfLike, toScale: boolean,
): Float32Array {
  if (expArg === null)
    throw new InconsistentTablesError('no argument column in the target output dataframe');

  const arg = expArg.name;
  const simArg = simDf.col(arg);

  if (simArg === null) {
    throw new InconsistentTablesError(
      `no "${arg}" column in the output dataframe "${simDf.name ?? ''}"`);
  }

  const indices = getIndices(expArg, simArg);
  const expColsCount = expFuncs.length;
  const errors = new Float32Array(expColsCount * indices.length);
  let errIdx = 0;

  for (let idx = 0; idx < expColsCount; ++idx) {
    const expCol = expFuncs[idx];
    const simCol = simDf.col(expCol.name);
    if (simCol === null) {
      throw new InconsistentTablesError(
        `no "${expCol.name}" column in the output dataframe "${simDf.name ?? ''}"`);
    }

    const simRaw = simCol.getRawData();
    const expRaw = expCol.getRawData();

    if (toScale) {
      const expScale = Math.max(Math.abs(expCol.stats.max), Math.abs(expCol.stats.min));
      const coef = (expScale > 0) ? expScale : 1;
      indices.forEach((simIdx, expIdx) => {
        errors[errIdx] = (simRaw[simIdx] - expRaw[expIdx]) / coef;
        ++errIdx;
      });
    } else {
      indices.forEach((simIdx, expIdx) => {
        errors[errIdx] = simRaw[simIdx] - expRaw[expIdx];
        ++errIdx;
      });
    }
  }
  return errors;
}
