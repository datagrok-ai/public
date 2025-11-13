// Pareto front computations

import {NumericArray, OPT_TYPE} from './defs';

export function getParetoMask(rawData: NumericArray[], sense: OPT_TYPE[], nPoints: number,
  nullIndices?: Set<number>): boolean[] {
  if (nPoints === 0)
    return [];

  const nDims = rawData.length;
  if (sense.length !== nDims)
    throw new Error('Sense array length must match number of dimensions');

  const pointIndices = new Uint32Array(nPoints);
  for (let i = 0; i < nPoints; i++)
    pointIndices[i] = i;

  pointIndices.sort((i1: number, i2: number) => {
    return sense[0] === OPT_TYPE.MIN ? rawData[0][i1] - rawData[0][i2] : rawData[0][i2] - rawData[0][i1];
  });

  const mask: boolean[] = Array(nPoints).fill(true);
  const paretoFrontIndices: number[] = [];

  // Set missing values to non-optimal
  nullIndices?.forEach((idx) => mask[idx] = false);

  for (const index of pointIndices) {
    if (!mask[index])
      continue;

    let dominated = false;

    for (const frontPointIndex of paretoFrontIndices) {
      let dominates = true;
      let strictlyBetter = false;

      for (let d = 0; d < nDims; d++) {
        const a = rawData[d][frontPointIndex];
        const b = rawData[d][index];
        const s = sense[d];

        if (s === OPT_TYPE.MIN) {
          if (a > b) dominates = false;
          if (a < b) strictlyBetter = true;
        } else {
          if (a < b) dominates = false;
          if (a > b) strictlyBetter = true;
        }
      } // for d

      if (dominates && strictlyBetter) {
        dominated = true;
        break;
      }
    } // for frontPointIndex

    if (dominated)
      mask[index] = false;
    else
      paretoFrontIndices.push(index);
  } // for index

  return mask;
} // getParetoMask
