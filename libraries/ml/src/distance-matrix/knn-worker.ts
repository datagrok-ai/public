import {insertSmaller, isNil} from './utils';
import {KnownMetrics, Measure, isBitArrayMetric} from '../typed-metrics';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

onmessage = async (event) => {
  const {values, startIdx, endIdx, fnName, opts, nNeighbours}:
    {values: any[], startIdx: number, endIdx: number, fnName: KnownMetrics, opts: any, nNeighbours: number} = event.data;
  try {
    const knnIndexes = new Array(values.length).fill(null).map(() => new Array<number>(nNeighbours).fill(-1));
    const knnDistances = new Array(values.length).fill(null).map(() => new Array<number>(nNeighbours).fill(999999));

    const chunkSize = endIdx - startIdx;

    if (isBitArrayMetric(fnName)) {
      for (let i = 0; i < values.length; ++i) {
        if (isNil(values[i])) continue;
        values[i] = new BitArray(values[i]._data, values[i]._length);
      }
    }
    let cnt = 0;
    const distanceFn = new Measure(fnName).getMeasure(opts);
    const startRow = values.length - 2 - Math.floor(
      Math.sqrt(-8 * startIdx + 4 * values.length * (values.length - 1) - 7) / 2 - 0.5);
    const startCol = startIdx - values.length * startRow + Math.floor((startRow + 1) * (startRow + 2) / 2);
    let mi = startRow;
    let mj = startCol;
    while (cnt < chunkSize) {
      //const value = seq1List[mi] && seq1List[mj] ? hamming(seq1List[mi], seq1List[mj]) : 0;
      const value = !isNil(values[mi]) && !isNil(values[mj]) ?
        distanceFn(values[mi], values[mj]) : 1; 
      insertSmaller(knnDistances[mi], knnIndexes[mi], value, mj);
      insertSmaller(knnDistances[mj], knnIndexes[mj], value, mi);

      cnt++;
      mj++;
      if (mj === values.length) {
        mi++;
        mj = mi + 1;
      }
    }

    postMessage({knnDistances, knnIndexes});
  } catch (e) {
    postMessage({error: e});
  }
};