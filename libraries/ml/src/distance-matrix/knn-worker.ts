import {getAggregationFunction, insertSmaller, isNil} from './utils';
import {KnownMetrics, Measure, isBitArrayMetric} from '../typed-metrics';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {DistanceAggregationMethod} from './types';

onmessage = async (event) => {
  const {values, startIdx, endIdx, fnNames, opts, nNeighbours, weights, aggregationMethod}:
    {values: Array<any[]>, startIdx: number, endIdx: number, fnNames: KnownMetrics[],
      opts: any[], nNeighbours: number, weights: number[], aggregationMethod: DistanceAggregationMethod} = event.data;
  try {
    const knnIndexes = new Array(values[0].length).fill(null).map(() => new Array<number>(nNeighbours).fill(-1));
    const knnDistances = new Array(values[0].length).fill(null).map(() => new Array<number>(nNeighbours).fill(999999));
    const aggregate = getAggregationFunction(aggregationMethod, weights);
    const chunkSize = endIdx - startIdx;
    values.forEach((v, colIdx) => {
      if (isBitArrayMetric(fnNames[colIdx])) {
        for (let i = 0; i < v.length; ++i) {
          if (isNil(v[i])) continue;
          values[colIdx][i] = new BitArray(values[colIdx][i]._data, values[colIdx][i]._length);
        }
      }
    });
    let cnt = 0;
    const distanceFns = new Array(fnNames.length).fill(null).map((_, i) => new Measure(fnNames[i]).getMeasure(opts[i]));
    const startRow = values[0].length - 2 - Math.floor(
      Math.sqrt(-8 * startIdx + 4 * values[0].length * (values[0].length - 1) - 7) / 2 - 0.5);
    const startCol = startIdx - values[0].length * startRow + Math.floor((startRow + 1) * (startRow + 2) / 2);
    let mi = startRow;
    let mj = startCol;
    while (cnt < chunkSize) {
      //const value = seq1List[mi] && seq1List[mj] ? hamming(seq1List[mi], seq1List[mj]) : 0;
      const distanceValues = distanceFns.map((fn, idx) => !isNil(values[idx][mi]) && !isNil(values[idx][mj]) ?
        fn(values[idx][mi], values[idx][mj]) : 1);
      // const value = !isNil(values[mi]) && !isNil(values[mj]) ?
      //   distanceFn(values[mi], values[mj]) : 1;
      const value = distanceValues.length === 1 ? distanceValues[0] : aggregate(distanceValues);
      insertSmaller(knnDistances[mi], knnIndexes[mi], value, mj);
      insertSmaller(knnDistances[mj], knnIndexes[mj], value, mi);

      cnt++;
      mj++;
      if (mj === values[0].length) {
        mi++;
        mj = mi + 1;
      }
    }

    postMessage({knnDistances, knnIndexes});
  } catch (e) {
    postMessage({error: e});
  }
};
