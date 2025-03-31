import {getAggregationFunction, insertSmaller, isNil} from './utils';
import {KnownMetrics, Measure, isBitArrayMetric} from '../typed-metrics';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {DistanceAggregationMethod} from './types';

onmessage = async (event) => {
  const {values, fnNames, opts, nNeighbours, weights, aggregationMethod, target, startIdx}:
    {values: Array<any[]>, fnNames: KnownMetrics[], startIdx: number,
      opts: any[], nNeighbours: number, weights: number[], aggregationMethod: DistanceAggregationMethod, target: any[]} = event.data;
  try {
    const knnIndexes = new Array<number>(nNeighbours).fill(-1);
    const knnDistances = new Array<number>(nNeighbours).fill(999999);
    const aggregate = getAggregationFunction(aggregationMethod, weights);
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

    for (let i = 0; i < values[0].length; i++) {
      const distanceValues = distanceFns.map((fn, idx) => !isNil(target[idx]) && !isNil(values[idx][i]) ?
        fn(target[idx], values[idx][i]) : 1);
      const value = distanceValues.length === 1 ? distanceValues[0] : aggregate(distanceValues);
      insertSmaller(knnDistances, knnIndexes, value, i + startIdx);
    }
    postMessage({knnDistances, knnIndexes});
  } catch (e) {
    postMessage({error: e});
  }
};
