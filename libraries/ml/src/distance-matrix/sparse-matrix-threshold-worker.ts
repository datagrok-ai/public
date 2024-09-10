import {getAggregationFunction, isNil} from './utils';
import {KnownMetrics, Measure, isBitArrayMetric} from '../typed-metrics';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {DistanceAggregationMethod} from './types';
onmessage = async (event) => {
  const {values, startIdx, endIdx, sampleLength, fnNames, opts, weights, aggregationMethod}:
    {values: Array<any[]>, startIdx: number, endIdx: number,
      sampleLength: number, fnNames: KnownMetrics[], opts: any[],
      weights: number[], aggregationMethod: DistanceAggregationMethod} = event.data;
  try {
    let distances: Float32Array = new Float32Array(sampleLength);
    const chunkSize = endIdx - startIdx;
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
    const increment = Math.floor(chunkSize / sampleLength);
    const distanceFns = new Array(fnNames.length).fill(null).map((_, i) => new Measure(fnNames[i]).getMeasure(opts[i]));
    const startRow = values[0].length - 2 - Math.floor(
      Math.sqrt(-8 * startIdx + 4 * values[0].length * (values[0].length - 1) - 7) / 2 - 0.5);
    const startCol = startIdx - values[0].length * startRow + Math.floor((startRow + 1) * (startRow + 2) / 2);
    let mi = startRow;
    let mj = startCol;
    let distanceArrayCounter = 0;
    while (cnt < chunkSize && distanceArrayCounter < sampleLength) {
      //const value = seq1List[mi] && seq1List[mj] ? hamming(seq1List[mi], seq1List[mj]) : 0;
      const distanceValues = distanceFns.map((fn, idx) => !isNil(values[idx][mi]) && !isNil(values[idx][mj]) ?
        fn(values[idx][mi], values[idx][mj]) : 1);
      const value = distanceValues.length === 1 ? distanceValues[0] : aggregate(distanceValues);
      distances[distanceArrayCounter] = value;
      distanceArrayCounter++;
      // const currentIncrement = Math.floor(Math.random() * increment) + 1
      cnt+=increment;
      mj+=increment;
      while (mj >= values[0].length && cnt < chunkSize) {
        mi++;
        mj = mi + 1 + (mj - values[0].length);
      }
    }

    if (distanceArrayCounter < sampleLength)
      distances = distances.slice(0, distanceArrayCounter);

    postMessage({distance: distances});
  } catch (e) {
    postMessage({error: e});
  }
};
