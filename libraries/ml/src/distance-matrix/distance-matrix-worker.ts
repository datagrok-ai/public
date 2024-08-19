import {getAggregationFunction, isNil} from './utils';
import {KnownMetrics, Measure, isBitArrayMetric} from '../typed-metrics';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {DistanceAggregationMethod} from './types';

onmessage = (event) => {
  const {values, fnNames, startRow, startCol, chunckSize, opts, weights, aggregationMethod}: {
    values: any[][], fnNames: KnownMetrics[], startRow: number, startCol: number, chunckSize: number,
    opts: {[_: string]: any}[], weights: number[], aggregationMethod: DistanceAggregationMethod,
  } = event.data;
  const data: { error?: any, distanceMatrixData?: Float32Array, min?: number, max?: number} = {};
  try {
    let i = startRow;
    let j = startCol;
    let cnt = 0;
    let lmin = 0;
    let lmax = Number.MIN_VALUE;
    const aggregate = getAggregationFunction(aggregationMethod, weights);
    values.forEach((v, colIdx) => {
      if (isBitArrayMetric(fnNames[colIdx])) {
        for (let mi = 0; mi < v.length; ++mi) {
          if (isNil(v[mi])) continue;
          values[colIdx][mi] = new BitArray(values[colIdx][mi]._data, values[colIdx][mi]._length);
        }
      }
    });
    const distanceFns = new Array(fnNames.length)
      .fill(null)
      .map((_, di) => new Measure(fnNames[di]).getMeasure(opts[di]));
    const retVal = new Float32Array(chunckSize);
    while (cnt < chunckSize) {
      const distanceValues = distanceFns.map((fn, idx) => !isNil(values[idx][i]) && !isNil(values[idx][j]) ?
        fn(values[idx][i], values[idx][j]) : 1);
      const value = distanceValues.length === 1 ? distanceValues[0] : aggregate(distanceValues);
      retVal[cnt] = value;
      if (value < lmin)
        lmin = value;
      if (value > lmax)
        lmax = value;
      cnt++;
      j++;
      if (j === values[0].length) {
        i++;
        j = i + 1;
      }
    }
    data.distanceMatrixData = retVal;
    data.min = lmin;
    data.max = lmax;
  } catch (e) {
    data.error = e;
  }
  postMessage(data);
};
