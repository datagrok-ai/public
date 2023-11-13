import {isNil} from './utils';
import {KnownMetrics, Measure, isBitArrayMetric} from '../typed-metrics';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
onmessage = async (event) => {
  const {values, startIdx, endIdx, threshold, fnName, opts}:
    {values: any[], startIdx: number, endIdx: number, threshold: number, fnName: KnownMetrics, opts: any} = event.data;
  try {
    // if (startIdx != -1)
    //   throw new Error('Error in sparse matrix worker'); // TODO: remove this line
    const i: number[] = [];
    const j: number[] = [];
    const distances: number[] = [];
    const chunkSize = endIdx - startIdx;
    //const mi = startRow;
    //const mj = startCol;

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
      if (1 - value >= threshold) {
        i.push(mi);
        j.push(mj);
        distances.push(value);
      }
      cnt++;
      mj++;
      if (mj === values.length) {
        mi++;
        mj = mi + 1;
      }
    }

    const iArray = new Int32Array(i);
    const jArray = new Int32Array(j);
    const distanceArray = new Float32Array(distances);
    postMessage({i: iArray, j: jArray, distance: distanceArray});
  } catch (e) {
    postMessage({error: e});
  }
};