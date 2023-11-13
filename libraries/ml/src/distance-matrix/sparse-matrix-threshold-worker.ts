import {isNil} from './utils';
import {KnownMetrics, Measure, isBitArrayMetric} from '../typed-metrics';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
onmessage = async (event) => {
  const {values, startIdx, endIdx, sampleLength, fnName, opts}:
    {values: any[], startIdx: number, endIdx: number, sampleLength: number, fnName: KnownMetrics, opts: any} = event.data;
  try {
    // if (startIdx != -1)
    //   throw new Error('Error in sparse threshold worker'); // TODO: remove this line
    let distances: Float32Array = new Float32Array(sampleLength);
    const chunkSize = endIdx - startIdx;

    if (isBitArrayMetric(fnName)) {
      for (let i = 0; i < values.length; ++i) {
        if (isNil(values[i])) continue;
        values[i] = new BitArray(values[i]._data, values[i]._length);
      }
    }
    let cnt = 0;
    const increment = Math.floor(chunkSize / sampleLength);
    const distanceFn = new Measure(fnName).getMeasure(opts);
    const startRow = values.length - 2 - Math.floor(
      Math.sqrt(-8 * startIdx + 4 * values.length * (values.length - 1) - 7) / 2 - 0.5);
    const startCol = startIdx - values.length * startRow + Math.floor((startRow + 1) * (startRow + 2) / 2);
    let mi = startRow;
    let mj = startCol;
    let distanceArrayCounter = 0;
    while (cnt < chunkSize && distanceArrayCounter < sampleLength) {
      //const value = seq1List[mi] && seq1List[mj] ? hamming(seq1List[mi], seq1List[mj]) : 0;
      const value = !isNil(values[mi]) && !isNil(values[mj]) ?
        distanceFn(values[mi], values[mj]) : 1;

      distances[distanceArrayCounter] = value;
      distanceArrayCounter++;
      // const currentIncrement = Math.floor(Math.random() * increment) + 1
      cnt+=increment;
      mj+=increment;
      while (mj >= values.length && cnt < chunkSize) {
        mi++;
        mj = mi + 1 + (mj - values.length);
      }
    }

    if (distanceArrayCounter < sampleLength) {
      distances = distances.slice(0, distanceArrayCounter);
    }
    postMessage({distance: distances});
  } catch (e) {
    postMessage({error: e});
  }
};