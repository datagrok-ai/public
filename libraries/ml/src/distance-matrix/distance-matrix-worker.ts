import {isNil} from './utils';
import {Measure, isBitArrayMetric} from '../typed-metrics';
import BitArray from '@datagrok-libraries/utils/src/bit-array';

onmessage = (event) => {
  const {values, fnName, startRow, startCol, chunckSize, opts} = event.data;
  const data: { error?: any, distanceMatrixData?: Float32Array, min?: number, max?: number} = {};
  try {
    // if (startRow != -1)
    //   throw new Error('Error in distance matrix worker'); // TODO: remove this line
    let i = startRow;
    let j = startCol;
    let cnt = 0;
    let lmin = 0;
    let lmax = Number.MIN_VALUE;
    if (isBitArrayMetric(fnName)) {
      for (let i = 0; i < values.length; ++i) {
        if (isNil(values[i])) continue;
        values[i] = new BitArray(values[i]._data, values[i]._length);
      }
    }

    const retVal = new Float32Array(chunckSize);
    const distanceFn = new Measure(fnName).getMeasure(opts);
    while (cnt < chunckSize) {
      const value = !isNil(values[i]) && !isNil(values[j]) ?
        distanceFn(values[i], values[j]) : 1;
      retVal[cnt] = value;
      if (value < lmin)
        lmin = value;
      if (value > lmax)
        lmax = value;
      cnt++;
      j++;
      if (j === values.length) {
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
