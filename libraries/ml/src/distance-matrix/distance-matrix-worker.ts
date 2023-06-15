import {isNil} from './utils';
import {Measure} from '../typed-metrics';

onmessage = (event) => {
  const {values, fnName, startRow, startCol, chunckSize, opts} = event.data;
  const data: { error?: any, distanceMatrixData?: Float32Array, min?: number, max?: number} = {};
  try {
    let i = startRow;
    let j = startCol;
    let cnt = 0;
    let lmin = 0;
    let lmax = Number.MIN_VALUE;
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
