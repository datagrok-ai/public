import {isNil} from '@datagrok-libraries/ml/src/distance-matrix';
import {mmDistanceFunctions} from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

onmessage = (event) => {
  const {fnName, values, templateIdx, start, end} = event.data;
  const data: { error?: any, distanceArrayData?: Float32Array, min?: number, max?: number} = {};
  try {
    let lmin = 0;
    let lmax = Number.MIN_VALUE;
    const retVal = new Float32Array(end - start);
    const distanceFn = mmDistanceFunctions[fnName as keyof typeof mmDistanceFunctions]();

    for (let i = start; i < end; i++) {
      const value = !isNil(values[i]) && !isNil(values[templateIdx]) ?
        distanceFn(values[i], values[templateIdx]) : 1;
      retVal[i - start] = value;
      if (value < lmin)
        lmin = value;
      if (value > lmax)
        lmax = value;
    }
    data.distanceArrayData = retVal;
    data.min = lmin;
    data.max = lmax;
  } catch (e) {
    data.error = e;
  }
  postMessage(data);
};
