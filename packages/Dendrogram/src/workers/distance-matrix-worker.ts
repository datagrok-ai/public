import {DistanceMatrix} from '@datagrok-libraries/bio/src/trees/distance-matrix';
import {distanceFunctions} from './distance-functions';

onmessage = (event) => {
  const {values, fnName} = event.data;
  const data: { error?: any; distanceMatrixData?: Float32Array } = {};
  try {
    const distanceMatrix = DistanceMatrix.calc(values, distanceFunctions[fnName as keyof typeof distanceFunctions]);
    data.distanceMatrixData = distanceMatrix.data;
  } catch (e) {
    data.error = e;
  }
  postMessage(data);
};
