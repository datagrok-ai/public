import {DistanceMatrix} from '@datagrok-libraries/bio/src/trees/distance-matrix';
import {mmDistanceFunctions, MmDistanceFunctionsNames}
  from '@datagrok-libraries/ml/src/macromolecule-distance-functions';

onmessage = (event) => {
  const {values, fnName} = event.data;
  const data: { error?: any; distanceMatrixData?: Float32Array } = {};
  try {
    const distanceMatrix = DistanceMatrix.calc(values, mmDistanceFunctions[fnName as MmDistanceFunctionsNames]());
    distanceMatrix.normalize();
    data.distanceMatrixData = distanceMatrix.data;
  } catch (e) {
    data.error = e;
  }
  postMessage(data);
};
