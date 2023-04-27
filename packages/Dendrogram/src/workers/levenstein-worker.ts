import { DistanceMatrix } from '@datagrok-libraries/bio/src/trees/distance-matrix';
import { distance } from 'fastest-levenshtein';


onmessage = (event) => {
  const { sequences } = event.data;
  const data: { error?: any; distanceMatrixData?: Float32Array } = {};
  try {
    const distanceMatrix = DistanceMatrix.calc(sequences, (a: string,b: string) => distance(a,b));
    data.distanceMatrixData = distanceMatrix.data;
  } catch (e) {
    data.error = e;
  }
  postMessage(data);
};
