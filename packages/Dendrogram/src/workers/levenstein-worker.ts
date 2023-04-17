import { DistanceMatrix } from '@datagrok-libraries/bio/src/trees/distance-matrix';
import { distance } from 'fastest-levenshtein';


onmessage = (event) => {
  const { sequences } = event.data;

  const data: { error?: any; distanceMatrix?: DistanceMatrix } = {};
  try {
    const distanceMatrix = DistanceMatrix.calc(sequences, distance);
    data.distanceMatrix = distanceMatrix;
  } catch (e) {
    data.error = e;
  }
  postMessage(data);
};