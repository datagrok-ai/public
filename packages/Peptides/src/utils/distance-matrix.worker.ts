import {DistanceMatrix} from '@datagrok-libraries/bio/src/trees/distance-matrix';
import {stringDistanceMetricsMethods} from '@datagrok-libraries/ml/src/typed-metrics';

const ctx: Worker = self as any;

ctx.addEventListener('message', async ({data: {peptidesList, metric}}) => {
  const data: {error?: any, distances?: Float32Array} = {};
  try {
    const dm = DistanceMatrix.calc(peptidesList,
      (a: string, b: string) => stringDistanceMetricsMethods[metric](a, b));
    data.distances = dm.data;
  } catch (e) {
    data.error = e;
  }
  self.postMessage(data);
});
