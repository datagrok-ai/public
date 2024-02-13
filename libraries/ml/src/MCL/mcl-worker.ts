import {SparseMatrixService} from '../distance-matrix/sparse-matrix-service';
import {DistanceAggregationMethod} from '../distance-matrix/types';
import {KnownMetrics} from '../typed-metrics';
import {MCLSparseReducer} from './marcov-cluster';

onmessage = async (event) => {
  const {data, threshold, weights, aggregationMethod, distanceFnArgs, distanceFns, maxIterations}:
   {data: any[][], threshold: number,
    weights: number[], aggregationMethod: DistanceAggregationMethod,
    distanceFns: KnownMetrics[], distanceFnArgs: any[], maxIterations: number} = event.data;
  const sparse = await new SparseMatrixService()
    .calcMultiColumn(data, distanceFns, threshold / 100, distanceFnArgs, weights, aggregationMethod);
  const res = await new MCLSparseReducer({maxIterations: maxIterations ?? 5}).transform(sparse, data[0].length);
  postMessage({res});
};
