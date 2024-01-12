import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {KnownMetrics} from '../typed-metrics';
import {createMultiDimRedWorker} from './multi-dim-red-worker-creator';
import {DistanceAggregationMethod} from '../distance-matrix/types';
import {normalize} from '@datagrok-libraries/utils/src/vector-operations';
import {DimReductionMethods} from './types';


export async function getNormalizedEmbeddings(
  dataCols: Array<any[]>,
  methodName: DimReductionMethods,
  distanceMetrics: KnownMetrics[],
  weights: number[],
  distanceAggregation: DistanceAggregationMethod,
  options: any, progressFunc?: (epoch: number, epochsLength: number, embedding: number[][]) => void
): Promise<Matrix> {
  let dimensionalityReduceRes: Matrix =
          await createMultiDimRedWorker(
            dataCols, distanceMetrics, methodName, weights, distanceAggregation, options, progressFunc
          );

  dimensionalityReduceRes = dimensionalityReduceRes.map((it) => normalize(it));
  return dimensionalityReduceRes;
}

