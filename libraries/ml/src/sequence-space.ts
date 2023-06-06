import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {Vector} from '@datagrok-libraries/utils/src/type-declarations';
import {normalize} from '@datagrok-libraries/utils/src/vector-operations';
import * as DG from 'datagrok-api/dg';
import {BitArrayMetrics, StringMetrics, ValidTypes, VectorMetrics} from './typed-metrics/typed-metrics';
import {createDimensinalityReducingWorker} from './workers/dimensionality-reducing-worker-creator';
import {MmDistanceFunctionsNames} from './macromolecule-distance-functions';
import {IReduceDimensionalityResult} from './reduce-dimensionality';

export async function reduceDimensinalityWithNormalization(
  dataCol: BitArray[]|Vector[]|string[],
  methodName: string,
  similarityMetric: BitArrayMetrics | VectorMetrics | StringMetrics | MmDistanceFunctionsNames,
  options?: any, parallelDistanceWorkers?: boolean): Promise<IReduceDimensionalityResult> {
  const dimensionalityReduceRes: IReduceDimensionalityResult =
        await createDimensinalityReducingWorker(
            {data: dataCol, metric: similarityMetric} as ValidTypes,
            methodName, options, parallelDistanceWorkers);

  dimensionalityReduceRes.embedding = dimensionalityReduceRes.embedding.map((it) => normalize(it));
  return dimensionalityReduceRes;
}
