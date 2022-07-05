import BitArray from '@datagrok-libraries/utils/src/bit-array';
import { Matrix, Vector } from '@datagrok-libraries/utils/src/type-declarations';
import { normalize } from '@datagrok-libraries/utils/src/vector-operations';
import * as DG from 'datagrok-api/dg';
import { BitArrayMetrics, StringMetrics, ValidTypes, VectorMetrics } from './typed-metrics';
import { createDimensinalityReducingWorker, IReduceDimensionalityResult } from './workers/dimensionality-reducing-worker-creator';


export async function reduceDimensinalityWithNormalization(
    dataCol: BitArray[]|Vector[]|string[], 
    methodName: string, 
    similarityMetric: BitArrayMetrics|VectorMetrics|StringMetrics,
    options?: any): Promise<IReduceDimensionalityResult> {

    const dimensionalityReduceRes: IReduceDimensionalityResult =
        await createDimensinalityReducingWorker(
            { data: dataCol, metric: similarityMetric } as ValidTypes,
            methodName, options);

    dimensionalityReduceRes.embedding = dimensionalityReduceRes.embedding.map(it => normalize(it));
    return dimensionalityReduceRes;
}