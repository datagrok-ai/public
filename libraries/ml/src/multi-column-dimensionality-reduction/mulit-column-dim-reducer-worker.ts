import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {KnownMetrics} from '../typed-metrics/typed-metrics';
import {DimReductionMethods} from './types';
import {DistanceAggregationMethod} from '../distance-matrix/types';
import {MultiColDimReducer} from './multi-column-dim-reducer';

/**
 * Worker thread receiving data function.
 *
 * @param {any[]} columnData Samples to process.
 * @param {KnownMethods} method Embedding method.
 * @param {KnownMetrics} measure Distance metric.
 * @param {any} options Options to pass to algorithm.
 * @param {boolean} parallelDistanceWorkers Whether to use parallel distance workers.
 * @return {any} Embedding (and distance matrix where applicable).
 */
async function onMessage(columnsData: Array<any[]>, method: DimReductionMethods, metrics: KnownMetrics[],
  weights: number[], aggregationMethod: DistanceAggregationMethod,
  options: any): Promise<Matrix> {
  const reducer = new MultiColDimReducer(
    columnsData,
    method,
    metrics,
    weights,
    aggregationMethod,
    {...options, progressFunc}
  );
  return await reducer.transform(true);
}

async function progressFunc(epochNum: number, epochsLength: number, embedding: number[][]) {
  if (epochNum % 5 === 0)
    self.postMessage({epochNum, epochsLength, embedding});
}

self.onmessage = async ({data: {columnsData, method, distanceMetrics, options, weights, aggregationMethod}}) => {
  let data: {error?: any, embedding?: any};
  try {
    const embedding = await onMessage(columnsData, method, distanceMetrics, weights, aggregationMethod, options);
    data = {embedding};
  } catch (e: any) {
    data = {error: e};
  }
  self.postMessage({
    error: data.error,
    embedding: data.embedding,
  });
};
