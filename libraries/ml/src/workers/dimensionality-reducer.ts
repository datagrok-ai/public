import {DimensionalityReducer, KnownMethods} from '../reduce-dimensionality';
import {KnownMetrics} from '../typed-metrics/typed-metrics';

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
async function onMessage(columnData: any[], method: KnownMethods, measure: KnownMetrics,
  options?: any, parallelDistanceWorkers?: boolean): Promise<{distance?: any, embedding?: any}> {
  const reducer = new DimensionalityReducer(
    columnData,
    method,
    measure,
    options,
  );
  return await reducer.transform(true, parallelDistanceWorkers);
}

self.onmessage = async ({data: {columnData, method, measure, options, parallelDistanceWorkers}}) => {
  let data: {error?: any, distance?: any, embedding?: any};
  try {
    data = await onMessage(columnData, method, measure, options, parallelDistanceWorkers);
  } catch (e: any) {
    data = {error: e};
  }
  self.postMessage({
    error: data.error,
    distance: data.distance,
    embedding: data.embedding,
  });
};
