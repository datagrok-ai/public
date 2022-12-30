import {DimensionalityReducer, KnownMethods} from '../reduce-dimensionality';
import {KnownMetrics} from '../typed-metrics';

/**
 * Worker thread receiving data function.
 *
 * @param {any[]} columnData Samples to process.
 * @param {KnownMethods} method Embedding method.
 * @param {KnownMetrics} measure Distance metric.
 * @param {any} options Options to pass to algorithm.
 * @return {any} Embedding (and distance matrix where applicable).
 */
function onMessage(columnData: any[], method: KnownMethods, measure: KnownMetrics, options?: any): {distance?: any, embedding?: any} {
  const reducer = new DimensionalityReducer(
    columnData,
    method,
    measure,
    options,
  );
  return reducer.transform(true);
}

self.onmessage = ({data: {columnData, method, measure, options}}) => {
  let data: {error?: any, distance?: any, embedding?: any};
  try{
    data = onMessage(columnData, method, measure, options);
  } catch (e: any) {
    data = {error: e};
  }
  self.postMessage({
    error: data.error,
    distance: data.distance,
    embedding: data.embedding,
  });
};
