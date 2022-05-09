import {DimensionalityReducer, KnownMethods} from '../reduce-dimensionality';
import {KnownMetrics} from '../typed-metrics';

/**
 * Worker thread receiving data function.
 *
 * @param {any[]} columnData Samples to process.
 * @param {KnownMethods} method Embedding method.
 * @param {KnownMetrics} measure Distance metric.
 * @param {number} cyclesCount Number of cycles to repeat.
 * @return {any} Embedding (and distance matrix where applicable).
 */
function onMessage(columnData: any[], method: KnownMethods, measure: KnownMetrics, options?: any) {
  const reducer = new DimensionalityReducer(
    columnData,
    method,
    measure,
    options,
  );
  return reducer.transform(true);
}

self.onmessage = ({data: {columnData, method, measure, options}}) => {
  const embedding = onMessage(columnData, method, measure, options);
  self.postMessage({
    distance: embedding.distance,
    embedding: embedding.embedding,
  });
};
