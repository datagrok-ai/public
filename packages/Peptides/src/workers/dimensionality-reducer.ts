import {DimensionalityReducer, KnownMethods} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {KnownMetrics} from '@datagrok-libraries/ml/src/typed-metrics';

/**
 * Worker thread receiving data function.
 *
 * @param {any[]} columnData Samples to process.
 * @param {string} method Embedding method.
 * @param {string} measure Distance metric.
 * @param {any} options Options to pass to algorithm.
 * @return {any} Embedding (and distance matrix where applicable).
 */
function onMessage(columnData: any[], method: KnownMethods, measure: KnownMetrics, options: any): any {
  const reducer = new DimensionalityReducer(
    columnData,
    method,
    measure,
    options,
  );
  return reducer.transform(true);
}

self.onmessage = ({data: {columnData, method, measure, options}}): void => {
  const embedding = onMessage(columnData, method, measure, options);
  self.postMessage({
    distance: embedding.distance,
    embedding: embedding.embedding,
  });
};

