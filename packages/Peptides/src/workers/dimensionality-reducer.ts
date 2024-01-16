import {DimensionalityReducer, KnownMethods} from '@datagrok-libraries/ml/src/reduce-dimensionality';
import {KnownMetrics} from '@datagrok-libraries/ml/src/typed-metrics';

/**
 * Worker thread receiving data function.
 * @param columnData Samples to process.
 * @param method Embedding method.
 * @param measure Distance metric.
 * @param options Options to pass to algorithm.
 * @return Embedding (and distance matrix where applicable).
 */
function onMessage(columnData: any[], method: KnownMethods, measure: KnownMetrics, options: any): any {
  const reducer = new DimensionalityReducer(columnData, method, measure, options);
  return reducer.transform(true);
}

self.onmessage = ({data: {columnData, method, measure, options}}): void => {
  const embedding = onMessage(columnData, method, measure, options);
  self.postMessage({
    distance: embedding.distance,
    embedding: embedding.embedding,
  });
};

