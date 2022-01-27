import {DimensionalityReducer, KnownMethods} from '../reduce-dimensionality';
import {KnownMetrics} from '../typed-metrics';

/**
 * Worker thread receiving data function.
 *
 * @param {any[]} columnData Samples to process.
 * @param {KnownMethods} method Embedding method.
 * @param {KnownMetrics} measure Distance metric.
 * @param {number} cyclesCount Number of cycles to repeat.
 * @return {Coordinates} Embedding.
 */
function onMessage(columnData: any[], method: KnownMethods, measure: KnownMetrics, cyclesCount?: number) {
  const reducer = new DimensionalityReducer(
    columnData,
    method,
    measure,
    cyclesCount ? {cycles: cyclesCount} : undefined,
  );
  return reducer.transform(true);
}

self.onmessage = ({data: {columnData, method, measure, cyclesCount}}) => {
  const embedding = onMessage(columnData, method, measure, cyclesCount);
  self.postMessage({
    embedding: embedding,
  });
};
