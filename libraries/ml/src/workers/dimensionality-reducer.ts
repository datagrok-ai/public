import {DimensionalityReducer, KnownMethods} from '../reduce-dimensionality';
import {Coordinates} from '@datagrok-libraries/utils/src/type-declarations';

/**
 * Worker thread receiving data function.
 *
 * @param {any[]} columnData Samples to process.
 * @param {string} method Embedding method.
 * @param {string} measure Distance metric.
 * @param {number} cyclesCount Number of cycles to repeat.
 * @return {Coordinates} Embedding.
 */
function onMessage(columnData: [], method: KnownMethods, measure?: string, cyclesCount?: number): Coordinates {
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
