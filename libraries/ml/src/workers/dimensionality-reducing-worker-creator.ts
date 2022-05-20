import {ValidTypes} from '../typed-metrics';

/**
 * A worker to perform dimensionality reduction.
 *
 * @param {ValidTypes} dataMetric The data to process.
 * @param {string} method A method of dimensionality reduction.
 * @param options
 * @param returnDistanceMatrix
 * @return {Promise<unknown>} Resulting embedding.
 */
export function createDimensinalityReducingWorker(dataMetric: ValidTypes, method: string,
      options?: any, returnDistanceMatrix?: boolean): Promise<unknown> {

  return new Promise(function(resolve) {
    const worker = new Worker(new URL('./dimensionality-reducer', import.meta.url));
    worker.postMessage({
      columnData: dataMetric.data,
      method: method,
      measure: dataMetric.metric,
      options: options,
    });
    worker.onmessage = ({data: {distance, embedding}}) => {
      returnDistanceMatrix ? resolve({distance: distance, embedding: embedding}) : resolve(embedding);
    };
  });
}
