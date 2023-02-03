import { Matrix } from '@datagrok-libraries/utils/src/type-declarations';
import {ValidTypes} from '../typed-metrics';

/**
 * A worker to perform dimensionality reduction.
 *
 * @param {ValidTypes} dataMetric The data to process.
 * @param {string} method A method of dimensionality reduction.
 * @param options - key-value pairs
 * @param returnDistanceMatrix
 * @return {Promise<IReduceDimensionalityResult>} Resulting embedding and distance matrix.
 */
export interface IReduceDimensionalityResult {
  distance: Matrix;
  embedding: Matrix;
}

export function createDimensinalityReducingWorker(dataMetric: ValidTypes, method: string,
      options?: any): Promise<IReduceDimensionalityResult> {

  return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./dimensionality-reducer', import.meta.url));
    worker.postMessage({
      columnData: dataMetric.data,
      method: method,
      measure: dataMetric.metric,
      options: options,
    });
    worker.onmessage = ({data: {error, distance, embedding}}) => {
      if (error)
        reject(error);
      else
        resolve({distance: distance, embedding: embedding});
    };
  });
}
