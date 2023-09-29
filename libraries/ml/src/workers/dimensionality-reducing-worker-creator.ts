import {ValidTypes} from '../typed-metrics/typed-metrics';
import {IReduceDimensionalityResult} from '../reduce-dimensionality';

/**
 * A worker to perform dimensionality reduction.
 *
 * @param {ValidTypes} dataMetric The data to process.
 * @param {string} method A method of dimensionality reduction.
 * @param {any}options - key-value pairs
 * @param {boolean}parallelDistanceWorkers - whether to use parallel distance matrix workers
 * @return {Promise<IReduceDimensionalityResult>} Resulting embedding and distance matrix.
 */
export async function createDimensinalityReducingWorker(dataMetric: ValidTypes, method: string,
  options?: any, parallelDistanceWorkers?: boolean,
  progressFunc?: (epoch: number, epochsLength: number, embedding: number[][]) => void
): Promise<IReduceDimensionalityResult> {
  return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./dimensionality-reducer', import.meta.url));
    worker.postMessage({
      columnData: dataMetric.data,
      method: method,
      measure: dataMetric.metric,
      options: options,
      parallelDistanceWorkers: parallelDistanceWorkers,
    });
    worker.onmessage = ({data: {error, distance, embedding, epochNum, epochsLength}}) => {
      if (epochNum && epochsLength) {
        progressFunc && progressFunc(epochNum, epochsLength, embedding);
        return;
      }
      if (error)
        reject(error);
      else
        resolve({distance: distance, embedding: embedding});
      // terminate the worker after some time. immidiate termination causes crashes.
      setTimeout(() => worker.terminate(), 0);
    };
  });
}
