import { ValidTypes } from "../distance-measures";

/**
 * A worker to perform dimensionality reduction.
 *
 * @param {any[]} columnData The data to process.
 * @param {string} method A method of dimensionality reduction.
 * @param {string} measure A distance metrics.
 * @param {number} cyclesCount Number of iterations to run.
 * @return {Promise<unknown>} Resulting embedding.
 */
export function createDimensinalityReducingWorker(
  dataMetric: ValidTypes,
  method: string,
  cyclesCount?: number,
): Promise<unknown> {
  return new Promise(function(resolve) {
    const worker = new Worker(new URL('./dimensionality-reducer', import.meta.url));
    worker.postMessage({
      columnData: dataMetric.data,
      method: method,
      measure: dataMetric.metric,
      cyclesCount: cyclesCount,
    });
    worker.onmessage = ({data: {embedding}}) => {
      resolve(embedding);
    };
  });
}
