import {ClusterMatrix} from '../types';

/** Calculates Cluster matrix based on distance matrix reduces array. i.e. the top right triangle of the matrix
 * turned into the array.
 * @param distMatArray - distance matrix reduces array.
 * @param n - dimension of the distance matrix, i.e. side length.
 * @param methodCode - method code for the clustering algorithm.
 */
export async function getClusterMatrixWorker(
  distMatArray: Float32Array, n: number, methodCode: number
): Promise<ClusterMatrix> {
  return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./clustering-worker.ts', import.meta.url));
    worker.postMessage({distMatArray, n, methodCode});
    worker.onmessage = ({data: {error, clusterMatrix}}): void => {
      worker.terminate();
      error ? reject(error) : resolve(clusterMatrix);
    };
  });
}
