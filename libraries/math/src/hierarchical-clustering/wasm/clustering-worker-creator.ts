import {ClusterMatrix} from '../types';
import {LinkageMethodCode} from '../consts';

/** Calculates Cluster matrix based on distance matrix reduces array. i.e. the top right triangle of the matrix
 * turned into the array.
 * @param distMatArray - reduced distance matrix array.
 * @param n - dimension of the distance matrix, i.e. side length.
 * @param methodCode - method code for the clustering algorithm. Method codes can be found in {@link LinkageMethodCode} object.
 * @returns Cluster matrix.
 */
export async function getClusterMatrixWorker(
  distMatArray: Float32Array, n: number, methodCode: number
): Promise<ClusterMatrix> {
  return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./clustering-worker', import.meta.url));
    worker.postMessage({distMatArray, n, methodCode});
    worker.onmessage = ({data: {error, clusterMatrix}}): void => {
      worker.terminate();
      error ? reject(error) : resolve(clusterMatrix);
    };
  });
}
