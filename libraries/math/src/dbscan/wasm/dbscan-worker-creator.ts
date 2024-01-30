
/** Calculates Clusters based on dbscan algorithm in wasm.
 * @param embedX - embeddings x location.
 * @param embedY - embeddings y location.
 * @param minPts - minimum number of points in a cluster.
 * @param epsilon - epsilon.
 * @returns array with cluster indexes.
 */
export async function getDbscanWorker(
  embedX: Float32Array, embedY: Float32Array, epsilon: number, minPts: number
): Promise<Int32Array> {
  return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./clustering-worker', import.meta.url));
    worker.postMessage({embedX, embedY, minPts, epsilon});
    worker.onmessage = ({data: {error, clusters}}): void => {
      worker.terminate();
      error ? reject(error) : resolve(clusters);
    };
  });
}
