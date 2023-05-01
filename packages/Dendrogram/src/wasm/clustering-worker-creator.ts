import {ClusterMatrix} from '@datagrok-libraries/bio/src/trees';

export async function getClusterMatrixWorker(
  distMatArray: Float32Array, n: number, methodCode: number
): Promise<ClusterMatrix> {
  return new Promise(function(resolve, reject) {
    // if (useWorker) {
    const worker = new Worker(new URL('./clustering-worker.ts', import.meta.url));
    worker.postMessage({distMatArray, n, methodCode});
    worker.onmessage = ({data: {error, clusterMatrix}}): void => {
      worker.terminate();
      error ? reject(error) : resolve(clusterMatrix);
    };
    // } else {
    //   try {
    //     getClustersFromDistMatWasm(distMatArray, n, methodCode, false).then((clusterMatrix) => {
    //       resolve(clusterMatrix);
    //     });
    //   } catch (e) {
    //     reject(e);
    //   }
    // }
  });
}
