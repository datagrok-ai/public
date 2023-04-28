import {ClusterMatrix} from '@datagrok-libraries/bio/src/trees';

export type ClusteringWorkerInput = {
    distMatArray: Float32Array;
    n: number;
    methodCode: number;
}

export async function getClusterMatrixWorker(
  {distMatArray, n, methodCode}: ClusteringWorkerInput): Promise<ClusterMatrix> {
  return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./clustering-worker.ts', import.meta.url));
    worker.postMessage({distMatArray, n, methodCode});
    worker.onmessage = ({data: {error, clusterMatrix}}): void => {
      worker.terminate();
      error ? reject(error) : resolve(clusterMatrix);
    };
  });
}
