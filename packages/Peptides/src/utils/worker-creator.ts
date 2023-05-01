import {StringMetrics} from '@datagrok-libraries/ml/src/typed-metrics';

export function createDistanceMatrixWorker(peptidesList: string[], metric: StringMetrics): Promise<Float32Array> {
  return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./distance-matrix.worker', import.meta.url));
    worker.postMessage({peptidesList, metric});
    worker.onmessage = ({data: {error, distances}}): void => {
      error ? reject(error) : resolve(distances);
    };
  });
}
