import {DistanceAggregationMethod} from '../distance-matrix/types';
import {KnownMetrics} from '../typed-metrics';

export * from './marcov-cluster';
export * from './types';

export async function createMCLWorker(data: any[][], threshold: number,
  weights: number[], aggregationMethod: DistanceAggregationMethod,
  distanceFns: KnownMetrics[], distanceFnArgs: any[], maxIterations: number = 10) {
  const worker = new Worker(new URL('mcl-worker', import.meta.url));
  worker.postMessage({data, threshold, weights, aggregationMethod, distanceFns, distanceFnArgs, maxIterations});
  return new Promise<{
    clusters: number[], embedX: Float32Array, embedY: Float32Array, is: number[], js: number[]}>((resolve, reject) => {
      worker.onmessage = (event) => {
        setTimeout(() => worker.terminate(), 100);
        resolve(event.data.res);
      };
      worker.onerror = (event) => {
        setTimeout(() => worker.terminate(), 100);
        reject(event);
      };
    });
}
