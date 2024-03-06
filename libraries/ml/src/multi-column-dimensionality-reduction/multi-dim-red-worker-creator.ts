import {KnownMetrics} from '../typed-metrics/typed-metrics';
import * as grok from 'datagrok-api/grok';
import {DIMENSIONALITY_REDUCER_TERMINATE_EVENT} from './consts';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';
import {DimReductionMethods} from './types';
import {DistanceAggregationMethod} from '../distance-matrix/types';
import {isNil} from '../distance-matrix/utils';

export async function createMultiDimRedWorker(data: Array<any[]>, metrics: KnownMetrics[], method: DimReductionMethods,
  weights: number[], aggregationMethod: DistanceAggregationMethod,
  options: any, progressFunc?: (epoch: number, epochsLength: number, embedding: number[][]) => void
): Promise<Matrix> {
  if (!options.distanceFnArgs)
    throw new Error('options.distanceFnArgs must be defined');
  if (data.length !== metrics.length || data.length !== options.distanceFnArgs.length || data.length !== weights.length)
    throw new Error('data, metrics and options and weights must have the same length');
  return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./mulit-column-dim-reducer-worker', import.meta.url));
    worker.postMessage({
      columnsData: data,
      distanceMetrics: metrics,
      method: method,
      options: options,
      weights: weights,
      aggregationMethod: aggregationMethod,
    });
    const terminateSub = grok.events.onCustomEvent(DIMENSIONALITY_REDUCER_TERMINATE_EVENT).subscribe(() => {
      try {
        worker?.terminate();
      } finally {
        terminateSub.unsubscribe();
      }
    });
    worker.onmessage = ({data: {error, embedding, epochNum, epochsLength}}) => {
      if (!isNil(epochNum) && !isNil(epochsLength)) {
        progressFunc && progressFunc(epochNum, epochsLength, embedding);
        return;
      }
      terminateSub.unsubscribe();
      if (error)
        reject(error);
      else
        resolve(embedding);
      // terminate the worker after some time. immidiate termination causes crashes.
      setTimeout(() => worker.terminate(), 100);
    };
  });
}
