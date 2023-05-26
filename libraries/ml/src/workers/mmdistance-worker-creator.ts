import * as DG from 'datagrok-api/dg';
import {MmDistanceFunctionsNames} from '../macromolecule-distance-functions';
import {Matrix} from '@datagrok-libraries/utils/src/type-declarations';

export function createMMDistanceWorker(valCol: DG.Column<any>,
  dataMetric: MmDistanceFunctionsNames, options?: any,
): Promise<Matrix> {
  return new Promise(function(resolve, reject) {
    const worker = new Worker(new URL('./mmdistance-worker', import.meta.url));
    worker.postMessage({
      values: valCol.toList(),
      fnName: dataMetric,
      options,
    });
    worker.onmessage = ({data: {distanceMatrixData, error}}) => {
      if (error)
        reject(error);
      else
        resolve(distanceMatrixData);
      worker.terminate();
    };
  });
}
