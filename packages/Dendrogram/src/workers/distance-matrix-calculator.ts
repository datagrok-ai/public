/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {DistanceFunctionNames} from './distance-functions';


export class DistanceMatrixWorker {
  private _worker: Worker;

  public constructor() {
    this._worker = new Worker(new URL('./distance-matrix-worker.ts', import.meta.url));
  };

  public async calc<T>(values: Array<T> | ArrayLike<T>, fnName: DistanceFunctionNames<T>): Promise<Float32Array> {
    this._worker.postMessage({values, fnName});
    return new Promise((resolve, reject) => {
      this._worker.onmessage = ({data: {error, distanceMatrixData}}): void => {
        error ? reject(error) : resolve(distanceMatrixData);
      };
    });
  }

  public terminate(): void {
    this._worker.terminate();
  }
}
