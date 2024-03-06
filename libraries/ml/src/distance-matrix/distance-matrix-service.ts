import {KnownMetrics} from '../typed-metrics';
import {DistanceAggregationMethod, DistanceAggregationMethods} from './types';

export class DistanceMatrixService {
    private _workers: Worker[];
    private _workerCount: number;
    private _terminateOnComplete: boolean;
    public constructor(useConcurrentWorkers = true, terminateOnComplete = true) {
      const threadCount = navigator.hardwareConcurrency;
      this._workerCount = useConcurrentWorkers ? Math.max(threadCount - 2, 1) : 1;
      this._workers = new Array(this._workerCount).fill(null)
        .map(() => new Worker(new URL('./distance-matrix-worker', import.meta.url)));
      this._terminateOnComplete = terminateOnComplete;
    };

    public async calc(values: ArrayLike<any>, fnName: KnownMetrics,
      normalize = true, opts?: {[_: string]: any}): Promise<Float32Array> {
      return await this.calcMulti([values], [fnName], normalize,
        [opts ?? {}], [1], DistanceAggregationMethods.MANHATTAN);
    }

    public async calcMulti(
      values: ArrayLike<any>[], fnNames: KnownMetrics[],
      normalize = true, opts: {[_: string]: any}[] = [{}],
      weights: number[] = [1], aggregationMethod: DistanceAggregationMethod = DistanceAggregationMethods.MANHATTAN
    ): Promise<Float32Array> {
      if (values.length < 1)
        throw new Error('values must contain at least one array');
      if (fnNames.length !== values.length || opts.length !== values.length || weights.length !== values.length)
        throw new Error('values, fnNames, weights and opts must have the same length');
      return new Promise(async (resolve, reject) => {
        try {
          const len = values[0].length;
          const promises = new Array<Promise<void>>(this._workerCount);
          const totalLength = len * (len - 1) / 2; // size of reduced distance matrix
          this._workerCount = Math.min(this._workerCount, totalLength);
          const chunkSize = totalLength / this._workerCount;
          const distanceMatrix = new Float32Array(totalLength);
          let endRow = 0;
          let endCol = 1;
          // minmax for normalization
          let lmin = 0;
          let lmax = Number.MIN_VALUE;
          for (let i = 0; i < this._workerCount; i++) {
            const start = Math.floor(i * chunkSize);
            const end = (i === this._workerCount - 1) ? totalLength : Math.floor((i + 1) * chunkSize);
            const startRow = endRow;
            const startCol = endCol;
            if (i !== this._workerCount - 1) {
              // These formulas map the linear index to the upper triangular matrix indices
              endRow = len - 2 - Math.floor(Math.sqrt(-8 * end + 4 * len * (len - 1) - 7) / 2 - 0.5);
              endCol = end - len * endRow + Math.floor((endRow + 1) * (endRow + 2) / 2);
            }
            this._workers[i].postMessage(
              {values, fnNames, startRow, startCol, chunckSize: end - start, opts, weights, aggregationMethod}
            );
            promises[i] = new Promise((resolveWorker, rejectWorker) => {
              this._workers[i].onmessage = ({data: {error, distanceMatrixData, min, max}}): void => {
                this._terminateOnComplete && setTimeout(() => this._workers[i].terminate());
                if (error) {
                  rejectWorker(error);
                } else {
                  distanceMatrix.set(distanceMatrixData, start);
                  if (min < lmin)
                    lmin = min;
                  if (max > lmax)
                    lmax = max;
                  resolveWorker();
                }
              };
            });
          }
          await Promise.all(promises);
          if (normalize)
            distanceMatrix.forEach((value, index) => { distanceMatrix[index] = (value - lmin) / (lmax - lmin); });
          resolve(distanceMatrix);
        } catch (e) {
          reject(e);
        }
      });
    }

    public terminate(): void {
      this._workers.forEach((worker) => worker.terminate());
    }
}
