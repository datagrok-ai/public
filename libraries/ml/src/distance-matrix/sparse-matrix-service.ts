import { getSimilarityFromDistance } from "../distance-metrics-methods";
import { BitArrayMetricsNames, KnownMetrics, Measure } from "../typed-metrics";
import { isNil } from "./utils";

export class SparseMatrixService {
    private _workerCount: number;
    private _workers: Worker[];
    constructor() {
      this._workerCount = Math.max(navigator.hardwareConcurrency - 2, 1);
      this._workers = new Array(this._workerCount).fill(null)
        .map(() => new Worker(new URL('./sparse-matrix-worker', import.meta.url)));
    }

    public async calc<T>(values: Array<T> | ArrayLike<T>, fnName: KnownMetrics, threshold: number) {
      const matSize = values.length * (values.length - 1) / 2;
      const chunkSize = Math.floor(matSize / this._workerCount);
      const promises =
        new Array<Promise<{i: Int32Array, j: Int32Array, distance: Float32Array, idx: number}>>(this._workerCount);
      for (let idx = 0; idx < this._workerCount; idx++) {
        promises[idx] = new Promise((resolveWorker, rejectWorker) => {
          const startIdx = idx * chunkSize;
          const endIdx = idx === this._workerCount - 1 ? matSize : (idx + 1) * chunkSize;
          this._workers[idx].postMessage({values, startIdx, endIdx, threshold, fnName});
          this._workers[idx].onmessage = ({data: {error, i, j, distance}}): void => {
            if (error) { rejectWorker(error); } else {
              this._workers[idx].terminate();
              resolveWorker({i, j, distance, idx});
            }
          };
        });
      }

      const results = await Promise.all(promises);
      const fullSize = results.reduce((acc, val) => acc + val.i.length, 0);
      const i = new Int32Array(fullSize);
      const j = new Int32Array(fullSize);
      const distance = new Float32Array(fullSize);
      let offset = 0;
      for (const res of results) {
        i.set(res.i, offset);
        j.set(res.j, offset);
        distance.set(res.distance, offset);
        offset += res.i.length;
      }
      return {i, j, distance};
    }

    public static calcSync<T> (values: Array<T> | ArrayLike<T>, fnName: KnownMetrics, distanceFn: Function, threshold: number) {
      const i: number[] = [];
      const j: number[] = [];
      const distances: number[] = [];
      let cnt = 0;
      let mi = 0;
      let mj = 0;
      const fullSize = values.length * (values.length - 1) / 2;
      while (cnt < fullSize) {
        //const value = seq1List[mi] && seq1List[mj] ? hamming(seq1List[mi], seq1List[mj]) : 0;
        const value = !isNil(values[mi]) && !isNil(values[mj]) ?
          distanceFn(values[mi], values[mj]) : 1;
        const similarity = Object.values(BitArrayMetricsNames).some((a) => a === fnName) ? getSimilarityFromDistance(value) : 1 - value;
        if (similarity >= threshold) {
          i.push(mi);
          j.push(mj);
          distances.push(value);
        }
        cnt++;
        mj++;
        if (mj === values.length) {
          mi++;
          mj = mi + 1;
        }
      }
    
      const iArray = new Int32Array(i);
      const jArray = new Int32Array(j);
      const distanceArray = new Float32Array(distances);

      return {i: iArray, j: jArray, distance: distanceArray};
    }
}