import { KnownMetrics } from "../typed-metrics";

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
        new Array<Promise<{i: number[], j: number[], similarity: number[], idx: number}>>(this._workerCount);
      for (let idx = 0; idx < this._workerCount; idx++) {
        promises[idx] = new Promise((resolveWorker, rejectWorker) => {
          const startIdx = idx * chunkSize;
          const endIdx = idx === this._workerCount - 1 ? matSize : (idx + 1) * chunkSize;
          this._workers[idx].postMessage({values, startIdx, endIdx, threshold, fnName});
          this._workers[idx].onmessage = ({data: {error, i, j, similarity}}): void => {
            if (error) { rejectWorker(error); } else {
              this._workers[idx].terminate();
              resolveWorker({i, j, similarity, idx});
            }
          };
        });
      }

      const results = await Promise.all(promises);
      const fullSize = results.reduce((acc, val) => acc + val.i.length, 0);
      const i = new Int32Array(fullSize);
      const j = new Int32Array(fullSize);
      const similarity = new Float32Array(fullSize);
      let offset = 0;
      for (const res of results) {
        i.set(res.i, offset);
        j.set(res.j, offset);
        similarity.set(res.similarity, offset);
        offset += res.i.length;
      }
      return {i, j, similarity};
    }
}