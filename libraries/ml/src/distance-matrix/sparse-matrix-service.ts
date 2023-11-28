import { getSimilarityFromDistance } from "../distance-metrics-methods";
import { BitArrayMetricsNames, KnownMetrics, Measure } from "../typed-metrics";
import { insertSmaller, isNil } from "./utils";

export type SparseMatrixResult = {
  i: Int32Array,
  j: Int32Array,
  distance: Float32Array,
  idx?: number
};

export type KnnResult = {
  knnDistances: number[][],
  knnIndexes: number[][]
}
export class SparseMatrixService {
    private _workerCount: number;
    constructor() {
      this._workerCount = Math.max(navigator.hardwareConcurrency - 2, 1);
    }

    public async calc<T>(values: Array<T>, fnName: KnownMetrics, threshold: number, opts: {[_: string]: any} = {}) {


      //size of full matrix
      const matSize = values.length * (values.length - 1) / 2;
      const chunkSize = Math.floor(matSize / this._workerCount);

      const minThreshold = values.length > 20_000 ? await this.getMinimalThreshold(values, fnName, opts) : 0;
      if (threshold < minThreshold) {
        console.log(`using threshold ${minThreshold}`);
        threshold = minThreshold;
      }
      opts['threshold'] = threshold;
      const promises =
        new Array<Promise<SparseMatrixResult>>(this._workerCount);

      const workers = new Array(this._workerCount).fill(null).map(() => new Worker(new URL('./sparse-matrix-worker', import.meta.url)));
      for (let idx = 0; idx < this._workerCount; idx++) {
        promises[idx] = new Promise((resolveWorker, rejectWorker) => {
          const startIdx = idx * chunkSize;
          const endIdx = idx === this._workerCount - 1 ? matSize : (idx + 1) * chunkSize;
          if (endIdx <= startIdx) 
            resolveWorker({i: new Int32Array(0), j: new Int32Array(0), distance: new Float32Array(0), idx});
          workers[idx].postMessage({values, startIdx, endIdx, threshold, fnName, opts});
          workers[idx].onmessage = ({data: {error, i, j, distance}}): void => {
            if (error) { 
              workers[idx].terminate();
              rejectWorker(error); 
            } else {
              workers[idx].terminate();
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
      // setting the results
      for (const res of results) {
        i.set(res.i, offset);
        j.set(res.j, offset);
        distance.set(res.distance, offset);
        offset += res.i.length;
      }
      return {i, j, distance};
    }

    public async getKNN<T>(values: Array<T>, fnName: KnownMetrics, nNeighbours: number = 15, opts: {[_: string]: any} = {}) {
      const matSize = values.length * (values.length - 1) / 2;
      const chunkSize = Math.floor(matSize / this._workerCount);
      const promises =
        new Array<Promise<KnnResult>>(this._workerCount);
      const workers = new Array(this._workerCount).fill(null).map(() => new Worker(new URL('./knn-worker', import.meta.url)));
      for (let idx = 0; idx < this._workerCount; idx++) {
        promises[idx] = new Promise((resolveWorker, rejectWorker) => {
          const startIdx = idx * chunkSize;
          const endIdx = idx === this._workerCount - 1 ? matSize : (idx + 1) * chunkSize;
          if (endIdx <= startIdx) 
            resolveWorker({knnDistances: new Array(0), knnIndexes: new Array(0)});
          workers[idx].postMessage({values, startIdx, endIdx, fnName, opts, nNeighbours});
          workers[idx].onmessage = ({data: {error, knnDistances, knnIndexes}}): void => {
            if (error) { 
              workers[idx].terminate();
              rejectWorker(error); 
            } else {
              workers[idx].terminate();
              resolveWorker({knnDistances, knnIndexes});
            }
          };
        });
      }
  
      const results = await Promise.all(promises);
      const knnRes: KnnResult = {knnDistances: new Array(values.length).fill(null).map(() => new Array<number>(nNeighbours).fill(99999)),
        knnIndexes: new Array(values.length).fill(null).map(() => new Array<number>(nNeighbours).fill(-1))};
      for (const res of results) {
        for (let i = 0; i < values.length; ++i) {
          for (let j = 0; j < res.knnDistances[i]?.length ?? 0; ++j) {
            insertSmaller(knnRes.knnDistances[i], knnRes.knnIndexes[i], res.knnDistances[i][j], res.knnIndexes[i][j]);
          }
        }
      }
      return knnRes;
    }

    public async getSampleDistances<T>(values: Array<T>, fnName: KnownMetrics, opts: {[_: string]: any} = {}): Promise<Float32Array> {
      const thresholdWorkers = new Array(this._workerCount).fill(null)
        .map(() => new Worker(new URL('./sparse-matrix-threshold-worker', import.meta.url)));
      // data may be sorted by clusters, which will hinder the random sampling
      // so we shuffle it first
      const shuffledValues = values.slice();
      for (let i = shuffledValues.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [shuffledValues[i], shuffledValues[j]] = [shuffledValues[j], shuffledValues[i]];
      }
  
      try {
        const matSize = values.length * (values.length - 1) / 2;
        const chunkSize = Math.floor(matSize / this._workerCount);
        const maxSampleSize = 1_000_000;
        const sampleSise = Math.max(Math.min(matSize / 1000, maxSampleSize), Math.min(matSize, maxSampleSize));
        const testSetSizePerWorker = Math.floor(sampleSise / this._workerCount);
        const tPromises = new Array<Promise<{distance: Float32Array}>>(this._workerCount);

        for (let idx = 0; idx < this._workerCount; idx++) {
          tPromises[idx] = new Promise((resolveWorker, rejectWorker) => {
            const startIdx = idx * chunkSize;
            const endIdx = idx === this._workerCount - 1 ? matSize : (idx + 1) * chunkSize;
            thresholdWorkers[idx].postMessage({values: shuffledValues, startIdx, endIdx, sampleLength: testSetSizePerWorker, fnName, opts});
            thresholdWorkers[idx].onmessage = ({data: {error, distance}}): void => {
              thresholdWorkers[idx].terminate();
              if (error) { rejectWorker(error); } else {
                resolveWorker({distance});
              }
            };
          });
        }

        const results = await Promise.all(tPromises);
        const fullSize = results.reduce((acc, val) => acc + val.distance.length, 0);
        const distance = new Float32Array(fullSize);
        let offset = 0;
        for (const res of results) {
          distance.set(res.distance, offset);
          offset += res.distance.length;
        }
        distance.sort();

        return distance;
      } catch (e) {
        thresholdWorkers?.forEach((w) => w?.terminate());
        console.error(e);
        return new Float32Array(1).fill(0.5);
      }
    }

    private async getMinimalThreshold<T>(values: Array<T>, fnName: KnownMetrics, opts: {[_: string]: any} = {}) {

      //We need to calculate the minimal threshold first,
      //in order to get matrix such that it does not exceed the maximum size of 1GB
      //we have 3 return arrays, each 4 bites per element, so if the maximum size of the matrix is 1GB,
      const max_Sparse_matrix_size = 70_000_000;
      try {
        const matSize = values.length * (values.length - 1) / 2;
        const distance = await this.getSampleDistances(values, fnName, opts);
        const fractionIndex = Math.floor(max_Sparse_matrix_size / matSize * distance.length);
        let threshold = 1 - distance[fractionIndex];
        // threshold = Math.max(threshold, 0.3);
        return threshold;
      } catch (e) {
        console.error(e);
        return 0.5;
      }
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