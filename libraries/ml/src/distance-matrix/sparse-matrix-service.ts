import {KnownMetrics} from '../typed-metrics';
import {DistanceAggregationMethod, DistanceAggregationMethods} from './types';
import {insertSmaller, isNil} from './utils';

export type SparseMatrixResult = {
  i: Int32Array | Uint32Array,
  j: Int32Array | Uint32Array,
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

    public async calcMultiColumn(values: Array<any[]>, fnNames: KnownMetrics[],
      threshold: number, opts: {[_: string]: any}[] = [{}], weights: number[] = [1],
      aggregationMethod: DistanceAggregationMethod = DistanceAggregationMethods.EUCLIDEAN
    ) {
      const matSize = values[0].length * (values[0].length - 1) / 2;
      const chunkSize = Math.floor(matSize / this._workerCount);

      const minThreshold = values[0].length > 20_000 ?
        await this.getMinimalThreshold(values, fnNames, opts, weights, aggregationMethod) : 0;
      if (threshold < minThreshold) {
        console.log(`using threshold ${minThreshold}`);
        threshold = minThreshold;
      }
      opts.forEach((_, i) => opts[i]['threshold'] = threshold);
      const promises =
        new Array<Promise<SparseMatrixResult>>(this._workerCount);

      const workers = new Array(this._workerCount)
        .fill(null).map(() => new Worker(new URL('./sparse-matrix-worker', import.meta.url)));
      for (let idx = 0; idx < this._workerCount; idx++) {
        promises[idx] = new Promise((resolveWorker, rejectWorker) => {
          const startIdx = idx * chunkSize;
          const endIdx = idx === this._workerCount - 1 ? matSize : (idx + 1) * chunkSize;
          if (endIdx <= startIdx)
            resolveWorker({i: new Int32Array(0), j: new Int32Array(0), distance: new Float32Array(0), idx});
          workers[idx].postMessage({values, startIdx, endIdx, threshold, fnNames, opts, weights, aggregationMethod});
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

    public async calc<T>(values: Array<T>, fnName: KnownMetrics, threshold: number, opts: {[_: string]: any} = {}) {
      //size of full matrix
      return await this.calcMultiColumn([values], [fnName], threshold, [opts], [1]);
    }

    public async getKNN(
      values: Array<any>, fnName: KnownMetrics, nNeighbours: number = 15, opts: {[_: string]: any} = {}
    ) {
      return await this.multiColumnKNN([values], [fnName], nNeighbours, [opts], [1]);
    }

    public async getThresholdKNN(
      values: Array<any>, fnName: KnownMetrics, threshold: number = 0.8, opts: {[_: string]: any} = {}
    ) {
      return await this.multiColumnThresholdKnn([values], [fnName], threshold, [opts], [1]);
    }

    public async multiColumnThresholdKnn(values: Array<Array<any>>, fnNames: KnownMetrics[], threshold: number = 0.8,
      opts: {[_: string]: any}[], weights: number[],
      aggregationMethod: DistanceAggregationMethod = DistanceAggregationMethods.EUCLIDEAN
    ) {
      if (values.length !== fnNames.length || values.length !== opts.length || values.length !== weights.length)
        throw new Error('values, distance functions, options and weights arrays should have the same length');

      if (values.some((v) => v.length !== values[0].length))
        throw new Error('all values arrays should have the same length');

      const matSize = values[0].length * (values[0].length - 1) / 2;
      const chunkSize = Math.floor(matSize / this._workerCount);
      const promises =
        new Array<Promise<KnnResult>>(this._workerCount);
      const workers = new Array(this._workerCount)
        .fill(null).map(() => new Worker(new URL('./knn-threshold-worker', import.meta.url)));
      for (let idx = 0; idx < this._workerCount; idx++) {
        promises[idx] = new Promise((resolveWorker, rejectWorker) => {
          const startIdx = idx * chunkSize;
          const endIdx = idx === this._workerCount - 1 ? matSize : (idx + 1) * chunkSize;
          if (endIdx <= startIdx)
            resolveWorker({knnDistances: new Array(0), knnIndexes: new Array(0)});
          workers[idx].postMessage({values, startIdx, endIdx, fnNames, opts, threshold, weights, aggregationMethod});
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
      const knnSizes = new Int32Array(values[0].length);
      for (const res of results) {
        for (let i = 0; i < values[0].length; ++i)
          knnSizes[i] += res.knnIndexes[i]?.length ?? 0;
      }
      const knnRes: KnnResult = {
        knnDistances: new Array(values[0].length).fill(null).map((_, i) => new Array<number>(knnSizes[i])),
        knnIndexes: new Array(values[0].length).fill(null).map((_, i) => new Array<number>(knnSizes[i]))};
      for (const res of results) {
        for (let i = 0; i < values[0].length; ++i) {
          for (let j = 0; j < (res.knnDistances[i]?.length ?? 0); ++j) {
            knnRes.knnDistances[i][knnSizes[i] - 1] = res.knnDistances[i][j];
            knnRes.knnIndexes[i][knnSizes[i] - 1] = res.knnIndexes[i][j];
            knnSizes[i] -= 1;
          }
        }
      }
      return knnRes;
    }

    public async multiColumnKNN(values: Array<Array<any>>, fnNames: KnownMetrics[], nNeighbours: number = 15,
      opts: {[_: string]: any}[], weights: number[],
      aggregationMethod: DistanceAggregationMethod = DistanceAggregationMethods.EUCLIDEAN
    ) {
      if (values.length !== fnNames.length || values.length !== opts.length || values.length !== weights.length)
        throw new Error('values, distance functions, options and weights arrays should have the same length');

      if (values.some((v) => v.length !== values[0].length))
        throw new Error('all values arrays should have the same length');

      const matSize = values[0].length * (values[0].length - 1) / 2;
      const chunkSize = Math.floor(matSize / this._workerCount);
      const promises =
        new Array<Promise<KnnResult>>(this._workerCount);
      const workers = new Array(this._workerCount)
        .fill(null).map(() => new Worker(new URL('./knn-worker', import.meta.url)));
      for (let idx = 0; idx < this._workerCount; idx++) {
        promises[idx] = new Promise((resolveWorker, rejectWorker) => {
          const startIdx = idx * chunkSize;
          const endIdx = idx === this._workerCount - 1 ? matSize : (idx + 1) * chunkSize;
          if (endIdx <= startIdx)
            resolveWorker({knnDistances: new Array(0), knnIndexes: new Array(0)});
          workers[idx].postMessage({values, startIdx, endIdx, fnNames, opts, nNeighbours, weights, aggregationMethod});
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
      const knnRes: KnnResult = {
        knnDistances: new Array(values[0].length).fill(null).map(() => new Array<number>(nNeighbours).fill(99999)),
        knnIndexes: new Array(values[0].length).fill(null).map(() => new Array<number>(nNeighbours).fill(-1))};
      for (const res of results) {
        for (let i = 0; i < values[0].length; ++i) {
          for (let j = 0; j < (res.knnDistances[i]?.length ?? 0); ++j)
            insertSmaller(knnRes.knnDistances[i], knnRes.knnIndexes[i], res.knnDistances[i][j], res.knnIndexes[i][j]);
        }
      }
      return knnRes;
    }

    public async getSampleDistances(values: Array<any[]>,
      fnNames: KnownMetrics[], opts: {[_: string]: any}[] = [], weights: number[],
      aggregationMethod: DistanceAggregationMethod = DistanceAggregationMethods.EUCLIDEAN): Promise<Float32Array> {
      const thresholdWorkers = new Array(this._workerCount).fill(null)
        .map(() => new Worker(new URL('./sparse-matrix-threshold-worker', import.meta.url)));

      try {
        const matSize = values[0].length * (values[0].length - 1) / 2;
        const chunkSize = Math.floor(matSize / this._workerCount);
        const maxSampleSize = 1_000_000;
        const sampleSise = Math.max(Math.min(matSize / 1000, maxSampleSize), Math.min(matSize, maxSampleSize));
        const testSetSizePerWorker = Math.floor(sampleSise / this._workerCount);
        const tPromises = new Array<Promise<{distance: Float32Array}>>(this._workerCount);

        for (let idx = 0; idx < this._workerCount; idx++) {
          tPromises[idx] = new Promise((resolveWorker, rejectWorker) => {
            const startIdx = idx * chunkSize;
            const endIdx = idx === this._workerCount - 1 ? matSize : (idx + 1) * chunkSize;
            thresholdWorkers[idx].postMessage({
              values: values, startIdx, endIdx, sampleLength: testSetSizePerWorker,
              fnNames, opts, weights, aggregationMethod
            });
            thresholdWorkers[idx].onmessage = ({data: {error, distance}}): void => {
              thresholdWorkers[idx].terminate();
              if (error) rejectWorker(error); else
                resolveWorker({distance});
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

    private async getMinimalThreshold(values: Array<any[]>,
      fnNames: KnownMetrics[], opts: {[_: string]: any}[] = [], weights: number[],
      aggregationMethod: DistanceAggregationMethod = DistanceAggregationMethods.EUCLIDEAN) {
      //We need to calculate the minimal threshold first,
      //in order to get matrix such that it does not exceed the maximum size of 1GB
      //we have 3 return arrays, each 4 bites per element, so if the maximum size of the matrix is 1GB,
      const maxSparseMatrixSize = 70_000_000;
      try {
        const matSize = values.length * (values.length - 1) / 2;
        const distance = await this.getSampleDistances(values, fnNames, opts, weights, aggregationMethod);
        const fractionIndex = Math.floor(maxSparseMatrixSize / matSize * distance.length);
        const threshold = 1 - distance[fractionIndex];
        // threshold = Math.max(threshold, 0.3);
        return threshold;
      } catch (e) {
        console.error(e);
        return 0.5;
      }
    }

    public static calcSync<T>(
      values: Array<T> | ArrayLike<T>, fnName: KnownMetrics, distanceFn: Function, threshold: number
    ) {
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
        const similarity = 1 - value;
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
