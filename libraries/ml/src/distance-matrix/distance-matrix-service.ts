import {KnownMetrics} from '../typed-metrics';
import {DistanceAggregationMethod, DistanceAggregationMethods, ClusterRepresentatives} from './types';

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

    /** For each cluster, computes every member's mean distance to the other members and ranks them
     * (rank 1 = medoid, the most representative; higher ranks = progressively less central, enabling
     * Top-N extraction). Each cluster is processed independently by a worker that reduces distances
     * to row sums on the fly, so the full distance matrix is never materialized (memory stays
     * proportional to cluster size, not to the dataset). Clusters run across a dedicated worker pool.
     *
     * @param values      Per-column value arrays (same encoding as `calcMulti`).
     * @param fnNames     Distance metric per column.
     * @param clusters    Each entry lists the row indices (into `values`) of one cluster's members.
     * @param opts        Per-column metric options.
     * @param weights     Per-column weights for multi-column aggregation.
     * @returns           One `ClusterRepresentatives` per cluster, in the same order as `clusters`;
     *                    its arrays are aligned with that cluster's member list. */
    public async calcMedoids(
      values: ArrayLike<any>[], fnNames: KnownMetrics[], clusters: number[][],
      opts: {[_: string]: any}[] = [{}], weights: number[] = [1],
      aggregationMethod: DistanceAggregationMethod = DistanceAggregationMethods.MANHATTAN
    ): Promise<ClusterRepresentatives[]> {
      if (values.length < 1)
        throw new Error('values must contain at least one array');
      if (fnNames.length !== values.length || opts.length !== values.length || weights.length !== values.length)
        throw new Error('values, fnNames, weights and opts must have the same length');

      // ranks members by ascending mean distance; ties broken by ascending original member index
      const rankMembers = (members: number[], meanDistances: number[]): number[] => {
        const order = members.map((_, k) => k).sort((a, b) => {
          const d = meanDistances[a] - meanDistances[b];
          return d !== 0 ? d : members[a] - members[b];
        });
        const ranks = new Array<number>(members.length);
        for (let r = 0; r < order.length; r++)
          ranks[order[r]] = r + 1;
        return ranks;
      };

      const results = new Array<ClusterRepresentatives>(clusters.length);
      const queue: number[] = []; // cluster indices that need a worker (size >= 2)
      for (let c = 0; c < clusters.length; c++) {
        const members = clusters[c];
        if (!members || members.length === 0)
          results[c] = {members: [], meanDistances: [], ranks: []};
        else if (members.length === 1)
          results[c] = {members: [members[0]], meanDistances: [0], ranks: [1]};
        else
          queue.push(c);
      }
      if (queue.length === 0)
        return results;

      const poolSize = Math.min(this._workerCount, queue.length);
      const workers = new Array(poolSize).fill(null)
        .map(() => new Worker(new URL('./medoid-worker', import.meta.url)));
      let nextQueueIdx = 0;

      const subsetFor = (clusterIdx: number): any[][] => {
        const members = clusters[clusterIdx];
        return values.map((col) => members.map((idx) => col[idx]));
      };

      const runWorker = (worker: Worker): Promise<void> => new Promise<void>((resolve, reject) => {
        const dispatch = (): void => {
          if (nextQueueIdx >= queue.length) {
            resolve();
            return;
          }
          const clusterIdx = queue[nextQueueIdx++];
          worker.onmessage = ({data: {error, meanDistances}}): void => {
            if (error) {
              reject(error);
              return;
            }
            const members = clusters[clusterIdx];
            const means = Array.from(meanDistances as Float64Array);
            results[clusterIdx] = {members, meanDistances: means, ranks: rankMembers(members, means)};
            dispatch();
          };
          worker.postMessage({values: subsetFor(clusterIdx), fnNames, opts, weights, aggregationMethod});
        };
        dispatch();
      });

      try {
        await Promise.all(workers.map((w) => runWorker(w)));
      } finally {
        workers.forEach((w) => w.terminate());
      }
      return results;
    }

    public terminate(): void {
      this._workers.forEach((worker) => worker.terminate());
    }
}
