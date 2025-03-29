import {Coordinates, DistanceMetric, Matrix, Options, Vector, Vectors}
  from '@datagrok-libraries/utils/src/type-declarations';
import {AvailableDataTypes, AvailableMetrics, KnownMetrics, Measure, isBitArrayMetric} from '../typed-metrics';
import {DistanceMatrixService} from '../distance-matrix';
import {DistanceAggregationMethod} from '../distance-matrix/types';
import {UMAP} from '../umap';
import {assert, transposeMatrix} from '@datagrok-libraries/utils/src/vector-operations';
import {KnnResult, SparseMatrixService} from '../distance-matrix/sparse-matrix-service';
import {DimReductionMethods} from './types';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import seedRandom from 'seedrandom';
//import {TSNE} from '@keckelt/tsne';
import {TSNE} from '../t-sne/t-sne';
import {multiColWebGPUKNN} from '@datagrok-libraries/math';
import {WebGPUUMAP} from '@datagrok-libraries/math/src/webGPU/umap';

export interface IUMAPOptions {
    learningRate?: number;
    nComponents?: number;
    nEpochs?: number;
    nNeighbors?: number;
    spread?: number;
    minDist?: number;
    progressFunc?: (epoc: number, epochsLength: number, embeddings: number[][]) => void;
    randomSeed?: string;
    useWebGPU?: boolean;
}

export interface ITSNEOptions {
    epsilon?: number;
    perplexity?: number;
    dim?: number;
}

export interface IDimReductionParam<T extends number | string | boolean = number> {
    uiName: string;
    value: (T extends number ? number : T extends string ? string : boolean) | null;
    tooltip: string;
    type?: T extends number ? 'number' : T extends string ? 'string' : 'boolean';
    placeholder?: string;
    min?: number;
    max?: number;
    step?: number;
    disable?: boolean;
    disableTooltip?: string;
}

abstract class MultiColumnReducer {
    protected data: Vectors[];
    protected weights: number[];
    protected aggregationMethod: DistanceAggregationMethod;
    constructor(options: Options) {
      this.data = options.data;
      this.weights = options.weights;
      this.aggregationMethod = options.aggregationMethod;
    }

    /** Embeds the data given into the two-dimensional space.
     * @return {any} Cartesian coordinate of this embedding and distance matrix where applicable. */
    abstract transform(): Promise<Matrix>;
}

class TSNEReducer extends MultiColumnReducer {
    protected reducer: TSNE;
    protected iterations: number;
    protected distanceFnames: KnownMetrics[];
    protected distanceFns: ((a: any, b: any) => number)[];
    protected distanceFnArgs: {[_: string]: any}[];
    protected progressFunc?: (epoc: number, epochsLength: number, embeddings?: number[][]) => void;
    /**
     * Creates an instance of TSNEReducer.
     * @param {Options} options Options to pass to the constructor.
     * @memberof TSNEReducer
     */
    constructor(options: Options) {
      super(options);
      const randomSeed: string = options.randomSeed ?? Date();
      const randomFn = seedRandom(randomSeed);
      options.dim = 2; // TODO: make it configurable
      options.random = randomFn;
      this.reducer = new TSNE(options);
      this.iterations = options?.iterations ?? this.reducer.getIterSize(this.data[0].length);
      this.distanceFnames = options.distanceFnames;
      this.distanceFns = options.distanceFns;
      this.distanceFnArgs = options.distanceFnArgs;
      this.progressFunc = options.progressFunc;
    }

    /**
     * Embeds the data given into the two-dimensional space using t-SNE method.\
     * @param {boolean} [parallelDistanceWorkers] Whether to use parallel distance workers.
     * @return {any} Cartesian coordinate of this embedding and distance matrix where applicable.
     */
    public async transform(): Promise<Matrix> {
      if (this.data[0].length > 10000)
        throw new Error('Maximum number of samples for T-SNE is 10000');
      const matrixService = new DistanceMatrixService(true, false);
      try {
        // const aggregate = getAggregationFunction(this.aggregationMethod, this.weights);
        // const distances: Array<Float32Array> = [];
        // for (let i = 0; i < this.data.length; ++i) {
        //   const dist = await matrixService.calc(this.data[i], this.distanceFnames[i], false, this.distanceFnArgs[i]);
        //   distances.push(dist);
        // }
        // const distance = new Float32Array(distances[0].length).fill(0);
        // for (let i = 0; i < distances[0].length; ++i)
        //   distance[i] = aggregate(distances.map((d) => d[i]));
        const distance = await matrixService.calcMulti(this.data, this.distanceFnames, false, this.distanceFnArgs,
          this.weights, this.aggregationMethod);
        matrixService.terminate();
        // const matrixProxy = distanceMatrixProxy(distance, this.data[0].length);
        this.reducer.initDataDist(distance, this.data[0].length);

        for (let i = 0; i < this.iterations; ++i) {
          this.reducer.step(); // every time you call this, solution gets better
          if (this.progressFunc)
            this.progressFunc(i, this.iterations, []);
        }
        return this.reducer.getSolution();
      } catch (e) {
        matrixService.terminate();
        throw e;
      }
    }
}

class UMAPReducer extends MultiColumnReducer {
    protected reducer: UMAP;
    protected distanceFnames: KnownMetrics[];
    protected distanceFns: Function[];
    protected vectors: number[];
    protected progressFunc?: (epoc: number, epochsLength: number, embedding: number[][]) => void;
    protected distanceFnArgs: {[_: string]: any}[];
    protected useWebGPU: boolean = false;
    protected webGPUReducer: WebGPUUMAP | undefined;
    /**
     * Creates an instance of UMAPReducer.
     * @param {Options} options Options to pass to the constructor.
     * @memberof UMAPReducer
     */
    constructor(options: Options) {
      const randomSeed: string = options.randomSeed ?? Date();
      const randomFn = seedRandom(randomSeed);
      super(options);
      assert('distanceFnames' in options);
      assert('distanceFns' in options);
      this.distanceFnArgs = options.distanceFnArgs;
      this.distanceFns = options.distanceFns!;
      this.progressFunc = options.progressFunc;
      this.useWebGPU = options.useWebGPU ?? false;
      this.distanceFnames = options.distanceFnames!;
      options.nComponents = 2; // TODO: make it configurable
      //Umap uses vector indexing, so we need to create an array of vectors as indeces.
      this.vectors = new Array(this.data[0].length).fill(0).map((_, i) => i);

      if (this.data[0].length <= (options.nNeighbors ?? 15))
        options.nNeighbors = this.data[0].length - 1;
      options.random = randomFn;
      this.reducer = new UMAP(options);
      // this.reducer.distanceFn = this._encodedDistance.bind(this);
    }

    /**
     * Embeds the data given into the two-dimensional space using UMAP method.
     * @param {boolean} [_parallelDistanceWorkers] Whether to use parallel distance matrix workers.
     * @return {any} Cartesian coordinate of this embedding.
     */
    public async transform(_parallelDistanceWorkers?: boolean): Promise<Matrix> {
      console.time('knn graph');
      let knnGraph: KnnResult | undefined = undefined;
      if (this.useWebGPU) {
        try {
          knnGraph = await multiColWebGPUKNN(
            this.data, this.reducer.neighbors, this.distanceFnames as any, this.aggregationMethod as any,
            this.weights, this.distanceFnArgs
          ) as unknown as KnnResult;
        } catch (e) {
          console.error(e);
        }
      }
      if (!knnGraph) {
        if (this.useWebGPU)
          console.error('WEBGPU KNN failed, falling back to multithreaded CPU implementation');
        knnGraph = await new SparseMatrixService()
          .multiColumnKNN(this.data, this.distanceFnames, this.reducer.neighbors,
            this.distanceFnArgs, this.weights, this.aggregationMethod);
      }
      console.timeEnd('knn graph');
      if (!knnGraph)
        throw new Error('Failed to compute KNN graph');

      // go through KNN graph and remove -99999 or -1 and such values
      const maxIndex = this.data[0].length;
      for (let i = 0; i < knnGraph.knnIndexes.length; ++i) {
        for (let j = 0; j < knnGraph.knnIndexes[i].length; ++j) {
          if (knnGraph.knnDistances[i][j] < 0)
            knnGraph.knnDistances[i][j] = 0;
          if (knnGraph.knnIndexes[i][j] < 0 || knnGraph.knnIndexes[i][j] >= maxIndex || knnGraph.knnDistances[i][j] > 10) {
            knnGraph.knnIndexes[i][j] = i == j ? j + 1 : j;
            knnGraph.knnDistances[i][j] = 10;
          }
        }
      }

      if (this.useWebGPU) {
        this.webGPUReducer = new WebGPUUMAP(this.vectors.length, {
          nComponents: this.reducer.nComponents,
          gamma: this.reducer.repulsionStrength,
          alpha: this.reducer.learningRate,
          nNeighbours: this.reducer.neighbors,
          spread: this.reducer.spread,
          minDist: this.reducer.minDist,
          negativeSampleRate: this.reducer.negativeSampleRate,
          localConnectivity: this.reducer.localConnectivity,
          setOpMixRatio: this.reducer.setOpMixRatio,
        });
        this.webGPUReducer.setPrecomputedKNN(knnGraph.knnIndexes, knnGraph.knnDistances);
      } else {
        this.reducer.setPrecomputedKNN(knnGraph.knnIndexes, knnGraph.knnDistances);
      }

      // needed so that garbage collector can free memory from distance matrix
      await new Promise<void>((resolve) => {
        setTimeout(() => {
          resolve();
        }, 300);
      });

      let embedding: number[][] | Float32Array[] | undefined |null = null;
      console.time('fit');
      try {
        if (this.useWebGPU && this.webGPUReducer) {
          embedding = await this.webGPUReducer.fit();
          if (!embedding)
            throw new Error('Failed to compute embedding');
          embedding = // transpose TODO: do a better job here
            new Array(embedding[0].length).fill(null).map((_, i) => new Float32Array(embedding!.map((v) => v[i])));
        }
      } catch (e) {
        console.error(e);
      }
      if (!embedding) {
        if (this.useWebGPU) {
          console.error('WEBGPU UMAP failed, falling back to CPU implementation');
          this.reducer.setPrecomputedKNN(knnGraph.knnIndexes, knnGraph.knnDistances);
        }
        embedding = await this.reducer.fitAsync(this.vectors, (epoc) => {
          if (this.progressFunc)
            this.progressFunc(epoc, this.reducer.getNEpochs(), this.reducer.getEmbedding());
        });
      }
      console.timeEnd('fit');
      return arrayCast2Coordinates(embedding);
      function arrayCast2Coordinates(data: number[][] | Float32Array[]): Coordinates {
        return new Array(data.length).fill(0).map((_, i) => (Vector.from(data[i])));
      }
    }
}

const AvailableReducers = {
  'UMAP': UMAPReducer,
  't-SNE': TSNEReducer,
};

export type KnownMethods = keyof typeof AvailableReducers;

export class MultiColDimReducer {
    private reducer: MultiColumnReducer | undefined;
    /**
   * Creates an instance of DimensionalityReducer.
   * @param {any[]} data Vectors to embed.
   * @param {KnownMethods} method Embedding method to be applied
   * @param {KnownMetrics} metric Distance metric to be computed between each of the vectors.
   * @param {Options} [options] Options to pass to the implementing embedders.
   * @memberof DimensionalityReducer
   */
    constructor(data: Array<any[]>, method: DimReductionMethods, metrics: KnownMetrics[],
      weights: number[], distanceAggregation: DistanceAggregationMethod, options: Options) {
      const measures: DistanceMetric[] = [];
      for (let idx = 0; idx < metrics.length; ++idx) {
        const measure = new Measure(metrics[idx]).getMeasure(options.distanceFnArgs[idx]);
        measures.push(measure);
        let bitArrayLength = 2048;
        for (let i = 0; i < data[idx].length; ++i) {
          if (data[idx][i] && data[idx][i]._length) {
            bitArrayLength = data[idx][i]._length;
            break;
          }
        }
        if (isBitArrayMetric(metrics[idx])) {
          for (let i = 0; i < data[idx].length; ++i) {
            if (data[idx][i] && data[idx][i]._data)
              data[idx][i] = new BitArray(data[idx][i]._data, data[idx][i]._length);
            else
              data[idx][i] = new BitArray(bitArrayLength);
          }
        }
      }

      if (method == 'UMAP') {
        this.reducer = new UMAPReducer({
          data: data,
          distanceFnames: metrics,
          distanceFns: measures,
          distanceFnArgs: options.distanceFnArgs,
          weights: weights,
          aggregationMethod: distanceAggregation,
          ...options
        });
      } else if (method == 't-SNE') {
        this.reducer = new TSNEReducer({
          data: data,
          distanceFnames: metrics,
          distanceFns: measures,
          distanceFnArgs: options.distanceFnArgs,
          weights: weights,
          aggregationMethod: distanceAggregation,
          ...options
        });
      }
    }

    /**
   * Embeds the data given into the two-dimensional space using the chosen method.
   *
   * @param {boolean} transpose Whether to transform coordinates to have columns-first orientation.
   * @param {boolean} parallelDistanceWorkers Whether to use parallel distance computation.
   * @throws {Error} If the embedding method was not found.
   * @return {any} Cartesian coordinate of this embedding and distance matrix where applicable.
   * @memberof DimensionalityReducer
   */
    public async transform(transpose: boolean = false): Promise<Matrix> {
      if (this.reducer === undefined)
        throw new Error('Reducer was not defined.');

      let embedding = await this.reducer.transform();

      if (transpose)
        embedding = transposeMatrix(embedding);

      return embedding;
    }

    /**
   * Returns metrics available by type.
   *
   * @param {AvailableDataTypes} typeName type name
   * @return {string[]} Metric names which expects the given data type
   * @memberof DimensionalityReducer
   */
    static availableMetricsByType(typeName: AvailableDataTypes) {
      return Object.keys(AvailableMetrics[typeName]);
    }

    /**
   * Returns dimensionality reduction methods available.
   *
   * @readonly
   * @memberof DimensionalityReducer
   */
    static get availableMethods() {
      return Object.keys(AvailableReducers);
    }

    /**
   * Returns metrics available.
   *
   * @readonly
   * @memberof DimensionalityReducer
   */
    static get availableMetrics() {
      let ans: string[] = [];
      Object.values(AvailableMetrics).forEach((obj) => {
        const array = Object.values(obj);
        ans = [...ans, ...array];
      });
      return ans;
    }
}
