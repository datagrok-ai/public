/* eslint-disable max-len */
import {TSNE} from '@keckelt/tsne';
import {
  Options,
  Coordinates,
  Vector,
  Vectors,
  Matrix,
} from '@datagrok-libraries/utils/src/type-declarations';
import {
  transposeMatrix,
  assert,
} from '@datagrok-libraries/utils/src/vector-operations';
import {SPEBase, PSPEBase, OriginalSPE} from './spe';
import {Measure, KnownMetrics, AvailableMetrics,
  isBitArrayMetric, AvailableDataTypes} from './typed-metrics/typed-metrics';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {UMAPParameters, UMAP} from './umap';
import {DistanceMatrix, DistanceMatrixService, distanceMatrixProxy, dmLinearIndex} from './distance-matrix';
import {SparseMatrixService} from './distance-matrix/sparse-matrix-service';
import {getKnnGraph, getKnnGraphFromDM} from './umap/knnGraph';

export type SparseMatrixTransferType = {
  i: Int32Array,
  j: Int32Array,
  distance: Float32Array,
}

export interface IReduceDimensionalityResult {
  distance?: Float32Array;
  sparseMatrix?: Map<number, Map<number, number>>;
  embedding: Matrix;
}

export enum DimReductionMethods{
  UMAP = 'UMAP',
  T_SNE = 't-SNE'
}

export interface IUMAPOptions {
  learningRate?: number;
  nComponents?: number;
  nEpochs?: number;
  nNeighbors?: number;
  spread?: number;
  minDist?: number;
  sparseMatrixThreshold?: number;
  preCalculateDistanceMatrix?: boolean;
  usingSparseMatrix?: boolean;
  sparseMatrix?: SparseMatrixTransferType;
  progressFunc?: (epoc: number, epochsLength: number, embeddings: number[][]) => void;
}

export interface ITSNEOptions {
  epsilon?: number;
  perplexity?: number;
  dim?: number;
}

export interface IDimReductionParam {
  uiName: string;
  value: number | null;
  tooltip: string;
  placeholder?: string;
}

/** Umap uses precalculated distance matrix to save time. though for too much data, memory becomes constraint.
 * if we have 100 000 rows, distance matrix will take ~10gb of memory and probably overflow.
 */
export const MAX_DISTANCE_MATRIX_ROWS = 20000;

export class UMAPOptions {
  learningRate: IDimReductionParam = {uiName: 'Learinig rate', value: 1, tooltip: 'The initial learning rate for the embedding optimization'};
  nComponents: IDimReductionParam = {uiName: 'Components', value: 2, tooltip: 'The number of components (dimensions) to project the data to'};
  nEpochs: IDimReductionParam = {uiName: 'Epochs', value: 0, tooltip: 'The number of epochs to optimize embeddings via SGD. Computed automatically if set to 0'};
  nNeighbors: IDimReductionParam = {uiName: 'Neighbors', value: 15, tooltip: 'The number of nearest neighbors to construct the fuzzy manifold'};
  spread: IDimReductionParam = {uiName: 'Spread', value: 1, tooltip: 'The effective scale of embedded points, used with min distance to control the clumped/dispersed nature of the embedding'};
  minDist: IDimReductionParam = {uiName: 'Min distance', value: 0.1, tooltip: 'The effective minimum distance between embedded points, used with spread to control the clumped/dispersed nature of the embedding'};

  constructor() {};
}

export class TSNEOptions {
  epsilon: IDimReductionParam = {uiName: 'Epsilon', value: 10, tooltip: 'Epsilon is learning rate'};
  perplexity: IDimReductionParam = {uiName: 'Perplexity', value: 30, tooltip: 'Roughly how many neighbors each point influences'};
  dim: IDimReductionParam = {uiName: 'Dimensionality', value: 2, tooltip: 'Dimensionality of the embedding'};

  constructor() {};
}

/** Abstract dimensionality reducer */
abstract class Reducer {
  protected data: Vectors;

  constructor(options: Options) {
    this.data = options.data;
  }

  /** Embeds the data given into the two-dimensional space.
   * @return {any} Cartesian coordinate of this embedding and distance matrix where applicable. */
  abstract transform(parallelDistanceWorkers?: boolean): Promise<IReduceDimensionalityResult>;
}

/** t-SNE dimensionality reduction. */
class TSNEReducer extends Reducer {
  protected reducer: TSNE;
  protected iterations: number;
  protected distanceFname: KnownMetrics;
  protected distanceFn: (a: any, b: any) => number;

  /**
   * Creates an instance of TSNEReducer.
   * @param {Options} options Options to pass to the constructor.
   * @memberof TSNEReducer
   */
  constructor(options: Options) {
    super(options);
    this.reducer = new TSNE(options);
    this.iterations = options?.iterations ?? 100;
    this.distanceFname = options.distanceFname;
    this.distanceFn = options.distanceFn;
  }

  /**
   * Embeds the data given into the two-dimensional space using t-SNE method.\
   * @param {boolean} [parallelDistanceWorkers] Whether to use parallel distance workers.
   * @return {any} Cartesian coordinate of this embedding and distance matrix where applicable.
   */
  public async transform(parallelDistanceWorkers?: boolean): Promise<IReduceDimensionalityResult> {
    const distance = parallelDistanceWorkers ? await (async () => {
      const matrixService = new DistanceMatrixService(true, false);
      try {
        const dist = await matrixService.calc(this.data, this.distanceFname);
        matrixService.terminate();
        return dist;
      } catch (e) {
        matrixService.terminate();
        throw e;
      }
    })() :
      (() => { const ret = DistanceMatrix.calc(this.data, (a, b) => this.distanceFn(a, b)); ret.normalize(); return ret.data; })();

    const matrixProxy = distanceMatrixProxy(distance, this.data.length);
    this.reducer.initDataDist(matrixProxy);

    for (let i = 0; i < this.iterations; ++i)
      this.reducer.step(); // every time you call this, solution gets better

    return {distance: distance, embedding: this.reducer.getSolution()};
  }
}

export type UmapOptions = Options & UMAPParameters & {
  preCalculateDistanceMatrix?: boolean,
  usingSparseMatrix?: boolean,
  sparseMatrixThreshold?: number,
  sparseMatrix?: SparseMatrixTransferType,
  progressFunc?: (epoc: number, epochsLength: number, embeddings: number[][]) => void,
};

/**
 * Implements UMAP dimensionality reduction.
 *
 * @class UMAPReducer
 * @extends {Reducer}
 */
class UMAPReducer extends Reducer {
  protected reducer: UMAP;
  protected distanceFname: KnownMetrics;
  protected distanceFn: Function;
  protected vectors: number[];
  protected distanceMatrix?: Float32Array;
  protected usingDistanceMatrix: boolean;
  protected sparseMatrix?: Map<number, Map<number, number>>;
  protected dmIndexFunc: (i: number, j: number) => number;
  protected usingSparseMatrix: boolean;
  protected sparseMatrixThreshold: number;
  protected transferedSparseMatrix?: SparseMatrixTransferType;
  protected progressFunc?: (epoc: number, epochsLength: number, embeddings: number[][]) => void;
  protected distanceFnArgs?: {[_: string]: any};
  /**
   * Creates an instance of UMAPReducer.
   * @param {Options} options Options to pass to the constructor.
   * @memberof UMAPReducer
   */
  constructor(options: UmapOptions) {
    super(options);
    assert('distanceFname' in options);
    assert('distanceFn' in options);
    this.distanceFnArgs = options?.distanceFnArgs;
    this.distanceFn = options.distanceFn!;
    this.usingSparseMatrix = !!options.usingSparseMatrix || !!options.sparseMatrix;
    this.sparseMatrixThreshold = options.sparseMatrixThreshold ?? 0.8;
    this.transferedSparseMatrix = options.sparseMatrix;
    this.progressFunc = options.progressFunc;

    this.distanceFname = options.distanceFname!;
    this.dmIndexFunc = dmLinearIndex(this.data.length);
    //Umap uses vector indexing, so we need to create an array of vectors as indeces.
    this.vectors = new Array(this.data.length).fill(0).map((_, i) => i);
    this.usingDistanceMatrix = !((!options.preCalculateDistanceMatrix && this.data.length > MAX_DISTANCE_MATRIX_ROWS)
      || this.usingSparseMatrix);
    if (this.usingDistanceMatrix)
      options.distanceFn = this._encodedDistanceMatrix.bind(this);
    else if (this.usingSparseMatrix)
      options.distanceFn = this._encodedSparseMatrix.bind(this);
    else
      options.distanceFn = this._encodedDistance.bind(this);

    if (this.data.length < 15)
      options.nNeighbors = this.data.length - 1;
    this.reducer = new UMAP(options);
    // this.reducer.distanceFn = this._encodedDistance.bind(this);
  }

  /**
   * Custom distance wrapper to have numeric inputs instead of string ones.
   *
   * @protected
   * @param {number[]} a The first item.
   * @param {number[]} b The first item.
   * @return {number} Distance metric.
   * @memberof UMAPReducer
   */
  protected _encodedDistanceMatrix(a: number, b: number): number {
    if (a === b)
      return 0;
    if (a > b)
      return this.distanceMatrix![this.dmIndexFunc(b, a)];
    return this.distanceMatrix![this.dmIndexFunc(a, b)];
  }

  protected _encodedSparseMatrix(a: number, b: number): number {
    return this.sparseMatrix!.get(a)?.get(b) ?? this.sparseMatrix!.get(b)?.get(a) ?? 1;
  }

  protected _encodedDistance(a: number, b: number): number {
    return this.distanceFn(this.data[a], this.data[b]);
  }

  /**
   * Embeds the data given into the two-dimensional space using UMAP method.
   * @param {boolean} [parallelDistanceWorkers] Whether to use parallel distance matrix workers.
   * @return {any} Cartesian coordinate of this embedding.
   */
  public async transform(parallelDistanceWorkers?: boolean): Promise<IReduceDimensionalityResult> {
    if (this.usingDistanceMatrix) {
      this.distanceMatrix = parallelDistanceWorkers ? await (async () => {
        const matrixService = new DistanceMatrixService(true, false);
        try {
          const dist = await matrixService.calc(this.data, this.distanceFname, false, this.distanceFnArgs);
          matrixService.terminate();
          return dist;
        } catch (e) {
          matrixService.terminate();
          throw e;
        }
      })() :
        (() => { 
          const ret = DistanceMatrix.calc(this.data, (a, b) => {
            const d = this.distanceFn(a, b);
            return d;}); 
          return ret.data; 
        })();
        if (this.distanceMatrix && !isBitArrayMetric(this.distanceFname)) {
          const knnRes = getKnnGraphFromDM(this.distanceMatrix, this.reducer.neighbors, this.data.length);
          this.reducer.setPrecomputedKNN(knnRes.knnIndexes, knnRes.knnDistances);
        }
        
    } else if (this.usingSparseMatrix) {
          console.time('sparse matrix');
          // let res: {[K in keyof SparseMatrixTransferType]: SparseMatrixTransferType[K] | null} | null
          //   = this.transferedSparseMatrix ??
          //     await new SparseMatrixService().calc(this.data, this.distanceFname, this.sparseMatrixThreshold, this.distanceFnArgs);
          // console.timeEnd('sparse matrix');
          //if (!isBitArrayMetric(this.distanceFname)) {
          //const knnRes = getKnnGraph(res.i!, res.j!, res.distance!, this.reducer.neighbors, this.data.length);
          const knnRes = await new SparseMatrixService().getKNN(this.data, this.distanceFname, this.reducer.neighbors, this.distanceFnArgs);
          console.timeEnd('sparse matrix');
          this.reducer.setPrecomputedKNN(knnRes.knnIndexes, knnRes.knnDistances);
          //}
          console.time('sparse matrix to map')
          this.sparseMatrix = new Map<number, Map<number, number>>();
          // for (let i = 0; i < res.i!.length; ++i) {
          //   const first = res.i![i];
          //   const second = res.j![i];
          //   const distance = res.distance![i];
          //   if (!this.sparseMatrix.has(first))
          //     this.sparseMatrix.set(first, new Map<number, number>());
          //   this.sparseMatrix.get(first)!.set(second, distance);
          // }
          for (let i = 0; i < knnRes.knnIndexes.length; ++i) {
            const first = i;
            if (!this.sparseMatrix.has(first))
              this.sparseMatrix.set(first, new Map<number, number>());
            for (let j = 0; j < knnRes.knnIndexes[i].length; ++j) {
              const second = knnRes.knnIndexes[i][j];
              if (!this.sparseMatrix.has(second))
                this.sparseMatrix.set(second, new Map<number, number>());
              const distance = knnRes.knnDistances[i][j];
              this.sparseMatrix.get(first)!.set(second, distance);
              this.sparseMatrix.get(second)!.set(first, distance);
            }
          }
          console.timeEnd('sparse matrix to map');
          // res.distance = null;
          // res.i = null;
          // res.j = null;
          // res = null;
          // needed so that garbage collector can free memory from distance matrix
          await new Promise<void>((resolve) => {
            setTimeout(() => {
              resolve();
            }, 500);
          })
        }

    
    const embedding = await this.reducer.fitAsync(this.vectors, (epoc) => {
      if (this.progressFunc)
        this.progressFunc(epoc, this.reducer.getNEpochs(), this.reducer.getEmbedding());
    });

    function arrayCast2Coordinates(data: number[][]): Coordinates {
      return new Array(data.length).fill(0).map((_, i) => (Vector.from(data[i])));
    }

    return {embedding: arrayCast2Coordinates(embedding), ...(this.distanceMatrix ? {distance: this.distanceMatrix} : {})};
  }
}

/**
 * Implements original SPE dimensionality reduction.
 *
 * @class SPEReducer
 * @extends {Reducer}
 */
class SPEReducer extends Reducer {
  protected reducer: SPEBase;

  /**
   * Creates an instance of SPEReducer.
   * @param {Options} options Options to pass to the constructor.
   * @memberof SPEReducer
   */
  constructor(options: Options) {
    super(options);
    this.reducer = new SPEBase(options);
  }

  /**
   * Embeds the data given into the two-dimensional space using the original SPE method.
   * @return {any} Cartesian coordinate of this embedding and distance matrix where applicable.
   */
  public async transform(): Promise<IReduceDimensionalityResult> {
    const emb = await this.reducer.embed(this.data);
    return {distance: this.reducer.distance, embedding: emb};
  }
}

/**
 * Implements modified SPE dimensionality reduction.
 *
 * @class PSPEReducer
 * @extends {Reducer}
 */
class PSPEReducer extends Reducer {
  protected reducer: PSPEBase;

  /**
   * Creates an instance of PSPEReducer.
   * @param {Options} options Options to pass to the constructor.
   * @memberof PSPEReducer
   */
  constructor(options: Options) {
    super(options);
    this.reducer = new PSPEBase(options);
  }

  /**
   * Embeds the data given into the two-dimensional space using the modified SPE method.
   * @return {any} Cartesian coordinate of this embedding and distance matrix where applicable.
   */
  public async transform(): Promise<IReduceDimensionalityResult> {
    const emb = await this.reducer.embed(this.data);
    return {distance: this.reducer.distance, embedding: emb};
  }
}

/**
 * Implements original SPE dimensionality reduction.
 *
 * @class OriginalSPEReducer
 * @extends {Reducer}
 */
class OriginalSPEReducer extends Reducer {
  protected reducer: OriginalSPE;

  /**
   * Creates an instance of OriginalSPEReducer.
   * @param {Options} options Options to pass to the constructor.
   * @memberof OriginalSPEReducer
   */
  constructor(options: Options) {
    super(options);
    this.reducer = new OriginalSPE(options);
  }

  /**
   * Embeds the data given into the two-dimensional space using the original SPE method.
   * @return {any} Cartesian coordinate of this embedding and distance matrix where applicable.
   */
  public async transform(): Promise<IReduceDimensionalityResult> {
    const emb = await this.reducer.embed(this.data);
    return {distance: this.reducer.distance, embedding: emb};
  }
}

const AvailableReducers = {
  'UMAP': UMAPReducer,
  't-SNE': TSNEReducer,
  'SPE': SPEReducer,
  'pSPE': PSPEReducer,
  'OriginalSPE': OriginalSPEReducer,
};

export type KnownMethods = keyof typeof AvailableReducers;

/**
 * Unified class implementing different dimensionality reduction methods.
 *
 * @export
 * @class DimensionalityReducer
 */
export class DimensionalityReducer {
  private reducer: Reducer | undefined;

  /**
   * Creates an instance of DimensionalityReducer.
   * @param {any[]} data Vectors to embed.
   * @param {KnownMethods} method Embedding method to be applied
   * @param {KnownMetrics} metric Distance metric to be computed between each of the vectors.
   * @param {Options} [options] Options to pass to the implementing embedders.
   * @memberof DimensionalityReducer
   */
  constructor(data: any[], method: KnownMethods, metric: KnownMetrics, options?: Options) {
    const measure = new Measure(metric).getMeasure(options?.distanceFnArgs);
    let specOptions = {};

    if (isBitArrayMetric(metric)) {
      for (let i = 0; i < data.length; ++i)
        data[i] = new BitArray(data[i]._data, data[i]._length);
    }

    if (method == 'UMAP') {
      specOptions = {
        ...{data: data},
        ...{distanceFn: measure},
        ...{distanceFname: metric},
        ...{nEpochs: options?.cycles},
        ...options,
      };
    } else if (method == 't-SNE') {
      specOptions = {
        ...{data: data},
        ...{distanceFn: measure},
        ...{distanceFname: metric},
        ...{iterations: options?.cycles ?? undefined},
        ...options,
      };
    } else if (method == 'SPE') {
      specOptions = {...{data: data}, ...{distance: measure}, distanceFunctionName: metric, ...options};
    } else {
      specOptions = {...{data: data}, ...{distance: measure}, distanceFunctionName: metric, ...options};
    }
    this.reducer = new AvailableReducers[method](specOptions);
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
  public async transform(transpose: boolean = false, parallelDistanceWorkers?: boolean): Promise<IReduceDimensionalityResult> {
    if (this.reducer === undefined)
      throw new Error('Reducer was not defined.');

    let {embedding, distance} = await this.reducer.transform(parallelDistanceWorkers);

    if (transpose)
      embedding = transposeMatrix(embedding);

    return {distance: distance, embedding: embedding};
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
