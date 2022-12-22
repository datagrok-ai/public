import * as umj from 'umap-js';
import {TSNE} from '@keckelt/tsne';
import {
  Options,
  DistanceMetric,
  Coordinates,
  Vector,
  Vectors,
  Matrix,
} from '@datagrok-libraries/utils/src/type-declarations';
import {
  calcDistanceMatrix,
  transposeMatrix,
  assert,
} from '@datagrok-libraries/utils/src/vector-operations';
import {SPEBase, PSPEBase, OriginalSPE} from './spe';
import {Measure, KnownMetrics, AvailableMetrics, isBitArrayMetric, AvailableDataTypes} from './typed-metrics';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {UMAPParameters} from 'umap-js';

/** Abstract dimensionality reducer */
abstract class Reducer {
  protected data: Vectors;

  constructor(options: Options) {
    this.data = options.data;
  }

  /** Embeds the data given into the two-dimensional space.
   * @return {any} Cartesian coordinate of this embedding and distance matrix where applicable. */
  abstract transform(): { [key: string]: Matrix };
}

/** t-SNE dimensionality reduction. */
class TSNEReducer extends Reducer {
  protected reducer: TSNE;
  protected iterations: number;
  protected distanceFn: DistanceMetric;

  /**
   * Creates an instance of TSNEReducer.
   * @param {Options} options Options to pass to the constructor.
   * @memberof TSNEReducer
   */
  constructor(options: Options) {
    super(options);
    this.reducer = new TSNE(options);
    this.iterations = options?.iterations ?? 100;
    this.distanceFn = options.distance;
  }

  /**
   * Embeds the data given into the two-dimensional space using t-SNE method.
   * @return {any} Cartesian coordinate of this embedding and distance matrix where applicable.
   */
  public transform(): { [key: string]: Matrix } {
    const distance = calcDistanceMatrix(this.data, this.distanceFn);
    this.reducer.initDataDist(distance);

    for (let i = 0; i < this.iterations; ++i) {
      this.reducer.step(); // every time you call this, solution gets better
    }
    return {distance: distance, embedding: this.reducer.getSolution()};
  }
}

export type UmapOptions = Options & UMAPParameters;

/**
 * Implements UMAP dimensionality reduction.
 *
 * @class UMAPReducer
 * @extends {Reducer}
 */
class UMAPReducer extends Reducer {
  protected reducer: umj.UMAP;
  protected distanceFn: Function;
  protected vectors: number[][];

  /**
   * Creates an instance of UMAPReducer.
   * @param {Options} options Options to pass to the constructor.
   * @memberof UMAPReducer
   */
  constructor(options: UmapOptions) {
    super(options);

    assert('distanceFn' in options);

    this.distanceFn = options.distanceFn!;
    this.vectors = [];
    options.distanceFn = this._encodedDistance.bind(this);
    if (this.data.length < 15)
      options.nNeighbors = this.data.length - 1;
    this.reducer = new umj.UMAP(options);
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
  protected _encodedDistance(a: number[], b: number[]): number {
    return this.distanceFn(this.data[a[0]], this.data[b[0]]);
  }

  /**
   * Encodes the input data as a vector of indices.
   *
   * @protected
   * @memberof UMAPReducer
   */
  protected _encode() {
    for (let i = 0; i < this.data.length; ++i) {
      this.vectors.push([i]);
    }
  }

  /**
   * Embeds the data given into the two-dimensional space using UMAP method.
   * @return {any} Cartesian coordinate of this embedding.
   */
  public transform(): { [key: string]: Matrix } {
    this._encode();

    const embedding = this.reducer.fit(this.vectors);

    function arrayCast2Coordinates(data: number[][]): Coordinates {
      return new Array(data.length).fill(0).map((_, i) => (Vector.from(data[i])));
    }

    return {embedding: arrayCast2Coordinates(embedding)};
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
  public transform(): { [key: string]: Matrix } {
    const emb = this.reducer.embed(this.data);
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
  public transform(): { [key: string]: Matrix } {
    const emb = this.reducer.embed(this.data);
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
  public transform(): { [key: string]: Matrix } {
    const emb = this.reducer.embed(this.data);
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
    const measure = new Measure(metric).getMeasure();
    let specOptions = {};

    if (isBitArrayMetric(metric)) {
      for (let i = 0; i < data.length; ++i) {
        data[i] = new BitArray(data[i]._data, data[i]._length);
      }
    }

    if (method == 'UMAP') {
      specOptions = {
        ...{data: data},
        ...{distanceFn: measure},
        ...{nEpochs: options?.cycles},
        ...options,
      };
    } else if (method == 't-SNE') {
      specOptions = {
        ...{data: data},
        ...{distance: measure},
        ...{iterations: options?.cycles ?? undefined},
        ...options,
      };
    } else if (method == 'SPE') {
      specOptions = {...{data: data}, ...{distance: measure}, ...options};
    } else {
      specOptions = {...{data: data}, ...{distance: measure}, ...options};
    }
    this.reducer = new AvailableReducers[method](specOptions);
  }

  /**
   * Embeds the data given into the two-dimensional space using the chosen method.
   *
   * @param {boolean} transpose Whether to transform coordinates to have columns-first orientation.
   * @throws {Error} If the embedding method was not found.
   * @return {any} Cartesian coordinate of this embedding and distance matrix where applicable.
   * @memberof DimensionalityReducer
   */
  public transform(transpose: boolean = false): { [key: string]: Matrix } {
    if (this.reducer == undefined) {
      throw new Error('Reducer was not defined.');
    }

    let {embedding, distance} = this.reducer.transform();

    if (transpose) {
      embedding = transposeMatrix(embedding);
    }
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
