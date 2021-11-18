import * as umj from 'umap-js';
import {TSNE} from '@keckelt/tsne';

import {Options, DistanceMetric, Coordinates, Vector, Vectors} from './type_declarations';
import {calcDistanceMatrix} from './operations';
import {SPEBase, PSPEBase} from './spe';
import {Measurer} from '../../../libraries/utils/src/string_measure';

/**
 * Abstract dimensionality reducer.
 *
 * @abstract
 * @class Reducer
 */
abstract class Reducer {
  protected data: Vectors;

  /**
   * Creates an instance of Reducer.
   * @param {Options} options Options to pass to the constructor.
   * @memberof Reducer
   */
  constructor(options: Options) {
    this.data = options.data;
  }

  /**
   * Is to embed the data given into the two-dimensional space.
   *
   * @abstract
   * @return {Coordinates} Cartesian coordinate of this embedding.
   * @memberof Reducer
   */
  abstract transform(): Coordinates;
}

/**
 * Implements t-SNE dimensionality reduction.
 *
 * @class TSNEReducer
 * @extends {Reducer}
 */
class TSNEReducer extends Reducer {
  protected reducer: TSNE;
  protected iterations: number;
  protected distance: DistanceMetric;

  /**
   * Creates an instance of TSNEReducer.
   * @param {Options} options Options to pass to the constructor.
   * @memberof TSNEReducer
   */
  constructor(options: Options) {
    super(options);
    this.reducer = new TSNE(options);
    this.iterations = options?.iterations ?? 100;
    this.distance = options.distance;
  }

  /**
   * Embeds the data given into the two-dimensional space using t-SNE method.
   * @return {Coordinates} Cartesian coordinate of this embedding.
   */
  public transform(): Coordinates {
    this.reducer.initDataDist(calcDistanceMatrix(this.data, this.distance));

    for (let i = 0; i < this.iterations; ++i) {
      this.reducer.step(); // every time you call this, solution gets better
    }
    return this.reducer.getSolution();
  }
}

/**
 * Implements UMAP dimensionality reduction.
 *
 * @class UMAPReducer
 * @extends {Reducer}
 */
class UMAPReducer extends Reducer {
  protected reducer: umj.UMAP;

  /**
   * Creates an instance of UMAPReducer.
   * @param {Options} options Options to pass to the constructor.
   * @memberof UMAPReducer
   */
  constructor(options: Options) {
    super(options);
    this.reducer = new umj.UMAP(options);
  }

  /**
   * Embeds the data given into the two-dimensional space using UMAP method.
   * @return {Coordinates} Cartesian coordinate of this embedding.
   */
  public transform(): Coordinates {
    const embedding = this.reducer.fit(this.data);

    function arrayCast2Coordinates(data: number[][]): Coordinates {
      return new Array(data.length).map((_, i) => (Vector.from(data[i])));
    }
    
    return arrayCast2Coordinates(embedding);
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
   * @return {Coordinates} Cartesian coordinate of this embedding.
   */
  public transform(): Coordinates {
    return this.reducer.embed(this.data);
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
   * @return {Coordinates} Cartesian coordinate of this embedding.
   */
  public transform(): Coordinates {
    return this.reducer.embed(this.data);
  }
}

/**
 * Unified class implementing different dimensionality reduction methods.
 *
 * @export
 * @class DimensionalityReducer
 */
export class DimensionalityReducer {
  private reducer: Reducer | undefined;
  private methods: string[];
  private measurer: Measurer;

  /**
   * Creates an instance of DimensionalityReducer.
   * @param {any[]} data Vectors to embed.
   * @param {string} method Embedding method to be applied
   * @param {string} metric Distance metric to be computed between each of the vectors.
   * @param {Options} [options] Options to pass to the implementing embedders.
   * @memberof DimensionalityReducer
   */
  constructor(data: any[], method: string, metric: string, options?: Options) {
    this.methods = ['UMAP', 'TSNE', 'SPE', 'PSPE'];

    if (!this.availableMethods.includes(method)) {
      throw new Error('The method "'+method+'" is not supported');
    }

    this.measurer = new Measurer(metric);
    const measure = this.measurer.getMeasure()

    if (method == 'UMAP') {
      this.reducer = new UMAPReducer({...{data: data}, ...{distanceFn: measure}, ...options});
    } else if (method == 'TSNE') {
      this.reducer = new TSNEReducer({
        ...{data: data},
        ...{distance: measure},
        ...{iterations: options?.cycles ?? undefined},
        ...options,
      });
    } else if (method == 'SPE') {
      this.reducer = new SPEReducer({...{data: data}, ...{distance: measure}, ...options});
    } else {
      this.reducer = new PSPEReducer({...{data: data}, ...{distance: measure}, ...options});
    }
  }

  /**
   * Embeds the data given into the two-dimensional space using the chosen method.
   *
   * @throws {Error} If the embedding method was not found.
   * @return {Coordinates} Cartesian coordinate of this embedding.
   * @memberof DimensionalityReducer
   */
  public transform(): Coordinates {
    if (this.reducer == undefined) {
      throw new Error('Reducer was not defined.');
    }
    return this.reducer.transform();
  }

  /**
   * Returns dimensionality reduction methods available.
   *
   * @readonly
   * @memberof DimensionalityReducer
   */
  get availableMethods() {
    return this.methods;
  }
}
