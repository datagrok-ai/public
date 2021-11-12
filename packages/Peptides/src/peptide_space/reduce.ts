import * as umj from 'umap-js';
import {SPEBase, PSPEBase, Vectors} from './spe';

type distanceMetric = (v1: any, v2: any) => (number);

abstract class Reducer {
  protected data: any[];

  constructor(options: {[name: string]: any}) {
    this.data = options.data;
  }

  abstract transform(): Vectors;
}

class UMAPReducer extends Reducer {
  protected reducer: umj.UMAP;

  constructor(options: {[name: string]: any}) {
    super(options);
    this.reducer = new umj.UMAP(options);
  }

  /**
   * transform
   * @return {Vectors} Dimensionality reduced vectors.
   */
  public transform(): Vectors {
    return this.reducer.fit(this.data);
  }
}

class SPEReducer extends Reducer {
  protected reducer: SPEBase;

  constructor(options: {[name: string]: any}) {
    super(options);
    this.reducer = new SPEBase(options);
  }

  /**
   * transform
   * @return {Vectors} Dimensionality reduced vectors.
   */
  public transform(): Vectors {
    return this.reducer.embed(this.data);
  }
}

class PSPEReducer extends Reducer {
  protected reducer: PSPEBase;

  constructor(options: {[name: string]: any}) {
    super(options);
    this.reducer = new PSPEBase(options);
  }

  /**
   * transform
   * @return {Vectors} Dimensionality reduced vectors.
   */
  public transform(): Vectors {
    return this.reducer.embed(this.data);
  }
}

export class DimensionalityReducer {
  private reducer: Reducer | undefined;
  private static methods: {[key: string]: any} = {
    UMAP: UMAPReducer,
    SPE: SPEReducer,
    PSPE: PSPEReducer,
  };

  constructor(data: any[], method: string, metric: distanceMetric, options?: {[name: string]: any}) {
    if (!(method in DimensionalityReducer.availableMethods())) {
      throw new Error('The method "'+method+'" is not supported');
    }

    if (method == 'UMAP') {
      this.reducer = new UMAPReducer({...{data: data}, ...{distanceFn: metric}, ...options});
    } else if (method == 'SPE') {
      this.reducer = new SPEReducer({...{data: data}, ...options});
    } else if (method == 'PSPE') {
      this.reducer = new PSPEReducer({...{data: data}, ...options});
    }
  }

  public transform(): Vectors {
    if (this.reducer == undefined) {
      throw new Error('Reducer was not defined.');
    }
    return this.reducer.transform();
  }

  static get availableMethods() {
    return DimensionalityReducer.methods.keys();
  }
}
