import * as fl from 'fastest-levenshtein';
import {jaroWinkler} from 'jaro-winkler-typescript';

import {DistanceMetric} from './type-declarations';
import {assert} from './operations';

/**
 * Unified class implementing the different string measures.
 *
 * @export
 * @class Measurer
 */
export class Measurer {
    protected method: string;
    public static receipt: {[name: string]: DistanceMetric} = {
      'Levenshtein': fl.distance,
      'Jaro-Winkler': jaroWinkler,
    };

    /**
     * Creates an instance of Measurer.
     * @param {string} method Method to calculate distance between strings.
     * @memberof Measurer
     */
    constructor(method: string) {
      assert(Measurer.availableMeasures.includes(method), 'The ${method} was not found.')
      this.method = method;
    }

    /**
     * Returns custom string distance function specified.
     * @return {DistanceMetric} Callback of the measure chosen.
     */
    public getMeasure(): DistanceMetric {
      return Measurer.receipt[this.method];
    }

    /**
     * Returns available string distance metrics.
     *
     * @readonly
     * @type {string[]}
     * @memberof Measurer
     */
    public static get availableMeasures() : string[] {
      let keys = [];

      for (let key in Measurer.receipt) {
        keys.push(key);
      }
      return keys;
    }
    
}
