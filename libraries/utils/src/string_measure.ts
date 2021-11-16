import * as fl from 'fastest-levenshtein';
import {jaroWinkler} from 'jaro-winkler-typescript';

import {DistanceMetric} from './type_declarations';

export class Measurer {
    protected method: string;
    protected receipt: {[name: string]: DistanceMetric} = {
      'Levenshtein': fl.distance,
      'Jaro-Winkler': jaroWinkler,
    };

    constructor(method: string) {
      this.method = method;
    }

    /**
     * getMeasure
     * @return {DistanceMetric} Callback of the measure chosen.
     */
    public getMeasure(): DistanceMetric {
      return this.receipt[this.method];
    }
}
