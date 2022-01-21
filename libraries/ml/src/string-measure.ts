import * as fl from 'fastest-levenshtein';
import {jaroWinkler} from 'jaro-winkler-typescript';

import {DistanceMetric} from '@datagrok-libraries/utils/src/type-declarations';
import {similarityMetric} from '@datagrok-libraries/utils/src/similarity-metrics';
import {calculateEuclideanDistance} from '@datagrok-libraries/utils/src/operations';

const AvailableMetrics: {[name: string]: DistanceMetric} = {
  'EuclideanDistance': calculateEuclideanDistance,
  'Levenshtein': fl.distance,
  'Jaro-Winkler': jaroWinkler,
  'Tanimoto': similarityMetric['Tanimoto'],
  'Dice': similarityMetric['Dice'],
  'Asymmetric': similarityMetric['Asymmetric'],
  'Braun-Blanquet': similarityMetric['Braun-Blanquet'],
  'Cosine': similarityMetric['Cosine'],
  'Kulczynski': similarityMetric['Kulczynski'],
  'Mc-Connaughey': similarityMetric['Mc-Connaughey'],
  'Rogot-Goldberg': similarityMetric['Rogot-Goldberg'],
  'Russel': similarityMetric['Russel'],
  'Sokal': similarityMetric['Sokal'],
};

export const MetricDataTypes: {[name: string]: string[]} = {
  'String': ['Levenshtein', 'Jaro-Winkler'],
  'BitArray': Object.keys(similarityMetric),
  'Vector': ['EuclideanDistance'],
  'Number': [],
  'Object': [],
};

export type KnownMetrics = keyof typeof AvailableMetrics;

/** Unified class implementing different string measures. */
export class StringMeasure {
  protected method: string;

  /**
   * Creates an instance of StringMeasure with .
   * @param {string} method Method to calculate distance between strings.
   * @memberof Measurer
   */
  constructor(method: KnownMetrics) {
    this.method = method as string;
  }

  /**
   * Returns custom string distance function specified.
   * @return {DistanceMetric} Callback of the measure chosen.
   */
  public getMeasure(): DistanceMetric {
    return AvailableMetrics[this.method];
  }

  /** Returns metric names available. */
  static get availableMeasures(): string[] {
    return Object.keys(AvailableMetrics);
  }
}
