import * as fl from 'fastest-levenshtein';
import {jaroWinkler} from 'jaro-winkler-typescript';

import {DistanceMetric} from '@datagrok-libraries/utils/src/type-declarations';
import {similarityMetric} from '@datagrok-libraries/utils/src/similarity-metrics';
import {calculateEuclideanDistance} from '@datagrok-libraries/utils/src/operations';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {Vector} from '@datagrok-libraries/utils/src/type-declarations';

export const AvailableMetrics = {
  'Vector': {
    'EuclideanDistance': calculateEuclideanDistance,
  },
  'String': {
    'Levenshtein': fl.distance,
    'Jaro-Winkler': jaroWinkler,
  },
  'BitArray': {
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
  }
};

export type AvailableDataTypes = keyof typeof AvailableMetrics;
export type StringMetrics = keyof typeof AvailableMetrics['String'];
export type BitArrayMetrics = keyof typeof AvailableMetrics['BitArray'];
export type VectorMetrics = keyof typeof AvailableMetrics['Vector'];
export type KnownMetrics = StringMetrics | BitArrayMetrics | VectorMetrics;

export type ValidTypes = {data: string[], metric: StringMetrics} | {data: Vector[], metric: VectorMetrics} | 
                         {data: BitArray[], metric: BitArrayMetrics};

export function isStringMetric(name: string) {
  return Object.keys(AvailableMetrics['String']).some(metricName => metricName == name);
}

export function isBitArrayMetric(name: string) {
  return Object.keys(AvailableMetrics['BitArray']).some(metricName => metricName == name);
}

export function isVectorMetric(name: string) {
  return Object.keys(AvailableMetrics['Vector']).some(metricName => metricName == name);
}

/** Unified class implementing different string measures. */
export class StringMeasure {
  protected method: KnownMetrics;

  /**
   * Creates an instance of StringMeasure with .
   * @param {string} method Method to calculate distance between strings.
   * @memberof Measurer
   */
  constructor(method: KnownMetrics) {
    this.method = method;
  }

  /**
   * Returns custom string distance function specified.
   * @return {DistanceMetric} Callback of the measure chosen.
   */
  public getMeasure(): DistanceMetric {
    for (const key of Object.keys(AvailableMetrics)) {
      const dict: {[name: string]: (v1: any, v2: any) => (number)} = AvailableMetrics[key as AvailableDataTypes];
      if (this.method in dict) {
        return dict[this.method];
      }
    }
    return calculateEuclideanDistance;
  }

  /**
   * Returns custom string distance function specified.
   * @return {string[]} Callback of the measure chosen.
   */
  public static getMetricByDataType(dataType: AvailableDataTypes): string[] {
    return Object.keys(AvailableMetrics[dataType]);
  }

  /** Returns metric names available. */
  static get availableMeasures(): string[] {
    return Object.keys(AvailableMetrics);
  }
}
