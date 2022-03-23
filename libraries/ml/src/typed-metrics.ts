import * as fl from 'fastest-levenshtein';
import {jaroWinkler} from 'jaro-winkler-typescript';

import {DistanceMetric} from '@datagrok-libraries/utils/src/type-declarations';
import {similarityMetric} from '@datagrok-libraries/utils/src/similarity-metrics';
import {calculateEuclideanDistance} from '@datagrok-libraries/utils/src/vector-operations';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {Vector, StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';

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
  },
};

export const MetricToDataType: StringDictionary = Object.keys(AvailableMetrics)
  .reduce((ret: StringDictionary, key) => {
    for (const val of Object.keys(AvailableMetrics[key as AvailableDataTypes])) {
      ret[val as AvailableDataTypes] = key;
    }
    return ret;
  }, {});

export type AvailableDataTypes = keyof typeof AvailableMetrics;
export type StringMetrics = keyof typeof AvailableMetrics['String'];
export type BitArrayMetrics = keyof typeof AvailableMetrics['BitArray'];
export type VectorMetrics = keyof typeof AvailableMetrics['Vector'];
export type KnownMetrics = StringMetrics | BitArrayMetrics | VectorMetrics;

export type ValidTypes = {data: string[], metric: StringMetrics} | {data: Vector[], metric: VectorMetrics} |
                         {data: BitArray[], metric: BitArrayMetrics};

export function isStringMetric(name: KnownMetrics) {
  return MetricToDataType[name] == 'String';
}

export function isBitArrayMetric(name: KnownMetrics) {
  return MetricToDataType[name] == 'BitArray';
}

export function isVectorMetric(name: KnownMetrics) {
  return MetricToDataType[name] == 'Vector';
}

/** Unified class implementing different string measures. */
export class Measure {
  protected method: KnownMetrics;
  protected dataType: AvailableDataTypes;

  /**
   * Creates an instance of Measure with .
   * @param {string} method Method to calculate distance between strings.
   * @memberof Measurer
   */
  constructor(method: KnownMetrics) {
    this.method = method;
    this.dataType = MetricToDataType[method] as AvailableDataTypes;
  }

  /**
   * Returns custom string distance function specified.
   * @return {DistanceMetric} Callback of the measure chosen.
   * @memberof Measurer
   */
  public getMeasure(): DistanceMetric {
    const dict: {[key: string]: {[key2: string]: DistanceMetric}} = AvailableMetrics;
    return dict[this.dataType][this.method];
  }

  /**
   * Returns custom string distance by the given data type.
   * @param {AvailableDataTypes} dataType Metric's data type
   * @return {string[]} Metric names which expects the given data type
   * @memberof Measurer
   */
  public static getMetricByDataType(dataType: AvailableDataTypes): string[] {
    return Object.keys(AvailableMetrics[dataType]);
  }

  /** Returns metric names available.
   * @memberof Measurer
  */
  static get availableMeasures(): string[] {
    return Object.keys(AvailableMetrics);
  }
}
