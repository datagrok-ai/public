import * as fl from 'fastest-levenshtein';
import {jaroWinkler} from 'jaro-winkler-typescript';
import {DistanceMetric} from '@datagrok-libraries/utils/src/type-declarations';
import {
  asymmetricDistance,
  braunBlanquetDistance,
  cosineDistance,
  diceDistance,
  euclideanDistance,
  hammingDistance,
  kulczynskiDistance,
  mcConnaugheyDistance,
  rogotGoldbergDistance,
  russelDistance,
  sokalDistance,
  tanimotoDistance,
} from './distance-metrics-methods';

import {calculateEuclideanDistance} from '@datagrok-libraries/utils/src/vector-operations';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {Vector, StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';

export enum StringMetricsNames {
  Levenshtein = 'Levenshtein',
  JaroWinkler = 'Jaro-Winkler',
  Manhattan = 'Manhattan',
}

export enum VectorMetricsNames {
  Euclidean = 'Euclidean',
}

export enum BitArrayMetricsNames {
  Tanimoto = 'Tanimoto',
  Dice = 'Dice',
  Asymmetric = 'Asymmetric',
  BraunBlanquet = 'Braun-Blanquet',
  Cosine = 'Cosine',
  Kulczynski = 'Kulczynski',
  McConnaughey = 'Mc-Connaughey',
  RogotGoldberg = 'Rogot-Goldberg',
  Russel = 'Russel',
  Sokal = 'Sokal',
  Hamming = 'Hamming',
  Euclidean = 'Euclidean',
}

export const vectorDistanceMetricsMethods: { [name: string]: (x: Vector, y: Vector) => number } = {
  [VectorMetricsNames.Euclidean]: calculateEuclideanDistance,
};

export const stringDistanceMetricsMethods: { [name: string]: (x: string, y: string) => number } = {
  [StringMetricsNames.Levenshtein]: fl.distance,
  [StringMetricsNames.JaroWinkler]: jaroWinkler,
  [StringMetricsNames.Manhattan]: manhattanDistance,
};

export const bitArrayDistanceMetricsMethods: { [name: string]: (x: BitArray, y: BitArray) => number } = {
  [BitArrayMetricsNames.Tanimoto]: tanimotoDistance,
  [BitArrayMetricsNames.Dice]: diceDistance,
  [BitArrayMetricsNames.Asymmetric]: asymmetricDistance,
  [BitArrayMetricsNames.BraunBlanquet]: braunBlanquetDistance,
  [BitArrayMetricsNames.Cosine]: cosineDistance,
  [BitArrayMetricsNames.Kulczynski]: kulczynskiDistance,
  [BitArrayMetricsNames.McConnaughey]: mcConnaugheyDistance,
  [BitArrayMetricsNames.RogotGoldberg]: rogotGoldbergDistance,
  [BitArrayMetricsNames.Russel]: russelDistance,
  [BitArrayMetricsNames.Sokal]: sokalDistance,
  [BitArrayMetricsNames.Hamming]: hammingDistance,
  [BitArrayMetricsNames.Euclidean]: euclideanDistance,
};

export enum AvailableMetricsTypes {
  Vector = 'Vector',
  String = 'String',
  BitArray = 'BitArray'
}

export const AvailableMetrics = {
  [AvailableMetricsTypes.Vector]: {
    [VectorMetricsNames.Euclidean]: vectorDistanceMetricsMethods[VectorMetricsNames.Euclidean],
  },
  [AvailableMetricsTypes.String]: {
    [StringMetricsNames.Levenshtein]: stringDistanceMetricsMethods[StringMetricsNames.Levenshtein],
    [StringMetricsNames.JaroWinkler]: stringDistanceMetricsMethods[StringMetricsNames.JaroWinkler],
    [StringMetricsNames.Manhattan]: stringDistanceMetricsMethods[StringMetricsNames.Manhattan],
  },
  [AvailableMetricsTypes.BitArray]: {
    [BitArrayMetricsNames.Tanimoto]: bitArrayDistanceMetricsMethods[BitArrayMetricsNames.Tanimoto],
    [BitArrayMetricsNames.Dice]: bitArrayDistanceMetricsMethods[BitArrayMetricsNames.Dice],
    [BitArrayMetricsNames.Asymmetric]: bitArrayDistanceMetricsMethods[BitArrayMetricsNames.Asymmetric],
    [BitArrayMetricsNames.BraunBlanquet]: bitArrayDistanceMetricsMethods[BitArrayMetricsNames.BraunBlanquet],
    [BitArrayMetricsNames.Cosine]: bitArrayDistanceMetricsMethods[BitArrayMetricsNames.Cosine],
    [BitArrayMetricsNames.Kulczynski]: bitArrayDistanceMetricsMethods[BitArrayMetricsNames.Kulczynski],
    [BitArrayMetricsNames.McConnaughey]: bitArrayDistanceMetricsMethods[BitArrayMetricsNames.McConnaughey],
    [BitArrayMetricsNames.RogotGoldberg]: bitArrayDistanceMetricsMethods[BitArrayMetricsNames.RogotGoldberg],
    [BitArrayMetricsNames.Russel]: bitArrayDistanceMetricsMethods[BitArrayMetricsNames.Russel],
    [BitArrayMetricsNames.Sokal]: bitArrayDistanceMetricsMethods[BitArrayMetricsNames.Sokal],
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
export type VectorMetrics = keyof typeof AvailableMetrics[AvailableMetricsTypes.Vector];
export type StringMetrics = keyof typeof AvailableMetrics[AvailableMetricsTypes.String];
export type BitArrayMetrics = keyof typeof AvailableMetrics[AvailableMetricsTypes.BitArray];
export type KnownMetrics = StringMetrics | BitArrayMetrics | VectorMetrics;

export type ValidTypes = { data: string[], metric: StringMetrics } | { data: Vector[], metric: VectorMetrics } |
  { data: BitArray[], metric: BitArrayMetrics };

export function isStringMetric(name: KnownMetrics) {
  return MetricToDataType[name] == 'String';
}

export function isBitArrayMetric(name: KnownMetrics) {
  return MetricToDataType[name] == 'BitArray';
}

export function isVectorMetric(name: KnownMetrics) {
  return MetricToDataType[name] == 'Vector';
}

/** Manhattan distance between two sequences (match - 0, mismatch - 1) normalized for length. */
export function manhattanDistance(s1: string, s2: string): number {
  if (s1.length !== s2.length) {
    return 1;
  } else {
    let dist: number = 0;
    for (let i = 1; i < s1.length; i++) {
      dist += s1[i] == s2[i] ? 0 : 1;
    }
    return dist / s1.length;
  }
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
    const dict: { [key: string]: { [key2: string]: DistanceMetric } } = AvailableMetrics;
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
