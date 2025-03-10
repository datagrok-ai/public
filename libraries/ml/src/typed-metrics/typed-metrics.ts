import * as fl from 'fastest-levenshtein';
import {jaroWinkler} from 'jaro-winkler-typescript';
import {DistanceMetric} from '@datagrok-libraries/utils/src/type-declarations';
import {
  asymmetricDistance,
  braunBlanquetDistance,
  cosineDistance,
  diceDistance,
  euclideanDistanceBitArray,
  hammingDistance,
  kulczynskiDistance,
  mcConnaugheyDistance,
  rogotGoldbergDistance,
  russelDistance,
  sokalDistance,
  tanimotoDistance,
  numericDistance,
  tanimotoDistanceIntArray,
  inverseCommonItemsCount,
  vectorEuclideanDistance,
  vectorManhattenDistance,
  vectorCosineDistance,
} from '../distance-metrics-methods';

import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {Vector, StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {mmDistanceFunctions, MmDistanceFunctionsNames} from '../macromolecule-distance-functions';
import {DistanceMetricsSubjects, BitArrayMetricsNames,
  StringMetricsNames, VectorMetricsNames, NumberMetricsNames, IntArrayMetricsNames,
  NumberArrayMetricsNames} from './consts';


export const vectorDistanceMetricsMethods: { [name: string]: (x: Vector, y: Vector) => number } = {
  [VectorMetricsNames.Euclidean]: vectorEuclideanDistance,
  [VectorMetricsNames.Manhattan]: vectorManhattenDistance,
  [VectorMetricsNames.Cosine]: vectorCosineDistance,
};

export const stringDistanceMetricsMethods: { [name: string]: (x: string, y: string) => number } = {
  [StringMetricsNames.Levenshtein]: fl.distance,
  [StringMetricsNames.JaroWinkler]: jaroWinkler,
  [StringMetricsNames.Manhattan]: manhattanDistance,
  [StringMetricsNames.Onehot]: categoricalDistance,
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
  [BitArrayMetricsNames.Euclidean]: euclideanDistanceBitArray,
};

export const intArrayDistanceMetricsMethods: { [name: string]: (x: Uint32Array, y: Uint32Array) => number } = {
  [IntArrayMetricsNames.TanimotoIntArray]: tanimotoDistanceIntArray,
};

export const numberDistanceMetricsMethods: { [name: string]: (args: any) => (x: number, y: number) => number } = {
  [NumberMetricsNames.Difference]: numericDistance,
};

export const numberArrayDistanceMetrics:
{ [name: string]: (args: any) => (x: ArrayLike<number>, y: ArrayLike<number>) => number } = {
  [NumberArrayMetricsNames.CommonItems]: inverseCommonItemsCount,
};

export const AvailableMetrics = {
  [DistanceMetricsSubjects.Vector]: {
    [VectorMetricsNames.Euclidean]: vectorDistanceMetricsMethods[VectorMetricsNames.Euclidean],
    [VectorMetricsNames.Manhattan]: vectorDistanceMetricsMethods[VectorMetricsNames.Manhattan],
    [VectorMetricsNames.Cosine]: vectorDistanceMetricsMethods[VectorMetricsNames.Cosine],
  },
  [DistanceMetricsSubjects.String]: {
    [StringMetricsNames.Levenshtein]: stringDistanceMetricsMethods[StringMetricsNames.Levenshtein],
    [StringMetricsNames.JaroWinkler]: stringDistanceMetricsMethods[StringMetricsNames.JaroWinkler],
    [StringMetricsNames.Manhattan]: stringDistanceMetricsMethods[StringMetricsNames.Manhattan],
    [StringMetricsNames.Onehot]: stringDistanceMetricsMethods[StringMetricsNames.Onehot],
  },
  [DistanceMetricsSubjects.BitArray]: {
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
  [DistanceMetricsSubjects.MacroMolecule]: { // optional args needed for macromolecule functions which initialize them
    [MmDistanceFunctionsNames.HAMMING]: mmDistanceFunctions[MmDistanceFunctionsNames.HAMMING],
    [MmDistanceFunctionsNames.LEVENSHTEIN]: mmDistanceFunctions[MmDistanceFunctionsNames.LEVENSHTEIN],
    [MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH]: mmDistanceFunctions[MmDistanceFunctionsNames.NEEDLEMANN_WUNSCH],
    [MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE]:
      mmDistanceFunctions[MmDistanceFunctionsNames.MONOMER_CHEMICAL_DISTANCE],
  },
  [DistanceMetricsSubjects.Number]: {
    [NumberMetricsNames.Difference]: numberDistanceMetricsMethods[NumberMetricsNames.Difference],
  },
  [DistanceMetricsSubjects.IntArray]: {
    [IntArrayMetricsNames.TanimotoIntArray]: intArrayDistanceMetricsMethods[IntArrayMetricsNames.TanimotoIntArray],
  },
  [DistanceMetricsSubjects.NumberArray]: {
    [NumberArrayMetricsNames.CommonItems]: numberArrayDistanceMetrics[NumberArrayMetricsNames.CommonItems],
  },
};

export const MetricToDataType: StringDictionary = Object.keys(AvailableMetrics)
  .reduce((ret: StringDictionary, key) => {
    for (const val of Object.keys(AvailableMetrics[key as AvailableDataTypes]))
      ret[val as AvailableDataTypes] = key;

    return ret;
  }, {});

export type AvailableDataTypes = keyof typeof AvailableMetrics;
export type VectorMetrics = keyof typeof AvailableMetrics[DistanceMetricsSubjects.Vector];
export type StringMetrics = keyof typeof AvailableMetrics[DistanceMetricsSubjects.String];
export type BitArrayMetrics = keyof typeof AvailableMetrics[DistanceMetricsSubjects.BitArray];
export type NumberMetrics = keyof typeof AvailableMetrics[DistanceMetricsSubjects.Number];
export type IntArrayMetrics = keyof typeof AvailableMetrics[DistanceMetricsSubjects.IntArray];
export type NumberArrayMetrics = keyof typeof AvailableMetrics[DistanceMetricsSubjects.NumberArray];
export type KnownMetrics = StringMetrics | BitArrayMetrics | VectorMetrics |
  MmDistanceFunctionsNames | NumberMetricsNames | IntArrayMetricsNames | NumberArrayMetricsNames;

export type ValidTypes =
  { data: string[], metric: StringMetrics | MmDistanceFunctionsNames } |
  { data: Vector[], metric: VectorMetrics } |
  { data: BitArray[], metric: BitArrayMetrics } |
  { data: number[], metric: NumberMetricsNames };

export function isStringMetric(name: KnownMetrics) {
  return MetricToDataType[name] == 'String';
}

export function isBitArrayMetric(name: KnownMetrics) {
  return MetricToDataType[name] == 'BitArray';
}

export function isVectorMetric(name: KnownMetrics) {
  return MetricToDataType[name] == 'Vector';
}

export function isMacroMoleculeMetric(name: KnownMetrics) {
  return MetricToDataType[name] == DistanceMetricsSubjects.MacroMolecule.toString();
}

export function isNumericMetric(name: KnownMetrics) {
  return MetricToDataType[name] == DistanceMetricsSubjects.Number.toString();
}

export function isNumberArrayMetric(name: KnownMetrics) {
  return MetricToDataType[name] == DistanceMetricsSubjects.NumberArray.toString();
}

/** Manhattan distance between two sequences (match - 0, mismatch - 1) normalized for length. */
export function manhattanDistance(s1: string, s2: string): number {
  if (s1.length !== s2.length) {
    return 1;
  } else {
    let dist: number = 0;
    for (let i = 1; i < s1.length; i++)
      dist += s1[i] == s2[i] ? 0 : 1;
    return dist / s1.length;
  }
}

export function categoricalDistance(s1: string, s2: string): number {
  return s1 === s2 ? 0 : 1;
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
   * Returns true if the metric needs arguments to be calculated.
   * @param {KnownMetrics} method Metric to check if it needs arguments.
   * @return {boolean} True if the metric needs arguments.
   * @memberof Measure
   */
  public metricNeedsArgs(method: KnownMetrics): boolean {
    return isMacroMoleculeMetric(method) || isNumericMetric(method) || isNumberArrayMetric(method);
  }
  /**
   * Returns custom string distance function specified.
   * @param {opts} opts Options for the measure. used for macromolecule distances
   * @return {DistanceMetric} Callback of the measure chosen.
   * @memberof Measurer
   */
  public getMeasure(opts?: any): DistanceMetric {
    const dict: { [key: string]:
      {[key2: string]: DistanceMetric | ((opts: any) => DistanceMetric)}
    } = AvailableMetrics;
    if (!dict.hasOwnProperty(this.dataType) || !dict[this.dataType].hasOwnProperty(this.method))
      throw new Error(`Unknown measure ${this.method} for data type ${this.dataType}`);
    return this.metricNeedsArgs(this.method) ?
      (dict[this.dataType][this.method] as ((opts: any) => DistanceMetric))(opts) :
      dict[this.dataType][this.method] as DistanceMetric;
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
