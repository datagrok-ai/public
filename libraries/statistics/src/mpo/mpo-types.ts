import * as DG from 'datagrok-api/dg';
import {TemplateFunction} from '../compute-functions/types';

/// An array of [x, y] points representing the desirability line, sorted by x ascending.
export type DesirabilityLine = number[][];

export enum DesirabilityMode {
  Freeform = 'freeform',
  Gaussian = 'gaussian',
  Sigmoid = 'sigmoid',
}

export type MissingValueConfig =
  | { strategy: 'exclude' }
  | { strategy: 'default'; score: number }
  | { strategy: 'skip' };

type BasePropertyDesirability = {
  weight: number; /// 0-1
  weightColumn?: string; /// optional column name to read per-row weights from (0-1)
  missingValues?: MissingValueConfig;
  computeFunction?: TemplateFunction;
}

export enum MpoScale {
  Linear = 'linear',
  Log = 'log',
}

export type NumericalDesirability = BasePropertyDesirability & {
  functionType: 'numerical';
  line: DesirabilityLine;
  inverted?: boolean; /// when true, desirability = 1 - d(x): lower-is-better / avoid-target
  min?: number; /// min value of the property (optional; used for editing the line)
  max?: number; /// max value of the property (optional; used for editing the line)

  scale?: MpoScale;

  mode?: DesirabilityMode;

  /// Gaussian mode parameters
  mean?: number;
  sigma?: number;

  /// Sigmoid mode parameters
  x0?: number;
  k?: number;

  freeformLine?: DesirabilityLine;
}

export type CategoricalDesirability = BasePropertyDesirability & {
  functionType: 'categorical';
  categories: { name: string; desirability: number }[];
}

export type PropertyDesirability = NumericalDesirability | CategoricalDesirability;

export const WEIGHTED_AGGREGATIONS = ['Average', 'Sum', 'Product', 'Geomean', 'Min', 'Max'] as const;
export const WEIGHTED_AGGREGATIONS_LIST: WeightedAggregation[] = [...WEIGHTED_AGGREGATIONS];
export type WeightedAggregation = typeof WEIGHTED_AGGREGATIONS[number];
export const DEFAULT_AGGREGATION: WeightedAggregation = 'Average';

export const DESIRABILITY_PROFILE_TYPE = 'MPO Desirability Profile';
export const CURRENT_MPO_VERSION = 1;

/// A map of desirability lines with their weights
export type DesirabilityProfile = {
  type: typeof DESIRABILITY_PROFILE_TYPE;
  version?: number;
  name: string;
  description: string;
  aggregation?: WeightedAggregation;
  properties: { [key: string]: PropertyDesirability };
}

export type MpoResult = {
  scoreColumn: DG.Column;
  desirabilityColumns?: DG.Column[];
}

/// The raw-data array kinds DG.Column.getRawData() can return (int/categorical → Int32Array, float → Float32/64, etc).
export type MpoTypedArray = ReturnType<DG.Column['getRawData']>;

/// How a row participates in the aggregation (stored in ColumnDesirability.state): Contribute (folded in
/// normally), Skip (dropped from this row's aggregation), Bail (forces the whole row's score to null).
export enum RowState { Contribute = 0, Skip = 1, Bail = 2 }

/// One property's column with everything the per-row loops need, hoisted out of the loop:
/// raw data, null sentinel, category lookup, and weight info. Built once by hoistColumns().
/// For categorical columns valNull is the index of the '' category (DG.Column.isNone(i) ⇔ raw[i] === valNull).
export interface HoistedColumn {
  col: DG.Column;
  template: PropertyDesirability;
  raw: MpoTypedArray;
  valNull: number;
  isCat: boolean;
  cats: string[] | null;
  catScore: Map<string, number> | null;
  staticW: number;
  rawW: MpoTypedArray | null;
  wNull: number;
}

/// Per-column desirability mapping result: D[i] is the desirability score (NaN where the column does not
/// contribute), and state marks each row with a RowState. state is null when every row contributes (the common case).
export interface ColumnDesirability {
  D: Float32Array;
  state: Uint8Array | null;
}

/// One cached column mapping in MpoCalculator: the result plus the line/missing-value fingerprint that validates it.
export interface CacheEntry {
  lineKey: string;
  D: Float32Array;
  state: Uint8Array | null;
}

// Numeric codes for the per-cell aggregation switch in reduceMpo() — cheaper to branch on than strings.
export enum AggregationCode { Sum, Average, Product, Geomean, Min, Max }

export const AGG_CODE: Record<WeightedAggregation, AggregationCode> = {
  Sum: AggregationCode.Sum, Average: AggregationCode.Average, Product: AggregationCode.Product,
  Geomean: AggregationCode.Geomean, Min: AggregationCode.Min, Max: AggregationCode.Max,
};
