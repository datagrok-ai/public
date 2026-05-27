import * as DG from 'datagrok-api/dg';
import type {PropertyDesirability, WeightedAggregation} from './mpo';

export type MpoTypedArray = Int32Array | Float32Array | Float64Array | Uint32Array;

/// One property's column with everything the per-row loops need, hoisted out of the loop:
/// raw data, null sentinel, category lookup, and weight info. Built once by hoistColumns().
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

/// Per-column desirability mapping result: D[i] is the desirability score (NaN where the column does
/// not contribute), and state marks each row — 0 contributes, 1 skip (drop from the row's aggregation),
/// 2 bail (force the row's score to null). state is null when every row contributes (the common case).
export interface ColumnDesirability {
  D: Float32Array;
  state: Uint8Array | null;
}

// Numeric codes for the per-cell aggregation switch in reduceMpo() — cheaper to branch on than strings.
export enum AggregationCode { Sum, Average, Product, Geomean, Min, Max }

export const AGG_CODE: Record<WeightedAggregation, AggregationCode> = {
  Sum: AggregationCode.Sum, Average: AggregationCode.Average, Product: AggregationCode.Product,
  Geomean: AggregationCode.Geomean, Min: AggregationCode.Min, Max: AggregationCode.Max,
};
