/**
 * Worker-safe DataFrame and Column shapes.
 *
 * Mirrors the subset of DG.DataFrame / DG.Column the cost function reads.
 * Built so a LiteColumn is drop-in compatible with DG.Column for the
 * surface fitting-utils.ts:getErrors uses (`.col(name)`, `.getRawData()`,
 * `.stats.{min, max}`).
 */

export type LiteColumnType =
  | 'int' | 'float' | 'double' | 'qnum' | 'string' | 'datetime' | 'bool' | 'bigint';

export type LiteRawData =
  Float64Array | Float32Array | Int32Array | BigInt64Array | string[] | boolean[] | unknown[];

export interface LiteColumnStats {
  readonly min: number;
  readonly max: number;
}

export interface LiteColumn {
  readonly name: string;
  readonly type: LiteColumnType;
  readonly length: number;
  getRawData(): LiteRawData;
  toList(): unknown[];
  get(index: number): unknown;
  readonly stats: LiteColumnStats;
}

export interface LiteColumnList {
  readonly length: number;
  byName(name: string): LiteColumn;
  names(): string[];
  [Symbol.iterator](): Iterator<LiteColumn>;
}

export interface LiteDataFrame {
  readonly columns: LiteColumnList;
  readonly rowCount: number;
  col(name: string): LiteColumn | null;
}
