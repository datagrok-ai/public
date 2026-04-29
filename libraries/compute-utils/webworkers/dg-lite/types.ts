/**
 * Worker-safe DataFrame and Column shapes.
 *
 * A minimal structural subset of DG.DataFrame / DG.Column — `.col(name)`,
 * `.getRawData()`, `.stats.{min, max}`, etc. — small enough for a worker
 * to construct and consume without `datagrok-api` and structurally
 * compatible with DG.Column for code that operates on either.
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
