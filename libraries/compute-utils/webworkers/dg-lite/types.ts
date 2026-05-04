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

export type LiteColumnFilter = 'numerical' | 'categorical' | LiteColumnType | null;

export interface LiteColumn {
  readonly name: string;
  readonly type: LiteColumnType;
  readonly length: number;
  readonly dataFrame: LiteDataFrame | null;
  readonly isNumerical: boolean;
  readonly isCategorical: boolean;
  readonly isEmpty: boolean;
  readonly min: number;
  readonly max: number;
  getRawData(): LiteRawData;
  toList(): unknown[];
  get(index: number): unknown;
  getString(index: number): string;
  getNumber(index: number): number;
  isNone(index: number): boolean;
  matches(filter: LiteColumnFilter): boolean;
  toString(): string;
  readonly stats: LiteColumnStats;
}

export interface LiteColumnList {
  readonly length: number;
  readonly dataFrame: LiteDataFrame | null;
  byName(name: string): LiteColumn;
  byIndex(index: number): LiteColumn;
  byNames(names: string[]): LiteColumn[];
  names(): string[];
  toList(): LiteColumn[];
  toMap(): Map<string, LiteColumn>;
  contains(name: string): boolean;
  firstWhere(predicate: (col: LiteColumn) => boolean): LiteColumn | undefined;
  readonly all: Iterable<LiteColumn>;
  readonly numerical: Iterable<LiteColumn>;
  readonly categorical: Iterable<LiteColumn>;
  readonly boolean: Iterable<LiteColumn>;
  readonly dateTime: Iterable<LiteColumn>;
  readonly numericalNoDateTime: Iterable<LiteColumn>;
  [Symbol.iterator](): Iterator<LiteColumn>;
}

export interface LiteDataFrame {
  id: string | undefined;
  readonly columns: LiteColumnList;
  readonly rowCount: number;
  col(nameOrIndex: string | number): LiteColumn | null;
  getCol(name: string): LiteColumn;
  get(name: string, idx: number): unknown;
  newId(): void;
  toJson(): Record<string, unknown>[];
}
