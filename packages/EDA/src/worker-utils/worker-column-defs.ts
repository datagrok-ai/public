/** Raw typed array as returned by DG.Column.getRawData(). */
export type RawData = Int32Array | Float32Array | Float64Array | Uint32Array;

/** Worker-side column type. Splits DG 'double' into 'float32'/'float64' by actual storage. */
export type ColumnType = 'int' | 'float32' | 'float64' | 'string' | 'bool' | 'datetime' | 'qnum' | 'bigint';

/** Plain object with all numerical stats from DG.Stats, plus the null sentinel. */
export interface WorkerColumnStats {
  totalCount: number;
  missingValueCount: number;
  uniqueCount: number;
  valueCount: number;
  min: number;
  max: number;
  sum: number;
  avg: number;
  stdev: number;
  variance: number;
  skew: number;
  kurt: number;
  med: number;
  q1: number;
  q2: number;
  q3: number;
  /** Null sentinel value used in rawData for this column type. */
  nullValue: number;
}

/** Structured-cloneable representation of a DG.Column for web worker transfer. */
export interface WorkerColumn {
  name: string;
  type: ColumnType;
  length: number;
  rawData: RawData;
  stats: WorkerColumnStats;
  /** String categories — present only for string columns. */
  categories?: string[];
}
