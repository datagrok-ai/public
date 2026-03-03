/**
 * Raw typed array as returned by {@link DG.Column.getRawData}.
 *
 * The concrete type depends on the column:
 * - `Int32Array` — int and string (category indices) columns
 * - `Float32Array` — single-precision float columns
 * - `Float64Array` — double-precision float, datetime, and qnum columns
 * - `Uint32Array` — bool columns (bit array)
 */
export type RawData = Int32Array | Float32Array | Float64Array | Uint32Array;

/**
 * Worker-side column type identifier.
 *
 * Unlike DG.COLUMN_TYPE, this splits the DG `'double'` type into `'float32'` and `'float64'`
 * based on the actual typed array storage returned by {@link DG.Column.getRawData}.
 *
 * @example
 * // A DG column of type 'double' stored as Float32Array becomes 'float32':
 * const wc = toWorkerColumn(floatCol);
 * console.log(wc.type); // 'float32'
 */
export type ColumnType = 'int' | 'float32' | 'float64' | 'string' | 'bool' | 'datetime' | 'qnum' | 'bigint';

/**
 * Plain object with all numerical fields from {@link DG.Stats}, plus the null sentinel.
 *
 * All values are eagerly extracted from the Dart-backed Stats object,
 * making this safe to transfer via `postMessage`.
 *
 * @example
 * // Access stats inside a web worker:
 * onmessage = (e) => {
 *   const col: WorkerColumn = e.data;
 *   if (col.stats.missingValueCount > 0)
 *     console.log(`Column has ${col.stats.missingValueCount} missing values`);
 *   // Check for nulls in raw data:
 *   const nullVal = col.stats.nullValue;
 *   for (let i = 0; i < col.length; i++)
 *     if (col.rawData[i] === nullVal) { // handle null }
 * };
 */
export interface WorkerColumnStats {
  /** Total number of values (including missing values). */
  totalCount: number;
  /** Number of missing (empty) values. */
  missingValueCount: number;
  /** Number of unique values. */
  uniqueCount: number;
  /** Number of non-empty values. */
  valueCount: number;
  /** Minimum value. */
  min: number;
  /** Maximum value. */
  max: number;
  /** Sum of all values. */
  sum: number;
  /** Average (mean). */
  avg: number;
  /** Standard deviation. */
  stdev: number;
  /** Variance. */
  variance: number;
  /** Skewness. */
  skew: number;
  /** Kurtosis. */
  kurt: number;
  /** Median value. */
  med: number;
  /** First quartile. */
  q1: number;
  /** Second quartile. */
  q2: number;
  /** Third quartile. */
  q3: number;
  /**
   * Null sentinel value used in {@link WorkerColumn.rawData} for this column type.
   *
   * - `INT_NULL` (-2147483648) for int, string, and bool columns
   * - `FLOAT_NULL` (2.6789344063684636e-34) for float32, float64, datetime, and qnum columns
   */
  nullValue: number;
}

/**
 * Structured-cloneable representation of a {@link DG.Column} for web worker transfer.
 *
 * Contains the raw typed array, pre-computed stats, column metadata, and (for string columns)
 * the categories array. All fields are plain JS values safe for `postMessage`.
 *
 * @example
 * // Main thread — send column to worker:
 * const wc: WorkerColumn = toWorkerColumn(df.col('age'));
 * worker.postMessage(wc);
 *
 * // Worker — receive and use:
 * onmessage = (e) => {
 *   const wc: WorkerColumn = e.data;
 *   const raw = wc.rawData as Float32Array;
 *   const mean = wc.stats.avg;
 *   const nullVal = wc.stats.nullValue;
 *   for (let i = 0; i < wc.length; i++) {
 *     if (raw[i] !== nullVal)
 *       raw[i] -= mean; // center the data
 *   }
 *   postMessage(wc);
 * };
 */
export interface WorkerColumn {
  /** Column name. */
  name: string;
  /** Worker-side column type (see {@link ColumnType}). */
  type: ColumnType;
  /** Number of rows. */
  length: number;
  /** Raw typed array extracted via {@link DG.Column.getRawData}. */
  rawData: RawData;
  /** Pre-computed descriptive statistics. */
  stats: WorkerColumnStats;
  /** String categories — present only for string columns (type === 'string'). */
  categories?: string[];
}

/**
 * Structured-cloneable representation of a {@link DG.DataFrame} for web worker transfer.
 *
 * @example
 * // Main thread — send entire table to worker:
 * const wdf: WorkerDataFrame = toWorkerDataFrame(table);
 * worker.postMessage(wdf);
 *
 * // Worker — iterate columns:
 * onmessage = (e) => {
 *   const wdf: WorkerDataFrame = e.data;
 *   for (const col of wdf.columns) {
 *     console.log(`${col.name}: ${col.type}, ${wdf.rowCount} rows`);
 *   }
 * };
 */
export interface WorkerDataFrame {
  /** DataFrame name. */
  name: string;
  /** Number of rows. */
  rowCount: number;
  /** Columns in order. */
  columns: WorkerColumn[];
}
