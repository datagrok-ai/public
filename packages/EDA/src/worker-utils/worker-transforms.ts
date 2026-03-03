import * as DG from 'datagrok-api/dg';
import {WorkerColumn, WorkerDataFrame, ColumnType} from './worker-defs';
import {getNullValue} from '../missing-values-imputation/knn-imputer';

/**
 * Resolves the worker {@link ColumnType} from a {@link DG.Column},
 * splitting DG `'double'` into `'float32'` or `'float64'` based on the actual typed array.
 *
 * @param col - source Datagrok column
 * @returns worker-side column type string
 */
function resolveColumnType(col: DG.Column): ColumnType {
  if (col.type === DG.COLUMN_TYPE.FLOAT)
    return col.getRawData() instanceof Float64Array ? 'float64' : 'float32';
  return col.type as ColumnType;
}

/**
 * Extracts a plain {@link WorkerColumn} from a {@link DG.Column}.
 *
 * Copies all stat values eagerly and retrieves the raw typed array via `getRawData()`.
 * The result is structured-cloneable and can be sent to a web worker via `postMessage`.
 *
 * @param col - source Datagrok column
 * @returns plain object safe for worker transfer
 *
 * @example
 * // Send a single column to a web worker:
 * const worker = new Worker('my-worker.ts');
 * const wc = toWorkerColumn(df.col('age'));
 * worker.postMessage(wc);
 */
export function toWorkerColumn(col: DG.Column): WorkerColumn {
  const type = resolveColumnType(col);
  const s = col.stats;

  const result: WorkerColumn = {
    name: col.name,
    type: type,
    length: col.length,
    rawData: col.getRawData(),
    stats: {
      totalCount: s.totalCount,
      missingValueCount: s.missingValueCount,
      uniqueCount: s.uniqueCount,
      valueCount: s.valueCount,
      min: s.min,
      max: s.max,
      sum: s.sum,
      avg: s.avg,
      stdev: s.stdev,
      variance: s.variance,
      skew: s.skew,
      kurt: s.kurt,
      med: s.med,
      q1: s.q1,
      q2: s.q2,
      q3: s.q3,
      nullValue: getNullValue(col),
    },
  };

  if (type === 'string')
    result.categories = col.categories;

  return result;
}

/**
 * Extracts a plain {@link WorkerDataFrame} from a {@link DG.DataFrame}.
 *
 * Converts each column via {@link toWorkerColumn}. The result is structured-cloneable
 * and can be sent to a web worker via `postMessage`.
 *
 * @param df - source Datagrok dataframe
 * @returns plain object safe for worker transfer
 *
 * @example
 * // Send an entire table to a web worker:
 * const worker = new Worker('my-worker.ts');
 * const wdf = toWorkerDataFrame(grok.shell.t);
 * worker.postMessage(wdf);
 */
export function toWorkerDataFrame(df: DG.DataFrame): WorkerDataFrame {
  const columns: WorkerColumn[] = [];
  for (const col of df.columns)
    columns.push(toWorkerColumn(col));
  return {name: df.name, rowCount: df.rowCount, columns: columns};
}

/**
 * Reconstructs a {@link DG.DataFrame} from a {@link WorkerDataFrame} received from a web worker.
 *
 * Creates columns via {@link fromWorkerColumn} and restores the dataframe name.
 *
 * @param wdf - worker dataframe received via `postMessage`
 * @returns Datagrok DataFrame
 *
 * @example
 * // Receive results from a web worker:
 * worker.onmessage = (e) => {
 *   const df = fromWorkerDataFrame(e.data as WorkerDataFrame);
 *   grok.shell.addTableView(df);
 * };
 */
export function fromWorkerDataFrame(wdf: WorkerDataFrame): DG.DataFrame {
  const df = DG.DataFrame.fromColumns(wdf.columns.map(fromWorkerColumn));
  df.name = wdf.name;
  return df;
}

/**
 * Reconstructs a {@link DG.Column} from a {@link WorkerColumn} received from a web worker.
 *
 * Selects the appropriate `DG.Column.fromXxxArray` factory based on {@link WorkerColumn.type}.
 *
 * @param wc - worker column received via `postMessage`
 * @returns Datagrok Column
 *
 * @example
 * // Receive a single column from a web worker:
 * worker.onmessage = (e) => {
 *   const col = fromWorkerColumn(e.data as WorkerColumn);
 *   table.columns.add(col);
 * };
 */
export function fromWorkerColumn(wc: WorkerColumn): DG.Column {
  switch (wc.type) {
  case 'int':
    return DG.Column.fromInt32Array(wc.name, wc.rawData as Int32Array, wc.length);
  case 'float32':
    return DG.Column.fromFloat32Array(wc.name, wc.rawData as Float32Array, wc.length);
  case 'float64':
    return DG.Column.fromFloat64Array(wc.name, wc.rawData as Float64Array, wc.length);
  case 'qnum':
  case 'datetime':
    return DG.Column.fromFloat64Array(wc.name, wc.rawData as Float64Array, wc.length);
  case 'string':
    return DG.Column.fromIndexes(wc.name, wc.categories!, wc.rawData as Int32Array);
  default:
    return DG.Column.fromFloat32Array(wc.name, wc.rawData as Float32Array, wc.length);
  }
}
