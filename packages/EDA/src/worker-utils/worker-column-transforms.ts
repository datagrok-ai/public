import * as DG from 'datagrok-api/dg';
import {WorkerColumn, ColumnType} from './worker-column-defs';
import {getNullValue} from '../missing-values-imputation/knn-imputer';

/** Resolves the worker ColumnType from a DG.Column, splitting 'double' into 'float32'/'float64'. */
function resolveColumnType(col: DG.Column): ColumnType {
  if (col.type === DG.COLUMN_TYPE.FLOAT)
    return col.getRawData() instanceof Float64Array ? 'float64' : 'float32';
  return col.type as ColumnType;
}

/** Extracts a plain WorkerColumn from a DG.Column. */
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

/** Reconstructs a DG.Column from a WorkerColumn received from a web worker. */
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
