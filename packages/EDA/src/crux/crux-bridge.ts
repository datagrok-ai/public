import * as DG from 'datagrok-api/dg';

import {Bool, Field, Float64, List, Table, Utf8, Vector, makeData, makeVector, tableToIPC,
  vectorFromArray} from 'apache-arrow';

function nullSentinel(col: DG.Column): number {
  if (col.type === DG.COLUMN_TYPE.INT)
    return DG.INT_NULL;
  if (col.type === DG.COLUMN_TYPE.FLOAT)
    return DG.FLOAT_NULL;
  throw new Error(`Crux: feature "${col.name}" is ${col.type}; only numeric features are supported`);
}

/** Row-major feature matrix of numeric columns; DG nulls become NaN (a missing value in Crux) */
function featureMatrix(columns: DG.Column[]): Float64Array {
  const n = columns[0].length;
  const m = columns.length;
  const data = new Float64Array(n * m);

  for (let j = 0; j < m; ++j) {
    const col = columns[j];
    const raw = col.getRawData().subarray(0, n);
    const none = nullSentinel(col);
    for (let i = 0; i < n; ++i)
      data[i * m + j] = raw[i] === none ? NaN : raw[i];
  }

  return data;
}

function targetVector(target: DG.Column): Vector {
  switch (target.type) {
  case DG.COLUMN_TYPE.INT:
  case DG.COLUMN_TYPE.FLOAT:
    return makeVector(featureMatrix([target]));
  case DG.COLUMN_TYPE.STRING: {
    const categories = target.categories;
    const raw = target.getRawData();
    const none = categories.indexOf('');
    return vectorFromArray(Array.from({length: target.length}, (_, i) => raw[i] === none ? null : categories[raw[i]]),
      new Utf8());
  }
  case DG.COLUMN_TYPE.BOOL:
    return vectorFromArray(target.toList(), new Bool());
  default:
    throw new Error(`Crux: target "${target.name}" is ${target.type}; numeric, string or boolean expected`);
  }
}

/** Arrow IPC stream with a row-major `rows` x `cols` matrix as `x: List<Float64>` and the target as `y` */
export function matrixIpc(data: Float64Array, rows: number, cols: number, target?: DG.Column): Uint8Array {
  const valueOffsets = Int32Array.from({length: rows + 1}, (_, i) => i * cols);
  const x = makeVector(makeData({
    type: new List(new Field('item', new Float64(), true)),
    length: rows,
    nullCount: 0,
    valueOffsets,
    child: makeData({type: new Float64(), length: rows * cols, nullCount: 0, data}),
  }));

  const columns: Record<string, Vector> = {x};
  if (target !== undefined)
    columns.y = targetVector(target);
  return tableToIPC(new Table(columns), 'stream');
}

/** Arrow IPC stream with the features as `x: List<Float64>` and the target as `y` */
export function toArrowIpc(features: DG.Column[], target?: DG.Column): Uint8Array {
  return matrixIpc(featureMatrix(features), features[0].length, features.length, target);
}

function asFloat64(values: ArrayLike<number>): Float64Array {
  return values instanceof Float64Array ? values : Float64Array.from(values);
}

/** Values of a one-row `List<Float64>` column */
export function flatList(table: Table, arrowName: string): Float64Array {
  return asFloat64(table.getChild(arrowName)!.get(0).toArray());
}

/** Values of a `List<Float64>` column with one list per row, concatenated */
export function flatBatch(table: Table, arrowName: string): Float64Array {
  return asFloat64(table.getChild(arrowName)!.getChildAt(0)!.toArray());
}

/** One FLOAT column per name from a `List<Float64>` column with one list per row (the lists as wide as `names`);
 * a null row is null in every column */
export function listColumns(table: Table, arrowName: string, names: string[]): DG.Column[] {
  const vector = table.getChild(arrowName)!;
  const n = table.numRows;
  const k = names.length;
  const out = names.map(() => new Float32Array(n));
  const flat = vector.nullCount === 0 ? flatBatch(table, arrowName) : null;

  for (let i = 0; i < n; ++i) {
    if (flat === null && !vector.isValid(i)) {
      for (const column of out)
        column[i] = DG.FLOAT_NULL;
      continue;
    }
    const row = flat === null ? vector.get(i)! : null;
    for (let j = 0; j < k; ++j) {
      const v = row === null ? flat![i * k + j] : Number(row.get(j));
      out[j][i] = Number.isFinite(v) ? v : DG.FLOAT_NULL;
    }
  }

  return names.map((name, j) => DG.Column.fromFloat32Array(name, out[j]));
}

/** Float column from a numeric Arrow column; nulls and non-finite values become DG nulls */
export function floatColumn(table: Table, arrowName: string, name: string): DG.Column {
  const vector = table.getChild(arrowName)!;
  const n = table.numRows;
  const out = new Float32Array(n);
  const values = vector.nullCount === 0 ? vector.toArray() : null;

  for (let i = 0; i < n; ++i) {
    const v = values !== null ? Number(values[i]) : vector.isValid(i) ? Number(vector.get(i)) : NaN;
    out[i] = Number.isFinite(v) ? v : DG.FLOAT_NULL;
  }

  return DG.Column.fromFloat32Array(name, out);
}

export function stringValues(table: Table, arrowName: string): (string | null)[] {
  const vector = table.getChild(arrowName)!;
  const out = new Array<string | null>(table.numRows);
  for (let i = 0; i < out.length; ++i)
    out[i] = vector.isValid(i) ? vector.get(i) : null;
  return out;
}
