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
  case DG.COLUMN_TYPE.STRING:
    return vectorFromArray(target.toList().map((v) => v === '' ? null : v), new Utf8());
  case DG.COLUMN_TYPE.BOOL:
    return vectorFromArray(target.toList(), new Bool());
  default:
    throw new Error(`Crux: target "${target.name}" is ${target.type}; numeric, string or boolean expected`);
  }
}

/** Arrow IPC stream with the features as `x: List<Float64>` and the target as `y` */
export function toArrowIpc(features: DG.Column[], target?: DG.Column): Uint8Array {
  const n = features[0].length;
  const m = features.length;
  const valueOffsets = Int32Array.from({length: n + 1}, (_, i) => i * m);
  const x = makeVector(makeData({
    type: new List(new Field('item', new Float64(), true)),
    length: n,
    nullCount: 0,
    valueOffsets,
    child: makeData({type: new Float64(), length: n * m, nullCount: 0, data: featureMatrix(features)}),
  }));

  const columns: Record<string, Vector> = {x};
  if (target !== undefined)
    columns.y = targetVector(target);
  return tableToIPC(new Table(columns), 'stream');
}

/** Float column from a numeric Arrow column; nulls and non-finite values become DG nulls */
export function floatColumn(table: Table, arrowName: string, name: string): DG.Column {
  const vector = table.getChild(arrowName)!;
  const n = table.numRows;
  const out = new Float32Array(n);

  for (let i = 0; i < n; ++i) {
    const v = vector.isValid(i) ? Number(vector.get(i)) : NaN;
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
