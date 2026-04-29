/**
 * Arrow IPC bytes → LiteDataFrame deserializer (worker-safe).
 *
 * Parallel implementation of `libraries/arrow/src/index.ts:fromFeather` that
 * builds LiteColumns instead of DG.Columns. Consumes apache-arrow directly;
 * has no `datagrok-api` dependency, so it can run inside a web worker.
 */

import * as arrow from 'apache-arrow';
import {createWorkerDG, FLOAT_NULL} from './dg-shim';
import {LiteColumn, LiteDataFrame} from './types';

const dg = createWorkerDG();

function unpackDictionary(vector: arrow.Vector): unknown[] {
  const codes = new Array(vector.length);
  const ks = vector?.data[vector.data.length - 1]?.dictionary?.toArray();
  let i = 0;
  for (const chunk of vector.data) {
    const nullMap = (chunk.nullBitmap || []) as Uint8Array;
    for (let j = 0; j < chunk.values.length; j++) {
      const ix = chunk.values[j];
      codes[i] = nullMap.length && !(nullMap[j >> 3] & (1 << (j % 8))) ? null : ks[ix];
      i++;
    }
  }
  return codes;
}

function bigIntToInt32OrBigInt(name: string, arr: BigInt64Array | BigUint64Array): LiteColumn {
  const overflow = BigInt(2 ** 31 - 1);
  for (const v of arr) {
    if (v > overflow)
      return dg.Column.fromBigInt64Array(name, arr instanceof BigInt64Array ? arr : new BigInt64Array(arr));
  }
  const out = new Int32Array(arr.length);
  for (let i = 0; i < arr.length; i++) out[i] = Number(arr[i]);
  return dg.Column.fromInt32Array(name, out);
}

function timestampBigIntToMs(value: bigint, unit: number): number {
  switch (unit) {
  case 0: return Number(value) * 1000;
  case 1: return Number(value);
  case 2: return Number(value / BigInt(1000));
  case 3: return Number(value / BigInt(1000000));
  default: throw new Error(`Unsupported time unit: ${unit}`);
  }
}

function arrowVectorToLiteColumn(name: string, vector: arrow.Vector, rowCount: number): LiteColumn {
  let type = vector.type;
  let values: any;

  if (arrow.DataType.isDictionary(type)) {
    type = vector.data[vector.data.length - 1].dictionary?.type;
    values = unpackDictionary(vector);
  } else
    values = vector.toArray();


  switch (type.typeId) {
  case arrow.Type.Int8:
  case arrow.Type.Int16:
  case arrow.Type.Int32:
  case arrow.Type.Int:
    if (ArrayBuffer.isView(values)) {
      if ((type as any).bitWidth < 64)
        return dg.Column.fromInt32Array(name, values as Int32Array);
      return bigIntToInt32OrBigInt(name, values as BigInt64Array);
    }
    return dg.Column.fromList('int', name, values);
  case arrow.Type.Uint32:
  case arrow.Type.Int64:
  case arrow.Type.Uint64:
    return bigIntToInt32OrBigInt(name, values);
  case arrow.Type.Float:
  case arrow.Type.Decimal:
    if (ArrayBuffer.isView(values)) {
      if ((type as any).bitWidth < 64)
        return dg.Column.fromFloat32Array(name, values as Float32Array);
      return dg.Column.fromFloat64Array(name, values as Float64Array);
    }
    return dg.Column.fromList('float', name, values);
  case arrow.Type.Utf8:
  case arrow.Type.Interval:
    return dg.Column.fromList('string', name, values);
  case arrow.Type.Bool: {
    // Arrow's Bool vector packs booleans into bytes. toArray() returns a
    // Uint8Array bitset; expand to boolean[] for the LiteColumn.
    const bools: boolean[] = new Array(rowCount);
    if (ArrayBuffer.isView(values)) {
      const bytes = values as Uint8Array;
      for (let i = 0; i < rowCount; ++i)
        bools[i] = (bytes[i >> 3] & (1 << (i & 7))) !== 0;
    } else
      for (let i = 0; i < rowCount; ++i) bools[i] = Boolean((values as ArrayLike<unknown>)[i]);

    return dg.Column.fromList('bool', name, bools);
  }
  case arrow.Type.Date:
  case arrow.Type.Timestamp: {
    // Build directly in DG raw layout (microseconds since epoch) and hand to
    // the dt factory — bypasses fromList's ms→μs conversion so we don't
    // have to round-trip through milliseconds. Preserves FLOAT_NULL too.
    const buf = new Float64Array(rowCount);
    if (values instanceof BigInt64Array)
      for (let i = 0; i < rowCount; ++i) buf[i] = timestampBigIntToMs(values[i], (type as any).unit) * 1000;
    else {
      for (let i = 0; i < rowCount; ++i) {
        const v = (values as any)[i];
        if (v == null) {buf[i] = FLOAT_NULL; continue;}
        const ms = typeof v === 'number' ? v : Number(v?.valueOf?.() ?? v);
        buf[i] = ms * 1000;
      }
    }
    return dg.Column.fromDateTimeMicros(name, buf);
  }
  case arrow.Type.Time:
    if ((type as any)?.bitWidth < 64)
      return dg.Column.fromInt32Array(name, new Int32Array((values as ArrayBufferView).buffer));
    return bigIntToInt32OrBigInt(name, values);
  default:
    return dg.Column.fromList('string', name, Array.from(values as ArrayLike<unknown>));
  }
}

export function arrowIpcToLite(bytes: Uint8Array): LiteDataFrame {
  const table = arrow.tableFromIPC(bytes);
  const cols: LiteColumn[] = [];
  for (let i = 0; i < table.numCols; i++) {
    const v = table.getChildAt(i)!;
    const name = table.schema.fields[i].name;
    cols.push(arrowVectorToLiteColumn(name, v, table.numRows));
  }
  return dg.DataFrame.fromColumns(cols);
}
