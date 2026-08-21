// DataFrame ↔ Arrow IPC (Feather) conversion.
//
// Extracted verbatim from packages/Arrow/src/api/api.ts. Parquet handling
// stays in the @datagrok/arrow package because parquet-wasm's WASM init is
// tied to that package's webRoot.

import * as DG from 'datagrok-api/dg';
import * as arrow from 'apache-arrow';

export function toFeather(table: DG.DataFrame, asStream: boolean = true): Uint8Array | null {
  const arrays: { [key: string]: any[] } = {};
  const fields: arrow.Field[] = [];

  for (const name of table.columns.names()) {
    const column = table.columns.byName(name);
    const columnType = column.type;

    let values: any[] = [];
    let type: arrow.DataType;

    if (columnType === 'double' || columnType === 'qnum') {
      const raw = (column.getRawData() as Float64Array).subarray(0, column.length);
      values = Array.from(raw, (v) => v === DG.FLOAT_NULL ? null : v);
      // always Float64: inferring Int32 from integral values made the written schema drift with data
      type = new arrow.Float64();
    } else if (columnType === 'int') {
      const raw = (column.getRawData() as Int32Array).subarray(0, column.length);
      values = Array.from(raw, (v) => v === DG.INT_NULL ? null : v);
      type = new arrow.Int32();
    } else if (columnType === 'datetime') {
      const raw = (column.getRawData() as Float64Array).subarray(0, column.length);
      values = Array.from(raw, (v) => v === DG.FLOAT_NULL ? null : new Date(v / 1000));
      type = new arrow.TimestampMillisecond();
    } else if (columnType === 'string') {
      values = column.toList();
      type = new arrow.Utf8();
    } else if (columnType === 'bool') {
      values = column.toList();
      type = new arrow.Bool();
    } else {
      const columnLength = column.length;
      const array = new Array(columnLength);
      for (let i = 0; i < columnLength; i++)
        array[i] = column.get(i);
      values = array;
      type = columnType === 'bigint' ? new arrow.Int64() : new arrow.Utf8();
    }

    arrays[name] = values;
    fields.push(new arrow.Field(name, type, true));
  }

  // Build column data positionally, not via an object map. JS object keys
  // that look like non-negative integers (e.g. "0", "25", "100") are
  // iterated by V8 in numeric order ahead of non-numeric keys, which would
  // misalign the table columns against the schema for any DataFrame whose
  // column names are numeric strings.
  const childData: arrow.Data<any>[] = [];
  for (const field of fields) {
    const type = field.type;
    let vec: arrow.Vector<any>;
    if (type instanceof arrow.Float64)
      vec = arrow.vectorFromArray(arrays[field.name] as number[], type);
    else if (type instanceof arrow.Int32)
      vec = arrow.vectorFromArray(arrays[field.name] as number[], type);
    else if (type instanceof arrow.TimestampMillisecond)
      vec = arrow.vectorFromArray(arrays[field.name] as Date[], type);
    else if (type instanceof arrow.Utf8)
      vec = arrow.vectorFromArray(arrays[field.name] as string[], type);
    else if (type instanceof arrow.Int64)
      vec = arrow.vectorFromArray(arrays[field.name] as (bigint | null)[], type);
    else if (type instanceof arrow.Bool)
      vec = arrow.vectorFromArray(arrays[field.name] as boolean[], type);
    else
      throw new Error(`Unsupported type for field ${field.name}`);
    childData.push(vec.data[0]);
  }

  const schema = new arrow.Schema(fields);
  const length = childData.length === 0 ? 0 : childData[0].length;
  const structData = arrow.makeData({
    type: new arrow.Struct(fields),
    length,
    nullCount: 0,
    children: childData,
  });
  const recordBatch = new arrow.RecordBatch(schema, structData);
  const tableArrow = new arrow.Table([recordBatch]);
  return arrow.tableToIPC(tableArrow, asStream ? 'stream' : 'file');
}

export function fromFeather(bytes: Uint8Array): DG.DataFrame | null {
  if (!bytes) return null;
  const table = arrow.tableFromIPC(bytes);
  const columns = [];
  for (let i = 0; i < table.numCols; i++) {
    const vector = table.getChildAt(i)!;
    let values: any;
    let type = vector.type;
    const name = table.schema.fields[i].name;
    const hasNulls = vector.nullCount > 0;
    if (arrow.DataType.isDictionary(type)) {
      type = vector.data[vector.data.length - 1].dictionary?.type;
      if (type.typeId === arrow.Type.Utf8) {
        columns.push(stringColumnFromDictionary(name, vector));
        continue;
      } else
        values = unpackDictionaryColumn(vector);
    } else {
      values = vector.toArray();
    }

    switch (type.typeId) {
    case arrow.Type.Int8:
    case arrow.Type.Int16:
    case arrow.Type.Int32:
    case arrow.Type.Int:
      if (ArrayBuffer.isView(values)) {
        if (type.bitWidth < 64) {
          if (hasNulls)
            columns.push(DG.Column.fromInt32Array(name,
              int32ArrayWithNulls(extractWithNulls(vector, (v) => Number(v)))));
          else
            columns.push(DG.Column.fromInt32Array(name, values as Int32Array));
        } else
          columns.push(convertInt64Column(vector, name));
      } else
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.INT as DG.ColumnType, name, values));
      break;
    case arrow.Type.Uint32:
    case arrow.Type.Int64:
    case arrow.Type.Uint64:
      columns.push(convertInt64Column(vector, name));
      break;
    case arrow.Type.Float:
    case arrow.Type.Decimal:
      if (ArrayBuffer.isView(values)) {
        if (hasNulls) {
          // float32 -> float64 widening is exact, so one Float64 path covers both bit widths
          columns.push(DG.Column.fromFloat64Array(name,
            float64ArrayWithNulls(extractWithNulls(vector, (v) => Number(v)))));
        } else if (type.bitWidth < 64)
          columns.push(DG.Column.fromFloat32Array(name, values as Float32Array));
        else
          columns.push(DG.Column.fromFloat64Array(name, values as Float64Array));
      } else
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT as DG.ColumnType, name, values));
      break;
    case arrow.Type.Utf8:
    case arrow.Type.Interval:
      columns.push(DG.Column.fromList(DG.COLUMN_TYPE.STRING as DG.ColumnType, name, values));
      break;
    case arrow.Type.Bool:
      if (hasNulls)
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.BOOL as DG.ColumnType, name,
          Array.from({length: vector.length}, (_, j) => vector.get(j))));
      else if (ArrayBuffer.isView(values))
        columns.push(DG.Column.fromBitSet(name, DG.BitSet.fromBytes(values.buffer as ArrayBuffer, table.numRows)));
      else
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.BOOL as DG.ColumnType, name, values));
      break;
    case arrow.Type.Date:
      // get() normalizes Date32 day counts and Date64 to epoch ms and honors the null bitmap
      columns.push(DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME as DG.ColumnType, name,
        Array.from({length: vector.length}, (_, j) => {
          const v = vector.get(j);
          if (v == null) return null;
          return v instanceof Date ? v : new Date(Number(v));
        })));
      break;
    case arrow.Type.Timestamp:
      if (hasNulls)
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME as DG.ColumnType, name,
          extractWithNulls(vector, (v) => timestampBigIntToDate(BigInt(v), (type as any).unit))));
      else
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME as DG.ColumnType, name,
          values instanceof BigInt64Array ?
            Array.from(values, (b) => timestampBigIntToDate(b, type.unit)) : values));
      break;
    case arrow.Type.Time:
      if (type?.bitWidth < 64)
        columns.push(DG.Column.fromInt32Array(name, new Int32Array(values.buffer)));
      else
        columns.push(convertInt64Column(vector, name));
      break;
    default:
      columns.push(DG.Column.fromStrings(name, values));
      break;
    }
  }
  return DG.DataFrame.fromColumns(columns);
}

/** Extracts values honoring the null bitmap, which toArray() discards for primitive vectors. */
function extractWithNulls<T>(vector: arrow.Vector, convert: (raw: any) => T): (T | null)[] {
  const out = new Array(vector.length);
  let i = 0;
  for (const chunk of vector.data) {
    const nullMap = chunk.nullBitmap;
    for (let j = 0; j < chunk.length; j++) {
      // the null bitmap is indexed from the buffer start (bit granularity prevents
      // slicing), while the values buffer is already offset-adjusted by Arrow JS
      const k = chunk.offset + j;
      if (nullMap && nullMap.length && !(nullMap[k >> 3] & (1 << (k % 8))))
        out[i] = null;
      else
        out[i] = convert(chunk.values[j]);
      i++;
    }
  }
  return out;
}

function int32ArrayWithNulls(values: (number | null)[]): Int32Array {
  const result = new Int32Array(values.length);
  for (let i = 0; i < values.length; i++)
    result[i] = values[i] == null ? DG.INT_NULL : (values[i] as number);
  return result;
}

function float64ArrayWithNulls(values: (number | null)[]): Float64Array {
  const result = new Float64Array(values.length);
  for (let i = 0; i < values.length; i++)
    result[i] = values[i] == null ? DG.FLOAT_NULL : (values[i] as number);
  return result;
}

function unpackDictionaryColumn(vector: arrow.Vector) {
  const codes = new Array(vector.length);
  const ks = vector?.data[vector.data.length - 1]?.dictionary?.toArray();
  let i = 0;
  for (const chunk of vector.data) {
    const nullMap = chunk.nullBitmap || [];
    for (let j = 0; j < chunk.values.length; j++) {
      const ix = chunk.values[j];
      if (nullMap.length && !(nullMap[j >> 3] & (1 << (j % 8))))
        codes[i] = null;
      else
        codes[i] = ks[ix];
      i++;
    }
  }
  return codes;
}

function stringColumnFromDictionary(name: string, vector: arrow.Vector): DG.Column {
  const data = vector?.data[vector.data.length - 1]?.dictionary?.toArray();
  const indexes = new Int32Array(vector.length);
  let i = 0;
  for (const chunk of vector.data) {
    const nullMap = chunk.nullBitmap || [];
    for (let j = 0; j < chunk.values.length; j++) {
      if (nullMap.length && !(nullMap[j >> 3] & (1 << (j % 8))))
        indexes[i] = -1;
      else
        indexes[i] = chunk.values[j];
      i++;
    }
  }
  return DG.Column.fromIndexes(name, data, indexes);
}

/** Converts a 64-bit (or Uint32) integer vector, honoring nulls: int32 + DG.INT_NULL when values fit, bigint otherwise. */
function convertInt64Column(vector: arrow.Vector, name: string): DG.Column {
  const values = extractWithNulls(vector, (v) => BigInt(v));
  let fitsInt32 = true;
  for (const v of values) {
    if (v != null && (v > BigInt(2 ** 31 - 1) || v < BigInt(-(2 ** 31) + 1))) {
      fitsInt32 = false;
      break;
    }
  }
  if (fitsInt32) {
    const result = new Int32Array(values.length);
    for (let i = 0; i < values.length; i++)
      result[i] = values[i] == null ? DG.INT_NULL : Number(values[i]);
    return DG.Column.fromInt32Array(name, result);
  }
  if (vector.nullCount > 0)
    // DG.toDart bridges each JS bigint to a Dart BigInt, as fromBigInt64Array does
    return DG.Column.fromList(DG.COLUMN_TYPE.BIG_INT as DG.ColumnType, name,
      values.map((v) => v == null ? null : DG.toDart(v)));
  const result = new BigInt64Array(values.length);
  for (let i = 0; i < values.length; i++)
    result[i] = values[i] as bigint;
  return DG.Column.fromBigInt64Array(name, result);
}

function timestampBigIntToDate(value: bigint, unit: number): Date {
  let ms: number;
  switch (unit) {
  case 0:
    ms = Number(value) * 1000;
    break;
  case 1:
    ms = Number(value);
    break;
  case 2:
    ms = Number(value / BigInt(1000));
    break;
  case 3:
    ms = Number(value / BigInt(1000000));
    break;
  default:
    throw new Error(`Unsupported time unit: ${unit}`);
  }
  return new Date(ms);
}
