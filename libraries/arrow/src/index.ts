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
      type = inferNumberType(values as number[]);
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
        if (type.bitWidth < 64)
          columns.push(DG.Column.fromInt32Array(name, values as Int32Array));
        else
          columns.push(convertInt64Column(values as BigInt64Array, name, true));
      } else
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.INT as DG.ColumnType, name, values));
      break;
    case arrow.Type.Uint32:
      columns.push(convertInt64Column(values, name));
      break;
    // 64-bit ints map to `bigint` regardless of the actual values — schema-oriented, to
    // mirror the Java (grok_connect) type map (int64/uint64 → bigint) so ADBC results
    // match the stored d42 baselines. Downcasting to int32 when values happen to fit
    // would produce a `int` column that fails strict type comparison against Java's.
    case arrow.Type.Int64:
    case arrow.Type.Uint64:
      columns.push(convertInt64Column(values, name, true));
      break;
    case arrow.Type.Float:
      if (!ArrayBuffer.isView(values))
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT as DG.ColumnType, name, values));
      else if (type.bitWidth < 64)
        columns.push(DG.Column.fromFloat32Array(name, values as Float32Array));
      else
        columns.push(DG.Column.fromFloat64Array(name, values as Float64Array));
      break;
    // A 128/256-bit decimal arrives as 4/8 raw `Uint32` words per value, so it cannot share
    // the `Float` path — that would read the words themselves as the column's values.
    case arrow.Type.Decimal:
      columns.push(DG.Column.fromFloat64Array(name, decimalColumnToDoubles(vector)));
      break;
    case arrow.Type.Utf8:
    case arrow.Type.Interval:
      columns.push(DG.Column.fromList(DG.COLUMN_TYPE.STRING as DG.ColumnType, name, values));
      break;
    case arrow.Type.Bool:
      if (ArrayBuffer.isView(values))
        columns.push(DG.Column.fromBitSet(name, DG.BitSet.fromBytes(values.buffer as ArrayBuffer, table.numRows)));
      else
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.BOOL as DG.ColumnType, name, values));
      break;
    case arrow.Type.Date:
    case arrow.Type.Timestamp:
      columns.push(DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME as DG.ColumnType, name, values instanceof BigInt64Array ?
        Array.from(values, (b) => timestampBigIntToDate(b, type.unit)) : values));
      break;
    case arrow.Type.Time:
      if (type?.bitWidth < 64)
        columns.push(DG.Column.fromInt32Array(name, new Int32Array(values.buffer)));
      else
        columns.push(convertInt64Column(values, name));
      break;
    // Nested types have no DataFrame equivalent, so they collapse into a string column,
    // rendered exactly the way the Java connector's `toString()` does — see `javaText`.
    case arrow.Type.List:
    case arrow.Type.FixedSizeList:
    case arrow.Type.Struct:
    case arrow.Type.Map:
      columns.push(DG.Column.fromStrings(name, nestedColumnToStrings(vector)));
      break;
    default:
      columns.push(DG.Column.fromStrings(name, values));
      break;
    }
  }
  return DG.DataFrame.fromColumns(columns);
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

function convertInt64Column(array: BigInt64Array | BigUint64Array, name: string, forceBigInt: boolean = false): DG.Column {
  if (!forceBigInt) {
    let fitsInt32 = true;
    for (const i of array)
      if (i > BigInt(2 ** 31 - 1)) {
        fitsInt32 = false;
        break;
      }
    if (fitsInt32) {
      const result: Int32Array = new Int32Array(new ArrayBuffer(array.length * 4));
      for (let i = 0; i < array.length; i++)
        result[i] = Number(array[i]);
      return DG.Column.fromInt32Array(name, result);
    }
  }
  return DG.Column.fromBigInt64Array(name, array);
}

// Decimals are decoded from their raw words rather than through `vector.get()`, which
// apache-arrow renders incorrectly for 256-bit values. Like the Java connector, which maps
// every `decimal` to `double`, precision beyond a float64 mantissa is deliberately dropped.
function decimalColumnToDoubles(vector: arrow.Vector): Float64Array {
  const type = vector.type as any;
  const wordsPerValue = type.bitWidth / 32;
  const words = vector.toArray() as Uint32Array;
  const doubles = new Float64Array(vector.length);
  for (let i = 0; i < vector.length; i++) {
    doubles[i] = vector.isValid(i)
      ? decimalToDouble(words, i * wordsPerValue, wordsPerValue, type.scale)
      : DG.FLOAT_NULL;
  }
  return doubles;
}

function decimalToDouble(words: Uint32Array, offset: number, wordsPerValue: number, scale: number): number {
  let unscaled = BigInt(0);
  for (let word = wordsPerValue - 1; word >= 0; word--)
    unscaled = (unscaled << BigInt(32)) | BigInt(words[offset + word]);

  const bits = BigInt(wordsPerValue * 32);
  if (unscaled >> (bits - BigInt(1)))
    unscaled -= BigInt(1) << bits;

  // Building the literal and letting the JS parser round it keeps the result correctly
  // rounded, matching Java's `BigDecimal.doubleValue()`.
  return scale === 0 ? Number(unscaled) : parseFloat(scaledLiteral(unscaled, scale));
}

function scaledLiteral(unscaled: bigint, scale: number): string {
  const negative = unscaled < BigInt(0);
  const digits = (negative ? -unscaled : unscaled).toString().padStart(scale + 1, '0');
  const point = digits.length - scale;
  return `${negative ? '-' : ''}${digits.slice(0, point)}.${digits.slice(point)}`;
}

function nestedColumnToStrings(vector: arrow.Vector): string[] {
  const strings = new Array<string>(vector.length);
  for (let i = 0; i < vector.length; i++)
    strings[i] = javaText(vector.type, vector.get(i));
  return strings;
}

// The `.d42` baselines these results are compared against were captured from the Java
// connector, which stores nested values as their Java `toString()`: lists and tuples as
// `[a, b]`, maps as `{k=v}`. Rendering is driven by the Arrow type rather than by the
// JavaScript value, because a `Float64` of `10` must still print as Java's `10.0`.
function javaText(type: any, value: any): string {
  if (value === null || value === undefined)
    return 'null';

  switch (type.typeId) {
  case arrow.Type.List:
  case arrow.Type.FixedSizeList: {
    const elementType = type.children[0].type;
    const items: string[] = [];
    for (const item of value as arrow.Vector)
      items.push(javaText(elementType, item));
    return `[${items.join(', ')}]`;
  }
  case arrow.Type.Struct:
    return `[${type.children.map((field: arrow.Field) => javaText(field.type, value[field.name])).join(', ')}]`;
  case arrow.Type.Map: {
    const [keyField, valueField] = type.children[0].type.children;
    const entries: string[] = [];
    for (const [key, item] of value as Iterable<[any, any]>)
      entries.push(`${javaText(keyField.type, key)}=${javaText(valueField.type, item)}`);
    return `{${entries.join(', ')}}`;
  }
  case arrow.Type.Float:
    // Java's `Double.toString` always renders a decimal point: `10` -> `10.0`.
    return Number.isInteger(value) ? (value as number).toFixed(1) : String(value);
  case arrow.Type.Bool:
    return value ? 'true' : 'false';
  default:
    return String(value);
  }
}

function inferNumberType(values: number[]): arrow.DataType {
  return values.every((v) => Number.isInteger(v)) ? new arrow.Int32() : new arrow.Float64();
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
