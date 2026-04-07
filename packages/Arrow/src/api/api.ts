import * as DG from 'datagrok-api/dg';
import * as arrow from 'apache-arrow';
import {Compression,  readParquet, Table, writeParquet, WriterPropertiesBuilder} from "parquet-wasm";

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
      values = Array.from(raw, v => v === DG.FLOAT_NULL ? null : v);
      type = inferNumberType(values as number[]);
    }
    else if (columnType === 'int') {
      const raw = (column.getRawData() as Int32Array).subarray(0, column.length);
      values = Array.from(raw, v => v === DG.INT_NULL ? null : v);
      type = new arrow.Int32();
    }
    else if (columnType === 'datetime') {
      const raw = (column.getRawData() as Float64Array).subarray(0, column.length);
      values = Array.from(raw, v => v === DG.FLOAT_NULL ? null : new Date(v / 1000));
      type = new arrow.TimestampMillisecond();
    }
    else if (columnType === 'string') {
      values = column.toList();
      type = new arrow.Utf8();  // plain string, no dictionary
    }
    else {
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

  // Build Arrow vectors
  const vectors: { [key: string]: arrow.Vector<any> } = {};
  for (const field of fields) {
    const name = field.name;
    const type = field.type;

    if (type instanceof arrow.Float64)
      vectors[name] = arrow.vectorFromArray(arrays[name] as number[], type);
    else if (type instanceof arrow.Int32)
      vectors[name] = arrow.vectorFromArray(arrays[name] as number[], type);
    else if (type instanceof arrow.TimestampMillisecond)
      vectors[name] = arrow.vectorFromArray(arrays[name] as Date[], type);
    else if (type instanceof arrow.Utf8)
      vectors[name] = arrow.vectorFromArray(arrays[name] as string[], type);
    else if (type instanceof arrow.Int64)
      vectors[name] = arrow.vectorFromArray(arrays[name] as (bigint | null)[], type);
    else
      throw new Error(`Unsupported type for field ${name}`);
  }

  const schema = new arrow.Schema(fields);
  const tableArrow = new arrow.Table(schema, vectors);
  return arrow.tableToIPC(tableArrow, asStream ? "stream" : "file");
}

export function fromFeather(bytes: Uint8Array): DG.DataFrame | null {
  if (!bytes) return null;
  const table = arrow.tableFromIPC(bytes);
  let columns = [];
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
      }
      else
        values = unpackDictionaryColumn(vector);
    }
    else {
      values = vector.toArray();
      // values = new Array(table.numRows);
      // for (let i = 0; i < table.numRows; i++)
      //   values[i] = vector.get(i);
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
            columns.push(convertInt64Column(values as BigInt64Array, name));
        }
        else
          columns.push(DG.Column.fromList(DG.COLUMN_TYPE.INT as DG.ColumnType, name, values));
        break;
      case arrow.Type.Uint32:
      case arrow.Type.Int64:
      case arrow.Type.Uint64:
        columns.push(convertInt64Column(values, name));
        break;
      case arrow.Type.Float:
      case arrow.Type.Decimal:
        if (ArrayBuffer.isView(values)) {
          if (type.bitWidth < 64)
            columns.push(DG.Column.fromFloat32Array(name, values as Float32Array));
          else
            columns.push(DG.Column.fromFloat64Array(name, values as Float64Array));
        }
        else
          columns.push(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT as DG.ColumnType, name, values));
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
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME as DG.ColumnType, name, values instanceof BigInt64Array
            ? Array.from(values, b => timestampBigIntToDate(b, type.unit)) : values));
        break;
      case arrow.Type.Time:
        if (type?.bitWidth < 64)
          columns.push(DG.Column.fromInt32Array(name, new Int32Array(values.buffer)));
        else
          columns.push(convertInt64Column(values, name));
        break;
      default:
        columns.push(DG.Column.fromStrings(name, values));
        break;
    }
  }
  return DG.DataFrame.fromColumns(columns);
}

export function toParquet(table: DG.DataFrame, compression?: Compression): Uint8Array | null {
  const arrowUint8Array = toFeather(table);
  if (arrowUint8Array == null) return null;
  const writerProperties = new WriterPropertiesBuilder().setCompression(compression ?? Compression.SNAPPY).build();
  return writeParquet(Table.fromIPCStream(arrowUint8Array), writerProperties);
}

export function fromParquet(bytes: Uint8Array): DG.DataFrame | null {
  if (!bytes) return null;
  const arrowUint8Array = readParquet(bytes).intoIPCStream();
  return fromFeather(arrowUint8Array);
}

// Unpacks dictionary vector to array. This method is much faster than calling directly toArray() on vector
function unpackDictionaryColumn(vector: arrow.Vector) {
  const codes = new Array(vector.length);
  // Gets unique values, a.k.a categories.
  // TODO: This call .toArray() is slow
  const ks = vector?.data[vector.data.length - 1]?.dictionary?.toArray();
  let i = 0;
  // Iterate over chunks of indexes and insert data to codes array
  for (const chunk of vector.data) {
    const nullMap = chunk.nullBitmap || []; // nullBitMap can be null if no null values
    for (let j = 0; j < chunk.values.length; j++) {
      const ix = chunk.values[j];
      // Fancy bit operations because the null masks pack 8 observations into each bit.
      // ix >> 3 advances the byte every 8 bits;
      // (1 << (ix % 8) checks if the bit is set for the particular position inside the byte.

      // Check nullMap.length because it can be null if there are no null values
      if (nullMap.length && !(nullMap[j >> 3] & (1 << (j % 8))))
        codes[i] = null;
      else
        codes[i] = ks[ix];
      i++;
    }
  }
  return codes
}

// Unpacks dictionary vector to DG.Column of string type. Uses direct creation of column from indexes
// and categories without copying a data
function stringColumnFromDictionary(name: string, vector: arrow.Vector): DG.Column {
  // Gets unique values, a.k.a categories
  const data = vector?.data[vector.data.length - 1]?.dictionary?.toArray();
  // Creates array for indexes
  const indexes = new Int32Array(vector.length);
  let i = 0;
  // Iterate over chunks of indexes and insert into array index. If category value is null inserts -1.
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

function convertInt64Column(array: BigInt64Array | BigUint64Array, name: string): DG.Column {
  for (const i of array) {
    if (i > 2**31 - 1)
      return DG.Column.fromBigInt64Array(name, array);
  }
  const result: Int32Array = new Int32Array(new ArrayBuffer(array.length * 4));
  for (let i = 0; i < array.length; i++)
    result[i] = Number(array[i]);
  return DG.Column.fromInt32Array(name, result);
}

function inferNumberType(values: number[]): arrow.DataType {
  return values.every(v => Number.isInteger(v)) ? new arrow.Int32() : new arrow.Float64();
}

function timestampBigIntToDate(value: bigint, unit: number): Date {
  let ms: number;
  switch (unit) {
    case 0: // SECOND
      ms = Number(value) * 1000;
      break;
    case 1: // MILLISECOND
      ms = Number(value);
      break;
    case 2: // MICROSECOND
      ms = Number(value / 1000n);
      break;
    case 3: // NANOSECOND
      ms = Number(value / 1000000n);
      break;
    default:
      throw new Error(`Unsupported time unit: ${unit}`);
  }
  return new Date(ms);
}
