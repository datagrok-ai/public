import * as DG from 'datagrok-api/dg';
import * as arrow from 'apache-arrow';
import {Compression,  readParquet, Table, writeParquet, WriterPropertiesBuilder} from "parquet-wasm";

export function toFeather(table: DG.DataFrame, asStream: boolean = true): Uint8Array | null {
  //todo: use direct creation of vectors from typed arrays
  if (!table) return null;
  let column_names = table.columns.names();
  const t: { [_: string]: any } = {};
  for (let i = 0; i < column_names.length; i++) {
    let column = table.columns.byName(column_names[i]);
    let columnType = column.type;
    if (columnType === 'double' || columnType === 'qnum' || columnType === 'int') {
      const rawData = (column.getRawData()).subarray(0, column.length);
      const nan = columnType !== 'int' ? DG.FLOAT_NULL : DG.INT_NULL;
      t[column_names[i]] = Array.from(rawData, (v, _) => v === nan ? null : v);
    }
    else if (columnType === 'datetime') {
      const rawData: Float64Array = (column.getRawData() as Float64Array).subarray(0, column.length);
      t[column_names[i]] = Array.from(rawData, (v, _) => v === DG.FLOAT_NULL ? null : new Date(v / 1000));
    }
    else if (columnType === 'string')
      t[column_names[i]] = column.toList();
    else {
      const columnLength = column.length;
      const array = new Array(columnLength);
      for (let i = 0; i < columnLength; i++)
        array[i] = column.get(i);
      t[column_names[i]] = array;
    }
  }
  const res = arrow.tableFromArrays(t);
  return arrow.tableToIPC(res, asStream ? "stream" : "file");
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
    else
      values = vector.toArray();

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
          columns.push(DG.Column.fromBitSet(name, DG.BitSet.fromBytes(values.buffer, table.numRows)));
        else
          columns.push(DG.Column.fromList(DG.COLUMN_TYPE.BOOL as DG.ColumnType, name, values));
        break;
      case arrow.Type.Date:
      case arrow.Type.Timestamp:
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME as DG.ColumnType, name, values));
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
