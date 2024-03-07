/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as arrow from 'apache-arrow';
import {Compression, Table, WriterPropertiesBuilder, default as init, readParquet, writeParquet} from "parquet-wasm/esm/arrow1";

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: autostart
export async function parquetInit() {
  await init(_package.webRoot + 'dist/arrow1_bg.wasm');
}

//name: toFeather
//description: Converts DG.DataFrame to arrow ipc stream
//input: dataframe table
//output: blob bytes
export function toFeather(table: DG.DataFrame): Uint8Array | null {
  if (table == null) return null;
  let column_names = table.columns.names();
  const t: { [_: string]: any } = {};
  for (let i = 0; i < column_names.length; i++) {
    let column = table.columns.byName(column_names[i]);
    if (['int', 'float', 'qnum'].includes(column.type))
      t[column_names[i]] = column.getRawData();
    else if (column.type === 'datetime') {
      const rawData: Float64Array = (column.getRawData() as Float64Array);
      t[column_names[i]] = Array.from(rawData, (v, _) => v === DG.FLOAT_NULL ? null : new Date(v / 1000));
    }
    else if (column.type === 'string') {
      const indexes = column.getRawData();
      t[column_names[i]] = Array.from(indexes, (v, _) => column.get(v));
    }
    else {
      const array = new Array(column.length);
      for (let i = 0; i < column.length; i++)
        array[i] = column.get(i);
      t[column_names[i]] = array;
    }
  }
  const res = arrow.tableFromArrays(t);
  return arrow.tableToIPC(res, "stream");
}

//name: fromFeather
//description: Converts arrow ipc stream to DG.DataFrame
//input: blob bytes
//output: dataframe table
export function fromFeather(bytes: Uint8Array): DG.DataFrame | null {
  if (bytes == null) return null;
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
        const col = stringColumnFromDictionary(name, vector);
        col.name = name;
        columns.push(col);
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
        if (ArrayBuffer.isView(values) && type.bitWidth < 64)
          columns.push(DG.Column.fromInt32Array(name, values as Int32Array));
        else if (type.bitWidth === 64)
          columns.push(DG.Column.fromBigInt64Array(name, values));
        else
          columns.push(DG.Column.fromList(DG.COLUMN_TYPE.INT as DG.ColumnType, name, values));
        break;
      case arrow.Type.Uint32:
      case arrow.Type.Int64:
      case arrow.Type.Uint64:
        columns.push(DG.Column.fromBigInt64Array(name, values));
        break;
      case arrow.Type.Float:
      case arrow.Type.Decimal:
        if (ArrayBuffer.isView(values))
          columns.push(DG.Column.fromFloat32Array(name, values as Float32Array));
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
        if (type?.bitWidth < 64 && ArrayBuffer.isView(values))
          columns.push(DG.Column.fromInt32Array(name, new Int32Array(values.buffer)));
        else
          columns.push(DG.Column.fromBigInt64Array(name, values));
        break;
      default:
        columns.push(DG.Column.fromStrings(name, values));
        break;
    }
  }
  return DG.DataFrame.fromColumns(columns);
}

//name: toParquet
//description: Converts DG.DataFrame to parquet
//input: dataframe table
//output: blob bytes
export function toParquet(table: DG.DataFrame, compression?: Compression): Uint8Array | null {
  const arrowUint8Array = toFeather(table);
  if (arrowUint8Array == null) return null;
  const writerProperties = new WriterPropertiesBuilder().setCompression(compression ?? Compression.SNAPPY).build();
  return writeParquet(Table.fromIPCStream(arrowUint8Array), writerProperties);
}

//name: fromParquet
//description: Converts binary data in parquet format to DG.DataFrame
//input: blob bytes
//output: dataframe table
export function fromParquet(bytes: Uint8Array): DG.DataFrame | null {
  if (bytes == null) return null;
  let date = new Date();
  const arrowUint8Array = readParquet(bytes).intoIPCStream();
  console.log(new Date().getTime() - date.getTime());
  date = new Date();
  let df = fromFeather(arrowUint8Array);
  console.log(new Date().getTime() - date.getTime());
  return df;
}

//input: list bytes
//output: list tables
//tags: file-handler
//meta.ext: parquet
export function parquetFileHandler(bytes: WithImplicitCoercion<ArrayBuffer | SharedArrayBuffer>) {
  return [fromParquet(bytes as Uint8Array)];
}

//input: list bytes
//output: list tables
//tags: file-handler
//meta.ext: feather
export function featherFileHandler(bytes: WithImplicitCoercion<ArrayBuffer | SharedArrayBuffer>) {
  const df = fromFeather(bytes as Uint8Array);
  if (df)
    return [df];
}

//name: saveAsParquet
//description: Save as Parquet
//tags: fileExporter
export function saveAsParquet() {
  let table = grok.shell.t;
  const parquetUint8Array = toParquet(table);
  DG.Utils.download(table.name + '.parquet', parquetUint8Array ?? new Uint8Array(0));
}

//name: saveAsFeather
//description: Save as Feather
//tags: fileExporter
export function saveAsFeather() {
  let table = grok.shell.t;
  const arrowUint8Array = toFeather(table);
  DG.Utils.download(table.name + '.feather', arrowUint8Array ?? new Uint8Array(0));
}

// Unpacks dictionary vector to array. This method is much faster than calling directly toArray() on vector
function unpackDictionaryColumn(vector: arrow.Vector) {
  const codes = new Array(vector.length);
  // Gets unique values, a.k.a categories
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
