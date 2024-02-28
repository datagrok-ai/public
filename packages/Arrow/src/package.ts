/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {tableFromIPC, tableFromArrays, tableToIPC, DataType, Type, Vector} from 'apache-arrow';
//@ts-ignore
import { default as init, readParquet, writeParquet, WriterPropertiesBuilder, Compression } from './arrow1';
import { Buffer } from 'buffer';
export const _package = new DG.Package();



//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}


//tags: autostart
export async function parquetInit() {
  await init(_package.webRoot + 'dist/arrow1_bg.wasm');
}

//input: list bytes
//output: list tables
//tags: file-handler
//meta.ext: parquet
export function parquetFileHandler(bytes: WithImplicitCoercion<ArrayBuffer | SharedArrayBuffer>) {
  const table = tableFromIPC(readParquet(bytes)).toArray();
  const df = DG.DataFrame.fromObjects(table);
  if (df)
    return [df];
}

//input: list bytes
//output: list tables
//tags: file-handler
//meta.ext: feather
export function featherFileHandler(bytes: WithImplicitCoercion<ArrayBuffer | SharedArrayBuffer>) {
  const arrow = Buffer.from(bytes);
  const table = tableFromIPC(arrow).toArray();
  const df = DG.DataFrame.fromObjects(table);
  if (df) 
    return [df];
}

//name: saveAsParquet
//description: Save as Parquet
//tags: fileExporter
export function saveAsParquet(){
  let table = grok.shell.t;
  let column_names = table.columns.names();
  const t: { [_: string]: any }= {};
  for(var i = 0; i < column_names.length; i++){
    if(table.col(column_names[i])?.type === 'int'){
      t[column_names[i]] = new Int32Array(table.columns.byName(column_names[i]).toList());
    }
    else{
      t[column_names[i]] = table.columns.byName(column_names[i]).toList();
    }
  }
  const res = tableFromArrays(t);
  const arrowUint8Array = tableToIPC(res, "stream");
  const writerProperties = new WriterPropertiesBuilder().setCompression(Compression.SNAPPY).build();
  const parquetUint8Array = writeParquet(arrowUint8Array, writerProperties);
  DG.Utils.download(table.name + '.parquet', parquetUint8Array);
}


//name: saveAsFeather
//description: Save as Feather
//tags: fileExporter
export function saveAsFeather(){
  let table = grok.shell.t;
  let column_names = table.columns.names();
  const t: { [_: string]: any } = {};
  for(var i = 0; i < column_names.length; i++){
    if(table.col(column_names[i])?.type === 'int'){
      t[column_names[i]] = new Int32Array(table.columns.byName(column_names[i]).toList());
    }
    else{
      t[column_names[i]] = table.columns.byName(column_names[i]).toList();
    }
  }
  const res = tableFromArrays(t);
  const arrowUint8Array = tableToIPC(res, "stream");
  DG.Utils.download(table.name + '.feather', arrowUint8Array);
}

//name: toParquet
//description: Convert DataFrame to parquet
//input: dataframe table
//output: blob bytes
export function toParquet(table: DG.DataFrame) {
  if (table === null) return null;
  let column_names = table.columns.names();
  const t: { [_: string]: any }= {};
  for(let i = 0; i < column_names.length; i++) {
    let column = table.columns.byName(column_names[i]);
    if(['int', 'float', 'qnum'].includes(column.type)) {
      t[column_names[i]] = column.getRawData();
    }
    else if (table.col(column_names[i])!.type === 'datetime') {
      const rawData = column.getRawData();
      t[column_names[i]] = Array.from(rawData, (v, _) =>
        v === DG.FLOAT_NULL ? null : new Date(v / 1000));
    }
    else if (['string', 'bigint'].includes(column.type)) {
      const indexes = column.getRawData();
      t[column_names[i]] = Array.from(indexes, (v, _) => column.get(v));
    }
    else
      t[column_names[i]] = column.toList();
  }
  const res = tableFromArrays(t);
  const arrowUint8Array = tableToIPC(res, "stream");
  const writerProperties = new WriterPropertiesBuilder().setCompression(Compression.SNAPPY).build();
  return writeParquet(arrowUint8Array, writerProperties);
}

//name: fromParquet
//description: Convert parquet to dataframe
//input: blob bytes
//output: dataframe table
export function fromParquet(bytes: Uint8Array | null) {
  if (bytes === null) return null;
  const table = tableFromIPC(readParquet(bytes))
  let columns = [];
  for (let i = 0; i < table.numCols; i++) {
    const vector = table.getChildAt(i)!;
    let values: any;
    let type = vector.type;
    if (DataType.isDictionary(type)) {
      values = unpackDictionaryColumn(vector);
      type = vector.data[vector.data.length - 1].dictionary?.type;
    }
    else
      values = vector.toArray();
    const name = table.schema.fields[i].name;
    switch (type.typeId) {
      case Type.Int8:
      case Type.Int16:
      case Type.Int32:
      case Type.Int:
        if (ArrayBuffer.isView(values) && type.bitWidth < 64)
          columns.push(DG.Column.fromInt32Array(name, values as Int32Array));
        else if (type.bitWidth === 64)
          columns.push(DG.Column.fromList(DG.COLUMN_TYPE.BIG_INT, name, values));
        else
          columns.push(DG.Column.fromList(DG.COLUMN_TYPE.INT, name, values));
        break;
      case Type.Uint32:
      case Type.Int64:
      case Type.Uint64:
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.BIG_INT, name, values));
        break;
      case Type.Float:
      case Type.Decimal:
        if (ArrayBuffer.isView(values))
          columns.push(DG.Column.fromFloat32Array(name, values as Float32Array));
        else
          columns.push(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, name, values));
        break;
      case Type.Utf8:
      case Type.Interval:
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.STRING, name, values));
        break;
      case Type.Bool:
        if (ArrayBuffer.isView(values))
          columns.push(DG.Column.fromBitSet(name, DG.BitSet.fromBytes(values.buffer, table.numRows)));
        else
          columns.push(DG.Column.fromList(DG.COLUMN_TYPE.BOOL, name, values));
        break;
      case Type.Date:
      case Type.Timestamp:
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, name, values));
        break;
      case Type.Time:
        if (type?.bitWidth < 64 && ArrayBuffer.isView(values))
          columns.push(DG.Column.fromInt32Array(name, new Int32Array(values.buffer)));
        else
          columns.push(DG.Column.fromList(DG.COLUMN_TYPE.BIG_INT, name, values));
        break;
      default:
        columns.push(DG.Column.fromStrings(name, values));
        break;
    }
  }
  return DG.DataFrame.fromColumns(columns);
}

function unpackDictionaryColumn(vector: Vector) {
  const codes = new Array(vector.length);
  const ks = vector?.data[vector.data.length - 1]?.dictionary?.toArray();
  let i = 0;
  for (const chunk of vector.data) {
    const nullMap = chunk.nullBitmap || [];
    for (let j = 0; j < chunk.values.length; j++) {
      const ix = chunk.values[j];
      // Fancy bit operations because the null masks pack 8 observations into each bit.
      // ix >> 3 advances the byte every 8 bits;
      // (1 << (ix % 8) checks if the bit is set for the particular position inside the byte.

      // You must check nullmap.length because if there are no null values in a chunk,
      // the nullmap doesn't exist.
      if (nullMap.length && !(nullMap[j >> 3] & (1 << (j % 8))))
        codes[i] = null;
      else
        codes[i] = ks[ix];
      i++;
    }
  }
  return codes
}
