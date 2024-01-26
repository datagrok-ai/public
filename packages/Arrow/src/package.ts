/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {tableFromArrays, tableFromIPC, tableToIPC} from 'apache-arrow';
//@ts-ignore
import {Compression, default as init, readParquet, writeParquet, WriterPropertiesBuilder} from './arrow1';
import {Buffer} from 'buffer';
import {FLOAT_NULL} from "datagrok-api/dg";

export const _package = new DG.Package();



//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}


//tags: init
export async function parquetInit() {
  await init(_package.webRoot + 'dist/arrow1_bg.wasm');
}

//name: fromParquet
//input: blob bytes
//output: dataframe table
export function fromParquet(bytes: Uint8Array) {
  let d1 = new Date();
  const table = tableFromIPC(readParquet(bytes));
  console.log(`tableFromIPC ${new Date().getTime() - d1.getTime()} ms`)
  d1 = new Date();
  const array = table.toArray();
  console.log(`table.toArray() ${new Date().getTime() - d1.getTime()} ms`)
  d1 = new Date();
  const df = DG.DataFrame.fromObjects(array);
  console.log(`.fromObjects(array) ${new Date().getTime() - d1.getTime()} ms`)
  return df;
}


//input: list bytes
//output: list tables
//tags: file-handler
//meta.ext: parquet
export function parquetFileHandler(bytes: WithImplicitCoercion<ArrayBuffer | SharedArrayBuffer>) {
  const table = tableFromIPC(readParquet(bytes));
  const array = table.toArray();
  const df = DG.DataFrame.fromObjects(array);
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
  const parquetUint8Array = toParquet(table);
  DG.Utils.download(table.name + '.parquet', parquetUint8Array);
}


//name: toParquet
//description: Convert DataFrame to parquet
//input: dataframe table
//output: blob bytes
export function toParquet(table: DG.DataFrame) {
  if (table == null) return null;
  try {
    let d1 = new Date();
    let column_names = table.columns.names();
    const t: { [_: string]: any }= {};
    for(let i = 0; i < column_names.length; i++) {
      let column = table.columns.byName(column_names[i]);
      if(['int', 'float', 'qnum'].includes(column.type)) {
        t[column_names[i]] = column.getRawData();
      }
      else if (table.col(column_names[i])?.type === 'datetime') {
        const rawData: Float64Array = (column.getRawData() as Float64Array);
        t[column_names[i]] = Array.from(rawData, (v, _) =>
            v === FLOAT_NULL ? null : new Date(v / 1000));
      }
      else if (table.col(column_names[i])?.type === 'string') {
        const indexes = column.getRawData();
        t[column_names[i]] = Array.from(indexes, (v, _) => column.get(v));
      }
      else {
        t[column_names[i]] = column.toList();
      }
    }
    console.log(`Converted columns ${new Date().getTime() - d1.getTime()} ms`);
    d1 = new Date();
    const res = tableFromArrays(t);
    console.log(`tableFromArrays ${new Date().getTime() - d1.getTime()} ms`);
    d1 = new Date();
    const arrowUint8Array = tableToIPC(res, "stream");
    console.log(`tableToIPC ${new Date().getTime() - d1.getTime()} ms`);
    d1 = new Date();
    const writerProperties = new WriterPropertiesBuilder().setCompression(Compression.SNAPPY).build();
    const writeParquet1 = writeParquet(arrowUint8Array, writerProperties);
    console.log(`writeParquet ${new Date().getTime() - d1.getTime()} ms`);
    return writeParquet1;
  } catch (e) {
    throw e;
  }
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
