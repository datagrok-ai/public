/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import { tableFromIPC, tableFromArrays, tableToIPC} from 'apache-arrow';
//@ts-ignore
import { default as init, readParquet, writeParquet, WriterPropertiesBuilder, Compression } from './arrow1';
import { Buffer } from 'buffer';
export const _package = new DG.Package();



//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}


//tags: init
export async function parquetInit() {
  await init(_package.webRoot + 'dist/arrow1_bg.wasm');
}

//input: list bytes
//output: list tables
//tags: file-handler
//meta.ext: parquet
export function parquetFileHandler(bytes: any) {
  const table = tableFromIPC(readParquet(bytes)).toArray();
  const df = DG.DataFrame.fromObjects(table);
  if (df)
    return grok.shell.addTableView(df);
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
    return grok.shell.addTableView(df);
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
