/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {default as init} from "parquet-wasm";
import {fromFeather as _fromFeather, toFeather as _toFeather, fromParquet as _fromParquet, toParquet as _toParquet} from "./api/api";

export const _package = new DG.Package();

//name: info
export function info() {``
  grok.shell.info(_package.webRoot);
}

//tags: autostart
export function parquetInit() {
}

//tags: init
export async function initPackage() {
  await init(_package.webRoot + 'dist/parquet_wasm_bg.wasm');
}

//name: toFeather
//description: Converts DG.DataFrame to arrow
//input: dataframe table
//input: bool asStream = true
//output: blob bytes
export function toFeather(table: DG.DataFrame, asStream: boolean = true): Uint8Array | null {
  return _toFeather(table, asStream);
}

//name: fromFeather
//description: Converts arrow ipc stream to DG.DataFrame
//input: blob bytes
//output: dataframe table
export function fromFeather(bytes: Uint8Array): DG.DataFrame | null {
  return _fromFeather(bytes);
}

//name: toParquet
//description: Converts DG.DataFrame to parquet
//input: dataframe table
//input: int compression {nullable: true}
//output: blob bytes
export function toParquet(table: DG.DataFrame, compression?: number): Uint8Array | null {
  return _toParquet(table, compression);
}

//name: fromParquet
//description: Converts binary data in parquet format to DG.DataFrame
//input: blob bytes
//output: dataframe table
export function fromParquet(bytes: Uint8Array): DG.DataFrame | null {
  return _fromParquet(bytes);
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
  const arrowUint8Array = toFeather(table, false);
  DG.Utils.download(table.name + '.feather', arrowUint8Array ?? new Uint8Array(0));
}
