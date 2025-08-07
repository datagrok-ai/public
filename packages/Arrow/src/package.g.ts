import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: parquetInit
//tags: autostart
export function parquetInit() {
  return PackageFunctions.parquetInit();
}

//name: initPackage
//tags: init
export async function initPackage() {
  return PackageFunctions.initPackage();
}

//name: toFeather
//description: Converts DG.DataFrame to arrow
//input: dataframe table 
//input: bool asStream { default: true }
//output: blob result
export function toFeather(table: DG.DataFrame, asStream: boolean) {
  return PackageFunctions.toFeather(table, asStream);
}

//name: fromFeather
//description: Converts arrow ipc stream to DG.DataFrame
//input: blob bytes 
//output: dataframe result
export function fromFeather(bytes: Uint8Array) {
  return PackageFunctions.fromFeather(bytes);
}

//name: toParquet
//description: Converts DG.DataFrame to parquet
//input: dataframe table 
//input: int compression { nullable: true }
//output: blob result
export function toParquet(table: DG.DataFrame, compression?: number) {
  return PackageFunctions.toParquet(table, compression);
}

//name: fromParquet
//description: Converts binary data in parquet format to DG.DataFrame
//input: blob bytes 
//output: dataframe result
export function fromParquet(bytes: Uint8Array) {
  return PackageFunctions.fromParquet(bytes);
}

//name: parquetFileHandler
//tags: file-handler
//input: list bytes 
//output: list<dataframe> result
//meta.ext: parquet
export function parquetFileHandler(bytes: any) {
  return PackageFunctions.parquetFileHandler(bytes);
}

//name: featherFileHandler
//tags: file-handler
//input: list bytes 
//output: list<dataframe> result
//meta.ext: feather
export function featherFileHandler(bytes: any) {
  return PackageFunctions.featherFileHandler(bytes);
}

//name: saveAsParquet
//description: Save as Parquet
//tags: fileExporter
export function saveAsParquet() {
  return PackageFunctions.saveAsParquet();
}

//name: saveAsFeather
//description: Save as Feather
//tags: fileExporter
export function saveAsFeather() {
  return PackageFunctions.saveAsFeather();
}
