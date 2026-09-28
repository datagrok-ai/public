import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: autostart
export function parquetInit() : void {
  PackageFunctions.parquetInit();
}

//meta.role: init
export async function initPackage() : Promise<void> {
  await PackageFunctions.initPackage();
}

//description: Converts DG.DataFrame to arrow
//input: dataframe table 
//input: bool asStream = true 
//output: blob result
export function toFeather(table: DG.DataFrame, asStream: boolean) : any {
  return PackageFunctions.toFeather(table, asStream);
}

//description: Converts arrow ipc stream to DG.DataFrame
//input: blob bytes 
//output: dataframe result
export function fromFeather(bytes: Uint8Array) : any {
  return PackageFunctions.fromFeather(bytes);
}

//description: Converts DG.DataFrame to parquet
//input: dataframe table 
//input: int compression { nullable: true }
//output: blob result
export function toParquet(table: DG.DataFrame, compression?: number) : any {
  return PackageFunctions.toParquet(table, compression);
}

//description: Converts binary data in parquet format to DG.DataFrame
//input: blob bytes 
//output: dataframe result
export function fromParquet(bytes: Uint8Array) : any {
  return PackageFunctions.fromParquet(bytes);
}

//input: list bytes 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: parquet
export function parquetFileHandler(bytes: any) : any {
  return PackageFunctions.parquetFileHandler(bytes);
}

//input: list bytes 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: feather
export function featherFileHandler(bytes: any) : any {
  return PackageFunctions.featherFileHandler(bytes);
}

//description: Save as Parquet
//meta.role: fileExporter
export function saveAsParquet() : void {
  PackageFunctions.saveAsParquet();
}

//description: Save as Feather
//meta.role: fileExporter
export function saveAsFeather() : void {
  PackageFunctions.saveAsFeather();
}
