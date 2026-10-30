import {JmpPackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//description: Opens JMP data table
//input: list bytes 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: jmp
export async function importJmp(bytes: Uint8Array) : Promise<any> {
  return await JmpPackageFunctions.importJmp(bytes);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: jmp
export async function previewJmp(file: DG.FileInfo) : Promise<any> {
  return await JmpPackageFunctions.previewJmp(file);
}
