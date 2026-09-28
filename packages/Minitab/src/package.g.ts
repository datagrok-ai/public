import {MinitabPackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//description: Opens Minitab Worksheet file
//input: list bytes 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: mwx
export async function importMwx(bytes: Uint8Array) : Promise<any> {
  return await MinitabPackageFunctions.importMwx(bytes);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: mwx
export async function previewMwx(file: DG.FileInfo) : Promise<any> {
  return await MinitabPackageFunctions.previewMwx(file);
}

//description: Opens Minitab Project file
//input: list bytes 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: mpx
export async function importMpx(bytes: Uint8Array) : Promise<any> {
  return await MinitabPackageFunctions.importMpx(bytes);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: mpx
export async function previewMpx(file: DG.FileInfo) : Promise<any> {
  return await MinitabPackageFunctions.previewMpx(file);
}
