import {PrismPackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//description: Opens GraphPad Prism file
//input: list bytes 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: prism
export async function importPrism(bytes: Uint8Array) : Promise<any> {
  return await PrismPackageFunctions.importPrism(bytes);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: prism
export async function previewPrism(file: DG.FileInfo) : Promise<any> {
  return await PrismPackageFunctions.previewPrism(file);
}
