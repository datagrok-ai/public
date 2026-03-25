import {FlowCytometryPackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//description: Opens FCS (Flow Cytometry Standard) file
//input: list bytes 
//output: list<dataframe> result
//meta.role: fileHandler
//meta.ext: fcs
export async function importFcs(bytes: Uint8Array) : Promise<any> {
  return await FlowCytometryPackageFunctions.importFcs(bytes);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: fcs
export async function previewFcs(file: DG.FileInfo) : Promise<any> {
  return await FlowCytometryPackageFunctions.previewFcs(file);
}
