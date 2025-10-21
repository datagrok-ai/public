import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: file-handler
//input: string bytes 
//output: list tables
//meta.fileViewerCheck: Nmrium:checkNmriumJdx
//meta.ext: jdx
export async function jdxFileHandler(bytes: string) {
  return await PackageFunctions.jdxFileHandler(bytes);
}

//name: jdxFileHandler
//tags: file-handler
//input: string bytes 
//output: list tables
//meta.fileViewerCheck: Nmrium:checkNmriumJdx
//meta.ext: dx
export async function dxFileHandler(bytes: string) {
  return await PackageFunctions.dxFileHandler(bytes);
}

//tags: file-handler
//input: string bytes 
//output: list tables
//meta.ext: nmrium
export async function nmriumFileHandler(bytes: string) {
  return await PackageFunctions.nmriumFileHandler(bytes);
}

//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: nmrium
export function previewNMRData(file: DG.FileInfo) : any {
  return PackageFunctions.previewNMRData(file);
}

//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: dx, jdx
//meta.fileViewerCheck: Nmrium:checkNmriumJdx
export function previewNMRFromDX(file: DG.FileInfo) : any {
  return PackageFunctions.previewNMRFromDX(file);
}

//input: string content 
//output: bool result
export function checkNmriumJdx(content: string) : boolean {
  return PackageFunctions.checkNmriumJdx(content);
}
