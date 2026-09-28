import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//input: string bytes 
//output: list tables
//meta.role: fileHandler
//meta.fileViewerCheck: Nmrium:checkNmriumJdx
//meta.ext: jdx
export async function jdxFileHandler(bytes: string) {
  return await PackageFunctions.jdxFileHandler(bytes);
}

//name: jdxFileHandler
//input: string bytes 
//output: list tables
//meta.role: fileHandler
//meta.fileViewerCheck: Nmrium:checkNmriumJdx
//meta.ext: dx
export async function dxFileHandler(bytes: string) {
  return await PackageFunctions.dxFileHandler(bytes);
}

//input: string bytes 
//output: list tables
//meta.role: fileHandler
//meta.ext: nmrium
export async function nmriumFileHandler(bytes: string) {
  return await PackageFunctions.nmriumFileHandler(bytes);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: nmrium
export function previewNMRData(file: DG.FileInfo) : any {
  return PackageFunctions.previewNMRData(file);
}

//input: file file 
//output: view result
//meta.role: fileViewer
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
