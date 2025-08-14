import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: jdxFileHandler
//tags: file-handler
//input: string bytes 
//output: list tables
//meta.fileViewerCheck: Nmrium:checkNmriumJdx
//meta.ext: jdx
export async function jdxFileHandler(bytes: string) {
  return PackageFunctions.jdxFileHandler(bytes);
}

//name: jdxFileHandler
//tags: file-handler
//input: string bytes 
//output: list tables
//meta.fileViewerCheck: Nmrium:checkNmriumJdx
//meta.ext: dx
export async function dxFileHandler(bytes: string) {
  return PackageFunctions.dxFileHandler(bytes);
}

//name: nmriumFileHandler
//tags: file-handler
//input: string bytes 
//output: list tables
//meta.ext: nmrium
export async function nmriumFileHandler(bytes: string) {
  return PackageFunctions.nmriumFileHandler(bytes);
}

//name: previewNMRData
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: nmrium
export function previewNMRData(file: DG.FileInfo) {
  return PackageFunctions.previewNMRData(file);
}

//name: previewNMRFromDX
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: dx, jdx
//meta.fileViewerCheck: Nmrium:checkNmriumJdx
export function previewNMRFromDX(file: DG.FileInfo) {
  return PackageFunctions.previewNMRFromDX(file);
}

//name: checkNmriumJdx
//input: string content 
//output: bool result
export function checkNmriumJdx(content: string) {
  return PackageFunctions.checkNmriumJdx(content);
}
