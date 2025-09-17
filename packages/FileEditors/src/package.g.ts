import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() {
  return PackageFunctions.info();
}

//name: previewPdf
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: pdf
export async function previewPdf(file: DG.FileInfo) {
  return PackageFunctions.previewPdf(file);
}

//name: viewPdf
//tags: file-handler
//input: list bytes 
//output: list<dataframe> result
//meta.ext: pdf
export async function viewPdf(bytes: Uint8Array) {
  return PackageFunctions.viewPdf(bytes);
}

//name: previewDocx
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: docx
export async function previewDocx(file: DG.FileInfo) {
  return PackageFunctions.previewDocx(file);
}

//name: viewDocx
//tags: file-handler
//input: list bytes 
//output: list<dataframe> result
//meta.ext: docx
export async function viewDocx(bytes: Uint8Array) {
  return PackageFunctions.viewDocx(bytes);
}

//name: previewRtf
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: rtf
export async function previewRtf(file: DG.FileInfo) {
  return PackageFunctions.previewRtf(file);
}

//name: previewTex
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: tex
export async function previewTex(file: DG.FileInfo) {
  return PackageFunctions.previewTex(file);
}
