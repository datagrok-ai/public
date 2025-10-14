import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: pdf
export async function previewPdf(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.previewPdf(file);
}

//tags: file-handler
//input: list bytes 
//output: list<dataframe> result
//meta.ext: pdf
export async function viewPdf(bytes: Uint8Array) : Promise<any> {
  return await PackageFunctions.viewPdf(bytes);
}

//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: docx
export async function previewDocx(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.previewDocx(file);
}

//tags: file-handler
//input: list bytes 
//output: list<dataframe> result
//meta.ext: docx
export async function viewDocx(bytes: Uint8Array) : Promise<any> {
  return await PackageFunctions.viewDocx(bytes);
}

//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: rtf
export async function previewRtf(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.previewRtf(file);
}

//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: tex
export async function previewTex(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.previewTex(file);
}
