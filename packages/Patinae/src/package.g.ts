import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//input: file file 
//output: view result
//meta.role: FileViewer
//meta.fileViewer: pse,prs,pml
export function previewPymol(file: DG.FileInfo) : any {
  return PackageFunctions.previewPymol(file);
}
