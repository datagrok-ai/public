import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: eps
export function previewEps(file: DG.FileInfo) : any {
  return PackageFunctions.previewEps(file);
}
