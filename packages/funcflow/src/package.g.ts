import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: FuncFlow
//description: Interactive function chain designer
//tags: app
//input: string path { meta.url: true; optional: true }
//output: view result
//meta.role: app
export function funcflowApp(path?: string) : any {
  return PackageFunctions.funcflowApp(path);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: ffjson
export function viewFuncFlow(file: DG.FileInfo) : any {
  return PackageFunctions.viewFuncFlow(file);
}
