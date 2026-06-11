import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Flow
//description: Interactive function chain designer
//tags: app
//input: string path { meta.url: true; optional: true }
//output: view result
//meta.role: app
export function funcflowApp(path?: string) : any {
  return PackageFunctions.funcflowApp(path);
}

//tags: autostart
//meta.role: autostart
export function autoS() : void {
  PackageFunctions.autoS();
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: ffjson
export function viewFuncFlow(file: DG.FileInfo) : any {
  return PackageFunctions.viewFuncFlow(file);
}

//description: Builds a flow diagram from a table creation script and opens it in the Flow editor
//input: string script 
//output: view result
export async function flowFromCreationScript(script: string) : Promise<any> {
  return await PackageFunctions.flowFromCreationScript(script);
}
