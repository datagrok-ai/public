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

//input: string script 
//input: list<string> tableIds 
//input: bool show 
//output: dynamic result
//meta.role: creationScriptEditor
//meta.includeInFlow: false
export async function openCreationScriptFlowDialog(script: string, tableIds: string[], show: boolean) : Promise<any> {
  return await PackageFunctions.openCreationScriptFlowDialog(script, tableIds, show);
}

//name: testDialog
export function testDialog() : void {
  PackageFunctions.testDialog();
}
