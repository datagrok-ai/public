import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: dev-tools
//input: object ent 
//output: widget result
export function renderDevPanel(ent: any) : any {
  return PackageFunctions.renderDevPanel(ent);
}

//output: widget result
//meta.inspectorPanel: true
//friendlyName: DevTools
export function _makeInspectorPanel() : any {
  return PackageFunctions._makeInspectorPanel();
}

//description: DevTools autostart function
//tags: autostart
export function autostartTools() : void {
  PackageFunctions.autostartTools();
}

//description: IconTool
export function _IconTool() : void {
  PackageFunctions._IconTool();
}

//name: Test Manager
//tags: app
//meta.browsePath: Admin
//top-menu: Tools | Dev | Test Manager
export async function testManager() : Promise<void> {
  await PackageFunctions.testManager();
}

//name: TestDetectors
//top-menu: Tools | Dev | Test | Detectors...
export function testDetectors() : void {
  PackageFunctions.testDetectors();
}

//top-menu: Tools | Dev | Test | Detectors Standard
export async function TestDetectorsStandard() : Promise<void> {
  await PackageFunctions.TestDetectorsStandard();
}

//input: map scope 
//output: dataframe result
export async function testFunctions(scope: any) : Promise<any> {
  return await PackageFunctions.testFunctions(scope);
}

//input: column<int> col1 
//input: column<string> col2 
//input: column<double> col3 
//input: list<string> out { optional: true }
//output: dataframe result
//meta.vectorFunc: true
export async function testFunction(col1: DG.Column, col2: DG.Column, col3: DG.Column, out?: string[]) : Promise<any> {
  return await PackageFunctions.testFunction(col1, col2, col3, out);
}

//name: ExceptionFunc
//input: int a 
//output: int result
export function exceptionFunc(a: number) : number {
  return PackageFunctions.exceptionFunc(a);
}
