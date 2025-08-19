import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: renderDevPanel
//tags: dev-tools
//input: object ent 
//output: widget result
export function renderDevPanel(ent: any) {
  return PackageFunctions.renderDevPanel(ent);
}

//name: _makeInspectorPanel
//output: widget result
//meta.inspectorPanel: true
//friendlyName: DevTools
export function _makeInspectorPanel() {
  return PackageFunctions._makeInspectorPanel();
}

//name: autostartTools
//description: DevTools autostart function
//tags: autostart
export function autostartTools() {
  return PackageFunctions.autostartTools();
}

//name: _IconTool
//description: IconTool
export function _IconTool() {
  return PackageFunctions._IconTool();
}

//name: Test Manager
//tags: app
//meta.browsePath: Admin
//top-menu: Tools | Dev | Test Manager
export async function testManager() {
  return PackageFunctions.testManager();
}

//name: TestDetectors
//top-menu: Tools | Dev | Test | Detectors...
export function testDetectors() {
  return PackageFunctions.testDetectors();
}

//name: TestDetectorsStandard
//top-menu: Tools | Dev | Test | Detectors Standard
export async function TestDetectorsStandard() {
  return PackageFunctions.TestDetectorsStandard();
}

//name: testFunctions
//input: map scope 
//output: dataframe result
export async function testFunctions(scope: any) {
  return PackageFunctions.testFunctions(scope);
}

//name: testFunction
//input: column<int> col1 
//input: column<string> col2 
//input: column<double> col3 
//input: list<string> out { optional: true }
//output: dataframe result
//meta.vectorFunc: true
export async function testFunction(col1: DG.Column, col2: DG.Column, col3: DG.Column, out?: string[]) {
  return PackageFunctions.testFunction(col1, col2, col3, out);
}

//name: ExceptionFunc
//input: int a 
//output: int result
export function exceptionFunc(a: number) {
  return PackageFunctions.exceptionFunc(a);
}
