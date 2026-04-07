import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

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
//meta.role: autostart
export function autostartTools() : void {
  PackageFunctions.autostartTools();
}

//description: IconTool
export function _IconTool() : void {
  PackageFunctions._IconTool();
}

//name: Test Manager
//meta.role: app
//meta.browsePath: Admin
//top-menu: Tools | Dev | Test Manager
export async function testManager() : Promise<void> {
  await PackageFunctions.testManager();
}

//input: dynamic treeNode 
//meta.role: appTreeBrowser
//meta.app: Test Manager
export async function testManagerAppTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.testManagerAppTreeBrowser(treeNode);
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
