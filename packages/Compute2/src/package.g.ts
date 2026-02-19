import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//name: renderRestPanel
//input: func func 
//output: widget result
export async function renderPanel(func: any) : Promise<any> {
  return await PackageFunctions.renderPanel(func);
}

//name: Model Hub
//tags: app
//output: view result
//meta.browsePath: Compute
export function modelCatalog() {
  return PackageFunctions.modelCatalog();
}

//input: dynamic treeNode 
//input: view browseView 
//meta.role:  
//meta.app:  
export function modelCatalogTreeBrowser(treeNode: any, browseView: DG.ViewBase) : void {
  PackageFunctions.modelCatalogTreeBrowser(treeNode, browseView);
}

//name: Custom Function View Editor
//tags: editor
//input: funccall call 
//output: view result
export async function CustomFunctionViewEditor(call: DG.FuncCall) : Promise<any> {
  return await PackageFunctions.CustomFunctionViewEditor(call);
}

//name: Rich Function View Editor
//tags: editor
//input: funccall call 
//output: view result
export async function RichFunctionViewEditor(call: DG.FuncCall) : Promise<any> {
  return await PackageFunctions.RichFunctionViewEditor(call);
}

//name: Tree Wizard Editor
//tags: editor
//input: funccall call 
//output: view result
export async function TreeWizardEditor(call: DG.FuncCall) : Promise<any> {
  return await PackageFunctions.TreeWizardEditor(call);
}

//input: string nqName 
//input: string version 
//input: object instanceConfig 
//output: object result
export async function StartWorkflow(nqName: string, version: string, instanceConfig?: any) {
  return await PackageFunctions.StartWorkflow(nqName, version, instanceConfig);
}

//input: object params 
//output: object result
export async function RunOptimizer(params: any) {
  return await PackageFunctions.RunOptimizer(params);
}

//name: ViewerTestApp
export async function ViewerTestApp() : Promise<void> {
  await PackageFunctions.ViewerTestApp();
}

//name: FormTestApp
export async function FormTestApp() : Promise<void> {
  await PackageFunctions.FormTestApp();
}

//name: HistoryTestApp
export async function HistoryTestApp() : Promise<void> {
  await PackageFunctions.HistoryTestApp();
}

//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockPipeline1(params: any) {
  return await PackageFunctions.MockPipeline1(params);
}

//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockPipeline2(params: any) {
  return await PackageFunctions.MockPipeline2(params);
}

//input: double a 
//input: double b 
//output: double result
export async function TestAdd2(a: number, b: number) : Promise<number> {
  return await PackageFunctions.TestAdd2(a, b);
}

//input: double a 
//input: double b 
//output: double result
export async function TestSub2(a: number, b: number) : Promise<number> {
  return await PackageFunctions.TestSub2(a, b);
}

//input: double a 
//input: double b 
//output: double result
export async function TestMul2(a: number, b: number) : Promise<number> {
  return await PackageFunctions.TestMul2(a, b);
}

//input: double a 
//input: double b 
//output: double result
export async function TestDiv2(a: number, b: number) : Promise<number> {
  return await PackageFunctions.TestDiv2(a, b);
}

//input: dataframe df 
//output: dataframe result
export async function TestDF1(df: DG.DataFrame) : Promise<any> {
  return await PackageFunctions.TestDF1(df);
}

//name: Custom View (Compute 2 Test)
//editor: Compute2:CustomFunctionViewEditor
export async function TestCustomView() : Promise<void> {
  await PackageFunctions.TestCustomView();
}

//description: Test for optimization: multiple scalars output
//input: double x1 = 1 { caption: param1; min: -3; max: 3 }
//input: double x2 = -1 { caption: param2; min: -3; max: 3 }
//input: dataframe y { caption: table }
//input: bool bool 
//output: int integer
//output: double float1
//output: double float2
//output: dataframe table1 { viewer: Line chart(block:60) | Grid(block:40) }
//output: dataframe table2 { viewer: Line chart(block:60) | Grid(block:40) }
//meta.features: {"fitting": true, "sens-analysis": true}
//meta.runOnOpen: true
//meta.runOnInput: true
//editor: Compute2:RichFunctionViewEditor
export function fitTestFunc(x1: number, x2: number, y: DG.DataFrame, bool: boolean) {
  return PackageFunctions.fitTestFunc(x1, x2, y, bool);
}

//description: Test for optimization: multiple scalars output
export async function testFittingOutputs() : Promise<void> {
  await PackageFunctions.testFittingOutputs();
}
