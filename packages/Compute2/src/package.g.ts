import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//name: Custom Function View Editor
//input: funccall call 
//output: view result
//meta.role: editor
export async function CustomFunctionViewEditor(call: DG.FuncCall) : Promise<any> {
  return await PackageFunctions.CustomFunctionViewEditor(call);
}

//name: Rich Function View Editor
//input: funccall call 
//output: view result
//meta.role: editor
export async function RichFunctionViewEditor(call: DG.FuncCall) : Promise<any> {
  return await PackageFunctions.RichFunctionViewEditor(call);
}

//name: Tree Wizard Editor
//input: funccall call 
//output: view result
//meta.role: editor
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

//tags: test, vue
export async function ViewerTestApp() : Promise<void> {
  await PackageFunctions.ViewerTestApp();
}

//tags: test, vue
export async function FormTestApp() : Promise<void> {
  await PackageFunctions.FormTestApp();
}

//tags: test, vue
export async function HistoryTestApp() : Promise<void> {
  await PackageFunctions.HistoryTestApp();
}

//tags: test, compute2
//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockPipeline1(params: any) {
  return await PackageFunctions.MockPipeline1(params);
}

//tags: test, compute2
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
//tags: test, compute2
//editor: Compute2:CustomFunctionViewEditor
export async function TestCustomView() : Promise<void> {
  await PackageFunctions.TestCustomView();
}
