import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
export async function init() {
  return PackageFunctions.init();
}

//name: Custom Function View Editor
//tags: editor
//input: funccall call 
//output: view result
export async function CustomFunctionViewEditor(call: DG.FuncCall) {
  return PackageFunctions.CustomFunctionViewEditor(call);
}

//name: Rich Function View Editor
//tags: editor
//input: funccall call 
//output: view result
export async function RichFunctionViewEditor(call: DG.FuncCall) {
  return PackageFunctions.RichFunctionViewEditor(call);
}

//name: Tree Wizard Editor
//tags: editor
//input: funccall call 
//output: view result
export async function TreeWizardEditor(call: DG.FuncCall) {
  return PackageFunctions.TreeWizardEditor(call);
}

//name: StartWorkflow
//input: string nqName 
//input: string version 
//input: object instanceConfig 
//output: object result
export async function StartWorkflow(nqName: string, version: string, instanceConfig?: any) {
  return PackageFunctions.StartWorkflow(nqName, version, instanceConfig);
}

//name: ViewerTestApp
//tags: test, vue
export async function ViewerTestApp() {
  return PackageFunctions.ViewerTestApp();
}

//name: FormTestApp
//tags: test, vue
export async function FormTestApp() {
  return PackageFunctions.FormTestApp();
}

//name: HistoryTestApp
//tags: test, vue
export async function HistoryTestApp() {
  return PackageFunctions.HistoryTestApp();
}

//name: MockPipeline1
//tags: test, compute2
//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockPipeline1(params: any) {
  return PackageFunctions.MockPipeline1(params);
}

//name: MockPipeline2
//tags: test, compute2
//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockPipeline2(params: any) {
  return PackageFunctions.MockPipeline2(params);
}

//name: TestAdd2
//input: double a 
//input: double b 
//output: double result
export async function TestAdd2(a: number, b: number) {
  return PackageFunctions.TestAdd2(a, b);
}

//name: TestSub2
//input: double a 
//input: double b 
//output: double result
export async function TestSub2(a: number, b: number) {
  return PackageFunctions.TestSub2(a, b);
}

//name: TestMul2
//input: string a { choices: ['0','1'] }
//input: double b 
//output: double result
export async function TestMul2(a: number, b: number) {
  return PackageFunctions.TestMul2(a, b);
}

//name: TestDiv2
//input: double a 
//input: double b 
//output: double result
export async function TestDiv2(a: number, b: number) {
  return PackageFunctions.TestDiv2(a, b);
}

//name: TestDF1
//input: dataframe df 
//output: dataframe result
export async function TestDF1(df: DG.DataFrame) {
  return PackageFunctions.TestDF1(df);
}

//name: Custom View (Compute 2 Test)
//tags: test, compute2
//editor: Compute2:CustomFunctionViewEditor
export async function TestCustomView() {
  return PackageFunctions.TestCustomView();
}
