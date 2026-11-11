import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
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
//output: view result
//meta.role: app
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

//name: Start Workflow
//description: Launch a compute workflow (pipeline) by its qualified name and open its editor.
//input: string nqName { description: Qualified name of the workflow function to launch }
//input: string version { description: Workflow version to run }
//input: object instanceConfig { description: Optional initial pipeline configuration }
//output: object result
export async function StartWorkflow(nqName: string, version: string, instanceConfig?: any) {
  return await PackageFunctions.StartWorkflow(nqName, version, instanceConfig);
}

//description: Run a compute workflow without UI and return the state of every step.
//input: string nqName { description: Qualified name of the workflow function to run }
//input: string version { description: Workflow version to run }
//input: object instanceConfig { description: Optional initial pipeline configuration }
//output: object result
export async function RunWorkflow(nqName: string, version: string, instanceConfig?: any) {
  return await PackageFunctions.RunWorkflow(nqName, version, instanceConfig);
}

//description: Save the state returned by RunWorkflow as a new workflow run in history.
//input: string nqName { description: Qualified name of the workflow function }
//input: string version { description: Workflow version }
//input: object state { description: The state of a RunWorkflow result }
//input: object metadata { description: Optional title, description and tags }
//output: funccall run
export async function SaveWorkflowRun(nqName: string, version: string, state: any, metadata?: any) {
  return await PackageFunctions.SaveWorkflowRun(nqName, version, state, metadata);
}

//description: The default Excel report of a run, as its export menu builds it: [Blob, ExcelJS.Workbook].
//input: funccall funcCall { description: The run to report }
//input: object validationState { description: Validation results by input; optional: true }
//input: object consistencyState { description: Consistency by input; optional: true }
//input: object viewersHook { description: The viewers hook of a workflow step, in place of the function's own; optional: true }
//input: object metaState { description: Meta values by output, for the hook; optional: true }
//output: object result
export async function ReportFuncCallExcel(funcCall: DG.FuncCall, validationState?: any, consistencyState?: any, viewersHook?: any, metaState?: any) {
  return await PackageFunctions.ReportFuncCallExcel(funcCall, validationState, consistencyState, viewersHook, metaState);
}

//description: Load a saved workflow run without UI and return the state of every step, as RunWorkflow does.
//input: string id { description: Meta FuncCall id of the saved workflow run }
//output: object result
export async function LoadWorkflowRun(id: string) {
  return await PackageFunctions.LoadWorkflowRun(id);
}

//description: Open a saved run by its FuncCall id — a workflow run in the Tree Wizard, or a single function run in its editor.
//input: string id { description: Meta FuncCall id of the saved workflow run }
export async function OpenWorkflowRun(id: string) : Promise<void> {
  await PackageFunctions.OpenWorkflowRun(id);
}

//name: Run Optimizer
//description: Run parameter optimization (fitting) for a model and return the resulting function calls.
//input: object params { description: Optimizer parameters: target function, variables, and objective }
//output: object result
export async function RunOptimizer(params: any) {
  return await PackageFunctions.RunOptimizer(params);
}

//name: Compare Runs
//description: Compare data across model runs: scalars or a single table column
export async function CompareRuns() : Promise<void> {
  await PackageFunctions.CompareRuns();
}

//name: Browse Runs
//description: Browse saved runs of any model in My runs and in spaces, and open or compare them
export async function BrowseRuns() : Promise<void> {
  await PackageFunctions.BrowseRuns();
}

//name: ViewerTestApp
export async function ViewerTestApp() : Promise<void> {
  await PackageFunctions.ViewerTestApp();
}

//name: FormTestApp
export async function FormTestApp() : Promise<void> {
  await PackageFunctions.FormTestApp();
}

//output: object hook
export async function FailingViewersHookMaker() {
  return await PackageFunctions.FailingViewersHookMaker();
}

//output: object hook
export async function SlowViewersHookMaker() {
  return await PackageFunctions.SlowViewersHookMaker();
}

//name: Mock Pipeline 1
//description: Sample static two-step workflow configuration used for testing the workflow engine.
//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockPipeline1(params: any) {
  return await PackageFunctions.MockPipeline1(params);
}

//name: Mock Pipeline 2
//description: Sample sequential workflow configuration with links, actions, and validators for testing.
//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockPipeline2(params: any) {
  return await PackageFunctions.MockPipeline2(params);
}

//name: Mock Keyboard Pipeline
//description: Workflow with fixed steps, a reorderable list and an invalid step, used for testing keyboard use.
//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockKeyboardPipeline(params: any) {
  return await PackageFunctions.MockKeyboardPipeline(params);
}

//name: Mock Tree Wizard Pipeline
//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockTreeWizardPipeline(params: any) {
  return await PackageFunctions.MockTreeWizardPipeline(params);
}

//name: Mock Tree Wizard Group
//description: Nested workflow of the Mock Tree Wizard Pipeline.
//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockTreeWizardGroup(params: any) {
  return await PackageFunctions.MockTreeWizardGroup(params);
}

//name: Mock Single Step Pipeline
//description: Single-step workflow used for testing the compact Tree Wizard view.
//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockSingleStepPipeline(params: any) {
  return await PackageFunctions.MockSingleStepPipeline(params);
}

//name: Mock Single Step Nested
//description: Single-step workflow reached through a ref, used for testing the compact Tree Wizard view.
//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function MockSingleStepNested(params: any) {
  return await PackageFunctions.MockSingleStepNested(params);
}

//tags: stress
//input: object params 
//output: object result
//editor: Compute2:TreeWizardEditor
export async function StressTestPipeline(params: any) {
  return await PackageFunctions.StressTestPipeline(params);
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
//output: object result
//meta.customExports: [{"name":"rec","function":"Compute2:TestCustomExportRecorder"}]
//editor: Compute2:CustomFunctionViewEditor
export async function TestCustomView() {
  return await PackageFunctions.TestCustomView();
}

//input: double a 
//output: double result
//meta.customExports: [{"name":"rec","function":"Compute2:TestCustomExportRecorder"},{"name":"none","function":"Compute2:TestCustomExportNoFile"},{"name":"hook","function":"Compute2:TestCustomExportHookArgs"}]
export async function TestCustomExportModel(a: number) : Promise<number> {
  return await PackageFunctions.TestCustomExportModel(a);
}

//input: funccall funcCall 
//output: file result
export async function TestCustomExportRecorder(funcCall: DG.FuncCall) : Promise<any> {
  return await PackageFunctions.TestCustomExportRecorder(funcCall);
}

//input: funccall funcCall 
//output: string result
export async function TestCustomExportNoFile(funcCall: DG.FuncCall) : Promise<string> {
  return await PackageFunctions.TestCustomExportNoFile(funcCall);
}

//input: funccall funcCall 
//input: object viewersHook { optional: true }
//input: object metaState { optional: true }
//output: file result
export async function TestCustomExportHookArgs(funcCall: DG.FuncCall, viewersHook?: any, metaState?: any) : Promise<any> {
  return await PackageFunctions.TestCustomExportHookArgs(funcCall, viewersHook, metaState);
}

//input: int rows 
//output: dataframe first
//output: dataframe middle
//output: dataframe last
export async function TestRfvTabs(rows: number) : Promise<any> {
  return await PackageFunctions.TestRfvTabs(rows);
}

//input: int x 
//output: dataframe res
export async function TestProjectExportSingleOut(x: number) : Promise<any> {
  return await PackageFunctions.TestProjectExportSingleOut(x);
}

//input: dataframe df 
//input: int x 
//output: dataframe res1
//output: dataframe res2
export async function TestProjectExportMultiOut(df: DG.DataFrame, x: number) : Promise<any> {
  return await PackageFunctions.TestProjectExportMultiOut(df, x);
}

//input: file file { optional: true; nullable: true }
//input: int x 
//output: dataframe res
export async function TestProjectExportFileIn(file: DG.FileInfo, x: number) : Promise<any> {
  return await PackageFunctions.TestProjectExportFileIn(file, x);
}
