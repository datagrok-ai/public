import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: openModelFromFuncall
//input: funccall funccall 
export function openModelFromFuncall(funccall: DG.FuncCall) {
  return PackageFunctions.openModelFromFuncall(funccall);
}

//name: OutliersSelectionViewer
//description: Creates an outliers selection viewer
//tags: viewer
//output: viewer result
export function OutliersSelection() {
  return PackageFunctions.OutliersSelection();
}

//name: RichFunctionViewEditor
//tags: editor
//input: funccall call 
//output: view result
export function RichFunctionViewEditor(call: DG.FuncCall) {
  return PackageFunctions.RichFunctionViewEditor(call);
}

//name: PipelineStepEditor
//tags: editor
//input: funccall call 
//output: view result
export function PipelineStepEditor(call: DG.FuncCall) {
  return PackageFunctions.PipelineStepEditor(call);
}

//name: renderRestPanel
//input: func func 
//output: widget result
export async function renderPanel(func: any) {
  return PackageFunctions.renderPanel(func);
}

//tags: init
export function init() {
  return PackageFunctions.init();
}

//name: Model Hub
//tags: app
//output: view result
//meta.browsePath: Compute
export function modelCatalog() {
  return PackageFunctions.modelCatalog();
}

//name: modelCatalogTreeBrowser
//input: dynamic treeNode 
//input: view browseView 
//meta.role:  
//meta.app:  
export function modelCatalogTreeBrowser(treeNode: any, browseView: DG.ViewBase) {
  return PackageFunctions.modelCatalogTreeBrowser(treeNode, browseView);
}

//name: CustomDataUploader
//input: dynamic func 
//output: object uploadedCalls
export async function CustomDataUploader(func: any) {
  return PackageFunctions.CustomDataUploader(func);
}

//name: CustomUploader
//input: object params 
//output: widget uploadWidget
//output: funccall uploadFuncCall
export async function CustomUploader(params: any) {
  return PackageFunctions.CustomUploader(params);
}

//name: CustomCustomizer
//input: object params 
export function CustomCustomizer(params: any) {
  return PackageFunctions.CustomCustomizer(params);
}

//name: SimTimeValidator
//input: object params 
//output: object validator
export function SimTimeValidator(params: any) {
  return PackageFunctions.SimTimeValidator(params);
}

//name: DesiredTempValidator
//input: object params 
//output: object validator
export function DesiredTempValidator(params: any) {
  return PackageFunctions.DesiredTempValidator(params);
}

//name: InitialTempValidator
//input: object params 
//output: object validator
export function InitialTempValidator(params: any) {
  return PackageFunctions.InitialTempValidator(params);
}

//name: AmbTempValidator
//input: object params 
//output: object validator
export function AmbTempValidator(params: any) {
  return PackageFunctions.AmbTempValidator(params);
}

//name: HeatCapValidator
//input: object params 
//output: object validator
export function HeatCapValidator(params: any) {
  return PackageFunctions.HeatCapValidator(params);
}

//name: CustomStringInput
//input: object params 
//output: object result
export function CustomStringInput(params: any) {
  return PackageFunctions.CustomStringInput(params);
}

//name: ObjectCoolingSelector
//input: object params 
//output: object result
export function ObjectCoolingSelector(params: any) {
  return PackageFunctions.ObjectCoolingSelector(params);
}

//name: fitTestFunc
//description: Test for optimization: multiple scalars output
//input: double x1 { caption: param1; min: -3; max: 3; default: 1 }
//input: double x2 { caption: param2; min: -3; max: 3; default: -1 }
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
//editor: Compute:RichFunctionViewEditor
export function fitTestFunc(x1: number, x2: number, y: DG.DataFrame, bool: boolean) {
  return PackageFunctions.fitTestFunc(x1, x2, y, bool);
}

//name: testFittingOutputs
//description: Test for optimization: multiple scalars output
export async function testFittingOutputs() {
  return PackageFunctions.testFittingOutputs();
}
