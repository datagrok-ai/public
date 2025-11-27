import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//input: funccall funccall 
export function openModelFromFuncall(funccall: DG.FuncCall) : void {
  PackageFunctions.openModelFromFuncall(funccall);
}

//name: OutliersSelectionViewer
//description: Creates an outliers selection viewer
//tags: viewer
//output: viewer result
export function OutliersSelection() {
  return PackageFunctions.OutliersSelection();
}

//tags: editor
//input: funccall call 
//output: view result
export function RichFunctionViewEditor(call: DG.FuncCall) {
  return PackageFunctions.RichFunctionViewEditor(call);
}

//tags: editor
//input: funccall call 
//output: view result
export function PipelineStepEditor(call: DG.FuncCall) {
  return PackageFunctions.PipelineStepEditor(call);
}

//input: dynamic func 
//output: object uploadedCalls
export async function CustomDataUploader(func: any) {
  return await PackageFunctions.CustomDataUploader(func);
}

//input: object params 
//output: widget uploadWidget
//output: funccall uploadFuncCall
export async function CustomUploader(params: any) {
  return await PackageFunctions.CustomUploader(params);
}

//input: object params 
export function CustomCustomizer(params: any) : void {
  PackageFunctions.CustomCustomizer(params);
}

//input: object params 
//output: object validator
export function SimTimeValidator(params: any) {
  return PackageFunctions.SimTimeValidator(params);
}

//input: object params 
//output: object validator
export function DesiredTempValidator(params: any) {
  return PackageFunctions.DesiredTempValidator(params);
}

//input: object params 
//output: object validator
export function InitialTempValidator(params: any) {
  return PackageFunctions.InitialTempValidator(params);
}

//input: object params 
//output: object validator
export function AmbTempValidator(params: any) {
  return PackageFunctions.AmbTempValidator(params);
}

//input: object params 
//output: object validator
export function HeatCapValidator(params: any) {
  return PackageFunctions.HeatCapValidator(params);
}

//input: object params 
//output: object result
export function CustomStringInput(params: any) {
  return PackageFunctions.CustomStringInput(params);
}

//input: object params 
//output: object result
export function ObjectCoolingSelector(params: any) {
  return PackageFunctions.ObjectCoolingSelector(params);
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
//editor: Compute:RichFunctionViewEditor
export function fitTestFunc(x1: number, x2: number, y: DG.DataFrame, bool: boolean) {
  return PackageFunctions.fitTestFunc(x1, x2, y, bool);
}

//description: Test for optimization: multiple scalars output
export async function testFittingOutputs() : Promise<void> {
  await PackageFunctions.testFittingOutputs();
}
