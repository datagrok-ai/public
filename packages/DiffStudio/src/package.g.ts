import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//name: dock
export function dock() : void {
  PackageFunctions.dock();
}

//input: object problem 
//output: dataframe result
export function solve(problem: any) : any {
  return PackageFunctions.solve(problem);
}

//input: object problem 
//input: object options 
//output: dataframe result
export function solveEquations(problem: any, options: any) : any {
  return PackageFunctions.solveEquations(problem, options);
}

//name: Diff Studio
//description: Solver of ordinary differential equations systems
//tags: app
//output: view result
//meta.browsePath: Compute
export async function runDiffStudio() : Promise<any> {
  return await PackageFunctions.runDiffStudio();
}

//name: Diff Studio Demo
//description: Interactive solver of ordinary differential equations (ODE)
//meta.demoPath: Compute | Diff Studio
//test: runDiffStudioDemo() //wait: 100 
export async function runDiffStudioDemo() : Promise<void> {
  await PackageFunctions.runDiffStudioDemo();
}

//tags: file-handler
//input: string content 
//meta.ext: ipv
export async function ivpFileHandler(content: string) : Promise<void> {
  await PackageFunctions.ivpFileHandler(content);
}

//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: ivp
export async function previewIvp(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.previewIvp(file);
}

//input: dynamic treeNode 
export async function runDiffStudioTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.runDiffStudioTreeBrowser(treeNode);
}

//name: Ball flight
//description: Ball flight simulation
//tags: model
//input: double dB = 0.01 { category: Ball; caption: Diameter; units: m; min: 0.01; max: 0.3; minFormula: roB / 20000; maxFormula: roB / 4000 }
//input: double roB = 200 { category: Ball; caption: Density; description: Material density; units: kg/m^3; min: 200; max: 1200 }
//input: double v = 50 { category: Throw parameters; caption: Velocity; min: 40; max: 60; units: m/sec }
//input: double a = 45 { category: Throw parameters; caption: Angle; min: 20; max: 70; units: deg }
//output: double maxDist { caption: Max distance }
//output: double maxHeight { caption: Max height }
//output: dataframe df { caption: Trajectory; viewer: Line chart(multiAxis: "false", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | Grid() }
//editor: Compute2:RichFunctionViewEditor
//meta.runOnOpen: true
//meta.runOnInput: true
//meta.features: {"sens-analysis": true, "fitting": true}
//meta.icon: files/icons/ball.png
//meta.dockSpawnConfig: {"Trajectory / Grid": {"dock-spawn-dock-ratio": 0.3, "dock-spawn-dock-type": "right", "dock-spawn-dock-to": "Trajectory / Line chart"}, "Output": {"dock-spawn-dock-ratio": 0.15, "dock-spawn-dock-type": "down", "dock-spawn-dock-to": "Trajectory / Line chart"}}
export function ballFlight(dB: number, roB: number, v: number, a: number) {
  return PackageFunctions.ballFlight(dB, roB, v, a);
}

//description: Return serialized initial value problem for ordinary differential equations
//input: string problem 
//output: object serialization
export function serializeEquations(problem: string) : any {
  return PackageFunctions.serializeEquations(problem);
}

//description: Perform ODEs serialization to JS-code
//input: dynamic serialization 
//output: string result
export function odesToCode(serialization: any) : string {
  return PackageFunctions.odesToCode(serialization);
}

//description: Solve initial value problem for ordinary differential equations
//input: string problem 
//output: dataframe result
export async function solveODE(problem: string) : Promise<any> {
  return await PackageFunctions.solveODE(problem);
}

//description: Run model with Diff Studio UI
//input: string model 
//input: int inputsTabDockRatio 
//input: int graphsDockRatio 
export async function runModel(model: string, inputsTabDockRatio: number, graphsDockRatio: number) : Promise<void> {
  await PackageFunctions.runModel(model, inputsTabDockRatio, graphsDockRatio);
}

//tags: scriptHandler
//input: funccall ivpCall 
//meta.scriptHandler.language: ivp
//meta.scriptHandler.extensions: ivp
//meta.scriptHandler.commentStart: #
//meta.scriptHandler.codeEditorMode: python
//meta.scriptHandler.parserFunction: DiffStudio:ivpLanguageParser
//meta.scriptHandler.templateScript: #name: Template\n#language: ivp\n#equations:\n  dy/dt = -y + sin(t) / t\n\n#argument: t\n  initial = 0.01 {min: 0.01; max: 10}\n  final = 15 {min: 15; max: 150}\n  step = 0.01 {min: 0.001; max: 0.1}\n\n#inits:\n  y = 0 {min: 0; max: 9}\n
//meta.icon: files/icons/package.png
export async function ivpLanguageHandler(ivpCall: DG.FuncCall) : Promise<void> {
  await PackageFunctions.ivpLanguageHandler(ivpCall);
}

//input: string code 
//output: dynamic result
export function ivpLanguageParser(code: string) : any {
  return PackageFunctions.ivpLanguageParser(code);
}
