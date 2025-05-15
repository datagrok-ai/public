import {PackageFunctions} from './package';

//name: init
//tags: init
//output: dynamic result
export async function init() {
  return PackageFunctions.init();
}

//name: solve
//input: object problem 
//output: dataframe result
export function solve(problem: any) {
  return PackageFunctions.solve(problem);
}

//name: solveEquations
//input: object problem 
//input: object options 
//output: dataframe result
export function solveEquations(problem: any, options: any) {
  return PackageFunctions.solveEquations(problem, options);
}

//name: Diff Studio
//description: Solver of ordinary differential equations systems
//tags: app
//output: view result
//meta.browsePath: Compute
export async function runDiffStudio() {
  return PackageFunctions.runDiffStudio();
}

//name: Diff Studio Demo
//description: Interactive solver of ordinary differential equations (ODE)
//meta.demoPath: Compute | Diff Studio
//test: runDiffStudioDemo() //wait: 100 
export async function runDiffStudioDemo() {
  return PackageFunctions.runDiffStudioDemo();
}

//name: ivpFileHandler
//tags: file-handler
//input: string content 
//output: dynamic result
//meta.ext: ipv
export async function ivpFileHandler(content: string) {
  return PackageFunctions.ivpFileHandler(content);
}

//name: previewIvp
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: ivp
export async function previewIvp(file: any) {
  return PackageFunctions.previewIvp(file);
}

//name: runDiffStudioTreeBrowser
//input: dynamic treeNode 
//input: dynamic browsePanel 
//output: dynamic result
export async function runDiffStudioTreeBrowser(treeNode: any, browsePanel: any) {
  return PackageFunctions.runDiffStudioTreeBrowser(treeNode, browsePanel);
}

//name: Ball flight
//description: Ball flight simulation
//tags: model
//input: double dB { default: 0.01; category: Ball; caption: Diameter; units: m; min: 0.01; max: 0.3 }
//input: double roB { default: 200; category: Ball; caption: Density; description: Material density; units: kg/m^3; min: 200; max: 1200 }
//input: double v { default: 50; category: Throw parameters; caption: Velocity; min: 40; max: 60; units: m/sec }
//input: double a { default: 45; category: Throw parameters; caption: Angle; min: 20; max: 70; units: deg }
//output: double maxDist { caption: Max distance }
//output: double maxHeight { caption: Max height }
//output: dataframe df { caption: Trajectory; viewer: Line chart(block: 60, multiAxis: "false", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | Grid(block: 40) }
//editor: Compute:RichFunctionViewEditor
//sidebar: @compute
//meta.runOnOpen: true
//meta.runOnInput: true
//meta.features: {"sens-analysis": true, "fitting": true}
//meta.icon: files/icons/ball.png
export function ballFlight(dB: number, roB: number, v: number, a: number) {
  return PackageFunctions.ballFlight(dB, roB, v, a);
}

//name: serializeEquations
//description: Return serialized initial value problem for ordinary differential equations
//input: string problem 
//output: object serialization
export function serializeEquations(problem: string) {
  return PackageFunctions.serializeEquations(problem);
}

//name: odesToCode
//description: Perform ODEs serialization to JS-code
//input: dynamic serialization 
//output: string result
export function odesToCode(serialization: any) {
  return PackageFunctions.odesToCode(serialization);
}

//name: solveODE
//description: Solve initial value problem for ordinary differential equations
//input: string problem 
//output: dataframe result
export async function solveODE(problem: string) {
  return PackageFunctions.solveODE(problem);
}

//name: pkPdNew
//tags: model
//meta.icon: files/icons/pkpd.png
export async function pkPdNew() {
  return PackageFunctions.pkPdNew();
}

//name: PK-PD Simulation Demo
//description: In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
//meta.demoPath: Compute | PK-PD Modeling
//test: demoSimPKPD() //wait: 100 
export async function demoSimPKPD() {
  return PackageFunctions.demoSimPKPD();
}

//name: Bioreactor
//description: Controlled fab-arm exchange mechanism simulation
//tags: model
//meta.icon: files/icons/bioreactor.png
export async function Bioreactor() {
  return PackageFunctions.Bioreactor();
}

//name: Bioreactor Demo
//description: In-browser simulation of controlled fab-arm exchange mechanism
//meta.demoPath: Compute | Bioreactor
//test: demoBioreactor() //wait: 100 
export async function demoBioreactor() {
  return PackageFunctions.demoBioreactor();
}

//name: runModel
//description: Run model with Diff Studio UI
//input: string model 
//input: int inputsTabDockRatio 
//input: int graphsDockRatio 
export async function runModel(model: string, inputsTabDockRatio: number, graphsDockRatio: number) {
  return PackageFunctions.runModel(model, inputsTabDockRatio, graphsDockRatio);
}

