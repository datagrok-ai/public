import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//tags: init
export async function init() : Promise<void> {
  PackageFunctions.init();
}

//name: solve
//input: object problem 
//output: dataframe result
export function solve(problem: any) : any {
  return PackageFunctions.solve(problem);
}

//name: solveEquations
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
  return PackageFunctions.runDiffStudio();
}

//name: Diff Studio Demo
//description: Interactive solver of ordinary differential equations (ODE)
//meta.demoPath: Compute | Diff Studio
//test: runDiffStudioDemo() //wait: 100 
export async function runDiffStudioDemo() : Promise<void> {
  PackageFunctions.runDiffStudioDemo();
}

//name: ivpFileHandler
//tags: file-handler
//input: string content 
//meta.ext: ipv
export async function ivpFileHandler(content: string) : Promise<void> {
  PackageFunctions.ivpFileHandler(content);
}

//name: previewIvp
//tags: fileViewer
//input: file file 
//output: view result
//meta.fileViewer: ivp
export async function previewIvp(file: DG.FileInfo) : Promise<any> {
  return PackageFunctions.previewIvp(file);
}

//name: runDiffStudioTreeBrowser
//input: dynamic treeNode 
//meta.role: appTreeBrowser
export async function runDiffStudioTreeBrowser(treeNode: any) : Promise<void> {
  PackageFunctions.runDiffStudioTreeBrowser(treeNode);
}

//name: Ball flight
//description: Ball flight simulation
//tags: model
//input: double dB { default: 0.01; category: Ball; caption: Diameter; units: m; min: 0.01; max: 0.3 }
//input: double roB { default: 200; category: Ball; caption: Density; description: Material density; units: kg/m^3; min: 200; max: 1200 }
//input: double v { default: 50; category: Throw parameters; caption: Velocity; min: 40; max: 60; units: m/sec }
//input: double a { default: 45; category: Throw parameters; caption: Angle; min: 20; max: 70; units: deg }
//output: double maxDist {  }
//output: double maxHeight {  }
//output: dataframe df {  }
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
export function serializeEquations(problem: string) : any {
  return PackageFunctions.serializeEquations(problem);
}

//name: odesToCode
//description: Perform ODEs serialization to JS-code
//input: dynamic serialization 
//output: string result
export function odesToCode(serialization: any) : string {
  return PackageFunctions.odesToCode(serialization);
}

//name: solveODE
//description: Solve initial value problem for ordinary differential equations
//input: string problem 
//output: dataframe result
export async function solveODE(problem: string) : Promise<any> {
  return PackageFunctions.solveODE(problem);
}

//name: PK-PD
//description: In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
//tags: model
//meta.icon: files/icons/pkpd.png
export async function pkPdNew() : Promise<void> {
  PackageFunctions.pkPdNew();
}

//name: PK-PD Simulation Demo
//description: In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
//output: dynamic result
//meta.demoPath: Compute | PK-PD Modeling
//test: demoSimPKPD() //wait: 100 
export async function demoSimPKPD() : Promise<any> {
  return PackageFunctions.demoSimPKPD();
}

//name: Bioreactor
//description: Controlled fab-arm exchange mechanism simulation
//tags: model
//meta.icon: files/icons/bioreactor.png
export async function Bioreactor() : Promise<void> {
  PackageFunctions.Bioreactor();
}

//name: Bioreactor Demo
//description: In-browser simulation of controlled fab-arm exchange mechanism
//output: dynamic result
//meta.demoPath: Compute | Bioreactor
//test: demoBioreactor() //wait: 100 
export async function demoBioreactor() : Promise<any> {
  return PackageFunctions.demoBioreactor();
}

//name: runModel
//description: Run model with Diff Studio UI
//input: string model 
//input: int inputsTabDockRatio 
//input: int graphsDockRatio 
export async function runModel(model: string, inputsTabDockRatio: number, graphsDockRatio: number) : Promise<void> {
  PackageFunctions.runModel(model, inputsTabDockRatio, graphsDockRatio);
}
