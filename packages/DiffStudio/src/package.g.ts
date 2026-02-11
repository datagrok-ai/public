import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//meta.role: init
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
//output: view result
//meta.role: app
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

//input: string content 
//meta.role: fileHandler
//meta.ext: ipv
export async function ivpFileHandler(content: string) : Promise<void> {
  await PackageFunctions.ivpFileHandler(content);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: ivp
export async function previewIvp(file: DG.FileInfo) : Promise<any> {
  return await PackageFunctions.previewIvp(file);
}

//input: dynamic treeNode 
export async function runDiffStudioTreeBrowser(treeNode: any) : Promise<void> {
  await PackageFunctions.runDiffStudioTreeBrowser(treeNode);
}

//name: Acid Production
//description: Gluconic acid (GA) production by Aspergillus niger modeling
//meta.role: model
//meta.icon: files/icons/ga-production.png
export async function acidProduction() : Promise<void> {
  await PackageFunctions.acidProduction();
}

//name: Ball flight
//description: Ball flight simulation
//input: double dB = 0.01 { category: Ball; caption: Diameter; units: m; min: 0.01; max: 0.3; minFormula: roB / 20000; maxFormula: roB / 4000 }
//input: double roB = 200 { category: Ball; caption: Density; description: Material density; units: kg/m^3; min: 200; max: 1200 }
//input: double v = 50 { category: Throw parameters; caption: Velocity; min: 40; max: 60; units: m/sec }
//input: double a = 45 { category: Throw parameters; caption: Angle; min: 20; max: 70; units: deg }
//output: double maxDist { caption: Max distance }
//output: double maxHeight { caption: Max height }
//output: dataframe df { caption: Trajectory; viewer: Line chart(block: 60, multiAxis: "false", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | Grid(block: 40) }
//meta.role: model
//editor: Compute:RichFunctionViewEditor
//sidebar: @compute
//meta.runOnOpen: true
//meta.runOnInput: true
//meta.features: {"sens-analysis": true, "fitting": true}
//meta.icon: files/icons/ball.png
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

//name: PK-PD
//description: In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
//meta.role: model
//meta.icon: files/icons/pkpd.png
export async function pkPdNew() : Promise<void> {
  await PackageFunctions.pkPdNew();
}

//name: PK-PD Simulation Demo
//description: In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
//output: dynamic result
//meta.demoPath: Compute | PK-PD Modeling
//test: demoSimPKPD() //wait: 100 
export async function demoSimPKPD() : Promise<any> {
  return await PackageFunctions.demoSimPKPD();
}

//name: Pollution
//description: The chemical reaction part of the air pollution model developed at The Dutch National Institute of Public Health and Environmental Protection
//meta.role: model
//meta.icon: files/icons/pollution.png
export async function pollution() : Promise<void> {
  await PackageFunctions.pollution();
}

//description: Controlled fab-arm exchange mechanism simulation
//meta.role: model
//meta.icon: files/icons/bioreactor.png
export async function Bioreactor() : Promise<void> {
  await PackageFunctions.Bioreactor();
}

//name: Bioreactor Demo
//description: In-browser simulation of controlled fab-arm exchange mechanism
//output: dynamic result
//meta.demoPath: Compute | Bioreactor
//test: demoBioreactor() //wait: 100 
export async function demoBioreactor() : Promise<any> {
  return await PackageFunctions.demoBioreactor();
}

//description: Run model with Diff Studio UI
//input: string model 
//input: int inputsTabDockRatio 
//input: int graphsDockRatio 
export async function runModel(model: string, inputsTabDockRatio: number, graphsDockRatio: number) : Promise<void> {
  await PackageFunctions.runModel(model, inputsTabDockRatio, graphsDockRatio);
}
