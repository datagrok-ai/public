/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {solveDefault, solveIVP} from './solver-tools';
import {DiffStudio} from './app';
import {getIVP, IVP, getScriptLines, getScriptParams} from './scripting-tools';

import {getBallFlightSim} from './demo/ball-flight';
import {PK_PD_DEMO} from './demo/pk-pd';
import {BIOREACTOR_DEMO} from './demo/bioreactor';
import {DF_NAME} from './constants';
import {UI_TIME} from './ui-constants';

import {ODEs, SolverOptions} from '@datagrok/diff-grok';
import {Model} from './model';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() {}

//name: solve
//input: object problem
//output: dataframe df
export function solve(problem: ODEs): DG.DataFrame {
  return solveDefault(problem);
}

//name: solveEquations
//input: object problem
//input: object options
//output: dataframe df
export function solveEquations(problem: ODEs, options: Partial<SolverOptions>): DG.DataFrame {
  return solveIVP(problem, options);
}

//name: Diff Studio
//description: Solver of ordinary differential equations systems
//tags: app
//output: view v
//meta.browsePath: Compute
export async function runDiffStudio(): Promise<DG.ViewBase> {
  const path = grok.shell.startUri;
  const toSetStartingPath = (path === window.location.href);

  const proxiView = DG.View.create();

  setTimeout(async () => {
    proxiView.close();

    const solver = new DiffStudio(false);

    if (toSetStartingPath)
      solver.setStartingPath(path);

    const view = await solver.runSolverApp();

    if (view !== null)
      grok.shell.addPreview(view);
  }, UI_TIME.APP_RUN_SOLVING);

  return proxiView;
}

//name: Diff Studio Demo
//description: Interactive solver of ordinary differential equations (ODE)
//meta.demoPath: Compute | Diff Studio
//test: runDiffStudioDemo() //wait: 100
export async function runDiffStudioDemo(): Promise<void> {
  const solver = new DiffStudio();
  await solver.runSolverDemoApp();
}

//name: ivpFileHandler
//tags: file-handler
//input: string content
//output: list tables
//meta.ext: ivp
export async function ivpFileHandler(content: string) {
  const solver = new DiffStudio();
  await solver.runSolverApp(content);

  return [];
}

//name: previewIvp
//tags: fileViewer
//meta.fileViewer: ivp
//input: file file
//output: view preview
export async function previewIvp(file: DG.FileInfo): Promise<DG.View> {
  let path: string;

  if (!DiffStudio.isStartingUriProcessed) {
    DiffStudio.isStartingUriProcessed = true;
    path = grok.shell.startUri;
  } else
    path = window.location.href;

  const proxiView = DG.View.create();

  setTimeout(async () => {
    proxiView.close();
    const solver = new DiffStudio(false, true, true);
    grok.shell.addPreview(await solver.getFilePreview(file, path));
  }, UI_TIME.PREVIEW_RUN_SOLVING);

  return proxiView;
}

//input: dynamic treeNode
//input: view browsePanel
export async function runDiffStudioTreeBrowser(treeNode: DG.TreeViewGroup, browsePanel: DG.BrowsePanel) {
  new DiffStudio(false, false, false, {treeNode: treeNode, browsePanel: browsePanel});
}

//name: Ball flight
//tags: model
//description: Ball flight simulation
//input: double dB = 0.01 {category: Ball; caption: Diameter; units: m; min: 0.01; max: 0.3}
//input: double roB = 200 {category: Ball; caption: Density; units: kg/m^3; min: 200; max: 1200} [Material density]
//input: double v = 50 {category: Throw parameters; caption: Velocity; min: 40; max: 60; units: m/sec}
//input: double a = 45 {category: Throw parameters; caption: Angle; min: 20; max: 70; units: deg}
//output: double maxDist {caption: Max distance}
//output: double maxHeight {caption: Max height}
//output: dataframe df {caption: Trajectory; viewer: Line chart(block: 60, multiAxis: "false", multiAxisLegendPosition: "RightCenter", autoLayout: "false", showAggrSelectors: "false") | Grid(block: 40)}
//editor: Compute:RichFunctionViewEditor
//sidebar: @compute
//meta.runOnOpen: true
//meta.runOnInput: true
//meta.features: {"sens-analysis": true, "fitting": true}
//meta.icon: files/icons/ball.png
export function ballFlight(dB: number, roB: number, v: number, a: number) {
  const simlulation = getBallFlightSim(v, Math.PI * a / 180, dB, roB);
  return {
    df: simlulation,
    maxDist: simlulation.col('Distance').stats.max,
    maxHeight: simlulation.col('Height').stats.max,
  };
}

//name: serializeEquations
//description: Return serialized initial value problem for ordinary differential equations
//input: string problem
//output: object serialization
export function serializeEquations(problem: string): IVP {
  return getIVP(problem);
}

//name: odesToCode
//description: Perform ODEs serialization to JS-code
//input: object serialization
//output: string code
export function odesToCode(serialization: IVP): string {
  return getScriptLines(serialization).join('\n');
}

//name: solveODE
//description: Solve initial value problem for ordinary differential equations
//input: string problem
//output: dataframe solution
export async function solveODE(problem: string): Promise<DG.DataFrame> {
  const ivp = getIVP(problem);
  const code = getScriptLines(ivp).join('\n');

  const script = DG.Script.create(code);
  const params = getScriptParams(ivp);
  const call = script.prepare(params);

  await call.call();

  return call.outputs[DF_NAME];
}

//name: PK-PD
//tags: model
//meta.icon: files/icons/pkpd.png
export async function pkPdNew(): Promise<void> {
  await PK_PD_DEMO.run();
}

//name: PK-PD Simulation Demo
//description: In-browser two-compartment pharmacokinetic-pharmacodynamic (PK-PD) simulation
//meta.demoPath: Compute | PK-PD Modeling
//test: demoSimPKPD() //wait: 100
export async function demoSimPKPD(): Promise<any> {
  await PK_PD_DEMO.runDemo();
}

//name: Bioreactor
//tags: model
//description: Controlled fab-arm exchange mechanism simulation
//meta.icon: files/icons/bioreactor.png
export async function Bioreactor(): Promise<void> {
  await BIOREACTOR_DEMO.run();
}

//name: Bioreactor Demo
//description: In-browser simulation of controlled fab-arm exchange mechanism
//meta.demoPath: Compute | Bioreactor
//test: demoBioreactor() //wait: 100
export async function demoBioreactor(): Promise<any> {
  await BIOREACTOR_DEMO.runDemo();
}

//name: runModel
//description: Run model with Diff Studio UI
//input: string model
//input: int inputsTabDockRatio
//input: int graphsDockRatio
export async function runModel(model: string, inputsTabDockRatio: number, graphsDockRatio: number): Promise<void> {
  const diffStudioModel = new Model(model, {
    inputsTabDockRatio: inputsTabDockRatio,
    graphsDockRatio: graphsDockRatio,
  }, '');

  await diffStudioModel.run();
}
