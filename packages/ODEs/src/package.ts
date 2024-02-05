/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {initMatrOperApi} from '../wasm/matrix-operations-api';
import {ODEs, solveODEs} from './solver';
import {Solver} from './app';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() {
  await initMatrOperApi();
}

//name: solve
//input: object problem
//output: dataframe df
export function solve(problem: ODEs): DG.DataFrame {  
  return solveODEs(problem); 
}

//name: Diff Studio
//description: Solver of ordinary differential equations systems
//tags: app
export async function DiffStudio(): Promise<void> {  
  //await runSolverApp();
  const solver = new Solver();
  await solver.runSolverApp();
}

//name: Diff Studio Demo
//description: Interactive solver of ordinary differential equations (ODE)
//meta.demoPath: Compute | Diff Studio
//test: demoEquaSleekX() //wait: 100
export async function demoDiffStudio(): Promise<void> {
  const solver = new Solver();
  await solver.runSolverDemoApp();
}

//name: ivpFileHandler
//tags: file-handler
//input: string content
//output: list tables
//meta.ext: ivp
export async function ivpFileHandler(content: string) {
  const solver = new Solver();
  await solver.runSolverApp(content);

  return [];
}

//name: previewIvp
//tags: fileViewer, fileViewer-ivp
//input: file file
//output: view preview
export async function previewIvp(file: DG.FileInfo): Promise<DG.View> {
  console.log('Start prview!');
  const path = window.location.href;
  const solver = new Solver(false);
  return await solver.getFilePreview(file, path);
}

//name: foo
export async function foo() {
  const op1 = {
    "name": "stage",
    "type": "int",
    "defaultValue": 0,
    "caption": 'bla1'
  };
const op2 = {
    "name": "slider",
    "type": DG.TYPE.FLOAT,
    "defaultValue": 1,
    "caption": 'bla2'
  };

const input1 = ui.input.forProperty(DG.Property.fromOptions(op1));
const input2 = ui.input.forProperty(DG.Property.fromOptions(op2));

const form = ui.form([]);

const tabControl = ui.tabControl();
const modelPane = tabControl.addPane('1-st', () => form);
const runPane = tabControl.addPane('2-nd', () => ui.panel([form]));

form.append(input1.root);
form.append(input2.root);

const v = grok.shell.addTableView(grok.data.demo.demog(10));
v.dockManager.dock(tabControl.root, 'left');
}
