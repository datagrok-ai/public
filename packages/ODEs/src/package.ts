/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {initMatrOperApi} from '../wasm/matrix-operations-api';
import {ODEs, solveODEs} from './solver';
import {runSolverApp, runSolverDemoApp} from './app';

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
export async function DiffStudio() {
  await runSolverApp(); 
}

//name: Diff Studio Demo
//description: Interactive solver of ordinary differential equations (ODE)
//meta.demoPath: Compute | Diff Studio
//test: demoEquaSleekX() //wait: 100
export async function demoDiffStudio(): Promise<any>  {
  await runSolverDemoApp();
}
