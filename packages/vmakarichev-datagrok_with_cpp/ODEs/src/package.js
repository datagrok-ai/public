/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

import { callWasm } from '../wasm/callWasm';

//tags: init
export async function init() {
  await initODEsolver();
}

//name: stiffExample
//input: double t0 = 0
//input: double t1 = 1 
//input: double step = 0.01
export function stiffExample(t0, t1, step) {
  let timesCount = Math.trunc((t1 - t0) / step) + 1;
  let varsCount = 3;

  let df = callWasm(ODEsolver, 'stiffExample', [t0, t1, timesCount, varsCount]);   

  let cols = df.columns; 


  cols.byIndex(0).name = 't';
  cols.byIndex(1).name = 'y1(t)';
  cols.byIndex(2).name = 'y2(t)';  

  df.name = 'Solution';

  let view = grok.shell.addTableView(df);

  view.lineChart({ markerType: 'dot' });    
}

//name: stiffExampleRK
//input: double t0 = 0
//input: double t1 = 1 
//input: double step = 0.01
export function stiffExampleRK(t0, t1, step) {
  let timesCount = Math.trunc((t1 - t0) / step) + 1;
  let varsCount = 3;

  let df = callWasm(ODEsolver, 'stiffExampleRK', [t0, t1, timesCount, varsCount]);   

  let cols = df.columns; 


  cols.byIndex(0).name = 't';
  cols.byIndex(1).name = 'y1(t)';
  cols.byIndex(2).name = 'y2(t)';  

  df.name = 'SolutionRK';

  let view = grok.shell.addTableView(df);

  view.lineChart({ markerType: 'dot' });  
  
}

//name: jnjStiff
//input: double t0 = 0
//input: double t1 = 100 
//input: double step = 0.01 
export function jnjStiff(t0, t1, step) {
  let timesCount = Math.trunc((t1 - t0) / step) + 1;
  let varsCount = 8;

  let df = callWasm(ODEsolver, 'jnjStiff', [t0, t1, timesCount, varsCount]);   

  let cols = df.columns; 

  cols.byIndex(0).name = 't';
  cols.byIndex(1).name = 'FFox(t)';
  cols.byIndex(2).name = 'KKox(t)'; 
  cols.byIndex(3).name = 'FFred(t)';
  cols.byIndex(4).name = 'KKred(t)';
  cols.byIndex(5).name = 'Ffree(t)';
  cols.byIndex(6).name = 'Kfree(t)'; 
  cols.byIndex(7).name = 'MEAthiol(t)';  

  df.name = 'FAE (stiff solver)';

  let view = grok.shell.addTableView(df);

  view.lineChart({ markerType: 'dot' });    
}

//name: jnjRK4
//input: double t0 = 0
//input: double t1 = 100 
//input: double step = 0.01 
export function jnjRK4(t0, t1, step) {
  let timesCount = Math.trunc((t1 - t0) / step) + 1;
  let varsCount = 8;

  let df = callWasm(ODEsolver, 'jnjRK4', [t0, t1, timesCount, varsCount]);   

  let cols = df.columns; 

  cols.byIndex(0).name = 't';
  cols.byIndex(1).name = 'FFox(t)';
  cols.byIndex(2).name = 'KKox(t)'; 
  cols.byIndex(3).name = 'FFred(t)';
  cols.byIndex(4).name = 'KKred(t)';
  cols.byIndex(5).name = 'Ffree(t)';
  cols.byIndex(6).name = 'Kfree(t)'; 
  cols.byIndex(7).name = 'MEAthiol(t)';  

  df.name = 'FAE (RK4)';

  let view = grok.shell.addTableView(df);

  view.lineChart({ markerType: 'dot' }); 
}

