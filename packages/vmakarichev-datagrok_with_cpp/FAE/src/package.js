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

//name: solveFAE
//input: double t0 = 0
//input: double t1 = 10
//input: double step = 0.01
//output: dataframe df 
export function solveFAE(t0, t1, step) {
  let timesCount = Math.trunc((t1 - t0) / step) + 1;
  let varsCount = 14;
  let df = callWasm(ODEsolver, 'solveFAE', [t0, t1, step, timesCount, varsCount]);

  let cols = df.columns; 
  cols.byIndex(0).name = 't';
  cols.byIndex(1).name = 'FFox(t)';
  cols.byIndex(2).name = 'KKox(t)'; 
  cols.byIndex(3).name = 'FFred(t)';
  cols.byIndex(4).name = 'KKred(t)';
  cols.byIndex(5).name = 'Ffree(t)';
  cols.byIndex(6).name = 'Kfree(t)'; 
  cols.byIndex(7).name = 'FKred(t)';  
  cols.byIndex(8).name = 'FKox(t)'; 
  cols.byIndex(9).name = 'MEAthiol(t)';
  cols.byIndex(10).name = 'CO2(t)';
  cols.byIndex(11).name = 'yO2P(t)';
  cols.byIndex(12).name = 'Cystamine(t)';
  cols.byIndex(13).name = 'VL(t)';
  df.name = `FAE,${t0}, ${t1},step=${step}`;

  return df;
} // solveFAE


