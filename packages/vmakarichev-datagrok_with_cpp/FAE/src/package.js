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
//input: double t0 =  {units: minutes; caption: initial; category: Time}
//input: double t1 =  {units: minutes; caption: final; category: Time}
//input: double h = 0.01 {units: minutes; caption: step; category: Time}
//input: double FFox =  {caption: FFox; category: Initial values}
//input: double KKox =  {caption: KKox; category: Initial values}
//input: double FFred =  {caption: FFred; category: Initial values}
//input: double KKred =  {caption: KKred; category: Initial values}
//input: double Ffree = {caption: Ffree; category: Initial values}
//input: double Kfree =  {caption: Kfree; category: Initial values}
//input: double FKred =  {caption: FKred; category: Initial values}
//input: double FKox =  {caption: FKox; category: Initial values}
//input: double MEAthiol =  {caption: MEAthiol; category: Initial values}
//input: double CO2 =  {caption: CO2; category: Initial values}
//input: double yO2P =  {caption: yO2P; category: Initial values}
//input: double Cystamine =  {caption: Cystamine; category: Initial values}
//input: double VL =  {caption: VL; category: Initial values}
//output: dataframe solution {caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", yGlobalScale: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
export function solveFAE(t0, t1, h, FFox, KKox, FFred, KKred, Ffree, 
  Kfree, FKred, FKox, MEAthiol, CO2, yO2P, Cystamine, VL) {
    
  let timesCount = Math.trunc((t1 - t0) / h) + 1;
  let varsCount = 14;

  return callWasm(ODEsolver, 'solveFAE', 
    [t0, t1, h, FFox, KKox, FFred, KKred, Ffree, Kfree, FKred, 
     FKox, MEAthiol, CO2, yO2P, Cystamine, VL, 
     timesCount, varsCount]);
}
