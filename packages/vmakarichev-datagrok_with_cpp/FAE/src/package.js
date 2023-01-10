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
//input: double t0 =  {units: minutes; caption: Initial time}
//input: double t1 = {units: minutes; caption: Final time}
//input: double h =  {units: minutes; caption: Time step}
//input: double FFox =  {caption: FFox(initial)}
//input: double KKox =  {caption: KKox(initial)}
//input: double FFred =  {caption: FFred(initial)}
//input: double KKred =  {caption: KKred(initial)}
//input: double Ffree =  {caption: Ffree(initial)}
//input: double Kfree =  {caption: Kfree(initial)}
//input: double FKred =  {caption: FKred(initial)}
//input: double FKox =  {caption: FKox(initial)}
//input: double MEAthiol =  {caption: MEAthiol(initial)}
//input: double CO2 =  {caption: CO2(initial)}
//input: double yO2P =  {caption: yO2P(initial)
//input: double Cystamine =  {caption: Cystamine(initial)}
//input: double VL =  {caption: VL(initial)}
//output: dataframe solution {caption: Solution; viewer: Line chart(x: "t, time (minutes)", sharex: "true", multiAxis: "true", yGlobalScale: "true", multiAxisLegendPosition: "RightCenter") | Grid(block: 100) }
//editor: Compute:RichFunctionViewEditor
export function solveFAE(t0, t1, h, FFox, KKox, FFred, KKred, Ffree, 
  Kfree, FKred, FKox, MEAthiol, CO2, yO2P, Cystamine, VL) {
    
  let timesCount = Math.trunc((t1 - t0) / h) + 1;
  let varsCount = 14;

  let df = callWasm(ODEsolver, 'solveFAE', 
    [t0, t1, h, FFox, KKox, FFred, KKred, Ffree, Kfree, FKred, 
     FKox, MEAthiol, CO2, yO2P, Cystamine, VL, 
     timesCount, varsCount]);

  let cols = df.columns; 
  cols.byIndex(0).name = 't, time (minutes)';
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
  df.name = `FAE,${t0}, ${t1},step=${h}`;
   
  return df;
}  // solveFAE
