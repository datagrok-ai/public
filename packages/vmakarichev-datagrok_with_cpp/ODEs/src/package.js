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

//name: example
//input: double t0 = 0
//input: double t1 = 10 
//input: double step = 0.01
//input: double VLinitial = 6.6
export function example(t0, t1, step, VLinitial) {

  let timesCount = Math.trunc((t1 - t0) / step) + 1;

  let resultCount = 5;

  let callReturn = callWasm(ODEsolver, 'basicExample', [t0, t1, timesCount, VLinitial, resultCount]);

  let timeCol = callReturn[0];
  timeCol.name = 'time'; 

  let cols = [timeCol];

  let resultCols = callReturn[1];

  resultCols[0].name = 'FFred';
  resultCols[1].name = 'KKred';
  resultCols[2].name = 'FFox';
  resultCols[3].name = 'KKox';
  resultCols[4].name = 'VL';

  for(let col of resultCols)
    cols.push(col);

  let df = DG.DataFrame.fromColumns(cols);

  df.name = 'Solution';

  let view = grok.shell.addTableView(df);

  view.lineChart();
  
  for(let col of resultCols) {
    let df = DG.DataFrame.fromColumns([timeCol, col]);
    df.name = col.name;    
    let view = grok.shell.addTableView(df);
    view.lineChart();
  }
}

