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
import {getCppInput, getResult} from '../wasm/callWasmUtils';

//tags: init
export async function init() {
  await initLib();
}

//name: sum
//input: int a = 11
//input: int b = 89
//output: int res 
export function sum(a, b) {
  return callWasm(Lib, 'sum', [a, b]);
}

//name: sumInWebWorker
//input: int a = 1012
//input: int b = 99
export function sumInWebWorker(a, b) {
   var worker = new Worker(new URL('../wasm/workerSum.js', import.meta.url));              

   worker.postMessage(getCppInput(Lib['sum'].arguments, [a, b]));       

   worker.onmessage = function(e) { 

    let output = getResult(Lib['sum'], e.data);

    // Here must be some UI
    alert(`${a} + ${b} = ${output}`);
   }
}

//name: maxFloatCol
//input: dataframe df
//input: column col
//output: double max 
export function maxFloatCol(df, col) {
  return callWasm(Lib, 'maxFloatCol', [col]);
}

//name: maxFloatColInWebWorker
//input: dataframe df
//input: column col
//output: double max 
export function maxFloatColInWebWorker(df, col) {
  var worker = new Worker(new URL('../wasm/workerMaxFloatCol.js', import.meta.url));              

  worker.postMessage(getCppInput(Lib['maxFloatCol'].arguments, [col]));       

  worker.onmessage = function(e) { 

   let output = getResult(Lib['maxFloatCol'], e.data);

   // Here must be some UI
   alert(`Max of the column '${col.name}' of the dataframe '${df.name}' is ${output}`);
  }
}

//name: maxIntCol
//input: dataframe df
//input: column col
//output: int max 
export function maxIntCol(df, col) {
  return callWasm(Lib, 'maxIntCol', [col]);
}

//name: maxIntColInWebWorker
//input: dataframe df
//input: column col
//output: double max 
export function maxIntColInWebWorker(df, col) {
  var worker = new Worker(new URL('../wasm/workerMaxIntCol.js', import.meta.url));              

  worker.postMessage(getCppInput(Lib['maxIntCol'].arguments, [col]));       

  worker.onmessage = function(e) { 

   let output = getResult(Lib['maxIntCol'], e.data);

   // Here must be some UI
   alert(`Max of the column '${col.name}' of the dataframe '${df.name}' is ${output}`);
  }
}

//name: addFloatCols
//input: dataframe df
//input: column col1
//input: column col2
export function addFloatCols(df, col1, col2) {
  let res = callWasm(Lib, 'addFloatCols', [col1, col2]);  
  res.name = `sum of ${col1.name} and ${col2.name}`;
  df.columns.add(res);
}

//name: addIntCols
//input: dataframe df
//input: column col1
//input: column col2
export function addIntCols(df, col1, col2) {
  let res = callWasm(Lib, 'addIntCols', [col1, col2]);  
  res.name = `sum of ${col1.name} and ${col2.name}`;
  df.columns.add(res);
}

//name: doubledInts
//input: dataframe table
//input: column_list cols
//output: dataframe result 
export function doubledInts(table, cols) {
  let res = callWasm(Lib, 'doubledInts', [cols]);
  res.name = 'Doubled_ints';
  return res;
}

//name: doubledFloats
//input: dataframe table
//input: column_list cols
//output: dataframe result 
export function doubledFloats(table, cols) {
  let res = callWasm(Lib, 'doubledFloats', [cols]);
  res.name = 'Doubled_floats';
  return res;
}

