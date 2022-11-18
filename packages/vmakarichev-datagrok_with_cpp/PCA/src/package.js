/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

import { callWasm } from '../wasm/callWasm';

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() {
  await initEigenPCA();  
}

//top-menu: Tools | Data Science | Principal Component Analysis by Eigen
//name: pca
//tags: ml
//input: dataframe table
//input: column_list columns
//input: int componentsCount
//output: dataframe result
export function pca(table, columns, componentsCount) {
  return callWasm(EigenPCA, 'principalComponentAnalysis', [columns, componentsCount]);
}

//name: mad
//input: dataframe table
//input: column column1
//input: column column2
//output: double result
export function mad(table, column1, column2) {
  return callWasm(EigenPCA, 'error', [column1, column2]);
}
