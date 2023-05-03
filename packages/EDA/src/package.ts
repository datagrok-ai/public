/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_initEDAAPI} from '../wasm/EDAAPI';
import {computePCA, computePLS} from './EDAtools';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init(): Promise<void> {
  await _initEDAAPI();
}

//name: PCA
//description: Principal component analysis (PCA).
//input: dataframe table
//input: column_list features
//input: int components = 3
//output: dataframe result {action:join(table)}
export async function PCA(table: DG.DataFrame, features: DG.ColumnList, components: number): Promise<DG.DataFrame> {
  return await computePCA(table, features, components);
}

//name: PLS
//description: Partial least square regression (PLS).
//input: dataframe table
//input: column_list features
//input: column predict
//input: int components = 3
export async function PLS(table: DG.DataFrame, features: DG.ColumnList, predict: DG.Column, components: number): Promise<void> {
  await computePLS(table, features, predict, components);
}
