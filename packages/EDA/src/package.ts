/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

import {_initEDAAPI, _principalComponentAnalysis, _principalComponentAnalysisInWebWorker, _error, _errorInWebWorker, _partialLeastSquareRegression, _partialLeastSquareRegressionInWebWorker} from '../wasm/EDAAPI';

//tags: init
export async function init(): Promise<void> {
  await _initEDAAPI();
}

//name: principalComponentAnalysis
//input: dataframe table
//input: column_list columns
//input: int componentsCount
//output: dataframe result 
export function principalComponentAnalysis(table: DG.DataFrame, columns: DG.ColumnList, componentsCount: number): DG.DataFrame {
  return _principalComponentAnalysis(table, columns, componentsCount);
}

//name: principalComponentAnalysisInWebWorker
//input: dataframe table
//input: column_list columns
//input: int componentsCount
//output: dataframe result 
export async function principalComponentAnalysisInWebWorker(table: DG.DataFrame, columns: DG.ColumnList, componentsCount: number): Promise<DG.DataFrame> {
  let _output: any;
  let _promise = _principalComponentAnalysisInWebWorker(table, columns, componentsCount);

  await _promise.then(
    _result => {  _output = _result; },
    _error => {  throw new Error (`Error: ${_error}`); }
  );

  return _output;
}

//name: error
//input: dataframe df
//input: column col1
//input: column col2
//output: double mad 
export function error(df: DG.DataFrame, col1: DG.Column, col2: DG.Column): number {
  return _error(df, col1, col2);
}

//name: errorInWebWorker
//input: dataframe df
//input: column col1
//input: column col2
//output: double mad 
export async function errorInWebWorker(df: DG.DataFrame, col1: DG.Column, col2: DG.Column): Promise<number> {
  let _output: any;
  let _promise = _errorInWebWorker(df, col1, col2);

  await _promise.then(
    _result => {  _output = _result; },
    _error => {  throw new Error (`Error: ${_error}`); }
  );

  return _output;
}

//name: partialLeastSquareRegression
//input: dataframe table
//input: column_list features
//input: column predict
//input: int componentsCount
export function partialLeastSquareRegression(table: DG.DataFrame, features: DG.ColumnList, predict: DG.Column, componentsCount: number): any {
  return _partialLeastSquareRegression(table, features, predict, componentsCount);
}

//name: partialLeastSquareRegressionInWebWorker
//input: dataframe table
//input: column_list features
//input: column predict
//input: int componentsCount
export async function partialLeastSquareRegressionInWebWorker(table: DG.DataFrame, features: DG.ColumnList, predict: DG.Column, componentsCount: number): Promise<any> {
  let _output: any;
  let _promise = _partialLeastSquareRegressionInWebWorker(table, features, predict, componentsCount);

  await _promise.then(
    _result => {  _output = _result; },
    _error => {  throw new Error (`Error: ${_error}`); }
  );

  return _output;
}


