// Exploratory data analysis (EDA) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_principalComponentAnalysisInWebWorker,
  _partialLeastSquareRegressionInWebWorker} from '../wasm/EDAAPI';

import {checkWasmDimensionReducerInputs, checkGeneratorSVMinputs} from './utils';

// Principal components analysis (PCA)
export async function computePCA(table: DG.DataFrame, features: DG.ColumnList, components: number,
  center: boolean, scale: boolean): Promise<DG.DataFrame> 
{
  checkWasmDimensionReducerInputs(features, components);

  const centerNum = center ? 1 : 0;
  const scaleNum = scale ? 1 : 0;

  let _output: any;
  let _promise = _principalComponentAnalysisInWebWorker(table, features, components, centerNum, scaleNum);

  await _promise.then(
    _result => { _output = _result; },
    _error => {  throw new Error (`Error: ${_error}`); }
  );

  return _output;  
} 

// Partial least square regression (PLS)
export async function computePLS(table: DG.DataFrame, features: DG.ColumnList, predict: DG.Column, components: number): Promise<any> 
{
  checkWasmDimensionReducerInputs(features, components);

  let _output: any;
  let _promise = _partialLeastSquareRegressionInWebWorker(table, features, predict, components);

  await _promise.then(
    _result => { _output = _result; },    
    _error => {  throw new Error (`Error: ${_error}`); }
  );

  return _output;
} 
