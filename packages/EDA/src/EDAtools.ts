// Exploratory data analysis (EDA) tools

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_principalComponentAnalysisInWebWorker,
  _partialLeastSquareRegressionInWebWorker, _oneWayAnovaInWebWorker} from '../wasm/EDAAPI';

import {checkComponenets, checkGeneratorSVMinputs, checkColumns} from './utils';

// Principal components analysis (PCA)
export async function computePCA(table: DG.DataFrame, features: DG.ColumnList, components: number,
  center: boolean, scale: boolean): Promise<DG.DataFrame> 
{
  checkComponenets(features, components);

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
  checkComponenets(features, components);

  let _output: any;
  let _promise = _partialLeastSquareRegressionInWebWorker(table, features, predict, components);

  await _promise.then(
    _result => { _output = _result; },    
    _error => {  throw new Error (`Error: ${_error}`); }
  );

  return _output;
} 

// Analysis of Variance (ANOVA)
export async function computeANOVA(table: DG.DataFrame, columns: DG.ColumnList): Promise<number> 
{
  checkColumns(columns);

  let _output: any;
  let _promise = _oneWayAnovaInWebWorker(table, columns);

  await _promise.then(
    _result => { _output = _result; },    
    _error => {  throw new Error (`Error: ${_error}`); }
  );

  return _output;
}
