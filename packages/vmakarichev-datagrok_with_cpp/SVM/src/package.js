/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

// Imports for call wasm runtime-system: in the main stream and in webworkers
import { callWasm } from '../wasm/callWasm';
import { getCppInput, getResult } from '../wasm/callWasmForWebWorker';

// Support vector machine (SVM) tools imports
import {showTrainReport, LINEAR, RBF, POLYNOMIAL, SIGMOID, getTrainedModel} from './svm';

//tags: init
export async function init() {
  await initSVMlib();
}

//name: Generate test data (linear kernel case)
//description: Generates dataset for testing SVN with linear kernel.
//input: string name = 'Data' {caption: name; category: Dataset}
//input: int samplesCount = 1000 {caption: samples; category: Size}
//input: int featuresCount = 2 {caption: features; category: Size}
//input: double min = -39 {caption: min; category: Range}
//input: double max = 173 {caption: max; category: Range}
//input: double violatorsPercentage = 5 {caption: violators; units: %; category: Dataset}
export function generateDatasetLinear(name, samplesCount, featuresCount, 
  min, max, violatorsPercentage) 
{
  // linear kernel ID
  let kernel = 0;

  // linear kernel parameters (any values can be used)
  let kernelParams = DG.Column.fromList('double', 'kernelParams', [0, 0]);   
    
  let output = callWasm(SVMlib, 'generateDataset', 
    [kernel, kernelParams, samplesCount, featuresCount, min, max, violatorsPercentage]);

  let df = DG.DataFrame.fromColumns(output[0]);
  df.name = name;  
  output[1].name = 'labels';
  df.columns.add(output[1]);
  grok.shell.addTableView(df);
} // generateDatasetLinear

//name: Generate test data (RBF kernel case)
//description: Generates dataset for testing SVN with RBF-kernel.
//input: string name = 'Data' {caption: name; category: Dataset}
//input: double sigma = 90  {caption: sigma; category: Hyperparameters}
//input: int samplesCount = 1000 {caption: samples; category: Size}
//input: int featuresCount = 2 {caption: features; category: Size}
//input: double min = -39 {caption: min; category: Range}
//input: double max = 173 {caption: max; category: Range}
//input: double violatorsPercentage = 5 {caption: violators; units: %; category: Dataset}
export function generateDatasetRBF(name, sigma, samplesCount, featuresCount, 
  min, max, violatorsPercentage) 
{
  // RBF kernel ID
  let kernel = 2;

  // RBF kernel parameters
  let kernelParams = DG.Column.fromList('double', 'kernelParams', [sigma, 0]);   
    
  let output = callWasm(SVMlib, 'generateDataset', 
    [kernel, kernelParams, samplesCount, featuresCount, min, max, violatorsPercentage]);

  let df = DG.DataFrame.fromColumns(output[0]);
  df.name = name;  
  output[1].name = 'labels';
  df.columns.add(output[1]);
  grok.shell.addTableView(df);
} // generateDatasetRBF

//name: trainLinearKernelSVM
//meta.mlname: linear kernel LS-SVM
//meta.mlrole: train
//input: dataframe df
//input: string predict_column
//input: double gamma = 1.0 {category: Hyperparameters}
//input: bool toShowReport = false {caption: to show report; category: Report}
//output: dynamic model
export function trainLinearKernelSVM(df, predict_column, gamma, toShowReport) {  

  let model = getTrainedModel({gamma: gamma, kernel: LINEAR}, df, predict_column);   

  if(toShowReport)
    showTrainReport(df, model);

  // TODO: add model packing
  let result = new Uint8Array(1000);

  return result;
} // trainLinearKernelSVM

//name: applyLinearKernelSVM
//meta.mlname: linear kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export function applyLinearKernelSVM(df, model) {   
  // TODO: to be implemented

  return DG.DataFrame.fromColumns(
    [DG.Column.fromFloat32Array('prediction', new Float32Array(df.rowCount))]
    );
} // applyLinearKernelSVM`

//name: trainRBFkernelSVM
//meta.mlname: RBF-kernel LS-SVM
//meta.mlrole: train
//input: dataframe df
//input: string predict_column
//input: double gamma = 1.0 {category: Hyperparameters}
//input: double sigma = 1.5 {category: Hyperparameters}
//input: bool toShowReport = false {caption: to show report; category: Report}
//output: dynamic model
export function trainRBFkernelSVM(df, predict_column, gamma, sigma, toShowReport) {  

  let model = getTrainedModel(
    {gamma: gamma, kernel: RBF, sigma: sigma}, 
    df, predict_column);   

  if(toShowReport)
    showTrainReport(df, model);

  // TODO: add model packing
  let result = new Uint8Array(1000);

  return result;
} // trainRBFkernelSVM

//name: applyRBFkernelSVM
//meta.mlname: RBF-kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export function applyRBFkernelSVM(df, model) {   
  // TODO: to be implemented

  return DG.DataFrame.fromColumns(
    [DG.Column.fromFloat32Array('prediction', new Float32Array(df.rowCount))]
    );
} // applyRBFkernelSVM

//name: trainPolynomialKernelSVM
//meta.mlname: polynomial kernel LS-SVM
//meta.mlrole: train
//input: dataframe df
//input: string predict_column
//input: double gamma = 1.0 {category: Hyperparameters}
//input: double c = 1 {category: Hyperparameters}
//input: double d = 2 {category: Hyperparameters}
//input: bool toShowReport = false {caption: to show report; category: Report}
//output: dynamic model
export function trainPolynomialKernelSVM(df, predict_column, gamma, c, d, toShowReport) {  

  let model = getTrainedModel(
    {gamma: gamma, kernel: POLYNOMIAL, cParam: c, dParam: d}, 
    df, predict_column);   

  if(toShowReport)
    showTrainReport(df, model);

  // TODO: add model packing
  let result = new Uint8Array(1000);

  return result;
} // trainPolynomialKernelSVM

//name: applyPolynomialKernelSVM
//meta.mlname: polynomial kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export function applyPolynomialKernelSVM(df, model) {   
  // TODO: to be implemented

  return DG.DataFrame.fromColumns(
    [DG.Column.fromFloat32Array('prediction', new Float32Array(df.rowCount))]
    );
} // applyPolynomialKernelSVM

//name: trainSigmoidKernelSVM
//meta.mlname: sigmoid kernel LS-SVM
//meta.mlrole: train
//input: dataframe df
//input: string predict_column
//input: double gamma = 1.0 {category: Hyperparameters}
//input: double kappa = 1 {category: Hyperparameters}
//input: double theta = 1 {category: Hyperparameters}
//input: bool toShowReport = false {caption: to show report; category: Report}
//output: dynamic model
export function trainSigmoidKernelSVM(df, predict_column, gamma, kappa, theta, toShowReport) {  

  let model = getTrainedModel(
    {gamma: gamma, kernel: SIGMOID, kappa: kappa, theta: theta}, 
    df, predict_column);   

  if(toShowReport)
    showTrainReport(df, model);

  // TODO: add model packing
  let result = new Uint8Array(1000);

  return result;
} // trainSigmoidKernelSVM

//name: applySigmoidKernelSVM
//meta.mlname: sigmoid kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export function applySigmoidKernelSVM(df, model) {   
  // TODO: to be implemented

  return DG.DataFrame.fromColumns(
    [DG.Column.fromFloat32Array('prediction', new Float32Array(df.rowCount))]
    );
} // applySigmoidKernelSVM
