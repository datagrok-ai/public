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
import {trainModel, showModel, showModelFullInfo, LINIEAR} from './svm';

//tags: init
export async function init() {
  await initSVMlib();
}

//name: getWrongPredictions
//input: dataframe df
//input: column col1
//input: column col2
export function getWrongPredictions(df, col1, col2) {  
  let arr1 = col1.getRawData();
  let arr2 = col2.getRawData();
  let size = arr1.length;
  let res = new Float32Array(size); 
  
  for(let i = 0; i < size; i++)
    res[i] = arr1[i] * arr2[i];

  df.columns.add(DG.Column.fromFloat32Array('FAULTS', res));
} // getWrongPredictions

//name: Generate test data (linear kernel case)
//description: Generates dataset for testing SVN with linear kernel.
//input: string name = 'Data' {caption: name; category: Dataset}
//input: int samplesCount = 100 {caption: samples; category: Size}
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
//input: double sigma = 60  {caption: sigma; category: Hyperparameters}
//input: int samplesCount = 100 {caption: samples; category: Size}
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

//name: Demo LS-SVM (linear kernel)
//description: Demo of training LS-SVM model with linear kernel.
//input: double gamma = 1.0 {category: Hyperparameters}
//input: dataframe df {caption: Table; category: Training data}
//input: column_list dataset {caption: feature columns; category: Training data}
//input: column labels {caption: labels; category: Training data}
//input: bool toAddPredictions = true {caption: to show predictions; category: Results}
//input: bool toShowWrongPredictions = true {caption: to show mistakes; category: Results}
//input: string modelInfoReport { choices:["short", "full"] }
export function demoLinearKernelLSSVM(gamma, df, dataset, labels, 
  toAddPredictions, toShowWrongPredictions, modelInfo) {  

  let hyperparameters = {gamma: gamma, kernel: LINIEAR};

  let model = trainModel(SVMlib, 'trainLSSVM', 'predictByLSSVM',
    hyperparameters, dataset, labels);  

  if(modelInfo === 'short')
    showModel(model);
  else
    showModelFullInfo(model);

  if(toAddPredictions) {
    df.columns.add(model.predictedLabels);    

    if(toShowWrongPredictions)
      getWrongPredictions(df, labels, model.predictedLabels);
  }
} // demoLinearKernelLSSVM


