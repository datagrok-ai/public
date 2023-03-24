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

//tags: init
export async function init() {
  await initSVMlib();
}

//name: evaluateFaults
//input: dataframe df
//input: column col1
//input: column col2
export function evaluateFaults(df, col1, col2) {
  // TODO: add correctness check
  let arr1 = col1.getRawData();
  let arr2 = col2.getRawData();
  let size = arr1.length;
  let faults = 0;
  let res = new Float32Array(size); 
  
  for(let i = 0; i < size; i++) {
    res[i] = arr1[i] * arr2[i];

    if(res[i] < 0)
      faults++;
  }

  df.columns.add(DG.Column.fromFloat32Array('FAULTS', res));

  alert(`Error: ${100.0 * faults / size} %.`);
} // evaluateFaults

//name: generateDataset
//input: string name = 'Data'
//input: int featuresCount = 2
//input: int samplesCount = 100
//input: double min = -250
//input: double max = 300
//input: double violatorsPercentage = 10
export function generateDataset(name, featuresCount, samplesCount, min, max, violatorsPercentage) 
{
  // TODO: add correctness check

  let output = callWasm(SVMlib, 'generateDataset', 
    [featuresCount, samplesCount, min, max, violatorsPercentage]);

  let df = DG.DataFrame.fromColumns(output[0]);
  df.name = name;  
  output[1].name = 'labels';
  df.columns.add(output[1]);
  grok.shell.addTableView(df);
}

//name: generateDatasetInWebWorker
//input: string name = 'Data'
//input: int featuresCount = 2
//input: int samplesCount = 100
//input: double min = -250
//input: double max = 300
//input: double violatorsPercentage = 10
export function generateDatasetInWebWorker(name, featuresCount, samplesCount, 
  min, max, violatorsPercentage) 
{
  // TODO: add correctness check

  var worker = new Worker(new URL('../wasm/generateDatasetWorker.js', import.meta.url));
  worker.postMessage(getCppInput(SVMlib['generateDataset'].arguments,
    [featuresCount, samplesCount, min, max, violatorsPercentage]));
  worker.onmessage = function(e) {
    let output = getResult(SVMlib['generateDataset'], e.data);

    // Provide output usage!
    let df = DG.DataFrame.fromColumns(output[0]);
    df.name = name;
    output[1].name = 'labels';
    df.columns.add(output[1]);
    grok.shell.addTableView(df);
  }
}

//name: demoLinearKernel
//input: double gamma = 1 {category: Hyperparameters}
//input: dataframe dfTrain {caption: table; category: Training data}
//input: column_list datasetTrain {caption: features; category: Training data}
//input: column labelsTrain {caption: class labels; category: Training data}
//input: dataframe dfTest {caption: table; category: Test data}
//input: column_list datasetTest {caption: features; category: Test data}
export function demoLinearKernel(gamma, dfTrain, datasetTrain, labelsTrain, dfTest, datasetTest) { 

  // TODO: add correctness check

  let trainCols = datasetTrain.toList();
  let weightsCount = trainCols.length + 1;
  let paramsCount = trainCols[0].length + 1;

  let output = callWasm(SVMlib, 'demoLinearKernel', 
    [gamma, datasetTrain, labelsTrain, datasetTest, weightsCount, paramsCount]);

  output[0].name = 'prediction';
  output[1].name = 'prediction';

  dfTrain.columns.add(output[0]);
  dfTest.columns.add(output[1]);

  output[2].name = 'mu';
  output[3].name = 'sigma';

  let dfStat = DG.DataFrame.fromColumns([output[2], output[3]]);
  dfStat.name = 'MODEL_mu_sigma';
  grok.shell.addTableView(dfStat);

  output[4].name = 'weights';  

  let dfModelWeights = DG.DataFrame.fromColumns([output[4]]);
  dfModelWeights.name = 'MODEL_weights';
  grok.shell.addTableView(dfModelWeights);

  output[5].name = 'params';

  let dfModelParams = DG.DataFrame.fromColumns([output[5]]);
  dfModelParams.name = 'MODEL_params';
  grok.shell.addTableView(dfModelParams);  
} // demoLinearKernel


//name: demoLinearKernelInWebWorker
//input: double gamma = 1 {category: Hyperparameters}
//input: dataframe dfTrain {caption: table; category: Training data}
//input: column_list datasetTrain {caption: features; category: Training data}
//input: column labelsTrain {caption: class labels; category: Training data}
//input: dataframe dfTest {caption: table; category: Test data}
//input: column_list datasetTest {caption: features; category: Test data}
export function demoLinearKernelInWebWorker(gamma, 
  dfTrain, datasetTrain, labelsTrain, dfTest, datasetTest) 
{
  // TODO: add correctness check

  let trainCols = datasetTrain.toList();
  let weightsCount = trainCols.length + 1;
  let paramsCount = trainCols[0].length + 1;

  var worker = new Worker(new URL('../wasm/demoLinearKernelWorker.js', import.meta.url));

  worker.postMessage(getCppInput(SVMlib['demoLinearKernel'].arguments,
    [gamma, datasetTrain, labelsTrain, datasetTest, weightsCount, paramsCount]));

  worker.onmessage = function(e) {
    let output = getResult(SVMlib['demoLinearKernel'], e.data);

    // Provide output usage!
    output[0].name = 'prediction';
    output[1].name = 'prediction';

    dfTrain.columns.add(output[0]);
    dfTest.columns.add(output[1]);

    output[2].name = 'mu';
    output[3].name = 'sigma';

    let dfStat = DG.DataFrame.fromColumns([output[2], output[3]]);
    dfStat.name = 'MODEL_mu_sigma';
    grok.shell.addTableView(dfStat);

    output[4].name = 'weights';  

    let dfModelWeights = DG.DataFrame.fromColumns([output[4]]);
    dfModelWeights.name = 'MODEL_weights';
    grok.shell.addTableView(dfModelWeights);

    output[5].name = 'params';

    let dfModelParams = DG.DataFrame.fromColumns([output[5]]);
    dfModelParams.name = 'MODEL_params';
    grok.shell.addTableView(dfModelParams);
  }
} // demoLinearKernelInWebWorker
