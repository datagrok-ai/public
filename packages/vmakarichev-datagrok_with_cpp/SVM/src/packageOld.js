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

//name: generateDataset
//input: string name = 'Data'
//input: int featuresCount = 2
//input: int samplesCount = 100
//input: double min = -250
//input: double max = 300
//input: double violatorsPercentage = 10
export function generateDataset(name, featuresCount, samplesCount, min, max, violatorsPercentage) {
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
export function generateDatasetInWebWorker(name,featuresCount, samplesCount, 
  min, max, violatorsPercentage) 
{
  var worker = new Worker(new URL('../wasm/generateDatasetWorker.js', import.meta.url));
  worker.postMessage(getCppInput(SVMlib['generateDataset'].arguments,[featuresCount, samplesCount, min, max, violatorsPercentage]));
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

//name: demoSVM
//input: double gamma = 1
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
export function demoSVM(gamma, dfTrain, datasetTrain, labelsTrain) {
  let output = callWasm(SVMlib, 'demo', [gamma, datasetTrain, labelsTrain]);

  console.log(output);

  output.name = '_LABEL';

  dfTrain.columns.add(output);
}

//name: demoSVM
//input: double gamma
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
export function demoSVMInWebWorker(gamma, dfTrain, datasetTrain, labelsTrain) {
  var worker = new Worker(new URL('../wasm/demoWorker.js', import.meta.url));
  worker.postMessage(getCppInput(SVMlib['demo'].arguments,[gamma, datasetTrain, labelsTrain]));
  worker.onmessage = function(e) {
    let output = getResult(SVMlib['demo'], e.data);

    // Provide output usage!
  }
}

//name: demoSVMupd
//input: double gamma = 1
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
//input: dataframe dfTest
//input: column_list datasetTest
export function demoSVMupd(gamma, dfTrain, datasetTrain, labelsTrain, dfTest, datasetTest) {
  let output = callWasm(SVMlib, 'demoUPD', [gamma, datasetTrain, labelsTrain, datasetTest]);

  console.log(output);

  output[0].name = '_LABEL';
  output[1].name = '_LABEL';

  dfTrain.columns.add(output[0]);
  dfTest.columns.add(output[1]);
}

//name: demoSVMupd
//input: double gamma
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
//input: dataframe dfTest
//input: column_list datasetTest
export function demoSVMupdInWebWorker(gamma, dfTrain, datasetTrain, labelsTrain, dfTest, datasetTest) {
  var worker = new Worker(new URL('../wasm/demoUPDWorker.js', import.meta.url));
  worker.postMessage(getCppInput(SVMlib['demoUPD'].arguments,[gamma, datasetTrain, labelsTrain, datasetTest]));
  worker.onmessage = function(e) {
    let output = getResult(SVMlib['demoUPD'], e.data);

    // Provide output usage!
  }
}

//name: demoSVMupdUPD
//input: double gamma = 1
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
//input: dataframe dfTest
//input: column_list datasetTest
export function demoSVMupdUPD(gamma, dfTrain, datasetTrain, labelsTrain, dfTest, datasetTest) {
  let output = callWasm(SVMlib, 'demoUPDupdUPD', [gamma, datasetTrain, labelsTrain, datasetTest]);  

  console.log(output);

  output[0].name = '_LABEL';
  output[1].name = '_LABEL';

  dfTrain.columns.add(output[0]);
  dfTest.columns.add(output[1]);
}

//name: demoSVMupdUPD
//input: double gamma
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
//input: dataframe dfTest
//input: column_list datasetTest
export function demoSVMupdUPDInWebWorker(gamma, dfTrain, datasetTrain, labelsTrain, dfTest, datasetTest) {
  var worker = new Worker(new URL('../wasm/demoUPDupdUPDWorker.js', import.meta.url));
  worker.postMessage(getCppInput(SVMlib['demoUPDupdUPD'].arguments,[gamma, datasetTrain, labelsTrain, datasetTest]));
  worker.onmessage = function(e) {
    let output = getResult(SVMlib['demoUPDupdUPD'], e.data);

    // Provide output usage!
  }
}


//name: demoLinear
//input: double gamma = 1
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
//input: dataframe dfTest
//input: column_list datasetTest
export function demoLinear(gamma, dfTrain, datasetTrain, labelsTrain, dfTest, datasetTest) {

  let weightsCount = 3;
  let paramsCount = 1001;

  let output = callWasm(SVMlib, 'demoLinear', 
    [gamma, datasetTrain, labelsTrain, datasetTest, weightsCount, paramsCount]);

  console.log(output);

  output[0].name = '_LABEL';
  output[1].name = '_LABEL';

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

//name: demoLinear
//input: double gamma
//input: dataframe dfTrain
//input: column_list datasetTrain
//input: column labelsTrain
//input: dataframe dfTest
//input: column_list datasetTest
//input: int paramsCount
export function demoLinearInWebWorker(gamma, dfTrain, datasetTrain, labelsTrain, dfTest, datasetTest, paramsCount) {
  var worker = new Worker(new URL('../wasm/demoLinearWorker.js', import.meta.url));
  worker.postMessage(getCppInput(SVMlib['demoLinear'].arguments,[gamma, datasetTrain, labelsTrain, datasetTest, paramsCount]));
  worker.onmessage = function(e) {
    let output = getResult(SVMlib['demoLinear'], e.data);

    // Provide output usage!
  }
}

