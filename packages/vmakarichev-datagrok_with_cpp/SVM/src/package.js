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



//name: normalizeDataset
//input: dataframe df
//input: column_list data
export function normalizeDataset(df, data) {
  let output = callWasm(SVMlib, 'normalizeDataset', [data]);

  console.log(output);
}

//name: normalizeDataset
//input: column_list data
export function normalizeDatasetInWebWorker(data) {
  var worker = new Worker(new URL('../wasm/normalizeDatasetWorker.js', import.meta.url));
  worker.postMessage(getCppInput(SVMlib['normalizeDataset'].arguments,[data]));
  worker.onmessage = function(e) {
    let output = getResult(SVMlib['normalizeDataset'], e.data);

    // Provide output usage!
  }
}

