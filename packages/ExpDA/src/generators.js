// generators.js

// Datasets generation tools

// Imports for call wasm runtime-system: in the main stream and in webworkers
import { callWasm } from '../wasm/callWasm';
import { getCppInput, getResult } from '../wasm/callWasmForWebWorker';

// Inputs checker
import { checkGeneratorSVMinputs } from './utils';

const SVM_GEN_FEATURES_INDEX = 0;
const SVM_GEN_LABELS_INDEX = 1;
const SVM_LABELS_NAME = 'labels';

// Generate dataset for testing support vector machine (SVM) method
export function generateDatasetForTestingSVM(kernel, kernelParams, 
  name, samplesCount, featuresCount, min, max, violatorsPercentage) {

  checkGeneratorSVMinputs(samplesCount, featuresCount, min, max, violatorsPercentage);
  
  let kernelParamsCol = DG.Column.fromList('double', 'kernelParams', kernelParams);   
    
  let output = callWasm(EDALib, 'generateDataset', 
    [kernel, kernelParamsCol, samplesCount, featuresCount, min, max, violatorsPercentage]);
  
  let df = DG.DataFrame.fromColumns(output[SVM_GEN_FEATURES_INDEX]);
  df.name = name;  
  output[SVM_GEN_LABELS_INDEX].name = SVM_LABELS_NAME;
  df.columns.add(output[SVM_GEN_LABELS_INDEX]);

  return df;
} // generateDatasetForTestingSVM
