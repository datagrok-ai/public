// svm.js

// Support vector machine (SVM) tools.
// Training & predicting is provided by wasm-computations.
// Least square support vector machine (LS-SVM) is implemented:
//   [1] Suykens, J., Vandewalle, J. "Least Squares Support Vector Machine Classifiers",
//	     Neural Processing Letters 9, 293-300 (1999). https://doi.org/10.1023/A:1018628609742 

// Imports for call wasm runtime-system: in the main stream and in webworkers
import { callWasm } from '../wasm/callWasm';
import { getCppInput, getResult } from '../wasm/callWasmForWebWorker';

// CONSTANTS

// kernel types
export const LINIEAR = 0;
export const POLYNOMIAL = 1;
export const RBF = 2;
export const SIGMOID = 3;

// names
const GENERAL = 'General';
const NORMALIZED = 'Normalized';
const CHARACTERISTICS = 'mean, deviation';
const PARAMETERS = 'Parameters';
const WEIGHTS = 'Weights';
const LABELS = 'Labels';
const PREDICTED = 'predicted';
const MEAN = 'mean';
const STD_DEV = 'std dev';
const MODEL_PARAMS_NAME = 'alpha';
const MODEL_WEIGHTS_NAME = 'weight';
const GAMMA = 'gamma';
const KERNEL = 'kernel';
const KERNEL_PARAMS = 'kernel params'; 
const KERNEL_PARAM_1 = 'kernel param 1';
const KERNEL_PARAM_2 = 'kernel param 2';
const TRAIN_ERROR = 'train error,%';
const MODEL_INFO = 'Model info';
const MODEL_INFO_FULL = 'Model full info';
const KERNEL_TYPE_TO_NAME_MAP = ['linear', 'polynomial', 'RBF', 'sigmoid'];

// misc
const INIT_VALUE = 0; // any number can be used
const LS_SVM_ADD_CONST = 1; // see [1] for more details

// Returnes labels predicted by the model specified
export function predict(module, predictFuncName, model, dataset)
{
  return callWasm(module, predictFuncName, 
    [ model.kernelType,
      model.kernelParams,
      model.normalizedTrainData.columns,
      model.trainLabels,      
      model.means, 
      model.stdDevs, 
      model.modelParams, 
      model.modelWeights, 
      dataset]);
} // predict

// Returns wrong labels percentage
export function getError(col1, col2)
{  
  // get raw data
  let arr1 = col1.getRawData();
  let arr2 = col2.getRawData();

  let size = arr1.length;

  // check sizes
  if(arr2.length != size)
    throw new Error('Coulmns must be of the same length.');

  let faultsCount = 0;

  for(let i = 0; i < size; i++)
    if(arr1[i] * arr2[i] < 0) // wrong predicted label has a wrong sign
      faultsCount++;
  
  return 100.0 * faultsCount / size;
}

// Returns trained LS-SVM model.
export function trainModel(module, trainFuncName, predictFuncName, 
    hyperparameters, dataset, labels) 
{  
  // create default kernel params array
  let kernelParamsArray = [INIT_VALUE, INIT_VALUE];
  
  // fill kernelParams
  switch(hyperparameters.kernel)
  {
    case LINIEAR: // no kernel parameters in the case of linear kernel
        break;
    default:
        throw new Error('Incorrect kernel ID.');
  };

  // create kernel params column
  let kernelParams = DG.Column.fromList('double', KERNEL_PARAMS, kernelParamsArray);
  
  // compute size of model params & precomputed weigths
  let trainCols = dataset.toList();  
  let modelParamsCount = trainCols[0].length + LS_SVM_ADD_CONST;
  let precomputedWeightsCount = trainCols.length + LS_SVM_ADD_CONST;
  
  // call webassembly training function
  let output = callWasm(module, trainFuncName, 
    [hyperparameters.gamma, hyperparameters.kernel, kernelParams, 
     modelParamsCount, precomputedWeightsCount, 
     dataset, labels]);       
  
  // rename output columns
  output[1].name = MEAN;
  output[2].name = STD_DEV;
  output[3].name = MODEL_PARAMS_NAME;
  output[4].name = MODEL_WEIGHTS_NAME;

  // complete model
  let model = {
    trainGamma: hyperparameters.gamma,
    kernelType: hyperparameters.kernel,
    kernelParams: kernelParams,
    trainLabels: labels,
    normalizedTrainData: DG.DataFrame.fromColumns(output[0]),
    means: output[1],
    stdDevs: output[2],
    modelParams: output[3],
    modelWeights: output[4],
    predictedLabels: undefined,
    trainError: undefined 
  };

  // compute predicted labels and add then to model specification
  let prediction = predict(module, predictFuncName, model, dataset);
  prediction.name = PREDICTED;
  model.predictedLabels = prediction;

  // evaluate train error;
  model.trainError = getError(model.trainLabels, model.predictedLabels);

  return model;
} // trainModel

// Returns dataframe with short info about model
export function getModelInfo(model) {
  let kernelParams = model.kernelParams.getRawData();

  return DG.DataFrame.fromColumns([
    DG.Column.fromList('double', GAMMA, [model.trainGamma]),
    DG.Column.fromStrings(KERNEL, [KERNEL_TYPE_TO_NAME_MAP[model.kernelType]]),
    DG.Column.fromList('double', KERNEL_PARAM_1, [kernelParams[0]]),
    DG.Column.fromList('double', KERNEL_PARAM_2, [kernelParams[1]]),
    DG.Column.fromList('double', TRAIN_ERROR, [model.trainError]), 
  ]);     
} 

// Show trained SVM-model: short from
export function showModel(model) {
  let info = getModelInfo(model);
  info.name = MODEL_INFO;
  grok.shell.addTableView(info);
}

// Show trained SVM-model: full info
export function showModelFullInfo(model) {

  // get model info
  let info = getModelInfo(model);
  info.name = MODEL_INFO_FULL;
  grok.shell.addTableView(info);
    
  // show normalized data
  //let normalizedData = DG.DataFrame.fromColumns(model.normalizedTrainData);
  let normalizedData = model.normalizedTrainData;
  normalizedData.name = NORMALIZED;    
  grok.shell.addTableView(normalizedData);

  // show characteristics
  let characteristics = DG.DataFrame.fromColumns([model.means, model.stdDevs]);
  characteristics.name = CHARACTERISTICS;
  grok.shell.addTableView(characteristics);    

  //show model params  
  let modelParamsDF = DG.DataFrame.fromColumns([model.modelParams]);
  modelParamsDF.name = PARAMETERS;
  grok.shell.addTableView(modelParamsDF);  

  //show model weights  
  let modelWeightsDF = DG.DataFrame.fromColumns([model.modelWeights]);
  modelWeightsDF.name = WEIGHTS;
  grok.shell.addTableView(modelWeightsDF);  

  //show model labels & prediciotns 
  let modelLabelsDF = DG.DataFrame.fromColumns([model.trainLabels, 
    model.predictedLabels]);
  modelLabelsDF.name = LABELS;
  grok.shell.addTableView(modelLabelsDF);  
} // showModelFullInfo
