/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();

// Principal component analysis (PCA) imports
import { performPCA } from './pca'

// Partial least-squares regression (PLS) imports
import { performPLS} from './pls'

// Support vector machine (SVM) tools imports
import { showTrainReport, getPackedModel, getTrainedModel, getPrediction,
  LINEAR, RBF, POLYNOMIAL, SIGMOID } from './svm';

// Data generation tools
import { generateDatasetForTestingSVM } from './generators';

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init() {
  await initEDALib();
}

//top-menu: Tools | Data Science | PCA
//name: PCA
//input: dataframe table
//input: column_list features
//input: int components = 2
export function PCA(table, features, components) {  
  performPCA(table, features, components);   
}

//top-menu: Tools | Data Science | PLS
//name: PLS
//input: dataframe table
//input: column_list features
//input: column predict
//input: int components = 3
export function PLS(table, features, predict, components) {
  performPLS(table, features, predict, components);  
}

//name: Generate test data (linear kernel case)
//description: Generates dataset for testing SVN with linear kernel.
//input: string name = 'Data' {caption: name; category: Dataset}
//input: int samplesCount = 1000 {caption: samples; category: Size}
//input: int featuresCount = 2 {caption: features; category: Size}
//input: double min = -39 {caption: min; category: Range}
//input: double max = 173 {caption: max; category: Range}
//input: double violatorsPercentage = 5 {caption: violators; units: %; category: Dataset}
//output: dataframe df
export function generateDatasetLinear(name, samplesCount, featuresCount, 
  min, max, violatorsPercentage) 
{
  return generateDatasetForTestingSVM(LINEAR, [0, 0], name, samplesCount, featuresCount,
    min, max, violatorsPercentage);
}

//name: Generate test data (RBF kernel case)
//description: Generates dataset for testing SVN with RBF-kernel.
//input: string name = 'Data' {caption: name; category: Dataset}
//input: double sigma = 90  {caption: sigma; category: Hyperparameters}
//input: int samplesCount = 1000 {caption: samples; category: Size}
//input: int featuresCount = 2 {caption: features; category: Size}
//input: double min = -39 {caption: min; category: Range}
//input: double max = 173 {caption: max; category: Range}
//input: double violatorsPercentage = 5 {caption: violators; units: %; category: Dataset}
//output: dataframe df
export function generateDatasetRBF(name, sigma, samplesCount, featuresCount, 
  min, max, violatorsPercentage) 
{
  return generateDatasetForTestingSVM(RBF, [sigma, 0], name, samplesCount, featuresCount,
    min, max, violatorsPercentage);
} 

//name: trainLinearKernelSVM
//meta.mlname: linear kernel LS-SVM
//meta.mlrole: train
//input: dataframe df
//input: string predict_column
//input: double gamma = 1.0 {category: Hyperparameters}
//input: bool toShowReport = false {caption: to show report; category: Report}
//output: dynamic model
export function trainLinearKernelSVM(df, predict_column, gamma, toShowReport) {    

  let trainedModel = getTrainedModel({gamma: gamma, kernel: LINEAR}, df, predict_column);   

  if(toShowReport)
    showTrainReport(df, trainedModel);

  return getPackedModel(trainedModel);
} // trainLinearKernelSVM

//name: applyLinearKernelSVM
//meta.mlname: linear kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export function applyLinearKernelSVM(df, model) { return getPrediction(df, model); }

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

  let trainedModel = getTrainedModel(
    {gamma: gamma, kernel: RBF, sigma: sigma}, 
    df, predict_column);   

  if(toShowReport)
    showTrainReport(df, trainedModel);

  return getPackedModel(trainedModel);
} // trainRBFkernelSVM

//name: applyRBFkernelSVM
//meta.mlname: RBF-kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export function applyRBFkernelSVM(df, model) { return getPrediction(df, model); } 

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

  let trainedModel = getTrainedModel(
    {gamma: gamma, kernel: POLYNOMIAL, cParam: c, dParam: d}, 
    df, predict_column);   

  if(toShowReport)
    showTrainReport(df, trainedModel);

  return getPackedModel(trainedModel);
} // trainPolynomialKernelSVM

//name: applyPolynomialKernelSVM
//meta.mlname: polynomial kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export function applyPolynomialKernelSVM(df, model) { return getPrediction(df, model); }

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

  let trainedModel = getTrainedModel(
    {gamma: gamma, kernel: SIGMOID, kappa: kappa, theta: theta}, 
    df, predict_column);   

  if(toShowReport)
    showTrainReport(df, trainedModel);

  return getPackedModel(trainedModel);
} // trainSigmoidKernelSVM

//name: applySigmoidKernelSVM
//meta.mlname: sigmoid kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export function applySigmoidKernelSVM(df, model) { return getPrediction(df, model); }

