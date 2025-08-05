import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
//output: dynamic result
export function info() {
  return PackageFunctions.info();
}

//name: init
//tags: init
export async function init() {
  return PackageFunctions.init();
}

//name: DBSCAN
//description: Density-based spatial clustering of applications with noise (DBSCAN)
//input: dataframe df 
//input: column xCol { type: numerical }
//input: column yCol { type: numerical }
//input: double epsilon { caption: Epsilon; default: 0.02; description: The maximum distance between two samples for them to be considered as in the same neighborhood. }
//input: int minPts { caption: Minimum points; default: 4; description: The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. }
//output: column result
//top-menu: ML | Cluster | DBSCAN...
export async function dbScan(df: DG.DataFrame, xCol: DG.Column, yCol: DG.Column, epsilon: number, minPts: number) {
  return PackageFunctions.dbScan(df, xCol, yCol, epsilon, minPts);
}

//name: PCA
//description: Principal component analysis (PCA)
//input: dataframe table 
//input: column_list features { type: numerical; nullable: false }
//input: int components { caption: Components; nullable: false; min: 1; default: 2; description: Number of components. }
//input: bool center { default: false; description: Indicating whether the variables should be shifted to be zero centered. }
//input: bool scale { default: false; description: Indicating whether the variables should be scaled to have unit variance. }
//top-menu: ML | Analyze | PCA...
export async function PCA(table: DG.DataFrame, features: DG.ColumnList, components: number, center: boolean, scale: boolean) {
  return PackageFunctions.PCA(table, features, components, center, scale);
}

//name: DBSCAN clustering
//tags: dim-red-postprocessing-function
//input: column col1 
//input: column col2 
//input: double epsilon { default: 0.01; description: Minimum distance between two points to be considered as in the same neighborhood. }
//input: int minimumPoints { default: 5; description: Minimum number of points to form a dense region. }
//output: dynamic result
//meta.defaultPostProcessingFunction: true
export async function dbscanPostProcessingFunction(col1: DG.Column, col2: DG.Column, epsilon: number, minimumPoints: number) {
  return PackageFunctions.dbscanPostProcessingFunction(col1, col2, epsilon, minimumPoints);
}

//name: None (number)
//tags: dim-red-preprocessing-function
//input: column col 
//input: string _metric { optional: true }
//output: dynamic result
//meta.supportedTypes: int,float,double,qnum
//meta.supportedDistanceFunctions: Difference
export function numberPreprocessingFunction(col: DG.Column, _metric: string) {
  return PackageFunctions.numberPreprocessingFunction(col, _metric);
}

//name: None (string)
//tags: dim-red-preprocessing-function
//input: column col 
//input: string _metric { optional: true }
//output: dynamic result
//meta.supportedTypes: string
//meta.supportedDistanceFunctions: One-Hot,Levenshtein,Hamming
export function stringPreprocessingFunction(col: DG.Column, _metric: string) {
  return PackageFunctions.stringPreprocessingFunction(col, _metric);
}

//name: Multi Column Dimensionality Reduction
//top-menu: ML | Reduce Dimensionality...
export async function reduceDimensionality() {
  return PackageFunctions.reduceDimensionality();
}

//name: GetMCLEditor
//tags: editor
//input: funccall call 
export function GetMCLEditor(call: DG.FuncCall) {
  return PackageFunctions.GetMCLEditor(call);
}

//name: MCLClustering
//description: Markov clustering (MCL) is an unsupervised clustering algorithm for graphs based on simulation of stochastic flow.
//input: dataframe df 
//input: list<column> cols 
//input: list<string> metrics 
//input: list<double> weights 
//input: string aggregationMethod 
//input: list<func> preprocessingFuncs 
//input: object preprocessingFuncArgs 
//input: int threshold { default: 80 }
//input: int maxIterations { default: 10 }
//input: bool useWebGPU { default: false }
//input: double inflate { default: 2 }
//input: int minClusterSize { default: 5 }
//output: dynamic result
//top-menu: ML | Cluster | MCL...
//editor: EDA: GetMCLEditor
export async function MCLClustering(df: DG.DataFrame, cols: DG.Column[], metrics: any, weights: number[], aggregationMethod: any, preprocessingFuncs: any[], preprocessingFuncArgs: any[], threshold: number, maxIterations: number, useWebGPU: boolean, inflate: number, minClusterSize: number) {
  return PackageFunctions.MCLClustering(df, cols, metrics, weights, aggregationMethod, preprocessingFuncs, preprocessingFuncArgs, threshold, maxIterations, useWebGPU, inflate, minClusterSize);
}

//name: MCL
//description: Markov clustering viewer
//tags: viewer
//output: viewer result
export function markovClusteringViewer() {
  return PackageFunctions.markovClusteringViewer();
}

//name: PLS
//description: Compute partial least squares (PLS) regression analysis components: prediction, regression coefficients, T- & U-scores, X-loadings.
//input: dataframe table 
//input: column_list features { type: numerical }
//input: column predict { type: numerical }
//input: int components { default: 3 }
//input: column names { type: string }
//output: object plsResults
export async function PLS(table: DG.DataFrame, features: DG.ColumnList, predict: DG.Column, components: number, names: DG.Column) {
  return PackageFunctions.PLS(table, features, predict, components, names);
}

//name: topMenuPLS
//description: Compute partial least squares (PLS) regression components. They maximally summarize the variation of the predictors while maximizing correlation with the response variable.
//top-menu: ML | Analyze | PLS...
export async function topMenuPLS() {
  return PackageFunctions.topMenuPLS();
}

//name: multivariateAnalysis
//description: Multidimensional data analysis using partial least squares (PLS) regression.
//top-menu: ML | Analyze | Multivariate Analysis...
export async function MVA() {
  return PackageFunctions.MVA();
}

//name: MVA demo
//description: Multidimensional data analysis using partial least squares (PLS) regression. It identifies latent factors and constructs a linear model based on them.
//output: dynamic result
//meta.demoPath: Compute | Multivariate Analysis
export async function demoMultivariateAnalysis() {
  return PackageFunctions.demoMultivariateAnalysis();
}

//name: trainLinearKernelSVM
//input: dataframe df 
//input: column predictColumn 
//input: double gamma { category: Hyperparameters; default: 1.0 }
//output: dynamic result
//meta.mlname: linear kernel LS-SVM
//meta.mlrole: train
export async function trainLinearKernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number) {
  return PackageFunctions.trainLinearKernelSVM(df, predictColumn, gamma);
}

//name: applyLinearKernelSVM
//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: linear kernel LS-SVM
//meta.mlrole: apply
export async function applyLinearKernelSVM(df: DG.DataFrame, model: any) {
  return PackageFunctions.applyLinearKernelSVM(df, model);
}

//name: isApplicableLinearKernelSVM
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: linear kernel LS-SVM
//meta.mlrole: isApplicable
export async function isApplicableLinearKernelSVM(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isApplicableLinearKernelSVM(df, predictColumn);
}

//name: isInteractiveLinearKernelSVM
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: linear kernel LS-SVM
//meta.mlrole: isInteractive
export async function isInteractiveLinearKernelSVM(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isInteractiveLinearKernelSVM(df, predictColumn);
}

//name: visualizeLinearKernelSVM
//input: dataframe df 
//input: column targetColumn 
//input: column predictColumn 
//input: dynamic model 
//output: dynamic result
//meta.mlname: linear kernel LS-SVM
//meta.mlrole: visualize
export async function visualizeLinearKernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any) {
  return PackageFunctions.visualizeLinearKernelSVM(df, targetColumn, predictColumn, model);
}

//name: trainRBFkernelSVM
//input: dataframe df 
//input: column predictColumn 
//input: double gamma { category: Hyperparameters; default: 1.0 }
//input: double sigma { category: Hyperparameters; default: 1.5 }
//output: dynamic result
//meta.mlname: RBF-kernel LS-SVM
//meta.mlrole: train
export async function trainRBFkernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number, sigma: number) {
  return PackageFunctions.trainRBFkernelSVM(df, predictColumn, gamma, sigma);
}

//name: applyRBFkernelSVM
//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: RBF-kernel LS-SVM
//meta.mlrole: apply
export async function applyRBFkernelSVM(df: DG.DataFrame, model: any) {
  return PackageFunctions.applyRBFkernelSVM(df, model);
}

//name: isApplicableRBFkernelSVM
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: RBF-kernel LS-SVM
//meta.mlrole: isApplicable
export async function isApplicableRBFkernelSVM(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isApplicableRBFkernelSVM(df, predictColumn);
}

//name: isInteractiveRBFkernelSVM
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: RBF-kernel LS-SVM
//meta.mlrole: isInteractive
export async function isInteractiveRBFkernelSVM(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isInteractiveRBFkernelSVM(df, predictColumn);
}

//name: visualizeRBFkernelSVM
//input: dataframe df 
//input: column targetColumn 
//input: column predictColumn 
//input: dynamic model 
//output: dynamic result
//meta.mlname: RBF-kernel LS-SVM
//meta.mlrole: visualize
export async function visualizeRBFkernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any) {
  return PackageFunctions.visualizeRBFkernelSVM(df, targetColumn, predictColumn, model);
}

//name: trainPolynomialKernelSVM
//input: dataframe df 
//input: column predictColumn 
//input: double gamma { category: Hyperparameters; default: 1.0 }
//input: double c { category: Hyperparameters; default: 1 }
//input: double d { category: Hyperparameters; default: 2 }
//output: dynamic result
//meta.mlname: polynomial kernel LS-SVM
//meta.mlrole: train
export async function trainPolynomialKernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number, c: number, d: number) {
  return PackageFunctions.trainPolynomialKernelSVM(df, predictColumn, gamma, c, d);
}

//name: applyPolynomialKernelSVM
//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: polynomial kernel LS-SVM
//meta.mlrole: apply
export async function applyPolynomialKernelSVM(df: DG.DataFrame, model: any) {
  return PackageFunctions.applyPolynomialKernelSVM(df, model);
}

//name: isApplicablePolynomialKernelSVM
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: polynomial kernel LS-SVM
//meta.mlrole: isApplicable
export async function isApplicablePolynomialKernelSVM(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isApplicablePolynomialKernelSVM(df, predictColumn);
}

//name: isInteractivePolynomialKernelSVM
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: polynomial kernel LS-SVM
//meta.mlrole: isInteractive
export async function isInteractivePolynomialKernelSVM(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isInteractivePolynomialKernelSVM(df, predictColumn);
}

//name: visualizePolynomialKernelSVM
//input: dataframe df 
//input: column targetColumn 
//input: column predictColumn 
//input: dynamic model 
//output: dynamic widget
//meta.mlname: polynomial kernel LS-SVM
//meta.mlrole: visualize
export async function visualizePolynomialKernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any) {
  return PackageFunctions.visualizePolynomialKernelSVM(df, targetColumn, predictColumn, model);
}

//name: trainSigmoidKernelSVM
//input: dataframe df 
//input: column predictColumn 
//input: double gamma { category: Hyperparameters; default: 1.0 }
//input: double kappa { category: Hyperparameters; default: 1 }
//input: double theta { category: Hyperparameters; default: 1 }
//output: dynamic result
//meta.mlname: sigmoid kernel LS-SVM
//meta.mlrole: train
export async function trainSigmoidKernelSVM(df: DG.DataFrame, predictColumn: DG.Column, gamma: number, kappa: number, theta: number) {
  return PackageFunctions.trainSigmoidKernelSVM(df, predictColumn, gamma, kappa, theta);
}

//name: applySigmoidKernelSVM
//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: sigmoid kernel LS-SVM
//meta.mlrole: apply
export async function applySigmoidKernelSVM(df: DG.DataFrame, model: any) {
  return PackageFunctions.applySigmoidKernelSVM(df, model);
}

//name: isApplicableSigmoidKernelSVM
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: sigmoid kernel LS-SVM
//meta.mlrole: isApplicable
export async function isApplicableSigmoidKernelSVM(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isApplicableSigmoidKernelSVM(df, predictColumn);
}

//name: isInteractiveSigmoidKernelSVM
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: sigmoid kernel LS-SVM
//meta.mlrole: isInteractive
export async function isInteractiveSigmoidKernelSVM(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isInteractiveSigmoidKernelSVM(df, predictColumn);
}

//name: visualizeSigmoidKernelSVM
//input: dataframe df 
//input: column targetColumn 
//input: column predictColumn 
//input: dynamic model 
//output: dynamic result
//meta.mlname: sigmoid kernel LS-SVM
//meta.mlrole: visualize
export async function visualizeSigmoidKernelSVM(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any) {
  return PackageFunctions.visualizeSigmoidKernelSVM(df, targetColumn, predictColumn, model);
}

//name: ANOVA
//description: One-way analysis of variances (ANOVA) determines whether the examined factor has a significant impact on the explored feature.
//top-menu: ML | Analyze | ANOVA...
export function anova() {
  return PackageFunctions.anova();
}

//name: KNN impute
//description: Missing values imputation using the k-nearest neighbors method (KNN)
//output: dynamic result
//top-menu: ML | Impute Missing Values...
export function kNNImputation() {
  return PackageFunctions.kNNImputation();
}

//name: KNN imputation for a table
//description: Missing values imputation using the k-nearest neighbors method
//input: dataframe table 
//output: dynamic result
export async function kNNImputationForTable(table: DG.DataFrame) {
  return PackageFunctions.kNNImputationForTable(table);
}

//name: trainLinearRegression
//input: dataframe df 
//input: column predictColumn 
//output: blob result
//meta.mlname: Linear Regression
//meta.mlrole: train
export async function trainLinearRegression(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.trainLinearRegression(df, predictColumn);
}

//name: applyLinearRegression
//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: Linear Regression
//meta.mlrole: apply
export function applyLinearRegression(df: DG.DataFrame, model: any) {
  return PackageFunctions.applyLinearRegression(df, model);
}

//name: isApplicableLinearRegression
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Linear Regression
//meta.mlrole: isApplicable
export function isApplicableLinearRegression(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isApplicableLinearRegression(df, predictColumn);
}

//name: isInteractiveLinearRegression
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Linear Regression
//meta.mlrole: isInteractive
export function isInteractiveLinearRegression(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isInteractiveLinearRegression(df, predictColumn);
}

//name: trainSoftmax
//input: dataframe df 
//input: column predictColumn 
//input: double rate { category: Hyperparameters; default: 1.0; min: 0.001; max: 20; description: Learning rate. }
//input: double iterations { category: Hyperparameters; default: 100; min: 1; max: 10000; step: 10; description: Fitting iterations count }
//input: double penalty { category: Hyperparameters; default: 0.1; min: 0.0001; max: 1; description: Regularization rate. }
//input: double tolerance { category: Hyperparameters; default: 0.001; min: 0.00001; max: 0.1; description: Fitting tolerance. }
//output: blob result
//meta.mlname: Softmax
//meta.mlrole: train
export async function trainSoftmax(df: DG.DataFrame, predictColumn: DG.Column, rate: number, iterations: number, penalty: number, tolerance: number) {
  return PackageFunctions.trainSoftmax(df, predictColumn, rate, iterations, penalty, tolerance);
}

//name: applySoftmax
//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: Softmax
//meta.mlrole: apply
export function applySoftmax(df: DG.DataFrame, model: any) {
  return PackageFunctions.applySoftmax(df, model);
}

//name: isApplicableSoftmax
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Softmax
//meta.mlrole: isApplicable
export function isApplicableSoftmax(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isApplicableSoftmax(df, predictColumn);
}

//name: isInteractiveSoftmax
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Softmax
//meta.mlrole: isInteractive
export function isInteractiveSoftmax(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isInteractiveSoftmax(df, predictColumn);
}

//name: trainPLSRegression
//input: dataframe df 
//input: column predictColumn 
//input: int components { min: 1; max: 10; default: 3; description: Number of latent components. }
//output: blob result
//meta.mlname: PLS Regression
//meta.mlrole: train
export async function trainPLSRegression(df: DG.DataFrame, predictColumn: DG.Column, components: number) {
  return PackageFunctions.trainPLSRegression(df, predictColumn, components);
}

//name: applyPLSRegression
//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: PLS Regression
//meta.mlrole: apply
export function applyPLSRegression(df: DG.DataFrame, model: any) {
  return PackageFunctions.applyPLSRegression(df, model);
}

//name: isApplicablePLSRegression
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: PLS Regression
//meta.mlrole: isApplicable
export function isApplicablePLSRegression(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isApplicablePLSRegression(df, predictColumn);
}

//name: visualizePLSRegression
//input: dataframe df 
//input: column targetColumn 
//input: column predictColumn 
//input: dynamic model 
//output: dynamic result
//meta.mlname: PLS Regression
//meta.mlrole: visualize
export async function visualizePLSRegression(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any) {
  return PackageFunctions.visualizePLSRegression(df, targetColumn, predictColumn, model);
}

//name: isInteractivePLSRegression
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: PLS Regression
//meta.mlrole: isInteractive
export function isInteractivePLSRegression(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isInteractivePLSRegression(df, predictColumn);
}

//name: trainXGBooster
//input: dataframe df 
//input: column predictColumn 
//input: int iterations { min: 1; max: 100; default: 20; description: Number of training iterations. }
//input: double eta { caption: Rate; min: 0; max: 1; default: 0.3; description: Learning rate. }
//input: int maxDepth { min: 0; max: 20; default: 6; description: Maximum depth of a tree. }
//input: double lambda { min: 0; max: 100; default: 1; description: L2 regularization term. }
//input: double alpha { min: 0; max: 100; default: 0; description: L1 regularization term. }
//output: blob result
//meta.mlname: XGBoost
//meta.mlrole: train
export async function trainXGBooster(df: DG.DataFrame, predictColumn: DG.Column, iterations: number, eta: number, maxDepth: number, lambda: number, alpha: number) {
  return PackageFunctions.trainXGBooster(df, predictColumn, iterations, eta, maxDepth, lambda, alpha);
}

//name: applyXGBooster
//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: XGBoost
//meta.mlrole: apply
export function applyXGBooster(df: DG.DataFrame, model: any) {
  return PackageFunctions.applyXGBooster(df, model);
}

//name: isInteractiveXGBooster
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: XGBoost
//meta.mlrole: isInteractive
export function isInteractiveXGBooster(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isInteractiveXGBooster(df, predictColumn);
}

//name: isApplicableXGBooster
//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: XGBoost
//meta.mlrole: isApplicable
export function isApplicableXGBooster(df: DG.DataFrame, predictColumn: DG.Column) {
  return PackageFunctions.isApplicableXGBooster(df, predictColumn);
}
