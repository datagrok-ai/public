import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
}

//tags: init
//meta.role: init
export async function init() : Promise<void> {
  await PackageFunctions.init();
}

//name: DBSCAN
//description: Density-based spatial clustering of applications with noise (DBSCAN)
//input: dataframe df 
//input: column xCol { type: numerical }
//input: column yCol { type: numerical }
//input: double epsilon = 0.02 { caption: Epsilon; description: The maximum distance between two samples for them to be considered as in the same neighborhood. }
//input: int minPts = 4 { caption: Minimum points; description: The number of samples (or total weight) in a neighborhood for a point to be considered as a core point. }
//output: column result
//top-menu: ML | Cluster | DBSCAN...
export async function dbScan(df: DG.DataFrame, xCol: DG.Column, yCol: DG.Column, epsilon: number, minPts: number) : Promise<any> {
  return await PackageFunctions.dbScan(df, xCol, yCol, epsilon, minPts);
}

//description: Principal component analysis (PCA).
//input: dataframe table { caption: Table }
//input: column_list features { type: numerical; nullable: false }
//input: int components = 2 { showPlusMinus: true; caption: Components; nullable: false; min: 1; description: Number of components. }
//input: bool center = false { caption: Center; description: Indicating whether the variables should be shifted to be zero centered. }
//input: bool scale = false { caption: Scale; description: Indicating whether the variables should be scaled to have unit variance. }
//top-menu: ML | Analyze | PCA...
//help-url: /help/explore/dim-reduction#pca
export async function PCA(table: DG.DataFrame, features: DG.ColumnList, components: number, center: boolean, scale: boolean) : Promise<void> {
  await PackageFunctions.PCA(table, features, components, center, scale);
}

//name: DBSCAN clustering
//tags: dim-red-postprocessing-function
//input: column col1 
//input: column col2 
//input: double epsilon = 0.01 { description: Minimum distance between two points to be considered as in the same neighborhood. }
//input: int minimumPoints = 5 { description: Minimum number of points to form a dense region. }
//meta.defaultPostProcessingFunction: true
//meta.role: dim-red-postprocessing-function
export async function dbscanPostProcessingFunction(col1: DG.Column, col2: DG.Column, epsilon: number, minimumPoints: number) : Promise<void> {
  await PackageFunctions.dbscanPostProcessingFunction(col1, col2, epsilon, minimumPoints);
}

//name: None (number)
//tags: dim-red-preprocessing-function
//input: column col 
//input: string _metric { optional: true }
//output: object result
//meta.supportedTypes: int,float,double,qnum
//meta.supportedDistanceFunctions: Difference
//meta.role: dim-red-preprocessing-function
export function numberPreprocessingFunction(col: DG.Column, _metric: string) {
  return PackageFunctions.numberPreprocessingFunction(col, _metric);
}

//name: None (string)
//tags: dim-red-preprocessing-function
//input: column col 
//input: string _metric { optional: true }
//output: object result
//meta.supportedTypes: string
//meta.supportedDistanceFunctions: One-Hot,Levenshtein,Hamming
//meta.role: dim-red-preprocessing-function
export function stringPreprocessingFunction(col: DG.Column, _metric: string) {
  return PackageFunctions.stringPreprocessingFunction(col, _metric);
}

//name: Multi Column Dimensionality Reduction
//top-menu: ML | Reduce Dimensionality...
export async function reduceDimensionality() : Promise<void> {
  await PackageFunctions.reduceDimensionality();
}

//tags: editor
//input: funccall call 
//meta.role: editor
export function GetMCLEditor(call: DG.FuncCall) : void {
  PackageFunctions.GetMCLEditor(call);
}

//description: Markov clustering (MCL) is an unsupervised clustering algorithm for graphs based on simulation of stochastic flow.
//input: dataframe df 
//input: list<column> cols 
//input: list<string> metrics 
//input: list<double> weights 
//input: string aggregationMethod 
//input: list<func> preprocessingFuncs 
//input: object preprocessingFuncArgs 
//input: int threshold = 80 
//input: int maxIterations = 10 
//input: bool useWebGPU = false 
//input: double inflate = 2 
//input: int minClusterSize = 5 
//top-menu: ML | Cluster | MCL...
//editor: EDA:GetMCLEditor
export async function MCLClustering(df: DG.DataFrame, cols: DG.Column[], metrics: any, weights: number[], aggregationMethod: any, preprocessingFuncs: any[], preprocessingFuncArgs: any[], threshold: number, maxIterations: number, useWebGPU: boolean, inflate: number, minClusterSize: number) : Promise<any> {
  return await PackageFunctions.MCLClustering(df, cols, metrics, weights, aggregationMethod, preprocessingFuncs, preprocessingFuncArgs, threshold, maxIterations, useWebGPU, inflate, minClusterSize);
}

//name: MCL
//description: Markov clustering viewer
//tags: viewer
//output: viewer result
//meta.showInGallery: false
//meta.role: viewer
export function markovClusteringViewer() : any {
  return PackageFunctions.markovClusteringViewer();
}

//description: Compute partial least squares (PLS) regression analysis components: prediction, regression coefficients, T- & U-scores, X-loadings.
//input: dataframe table 
//input: column_list features { type: numerical }
//input: column predict { type: numerical }
//input: int components = 3 { description: Number of latent factors the model extracts from the predictors. }
//input: column names { type: string }
//output: object plsResults
export async function PLS(table: DG.DataFrame, features: DG.ColumnList, predict: DG.Column, components: number, names: DG.Column) : Promise<any> {
  return await PackageFunctions.PLS(table, features, predict, components, names);
}

//description: Compute partial least squares (PLS) regression components. They maximally summarize the variation of the predictors while maximizing correlation with the response variable.
//top-menu: ML | Analyze | PLS...
export async function topMenuPLS() : Promise<void> {
  await PackageFunctions.topMenuPLS();
}

//name: multivariateAnalysis
//description: Multidimensional data analysis using partial least squares (PLS) regression.
//top-menu: ML | Analyze | Multivariate Analysis...
export async function MVA() : Promise<void> {
  await PackageFunctions.MVA();
}

//name: MVA demo
//description: Multidimensional data analysis using partial least squares (PLS) regression. It identifies latent factors and constructs a linear model based on them.
//meta.demoPath: Compute | Multivariate Analysis
export async function demoMultivariateAnalysis() : Promise<void> {
  await PackageFunctions.demoMultivariateAnalysis();
}

//name: T-test
//description: Two-sample t-test (Welch or Student) compares the means of a feature between two groups.
//top-menu: ML | Analyze | Group Comparison | T-test...
export function tTest() : void {
  PackageFunctions.tTest();
}

//name: ANOVA
//description: One-way analysis of variances (ANOVA) determines whether the examined factor has a significant impact on the explored feature.
//top-menu: ML | Analyze | Group Comparison | ANOVA...
export function anova() : void {
  PackageFunctions.anova();
}

//name: Control comparisons
//description: Compare several groups against a single control (Dunnett's test or Holm-corrected Welch's t-tests).
//top-menu: ML | Analyze | Group Comparison | Control Comparisons...
export function controlComparisons() : void {
  PackageFunctions.controlComparisons();
}

//name: KNN impute
//description: Missing values imputation using the k-nearest neighbors method (KNN)
//top-menu: ML | Impute Missing Values...
export function kNNImputation() : void {
  PackageFunctions.kNNImputation();
}

//name: KNN imputation for a table
//description: Missing values imputation using the k-nearest neighbors method
//input: dataframe table 
export async function kNNImputationForTable(table: DG.DataFrame) : Promise<void> {
  await PackageFunctions.kNNImputationForTable(table);
}

//input: dataframe df 
//input: column predictColumn 
//input: double rate = 0.1 { caption: Rate; min: 0; max: 10; step: 0.01; description: Gradient descent learning rate. }
//input: int iterations = 1000 { caption: Iterations; min: 1; step: 50; max: 10000; description: Largest number of training steps before training stops. }
//input: double alpha = 0 { caption: L1; min: 0; max: 100; description: L1 (Lasso) regularization term. 0 means plain ordinary least squares. }
//input: double lambda = 0 { caption: L2; min: 0; max: 100; description: L2 (Ridge) regularization term. 0 means plain ordinary least squares. }
//output: dynamic model
//meta.mlname: Linear Regression
//meta.mlrole: train
export async function trainLinearRegression(df: DG.DataFrame, predictColumn: DG.Column, rate: number, iterations: number, alpha: number, lambda: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainLinearRegression(df, predictColumn, rate, iterations, alpha, lambda);
}

//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: Linear Regression
//meta.mlrole: apply
export function applyLinearRegression(df: DG.DataFrame, model: any) : any {
  return PackageFunctions.applyLinearRegression(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Linear Regression
//meta.mlrole: isApplicable
export function isApplicableLinearRegression(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableLinearRegression(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Linear Regression
//meta.mlrole: isInteractive
export function isInteractiveLinearRegression(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveLinearRegression(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: double rate = 1.0 { min: 0.001; max: 20; description: Gradient descent learning rate. }
//input: double iterations = 100 { min: 1; max: 10000; step: 10; description: Largest number of training steps before training stops. }
//input: double penalty = 0.1 { min: 0.0001; max: 1; description: Regularization rate. }
//input: double tolerance = 0.001 { min: 0.00001; max: 0.1; description: Smallest improvement worth continuing training for. }
//output: dynamic model
//meta.mlname: Softmax
//meta.mlrole: train
export async function trainSoftmax(df: DG.DataFrame, predictColumn: DG.Column, rate: number, iterations: number, penalty: number, tolerance: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainSoftmax(df, predictColumn, rate, iterations, penalty, tolerance);
}

//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: Softmax
//meta.mlrole: apply
export function applySoftmax(df: DG.DataFrame, model: any) : any {
  return PackageFunctions.applySoftmax(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Softmax
//meta.mlrole: isApplicable
export function isApplicableSoftmax(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableSoftmax(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Softmax
//meta.mlrole: isInteractive
export function isInteractiveSoftmax(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveSoftmax(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: int components = 3 { min: 1; max: 10; description: Number of latent components. }
//output: dynamic model
//meta.mlname: PLS Regression
//meta.mlrole: train
export async function trainPLSRegression(df: DG.DataFrame, predictColumn: DG.Column, components: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainPLSRegression(df, predictColumn, components);
}

//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: PLS Regression
//meta.mlrole: apply
export function applyPLSRegression(df: DG.DataFrame, model: any) : any {
  return PackageFunctions.applyPLSRegression(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: PLS Regression
//meta.mlrole: isApplicable
export function isApplicablePLSRegression(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicablePLSRegression(df, predictColumn);
}

//input: dataframe df 
//input: column targetColumn 
//input: column predictColumn 
//input: dynamic model 
//output: dynamic result
//meta.mlname: PLS Regression
//meta.mlrole: visualize
export async function visualizePLSRegression(df: DG.DataFrame, targetColumn: DG.Column, predictColumn: DG.Column, model: any) : Promise<any> {
  return await PackageFunctions.visualizePLSRegression(df, targetColumn, predictColumn, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: PLS Regression
//meta.mlrole: isInteractive
export function isInteractivePLSRegression(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractivePLSRegression(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: int iterations = 20 { min: 1; max: 100; description: Number of training iterations. }
//input: double eta = 0.3 { caption: Rate; min: 0; max: 1; description: Learning rate. }
//input: int maxDepth = 6 { min: 0; max: 20; description: Maximum depth of a tree. }
//input: double lambda = 1 { min: 0; max: 100; description: L2 regularization term. }
//input: double alpha = 0 { min: 0; max: 100; description: L1 regularization term. }
//output: dynamic model
//meta.mlname: XGBoost
//meta.mlrole: train
export async function trainXGBooster(df: DG.DataFrame, predictColumn: DG.Column, iterations: number, eta: number, maxDepth: number, lambda: number, alpha: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainXGBooster(df, predictColumn, iterations, eta, maxDepth, lambda, alpha);
}

//input: dataframe df 
//input: dynamic model 
//output: dataframe result
//meta.mlname: XGBoost
//meta.mlrole: apply
export function applyXGBooster(df: DG.DataFrame, model: any) : any {
  return PackageFunctions.applyXGBooster(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: XGBoost
//meta.mlrole: isInteractive
export function isInteractiveXGBooster(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveXGBooster(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: XGBoost
//meta.mlrole: isApplicable
export function isApplicableXGBooster(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableXGBooster(df, predictColumn);
}

//name: Pareto Front
//description: Perform optimization across multiple objectives: analyze trade-offs between conflicting objectives and identify Pareto-optimal points.
//top-menu: ML | Pareto Front...
export function paretoFront() : void {
  PackageFunctions.paretoFront();
}

//name: Pareto front
//description: Pareto front viewer
//tags: viewer
//output: viewer result
//meta.icon: icons/pareto-front-viewer.svg
//meta.role: viewer
export function paretoFrontViewer() : any {
  return PackageFunctions.paretoFrontViewer();
}

//description: Train probabilistic multi-parameter optimization (pMPO) model
export function trainPmpo() : void {
  PackageFunctions.trainPmpo();
}

//input: view view 
//output: object result
export async function getPmpoAppItems(view: any) : Promise<any> {
  return await PackageFunctions.getPmpoAppItems(view);
}

//description: Generates syntethetic dataset oriented on the pMPO modeling
//input: int samples 
//output: dataframe Synthetic
export async function generatePmpoDataset(samples: number) : Promise<any> {
  return await PackageFunctions.generatePmpoDataset(samples);
}
