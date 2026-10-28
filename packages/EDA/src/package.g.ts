import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: info
export function info() : void {
  PackageFunctions.info();
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

//name: MCL Clustering
//description: Markov clustering (MCL) is an unsupervised clustering algorithm for graphs based on simulation of stochastic flow.
//input: dataframe df 
//input: list<column> cols 
//input: list<string> metrics { description: Distance metric per column used to measure similarity between rows. }
//input: list<double> weights { description: Relative weight of each column when combining per-column distances. }
//input: string aggregationMethod { description: How per-column distances are aggregated into a single distance. }
//input: list<func> preprocessingFuncs { description: Preprocessing function applied to each column before distances are computed. }
//input: object preprocessingFuncArgs { description: Arguments passed to the preprocessing functions. }
//input: int threshold = 80 { description: Similarity threshold (percentile): edges below it are dropped before clustering. }
//input: int maxIterations = 10 { description: Maximum number of expansion/inflation iterations. }
//input: bool useWebGPU = false { description: Run the computation on the GPU via WebGPU when available. }
//input: double inflate = 2 { description: Inflation factor controlling cluster granularity: higher values yield more, smaller clusters. }
//input: int minClusterSize = 5 { description: Clusters smaller than this are merged into noise. }
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

//description: Add the PLS-based multivariate analysis results to the table.
//input: dataframe table 
//input: list<string> featureNames 
//input: column predict { type: numerical }
//input: int components { description: Number of latent factors the model extracts from the predictors. }
//input: bool isQuadratic { description: Include squared terms as additional predictors. }
//input: bool componentsOnly { description: Add just the PLS components. }
//input: list<string> xScoreNames 
//input: list<string> yScoreNames 
//input: string predictionName 
//input: string analysisTableName 
//input: string explVarTableName 
//meta.role: transform
export async function multivariateAnalysisTransform(table: DG.DataFrame, featureNames: string[], predict: DG.Column, components: number, isQuadratic: boolean, componentsOnly: boolean, xScoreNames: string[], yScoreNames: string[], predictionName: string, analysisTableName: string, explVarTableName: string) : Promise<void> {
  await PackageFunctions.multivariateAnalysisTransform(table, featureNames, predict, components, isQuadratic, componentsOnly, xScoreNames, yScoreNames, predictionName, analysisTableName, explVarTableName);
}

//description: Restore the multivariate analysis model viewer.
//input: viewer viewer 
export function mvaModelInitFunction(viewer: any) : void {
  PackageFunctions.mvaModelInitFunction(viewer);
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

//description: Fills missing values in place using the k-nearest neighbors method; cells with no usable neighbors stay empty
//input: dataframe table 
//input: list<string> columns { description: Columns to fill, an empty list fills every feature with missing values }
//input: list<string> features { description: Columns used to find the nearest rows }
//input: int neighbors = 4 { min: 1; description: Nearest rows used to fill a value. }
//input: string distance = 'Euclidean' { choices: ["Euclidean","Manhattan"]; description: How the distance between rows is measured. }
//output: dataframe result
export function knnImpute(table: DG.DataFrame, columns: string[], features: string[], neighbors: number, distance: string) : any {
  return PackageFunctions.knnImpute(table, columns, features, neighbors, distance);
}

//description: Classification quality of predicted vs actual classes from Crux (scikit-learn metrics): accuracy, balanced accuracy, MCC, precision, recall, F1, the confusion matrix, per-class values and, with class probabilities, AUC-ROC, average precision and the ROC and precision-recall curves of the positive class. Used by Forge
//input: column actual { description: The actual classes }
//input: column predicted { description: The predicted classes }
//input: string positiveClass { optional: true; description: The positive class of two classes - its precision, recall, F1, AUC-ROC, average precision and curves. Empty for none }
//input: list<string> classes { optional: true; description: The classes the probability columns stand for, in column order }
//input: dataframe probabilities { optional: true; description: One column per entry of classes - the probability of the class per row }
//input: bool curves = false { description: Also the ROC and precision-recall curves of the positive class }
//output: object report { description: n, labels, accuracy, balancedAccuracy, mcc, optional cohenKappa, precision, recall, auc and averagePrecision, f1, confusion, classes (per class: label, rows, optional precision and recall, f1, optional auc and averagePrecision), optional roc (fpr, tpr, threshold) and pr (precision, recall, threshold) }
export async function classificationReport(actual: DG.Column, predicted: DG.Column, positiveClass: string, classes: string[], probabilities: DG.DataFrame, curves: boolean) : Promise<any> {
  return await PackageFunctions.classificationReport(actual, predicted, positiveClass, classes, probabilities, curves);
}

//input: dataframe df 
//input: column predictColumn 
//input: int iterations = 1000 { caption: Iterations; min: 1; step: 50; max: 10000; description: Largest number of solver iterations (elastic net). }
//input: double alpha = 0 { caption: L1; min: 0; max: 100; description: L1 (Lasso) regularization term. 0 means plain ordinary least squares. }
//input: double lambda = 0 { caption: L2; min: 0; max: 100; description: L2 (Ridge) regularization term. 0 means plain ordinary least squares. }
//output: dynamic model
//meta.mlname: Linear Regression
//meta.mlrole: train
export async function trainLinearRegression(df: DG.DataFrame, predictColumn: DG.Column, iterations: number, alpha: number, lambda: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainLinearRegression(df, predictColumn, iterations, alpha, lambda);
}

//name: Apply Linear Regression
//description: Predict the target for a table using a trained linear regression model.
//input: dataframe df 
//input: dynamic model { description: Trained linear regression model to apply. }
//output: dataframe result
//meta.mlname: Linear Regression
//meta.mlrole: apply
export async function applyLinearRegression(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyLinearRegression(df, model);
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
//input: double iterations = 100 { min: 1; max: 10000; step: 10; description: Largest number of training steps before training stops. }
//input: double penalty = 0.1 { min: 0.0001; max: 1; description: Regularization rate. }
//input: double tolerance = 0.001 { min: 0.00001; max: 0.1; description: Smallest improvement worth continuing training for. }
//output: dynamic model
//meta.mlname: Softmax
//meta.mlrole: train
export async function trainSoftmax(df: DG.DataFrame, predictColumn: DG.Column, iterations: number, penalty: number, tolerance: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainSoftmax(df, predictColumn, iterations, penalty, tolerance);
}

//name: Apply Softmax
//description: Classify the rows of a table using a trained softmax (multinomial logistic regression) model.
//input: dataframe df 
//input: dynamic model { description: Trained softmax classifier model to apply. }
//output: dataframe result
//meta.mlname: Softmax
//meta.mlrole: apply
export async function applySoftmax(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applySoftmax(df, model);
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

//name: Apply PLS Regression
//description: Predict the target for a table using a trained partial least squares (PLS) regression model.
//input: dataframe df 
//input: dynamic model { description: Trained PLS regression model to apply. }
//output: dataframe result
//meta.mlname: PLS Regression
//meta.mlrole: apply
export async function applyPLSRegression(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyPLSRegression(df, model);
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
//meta.mlnativenans: true
export async function trainXGBooster(df: DG.DataFrame, predictColumn: DG.Column, iterations: number, eta: number, maxDepth: number, lambda: number, alpha: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainXGBooster(df, predictColumn, iterations, eta, maxDepth, lambda, alpha);
}

//name: Apply XGBoost
//description: Predict the target for a table using a trained XGBoost gradient-boosting model.
//input: dataframe df 
//input: dynamic model { description: Trained XGBoost model to apply. }
//output: dataframe result
//meta.mlname: XGBoost
//meta.mlrole: apply
export async function applyXGBooster(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyXGBooster(df, model);
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

//input: dataframe df 
//input: column predictColumn 
//input: string kernel = 'RBF' { choices: ["Linear","Polynomial","RBF","Sigmoid"]; description: Kernel function: linear, polynomial, RBF (radial basis function), or sigmoid. }
//input: double cost = 1 { caption: Penalty; min: 0.001; max: 1000; description: Penalty for misclassifying training points. Higher values fit the training data more tightly and risk overfitting, lower values give a smoother, more general model. }
//input: double gamma = 0 { min: 0; max: 100; description: Kernel coefficient. Used by the RBF, polynomial, and sigmoid kernels. Leave 0 to use 1/number-of-features. }
//input: int degree = 3 { min: 1; max: 10; description: Degree of the polynomial kernel. Higher values give a more flexible, curved decision boundary. Used by the polynomial kernel only. }
//input: double coef0 = 0 { caption: Offset; min: -100; max: 100; description: Constant added inside the kernel before the polynomial power or the sigmoid. Used by the polynomial and sigmoid kernels only. }
//input: double epsilon = 0.1 { min: 0; max: 10; description: Width of the error-insensitive tube in SVR regression. Training errors smaller than this are not penalized, so larger values give a smoother model with fewer support vectors. Used for regression targets only. }
//output: dynamic model
//meta.mlname: SVM
//meta.mlrole: train
//meta.mlhasrows: true
export async function trainSVM(df: DG.DataFrame, predictColumn: DG.Column, kernel: string, cost: number, gamma: number, degree: number, coef0: number, epsilon: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainSVM(df, predictColumn, kernel, cost, gamma, degree, coef0, epsilon);
}

//name: Apply SVM
//description: Predict the target for a table using a trained support vector machine (SVM) model.
//input: dataframe df 
//input: dynamic model { description: Trained SVM model to apply. }
//output: dataframe result
//meta.mlname: SVM
//meta.mlrole: apply
export async function applySVM(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applySVM(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: SVM
//meta.mlrole: isApplicable
export function isApplicableSVM(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableSVM(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: SVM
//meta.mlrole: isInteractive
export function isInteractiveSVM(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveSVM(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: int trees = 100 { caption: Trees; min: 1; max: 1000; description: Number of trees in the forest. }
//input: int maxDepth = 10 { caption: Max depth; min: 1; max: 50; description: Largest depth of a tree. Bounds the model size. }
//input: int minSamplesLeaf = 1 { caption: Min rows per leaf; min: 1; max: 1000; description: Smallest number of training rows a leaf may hold. }
//input: string maxFeatures = 'auto' { caption: Features per split; choices: ["auto","sqrt","log2","all"]; description: Features tried at each split: auto is sqrt for classification and all for regression. }
//input: int seed = 0 { caption: Seed; min: 0; max: 2147483647; description: Random seed. The same seed gives the same model. }
//output: dynamic model
//meta.mlname: Random Forest
//meta.mlrole: train
//meta.mlnativenans: true
export async function trainRandomForest(df: DG.DataFrame, predictColumn: DG.Column, trees: number, maxDepth: number, minSamplesLeaf: number, maxFeatures: string, seed: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainRandomForest(df, predictColumn, trees, maxDepth, minSamplesLeaf, maxFeatures, seed);
}

//name: Apply Random Forest
//description: Predict the target for a table using a trained random forest model.
//input: dataframe df 
//input: dynamic model { description: Trained random forest model to apply. }
//output: dataframe result
//meta.mlname: Random Forest
//meta.mlrole: apply
export async function applyRandomForest(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyRandomForest(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Random Forest
//meta.mlrole: isApplicable
export function isApplicableRandomForest(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableRandomForest(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Random Forest
//meta.mlrole: isInteractive
export function isInteractiveRandomForest(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveRandomForest(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: int trees = 100 { caption: Trees; min: 1; max: 1000; description: Number of trees in the forest. }
//input: int maxDepth = 10 { caption: Max depth; min: 1; max: 50; description: Largest depth of a tree. Bounds the model size. }
//input: int minSamplesLeaf = 1 { caption: Min rows per leaf; min: 1; max: 1000; description: Smallest number of training rows a leaf may hold. }
//input: string maxFeatures = 'auto' { caption: Features per split; choices: ["auto","sqrt","log2","all"]; description: Features tried at each split: auto is sqrt for classification and all for regression. }
//input: int seed = 0 { caption: Seed; min: 0; max: 2147483647; description: Random seed. The same seed gives the same model. }
//output: dynamic model
//meta.mlname: Extra Trees
//meta.mlrole: train
//meta.mlnativenans: true
export async function trainExtraTrees(df: DG.DataFrame, predictColumn: DG.Column, trees: number, maxDepth: number, minSamplesLeaf: number, maxFeatures: string, seed: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainExtraTrees(df, predictColumn, trees, maxDepth, minSamplesLeaf, maxFeatures, seed);
}

//name: Apply Extra Trees
//description: Predict the target for a table using a trained extra trees model.
//input: dataframe df 
//input: dynamic model { description: Trained extra trees model to apply. }
//output: dataframe result
//meta.mlname: Extra Trees
//meta.mlrole: apply
export async function applyExtraTrees(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyExtraTrees(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Extra Trees
//meta.mlrole: isApplicable
export function isApplicableExtraTrees(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableExtraTrees(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Extra Trees
//meta.mlrole: isInteractive
export function isInteractiveExtraTrees(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveExtraTrees(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: int iterations = 100 { caption: Iterations; min: 1; max: 1000; description: Boosting iterations (trees per class). }
//input: double learningRate = 0.1 { caption: Learning rate; min: 0.001; max: 1; description: Shrinkage of each tree's contribution. }
//input: int maxLeafNodes = 31 { caption: Max leaves; min: 2; max: 255; description: Largest number of leaves per tree. }
//input: int maxDepth = 0 { caption: Max depth; min: 0; max: 50; description: Largest depth of a tree. 0 means no limit. }
//input: int minSamplesLeaf = 20 { caption: Min rows per leaf; min: 1; max: 1000; description: Smallest number of training rows a leaf may hold. }
//input: double l2 = 0 { caption: L2; min: 0; max: 100; description: L2 penalty of the leaf values. }
//input: int seed = 0 { caption: Seed; min: 0; max: 2147483647; description: Random seed. The same seed gives the same model. }
//output: dynamic model
//meta.mlname: Hist Gradient Boosting
//meta.mlrole: train
//meta.mlnativenans: true
export async function trainHistGradientBoosting(df: DG.DataFrame, predictColumn: DG.Column, iterations: number, learningRate: number, maxLeafNodes: number, maxDepth: number, minSamplesLeaf: number, l2: number, seed: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainHistGradientBoosting(df, predictColumn, iterations, learningRate, maxLeafNodes, maxDepth, minSamplesLeaf, l2, seed);
}

//name: Apply Hist Gradient Boosting
//description: Predict the target for a table using a trained histogram gradient boosting model.
//input: dataframe df 
//input: dynamic model { description: Trained histogram gradient boosting model to apply. }
//output: dataframe result
//meta.mlname: Hist Gradient Boosting
//meta.mlrole: apply
export async function applyHistGradientBoosting(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyHistGradientBoosting(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Hist Gradient Boosting
//meta.mlrole: isApplicable
export function isApplicableHistGradientBoosting(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableHistGradientBoosting(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Hist Gradient Boosting
//meta.mlrole: isInteractive
export function isInteractiveHistGradientBoosting(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveHistGradientBoosting(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: int iterations = 100 { caption: Iterations; min: 1; max: 5000; description: Boosting iterations (trees). }
//input: double learningRate = 0.1 { caption: Learning rate; min: 0.001; max: 1; description: Shrinkage of each tree's contribution. }
//input: int depth = 6 { caption: Depth; min: 1; max: 16; description: Depth of the symmetric trees. }
//input: double l2LeafReg = 3 { caption: L2; min: 0; max: 100; description: L2 regularization of the leaf values. }
//input: int seed = 0 { caption: Seed; min: 0; max: 2147483647; description: Random seed. The same seed gives the same model. }
//output: dynamic model
//meta.mlname: CatBoost
//meta.mlrole: train
//meta.mlnativenans: true
export async function trainCatBoost(df: DG.DataFrame, predictColumn: DG.Column, iterations: number, learningRate: number, depth: number, l2LeafReg: number, seed: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainCatBoost(df, predictColumn, iterations, learningRate, depth, l2LeafReg, seed);
}

//name: Apply CatBoost
//description: Predict the target for a table using a trained CatBoost gradient-boosting model.
//input: dataframe df 
//input: dynamic model { description: Trained CatBoost model to apply. }
//output: dataframe result
//meta.mlname: CatBoost
//meta.mlrole: apply
export async function applyCatBoost(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyCatBoost(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: CatBoost
//meta.mlrole: isApplicable
export function isApplicableCatBoost(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableCatBoost(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: CatBoost
//meta.mlrole: isInteractive
export function isInteractiveCatBoost(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveCatBoost(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: int maxDepth = 10 { caption: Max depth; min: 1; max: 50; description: Largest depth of the tree. Bounds the model size. }
//input: int minSamplesLeaf = 1 { caption: Min rows per leaf; min: 1; max: 1000; description: Smallest number of training rows a leaf may hold. }
//input: int seed = 0 { caption: Seed; min: 0; max: 2147483647; description: Random seed. The same seed gives the same model. }
//output: dynamic model
//meta.mlname: Decision Tree
//meta.mlrole: train
//meta.mlnativenans: true
export async function trainDecisionTree(df: DG.DataFrame, predictColumn: DG.Column, maxDepth: number, minSamplesLeaf: number, seed: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainDecisionTree(df, predictColumn, maxDepth, minSamplesLeaf, seed);
}

//name: Apply Decision Tree
//description: Predict the target for a table using a trained decision tree model.
//input: dataframe df 
//input: dynamic model { description: Trained decision tree model to apply. }
//output: dataframe result
//meta.mlname: Decision Tree
//meta.mlrole: apply
export async function applyDecisionTree(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyDecisionTree(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Decision Tree
//meta.mlrole: isApplicable
export function isApplicableDecisionTree(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableDecisionTree(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Decision Tree
//meta.mlrole: isInteractive
export function isInteractiveDecisionTree(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveDecisionTree(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: int neighbors = 5 { caption: Neighbors; min: 1; max: 100; description: Training rows that vote for a prediction. }
//input: string weights = 'Uniform' { caption: Weights; choices: ["Uniform","Distance"]; description: Uniform: every neighbor counts the same. Distance: closer neighbors count more. }
//output: dynamic model
//meta.mlname: k-Nearest Neighbors
//meta.mlrole: train
//meta.mlhasrows: true
export async function trainKnn(df: DG.DataFrame, predictColumn: DG.Column, neighbors: number, weights: string) : Promise<Uint8Array> {
  return await PackageFunctions.trainKnn(df, predictColumn, neighbors, weights);
}

//name: Apply k-Nearest Neighbors
//description: Predict the target for a table using a trained k-nearest neighbors model.
//input: dataframe df 
//input: dynamic model { description: Trained k-nearest neighbors model to apply. }
//output: dataframe result
//meta.mlname: k-Nearest Neighbors
//meta.mlrole: apply
export async function applyKnn(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyKnn(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: k-Nearest Neighbors
//meta.mlrole: isApplicable
export function isApplicableKnn(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableKnn(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: k-Nearest Neighbors
//meta.mlrole: isInteractive
export function isInteractiveKnn(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveKnn(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: dynamic model
//meta.mlname: Naive Bayes
//meta.mlrole: train
export async function trainNaiveBayes(df: DG.DataFrame, predictColumn: DG.Column) : Promise<Uint8Array> {
  return await PackageFunctions.trainNaiveBayes(df, predictColumn);
}

//name: Apply Naive Bayes
//description: Classify the rows of a table using a trained Gaussian naive Bayes model.
//input: dataframe df 
//input: dynamic model { description: Trained naive Bayes classifier model to apply. }
//output: dataframe result
//meta.mlname: Naive Bayes
//meta.mlrole: apply
export async function applyNaiveBayes(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyNaiveBayes(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Naive Bayes
//meta.mlrole: isApplicable
export function isApplicableNaiveBayes(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableNaiveBayes(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Naive Bayes
//meta.mlrole: isInteractive
export function isInteractiveNaiveBayes(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveNaiveBayes(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: double shrinkage = 0 { caption: Shrinkage; min: 0; max: 1; description: Covariance shrinkage from 0 (none) to 1. Helps when features outnumber rows. }
//output: dynamic model
//meta.mlname: LDA
//meta.mlrole: train
export async function trainLda(df: DG.DataFrame, predictColumn: DG.Column, shrinkage: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainLda(df, predictColumn, shrinkage);
}

//name: Apply LDA
//description: Classify the rows of a table using a trained linear discriminant analysis (LDA) model.
//input: dataframe df 
//input: dynamic model { description: Trained LDA classifier model to apply. }
//output: dataframe result
//meta.mlname: LDA
//meta.mlrole: apply
export async function applyLda(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyLda(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: LDA
//meta.mlrole: isApplicable
export function isApplicableLda(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableLda(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: LDA
//meta.mlrole: isInteractive
export function isInteractiveLda(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveLda(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: string hiddenLayers = 100 { caption: Hidden layers; description: Comma-separated sizes of the hidden layers, e.g. 100 or 64,32. }
//input: string activation = 'ReLU' { caption: Activation; choices: ["ReLU","Tanh","Logistic","Identity"]; description: Activation function of the hidden layers. }
//input: int iterations = 200 { caption: Iterations; min: 1; max: 2000; description: Largest number of training epochs. }
//input: double learningRate = 0.001 { caption: Learning rate; min: 0.00001; max: 1; description: Step size of the weight updates. }
//input: double alpha = 0.0001 { caption: L2; min: 0; max: 10; description: L2 penalty of the weights. }
//input: int seed = 0 { caption: Seed; min: 0; max: 2147483647; description: Random seed. The same seed gives the same model. }
//output: dynamic model
//meta.mlname: Neural Network
//meta.mlrole: train
export async function trainMlp(df: DG.DataFrame, predictColumn: DG.Column, hiddenLayers: string, activation: string, iterations: number, learningRate: number, alpha: number, seed: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainMlp(df, predictColumn, hiddenLayers, activation, iterations, learningRate, alpha, seed);
}

//name: Apply Neural Network
//description: Predict the target for a table using a trained neural network (multilayer perceptron) model.
//input: dataframe df 
//input: dynamic model { description: Trained neural network model to apply. }
//output: dataframe result
//meta.mlname: Neural Network
//meta.mlrole: apply
export async function applyMlp(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyMlp(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Neural Network
//meta.mlrole: isApplicable
export function isApplicableMlp(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableMlp(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: Neural Network
//meta.mlrole: isInteractive
export function isInteractiveMlp(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveMlp(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: string variant = 'TD' { caption: Variant; choices: ["TD","TD-S"]; description: TD: the tuned defaults (numerical embeddings, parametric activations, dropout, weight decay). TD-S: the simplified variant. }
//input: int epochs = 128 { caption: Epochs; min: 1; max: 2000; description: Training epochs. The best epoch on the validation split is kept. }
//input: string hiddenLayers = 256,256,256 { caption: Hidden layers; description: Comma-separated sizes of the hidden layers, e.g. 256,256,256 or 128,128. }
//input: double learningRate = 0 { caption: Learning rate; min: 0; max: 1; description: Base learning rate. 0 keeps the RealMLP default (0.04 for classification, 0.2 for regression). }
//input: int batchSize = 256 { caption: Batch size; min: 1; max: 4096; description: Rows per optimizer step. }
//input: double validationFraction = 0.2 { caption: Validation fraction; min: 0.05; max: 0.5; description: Share of the training rows held out to pick the best epoch. }
//input: int seed = 0 { caption: Seed; min: 0; max: 2147483647; description: Random seed. The same seed gives the same model. }
//output: dynamic model
//meta.mlname: RealMLP
//meta.mlrole: train
export async function trainRealMlp(df: DG.DataFrame, predictColumn: DG.Column, variant: string, epochs: number, hiddenLayers: string, learningRate: number, batchSize: number, validationFraction: number, seed: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainRealMlp(df, predictColumn, variant, epochs, hiddenLayers, learningRate, batchSize, validationFraction, seed);
}

//name: Apply RealMLP
//description: Predict the target for a table using a trained RealMLP (tuned-default tabular MLP) model.
//input: dataframe df 
//input: dynamic model { description: Trained RealMLP model to apply. }
//output: dataframe result
//meta.mlname: RealMLP
//meta.mlrole: apply
export async function applyRealMlp(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyRealMlp(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: RealMLP
//meta.mlrole: isApplicable
export function isApplicableRealMlp(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableRealMlp(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: RealMLP
//meta.mlrole: isInteractive
export function isInteractiveRealMlp(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveRealMlp(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//input: int ensembleSize = 8 { caption: Ensemble size; min: 1; max: 128; description: Weight-sharing MLPs averaged at prediction (the TabM paper uses 32). }
//input: int width = 128 { caption: Width; min: 1; max: 2048; description: Units per hidden layer (the paper uses 512). }
//input: string embeddings = 'Piecewise linear' { caption: Embeddings; choices: ["Piecewise linear","Periodic","Linear-ReLU","None"]; description: Embeddings of the numerical features. Piecewise linear over quantile bins is the paper's best. }
//input: int epochs = 32 { caption: Epochs; min: 1; max: 1000; description: Largest number of training epochs. Early stopping may end training sooner. }
//input: int patience = 16 { caption: Patience; min: 1; max: 100; description: Epochs without a better validation score before training stops. }
//input: double learningRate = 0.002 { caption: Learning rate; min: 0.00001; max: 1; description: AdamW learning rate. }
//input: double dropout = 0.1 { caption: Dropout; min: 0; max: 0.9; description: Dropout rate of the hidden layers. }
//input: double validationFraction = 0.2 { caption: Validation fraction; min: 0.05; max: 0.5; description: Share of the training rows held out to pick the best epoch. }
//input: int seed = 0 { caption: Seed; min: 0; max: 2147483647; description: Random seed. The same seed gives the same model. }
//output: dynamic model
//meta.mlname: TabM
//meta.mlrole: train
//meta.mlnativenans: true
export async function trainTabm(df: DG.DataFrame, predictColumn: DG.Column, ensembleSize: number, width: number, embeddings: string, epochs: number, patience: number, learningRate: number, dropout: number, validationFraction: number, seed: number) : Promise<Uint8Array> {
  return await PackageFunctions.trainTabm(df, predictColumn, ensembleSize, width, embeddings, epochs, patience, learningRate, dropout, validationFraction, seed);
}

//name: Apply TabM
//description: Predict the target for a table using a trained TabM (ensemble tabular MLP) model.
//input: dataframe df 
//input: dynamic model { description: Trained TabM model to apply. }
//output: dataframe result
//meta.mlname: TabM
//meta.mlrole: apply
export async function applyTabm(df: DG.DataFrame, model: any) : Promise<any> {
  return await PackageFunctions.applyTabm(df, model);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: TabM
//meta.mlrole: isApplicable
export function isApplicableTabm(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isApplicableTabm(df, predictColumn);
}

//input: dataframe df 
//input: column predictColumn 
//output: bool result
//meta.mlname: TabM
//meta.mlrole: isInteractive
export function isInteractiveTabm(df: DG.DataFrame, predictColumn: DG.Column) : boolean {
  return PackageFunctions.isInteractiveTabm(df, predictColumn);
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

//name: Generate pMPO Dataset
//description: Generate a synthetic dataset oriented on probabilistic multi-parameter optimization (pMPO) modeling.
//input: int samples { description: Number of rows (samples) to generate. }
//output: dataframe Synthetic
export async function generatePmpoDataset(samples: number) : Promise<any> {
  return await PackageFunctions.generatePmpoDataset(samples);
}
