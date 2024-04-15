/* eslint-disable camelcase */
/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DemoScript} from '@datagrok-libraries/tutorials/src/demo-script';

import {_initEDAAPI} from '../wasm/EDAAPI';
import {computePCA, computePLS} from './eda-tools';
import {addPrefixToEachColumnName, addPLSvisualization, regressionCoefficientsBarChart,
  scoresScatterPlot, predictedVersusReferenceScatterPlot} from './eda-ui';
import {carsDataframe, testDataForBinaryClassification} from './data-generators';
import {LINEAR, RBF, POLYNOMIAL, SIGMOID,
  getTrainedModel, getPrediction, showTrainReport, getPackedModel} from './svm';

import {oneWayAnova} from './anova/anova-tools';
import {runOneWayAnova} from './anova/anova-ui';
import {getDbscanWorker} from '@datagrok-libraries/math';

import {DistanceAggregationMethod, DistanceAggregationMethods} from '@datagrok-libraries/ml/src/distance-matrix/types';
import {MultiColumnDimReductionEditor} from
  '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/multi-column-dim-reduction-editor';
import {multiColReduceDimensionality} from
  '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import {KnownMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';

import {runKNNImputer} from './missing-values-imputation/ui';
import {MCLEditor} from '@datagrok-libraries/ml/src/MCL/mcl-editor';
import {markovCluster} from '@datagrok-libraries/ml/src/MCL/clustering-view';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function init(): Promise<void> {
  await _initEDAAPI();
}

//top-menu: ML | Cluster | DBSCAN...
//name: DBSCAN
//description: Density-based spatial clustering of applications with noise (DBSCAN)
//input: dataframe df
//input: column xCol {type: numerical}
//input: column yCol {type: numerical}
//input: double epsilon = 0.02 {caption: Epsilon} [The maximum distance between two samples for them to be considered as in the same neighborhood.]
//input: int minPts = 4 {caption: Minimum points} [The number of samples (or total weight) in a neighborhood for a point to be considered as a core point.]
//output: column cluster
export async function dbScan(df: DG.DataFrame, xCol: DG.Column, yCol: DG.Column, epsilon: number, minPts: number) {
  const x = xCol.getRawData() as Float32Array;
  const y = yCol.getRawData() as Float32Array;
  const res = await getDbscanWorker(x, y, epsilon, minPts);
  const clusterColName = df.columns.getUnusedName('Cluster (DBSCAN)');
  const cluster = DG.Column.fromInt32Array(clusterColName, res);
  df.columns.add(cluster);
  return cluster;
}

//top-menu: ML | Analyze | PCA...
//name: PCA
//description: Principal component analysis (PCA)
//input: dataframe table
//input: column_list features {type: numerical}
//input: int components = 2 {caption: Components} [Number of components.]
//input: bool center = false [Indicating whether the variables should be shifted to be zero centered.]
//input: bool scale = false [Indicating whether the variables should be scaled to have unit variance.]
export async function PCA(table: DG.DataFrame, features: DG.ColumnList, components: number, center: boolean, scale: boolean): Promise<void> {
  const pcaTable = await computePCA(table, features, components, center, scale);
  addPrefixToEachColumnName('PCA', pcaTable.columns);

  if (table.id === null) // table is loaded from a local file
    grok.shell.addTableView(pcaTable);
  else {
    const cols = table.columns;

    for (const col of pcaTable.columns) {
      col.name = cols.getUnusedName(col.name);
      cols.add(col);
    }
  }
}

//name: DBSCAN clustering
//tags: dim-red-postprocessing-function
//meta.defaultPostProcessingFunction: true
//input: column col1
//input: column col2
//input: double epsilon = 0.01 {default: 0.01}[Minimum distance between two points to be considered as in the same neighborhood.]
//input: int minimumPoints = 5 {default: 5}[Minimum number of points to form a dense region.]
export async function dbscanPostProcessingFunction(col1: DG.Column, col2: DG.Column, epsilon: number, minimumPoints: number) {
  const df = col1.dataFrame;
  if (df === null)
    return;
  const resCol = await dbScan(df, col1, col2, epsilon, minimumPoints);
  df.changeColumnType(resCol, 'string');
  const colNames = [col1.name, col2.name];
  const tv = grok.shell.tableView(df.name);
  if (!tv)
    return;
  // find the correct scatterPlotViewer and set the colorColumnName
  for (const v of tv.viewers) {
    if (v instanceof DG.ScatterPlotViewer && colNames.includes(v.props.xColumnName) && colNames.includes(v.props.yColumnName)) {
      v.props.colorColumnName = resCol.name;
      return;
    }
  }
}

//name: None (number)
//tags: dim-red-preprocessing-function
//meta.supportedTypes: int,float,double,qnum
//meta.supportedDistanceFunctions: Difference
//input: column col
//input: string _metric {optional: true}
//output: object result
export function numberPreprocessingFunction(col: DG.Column, _metric: string) {
  const range = col.stats.max - col.stats.min;
  const entries = col.toList();
  return {entries, options: {range}};
}

//name: None (string)
//tags: dim-red-preprocessing-function
//meta.supportedTypes: string
//meta.supportedDistanceFunctions: Levenshtein,Hamming,One-Hot
//input: column col
//input: string _metric {optional: true}
//output: object result
export function stringPreprocessingFunction(col: DG.Column, _metric: string) {
  const entries = col.toList();
  return {entries, options: {}};
}

//top-menu: ML | Reduce Dimensionality...
//name: Multi Column Dimensionality Reduction
export async function reduceDimensionality(): Promise<void> {
  const editor = new MultiColumnDimReductionEditor();
  ui.dialog('Dimensionality reduction').add(editor.getEditor()).onOK(async () => {
    const params = editor.getParams();
    if (params.columns.length === 0)
      return;
    await multiColReduceDimensionality(params.table, params.columns, params.methodName as DimReductionMethods,
      params.distanceMetrics as KnownMetrics[],
      params.weights, params.preprocessingFunctions, params.aggreaggregationMethod as DistanceAggregationMethods,
      !!params.plotEmbeddings, !!params.clusterEmbeddings, params.options, {
        fastRowCount: 10000,
      }, params.postProcessingFunction, params.postProcessingFunctionArgs);
  }).show();
}

//name: GetMCLEditor
//tags: editor
//input: funccall call
export function GetMCLEditor(call: DG.FuncCall): void {
  try {
    const funcEditor = new MCLEditor();
    ui.dialog('Markov clustering')
      .add(funcEditor.getEditor())
      .onOK(async () => {
        const params = funcEditor.params;
        return call.func.prepare({
          df: params.table, cols: params.columns, metrics: params.distanceMetrics,
          weights: params.weights, aggregationMethod: params.aggreaggregationMethod, preprocessingFuncs: params.preprocessingFunctions,
          preprocessingFuncArgs: params.preprocessingFuncArgs, threshold: params.threshold, maxIterations: params.maxIterations,
        }).call(true);
      }).show();
  } catch (err: any) {
    const errMsg = err instanceof Error ? err.message : err.toString();
    const errStack = err instanceof Error ? err.stack : undefined;
    grok.shell.error(`Get region editor error: ${errMsg}`);
    _package.logger.error(errMsg, undefined, errStack);
  }
}


//top-menu: ML | Cluster | MCL...
//name: MCL
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
//editor: EDA: GetMCLEditor
export async function MCL(df: DG.DataFrame, cols: DG.Column[], metrics: KnownMetrics[],
  weights: number[], aggregationMethod: DistanceAggregationMethod, preprocessingFuncs: (DG.Func | null | undefined)[],
  preprocessingFuncArgs: any[], threshold: number = 80, maxIterations: number = 10) {
  const res = (await markovCluster(df, cols, metrics, weights,
    aggregationMethod, preprocessingFuncs, preprocessingFuncArgs, threshold, maxIterations));
  return res?.sc;
}

//top-menu: ML | Analyze | Multivariate Analysis...
//name: Multivariate Analysis (PLS)
//description: Multidimensional data analysis using partial least squares (PLS) regression. It reduces the predictors to a smaller set of uncorrelated components and performs least squares regression on them.
//input: dataframe table
//input: column names
//input: column_list features {type: numerical}
//input: column predict {type: numerical}
//input: int components = 3
export async function PLS(table: DG.DataFrame, names: DG.Column, features: DG.ColumnList,
  predict: DG.Column, components: number): Promise<void> {
  const plsResults = await computePLS(table, features, predict, components);
  addPLSvisualization(table, names, features, predict, plsResults);
}

//name: MVA demo
//description: Multidimensional data analysis using partial least squares (PLS) regression. It reduces the predictors to a smaller set of uncorrelated components and performs least squares regression on them.
//meta.demoPath: Compute | Multivariate analysis
//meta.isDemoScript: True
export async function demoMultivariateAnalysis(): Promise<any> {
  const demoScript = new DemoScript('Partial least squares regression',
    'Analysis of multidimensional data.');

  const cars = carsDataframe();

  const components = 3;
  const names = cars.columns.byName('model');
  const predict = cars.columns.byName('price');
  const features = cars.columns.remove('price').remove('model');
  const plsOutput = await computePLS(cars, features, predict, components);

  const sourceCars = carsDataframe();
  sourceCars.name = 'Cars';
  let view: any;
  let dialog: any;

  await demoScript
    .step('Data', async () => {
      grok.shell.addTableView(sourceCars);
      view = grok.shell.getTableView(sourceCars.name);
    }, {description: 'Each car has many features - patterns extraction is complicated.', delay: 0})
    .step('Model', async () => {
      dialog = ui.dialog({title: 'Multivariate Analysis (PLS)'})
        .add(ui.tableInput('Table', sourceCars))
        .add(ui.columnsInput('Features', cars, features.toList, {available: undefined, checked: features.names()}))
        .add(ui.columnInput('Names', cars, names, undefined))
        .add(ui.columnInput('Predict', cars, predict, undefined))
        .add(ui.intInput('Components', components, undefined))
        .onOK(() => {
          grok.shell.info('Multivariate analysis has been already performed.');
        })
        .show({x: 400, y: 140});
    }, {description: 'Predict car price by its other features.', delay: 0})
    .step('Regression coeffcicients', async () => {
      dialog.close();
      view.addViewer(regressionCoefficientsBarChart(features, plsOutput[1]));
    },
    {description: 'The feature "diesel" affects the price the most.', delay: 0})
    .step('Scores', async () => {view.addViewer(scoresScatterPlot(names, plsOutput[2], plsOutput[3]));},
      {description: 'Similarities & dissimilarities: alfaromeo and mercedes are different.', delay: 0})
    .step('Prediction', async () => {view.addViewer(predictedVersusReferenceScatterPlot(names, predict, plsOutput[0]));},
      {description: 'Closer to the line means better price prediction.', delay: 0})
    .start();
}

//name: Generate linear separable dataset
//description: Generates linear separble dataset for testing binary classificators
//input: string name = 'Data' {caption: name; category: Dataset}
//input: int samplesCount = 1000 {caption: samples; category: Size}
//input: int featuresCount = 2 {caption: features; category: Size}
//input: double min = -39 {caption: min; category: Range}
//input: double max = 173 {caption: max; category: Range}
//input: double violatorsPercentage = 5 {caption: violators; units: %; category: Dataset}
//output: dataframe df
export async function testDataLinearSeparable(name: string, samplesCount: number, featuresCount: number,
  min: number, max: number, violatorsPercentage: number): Promise<DG.DataFrame> {
  return await testDataForBinaryClassification(LINEAR, [0, 0], name, samplesCount, featuresCount,
    min, max, violatorsPercentage);
}

//name: Generate linear non-separable dataset
//description: Generates linear non-separble dataset for testing binary classificators
//input: string name = 'Data' {caption: name; category: Dataset}
//input: double sigma = 90  {caption: sigma; category: Hyperparameters} [RBF-kernel paramater]
//input: int samplesCount = 1000 {caption: samples; category: Size}
//input: int featuresCount = 2 {caption: features; category: Size}
//input: double min = -39 {caption: min; category: Range}
//input: double max = 173 {caption: max; category: Range}
//input: double violatorsPercentage = 5 {caption: violators; units: %; category: Dataset}
//output: dataframe df
export async function testDataLinearNonSeparable(name: string, sigma: number, samplesCount: number,
  featuresCount: number, min: number, max: number, violatorsPercentage: number): Promise<DG.DataFrame> {
  return await testDataForBinaryClassification(RBF, [sigma, 0], name, samplesCount, featuresCount,
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
export async function trainLinearKernelSVM(df: DG.DataFrame, predict_column: string,
  gamma: number, toShowReport: boolean): Promise<any> {
  const trainedModel = await getTrainedModel({gamma: gamma, kernel: LINEAR}, df, predict_column);

  if (toShowReport)
    showTrainReport(df, trainedModel);

  return getPackedModel(trainedModel);
}

//name: applyLinearKernelSVM
//meta.mlname: linear kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export async function applyLinearKernelSVM(df: DG.DataFrame, model: any): Promise<DG.DataFrame> {
  return await getPrediction(df, model);
}

//name: trainRBFkernelSVM
//meta.mlname: RBF-kernel LS-SVM
//meta.mlrole: train
//input: dataframe df
//input: string predict_column
//input: double gamma = 1.0 {category: Hyperparameters}
//input: double sigma = 1.5 {category: Hyperparameters}
//input: bool toShowReport = false {caption: to show report; category: Report}
//output: dynamic model
export async function trainRBFkernelSVM(df: DG.DataFrame, predict_column: string,
  gamma: number, sigma: number, toShowReport: boolean): Promise<any> {
  const trainedModel = await getTrainedModel(
    {gamma: gamma, kernel: RBF, sigma: sigma},
    df, predict_column);

  if (toShowReport)
    showTrainReport(df, trainedModel);

  return getPackedModel(trainedModel);
}

//name: applyRBFkernelSVM
//meta.mlname: RBF-kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export async function applyRBFkernelSVM(df: DG.DataFrame, model: any): Promise<DG.DataFrame> {
  return await getPrediction(df, model);
}

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
export async function trainPolynomialKernelSVM(df: DG.DataFrame, predict_column: string,
  gamma: number, c: number, d: number, toShowReport: boolean): Promise<any> {
  const trainedModel = await getTrainedModel(
    {gamma: gamma, kernel: POLYNOMIAL, cParam: c, dParam: d},
    df, predict_column);

  if (toShowReport)
    showTrainReport(df, trainedModel);

  return getPackedModel(trainedModel);
} // trainPolynomialKernelSVM

//name: applyPolynomialKernelSVM
//meta.mlname: polynomial kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export async function applyPolynomialKernelSVM(df: DG.DataFrame, model: any): Promise<DG.DataFrame> {
  return await getPrediction(df, model);
}

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
export async function trainSigmoidKernelSVM(df: DG.DataFrame, predict_column: string,
  gamma: number, kappa: number, theta: number, toShowReport: boolean): Promise<any> {
  const trainedModel = await getTrainedModel(
    {gamma: gamma, kernel: SIGMOID, kappa: kappa, theta: theta},
    df, predict_column);

  if (toShowReport)
    showTrainReport(df, trainedModel);

  return getPackedModel(trainedModel);
} // trainSigmoidKernelSVM

//name: applySigmoidKernelSVM
//meta.mlname: sigmoid kernel LS-SVM
//meta.mlrole: apply
//input: dataframe df
//input: dynamic model
//output: dataframe table
export async function applySigmoidKernelSVM(df: DG.DataFrame, model: any): Promise<DG.DataFrame> {
  return await getPrediction(df, model);
}

//top-menu: ML | Analyze | ANOVA old...
//name: ANOVAold
//description: One-way analysis of variances (ANOVA) determines whether the examined factor has a significant impact on the studied feature.
//input: dataframe table
//input: column factor {type: categorical}
//input: column feature {type: numerical}
//input: double significance = 0.05 [The significance level is a value from the interval (0, 1) specifying the criterion used for rejecting the null hypothesis.]
//input: bool validate = false [Indicates whether the normality of distribution and an eqaulity of varainces should be checked.]
export function anovaOld(table: DG.DataFrame, factor: DG.Column, feature: DG.Column, significance: number, validate: boolean) {
  const res = oneWayAnova(factor, feature, significance, validate);
  //addOneWayAnovaVizualization(table, factor, feature, res);
}

//top-menu: ML | Analyze | ANOVA...
//name: ANOVA
//description: One-way analysis of variances (ANOVA) determines whether the examined factor has a significant impact on the explored feature.
export function anova(): void {
  runOneWayAnova();
}

//top-menu: ML | Missing Values Imputation ...
//name: KNN impute
//desription: Missing values imputation using the k-nearest neighbors method
export function kNNImputation() {
  runKNNImputer();
}
