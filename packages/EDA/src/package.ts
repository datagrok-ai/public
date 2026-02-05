/* eslint-disable camelcase */
/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {_initEDAAPI} from '../wasm/EDAAPI';
import {computePCA} from './eda-tools';
import {addPrefixToEachColumnName} from './eda-ui';
import {LINEAR, RBF, POLYNOMIAL, SIGMOID,
  getTrainedModel, getPrediction, isApplicableSVM, isInteractiveSVM, showTrainReport, getPackedModel} from './svm';

import {PLS_ANALYSIS} from './pls/pls-constants';
import {runMVA, runDemoMVA, getPlsAnalysis, PlsOutput} from './pls/pls-tools';
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
import {MCLViewer} from '@datagrok-libraries/ml/src/MCL/mcl-viewer';
import {MCLSerializableOptions} from '@datagrok-libraries/ml/src/MCL';

import {getLinearRegressionParams, getPredictionByLinearRegression} from './regression';
import {PlsModel} from './pls/pls-ml';
import {SoftmaxClassifier} from './softmax-classifier';

import {initXgboost} from '../wasm/xgbooster';
import {XGBooster} from './xgbooster';

import {ParetoOptimizer} from './pareto-optimization/pareto-optimizer';
import {ParetoFrontViewer} from './pareto-optimization/pareto-front-viewer';

import {Pmpo} from './probabilistic-scoring/prob-scoring';
import {getSynteticPmpoData} from './probabilistic-scoring/data-generator';

export const _package = new DG.Package();
export * from './package.g';

export class PackageFunctions {
  @grok.decorators.func({
    'name': 'info',
  })
  static info() {
    grok.shell.info(_package.webRoot);
  }


  @grok.decorators.init({})
  static async init(): Promise<void> {
    await _initEDAAPI();
    await initXgboost();
  }


  @grok.decorators.func({
    'top-menu': 'ML | Cluster | DBSCAN...',
    'name': 'DBSCAN',
    'description': 'Density-based spatial clustering of applications with noise (DBSCAN)',
  })
  static async dbScan(
    df: DG.DataFrame,
    @grok.decorators.param({'type': 'column', 'options': {'type': 'numerical'}}) xCol: DG.Column,
    @grok.decorators.param({'type': 'column', 'options': {'type': 'numerical'}}) yCol: DG.Column,
    @grok.decorators.param({'options': {'caption': 'Epsilon', 'initialValue': '0.02', 'description': 'The maximum distance between two samples for them to be considered as in the same neighborhood.'}}) epsilon: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Minimum points', 'initialValue': '4', 'description': 'The number of samples (or total weight) in a neighborhood for a point to be considered as a core point.'}}) minPts: number) : Promise<DG.Column> {
    const x = xCol.getRawData() as Float32Array;
    const y = yCol.getRawData() as Float32Array;
    const res = await getDbscanWorker(x, y, epsilon, minPts);
    const clusterColName = df.columns.getUnusedName('Cluster (DBSCAN)');
    const cluster = DG.Column.fromInt32Array(clusterColName, res);
    df.columns.add(cluster);
    return cluster;
  }

  @grok.decorators.func({
    'top-menu': 'ML | Analyze | PCA...',
    'description': 'Principal component analysis (PCA)',
    'helpUrl': '/help/explore/dim-reduction#pca',
  })
  static async PCA(
    @grok.decorators.param({'type': 'dataframe', 'options': {'caption': 'Table'}}) table: DG.DataFrame,
    @grok.decorators.param({'type': 'column_list', 'options': {'type': 'numerical', 'nullable': false}}) features: DG.ColumnList,
    //@ts-ignore
    @grok.decorators.param({'type': 'int', 'options': {'showPlusMinus': true, 'caption': 'Components', 'nullable': false, 'min': '1', 'initialValue': '2', 'description': 'Number of components.'}}) components: number,
    @grok.decorators.param({'type': 'bool', 'options': {'caption': 'Center', 'initialValue': 'false', 'description': 'Indicating whether the variables should be shifted to be zero centered.'}}) center: boolean,
    @grok.decorators.param({'type': 'bool', 'options': {'caption': 'Scale', 'initialValue': 'false', 'description': 'Indicating whether the variables should be scaled to have unit variance.'}}) scale: boolean): Promise<void> {
    try {
      const pcaTable = await computePCA(table, features, components, center, scale);
      addPrefixToEachColumnName('PC', pcaTable.columns);

      if (table.id === null) // table is loaded from a local file
        grok.shell.addTableView(pcaTable);
      else {
        const cols = table.columns;
        const pcaTableCols = pcaTable.columns.toList();

        for (const col of pcaTableCols) {
          pcaTable.columns.remove(col.name);
          col.name = cols.getUnusedName(col.name);
          cols.add(col);
        }
      }
    } catch (error) {
      grok.shell.warning(`Failed to compute PCA: ${error instanceof Error ? error.message : 'platform issue'}`);
    }
  }


  @grok.decorators.func({
    'meta': {'defaultPostProcessingFunction': 'true', 'role': 'dimRedPostprocessingFunction'},
    'name': 'DBSCAN clustering',
  })
  static async dbscanPostProcessingFunction(
    col1: DG.Column,
    col2: DG.Column,
    @grok.decorators.param({options: {initialValue: '0.01', description: 'Minimum distance between two points to be considered as in the same neighborhood.'}}) epsilon: number,
    @grok.decorators.param({type: 'int', options: {initialValue: '5', description: 'Minimum number of points to form a dense region.'}}) minimumPoints: number) {
    const df = col1.dataFrame;
    if (df === null)
      return;
    const resCol = await PackageFunctions.dbScan(df, col1, col2, epsilon, minimumPoints);
    df.changeColumnType(resCol, 'string');
    const colNames = [col1.name, col2.name];
    const tv = grok.shell.tableView(df.name);
    if (!tv)
      return;
    for (const v of tv.viewers) {
      if (v instanceof DG.ScatterPlotViewer && colNames.includes(v.props.xColumnName) && colNames.includes(v.props.yColumnName)) {
        v.props.colorColumnName = resCol.name;
        return;
      }
    }
  }


  @grok.decorators.func({
    'meta': {
      'supportedTypes': 'int,float,double,qnum',
      'supportedDistanceFunctions': 'Difference',
      'role': 'dimRedPreprocessingFunction',
    },
    'name': 'None (number)',
    'outputs': [{name: 'result', type: 'object'}],
  })
  static numberPreprocessingFunction(
    col: DG.Column,
    @grok.decorators.param({'options': {'optional': true}}) _metric: string) {
    const range = col.stats.max - col.stats.min;
    const entries = col.toList();
    return {entries, options: {range}};
  }


  @grok.decorators.func({
    'meta': {
      'supportedTypes': 'string',
      'supportedDistanceFunctions': 'One-Hot,Levenshtein,Hamming',
      'role': 'dimRedPreprocessingFunction',
    },
    'name': 'None (string)',
    'outputs': [{name: 'result', type: 'object'}],
  })
  static stringPreprocessingFunction(
    col: DG.Column,
    @grok.decorators.param({'options': {'optional': true}}) _metric: string) {
    const entries = col.toList();
    return {entries, options: {}};
  }


  @grok.decorators.func({
    'top-menu': 'ML | Reduce Dimensionality...',
    'name': 'Multi Column Dimensionality Reduction',
  })
  static async reduceDimensionality(): Promise<void> {
    const editor = new MultiColumnDimReductionEditor();
    const dialog = ui.dialog('Reduce Dimensionality')
      .add(editor.getEditor())
      .onOK(async () => {
        const params = editor.getParams();
        if (params.columns.length === 0)
          return;
        await multiColReduceDimensionality(params.table, params.columns, params.methodName as DimReductionMethods,
        params.distanceMetrics as KnownMetrics[],
        params.weights, params.preprocessingFunctions, params.aggreaggregationMethod as DistanceAggregationMethods,
        !!params.plotEmbeddings, !!params.clusterEmbeddings, params.options, {
          fastRowCount: 10000,
        }, params.postProcessingFunction, params.postProcessingFunctionArgs, params.vectorDistanceMetric);
      }).show();
    dialog.helpUrl = 'https://datagrok.ai/help/explore/dim-reduction.md';
    const validate = () => {
      const cols = editor.columnsInput.value;
      const okButton = dialog.getButton('OK');
      if (!okButton)
        return;
      const isDisabled = !cols || cols.length === 0;
      if (isDisabled)
        okButton.classList.add('disabled');
      else
        okButton.classList.remove('disabled');
    };
    dialog.history(() => ({editorSettings: editor.getStringInput()}), (x: any) => editor.applyStringInput(x['editorSettings']));
    editor.onColumnsChanged.subscribe(() => {
      try {
        validate();
      } catch (e) {
        console.error(e);
      }
    });
    validate();
  }


  @grok.decorators.editor()
  static GetMCLEditor(
    call: DG.FuncCall): void {
    try {
      const funcEditor = new MCLEditor();
      const dialog = ui.dialog('Markov clustering')
        .add(funcEditor.getEditor())
        .onOK(async () => {
          const params = funcEditor.params;
          return call.func.prepare({
            df: params.table, cols: params.columns, metrics: params.distanceMetrics,
            weights: params.weights, aggregationMethod: params.aggreaggregationMethod, preprocessingFuncs: params.preprocessingFunctions,
            preprocessingFuncArgs: params.preprocessingFuncArgs, threshold: params.threshold, maxIterations: params.maxIterations,
            useWebGPU: params.useWebGPU, inflate: params.inflateFactor, minClusterSize: params.minClusterSize,
          }).call(true);
        }).show();
      dialog.history(() => ({editorSettings: funcEditor.getStringInput()}), (x: any) => funcEditor.applyStringInput(x['editorSettings']));
    } catch (err: any) {
      const errMsg = err instanceof Error ? err.message : err.toString();
      const errStack = err instanceof Error ? err.stack : undefined;
      grok.shell.error(`Get region editor error: ${errMsg}`);
      _package.logger.error(errMsg, undefined, errStack);
    }
  }


  @grok.decorators.func({
    'top-menu': 'ML | Cluster | MCL...',
    'name': 'MCLClustering',
    'description': 'Markov clustering (MCL) is an unsupervised clustering algorithm for graphs based on simulation of stochastic flow.',
    'editor': 'EDA:GetMCLEditor',
    'outputs': [],
  })
  static async MCLClustering(
    df: DG.DataFrame,
    cols: DG.Column[],
    @grok.decorators.param({'type': 'list<string>'}) metrics: KnownMetrics[],
    weights: number[],
    @grok.decorators.param({'type': 'string'}) aggregationMethod: DistanceAggregationMethod,
    @grok.decorators.param({'type': 'list<func>'}) preprocessingFuncs: any[],
    @grok.decorators.param({'type': 'object'}) preprocessingFuncArgs: any[],
    @grok.decorators.param({'type': 'int', 'options': {'initialValue': '80'}}) threshold: number = 80,
    @grok.decorators.param({'type': 'int', 'options': {'initialValue': '10'}}) maxIterations: number = 10,
    @grok.decorators.param({'type': 'bool', 'options': {'initialValue': 'false'}}) useWebGPU: boolean = false,
    @grok.decorators.param({'type': 'double', 'options': {'initialValue': '2'}}) inflate: number = 0,
    @grok.decorators.param({'type': 'int', 'options': {'initialValue': '5'}}) minClusterSize: number = 5): Promise<MCLViewer> {
    const tv = grok.shell.tableView(df.name) ?? grok.shell.addTableView(df);
    const serializedOptions: string = JSON.stringify({
      cols: cols.map((col) => col.name),
      metrics: metrics,
      weights: weights,
      aggregationMethod: aggregationMethod,
      preprocessingFuncs: preprocessingFuncs.map((func) => func?.name ?? null),
      preprocessingFuncArgs: preprocessingFuncArgs,
      threshold: threshold,
      maxIterations: maxIterations,
      useWebGPU: useWebGPU,
      inflate: inflate,
      minClusterSize: minClusterSize ?? 5,
    } satisfies MCLSerializableOptions);

    const viewer = tv.addViewer('MCL', {mclProps: serializedOptions}) as MCLViewer;
    return viewer;
  }

  @grok.decorators.func({
    'outputs': [{'name': 'result', 'type': 'viewer'}],
    'meta': {showInGallery: 'false', role: 'viewer'},
    'name': 'MCL',
    'description': 'Markov clustering viewer',
  })
  static markovClusteringViewer(): MCLViewer {
    return new MCLViewer();
  }

  @grok.decorators.func({
    'outputs': [{'name': 'plsResults', 'type': 'object'}],
    'description': 'Compute partial least squares (PLS) regression analysis components: prediction, regression coefficients, T- & U-scores, X-loadings.',
  })
  static async PLS(
    table: DG.DataFrame,
    @grok.decorators.param({'type': 'column_list', 'options': {'type': 'numerical'}}) features: DG.ColumnList,
    @grok.decorators.param({'type': 'column', 'options': {'type': 'numerical'}}) predict: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'initialValue': '3'}}) components: number,
    @grok.decorators.param({'type': 'column', 'options': {'type': 'string'}}) names: DG.Column): Promise<PlsOutput> {
    return await getPlsAnalysis({
      table: table,
      features: features,
      predict: predict,
      components: components,
      isQuadratic: false,
      names: names,
    });
  }


  @grok.decorators.func({
    'top-menu': 'ML | Analyze | PLS...',
    'description': 'Compute partial least squares (PLS) regression components. They maximally summarize the variation of the predictors while maximizing correlation with the response variable.',
  })
  static async topMenuPLS(): Promise<void> {
    await runMVA(PLS_ANALYSIS.COMPUTE_COMPONENTS);
  }


  @grok.decorators.func({
    'top-menu': 'ML | Analyze | Multivariate Analysis...',
    'name': 'multivariateAnalysis',
    'description': 'Multidimensional data analysis using partial least squares (PLS) regression.',
  })
  static async MVA(): Promise<void> {
    await runMVA(PLS_ANALYSIS.PERFORM_MVA);
  }


  @grok.decorators.func({
    'meta': {
      'demoPath': 'Compute | Multivariate Analysis',
    },
    'name': 'MVA demo',
    'description': 'Multidimensional data analysis using partial least squares (PLS) regression. It identifies latent factors and constructs a linear model based on them.',
  })
  static async demoMultivariateAnalysis(): Promise<void> {
    await runDemoMVA();
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'linear kernel LS-SVM',
      'mlrole': 'train',
    },
  })
  static async trainLinearKernelSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '1.0'}}) gamma: number): Promise<any> {
    const trainedModel = await getTrainedModel({gamma: gamma, kernel: LINEAR}, df, predictColumn);
    return getPackedModel(trainedModel);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'linear kernel LS-SVM',
      'mlrole': 'apply',
    },
  })
  static async applyLinearKernelSVM(
    df: DG.DataFrame,
    model: any): Promise<DG.DataFrame> {
    return await getPrediction(df, model);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'linear kernel LS-SVM',
      'mlrole': 'isApplicable',
    },
  })
  static async isApplicableLinearKernelSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column): Promise<boolean> {
    return isApplicableSVM(df, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'linear kernel LS-SVM',
      'mlrole': 'isInteractive',
    },
  })
  static async isInteractiveLinearKernelSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column): Promise<boolean> {
    return isInteractiveSVM(df, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'linear kernel LS-SVM',
      'mlrole': 'visualize',
    },
  })
  static async visualizeLinearKernelSVM(
    df: DG.DataFrame,
    targetColumn: DG.Column,
    predictColumn: DG.Column,
    model: any): Promise<any> {
    return showTrainReport(df, model);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'RBF-kernel LS-SVM',
      'mlrole': 'train',
    },
  })
  static async trainRBFkernelSVM(df: DG.DataFrame, predictColumn: DG.Column,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '1.0'}}) gamma: number,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '1.5'}}) sigma: number): Promise<any> {
    const trainedModel = await getTrainedModel(
      {gamma: gamma, kernel: RBF, sigma: sigma},
      df, predictColumn);

    return getPackedModel(trainedModel);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'RBF-kernel LS-SVM',
      'mlrole': 'apply',
    },
  })
  static async applyRBFkernelSVM(
    df: DG.DataFrame,
    model: any): Promise<DG.DataFrame> {
    return await getPrediction(df, model);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'RBF-kernel LS-SVM',
      'mlrole': 'isApplicable',
    },
  })
  static async isApplicableRBFkernelSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column): Promise<boolean> {
    return isApplicableSVM(df, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'RBF-kernel LS-SVM',
      'mlrole': 'isInteractive',
    },
  })
  static async isInteractiveRBFkernelSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column): Promise<boolean> {
    return isInteractiveSVM(df, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'RBF-kernel LS-SVM',
      'mlrole': 'visualize',
    },
  })
  static async visualizeRBFkernelSVM(
    df: DG.DataFrame,
    targetColumn: DG.Column,
    predictColumn: DG.Column,
    model: any): Promise<any> {
    return showTrainReport(df, model);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'polynomial kernel LS-SVM',
      'mlrole': 'train',
    },
  })
  static async trainPolynomialKernelSVM(df: DG.DataFrame, predictColumn: DG.Column,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '1.0'}}) gamma: number,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '1'}}) c: number,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '2'}}) d: number): Promise<any> {
    const trainedModel = await getTrainedModel(
      {gamma: gamma, kernel: POLYNOMIAL, cParam: c, dParam: d},
      df, predictColumn);

    return getPackedModel(trainedModel);
  } // trainPolynomialKernelSVM


  @grok.decorators.func({
    'meta': {
      'mlname': 'polynomial kernel LS-SVM',
      'mlrole': 'apply',
    },
  })
  static async applyPolynomialKernelSVM(
    df: DG.DataFrame,
    model: any): Promise<DG.DataFrame> {
    return await getPrediction(df, model);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'polynomial kernel LS-SVM',
      'mlrole': 'isApplicable',
    },
  })
  static async isApplicablePolynomialKernelSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column): Promise<boolean> {
    return isApplicableSVM(df, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'polynomial kernel LS-SVM',
      'mlrole': 'isInteractive',
    },
  })
  static async isInteractivePolynomialKernelSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column): Promise<boolean> {
    return isInteractiveSVM(df, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'polynomial kernel LS-SVM',
      'mlrole': 'visualize',
    },
    'outputs': [{'name': 'widget', 'type': 'dynamic'}],
    'name': 'visualizePolynomialKernelSVM',
  })
  static async visualizePolynomialKernelSVM(
    df: DG.DataFrame,
    targetColumn: DG.Column,
    predictColumn: DG.Column,
    model: any): Promise<any> {
    return showTrainReport(df, model);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'sigmoid kernel LS-SVM',
      'mlrole': 'train',
    },
    'name': 'trainSigmoidKernelSVM',
  })
  static async trainSigmoidKernelSVM(df: DG.DataFrame, predictColumn: DG.Column,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '1.0'}}) gamma: number,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '1'}}) kappa: number,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '1'}}) theta: number): Promise<any> {
    const trainedModel = await getTrainedModel(
      {gamma: gamma, kernel: SIGMOID, kappa: kappa, theta: theta},
      df, predictColumn);

    return getPackedModel(trainedModel);
  } // trainSigmoidKernelSVM


  @grok.decorators.func({
    'meta': {
      'mlname': 'sigmoid kernel LS-SVM',
      'mlrole': 'apply',
    },
    'name': 'applySigmoidKernelSVM',
  })
  static async applySigmoidKernelSVM(
    df: DG.DataFrame,
    model: any): Promise<DG.DataFrame> {
    return await getPrediction(df, model);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'sigmoid kernel LS-SVM',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableSigmoidKernelSVM',
  })
  static async isApplicableSigmoidKernelSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column): Promise<boolean> {
    return isApplicableSVM(df, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'sigmoid kernel LS-SVM',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveSigmoidKernelSVM',
  })
  static async isInteractiveSigmoidKernelSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column): Promise<boolean> {
    return isInteractiveSVM(df, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'sigmoid kernel LS-SVM',
      'mlrole': 'visualize',
    },
    'name': 'visualizeSigmoidKernelSVM',
  })
  static async visualizeSigmoidKernelSVM(
    df: DG.DataFrame,
    targetColumn: DG.Column,
    predictColumn: DG.Column,
    model: any): Promise<any> {
    return showTrainReport(df, model);
  }


  @grok.decorators.func({
    'top-menu': 'ML | Analyze | ANOVA...',
    'name': 'ANOVA',
    'description': 'One-way analysis of variances (ANOVA) determines whether the examined factor has a significant impact on the explored feature.',
  })
  static anova(): void {
    runOneWayAnova();
  }


  @grok.decorators.func({
    'top-menu': 'ML | Impute Missing Values...',
    'name': 'KNN impute',
    'description': 'Missing values imputation using the k-nearest neighbors method (KNN)',
  })
  static kNNImputation() {
    runKNNImputer();
  }


  @grok.decorators.func({
    'name': 'KNN imputation for a table',
    'description': 'Missing values imputation using the k-nearest neighbors method',
  })
  static async kNNImputationForTable(
    table: DG.DataFrame) {
    await runKNNImputer(table);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Linear Regression',
      'mlrole': 'train',
    },
    'name': 'trainLinearRegression',
    'outputs': [{'type': 'dynamic', 'name': 'model'}],
  })
  static async trainLinearRegression(
    df: DG.DataFrame,
    predictColumn: DG.Column): Promise<Uint8Array> {
    const features = df.columns;
    const params = await getLinearRegressionParams(features, predictColumn);

    return new Uint8Array(params.buffer);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Linear Regression',
      'mlrole': 'apply',
    },
    'name': 'applyLinearRegression',
  })
  static applyLinearRegression(
    df: DG.DataFrame,
    model: any): DG.DataFrame {
    const features = df.columns;
    const params = new Float32Array((model as Uint8Array).buffer);
    return DG.DataFrame.fromColumns([getPredictionByLinearRegression(features, params)]);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Linear Regression',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableLinearRegression',
  })
  static isApplicableLinearRegression(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    for (const col of df.columns) {
      if (!col.matches('numerical'))
        return false;
    }

    return predictColumn.matches('numerical');
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Linear Regression',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveLinearRegression',
  })
  static isInteractiveLinearRegression(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return df.rowCount <= 100000;
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Softmax',
      'mlrole': 'train',
    },
    'name': 'trainSoftmax',
    'outputs': [{'type': 'dynamic', 'name': 'model'}],
  })
  static async trainSoftmax(df: DG.DataFrame, predictColumn: DG.Column,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '1.0', 'min': '0.001', 'max': '20', 'description': 'Learning rate.'}}) rate: number,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '100', 'min': '1', 'max': '10000', 'step': '10', 'description': 'Fitting iterations count'}}) iterations: number,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '0.1', 'min': '0.0001', 'max': '1', 'description': 'Regularization rate.'}}) penalty: number,
    @grok.decorators.param({'options': {'category': 'Hyperparameters', 'initialValue': '0.001', 'min': '0.00001', 'max': '0.1', 'description': 'Fitting tolerance.'}}) tolerance: number): Promise<Uint8Array> {
    const features = df.columns;

    const model = new SoftmaxClassifier({
      classesCount: predictColumn.categories.length,
      featuresCount: features.length,
    });

    await model.fit(features, predictColumn, rate, iterations, penalty, tolerance);

    return model.toBytes();
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Softmax',
      'mlrole': 'apply',
    },
    'name': 'applySoftmax',
  })
  static applySoftmax(
    df: DG.DataFrame,
    model: any): DG.DataFrame {
    const features = df.columns;
    const unpackedModel = new SoftmaxClassifier(undefined, model);

    return DG.DataFrame.fromColumns([unpackedModel.predict(features)]);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Softmax',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableSoftmax',
  })
  static isApplicableSoftmax(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return SoftmaxClassifier.isApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Softmax',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveSoftmax',
  })
  static isInteractiveSoftmax(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return SoftmaxClassifier.isInteractive(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'PLS Regression',
      'mlrole': 'train',
    },
    'name': 'trainPLSRegression',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainPLSRegression(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'min': '1', 'max': '10', 'initialValue': '3', 'description': 'Number of latent components.'}}) components: number): Promise<Uint8Array> {
    const features = df.columns;

    const model = new PlsModel();
    await model.fit(
      features,
      predictColumn,
      Math.min(components, features.length),
    );

    return model.toBytes();
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'PLS Regression',
      'mlrole': 'apply',
    },
    'name': 'applyPLSRegression',
  })
  static applyPLSRegression(
    df: DG.DataFrame,
    model: any): DG.DataFrame {
    const unpackedModel = new PlsModel(model);
    return DG.DataFrame.fromColumns([unpackedModel.predict(df.columns)]);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'PLS Regression',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicablePLSRegression',
  })
  static isApplicablePLSRegression(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return PlsModel.isApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'PLS Regression',
      'mlrole': 'visualize',
    },
    'name': 'visualizePLSRegression',
  })
  static async visualizePLSRegression(
    df: DG.DataFrame,
    targetColumn: DG.Column,
    predictColumn: DG.Column,
    model: any): Promise<any> {
    const unpackedModel = new PlsModel(model);
    const viewers = unpackedModel.viewers();

    return viewers.map((v) => v.root);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'PLS Regression',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractivePLSRegression',
  })
  static isInteractivePLSRegression(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return PlsModel.isInteractive(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'XGBoost',
      'mlrole': 'train',
    },
    'name': 'trainXGBooster',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainXGBooster(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'min': '1', 'max': '100', 'initialValue': '20', 'description': 'Number of training iterations.'}}) iterations: number,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'Rate', 'min': '0', 'max': '1', 'initialValue': '0.3', 'description': 'Learning rate.'}}) eta: number,
    @grok.decorators.param({'type': 'int', 'options': {'min': '0', 'max': '20', 'initialValue': '6', 'description': 'Maximum depth of a tree.'}}) maxDepth: number,
    @grok.decorators.param({'type': 'double', 'options': {'min': '0', 'max': '100', 'initialValue': '1', 'description': 'L2 regularization term.'}}) lambda: number,
    @grok.decorators.param({'type': 'double', 'options': {'min': '0', 'max': '100', 'initialValue': '0', 'description': 'L1 regularization term.'}}) alpha: number): Promise<Uint8Array> {
    const features = df.columns;

    const booster = new XGBooster();
    await booster.fit(features, predictColumn, iterations, eta, maxDepth, lambda, alpha);

    return booster.toBytes();
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'XGBoost',
      'mlrole': 'apply',
    },
    'name': 'applyXGBooster',
  })
  static applyXGBooster(
    df: DG.DataFrame,
    model: any): DG.DataFrame {
    const unpackedModel = new XGBooster(model);
    return DG.DataFrame.fromColumns([unpackedModel.predict(df.columns)]);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'XGBoost',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveXGBooster',
  })
  static isInteractiveXGBooster(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return XGBooster.isInteractive(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'XGBoost',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableXGBooster',
  })
  static isApplicableXGBooster(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return XGBooster.isApplicable(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'top-menu': 'ML | Pareto Front...',
    'name': 'Pareto Front',
    'description': 'Perform optimization across multiple objectives: analyze trade-offs between conflicting objectives and identify Pareto-optimal points.',
  })
  static paretoFront(): void {
    const df = grok.shell.t;
    if (df === null) {
      grok.shell.warning('No dataframe is opened');
      return;
    }

    const optimizer = new ParetoOptimizer(df);
    optimizer.run();
    //runParetoOptimizer();
  }

  @grok.decorators.func({
    'name': 'Pareto front',
    'description': 'Pareto front viewer',
    'outputs': [{'name': 'result', 'type': 'viewer'}],
    'meta': {'icon': 'icons/pareto-front-viewer.svg', 'role': 'viewer'},
  })
  static paretoFrontViewer(): DG.Viewer {
    return new ParetoFrontViewer();
  }

  @grok.decorators.func({
    'top-menu': 'Chem | Calculate | Train pMPO...',
    'name': 'trainPmpo',
    'description': 'Train probabilistic multi-parameter optimization (pMPO) model',
  })
  static trainPmpo(): void {
    const df = grok.shell.t;
    if (df === null) {
      grok.shell.warning('No dataframe is opened');
      return;
    }

    if (!Pmpo.isTableValid(df))
      return;

    const pMPO = new Pmpo(df);
    pMPO.runTrainingApp();
  }

  @grok.decorators.func({'name': 'getPmpoAppItems', 'outputs': [{name: 'result', type: 'object'}]})
  static getPmpoAppItems(@grok.decorators.param({type: 'view'}) view: DG.TableView): any | null {
    const df = view.dataFrame;
    if (!Pmpo.isTableValid(df))
      return null;

    const pMPO = new Pmpo(df, view);

    return pMPO.getPmpoAppItems();
  }

  @grok.decorators.func({
    'name': 'generatePmpoDataset',
    'description': 'Generates syntethetic dataset oriented on the pMPO modeling',
    'outputs': [{name: 'Synthetic', type: 'dataframe'}],
  })
  static async generatePmpoDataset(@grok.decorators.param({'type': 'int'}) samples: number): Promise<DG.DataFrame> {
    const df = await getSynteticPmpoData(samples);
    df.name = 'Synthetic';
    return df;
  }
}
