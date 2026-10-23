/* eslint-disable camelcase */
/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {computePCA} from './eda-tools';
import {addPrefixToEachColumnName} from './eda-ui';

import {PLS_ANALYSIS} from './pls/pls-constants';
import {runMVA, runDemoMVA, getPlsAnalysis, PlsOutput, addMvaResults, initMvaModelViewer} from './pls/pls-tools';
import {runOneWayAnova} from './anova/anova-ui';
import {runTwoSampleTTest} from './ttest/ttest-ui';
import {runControlComparisons} from './control-comparisons/control-comparisons-ui';

import {getDbscanWorker} from '@datagrok-libraries/math';

import {DistanceAggregationMethod, DistanceAggregationMethods} from '@datagrok-libraries/ml/src/distance-matrix/types';
import {MultiColumnDimReductionEditor} from
  '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/multi-column-dim-reduction-editor';
import {multiColReduceDimensionality} from
  '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/reduce-dimensionality';
import {KnownMetrics} from '@datagrok-libraries/ml/src/typed-metrics';
import {DimReductionMethods} from '@datagrok-libraries/ml/src/multi-column-dimensionality-reduction/types';

import {runKNNImputer} from './missing-values-imputation/ui';
import {imputeColumns} from './missing-values-imputation/impute-columns';
import {DISTANCE_TYPE} from './missing-values-imputation/knn-imputer';
import {MCLEditor} from '@datagrok-libraries/ml/src/MCL/mcl-editor';
import {MCLViewer} from '@datagrok-libraries/ml/src/MCL/mcl-viewer';
import {MCLSerializableOptions} from '@datagrok-libraries/ml/src/MCL';

import {isLinearRegressionApplicable, isLinearRegressionInteractive, trainLinearRegression} from './regression';
import {isPlsApplicable, isPlsInteractive, plsRegressionViewers, trainPlsRegression} from './pls/pls-ml';
import {isSoftmaxApplicable, isSoftmaxInteractive, isSvmApplicable, isSvmInteractive, isXgboostApplicable,
  isXgboostInteractive, trainSoftmax, trainSvm, trainXgboost, isRandomForestApplicable, isRandomForestInteractive,
  trainRandomForest, isExtraTreesApplicable, isExtraTreesInteractive, trainExtraTrees,
  isHistGradientBoostingApplicable, isHistGradientBoostingInteractive, trainHistGradientBoosting,
  isCatBoostApplicable, isCatBoostInteractive, trainCatBoost, isDecisionTreeApplicable, isDecisionTreeInteractive,
  trainDecisionTree, isKnnApplicable, isKnnInteractive, trainKnn, isNaiveBayesApplicable, isNaiveBayesInteractive,
  trainNaiveBayes, isLdaApplicable, isLdaInteractive, trainLda, isMlpApplicable, isMlpInteractive,
  trainMlp} from './crux/methods';
import {cruxPredict, modelFeatures, readCruxManifest} from './crux/crux-model';

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
    'description': 'Principal component analysis (PCA).',
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
      // "column_list" might be passed as array of columns, not DG.ColumnList (slack: https://datagrok.slack.com/archives/C04BF3YM6CF/p1776968116891119)
      const featuresToPass = Array.isArray(features) ? DG.DataFrame.fromColumns(features).columns : features;

      const pcaTable = await computePCA(table, featuresToPass, components, center, scale);
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
    'meta': {'defaultPostProcessingFunction': 'true', 'role': 'dim-red-postprocessing-function'},
    'name': 'DBSCAN clustering',
    'tags': ['dim-red-postprocessing-function'],
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
      'role': 'dim-red-preprocessing-function',
    },
    'name': 'None (number)',
    'tags': ['dim-red-preprocessing-function'],
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
      'role': 'dim-red-preprocessing-function',
    },
    'tags': ['dim-red-preprocessing-function'],
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


  @grok.decorators.editor({tags: ['editor']})
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
    'name': 'MCL Clustering',
    'description': 'Markov clustering (MCL) is an unsupervised clustering algorithm for graphs based on simulation of stochastic flow.',
    'editor': 'EDA:GetMCLEditor',
    'outputs': [],
  })
  static async MCLClustering(
    df: DG.DataFrame,
    cols: DG.Column[],
    @grok.decorators.param({'type': 'list<string>', 'options': {'description': 'Distance metric per column used to measure similarity between rows.'}}) metrics: KnownMetrics[],
    @grok.decorators.param({'options': {'description': 'Relative weight of each column when combining per-column distances.'}}) weights: number[],
    @grok.decorators.param({'type': 'string', 'options': {'description': 'How per-column distances are aggregated into a single distance.'}}) aggregationMethod: DistanceAggregationMethod,
    @grok.decorators.param({'type': 'list<func>', 'options': {'description': 'Preprocessing function applied to each column before distances are computed.'}}) preprocessingFuncs: any[],
    @grok.decorators.param({'type': 'object', 'options': {'description': 'Arguments passed to the preprocessing functions.'}}) preprocessingFuncArgs: any[],
    @grok.decorators.param({'type': 'int', 'options': {'initialValue': '80', 'description': 'Similarity threshold (percentile): edges below it are dropped before clustering.'}}) threshold: number = 80,
    @grok.decorators.param({'type': 'int', 'options': {'initialValue': '10', 'description': 'Maximum number of expansion/inflation iterations.'}}) maxIterations: number = 10,
    @grok.decorators.param({'type': 'bool', 'options': {'initialValue': 'false', 'description': 'Run the computation on the GPU via WebGPU when available.'}}) useWebGPU: boolean = false,
    @grok.decorators.param({'type': 'double', 'options': {'initialValue': '2', 'description': 'Inflation factor controlling cluster granularity: higher values yield more, smaller clusters.'}}) inflate: number = 0,
    @grok.decorators.param({'type': 'int', 'options': {'initialValue': '5', 'description': 'Clusters smaller than this are merged into noise.'}}) minClusterSize: number = 5): Promise<MCLViewer> {
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
    'tags': ['viewer'],
    'name': 'MCL',
    'description': 'Markov clustering viewer',
  })
  static markovClusteringViewer(): MCLViewer {
    return new MCLViewer();
  }

  @grok.decorators.func({
    'name': 'PLS',
    'outputs': [{'name': 'plsResults', 'type': 'object'}],
    'description': 'Compute partial least squares (PLS) regression analysis components: prediction, regression coefficients, T- & U-scores, X-loadings.',
  })
  static async PLS(
    table: DG.DataFrame,
    @grok.decorators.param({'type': 'column_list', 'options': {'type': 'numerical'}}) features: DG.ColumnList,
    @grok.decorators.param({'type': 'column', 'options': {'type': 'numerical'}}) predict: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'initialValue': '3', 'description': 'Number of latent factors the model extracts from the predictors.'}}) components: number,
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
    'description': 'Add the PLS-based multivariate analysis results to the table.',
    'meta': {'role': 'transform'},
  })
  static async multivariateAnalysisTransform(
    table: DG.DataFrame,
    // names, not a column list: the platform passes column_list values of a replayed creation script as strings
    featureNames: string[],
    @grok.decorators.param({'type': 'column', 'options': {'type': 'numerical'}}) predict: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'description': 'Number of latent factors the model extracts from the predictors.'}}) components: number,
    @grok.decorators.param({'options': {'description': 'Include squared terms as additional predictors.'}}) isQuadratic: boolean,
    @grok.decorators.param({'options': {'description': 'Add just the PLS components.'}}) componentsOnly: boolean,
    xScoreNames: string[],
    yScoreNames: string[],
    predictionName: string,
    analysisTableName: string,
    explVarTableName: string): Promise<void> {
    await addMvaResults({
      table: table,
      features: DG.DataFrame.fromColumns(featureNames.map((name) => table.col(name)!)).columns,
      predict: predict,
      components: components,
      isQuadratic: isQuadratic,
      names: undefined,
    }, {
      xScores: xScoreNames,
      yScores: yScoreNames,
      prediction: predictionName,
      analysisTable: analysisTableName,
      explVarTable: explVarTableName,
    }, componentsOnly);
  }


  @grok.decorators.func({
    'description': 'Restore the multivariate analysis model viewer.',
  })
  static mvaModelInitFunction(
    @grok.decorators.param({'type': 'viewer'}) viewer: DG.Viewer): void {
    initMvaModelViewer(viewer);
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
    'top-menu': 'ML | Analyze | Group Comparison | T-test...',
    'name': 'T-test',
    'description': 'Two-sample t-test (Welch or Student) compares the means of a feature between two groups.',
  })
  static tTest(): void {
    runTwoSampleTTest();
  }


  @grok.decorators.func({
    'top-menu': 'ML | Analyze | Group Comparison | ANOVA...',
    'name': 'ANOVA',
    'description': 'One-way analysis of variances (ANOVA) determines whether the examined factor has a significant impact on the explored feature.',
  })
  static anova(): void {
    runOneWayAnova();
  }


  @grok.decorators.func({
    'top-menu': 'ML | Analyze | Group Comparison | Control Comparisons...',
    'name': 'Control comparisons',
    'description': 'Compare several groups against a single control (Dunnett\'s test or Holm-corrected Welch\'s t-tests).',
  })
  static controlComparisons(): void {
    runControlComparisons();
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
    'name': 'knnImpute',
    'description': 'Fills missing values in place using the k-nearest neighbors method; cells with no usable neighbors stay empty',
    'outputs': [{'type': 'dataframe', 'name': 'result'}],
  })
  static knnImpute(
    table: DG.DataFrame,
    @grok.decorators.param({'type': 'list<string>', 'options': {'description': 'Columns to fill, an empty list fills every feature with missing values'}}) columns: string[],
    @grok.decorators.param({'type': 'list<string>', 'options': {'description': 'Columns used to find the nearest rows'}}) features: string[],
    @grok.decorators.param({'type': 'int', 'options': {'initialValue': '4', 'min': '1', 'description': 'Nearest rows used to fill a value.'}}) neighbors: number,
    @grok.decorators.param({'type': 'string', 'options': {'choices': ['Euclidean', 'Manhattan'], 'initialValue': 'Euclidean', 'description': 'How the distance between rows is measured.'}}) distance: string): DG.DataFrame {
    const distances: string[] = Object.values(DISTANCE_TYPE);
    if (!distances.includes(distance))
      throw new Error(`knnImpute: unknown distance "${distance}"; use ${distances.join(' or ')}`);
    imputeColumns(table, columns?.length ? columns : features, features, neighbors, distance as DISTANCE_TYPE);
    return table;
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
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Iterations', 'min': '1', 'step': '50', 'max': '10000', 'initialValue': '1000', 'description': 'Largest number of solver iterations (elastic net).'}}) iterations: number,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'L1', 'min': '0', 'max': '100', 'initialValue': '0', 'description': 'L1 (Lasso) regularization term. 0 means plain ordinary least squares.'}}) alpha: number,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'L2', 'min': '0', 'max': '100', 'initialValue': '0', 'description': 'L2 (Ridge) regularization term. 0 means plain ordinary least squares.'}}) lambda: number): Promise<Uint8Array> {
    return trainLinearRegression(df, predictColumn, alpha, lambda, iterations);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Linear Regression',
      'mlrole': 'apply',
    },
    'name': 'Apply Linear Regression',
    'description': 'Predict the target for a table using a trained linear regression model.',
  })
  static async applyLinearRegression(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained linear regression model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
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
    return isLinearRegressionApplicable(df.columns, predictColumn);
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
    return isLinearRegressionInteractive(df.columns, predictColumn);
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
    @grok.decorators.param({'options': {'initialValue': '100', 'min': '1', 'max': '10000', 'step': '10', 'description': 'Largest number of training steps before training stops.'}}) iterations: number,
    @grok.decorators.param({'options': {'initialValue': '0.1', 'min': '0.0001', 'max': '1', 'description': 'Regularization rate.'}}) penalty: number,
    @grok.decorators.param({'options': {'initialValue': '0.001', 'min': '0.00001', 'max': '0.1', 'description': 'Smallest improvement worth continuing training for.'}}) tolerance: number): Promise<Uint8Array> {
    return trainSoftmax(df, predictColumn, iterations, penalty, tolerance);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Softmax',
      'mlrole': 'apply',
    },
    'name': 'Apply Softmax',
    'description': 'Classify the rows of a table using a trained softmax (multinomial logistic regression) model.',
  })
  static async applySoftmax(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained softmax classifier model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
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
    return isSoftmaxApplicable(df.columns, predictColumn);
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
    return isSoftmaxInteractive(df.columns, predictColumn);
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
    return trainPlsRegression(df, predictColumn, components);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'PLS Regression',
      'mlrole': 'apply',
    },
    'name': 'Apply PLS Regression',
    'description': 'Predict the target for a table using a trained partial least squares (PLS) regression model.',
  })
  static async applyPLSRegression(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained PLS regression model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
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
    return isPlsApplicable(df.columns, predictColumn);
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
    const manifest = readCruxManifest(model as Uint8Array);
    const features = df.clone(null, modelFeatures(df, manifest).map((col) => col.name));
    const components = manifest.params.n_components as number;
    return (await plsRegressionViewers(features, targetColumn, components)).map((v) => v.root);
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
    return isPlsInteractive(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'XGBoost',
      'mlrole': 'train',
      'mlnativenans': 'true',
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
    return trainXgboost(df, predictColumn, iterations, eta, maxDepth, lambda, alpha);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'XGBoost',
      'mlrole': 'apply',
    },
    'name': 'Apply XGBoost',
    'description': 'Predict the target for a table using a trained XGBoost gradient-boosting model.',
  })
  static async applyXGBooster(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained XGBoost model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
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
    return isXgboostInteractive(df.columns, predictColumn);
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
    return isXgboostApplicable(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'SVM',
      'mlrole': 'train',
      'mlhasrows': 'true',
    },
    'name': 'trainSVM',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'string', 'options': {'choices': ['Linear', 'Polynomial', 'RBF', 'Sigmoid'], 'initialValue': 'RBF', 'description': 'Kernel function: linear, polynomial, RBF (radial basis function), or sigmoid.'}}) kernel: string,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'Penalty', 'min': '0.001', 'max': '1000', 'initialValue': '1', 'description': 'Penalty for misclassifying training points. Higher values fit the training data more tightly and risk overfitting, lower values give a smoother, more general model.'}}) cost: number,
    @grok.decorators.param({'type': 'double', 'options': {'min': '0', 'max': '100', 'initialValue': '0', 'description': 'Kernel coefficient. Used by the RBF, polynomial, and sigmoid kernels. Leave 0 to use 1/number-of-features.'}}) gamma: number,
    @grok.decorators.param({'type': 'int', 'options': {'min': '1', 'max': '10', 'initialValue': '3', 'description': 'Degree of the polynomial kernel. Higher values give a more flexible, curved decision boundary. Used by the polynomial kernel only.'}}) degree: number,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'Offset', 'min': '-100', 'max': '100', 'initialValue': '0', 'description': 'Constant added inside the kernel before the polynomial power or the sigmoid. Used by the polynomial and sigmoid kernels only.'}}) coef0: number,
    @grok.decorators.param({'type': 'double', 'options': {'min': '0', 'max': '10', 'initialValue': '0.1', 'description': 'Width of the error-insensitive tube in SVR regression. Training errors smaller than this are not penalized, so larger values give a smoother model with fewer support vectors. Used for regression targets only.'}}) epsilon: number): Promise<Uint8Array> {
    return trainSvm(df, predictColumn, kernel, cost, gamma, degree, coef0, epsilon);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'SVM',
      'mlrole': 'apply',
    },
    'name': 'Apply SVM',
    'description': 'Predict the target for a table using a trained support vector machine (SVM) model.',
  })
  static async applySVM(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained SVM model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'SVM',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableSVM',
  })
  static isApplicableSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isSvmApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'SVM',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveSVM',
  })
  static isInteractiveSVM(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isSvmInteractive(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'Random Forest',
      'mlrole': 'train',
      'mlnativenans': 'true',
    },
    'name': 'trainRandomForest',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainRandomForest(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Trees', 'min': '1', 'max': '1000', 'initialValue': '100', 'description': 'Number of trees in the forest.'}}) trees: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Max depth', 'min': '1', 'max': '50', 'initialValue': '10', 'description': 'Largest depth of a tree. Bounds the model size.'}}) maxDepth: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Min rows per leaf', 'min': '1', 'max': '1000', 'initialValue': '1', 'description': 'Smallest number of training rows a leaf may hold.'}}) minSamplesLeaf: number,
    @grok.decorators.param({'type': 'string', 'options': {'caption': 'Features per split', 'choices': ['auto', 'sqrt', 'log2', 'all'], 'initialValue': 'auto', 'description': 'Features tried at each split: auto is sqrt for classification and all for regression.'}}) maxFeatures: string,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Seed', 'min': '0', 'max': '2147483647', 'initialValue': '0', 'description': 'Random seed. The same seed gives the same model.'}}) seed: number): Promise<Uint8Array> {
    return trainRandomForest(df, predictColumn, trees, maxDepth, minSamplesLeaf, maxFeatures, seed);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Random Forest',
      'mlrole': 'apply',
    },
    'name': 'Apply Random Forest',
    'description': 'Predict the target for a table using a trained random forest model.',
  })
  static async applyRandomForest(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained random forest model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Random Forest',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableRandomForest',
  })
  static isApplicableRandomForest(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isRandomForestApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Random Forest',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveRandomForest',
  })
  static isInteractiveRandomForest(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isRandomForestInteractive(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'Extra Trees',
      'mlrole': 'train',
      'mlnativenans': 'true',
    },
    'name': 'trainExtraTrees',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainExtraTrees(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Trees', 'min': '1', 'max': '1000', 'initialValue': '100', 'description': 'Number of trees in the forest.'}}) trees: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Max depth', 'min': '1', 'max': '50', 'initialValue': '10', 'description': 'Largest depth of a tree. Bounds the model size.'}}) maxDepth: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Min rows per leaf', 'min': '1', 'max': '1000', 'initialValue': '1', 'description': 'Smallest number of training rows a leaf may hold.'}}) minSamplesLeaf: number,
    @grok.decorators.param({'type': 'string', 'options': {'caption': 'Features per split', 'choices': ['auto', 'sqrt', 'log2', 'all'], 'initialValue': 'auto', 'description': 'Features tried at each split: auto is sqrt for classification and all for regression.'}}) maxFeatures: string,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Seed', 'min': '0', 'max': '2147483647', 'initialValue': '0', 'description': 'Random seed. The same seed gives the same model.'}}) seed: number): Promise<Uint8Array> {
    return trainExtraTrees(df, predictColumn, trees, maxDepth, minSamplesLeaf, maxFeatures, seed);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Extra Trees',
      'mlrole': 'apply',
    },
    'name': 'Apply Extra Trees',
    'description': 'Predict the target for a table using a trained extra trees model.',
  })
  static async applyExtraTrees(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained extra trees model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Extra Trees',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableExtraTrees',
  })
  static isApplicableExtraTrees(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isExtraTreesApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Extra Trees',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveExtraTrees',
  })
  static isInteractiveExtraTrees(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isExtraTreesInteractive(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'Hist Gradient Boosting',
      'mlrole': 'train',
      'mlnativenans': 'true',
    },
    'name': 'trainHistGradientBoosting',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainHistGradientBoosting(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Iterations', 'min': '1', 'max': '1000', 'initialValue': '100', 'description': 'Boosting iterations (trees per class).'}}) iterations: number,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'Learning rate', 'min': '0.001', 'max': '1', 'initialValue': '0.1', 'description': 'Shrinkage of each tree\'s contribution.'}}) learningRate: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Max leaves', 'min': '2', 'max': '255', 'initialValue': '31', 'description': 'Largest number of leaves per tree.'}}) maxLeafNodes: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Max depth', 'min': '0', 'max': '50', 'initialValue': '0', 'description': 'Largest depth of a tree. 0 means no limit.'}}) maxDepth: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Min rows per leaf', 'min': '1', 'max': '1000', 'initialValue': '20', 'description': 'Smallest number of training rows a leaf may hold.'}}) minSamplesLeaf: number,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'L2', 'min': '0', 'max': '100', 'initialValue': '0', 'description': 'L2 penalty of the leaf values.'}}) l2: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Seed', 'min': '0', 'max': '2147483647', 'initialValue': '0', 'description': 'Random seed. The same seed gives the same model.'}}) seed: number): Promise<Uint8Array> {
    return trainHistGradientBoosting(df, predictColumn, iterations, learningRate, maxLeafNodes, maxDepth, minSamplesLeaf,
      l2, seed);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Hist Gradient Boosting',
      'mlrole': 'apply',
    },
    'name': 'Apply Hist Gradient Boosting',
    'description': 'Predict the target for a table using a trained histogram gradient boosting model.',
  })
  static async applyHistGradientBoosting(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained histogram gradient boosting model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Hist Gradient Boosting',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableHistGradientBoosting',
  })
  static isApplicableHistGradientBoosting(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isHistGradientBoostingApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Hist Gradient Boosting',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveHistGradientBoosting',
  })
  static isInteractiveHistGradientBoosting(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isHistGradientBoostingInteractive(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'CatBoost',
      'mlrole': 'train',
      'mlnativenans': 'true',
    },
    'name': 'trainCatBoost',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainCatBoost(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Iterations', 'min': '1', 'max': '5000', 'initialValue': '100', 'description': 'Boosting iterations (trees).'}}) iterations: number,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'Learning rate', 'min': '0.001', 'max': '1', 'initialValue': '0.1', 'description': 'Shrinkage of each tree\'s contribution.'}}) learningRate: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Depth', 'min': '1', 'max': '16', 'initialValue': '6', 'description': 'Depth of the symmetric trees.'}}) depth: number,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'L2', 'min': '0', 'max': '100', 'initialValue': '3', 'description': 'L2 regularization of the leaf values.'}}) l2LeafReg: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Seed', 'min': '0', 'max': '2147483647', 'initialValue': '0', 'description': 'Random seed. The same seed gives the same model.'}}) seed: number): Promise<Uint8Array> {
    return trainCatBoost(df, predictColumn, iterations, learningRate, depth, l2LeafReg, seed);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'CatBoost',
      'mlrole': 'apply',
    },
    'name': 'Apply CatBoost',
    'description': 'Predict the target for a table using a trained CatBoost gradient-boosting model.',
  })
  static async applyCatBoost(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained CatBoost model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'CatBoost',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableCatBoost',
  })
  static isApplicableCatBoost(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isCatBoostApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'CatBoost',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveCatBoost',
  })
  static isInteractiveCatBoost(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isCatBoostInteractive(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'Decision Tree',
      'mlrole': 'train',
      'mlnativenans': 'true',
    },
    'name': 'trainDecisionTree',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainDecisionTree(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Max depth', 'min': '1', 'max': '50', 'initialValue': '10', 'description': 'Largest depth of the tree. Bounds the model size.'}}) maxDepth: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Min rows per leaf', 'min': '1', 'max': '1000', 'initialValue': '1', 'description': 'Smallest number of training rows a leaf may hold.'}}) minSamplesLeaf: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Seed', 'min': '0', 'max': '2147483647', 'initialValue': '0', 'description': 'Random seed. The same seed gives the same model.'}}) seed: number): Promise<Uint8Array> {
    return trainDecisionTree(df, predictColumn, maxDepth, minSamplesLeaf, seed);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Decision Tree',
      'mlrole': 'apply',
    },
    'name': 'Apply Decision Tree',
    'description': 'Predict the target for a table using a trained decision tree model.',
  })
  static async applyDecisionTree(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained decision tree model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Decision Tree',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableDecisionTree',
  })
  static isApplicableDecisionTree(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isDecisionTreeApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Decision Tree',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveDecisionTree',
  })
  static isInteractiveDecisionTree(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isDecisionTreeInteractive(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'k-Nearest Neighbors',
      'mlrole': 'train',
      'mlhasrows': 'true',
    },
    'name': 'trainKnn',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainKnn(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Neighbors', 'min': '1', 'max': '100', 'initialValue': '5', 'description': 'Training rows that vote for a prediction.'}}) neighbors: number,
    @grok.decorators.param({'type': 'string', 'options': {'caption': 'Weights', 'choices': ['Uniform', 'Distance'], 'initialValue': 'Uniform', 'description': 'Uniform: every neighbor counts the same. Distance: closer neighbors count more.'}}) weights: string): Promise<Uint8Array> {
    return trainKnn(df, predictColumn, neighbors, weights);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'k-Nearest Neighbors',
      'mlrole': 'apply',
    },
    'name': 'Apply k-Nearest Neighbors',
    'description': 'Predict the target for a table using a trained k-nearest neighbors model.',
  })
  static async applyKnn(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained k-nearest neighbors model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'k-Nearest Neighbors',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableKnn',
  })
  static isApplicableKnn(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isKnnApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'k-Nearest Neighbors',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveKnn',
  })
  static isInteractiveKnn(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isKnnInteractive(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'Naive Bayes',
      'mlrole': 'train',
    },
    'name': 'trainNaiveBayes',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainNaiveBayes(
    df: DG.DataFrame,
    predictColumn: DG.Column): Promise<Uint8Array> {
    return trainNaiveBayes(df, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Naive Bayes',
      'mlrole': 'apply',
    },
    'name': 'Apply Naive Bayes',
    'description': 'Classify the rows of a table using a trained Gaussian naive Bayes model.',
  })
  static async applyNaiveBayes(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained naive Bayes classifier model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Naive Bayes',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableNaiveBayes',
  })
  static isApplicableNaiveBayes(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isNaiveBayesApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Naive Bayes',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveNaiveBayes',
  })
  static isInteractiveNaiveBayes(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isNaiveBayesInteractive(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'LDA',
      'mlrole': 'train',
    },
    'name': 'trainLda',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainLda(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'Shrinkage', 'min': '0', 'max': '1', 'initialValue': '0', 'description': 'Covariance shrinkage from 0 (none) to 1. Helps when features outnumber rows.'}}) shrinkage: number): Promise<Uint8Array> {
    return trainLda(df, predictColumn, shrinkage);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'LDA',
      'mlrole': 'apply',
    },
    'name': 'Apply LDA',
    'description': 'Classify the rows of a table using a trained linear discriminant analysis (LDA) model.',
  })
  static async applyLda(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained LDA classifier model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'LDA',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableLda',
  })
  static isApplicableLda(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isLdaApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'LDA',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveLda',
  })
  static isInteractiveLda(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isLdaInteractive(df.columns, predictColumn);
  }

  @grok.decorators.func({
    'meta': {
      'mlname': 'Neural Network',
      'mlrole': 'train',
    },
    'name': 'trainMlp',
    'outputs': [{'name': 'model', 'type': 'dynamic'}],
  })
  static async trainMlp(
    df: DG.DataFrame,
    predictColumn: DG.Column,
    @grok.decorators.param({'type': 'string', 'options': {'caption': 'Hidden layers', 'initialValue': '100', 'description': 'Comma-separated sizes of the hidden layers, e.g. 100 or 64,32.'}}) hiddenLayers: string,
    @grok.decorators.param({'type': 'string', 'options': {'caption': 'Activation', 'choices': ['ReLU', 'Tanh', 'Logistic', 'Identity'], 'initialValue': 'ReLU', 'description': 'Activation function of the hidden layers.'}}) activation: string,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Iterations', 'min': '1', 'max': '2000', 'initialValue': '200', 'description': 'Largest number of training epochs.'}}) iterations: number,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'Learning rate', 'min': '0.00001', 'max': '1', 'initialValue': '0.001', 'description': 'Step size of the weight updates.'}}) learningRate: number,
    @grok.decorators.param({'type': 'double', 'options': {'caption': 'L2', 'min': '0', 'max': '10', 'initialValue': '0.0001', 'description': 'L2 penalty of the weights.'}}) alpha: number,
    @grok.decorators.param({'type': 'int', 'options': {'caption': 'Seed', 'min': '0', 'max': '2147483647', 'initialValue': '0', 'description': 'Random seed. The same seed gives the same model.'}}) seed: number): Promise<Uint8Array> {
    return trainMlp(df, predictColumn, hiddenLayers, activation, iterations, learningRate, alpha, seed);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Neural Network',
      'mlrole': 'apply',
    },
    'name': 'Apply Neural Network',
    'description': 'Predict the target for a table using a trained neural network (multilayer perceptron) model.',
  })
  static async applyMlp(
    df: DG.DataFrame,
    @grok.decorators.param({'options': {'description': 'Trained neural network model to apply.'}}) model: any): Promise<DG.DataFrame> {
    return cruxPredict(df, model as Uint8Array);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Neural Network',
      'mlrole': 'isApplicable',
    },
    'name': 'isApplicableMlp',
  })
  static isApplicableMlp(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isMlpApplicable(df.columns, predictColumn);
  }


  @grok.decorators.func({
    'meta': {
      'mlname': 'Neural Network',
      'mlrole': 'isInteractive',
    },
    'name': 'isInteractiveMlp',
  })
  static isInteractiveMlp(
    df: DG.DataFrame,
    predictColumn: DG.Column): boolean {
    return isMlpInteractive(df.columns, predictColumn);
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
    'tags': ['viewer'],
  })
  static paretoFrontViewer(): DG.Viewer {
    return new ParetoFrontViewer();
  }

  @grok.decorators.func({
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
  static async getPmpoAppItems(@grok.decorators.param({type: 'view'}) view: DG.TableView): Promise<any | null> {
    const df = view.dataFrame;
    if (!Pmpo.isTableValid(df))
      return null;

    const pMPO = new Pmpo(df, view);
    return pMPO.getPmpoAppItems();
  }

  @grok.decorators.func({
    'name': 'Generate pMPO Dataset',
    'description': 'Generate a synthetic dataset oriented on probabilistic multi-parameter optimization (pMPO) modeling.',
    'outputs': [{name: 'Synthetic', type: 'dataframe'}],
  })
  static async generatePmpoDataset(@grok.decorators.param({'type': 'int', 'options': {'description': 'Number of rows (samples) to generate.'}}) samples: number): Promise<DG.DataFrame> {
    const df = await getSynteticPmpoData(samples, false);
    df.name = 'Synthetic';
    return df;
  }
}
