// Predicitve tools based on the PLS method

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getPlsAnalysis} from './pls-tools';
import {LINK} from './pls-constants';
import {getPredictionByLinearRegression} from '../regression';

/** Acceptable column types */
const NUMERIC_TYPES = [
  DG.COLUMN_TYPE.INT,
  DG.COLUMN_TYPE.FLOAT,
  DG.COLUMN_TYPE.BIG_INT,
  DG.COLUMN_TYPE.QNUM,
];

// PLS ML specific constants
const EXTRA_ROWS = 1;
const COL_OFFSET = 2;
const BYTES_PER_MODEL_SIZE = 4;
const MIN_COLS_COUNT = 3;

enum TITLE {
  FEATURES = 'Feature names',
  PARAMS = 'Regression Coefficients',
  X_LOADING = 'x.loading.p',
  LOADINGS = 'Loadings',
  COMP = 'component',
  COMPS = 'components',
  EXPL_VAR = 'explained variance',
};

/** PLS regression modeling tool */
export class PLSmodel {
  /** Check applicability */
  static isApplicable(features: DG.ColumnList): boolean {
    for (const col of features) {
      if (!NUMERIC_TYPES.includes(col.type as DG.COLUMN_TYPE))
        return false;
    }

    return true;
  }

  private params: Float32Array | null = null;
  private featureNames: string[] | null = null;
  private loadings: Float32Array[] | null = null;

  constructor(packedModel?: Uint8Array) {
    if (packedModel) {
      try {
        // Extract model's bytes count
        const sizeArr = new Uint32Array(packedModel.buffer, 0, 1); // 1-st element is a size of model bytes
        const bytesCount = sizeArr[0];

        // Model's bytes
        const bytes = new Uint8Array(packedModel.buffer, BYTES_PER_MODEL_SIZE, bytesCount);

        // Model as dataframe
        const modelDf = DG.DataFrame.fromByteArray(bytes);
        const rowCount = modelDf.rowCount;
        const columns = modelDf.columns;
        const colsCount = columns.length;

        if (colsCount < MIN_COLS_COUNT)
          throw new Error('incorrect columns count');

        // Extract names of features
        this.featureNames = columns.byName(TITLE.FEATURES).toList();

        // Extract parameters of the linear model
        this.params = new Float32Array(rowCount);
        this.params.set(columns.byName(TITLE.PARAMS).getRawData());

        // Extract loadings
        const components = colsCount - COL_OFFSET;
        this.loadings = new Array<Float32Array>(components);

        for (let i = 0; i < components; ++i) {
          this.loadings[i] = new Float32Array(rowCount);
          this.loadings[i].set(columns.byIndex(i + COL_OFFSET).getRawData());
        }
      } catch (error) {
        throw new Error(`Failed to load model: ${(error instanceof Error ? error.message : 'the platform issue')}`);
      }
    }
  }

  /** Train model */
  public async fit(features: DG.ColumnList, target: DG.Column, components: number) {
    const analysis = await getPlsAnalysis({
      table: DG.DataFrame.fromColumns([target]),
      features: features,
      predict: target,
      components: components,
      names: null,
    });

    // 1. Names of features
    this.featureNames = features.names();
    this.featureNames.push('_'); // add extra item

    // 2. Regression coefficients
    this.params = this.getRegrCoeffs(features, target, analysis.regressionCoefficients);

    // 3. Loadings
    this.loadings = this.getLoadings(components, analysis.xLoadings);

    // 4. Compute explained variances
    this.computeExplVars(target.length, components, analysis.yLoadings);
  } // fit

  /** Return x-loadings with extra items reserved for explained variances */
  private getLoadings(components: number, loadingsCols: DG.Column[]): Float32Array[] {
    const res = Array<Float32Array>(components);
    const len = loadingsCols[0].length + EXTRA_ROWS;

    for (let i = 0; i < components; ++i) {
      res[i] = new Float32Array(len);
      res[i].set(loadingsCols[i].getRawData());
    }

    return res;
  }

  /** Return regression coefficients */
  private getRegrCoeffs(features: DG.ColumnList, target: DG.Column, regrCoefsCol: DG.Column): Float32Array {
    const featuresCount = features.length;
    const params = new Float32Array(featuresCount + EXTRA_ROWS);
    const paramsByPLS = regrCoefsCol.getRawData();

    let tmpSum = 0;

    for (let i = 0; i < featuresCount; ++i) {
      params[i] = paramsByPLS[i];
      tmpSum += paramsByPLS[i] * features.byIndex(i).stats.avg;
    }

    // compute bias
    params[featuresCount] = target.stats.avg - tmpSum;

    return params;
  }

  /** Return explained variances */
  private computeExplVars(samplesCount: number, components: number, yLoadings: DG.Column) {
    if (this.loadings === null)
      throw new Error('Failed to compute explained variances');

    const raw = yLoadings.getRawData();
    const featuresCount = this.loadings[0].length - 1;

    // Compute, source: the paper https://doi.org/10.1002/cem.2589
    let explVar = raw[0]**2 / samplesCount;

    this.loadings[0][featuresCount] = explVar;

    for (let comp = 1; comp < components; ++comp) {
      explVar += raw[comp]**2 / samplesCount;
      this.loadings[comp][featuresCount] = explVar;
    }
  }

  /** Return packed model */
  public toBytes(): Uint8Array {
    if ((this.featureNames === null) || (this.params === null) || (this.loadings === null))
      throw new Error('Failed to pack untrained model');

    // 1. Store model in dataframe
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings(TITLE.FEATURES, this.featureNames),
      DG.Column.fromFloat32Array(TITLE.PARAMS, this.params),
    ]);

    this.loadings.forEach((array, idx) => df.columns.add(DG.Column.fromFloat32Array(
      `${TITLE.X_LOADING} ${idx + 1}`,
      array,
    )));

    // 2. Pack model dataframe
    const modelBytes = df.toByteArray();
    const bytesCount = modelBytes.length;

    const packedModel = new Uint8Array(bytesCount + BYTES_PER_MODEL_SIZE);

    // 4 bytes for storing model's bytes count
    const sizeArr = new Uint32Array(packedModel.buffer, 0, 1);
    sizeArr[0] = bytesCount;

    // Store model's bytes
    packedModel.set(modelBytes, BYTES_PER_MODEL_SIZE);

    return packedModel;

    return new Uint8Array();
  }

  /** Return prediction */
  public predict(features: DG.ColumnList): DG.Column {
    if (this.params === null)
      throw new Error('Predicting failed: model is not trained');

    return getPredictionByLinearRegression(features, this.params);
  }

  /** Returns viewers */
  public viewers(): DG.Viewer[] {
    if ((this.featureNames === null) || (this.params === null) || (this.loadings === null))
      throw new Error('Failed to create viewers: untrained model');

    const viewers: DG.Viewer[] = [];

    const featuresCount = this.featureNames.length - 1;

    // Parameters and loadings dataframe
    const loadingsDf = DG.DataFrame.fromColumns([
      DG.Column.fromStrings(TITLE.FEATURES, this.featureNames.slice(0, -1)),
      DG.Column.fromFloat32Array(TITLE.PARAMS, this.params, featuresCount),
    ]);

    this.loadings.forEach((arr, idx) => loadingsDf.columns.add(
      DG.Column.fromFloat32Array(`${TITLE.X_LOADING} ${idx + 1}`, arr, featuresCount),
    ));

    const columns = loadingsDf.columns;
    const components = columns.length - COL_OFFSET;

    // Loading scatterplot
    viewers.push(DG.Viewer.scatterPlot(loadingsDf, {
      title: TITLE.LOADINGS,
      xColumnName: columns.byIndex(COL_OFFSET).name,
      yColumnName: columns.byIndex(COL_OFFSET + (components > 1 ? 1 : 0)).name,
      markerType: DG.MARKER_TYPE.CIRCLE,
      labels: TITLE.FEATURES,
      help: LINK.LOADINGS,
    }));

    // Regression coefficients barchart
    viewers.push(DG.Viewer.barChart(loadingsDf, {
      title: TITLE.PARAMS,
      splitColumnName: TITLE.FEATURES,
      valueColumnName: TITLE.PARAMS,
      valueAggrType: DG.AGG.AVG,
      help: LINK.COEFFS,
      showValueSelector: false,
      showStackSelector: false,
    }));

    // Explained variances dataframe
    const compNames = new Array<string>(components);
    const explVars = new Float32Array(components);

    compNames[0] = `${TITLE.COMP} 1`;
    explVars[0] = this.loadings[0][featuresCount];

    for (let i = 1; i < components; ++i) {
      compNames[i] = `${TITLE.COMPS} ${i + 1}`;
      explVars[i] = this.loadings[i][featuresCount];
    }

    // Explained variances barchart
    viewers.push(DG.Viewer.barChart(DG.DataFrame.fromColumns([
      DG.Column.fromStrings(TITLE.COMPS, compNames),
      DG.Column.fromFloat32Array(TITLE.EXPL_VAR, explVars),
    ]), {
      title: TITLE.EXPL_VAR,
      splitColumnName: TITLE.COMPS,
      valueColumnName: TITLE.EXPL_VAR,
      valueAggrType: DG.AGG.AVG,
      help: LINK.EXPL_VARS,
      showCategorySelector: false,
      showStackSelector: false,
      showValueSelector: false,
    }));

    return viewers;
  }
};
