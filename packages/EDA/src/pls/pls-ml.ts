// Predicitve tools based on the PLS method

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {getPlsAnalysis} from './pls-tools';
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

enum TITLES {
  FEATURES = 'Feature names',
  PARAMS = 'Parameters',
  LOADING = 'x-loading',
};

/** PLS regression modeling tool */
export class PLSmodel {
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
        this.featureNames = columns.byName(TITLES.FEATURES).toList();

        // Extract parameters of the linear model
        this.params = new Float32Array(rowCount);
        this.params.set(columns.byName(TITLES.PARAMS).getRawData());

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
      DG.Column.fromStrings(TITLES.FEATURES, this.featureNames),
      DG.Column.fromFloat32Array(TITLES.PARAMS, this.params),
    ]);

    this.loadings.forEach((array, idx) => df.columns.add(DG.Column.fromFloat32Array(
      `${TITLES.LOADING} ${idx + 1}`,
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

  /** Check applicability */
  static isApplicable(features: DG.ColumnList): boolean {
    for (const col of features) {
      if (!NUMERIC_TYPES.includes(col.type as DG.COLUMN_TYPE))
        return false;
    }

    return true;
  }
};
