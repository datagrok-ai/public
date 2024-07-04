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
const SHIFT = 4;
const MIN_LOADINGS = 1;
const BYTES_PER_MODEL_SIZE = 4;
const MIN_COLS_COUNT = SHIFT + MIN_LOADINGS;

/** Titles */
enum TITLE {
  FEATURES = 'Feature names',
  PARAMS = 'Regression Coefficients',
  X_LOADING = 'x.loading.p',
  LOADINGS = 'Loadings',
  COMP = 'component',
  COMPS = 'components',
  EXPL_VAR = 'explained variance',
  AVG = 'average',
  STDEV = 'sigma',
};

/** Model specification */
type PlsModelSpecification = {
  params: Float32Array,
  names: string[],
  loadings: Float32Array[],
  avgs: Float32Array,
  stdevs: Float32Array,
  dim: number,
  components: number,
}

/** PLS regression modeling tool */
export class PlsModel {
  /** Check applicability */
  static isApplicable(features: DG.ColumnList): boolean {
    for (const col of features) {
      if (!NUMERIC_TYPES.includes(col.type as DG.COLUMN_TYPE))
        return false;
    }

    return true;
  }

  /** Specification of the PLS model */
  private specn: PlsModelSpecification | null = null;

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
        const featureNames = columns.byName(TITLE.FEATURES).toList();

        // Extract parameters of the linear model
        const params = new Float32Array(rowCount);
        params.set(columns.byName(TITLE.PARAMS).getRawData());

        // Extract averages & stdevs
        const featuresAvgs = columns.byName(TITLE.AVG).getRawData() as Float32Array;
        const featuresStdevs = columns.byName(TITLE.STDEV).getRawData() as Float32Array;

        // Extract loadings
        const components = colsCount - SHIFT;
        const loadings = new Array<Float32Array>(components);

        for (let i = 0; i < components; ++i) {
          loadings[i] = new Float32Array(rowCount);
          loadings[i].set(columns.byIndex(i + SHIFT).getRawData());
        }

        this.specn = {
          params: params,
          loadings: loadings,
          names: featureNames,
          avgs: featuresAvgs,
          stdevs: featuresStdevs,
          dim: rowCount - EXTRA_ROWS,
          components: colsCount - SHIFT,
        };
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
    const featureNames = features.names();
    featureNames.push('_'); // add extra item

    // 2. Regression coefficients
    const params = this.getRegrCoeffs(features, target, analysis.regressionCoefficients);

    // 3. Loadings
    const loadings = this.getLoadings(components, analysis.xLoadings);

    // 4. Features stats
    const stats = this.extractStats(features);

    // 5. Model specification
    this.specn = {
      names: featureNames,
      params: params,
      avgs: stats.avgs,
      stdevs: stats.stdevs,
      loadings: loadings,
      components: components,
      dim: features.length,
    };

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
    const dim = features.length;
    const params = new Float32Array(dim + EXTRA_ROWS);
    const paramsByPLS = regrCoefsCol.getRawData();

    let tmpSum = 0;

    for (let i = 0; i < dim; ++i) {
      params[i] = paramsByPLS[i];
      tmpSum += paramsByPLS[i] * features.byIndex(i).stats.avg;
    }

    // compute bias
    params[dim] = target.stats.avg - tmpSum;

    return params;
  }

  /** Return explained variances */
  private computeExplVars(samplesCount: number, components: number, yLoadings: DG.Column) {
    if (this.specn === null)
      throw new Error('Failed to compute explained variances');

    const raw = yLoadings.getRawData();
    const dim = this.specn.loadings[0].length - EXTRA_ROWS;

    // Compute, source: the paper https://doi.org/10.1002/cem.2589
    let explVar = raw[0]**2 / samplesCount;

    this.specn.loadings[0][dim] = explVar;

    for (let comp = 1; comp < components; ++comp) {
      explVar += raw[comp]**2 / samplesCount;
      this.specn.loadings[comp][dim] = explVar;
    }
  }

  /** Return packed model */
  public toBytes(): Uint8Array {
    if (this.specn === null)
      throw new Error('Failed to pack untrained model');

    // 1. Store model in dataframe
    const df = DG.DataFrame.fromColumns([
      DG.Column.fromStrings(TITLE.FEATURES, this.specn.names),
      DG.Column.fromFloat32Array(TITLE.PARAMS, this.specn.params),
      DG.Column.fromFloat32Array(TITLE.AVG, this.specn.avgs),
      DG.Column.fromFloat32Array(TITLE.STDEV, this.specn.stdevs),
    ]);

    this.specn.loadings.forEach((array, idx) => df.columns.add(DG.Column.fromFloat32Array(
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
  }

  /** Return prediction */
  public predict(features: DG.ColumnList): DG.Column {
    if (this.specn === null)
      throw new Error('Predicting failed: model is not trained');

    return getPredictionByLinearRegression(features, this.specn.params);
  }

  /** Return loadings and regression coefficients viewers */
  private loadingsParamsViewers(): DG.Viewer[] {
    if (this.specn === null)
      throw new Error('Failed to create loadings and parameters viewers: untrained model');

    const viewers: DG.Viewer[] = [];

    const dim = this.specn.dim;

    // Parameters and loadings dataframe
    const loadingsDf = DG.DataFrame.fromColumns([
      DG.Column.fromStrings(TITLE.FEATURES, this.specn.names.slice(0, -1)),
      DG.Column.fromFloat32Array(TITLE.PARAMS, this.specn.params, dim),
    ]);

    const columns = loadingsDf.columns;
    const shift = columns.length;
    const components = this.specn.components;

    this.specn.loadings.forEach((arr, idx) => loadingsDf.columns.add(
      DG.Column.fromFloat32Array(`${TITLE.X_LOADING} ${idx + 1}`, arr, dim),
    ));

    // Loading scatterplot
    viewers.push(DG.Viewer.scatterPlot(loadingsDf, {
      title: TITLE.LOADINGS,
      xColumnName: columns.byIndex(shift).name,
      yColumnName: columns.byIndex(shift + (components > 1 ? 1 : 0)).name,
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

    return viewers;
  } // getLoadingsParamsViewers

  /** Return explained variances viewer */
  private explVarsViewer(): DG.Viewer {
    if (this.specn === null)
      throw new Error('Failed to create exaplained variances viewer: untrained model');

    const components = this.specn.components;
    const dim = this.specn.dim;

    const compNames = new Array<string>(components);
    const explVars = new Float32Array(components);

    compNames[0] = `${TITLE.COMP} 1`;
    explVars[0] = this.specn.loadings[0][dim];

    for (let i = 1; i < components; ++i) {
      compNames[i] = `${TITLE.COMPS} ${i + 1}`;
      explVars[i] = this.specn.loadings[i][dim];
    }

    return DG.Viewer.barChart(DG.DataFrame.fromColumns([
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
    });
  }

  /** Returns viewers */
  public viewers(): DG.Viewer[] {
    if (this.specn === null)
      throw new Error('Failed to create viewers: untrained model');

    const viewers = this.loadingsParamsViewers();
    viewers.push(this.explVarsViewer());

    return viewers;
  }

  /** Extract average & standard deviations of features */
  private extractStats(features: DG.ColumnList) {
    const lenght = features.length + EXTRA_ROWS;

    const avgs = new Float32Array(lenght);
    const stdevs = new Float32Array(lenght);

    let i = 0;
    let stats: DG.Stats;

    for (const col of features) {
      stats = col.stats;
      avgs[i] = stats.avg;
      stdevs[i] = stats.stdev;
      ++i;
    }

    return {
      avgs: avgs,
      stdevs: stdevs,
    };
  }
};
