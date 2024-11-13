// Predicitve tools based on the PLS method

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {TITLE, RESULT_NAMES} from './pls-constants';
import {getPlsAnalysis, PlsOutput, getLines} from './pls-tools';
import {LINK} from './pls-constants';
import {getPredictionByLinearRegression} from '../regression';

// PLS ML specific constants
const EXTRA_ROWS = 1;
const SHIFT = 2;
const MIN_LOADINGS = 1;
const MIN_COLS_COUNT = SHIFT + MIN_LOADINGS;
const SIZE_ARR_LEN = 2;
const MODEL_IDX = 0;
const SCORES_IDX = 1;
const BYTES_PER_SIZES = SIZE_ARR_LEN * 4;
const BLOCK_SIZE = 64;

/** Interactivity tresholds */
enum INTERACTIVITY {
  MAX_SAMLPES = 100000,
  MAX_FEATURES = 1000,
};

/** Model specification */
type PlsModelSpecification = {
  params: Float32Array,
  names: string[],
  loadings: Float32Array[],
  dim: number,
  components: number,
  scores: DG.DataFrame,
}

/** PLS regression modeling tool */
export class PlsModel {
  /** Check applicability */
  static isApplicable(features: DG.ColumnList, predictColumn: DG.Column): boolean {
    for (const col of features) {
      if (!col.matches('numerical'))
        return false;
    }
    if (!predictColumn.matches('numerical'))
      return false;

    return true;
  }

  /** Check interactivity */
  static isInteractive(features: DG.ColumnList, predictColumn: DG.Column): boolean {
    return (features.length <= INTERACTIVITY.MAX_FEATURES) &&
      (predictColumn.length <= INTERACTIVITY.MAX_SAMLPES);
  }

  /** Specification of the PLS model */
  private specn: PlsModelSpecification | null = null;

  constructor(packedModel?: Uint8Array) {
    if (packedModel) {
      try {
        // Extract model's bytes count
        const sizeArr = new Uint32Array(packedModel.buffer, 0, SIZE_ARR_LEN); // 1-st element is a size of model bytes
        const modelDfBytesCount = sizeArr[MODEL_IDX];
        const scoresDfBytesCount = sizeArr[SCORES_IDX];

        // Model's bytes
        const modelBytes = new Uint8Array(packedModel.buffer, BYTES_PER_SIZES, modelDfBytesCount);

        // Model as dataframe
        const modelDf = DG.DataFrame.fromByteArray(modelBytes);
        const rowCount = modelDf.rowCount;
        const columns = modelDf.columns;
        const colsCount = columns.length;

        // Scores
        const scoresBytes = new Uint8Array(packedModel.buffer, BYTES_PER_SIZES + modelDfBytesCount, scoresDfBytesCount);
        const scores = DG.DataFrame.fromByteArray(scoresBytes);

        if (colsCount < MIN_COLS_COUNT)
          throw new Error('incorrect columns count');

        // Extract names of features
        const featureNames = columns.byName(TITLE.FEATURES).toList();

        // Extract parameters of the linear model
        const params = new Float32Array(rowCount);
        params.set(columns.byName(TITLE.REGR_COEFS).getRawData());

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
          dim: rowCount - EXTRA_ROWS,
          components: colsCount - SHIFT,
          scores: scores,
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
      names: undefined,
    });

    // 1. Names of features
    const featureNames = features.names();
    featureNames.push('_'); // add extra item

    // 2. Regression coefficients
    const params = this.getRegrCoeffs(features, target, analysis.regressionCoefficients);

    // 3. Loadings
    const loadings = this.getLoadings(components, analysis.xLoadings);

    // 4. Model specification
    this.specn = {
      names: featureNames,
      params: params,
      loadings: loadings,
      components: components,
      dim: features.length,
      scores: this.getScoresDf(analysis),
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

    // 1. Store model params in dataframe
    const modelDf = DG.DataFrame.fromColumns([
      DG.Column.fromStrings(TITLE.FEATURES, this.specn.names),
      DG.Column.fromFloat32Array(TITLE.REGR_COEFS, this.specn.params),
    ]);

    this.specn.loadings.forEach((array, idx) => modelDf.columns.add(DG.Column.fromFloat32Array(
      `${TITLE.XLOADING}${idx + 1}`,
      array,
    )));

    // 2. Pack model dataframe
    const modelDfBytes = modelDf.toByteArray();
    const modelDfBytesCount = modelDfBytes.length;

    const scoresBytes = this.specn.scores.toByteArray();
    const scoresBytesCount = scoresBytes.length;

    const requiredBytes = modelDfBytesCount + scoresBytesCount + BYTES_PER_SIZES;

    const packedModel = new Uint8Array((Math.ceil(requiredBytes / BLOCK_SIZE) + 1) * BLOCK_SIZE);

    // 4 bytes for storing model's bytes count
    const sizeArr = new Uint32Array(packedModel.buffer, 0, SIZE_ARR_LEN);
    sizeArr[MODEL_IDX] = modelDfBytesCount;
    sizeArr[SCORES_IDX] = scoresBytesCount;

    // Store model's bytes
    packedModel.set(modelDfBytes, BYTES_PER_SIZES);

    // Store scores bytes
    packedModel.set(scoresBytes, BYTES_PER_SIZES + modelDfBytesCount);

    return packedModel;
  } // toBytes

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
      DG.Column.fromFloat32Array(TITLE.REGR_COEFS, this.specn.params, dim),
    ]);

    const columns = loadingsDf.columns;
    const shift = columns.length;
    const components = this.specn.components;

    this.specn.loadings.forEach((arr, idx) => loadingsDf.columns.add(
      DG.Column.fromFloat32Array(`${TITLE.XLOADING}${idx + 1}`, arr, dim),
    ));

    // Loading scatterplot
    viewers.push(DG.Viewer.scatterPlot(loadingsDf, {
      title: TITLE.LOADINGS,
      xColumnName: columns.byIndex(shift).name,
      yColumnName: columns.byIndex(shift + (components > 1 ? 1 : 0)).name,
      markerType: DG.MARKER_TYPE.CIRCLE,
      //@ts-ignore
      labelFormColumnNames: [TITLE.FEATURES],
      help: LINK.LOADINGS,
    }));

    // Regression coefficients barchart
    viewers.push(DG.Viewer.barChart(loadingsDf, {
      title: TITLE.REGR_COEFS,
      splitColumnName: TITLE.FEATURES,
      valueColumnName: TITLE.REGR_COEFS,
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

    compNames[0] = `${RESULT_NAMES.COMP} 1`;
    explVars[0] = this.specn.loadings[0][dim];

    for (let i = 1; i < components; ++i) {
      compNames[i] = `${RESULT_NAMES.COMPS} ${i + 1}`;
      explVars[i] = this.specn.loadings[i][dim];
    }

    return DG.Viewer.barChart(DG.DataFrame.fromColumns([
      DG.Column.fromStrings(RESULT_NAMES.COMPS, compNames),
      DG.Column.fromFloat32Array(TITLE.EXPL_VAR, explVars),
    ]), {
      title: TITLE.EXPL_VAR,
      splitColumnName: RESULT_NAMES.COMPS,
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
    viewers.push(
      this.explVarsViewer(),
      this.getScoresScatter(),
    );

    return viewers;
  }

  /** Return dataframe with scores */
  private getScoresDf(analysis: PlsOutput): DG.DataFrame {
    const tScores = analysis.tScores;
    const uScores = analysis.uScores;

    tScores.forEach((col, idx) => col.name = `${TITLE.XSCORE}${idx + 1}`);
    uScores.forEach((col, idx) => col.name = `${TITLE.YSCORE}${idx + 1}`);

    return DG.DataFrame.fromColumns(tScores.concat(uScores));
  }

  /** Return scores scatter */
  private getScoresScatter(): DG.Viewer {
    if (this.specn === null)
      throw new Error('Failed to create scores scatter: untrained model');

    const names = this.specn.scores.columns.names();

    const scatter = DG.Viewer.scatterPlot(this.specn.scores, {
      title: TITLE.SCORES,
      xColumnName: names[0],
      yColumnName: names[1],
      markerType: DG.MARKER_TYPE.CIRCLE,
      help: LINK.SCORES,
      showViewerFormulaLines: true,
    });

    scatter.meta.formulaLines.addAll(getLines(names));

    return scatter;
  }
};
