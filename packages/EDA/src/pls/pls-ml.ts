// Predicitve tools based on the PLS method

import * as DG from 'datagrok-api/dg';

import {TITLE, RESULT_NAMES, LINK} from './pls-constants';
import {getPlsAnalysis, getLines, PlsOutput} from './pls-tools';
import {cruxFit, runCrux} from '../crux/crux-model';
import {toArrowIpc} from '../crux/crux-bridge';
import {isNumeric, isRegressionTarget} from '../crux/methods';
import {checkColumnType, checkDimensionReducerInputs, nonNullRows} from '../utils';

/** Interactivity tresholds */
enum INTERACTIVITY {
  MAX_SAMLPES = 100000,
  MAX_FEATURES = 1000,
};

/** Applicability bound: design-matrix budget (Float32 cells, ≈ 200 MB) */
const APPLICABILITY_MAX_CELLS = 5e7;

export function isPlsApplicable(features: DG.ColumnList, target: DG.Column): boolean {
  return isNumeric(features) && isRegressionTarget(target) &&
    features.length * target.length <= APPLICABILITY_MAX_CELLS;
}

export function isPlsInteractive(features: DG.ColumnList, target: DG.Column): boolean {
  return features.length <= INTERACTIVITY.MAX_FEATURES && target.length <= INTERACTIVITY.MAX_SAMLPES;
}

/** The rows of `df` and `target` that have a target value: `ml.pls_analysis` rejects a missing y */
function targetRows(df: DG.DataFrame, target: DG.Column): {df: DG.DataFrame, target: DG.Column} {
  const missing = target.stats.missingValueCount;
  if (missing === 0)
    return {df, target};
  if (missing === target.length)
    throw new Error(`The column '${target.name}' has only missing values.`);
  const mask = nonNullRows([target]);
  return {df: df.clone(mask), target: target.clone(mask)};
}

/** Trains PLS regression; `ml.pls_analysis` checks first that the data supports `components` */
export async function trainPlsRegression(df: DG.DataFrame, target: DG.Column, components: number):
  Promise<Uint8Array> {
  const k = Math.min(components, df.columns.length);
  try {
    checkDimensionReducerInputs(df.columns, k);
    checkColumnType(target);
    const rows = targetRows(df, target);
    await runCrux('ml.pls_analysis', toArrowIpc(rows.df.columns.toList(), rows.target), {n_components: k});
  } catch (e) {
    const message = e instanceof Error ? e.message : String(e);
    throw new Error(`PLS Regression: ${message.replace(/^ml\.pls_analysis: /, '')}`);
  }
  return cruxFit('regression', 'pls', df, target, {estimator: {n_components: k}});
}

/** Loadings, regression coefficients, VIP, explained variances and scores of the PLS model of `df` */
export async function plsRegressionViewers(df: DG.DataFrame, target: DG.Column, components: number):
  Promise<DG.Viewer[]> {
  const rows = targetRows(df, target);
  const analysis = await getPlsAnalysis({table: rows.df, features: rows.df.columns, predict: rows.target, components,
    isQuadratic: false, names: undefined});
  const k = analysis.xLoadings.length;

  return [...loadingsParamsViewers(df.columns.names(), analysis, k),
    explVarsViewer(analysis.yLoadings, rows.target.length, k), scoresScatter(analysis)];
}

function loadingsParamsViewers(names: string[], analysis: PlsOutput, components: number): DG.Viewer[] {
  const coefs = analysis.regressionCoefficients;
  coefs.name = TITLE.REGR_COEFS;
  const loadingsDf = DG.DataFrame.fromColumns([DG.Column.fromStrings(TITLE.FEATURES, names), coefs]);

  const columns = loadingsDf.columns;
  for (const [idx, col] of analysis.xLoadings.entries()) {
    col.name = `${TITLE.XLOADING}${idx + 1}`;
    columns.add(col);
  }
  const vip = analysis.vip;
  vip.name = TITLE.VIP;
  columns.add(vip);

  const bars = (column: string, help: string) => DG.Viewer.barChart(loadingsDf, {
    title: column,
    splitColumnName: TITLE.FEATURES,
    valueColumnName: column,
    valueAggrType: DG.AGG.AVG,
    help,
    showValueSelector: false,
    showStackSelector: false,
  });

  return [
    DG.Viewer.scatterPlot(loadingsDf, {
      title: TITLE.LOADINGS,
      xColumnName: `${TITLE.XLOADING}1`,
      yColumnName: `${TITLE.XLOADING}${components > 1 ? 2 : 1}`,
      markerType: DG.MARKER_TYPE.CIRCLE,
      //@ts-ignore
      labelFormColumnNames: [TITLE.FEATURES],
      help: LINK.LOADINGS,
    }),
    bars(TITLE.REGR_COEFS, LINK.COEFFS),
    bars(TITLE.VIP, LINK.VIPS),
  ];
}

/** Cumulative explained variances of y, source: https://doi.org/10.1002/cem.2589 */
function explVarsViewer(yLoadings: DG.Column, samplesCount: number, components: number): DG.Viewer {
  const q = yLoadings.getRawData();
  const compNames = new Array<string>(components);
  const explVars = new Float32Array(components);

  let explVar = 0;
  for (let i = 0; i < components; ++i) {
    explVar += q[i] ** 2 / samplesCount;
    compNames[i] = `${i === 0 ? RESULT_NAMES.COMP : RESULT_NAMES.COMPS} ${i + 1}`;
    explVars[i] = explVar;
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

function scoresScatter(analysis: PlsOutput): DG.Viewer {
  for (const [idx, col] of analysis.tScores.entries())
    col.name = `${TITLE.XSCORE}${idx + 1}`;
  for (const [idx, col] of analysis.uScores.entries())
    col.name = `${TITLE.YSCORE}${idx + 1}`;
  const scores = DG.DataFrame.fromColumns(analysis.tScores.concat(analysis.uScores));
  const names = scores.columns.names();

  const scatter = DG.Viewer.scatterPlot(scores, {
    title: TITLE.SCORES,
    xColumnName: names[0],
    yColumnName: names[1],
    markerType: DG.MARKER_TYPE.CIRCLE,
    help: LINK.SCORES,
    showViewerFormulaLines: true,
  });

  scatter.meta.formulaLines.addAll(getLines(scores.columns.toList()));

  return scatter;
}
