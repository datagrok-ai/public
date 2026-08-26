// Tools for multivariate analysis by PLS

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as jStat from 'jstat';

import {PLS_ANALYSIS, ERROR_MSG, TITLE, HINT, LINK, COMPONENTS,
  RESULT_NAMES, WASM_OUTPUT_IDX, ELLIPSES, LINE_WIDTH, COLOR, X_COORD, Y_COORD,
  DEMO_INTRO_MD, DEMO_RESULTS, NUMS_AFTER_COMMA,
  MAX_ROWS_IN_PREDICTION_TOOLTIP, DELAY,
  MVA_MODEL_TAG, MVA_TRANSFORM_FUNC, MVA_INIT_FUNC} from './pls-constants';
import {checkWasmDimensionReducerInputs, checkColumnType, checkMissingVals, describeElements} from '../utils';
// PLS1 migrated to Rust + WASM (sci-comp-ml).
import {_partialLeastSquareRegressionInWebWorker} from '../../wasm/eda-api';
import {carsDataframe} from '../data-generators';

const min = Math.min;
const max = Math.max;

/** PLS analysis results */
export type PlsOutput = {
  prediction: DG.Column<DG.COLUMN_TYPE.FLOAT>,
  regressionCoefficients: DG.Column<DG.COLUMN_TYPE.FLOAT>,
  tScores: DG.Column<DG.COLUMN_TYPE.FLOAT>[],
  uScores: DG.Column<DG.COLUMN_TYPE.FLOAT>[],
  xLoadings: DG.Column<DG.COLUMN_TYPE.FLOAT>[],
  yLoadings: DG.Column<DG.COLUMN_TYPE.FLOAT>,
  vip: DG.Column<DG.COLUMN_TYPE.FLOAT>,
};

/** PLS analysis input */
export type PlsInput = {
  table: DG.DataFrame,
  features: DG.ColumnList,
  predict: DG.Column,
  components: number,
  isQuadratic: boolean,
  names : DG.Column | undefined,
};

/** Multivariate analysis input: features are columns of the source table, so that
 * the transform function call can be replayed against the restored table */
export type MvaInput = Omit<PlsInput, 'features'> & {features: DG.Column[]};

/** Names of the columns & tables that the multivariate analysis produces */
export type MvaNames = {
  xScores: string[],
  yScores: string[],
  prediction: string,
  analysisTable: string,
  explVarTable: string,
};

/** The fitted model: enough to restore the prediction tooltip once a project is reopened */
type MvaModel = {
  features: string[],
  coefficients: number[],
  bias: number,
  prediction: string,
};

type TypedArray = Int32Array | Float32Array | Uint32Array | Float64Array;

/** Set style for input element depending on the validity of the value */
function setStyle(valid: boolean, element: HTMLElement, tooltip: string, errorMsg: string) {
  if (valid) {
    element.style.color = COLOR.VALID_TEXT;
    element.style.borderBottomColor = COLOR.VALID_LINE;
    ui.tooltip.bind(element, tooltip);
  } else {
    element.style.color = COLOR.INVALID;
    element.style.borderBottomColor = COLOR.INVALID;
    ui.tooltip.bind(element, () => {
      const hint = ui.label(tooltip);
      const err = ui.label(errorMsg);
      err.style.color = COLOR.INVALID;
      return ui.divV([hint, err]);
    });
  }
};

function getModelFormulaTerms(model: MvaModel): Map<string, number> {
  const terms = new Map([[TITLE.BIAS as string, model.bias]]);

  model.features.forEach((name, idx) => {
    terms.set(name, model.coefficients[idx]);
  });

  return terms;
}

/** Hotelling's T² limit for a 2D score plot at the given confidence (p = 2). */
export function hotellingT2Limit(conf: number, n: number): number {
  const p = 2;
  return (p * (n - 1) / (n - p)) * jStat.centralF.inv(conf, p, n - p);
}

/** Axis-aligned Hotelling's T² ellipse for a score pair; DG variance is n-1, like R's var(). */
export function hotellingEllipseParams(xCol: DG.Column, yCol: DG.Column, conf: number):
  {cx: number, cy: number, a: number, b: number} {
  const sx = xCol.stats;
  const sy = yCol.stats;
  const t2 = hotellingT2Limit(conf, xCol.length);
  return {cx: sx.avg, cy: sy.avg, a: Math.sqrt(t2 * sx.variance), b: Math.sqrt(t2 * sy.variance)};
}

/** Decimal string for a formula-line expression — the parser rejects exponential notation. */
function fmt(v: number): string {
  return v.toLocaleString('en-US', {useGrouping: false, maximumFractionDigits: 20});
}

/** Scores-plot lines: an x=0 axis per component and the 95%/99% Hotelling's T² ellipses. */
export function getLines(cols: DG.Column[]): DG.FormulaLine[] {
  const lines: DG.FormulaLine[] = [];
  const n = cols.length > 0 ? cols[0].length : 0;
  const stats = new Map(cols.map((c) => [c.name, c.stats]));
  // T-squared is defined only for n - p > 0 (p = 2)
  const limits = n > 2 ? ELLIPSES.map((e) => ({...e, t2: hotellingT2Limit(e.conf, n)})) : [];

  const ref = (name: string) => '${' + name + '}';

  const addEllipseHalf = (formula: string, minX: number, maxX: number, color: string) =>
    lines.push({type: 'line', formula, width: LINE_WIDTH, visible: true, title: ' ', min: minX, max: maxX, color});

  cols.forEach((xCol) => {
    const x = ref(xCol.name);
    lines.push({type: 'line', formula: `${x} = 0`, width: LINE_WIDTH, visible: true, title: ' ', color: COLOR.AXIS});

    cols.forEach((yCol) => {
      if (yCol.name === xCol.name)
        return;

      const y = ref(yCol.name);
      const sx = stats.get(xCol.name)!;
      const sy = stats.get(yCol.name)!;
      if (sx.variance <= 0 || sy.variance <= 0)
        return;

      // both orderings emitted: the y=f(x) and x=f(y) arcs fill each other's gaps
      for (const {t2, color} of limits) {
        const a = Math.sqrt(t2 * sx.variance);
        const b = Math.sqrt(t2 * sy.variance);
        const cx = sx.avg;
        const cy = sy.avg;
        // (x - cx)²/a² + (y - cy)²/b² = 1  ->  y = cy ± b·sqrt(1 - (x - cx)²/a²)
        const dx = `(${x} - ${fmt(cx)})`;
        const inner = `sqrt(1 - ${dx} * ${dx} / ${fmt(a * a)})`;
        addEllipseHalf(`${y} = ${fmt(cy)} + ${fmt(b)} * ${inner}`, cx - a, cx + a, color);
        addEllipseHalf(`${y} = ${fmt(cy)} - ${fmt(b)} * ${inner}`, cx - a, cx + a, color);
      }
    });
  });

  return lines;
}

/** Partial least square regression (PLS) */
export async function getPlsAnalysis(input: PlsInput): Promise<PlsOutput> {
  checkWasmDimensionReducerInputs(input.features, input.components);

  // Check the responce column
  checkColumnType(input.predict);
  checkMissingVals(input.predict);

  const result = await _partialLeastSquareRegressionInWebWorker(
    input.table,
    input.features,
    input.predict,
    input.components,
  );

  return {
    prediction: result[WASM_OUTPUT_IDX.PREDICTION],
    regressionCoefficients: result[WASM_OUTPUT_IDX.REGR_COEFFS],
    tScores: result[WASM_OUTPUT_IDX.T_SCORES],
    uScores: result[WASM_OUTPUT_IDX.U_SCORES],
    xLoadings: result[WASM_OUTPUT_IDX.X_LOADINGS],
    yLoadings: result[WASM_OUTPUT_IDX.Y_LOADINGS],
    vip: result[WASM_OUTPUT_IDX.VIP],
  };
}

/** Return debiased predction by PLS regression */
function debiasedPrediction(features: DG.ColumnList, params: DG.Column,
  target: DG.Column, biasedPrediction: DG.Column): {debiased: DG.Column, bias: number} {
  const samples = target.length;
  const dim = features.length;
  const rawParams = params.getRawData();
  const debiased = new Float32Array(samples);
  const biased = biasedPrediction.getRawData();

  // Compute bias
  let bias = target.stats.avg;
  for (let i = 0; i < dim; ++i)
    bias -= rawParams[i] * features.byIndex(i).stats.avg;

  // Compute debiased prediction
  for (let i = 0; i < samples; ++i)
    debiased[i] = bias + biased[i];

  return {debiased: DG.Column.fromFloat32Array('Debiased', debiased, samples), bias: bias};
}

/** Return an input for the quadratic PLS regression */
function getQuadraticPlsInput(input: PlsInput): PlsInput {
  if (!input.isQuadratic)
    return input;

  const cols: DG.Column[] = input.features.toList();
  const colsCount = cols.length;
  const rowsCount = input.table.rowCount;
  const quadrCols: DG.Column[] = [];
  let col1: DG.Column;
  let raw1: TypedArray;
  let col2: DG.Column;
  let raw2: TypedArray;
  let qaudrRaw: Float32Array;

  for (let i = 0; i < colsCount; ++i) {
    col1 = cols[i];
    raw1 = col1.getRawData();

    for (let j = i; j < colsCount; ++j) {
      col2 = cols[j];
      raw2 = col2.getRawData();
      qaudrRaw = new Float32Array(rowsCount);

      for (let k = 0; k < rowsCount; ++k)
        qaudrRaw[k] = raw1[k] * raw2[k];

      const quadrCol = DG.Column.fromFloat32Array(`${col1.name} x ${col2.name}`, qaudrRaw);

      if (quadrCol.stats.stdev > 0)
        quadrCols.push(quadrCol);
    }
  }

  const extendedTable = DG.DataFrame.fromColumns(cols.concat(quadrCols));

  return {
    table: extendedTable,
    features: extendedTable.columns,
    isQuadratic: true,
    names: input.names,
    predict: input.predict,
    components: input.components,
  };
}

/** Return a table name that is not yet used in the workspace */
function getUnusedTableName(base: string): string {
  const names = grok.shell.tableNames;
  let name = base;
  let idx = 1;

  while (names.includes(name))
    name = `${base} (${++idx})`;

  return name;
}

/** Names of the columns & tables that the multivariate analysis produces. They are reserved
 * before the fit and passed to the transform function explicitly: replaying the table creation
 * script must reproduce exactly the names that the saved layout refers to. */
export function getMvaNames(input: MvaInput, componentsOnly: boolean): MvaNames {
  const cols = input.table.columns;
  const comps = [...Array(input.components).keys()];
  const xScores = comps.map((idx) =>
    cols.getUnusedName(`${componentsOnly ? RESULT_NAMES.PREFIX : TITLE.XSCORE}${idx + 1}`));

  if (componentsOnly)
    return {xScores: xScores, yScores: [], prediction: '', analysisTable: '', explVarTable: ''};

  return {
    xScores: xScores,
    yScores: comps.map((idx) => cols.getUnusedName(`${TITLE.YSCORE}${idx + 1}`)),
    prediction: cols.getUnusedName(`${input.predict.name} ${RESULT_NAMES.SUFFIX}`),
    analysisTable: getUnusedTableName(`${input.table.name}(${TITLE.ANALYSIS})`),
    explVarTable: getUnusedTableName(`${input.table.name}(${TITLE.EXPL_VAR})`),
  };
}

/** Add the multivariate analysis results: PLS components, prediction & the analysis tables */
export async function addMvaResults(input: PlsInput, names: MvaNames, componentsOnly: boolean): Promise<void> {
  const sourceTable = input.table;
  const plsInput = input.isQuadratic ? getQuadraticPlsInput(input) : input;
  const result = await getPlsAnalysis(plsInput);
  const cols = sourceTable.columns;

  // add PLS components to the table
  result.tScores.forEach((col, idx) => {
    col.name = cols.getUnusedName(names.xScores[idx]);
    cols.add(col);
  });

  if (componentsOnly)
    return;

  const features = plsInput.features;
  const featuresNames = features.names();

  // debias prediction (since PLS centers data)
  const debiased = debiasedPrediction(features, result.regressionCoefficients, input.predict, result.prediction);
  const pred = debiased.debiased;
  pred.name = cols.getUnusedName(names.prediction);
  cols.add(pred);

  result.uScores.forEach((col, idx) => {
    col.name = cols.getUnusedName(names.yScores[idx]);
    cols.add(col);
  });

  // features analysis table: regression coefficients, X-loadings & VIP
  result.regressionCoefficients.name = TITLE.REGR_COEFS;
  const loadingsRegrCoefsTable = DG.DataFrame.fromColumns([
    DG.Column.fromStrings(TITLE.FEATURE, featuresNames),
    result.regressionCoefficients,
  ]);

  loadingsRegrCoefsTable.name = names.analysisTable;
  grok.shell.addTable(loadingsRegrCoefsTable);

  result.xLoadings.forEach((col, idx) => {
    col.name = loadingsRegrCoefsTable.columns.getUnusedName(`${TITLE.XLOADING}${idx + 1}`);
    loadingsRegrCoefsTable.columns.add(col);
  });

  result.vip.name = loadingsRegrCoefsTable.columns.getUnusedName(TITLE.VIP);
  loadingsRegrCoefsTable.columns.add(result.vip);

  // explained variances, source: the paper https://doi.org/10.1002/cem.2589
  // here, we use notations from this paper
  const q = result.yLoadings.getRawData();
  const p = result.xLoadings.map((col) => col.getRawData());
  const n = sourceTable.rowCount;
  const m = featuresNames.length;
  const A = input.components;
  const yExplVars = new Float32Array(A);
  const compNames = [] as string[];
  const xExplVars: Float32Array[] = [];
  for (let i = 0; i < m; ++i)
    xExplVars.push(new Float32Array(A));

  yExplVars[0] = q[0]**2 / n;
  compNames.push(`1 ${RESULT_NAMES.COMP}`);
  xExplVars.forEach((arr, idx) => {arr[0] = p[0][idx]**2 / n;});

  for (let comp = 1; comp < A; ++comp) {
    yExplVars[comp] = yExplVars[comp - 1] + q[comp]**2 / n;
    xExplVars.forEach((arr, idx) => arr[comp] = arr[comp - 1] + p[comp][idx]**2 / n);
    compNames.push(`${comp + 1} ${RESULT_NAMES.COMPS}`);
  }

  const explVarsDF = DG.DataFrame.fromColumns([
    DG.Column.fromStrings(TITLE.COMPONENTS, compNames),
    DG.Column.fromFloat32Array(input.predict.name, yExplVars),
  ]);

  explVarsDF.name = names.explVarTable;
  grok.shell.addTable(explVarsDF);

  xExplVars.forEach((arr, idx) => explVarsDF.columns.add(DG.Column.fromFloat32Array(featuresNames[idx], arr)));

  sourceTable.setTag(MVA_MODEL_TAG, JSON.stringify({
    features: featuresNames,
    coefficients: result.regressionCoefficients.toList(),
    bias: debiased.bias,
    prediction: pred.name,
  } as MvaModel));
} // addMvaResults

/** Add the multivariate analysis viewers to the table view */
function addMvaViewers(input: MvaInput, names: MvaNames, analysisType: PLS_ANALYSIS): void {
  const sourceTable = input.table;
  const view = grok.shell.tableView(sourceTable.name);
  const loadingsRegrCoefsTable = grok.shell.tableByName(names.analysisTable);
  const explVarsDF = grok.shell.tableByName(names.explVarTable);
  const model: MvaModel = JSON.parse(sourceTable.getTag(MVA_MODEL_TAG)!);

  // 1. Predicted vs Reference scatter plot
  const predictVsReferScatter = view.addViewer(DG.Viewer.scatterPlot(sourceTable, {
    title: TITLE.MODEL,
    xColumnName: input.predict.name,
    yColumnName: model.prediction,
    showRegressionLine: true,
    markerType: DG.MARKER_TYPE.CIRCLE,
    showLabels: 'Always',
    help: LINK.MODEL,
    initializationFunction: MVA_INIT_FUNC,
  }));

  if ((input.names !== undefined) && (input.names !== null))
    predictVsReferScatter.setOptions({labelColumnNames: [input.names?.name]});

  // 2. Regression Coefficients Bar Chart
  const regrCoeffsBar = view.addViewer(DG.Viewer.barChart(loadingsRegrCoefsTable, {
    table: loadingsRegrCoefsTable.name,
    title: TITLE.REGR_COEFS,
    splitColumnName: TITLE.FEATURE,
    valueColumnName: TITLE.REGR_COEFS,
    valueAggrType: DG.AGG.AVG,
    help: LINK.COEFFS,
    showValueSelector: false,
    showStackSelector: false,
    description: `${TITLE.BIAS} = ${model.bias.toFixed(NUMS_AFTER_COMMA)}`,
    descriptionVisibilityMode: 'Always',
    descriptionPosition: 'Bottom',
  }));

  // 3. Loadings Scatter Plot
  const loadingsScatter = view.addViewer(DG.Viewer.scatterPlot(loadingsRegrCoefsTable, {
    table: loadingsRegrCoefsTable.name,
    title: TITLE.LOADINGS,
    xColumnName: `${TITLE.XLOADING}1`,
    yColumnName: `${TITLE.XLOADING}${names.xScores.length > 1 ? '2' : '1'}`,
    markerType: DG.MARKER_TYPE.CIRCLE,
    labelColumnNames: [TITLE.FEATURE],
    help: LINK.LOADINGS,
  }));

  // 4. Scores Scatter Plot
  const scoreNames = names.xScores.concat(names.yScores);
  const scoresScatter = DG.Viewer.scatterPlot(sourceTable, {
    title: TITLE.SCORES,
    xColumnName: names.xScores[0],
    yColumnName: (names.xScores.length > 1) ? names.xScores[1] : names.yScores[0],
    markerType: DG.MARKER_TYPE.CIRCLE,
    help: LINK.SCORES,
    showViewerFormulaLines: true,
    labelColumnNames: ((input.names !== undefined) && (input.names !== null)) ? [input.names?.name] : undefined,
  });

  // create axes & Hotelling's T-squared confidence ellipses
  view.addViewer(scoresScatter);
  const scoreCols = scoreNames.map((name) => sourceTable.col(name)).filter((c): c is DG.Column => c !== null);
  scoresScatter.meta.formulaLines.addAll(getLines(scoreCols));

  // 5. Explained Variance Bar Chart
  const explVarsBar = view.addViewer(DG.Viewer.barChart(explVarsDF, {
    table: explVarsDF.name,
    title: TITLE.EXPL_VAR,
    splitColumnName: TITLE.COMPONENTS,
    valueColumnName: input.predict.name,
    valueAggrType: DG.AGG.AVG,
    help: LINK.EXPL_VARS,
    showCategorySelector: false,
    showStackSelector: false,
  }));

  // 6. Variable Importance Bar Chart
  const vipViewer = DG.Viewer.barChart(loadingsRegrCoefsTable, {
    table: loadingsRegrCoefsTable.name,
    title: TITLE.VIP,
    splitColumnName: TITLE.FEATURE,
    valueColumnName: TITLE.VIP,
    valueAggrType: DG.AGG.AVG,
    help: LINK.VIPS,
    showValueSelector: false,
    showStackSelector: false,
  });

  const regrCoeffsNode = view.dockManager.findNode(regrCoeffsBar.root);

  if (regrCoeffsNode !== null) {
    view.dockManager.dock(
      vipViewer,
      DG.DOCK_TYPE.FILL,
      regrCoeffsNode,
    );
  }

  // emphasize viewers in the demo case
  if (analysisType === PLS_ANALYSIS.DEMO) {
    setTimeout(() => {
      describeElements(
        [predictVsReferScatter, scoresScatter, loadingsScatter, vipViewer, explVarsBar].map((v) => v.root),
        DEMO_RESULTS.map((info) => `<b>${info.caption}</b>\n\n${info.text}`),
        ['left', 'left', 'right', 'right', 'left'],
        view.root,
      );
    }, DELAY);
  }
} // addMvaViewers

/** Show the model formula in the prediction column tooltip. Called by the platform each time
 * the model viewer is created, including its restoration from a project or layout. */
export function initMvaModelViewer(viewer: DG.Viewer): void {
  const table = viewer.dataFrame;
  const modelJson = table?.getTag(MVA_MODEL_TAG);

  if (!modelJson)
    return;

  const model: MvaModel = JSON.parse(modelJson);
  const view = grok.shell.tableView(table.name);
  const pred = table.col(model.prediction);

  if ((view === null) || (pred === null))
    return;

  setPredictionTooltip(view, pred, getModelFormulaTerms(model));
}

/** Add the analysis results by calling the transform function: the call goes to the table creation
 * script, so a data-synced project restores the columns that the saved viewers refer to */
export async function callMvaTransform(input: MvaInput, names: MvaNames, componentsOnly: boolean): Promise<void> {
  await DG.Func.find({package: 'EDA', name: MVA_TRANSFORM_FUNC})[0].prepare({
    table: input.table,
    featureNames: input.features.map((col) => col.name),
    predict: input.predict,
    components: input.components,
    isQuadratic: input.isQuadratic,
    componentsOnly: componentsOnly,
    xScoreNames: names.xScores,
    yScoreNames: names.yScores,
    predictionName: names.prediction,
    analysisTableName: names.analysisTable,
    explVarTableName: names.explVarTable,
  }).call(undefined, undefined, {processed: false});
}

/** Perform multivariate analysis using the PLS regression */
async function performMVA(input: MvaInput, analysisType: PLS_ANALYSIS): Promise<void> {
  const componentsOnly = (analysisType === PLS_ANALYSIS.COMPUTE_COMPONENTS);
  const names = getMvaNames(input, componentsOnly);

  await callMvaTransform(input, names, componentsOnly);

  if (!componentsOnly)
    addMvaViewers(input, names, analysisType);
} // performMVA

/** Run multivariate analysis (PLS) */
export async function runMVA(analysisType: PLS_ANALYSIS): Promise<void> {
  const table = grok.shell.t;

  if (table === null) {
    grok.shell.warning(ERROR_MSG.NO_DF);
    return;
  }

  if (table.rowCount === 0) {
    grok.shell.warning(ERROR_MSG.EMPTY_DF);
    return;
  }

  const numColNames = [] as string[];
  const numCols = [] as DG.Column[];
  const strCols = [] as DG.Column[];

  const isValidNumeric = (col: DG.Column) =>
    ((col.type === DG.COLUMN_TYPE.INT) || (col.type === DG.COLUMN_TYPE.FLOAT)) &&
        (col.stats.missingValueCount === 0);

  table.columns.toList().forEach((col) => {
    if (isValidNumeric(col)) {
      numColNames.push(col.name);
      numCols.push(col);
    } else if (col.type === DG.COLUMN_TYPE.STRING)
      strCols.push(col);
  });

  if (numColNames.length === 0) {
    grok.shell.warning(ERROR_MSG.NO_COLS);
    return;
  }

  if (numColNames.length === 1) {
    grok.shell.warning(ERROR_MSG.ONE_COL);
    return;
  }

  const doFeaturesIncludePredict = () => {
    return featuresInput.value.some((col) => col.name === predictInput.value!.name);
  };

  const isCompConsistent = () => {
    if (componentsInput.value! < 1)
      return false;

    if (componentsInput.value! > table.rowCount)
      return false;

    const n = featuresInput.value.length;

    if (isQuadraticInput.value)
      return componentsInput.value! <= (n + 1) * n / 2 + n;

    return componentsInput.value! <= n;
  };

  // response (to predict)
  const predictInput = ui.input.column(TITLE.PREDICT, {
    table: table,
    value: numCols[numCols.length - 1],
    nullable: false,
    onValueChanged: (_) => {
      updateInputStyles();
    },
    filter: (col: DG.Column) => isValidNumeric(col),
    tooltipText: HINT.PREDICT,
  });

  // predictors (features)
  const featuresInput = ui.input.columns(TITLE.USING, {
    table: table,
    available: numColNames,
    value: numCols.slice(0, numCols.length - 1),
    onValueChanged: (_) => {
      updateInputStyles();
    },
    tooltipText: HINT.FEATURES,
    nullable: false,
  });

  // components count
  const componentsInput = ui.input.int(TITLE.COMPONENTS, {
    value: min(numColNames.length - 1, COMPONENTS.DEFAULT as number),
    showPlusMinus: true,
    nullable: false,
    onValueChanged: (_) => {
      updateInputStyles();
    },
    tooltipText: HINT.COMPONENTS,
  });

  let dlgTitle: string;
  let dlgHelpUrl: string;
  let dlgRunBtnTooltip: string;

  if (analysisType === PLS_ANALYSIS.COMPUTE_COMPONENTS) {
    dlgTitle = TITLE.PLS;
    dlgHelpUrl = LINK.PLS;
    dlgRunBtnTooltip = HINT.PLS;
  } else {
    dlgTitle = TITLE.MVA;
    dlgHelpUrl = LINK.MVA;
    dlgRunBtnTooltip = HINT.MVA;
  }

  const updateInputStyles = () => {
    const featuresValid = featuresInput.value.length >= 1;
    const predValid = featuresValid && !doFeaturesIncludePredict();
    let compValid: boolean;

    if (!featuresValid)
      setStyle(false, featuresInput.input, HINT.FEATURES, ERROR_MSG.ENOUGH);
    else if (predValid) {
      setStyle(true, predictInput.input, HINT.PREDICT, '');
      setStyle(true, featuresInput.input, HINT.FEATURES, '');
    } else {
      setStyle(false, predictInput.input, HINT.PREDICT, ERROR_MSG.PREDICT);
      setStyle(false, featuresInput.input, HINT.FEATURES, ERROR_MSG.PREDICT);
    }

    if (componentsInput.value == null) {
      setStyle(false, componentsInput.input, HINT.COMPONENTS, ERROR_MSG.NULL_COMPS);
      compValid = false;
    } else if (componentsInput.value < 1) {
      setStyle(false, componentsInput.input, HINT.COMPONENTS, ERROR_MSG.COMPONENTS);
      compValid = false;
    } else {
      compValid = isCompConsistent();

      if (compValid) {
        setStyle(true, componentsInput.input, HINT.COMPONENTS, '');
        if (predValid)
          setStyle(true, featuresInput.input, HINT.FEATURES, '');
      } else {
        const errMsg = componentsInput.value! > table.rowCount ?
          ERROR_MSG.COMP_ROWS :
          isQuadraticInput.value ? ERROR_MSG.COMP_QUA_PLS : ERROR_MSG.COMP_LIN_PLS;
        setStyle(false, componentsInput.input, HINT.COMPONENTS, errMsg);
        setStyle(false, featuresInput.input, HINT.FEATURES, ERROR_MSG.ENOUGH);
      }
    }

    const isValid = predValid && compValid;

    dlg.getButton(TITLE.RUN).disabled = !isValid;

    return isValid;
  }; // updateInputStyles

  const getStrColWithUniqueVals = () => {
    for (const col of strCols) {
      if (col.stats.uniqueCount === table.rowCount)
        return col;
    }
    return undefined;
  };

  // names of samples
  let names = getStrColWithUniqueVals();
  const namesInputs = ui.input.column(TITLE.NAMES, {
    table: table,
    value: names,
    onValueChanged: (value) => names = value ?? undefined,
    filter: (col: DG.Column) => col.type === DG.COLUMN_TYPE.STRING},
  );
  namesInputs.setTooltip(HINT.NAMES);
  namesInputs.root.hidden = (strCols.length === 0) || (analysisType === PLS_ANALYSIS.COMPUTE_COMPONENTS);

  // quadratic/linear model
  const isQuadraticInput = ui.input.bool(TITLE.QUADRATIC, {
    value: false,
    tooltipText: HINT.QUADRATIC,
    onValueChanged: (_) => {
      updateInputStyles();
    },
  });

  const dlg = ui.dialog({title: dlgTitle, helpUrl: dlgHelpUrl})
    .add(ui.form([predictInput, featuresInput, componentsInput, isQuadraticInput, namesInputs]))
    .addButton(TITLE.RUN, async () => {
      dlg.close();

      await performMVA({
        table: table,
        features: featuresInput.value,
        predict: predictInput.value!,
        components: componentsInput.value!,
        isQuadratic: isQuadraticInput.value,
        names: names,
      }, analysisType);
    }, undefined, dlgRunBtnTooltip)
    .show({x: X_COORD, y: Y_COORD});

  grok.shell.v.append(dlg.root);
} // runMVA

/** Run multivariate analysis demo */
export async function runDemoMVA(): Promise<void> {
  const table = carsDataframe();
  grok.shell.addTableView(table);
  grok.shell.windows.help.visible = true;
  grok.shell.windows.help.showHelp(ui.markdown(DEMO_INTRO_MD));
  grok.shell.windows.showContextPanel = false;
  grok.shell.windows.showProperties = false;

  const cols = table.columns.toList();
  const numCols = cols.filter((col) =>
    ((col.type === DG.COLUMN_TYPE.INT) || (col.type === DG.COLUMN_TYPE.FLOAT)) &&
    (col.stats.missingValueCount === 0));
  const names = cols.find((col) =>
    (col.type === DG.COLUMN_TYPE.STRING) && (col.stats.uniqueCount === table.rowCount));

  await performMVA({
    table: table,
    features: numCols.slice(0, -1),
    predict: numCols[numCols.length - 1],
    components: min(numCols.length - 1, COMPONENTS.DEFAULT as number),
    isQuadratic: false,
    names: names,
  }, PLS_ANALYSIS.DEMO);
}

function setPredictionTooltip(view: DG.TableView, predCol: DG.Column, modelTerms: Map<string, number>): void {
  view.grid.onCellTooltip((cell, x, y) => {
    if (cell.isColHeader) {
      const cellCol = cell.tableColumn;

      if (cellCol == null)
        return false;

      if (cellCol.name === predCol.name) {
        ui.tooltip.show(getPredictionTooltip(modelTerms, predCol), x, y);
        return true;
      }
    }
    return false;
  });
}

function getPredictionTooltip(modelTerms: Map<string, number>, predCol: DG.Column): HTMLElement {
  let idx = 0;
  const bias = modelTerms.get(TITLE.BIAS) ?? 0;
  const elements: HTMLElement[] = [];
  if (Math.abs(bias) > 0) {
    const biasEl = ui.divText(`${bias}`);
    biasEl.style.marginTop = '2px';
    biasEl.style.marginLeft = '4px';
    elements.push(biasEl);
    ++idx;
  }

  const sortedTerms = [...modelTerms.entries()]
    .filter(([key]) => key !== TITLE.BIAS)
    .sort((a, b) => Math.abs(b[1]) - Math.abs(a[1]));

  const maxFeatureRows = MAX_ROWS_IN_PREDICTION_TOOLTIP - elements.length;
  const hasOverflow = sortedTerms.length > maxFeatureRows;
  const visibleTerms = hasOverflow ? sortedTerms.slice(0, maxFeatureRows - 1) : sortedTerms;

  for (const [key, value] of visibleTerms) {
    const signEl = ui.divText(idx > 0 ? '+ ' : '');
    signEl.style.marginRight = '4px';
    signEl.style.marginLeft = '4px';

    const featureEl = ui.divText(`${key}`);
    featureEl.style.fontWeight = 'bold';

    const valueEl = ui.divText(` * ${value > 0 ? value : `(${value})`}`);
    valueEl.style.marginLeft = '4px';

    const rowEl = ui.divH([signEl, featureEl, valueEl]);
    rowEl.style.marginTop = '4px';
    elements.push(rowEl);

    ++idx;
  }

  if (hasOverflow) {
    const hidden = sortedTerms.length - visibleTerms.length;
    const ellipsisEl = ui.divText(`(${hidden} more term${hidden > 1 ? 's' : ''})`);
    ellipsisEl.style.marginTop = '4px';
    ellipsisEl.style.marginLeft = '4px';
    ellipsisEl.style.fontStyle = 'italic';
    elements.push(ellipsisEl);
  }

  const headerEl = ui.divText('Formula:');

  const leftEl = ui.divText(`${predCol.name} = `);
  leftEl.style.fontWeight = 'bold';
  leftEl.style.marginTop = '4px';

  const elementsContainer = ui.divV(elements);
  elementsContainer.style.marginTop = '4px';

  return ui.divV([headerEl, leftEl, elementsContainer]);
}
