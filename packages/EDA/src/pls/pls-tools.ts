// Tools for multivariate analysis by PLS

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PLS_ANALYSIS, ERROR_MSG, TITLE, HINT, LINK, COMPONENTS,
  INT, TIMEOUT, RESULT_NAMES, WASM_OUTPUT_IDX} from './pls-constants';
import {checkWasmDimensionReducerInputs, checkColumnType, checkMissingVals} from '../utils';
import {_partialLeastSquareRegressionInWebWorker} from '../../wasm/EDAAPI';

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
};

/** PLS analysis input */
export type PlsInput = {
  table: DG.DataFrame,
  features: DG.ColumnList,
  predict: DG.Column,
  components: number,
  names : DG.Column | null,
};

// Partial least square regression (PLS)
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
  };
}

/** Perform multivariate analysis using the PLS regression */
async function performMVA(input: PlsInput, type: PLS_ANALYSIS): Promise<void> {
  const res = await getPlsAnalysis(input);

  const plsCols = res.tScores;
  const cols = input.table.columns;
  const featuresNames = input.features.names();

  // add PLS components to the table
  plsCols.forEach((col, idx) => {
    col.name = cols.getUnusedName(`${RESULT_NAMES.PREFIX}${idx + 1}`);
    cols.add(col);
  });

  if (type === PLS_ANALYSIS.COMPUTE_COMPONENTS)
    return;

  const view = grok.shell.tableView(input.table.name);

  // 0.1 Buffer table
  const buffer = DG.DataFrame.fromColumns([
    DG.Column.fromStrings(TITLE.FEATURE, featuresNames),
    res.regressionCoefficients,
  ]);

  // 0.2. Add X-Loadings
  res.xLoadings.forEach((col, idx) => {
    col.name = buffer.columns.getUnusedName(`${TITLE.XLOADING}${idx + 1}`);
    buffer.columns.add(col);
  });

  // 1. Predicted vs Reference scatter plot
  const pred = res.prediction;
  pred.name = cols.getUnusedName(`${input.predict.name} ${RESULT_NAMES.SUFFIX}`);
  cols.add(pred);
  view.addViewer(DG.VIEWER.SCATTER_PLOT, {
    title: TITLE.MODEL,
    xColumnName: input.predict.name,
    yColumnName: pred.name,
    showRegressionLine: true,
    markerType: DG.MARKER_TYPE.CIRCLE,
    labels: input.names?.name,
    help: LINK.MODEL,
  });

  // 2. Regression Coefficients Bar Chart
  res.regressionCoefficients.name = TITLE.REGR_COEFS;
  view.addViewer(DG.Viewer.barChart(buffer, {
    title: TITLE.REGR_COEFS,
    splitColumnName: TITLE.FEATURE,
    valueColumnName: res.regressionCoefficients.name,
    valueAggrType: DG.AGG.AVG,
    help: LINK.COEFFS,
  }));

  // 3. Loading Scatter Plot
  res.xLoadings.forEach((col, idx) => col.name = `${TITLE.XLOADING}${idx + 1}`);
  view.addViewer(DG.Viewer.scatterPlot(buffer, {
    title: TITLE.LOADINGS,
    xColumnName: `${TITLE.XLOADING}1`,
    yColumnName: `${TITLE.XLOADING}${res.xLoadings.length > 1 ? '2' : '1'}`,
    markerType: DG.MARKER_TYPE.CIRCLE,
    labels: TITLE.FEATURE,
    help: LINK.LOADINGS,
  },
  ));

  // 4. Scores Scatter Plot
  plsCols.forEach((col, idx) => col.name = cols.getUnusedName(`${TITLE.XSCORE}${idx + 1}`));
  res.uScores.forEach((col, idx) => {
    col.name = cols.getUnusedName(`${TITLE.YSCORE}${idx + 1}`);
    cols.add(col);
  });
  view.addViewer(DG.VIEWER.SCATTER_PLOT, {
    title: TITLE.SCORES,
    xColumnName: plsCols[0].name,
    yColumnName: (plsCols.length > 1) ? plsCols[1].name : res.uScores[0],
    markerType: DG.MARKER_TYPE.CIRCLE,
    labels: input.names?.name,
    help: LINK.MODEL,
  });

  // 5. Explained Variances

  // 5.1) computation, source: the paper https://doi.org/10.1002/cem.2589
  //      here, we use notations from this paper
  const q = res.yLoadings.getRawData();
  const n = input.table.rowCount;
  const A = input.components;
  const explVars = new Float32Array(A);
  const compNames = [] as string[];

  explVars[0] = q[0]**2 / n;
  compNames.push(`1 ${RESULT_NAMES.COMP}`);

  for (let i = 1; i < A; ++i) {
    explVars[i] = explVars[i - 1] + q[i]**2 / n;
    compNames.push(`${i + 1} ${RESULT_NAMES.COMPS}`);
  }

  // 5.2) bar chart
  view.addViewer(DG.Viewer.barChart(DG.DataFrame.fromColumns([
    DG.Column.fromStrings(TITLE.COMPONENTS, compNames),
    DG.Column.fromFloat32Array(TITLE.EXPL_VAR, explVars),
  ]), {
    title: TITLE.EXPL_VAR,
    splitColumnName: TITLE.COMPONENTS,
    valueColumnName: TITLE.EXPL_VAR,
    valueAggrType: DG.AGG.AVG,
    help: LINK.EXPL_VARS,
  }));
} // performMVA

/** Run PLS */
export async function runMVA(type: PLS_ANALYSIS): Promise<void> {
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

  // responce (to predict)
  let predict = numCols[numCols.length - 1];
  const predictInput = ui.columnInput(TITLE.PREDICT, table, predict, () => {
    predict = predictInput.value!;
    updateIputs();
  },
  {filter: (col: DG.Column) => isValidNumeric(col)},
  );
  predictInput.setTooltip(HINT.PREDICT);

  // predictors (features)
  let features: DG.Column[];
  const featuresInput = ui.columnsInput(TITLE.USING, table, () => {}, {available: numColNames});
  featuresInput.onInput(() => updateIputs());
  featuresInput.setTooltip(HINT.FEATURES);

  // components count
  let components = min(numColNames.length - 1, COMPONENTS.DEFAULT as number);
  const componentsInput = ui.input.forProperty(DG.Property.fromOptions({
    name: TITLE.COMPONENTS,
    inputType: INT,
    defaultValue: components,
    //@ts-ignore
    showPlusMinus: true,
    min: COMPONENTS.MIN,
  }));
  componentsInput.onInput(() => updateIputs());
  componentsInput.setTooltip(HINT.COMPONENTS);

  let dlgTitle: string;
  let dlgHelpUrl: string;
  let dlgRunBtnTooltip: string;

  if (type === PLS_ANALYSIS.COMPUTE_COMPONENTS) {
    dlgTitle = TITLE.PLS;
    dlgHelpUrl = LINK.PLS;
    dlgRunBtnTooltip = HINT.PLS;
  } else {
    dlgTitle = TITLE.MVA;
    dlgHelpUrl = LINK.MVA;
    dlgRunBtnTooltip = HINT.MVA;
  }

  const updateIputs = () => {
    featuresInput.value = featuresInput.value.filter((col) => col !== predict);
    features = featuresInput.value;

    componentsInput.value = min(max(componentsInput.value ?? components, COMPONENTS.MIN), features.length);
    components = componentsInput.value;

    dlg.getButton(TITLE.RUN).disabled = (features.length === 0) || (components <= 0);
  };

  // names of samples
  let names = (strCols.length > 0) ? strCols[0] : null;
  const namesInputs = ui.columnInput(TITLE.NAMES, table, names, () => names = predictInput.value,
    {filter: (col: DG.Column) => col.type === DG.COLUMN_TYPE.STRING},
  );
  namesInputs.setTooltip(HINT.NAMES);
  namesInputs.root.hidden = (strCols.length === 0) || (type === PLS_ANALYSIS.COMPUTE_COMPONENTS);

  const dlg = ui.dialog({title: dlgTitle, helpUrl: dlgHelpUrl})
    .add(ui.form([predictInput, featuresInput, componentsInput, namesInputs]))
    .addButton(TITLE.RUN, async () => {
      dlg.close();

      await performMVA({
        table: table,
        features: DG.DataFrame.fromColumns(features).columns,
        predict: predict,
        components: components,
        names: names,
      }, type);
    }, undefined, dlgRunBtnTooltip)
    .show();

  // the following delay provides correct styles (see https://reddata.atlassian.net/browse/GROK-15196)
  setTimeout(() => {
    featuresInput.value = numCols.filter((col) => col !== predict);
    features = featuresInput.value;
  }, TIMEOUT);

  grok.shell.v.append(dlg.root);
} // addPLS
