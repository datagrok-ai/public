// Tools for multivariate analysis by PLS

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PLS_ANALYSIS, ERROR_MSG, TITLE, HINT, LINK, COMPONENTS, INT, TIMEOUT,
  RESULT_NAMES, WASM_OUTPUT_IDX, RADIUS, LINE_WIDTH, COLOR, X_COORD, Y_COORD,
  DEMO_INTRO_MD, DEMO_RESULTS_MD, DEMO_RESULTS} from './pls-constants';
import {checkWasmDimensionReducerInputs, checkColumnType, checkMissingVals, describeElements} from '../utils';
import {_partialLeastSquareRegressionInWebWorker} from '../../wasm/EDAAPI';
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

type TypedArray = Int32Array | Float32Array | Uint32Array | Float64Array;

/** Return lines */
export function getLines(names: string[]): DG.FormulaLine[] {
  const lines: DG.FormulaLine[] = [];

  const addLine = (formula: string, radius: number) => {
    lines.push({
      type: 'line',
      formula: formula,
      width: LINE_WIDTH,
      visible: true,
      title: ' ',
      min: -radius,
      max: radius,
      color: COLOR.CIRCLE,
    });
  };

  names.forEach((xName) => {
    const x = '${' + xName + '}';
    lines.push({type: 'line', formula: `${x} = 0`, width: LINE_WIDTH, visible: true, title: ' ', color: COLOR.AXIS});

    names.forEach((yName) => {
      const y = '${' + yName + '}';

      RADIUS.forEach((r) => {
        addLine(y + ` = sqrt(${r*r} - ${x} * ${x})`, r);
        addLine(y + ` = -sqrt(${r*r} - ${x} * ${x})`, r);
      });
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
  };
}

/** Return debiased predction by PLS regression */
function debiasedPrediction(features: DG.ColumnList, params: DG.Column,
  target: DG.Column, biasedPrediction: DG.Column): DG.Column {
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

  return DG.Column.fromFloat32Array('Debiased', debiased, samples);
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

/** Perform multivariate analysis using the PLS regression */
async function performMVA(input: PlsInput, analysisType: PLS_ANALYSIS): Promise<void> {
  const sourceTable = input.table;

  if (input.isQuadratic)
    input = getQuadraticPlsInput(input);

  const result = await getPlsAnalysis(input);

  const plsCols = result.tScores;
  const cols = sourceTable.columns;
  const features = input.features;
  const featuresNames = features.names();
  const prefix = (analysisType === PLS_ANALYSIS.COMPUTE_COMPONENTS) ? RESULT_NAMES.PREFIX : TITLE.XSCORE;

  // add PLS components to the table
  plsCols.forEach((col, idx) => {
    col.name = cols.getUnusedName(`${prefix}${idx + 1}`);
    cols.add(col);
  });

  if (analysisType === PLS_ANALYSIS.COMPUTE_COMPONENTS)
    return;

  const view = grok.shell.tableView(sourceTable.name);

  // 0.1 Buffer table
  const loadingsRegrCoefsTable = DG.DataFrame.fromColumns([
    DG.Column.fromStrings(TITLE.FEATURE, featuresNames),
    result.regressionCoefficients,
  ]);

  loadingsRegrCoefsTable.name = `${sourceTable.name}(${TITLE.ANALYSIS})`;
  grok.shell.addTable(loadingsRegrCoefsTable);

  // 0.2. Add X-Loadings
  result.xLoadings.forEach((col, idx) => {
    col.name = loadingsRegrCoefsTable.columns.getUnusedName(`${TITLE.XLOADING}${idx + 1}`);
    loadingsRegrCoefsTable.columns.add(col);
  });

  // 1. Predicted vs Reference scatter plot
  // Debias prediction (since PLS center data)
  const pred = debiasedPrediction(features, result.regressionCoefficients, input.predict, result.prediction);
  pred.name = cols.getUnusedName(`${input.predict.name} ${RESULT_NAMES.SUFFIX}`);
  cols.add(pred);
  const predictVsReferScatter = view.addViewer(DG.Viewer.scatterPlot(sourceTable, {
    title: TITLE.MODEL,
    xColumnName: input.predict.name,
    yColumnName: pred.name,
    showRegressionLine: true,
    markerType: DG.MARKER_TYPE.CIRCLE,
    showLabels: 'Always',
    help: LINK.MODEL,
  }));

  if ((input.names !== undefined) && (input.names !== null))
    predictVsReferScatter.setOptions({labelColumnNames: [input.names?.name]});

  // 2. Regression Coefficients Bar Chart
  result.regressionCoefficients.name = TITLE.REGR_COEFS;
  const regrCoeffsBar = view.addViewer(DG.Viewer.barChart(loadingsRegrCoefsTable, {
    table: loadingsRegrCoefsTable.name,
    title: TITLE.REGR_COEFS,
    splitColumnName: TITLE.FEATURE,
    valueColumnName: result.regressionCoefficients.name,
    valueAggrType: DG.AGG.AVG,
    help: LINK.COEFFS,
    showValueSelector: false,
    showStackSelector: false,
  }));

  // 3. Loadings Scatter Plot
  result.xLoadings.forEach((col, idx) => col.name = `${TITLE.XLOADING}${idx + 1}`);
  const loadingsScatter = view.addViewer(DG.Viewer.scatterPlot(loadingsRegrCoefsTable, {
    table: loadingsRegrCoefsTable.name,
    title: TITLE.LOADINGS,
    xColumnName: `${TITLE.XLOADING}1`,
    yColumnName: `${TITLE.XLOADING}${result.xLoadings.length > 1 ? '2' : '1'}`,
    markerType: DG.MARKER_TYPE.CIRCLE,
    labelColumnNames: [TITLE.FEATURE],
    help: LINK.LOADINGS,
  }));

  // 4. Scores Scatter Plot

  // 4.1) data
  const scoreNames = plsCols.map((col) => col.name);
  result.uScores.forEach((col, idx) => {
    col.name = cols.getUnusedName(`${TITLE.YSCORE}${idx + 1}`);
    cols.add(col);
    scoreNames.push(col.name);
  });

  // 4.2) create scatter
  const scoresScatter = DG.Viewer.scatterPlot(sourceTable, {
    title: TITLE.SCORES,
    xColumnName: plsCols[0].name,
    yColumnName: (plsCols.length > 1) ? plsCols[1].name : result.uScores[0].name,
    markerType: DG.MARKER_TYPE.CIRCLE,
    help: LINK.SCORES,
    showViewerFormulaLines: true,
    labelColumnNames: ((input.names !== undefined) && (input.names !== null)) ? [input.names?.name] : undefined,
  });


  // 4.3) create lines & circles  
  view.addViewer(scoresScatter);
  scoresScatter.meta.formulaLines.addAll(getLines(scoreNames));

  // 5. Explained Variances

  // 5.1) computation, source: the paper https://doi.org/10.1002/cem.2589
  //      here, we use notations from this paper
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

  // 5.2) create df
  const explVarsDF = DG.DataFrame.fromColumns([
    DG.Column.fromStrings(TITLE.COMPONENTS, compNames),
    DG.Column.fromFloat32Array(input.predict.name, yExplVars),
  ]);

  explVarsDF.name = `${sourceTable.name}(${TITLE.EXPL_VAR})`;
  grok.shell.addTable(explVarsDF);

  xExplVars.forEach((arr, idx) => explVarsDF.columns.add(DG.Column.fromFloat32Array(featuresNames[idx], arr)));

  // 5.3) bar chart
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

  // emphasize viewers in the demo case
  if (analysisType === PLS_ANALYSIS.DEMO) {    
    grok.shell.windows.help.showHelp(ui.markdown(DEMO_RESULTS_MD));

    describeElements(
      [predictVsReferScatter, scoresScatter, loadingsScatter, regrCoeffsBar, explVarsBar].map((v) => v.root),
      DEMO_RESULTS.map((info) => `<b>${info.caption}</b>\n\n${info.text}`),
      ['left', 'left', 'right', 'right', 'left'],
    );
  }
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

  let features: DG.Column[] = numCols.slice(0, numCols.length - 1);
  let predict = numCols[numCols.length - 1];
  let components = min(numColNames.length - 1, COMPONENTS.DEFAULT as number);
  let isQuadratic = false;

  const isPredictValid = () => {
    for (const col of features)
      if (col.name === predict.name)
        return false;
    return true;
  };

  const isCompConsistent = () => {
    if (components < 1)
      return false;

    const n = features.length;

    if (isQuadratic)      
      return components <= (n + 1) * n / 2 + n;

    return components <= n;
  }

  // response (to predict)
  const predictInput = ui.input.column(TITLE.PREDICT, {
    table: table,
    value: predict,
    nullable: false,
    onValueChanged: (value) => {
      predict = value;
      updateIputs();
    },
    filter: (col: DG.Column) => isValidNumeric(col),
    tooltipText: HINT.PREDICT,
  });

  // predictors (features)
  const featuresInput = ui.input.columns(TITLE.USING, {
    table: table,
    available: numColNames,
    value: features,
    onValueChanged: (val) => {
      features = val;
      updateIputs();      
    },
    tooltipText: HINT.FEATURES,
  });

  // components count
  const componentsInput = ui.input.int(TITLE.COMPONENTS, {
    value: components,
    showPlusMinus: true,
    onValueChanged: (val) => {
      components = val;
      updateIputs();
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

  const setStyle = (valid: boolean, element: HTMLElement, tooltip: string, errorMsg: string) => {
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

  const updateIputs = () => {
    const predValid = isPredictValid();
    let compValid: boolean;

    if (predValid) {
      setStyle(true, predictInput.input, HINT.PREDICT, '');
      setStyle(true, featuresInput.input, HINT.FEATURES, '');
    } else {
      setStyle(false, predictInput.input, HINT.PREDICT, ERROR_MSG.PREDICT);
      setStyle(false, featuresInput.input, HINT.FEATURES, ERROR_MSG.PREDICT);
    }

    if (components < 1) {
      setStyle(false, componentsInput.input, HINT.COMPONENTS, ERROR_MSG.COMPONENTS);
      compValid = false;  
    } else {
      compValid = isCompConsistent();

      if (compValid) {
        setStyle(true, componentsInput.input, HINT.COMPONENTS, '');
        if (predValid)
          setStyle(true, featuresInput.input, HINT.FEATURES, '');
      } else {
        const errMsg = isQuadratic ? ERROR_MSG.COMP_QUA_PLS : ERROR_MSG.COMP_LIN_PLS;
        setStyle(false, componentsInput.input, HINT.COMPONENTS, errMsg);
        setStyle(false, featuresInput.input, HINT.FEATURES, ERROR_MSG.ENOUGH);
      }
    }

    const isValid = predValid && compValid;

    dlg.getButton(TITLE.RUN).disabled = !isValid;

    return isValid;
  };

  // names of samples
  let names = (strCols.length > 0) ? strCols[0] : undefined;
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
    value: isQuadratic,
    tooltipText: HINT.QUADRATIC,
    onValueChanged: (val) => {
      isQuadratic = val;
      updateIputs();
    },
  });

  const dlg = ui.dialog({title: dlgTitle, helpUrl: dlgHelpUrl})
    .add(ui.form([predictInput, featuresInput, componentsInput, isQuadraticInput, namesInputs]))
    .addButton(TITLE.RUN, async () => {
      dlg.close();

      await performMVA({
        table: table,
        features: DG.DataFrame.fromColumns(features).columns,
        predict: predict,
        components: components,
        isQuadratic: isQuadratic,
        names: names,
      }, analysisType);
    }, undefined, dlgRunBtnTooltip)
    .show({x: X_COORD, y: Y_COORD});

  // the following delay provides correct styles (see https://reddata.atlassian.net/browse/GROK-15196)
  setTimeout(() => {
    featuresInput.value = numCols.filter((col) => col !== predict);
    features = featuresInput.value;
  }, TIMEOUT);

  grok.shell.v.append(dlg.root);
} // runMVA

/** Run multivariate analysis demo */
export async function runDemoMVA(): Promise<void> {
  grok.shell.addTableView(carsDataframe());
  grok.shell.windows.help.visible = true;
  grok.shell.windows.help.showHelp(ui.markdown(DEMO_INTRO_MD));
  grok.shell.windows.showContextPanel = false;
  grok.shell.windows.showProperties = false;

  await runMVA(PLS_ANALYSIS.DEMO);
}
