// PLS user interface

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PLS_ANALYSIS, ERROR_MSG, TITLE, HINT, LINK, COMPONENTS, INT, TIMEOUT, PREFIX} from './pls-constants';
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
    prediction: result[0],
    regressionCoefficients: result[1],
    tScores: result[2],
    uScores: result[3],
    xLoadings: result[4],
  };
}

/** Compute PLS components and add them to the table */
async function addPLS(input: PlsInput): Promise<void> {
  const res = await getPlsAnalysis(input);
  const plsCols = res.tScores;
  const cols = input.table.columns;

  plsCols.forEach((col, idx) => {
    col.name = cols.getUnusedName(`${PREFIX}${idx}`);
    cols.add(col);
  });
}

/** Run PLS */
export async function runPLS(type: PLS_ANALYSIS): Promise<void> {
  const table = grok.shell.t;

  if (table === null) {
    grok.shell.warning(ERROR_MSG.NO_DF);
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
      const plsInput = {
        table: table,
        features: DG.DataFrame.fromColumns(features).columns,
        predict: predict,
        components: components,
        names: names,
      };

      console.log(plsInput);
      dlg.close();

      if (type === PLS_ANALYSIS.COMPUTE_COMPONENTS)
        await addPLS(plsInput);
    }, undefined, dlgRunBtnTooltip)
    .show();

  // the following delay provides correct styles
  setTimeout(() => {
    featuresInput.value = numCols.filter((col) => col !== predict);
    features = featuresInput.value;
  }, TIMEOUT);

  grok.shell.v.append(dlg.root);
} // addPLS
