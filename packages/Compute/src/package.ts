/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelHandler} from './model-handler';
import {selectOutliersManually} from './outliers-selection';
import {exportFuncCall} from './export-funccall';
import {_functionParametersGrid} from "./function-views/function-parameters-grid";
import {ModelCatalogView} from "./model-catalog-view";

export const _package = new DG.Package();

//name: test
export function test() {
  grok.shell.info(_package.webRoot);
}

//tags: init, autostart
export function init() {
  console.log('init');
  DG.ObjectHandler.register(new ModelHandler());
}

//name: Model Catalog
//tags: app
export function modelCatalog() {
  grok.shell.addView(new ModelCatalogView());
}

//name: manualOutlierDetectionDialog
//input: dataframe inputData
//output: dataframe augmentedInput
//output: dataframe editedInput
export async function manualOutlierSelectionDialog(inputData: DG.DataFrame) {
  const call = grok.functions.getCurrentCall();

  const IS_OUTLIER_COL_LABEL = 'isOutlier';
  const OUTLIER_REASON_COL_LABEL = 'Reason';

  if (call.options['interactive']) {
    const {augmentedInput, editedInput} = await selectOutliersManually(inputData);
    return {augmentedInput, editedInput};
  }
  return new Promise<{augmentedInput: DG.DataFrame, editedInput: DG.DataFrame}>((resolve, reject) => {
    if (!inputData.columns.byName(IS_OUTLIER_COL_LABEL)) {
      inputData.columns
        .add(DG.Column.fromBitSet(IS_OUTLIER_COL_LABEL, DG.BitSet.create(inputData.rowCount, () => false)));
    }

    if (!inputData.columns.byName(OUTLIER_REASON_COL_LABEL)) {
      inputData.columns
        .add(DG.Column.fromStrings(OUTLIER_REASON_COL_LABEL, Array.from({length: inputData.rowCount}, () => '')));
    }
    resolve({augmentedInput: inputData, editedInput: inputData});
  });
}

//name: export To Excel
//input: funccall call
//tags: export
export function exportToExcel(call: DG.FuncCall) {
  exportFuncCall(call);
}

//description: A spreadsheet that lets you interactively edit parameters and evaluate functions
//tags: functionAnalysis
//input: func f
//output: view result
export function functionParametersGrid(f: DG.Func): DG.View {
  return _functionParametersGrid(f);
}

//name: hof
export function hof() {
  let f: DG.Func = DG.Func.byName('Sin');
  let v: DG.View = functionParametersGrid(f);
  grok.shell.addView(v);
}