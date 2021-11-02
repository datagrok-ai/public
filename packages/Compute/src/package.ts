/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelHandler} from './model_handler';
import {selectOutliersManually} from './outliers_selection';
import {exportFuncCall} from './export_funccall';

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
  const v = DG.CardView.create({dataSource: grok.dapi.scripts, permanentFilter: '#model'});
  v.meta = new ModelHandler();
  v.name = 'Models';
  v.permanentFilter = '#model';
  grok.shell.addView(v);
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
