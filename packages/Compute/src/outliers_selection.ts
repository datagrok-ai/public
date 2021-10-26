import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BitSet, DataFrame} from 'datagrok-api/dg';

export async function selectOutliersManually(inputData: DataFrame) {
  const IS_OUTLIER_COL_LABEL = 'isOutlier';

  const OUTLIER_REASON_COL_LABEL = 'Reason';

  if (!inputData.columns.byName(IS_OUTLIER_COL_LABEL)) {
    inputData.columns
      .add(DG.Column.fromBitSet(IS_OUTLIER_COL_LABEL, BitSet.create(inputData.rowCount, () => false)));
  }

  if (!inputData.columns.byName(OUTLIER_REASON_COL_LABEL)) {
    inputData.columns
      .add(DG.Column.fromStrings(OUTLIER_REASON_COL_LABEL, Array.from({length: inputData.rowCount}, () => '')));
  }

  const reasonInput = ui.textInput('Reason', '');
  const scatterPlot = DG.Viewer.scatterPlot(inputData, {
    'color': OUTLIER_REASON_COL_LABEL,
    'lassoTool': true,
  });

  let isInnerModalOpened = false;

  const addOutlierGroupBtn = {
    text: 'ADD OUTLIERS GROUP',
    action: () => {
      if (isInnerModalOpened) return;

      const innerDialog = ui.dialog('Add outliers group')
        .add(reasonInput)
        .onOK(
          () => {
            inputData.selection.getSelectedIndexes().forEach((selectedIndex: number) => {
              inputData.set(IS_OUTLIER_COL_LABEL, selectedIndex, true);
              inputData.set(OUTLIER_REASON_COL_LABEL, selectedIndex, reasonInput.value);
            });
            inputData.selection.setAll(false);
          },
        )
        .show();
      innerDialog.onClose.subscribe(() => isInnerModalOpened = false);
      isInnerModalOpened = true;
    },
  };

  const removeOutlierGroupBtn = {
    text: 'REMOVE OUTLIERS GROUP',
    action: () => {
      if (isInnerModalOpened) return;

      inputData.selection.getSelectedIndexes().forEach((selectedIndex: number) => {
        inputData.set(IS_OUTLIER_COL_LABEL, selectedIndex, false);
        inputData.set(OUTLIER_REASON_COL_LABEL, selectedIndex, '');
      });
      inputData.selection.setAll(false);
    },
  };

  const result = new Promise<{augmentedInput: DataFrame, editedInput: DataFrame}>((resolve, reject) => {
    ui.dialog('Manual outliers selection')
      .add(
        ui.divH([scatterPlot.root]),
      )
      .addButton(addOutlierGroupBtn.text, addOutlierGroupBtn.action)
      .addButton(removeOutlierGroupBtn.text, removeOutlierGroupBtn.action)
      .onOK(() => {
        const editedInput = inputData.clone();
        editedInput.rows.filter((row) => !inputData.get(IS_OUTLIER_COL_LABEL, row.idx));
        resolve({augmentedInput: inputData, editedInput});
      })
      .onCancel(() => {
        reject(new Error('Manual outliers selection is aborted'));
      })
      .show({width: 900, height: 400});
  });

  return result;
}
