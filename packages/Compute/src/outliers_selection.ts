import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BitSet, DataFrame} from 'datagrok-api/dg';

export async function selectOutliersManually(inputData: DataFrame) {
  const IS_OUTLIER_COL_LABEL = 'isOutlier';
  const OUTLIER_REASON_COL_LABEL = 'Reason';
  const OUTLIER_COUNT_COL_LABEL = 'Count';

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
    'legendVisibility': 'Never',
  });

  const clearTable = () => {
    return DG.DataFrame.fromColumns([
      DG.Column.fromStrings(OUTLIER_REASON_COL_LABEL, []),
      DG.Column.fromInt32Array(OUTLIER_COUNT_COL_LABEL, new Int32Array([])),
      DG.Column.fromStrings('', []),
    ]);
  };

  const groupsListGrid = DG.Viewer.grid(clearTable());

  groupsListGrid.onCellPrepare(function(gc) {
    const btn = (reason: string) => ui.icons.delete(() => {
      for (let i = 0; i < inputData.rowCount; i++) {
        if (inputData.columns.byName(OUTLIER_REASON_COL_LABEL).get(i) === reason) {
          inputData.columns.byName(OUTLIER_REASON_COL_LABEL).set(i, '');
        }
      }
      updateTable();
    });

    if (!gc.isTableCell) {
      return;
    }

    if (gc.gridColumn.name === '') {
      gc.gridColumn.cellType = 'html';
      gc.style.element = btn(gc.grid.dataFrame?.get(OUTLIER_REASON_COL_LABEL, gc.gridRow));
    }
  });

  let isInnerModalOpened = false;

  const updateTable = () => {
    const uniqueValues: {[key:string]: number} = {};
    for (let i = 0; i < inputData.rowCount; i++) {
      const record = inputData.columns.byName(OUTLIER_REASON_COL_LABEL).get(i);
      const count = uniqueValues[record];
      if (record != '') {
        count ?
          uniqueValues[record]++ :
          uniqueValues[record] = 1;
      }
    }
    groupsListGrid.dataFrame = clearTable();
    Object.keys(uniqueValues).map((key: string) => {
      groupsListGrid.dataFrame?.rows.addNew([key, uniqueValues[key], '']);
    });
  };

  const addOutlierGroupBtn = ui.button(
    'MARK',
    () => {
      if (isInnerModalOpened) return;

      const innerDialog = ui.dialog('Add outliers group')
        .add(reasonInput)
        .onOK(
          () => {
            inputData.selection.getSelectedIndexes().forEach((selectedIndex: number) => {
              inputData.set(IS_OUTLIER_COL_LABEL, selectedIndex, true);
              inputData.set(OUTLIER_REASON_COL_LABEL, selectedIndex, reasonInput.value);
            });
            updateTable();
            inputData.selection.setAll(false);
          },
        )
        .show();
      innerDialog.onClose.subscribe(() => isInnerModalOpened = false);
      isInnerModalOpened = true;
    },
  );

  const removeOutlierGroupBtn = ui.button(
    'UNMARK',
    () => {
      if (isInnerModalOpened) return;

      inputData.selection.getSelectedIndexes().forEach((selectedIndex: number) => {
        inputData.set(IS_OUTLIER_COL_LABEL, selectedIndex, false);
        inputData.set(OUTLIER_REASON_COL_LABEL, selectedIndex, '');
      });
      updateTable();
      inputData.selection.setAll(false);
    },
  );

  inputData.onSelectionChanged.subscribe(() => {
    if (inputData.selection.trueCount === 0) {
      addOutlierGroupBtn.classList.add('disabled');
      removeOutlierGroupBtn.classList.add('disabled');
    } else {
      addOutlierGroupBtn.classList.remove('disabled');
      removeOutlierGroupBtn.classList.remove('disabled');
    }
  });

  updateTable();
  addOutlierGroupBtn.classList.add('disabled');
  removeOutlierGroupBtn.classList.add('disabled');

  const result = new Promise<{augmentedInput: DataFrame, editedInput: DataFrame}>((resolve, reject) => {
    ui.dialog('Manual outliers selection')
      .add(
        ui.divH([
          ui.block75([scatterPlot.root]),
          ui.block25([
            ui.divV([
              groupsListGrid.root,
              ui.divH([
                ui.block50([addOutlierGroupBtn], {style: {'text-align': 'center'}}),
                ui.block50([removeOutlierGroupBtn], {style: {'text-align': 'center'}}),
              ]),
            ], {style: {height: '300px'}}),
          ]),
        ]),
      )
      .onOK(() => {
        const editedInput = inputData.clone();
        editedInput.rows.filter((row) => !inputData.get(IS_OUTLIER_COL_LABEL, row.idx));
        resolve({augmentedInput: inputData, editedInput});
      })
      .onCancel(() => {
        reject(new Error('Manual outliers selection is aborted'));
      })
      .show({width: 950, height: 400});
  });

  return result;
}
