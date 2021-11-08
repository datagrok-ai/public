import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {_package} from './package';
import {std, mean, abs} from 'mathjs';

export async function selectOutliersManually(inputData: DG.DataFrame) {
  const IS_OUTLIER_COL_LABEL = 'isOutlier';
  const OUTLIER_REASON_COL_LABEL = 'Reason';
  const OUTLIER_COUNT_COL_LABEL = 'Count';

  if (!inputData.columns.byName(IS_OUTLIER_COL_LABEL)) {
    inputData.columns
      .add(DG.Column.fromBitSet(IS_OUTLIER_COL_LABEL, DG.BitSet.create(inputData.rowCount, () => false)));
  }

  if (!inputData.columns.byName(OUTLIER_REASON_COL_LABEL)) {
    inputData.columns
      .add(DG.Column.fromStrings(OUTLIER_REASON_COL_LABEL, Array.from({length: inputData.rowCount}, () => '')));
  }

  const initialData = inputData.clone();

  const reasonInput = ui.textInput('Reason', '');
  const scatterPlot = DG.Viewer.scatterPlot(inputData, {
    'color': OUTLIER_REASON_COL_LABEL,
    'lassoTool': true,
    'legendVisibility': 'Never',
    'filterByZoom': false,
  });
  scatterPlot.root.style.height = '100%';
  scatterPlot.root.style.width = '100%';

  const clearTable = () => {
    return DG.DataFrame.fromColumns([
      DG.Column.fromStrings(OUTLIER_REASON_COL_LABEL, []),
      DG.Column.fromInt32Array(OUTLIER_COUNT_COL_LABEL, new Int32Array([])),
      DG.Column.fromStrings('', []),
    ]);
  };

  const groupsListGrid = DG.Viewer.grid(clearTable());
  groupsListGrid.root.style.width = '100%';

  groupsListGrid.onCellPrepare(function(gc) {
    const btn = (reason: string) => ui.div(
      ui.icons.delete(() => {
        for (let i = 0; i < inputData.rowCount; i++) {
          if (inputData.columns.byName(OUTLIER_REASON_COL_LABEL).get(i) === reason) {
            inputData.columns.byName(OUTLIER_REASON_COL_LABEL).set(i, '');
            inputData.columns.byName(IS_OUTLIER_COL_LABEL).set(i, false);
          }
        }
        updateTable();
      }, 'Remove the outliers group'), {style: {'text-align': 'center', 'margin': '6px'}},
    );

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

  let shouldCancel = true;

  const cancelAllChanges = () => {
    if (shouldCancel) {
      (inputData.columns as DG.ColumnList).byName(IS_OUTLIER_COL_LABEL).init(
        (index) => initialData.get(IS_OUTLIER_COL_LABEL, index),
      );
      (inputData.columns as DG.ColumnList).byName(OUTLIER_REASON_COL_LABEL).init(
        (index) => initialData.get(OUTLIER_REASON_COL_LABEL, index),
      );
    }
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
    'Mark the selected points as outliers',
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
    'Remove the selected points from oultiers list',
  );

  const autoOutlierGroupBtn = ui.button(
    'STDDEV RULE...',
    () => {
      if (isInnerModalOpened) return;

      const intInput = ui.intInput('N', 3);
      const columnInput = ui.columnInput('Values', inputData, null);
      const autoDetectionDialog = ui.dialog('Edit standard deviation rule')
        .add(intInput.root)
        .add(columnInput.root)
        .onOK(()=>{
          const values = [...columnInput.value.values()];
          const meanValue = mean(values);
          const stddev = std(values);
          inputData.selection.init((index)=> abs(meanValue - values[index]) > stddev*intInput.value);
          inputData.selection.getSelectedIndexes().forEach((selectedIndex: number) => {
            inputData.set(IS_OUTLIER_COL_LABEL, selectedIndex, true);
            inputData.set(OUTLIER_REASON_COL_LABEL, selectedIndex, `${intInput.value}x stddev rule`);
          });
          updateTable();
          inputData.selection.setAll(false);
        })
        .show();
      autoDetectionDialog.onClose.subscribe(() => isInnerModalOpened = false);
      isInnerModalOpened = true;
    },
    'Edit standard deviation rule',
  );

  inputData.selection.setAll(false);

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

  const result = new Promise<{augmentedInput: DG.DataFrame, editedInput: DG.DataFrame}>((resolve, reject) => {
    const selectionDialog = ui.dialog({title: 'Select outliers', helpUrl: `${_package.webRoot}/help/outliers_selection/main.md`})
      .add(
        ui.info(
          ui.div([
            ui.p('Hold the “SHIFT” key and start to draw a freehand selection on the scatterplot area'),
          ]),
        ),
      )
      .add(
        ui.divH([
          ui.block75([scatterPlot.root]),
          ui.block25([
            ui.divV([
              ui.divH([
                addOutlierGroupBtn, removeOutlierGroupBtn, autoOutlierGroupBtn,
              ], {style: {'text-align': 'center'}}),
              groupsListGrid.root,
            ], {style: {height: '75%'}}),
          ]),
        ], {style: {height: '100%'}}),
      )
      .onOK(() => {
        shouldCancel = false;
        const editedInput = inputData.clone();
        editedInput.rows.filter((row) => !inputData.get(IS_OUTLIER_COL_LABEL, row.idx));
        resolve({augmentedInput: inputData, editedInput});
      })
      .onCancel(() => {
        const editedInput = initialData.clone();
        editedInput.rows.filter((row) => !initialData.get(IS_OUTLIER_COL_LABEL, row.idx));
        cancelAllChanges();
        resolve({augmentedInput: inputData, editedInput});
      });

    // TODO: change the resolving strategy after API update
    // onClose is called before onOK or onCancel, thus, the timeout is requierd
    selectionDialog.onClose.subscribe(() =>
      setTimeout(()=>{
        const editedInput = initialData.clone();
        editedInput.rows.filter((row) => !initialData.get(IS_OUTLIER_COL_LABEL, row.idx));
        cancelAllChanges();
        resolve({augmentedInput: inputData, editedInput});
      }, 100));
    selectionDialog.show({width: 1000, height: 800});
  });

  return result;
}
