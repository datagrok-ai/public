import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import $ from 'cash-dom';

const getBalloonContainer = () => {
  return document.querySelector('.d4-balloon-container');
};

export class OutliersSelectionViewer extends DG.JsViewer {
  constructor() {
    super();
  }

  onFrameAttached(dataFrame: DG.DataFrame): void {
    this.render();
  }

  onTableAttached(): void {
    this.render();
  }

  render() {
    if (!this.dataFrame || this.dataFrame.rowCount <= 0) {
      this.root.innerHTML = '';
      return;
    }

    const inputData = this.dataFrame;
    const IS_OUTLIER_COL_LABEL = 'isOutlier';
    const IS_OUTLIER_SECONDARY_COL_LABEL = 'isOutlierSecondary';
    const OUTLIER_RATIONALE_COL_LABEL = 'Rationale';
    const OUTLIER_COUNT_COL_LABEL = 'Count';
    const IS_GROUP_CONFIRMED_LABEL = 'isConfirmed';

    if (!inputData.columns.byName(IS_OUTLIER_COL_LABEL)) {
      inputData.columns
        .add(DG.Column.fromList(DG.TYPE.STRING, IS_OUTLIER_COL_LABEL, new Array(inputData.rowCount).fill('Inlier')));
    }
    if (!inputData.columns.byName(IS_OUTLIER_SECONDARY_COL_LABEL)) {
      inputData.columns
        .add(DG.Column.fromList(DG.TYPE.STRING, IS_OUTLIER_SECONDARY_COL_LABEL, new Array(inputData.rowCount).fill('Inlier')));
    }

    if (!inputData.columns.byName(OUTLIER_RATIONALE_COL_LABEL)) {
      inputData.columns
        .add(DG.Column.fromStrings(OUTLIER_RATIONALE_COL_LABEL, Array.from({length: inputData.rowCount}, () => '')));
    }

    const clearTable = () => {
      return DG.DataFrame.fromColumns([
        DG.Column.fromStrings(OUTLIER_RATIONALE_COL_LABEL, []),
        DG.Column.fromInt32Array(OUTLIER_COUNT_COL_LABEL, new Int32Array([])),
        DG.Column.fromStrings('Actions', []),
        DG.Column.fromBitSet(IS_GROUP_CONFIRMED_LABEL, DG.BitSet.create(0, () => false)),
      ]);
    };

    const groupsListGrid = DG.Viewer.grid(clearTable());
    groupsListGrid.root.style.width = '100%';
    groupsListGrid.root.style.minWidth = '230px';
    groupsListGrid.columns.setVisible([OUTLIER_RATIONALE_COL_LABEL, OUTLIER_COUNT_COL_LABEL, 'Actions']);

    this.root.append(
      ui.divV([
        ui.divV([
          groupsListGrid.root,
        ], {style: {'height': '75%'}}),
      ], {style: {'height': '100%'}}),
    );

    const confirmOutliers = () => {
      if (!groupsListGrid.dataFrame) return;

      inputData.col(IS_OUTLIER_COL_LABEL)?.meta.markers.assign('Outlier', DG.MARKER_TYPE.OUTLIER);
      inputData.col(IS_OUTLIER_COL_LABEL)?.meta.markers.assign('Last Inlier', DG.MARKER_TYPE.STAR);
      inputData.col(IS_OUTLIER_COL_LABEL)?.meta.colors.setCategorical({'Outlier': DG.Color.scatterPlotMarker,
        'Inlier': DG.Color.scatterPlotMarker, 'Last Inlier': DG.Color.red});
      inputData.col(IS_OUTLIER_SECONDARY_COL_LABEL)?.meta.markers.assign('Outlier', DG.MARKER_TYPE.OUTLIER);
      inputData.col(IS_OUTLIER_SECONDARY_COL_LABEL)?.meta.markers.assign('Last Inlier', DG.MARKER_TYPE.STAR);
      inputData.col(IS_OUTLIER_SECONDARY_COL_LABEL)?.meta.colors.setCategorical({'Outlier': DG.Color.cyan,
        'Inlier': DG.Color.cyan, 'Last Inlier': DG.Color.red});

      const newRationale =
        groupsListGrid.dataFrame?.cell(groupsListGrid.dataFrame.rowCount-1, OUTLIER_RATIONALE_COL_LABEL).value;

      inputData.selection.getSelectedIndexes().forEach((selectedIndex: number) => {
        inputData.col(IS_OUTLIER_COL_LABEL)!.set(selectedIndex, 'Outlier', false);
        inputData.col(IS_OUTLIER_SECONDARY_COL_LABEL)!.set(selectedIndex, 'Outlier', false);
        inputData.col(OUTLIER_RATIONALE_COL_LABEL)!.set(
          selectedIndex,
          newRationale,
          false,
        );
      });
      inputData.selection.setAll(false);
      inputData.fireValuesChanged();

      $(getBalloonContainer()).empty();
    };

    const cancelOutliers = () => {
      if (!groupsListGrid.dataFrame) return;

      groupsListGrid.dataFrame.rows.removeAt(groupsListGrid.dataFrame.rowCount - 1);
      inputData.selection.setAll(false);

      $(getBalloonContainer()).empty();
    };

    groupsListGrid.onCellPrepare((gc) => {
      if (!gc.isTableCell)
        return;


      const confirmBtn = () => ui.div(
        ui.iconFA('check', confirmOutliers, 'Confirm the outliers'), {style: {'textAlign': 'center', 'margin': '6px'}},
      );

      const cancelBtn = () => ui.div(
        ui.iconFA('times', cancelOutliers, 'Cancel the outliers'), {style: {'textAlign': 'center', 'margin': '6px'}},
      );

      const deleteBtn = (rationale: string) => ui.div(
        ui.icons.delete(() => {
          for (let i = 0; i < inputData.rowCount; i++) {
            if (inputData.columns.byName(OUTLIER_RATIONALE_COL_LABEL).get(i) === rationale) {
              inputData.columns.byName(OUTLIER_RATIONALE_COL_LABEL).set(i, '', false);
              inputData.columns.byName(IS_OUTLIER_COL_LABEL).set(i, 'Inlier', false);
              inputData.columns.byName(IS_OUTLIER_SECONDARY_COL_LABEL).set(i, 'Inlier', false);
            }
          }
          inputData.fireValuesChanged();
        }, 'Remove the outliers group'), {style: {textAlign: 'center', margin: '6px'}},
      );

      if (gc.gridColumn.name === 'Actions') {
        gc.gridColumn.cellType = 'html';
        gc.style.element =
          ui.divH([
            ...!gc.grid.dataFrame?.get(IS_GROUP_CONFIRMED_LABEL, gc.gridRow) ?
              [confirmBtn(), cancelBtn()] :
              [deleteBtn(gc.grid.dataFrame?.get(OUTLIER_RATIONALE_COL_LABEL, gc.gridRow))],
          ]);
      }
    });

    groupsListGrid.onCellValueEdited.subscribe((editedCell) => {
      if (!groupsListGrid.dataFrame || editedCell.gridColumn.name !== OUTLIER_RATIONALE_COL_LABEL) return;

      if (!(groupsListGrid.dataFrame.columns as DG.ColumnList)
        .byName(IS_GROUP_CONFIRMED_LABEL).get(editedCell.gridRow)) return;

      const uniqueValues= new Set<string>();
      for (let i = 0; i < groupsListGrid.dataFrame?.rowCount; i++) {
        const record = groupsListGrid.dataFrame?.columns.byName(OUTLIER_RATIONALE_COL_LABEL).get(i);
        uniqueValues.add(record);
      }

      for (let i = 0; i < inputData.rowCount; i++) {
        const currentCellValue = inputData.columns.byName(OUTLIER_RATIONALE_COL_LABEL).get(i);
        if (currentCellValue != '' && !uniqueValues.has(currentCellValue))
          inputData.columns.byName(OUTLIER_RATIONALE_COL_LABEL).set(i, editedCell.cell.value);
      }
    });

    const updateGroupsDf = () => {
      if ((inputData.columns as DG.ColumnList).byName(OUTLIER_RATIONALE_COL_LABEL)) {
        for (let i = 0; i < inputData.rowCount; i++) {
          if (inputData.columns.byName(IS_OUTLIER_COL_LABEL).get(i) === 'Outlier') {
            if (inputData.columns.byName(OUTLIER_RATIONALE_COL_LABEL).get(i) === '')
              (inputData.columns as DG.ColumnList).byName(OUTLIER_RATIONALE_COL_LABEL).set(i, 'Manual', false);
          } else {
            if (inputData.columns.byName(OUTLIER_RATIONALE_COL_LABEL).get(i) !== '')
              (inputData.columns as DG.ColumnList).byName(OUTLIER_RATIONALE_COL_LABEL).set(i, '', false);
          }
        }
      }

      const uniqueValues: {[key:string]: number} = {};
      if ((inputData.columns as DG.ColumnList).byName(OUTLIER_RATIONALE_COL_LABEL)) {
        for (let i = 0; i < inputData.rowCount; i++) {
          const record = inputData.columns.byName(OUTLIER_RATIONALE_COL_LABEL).get(i);
          const count = uniqueValues[record];
          if (record != '') {
            count ?
              uniqueValues[record]++ :
              uniqueValues[record] = 1;
          }
        }
      }
      groupsListGrid.dataFrame = clearTable();
      groupsListGrid.columns.setVisible([OUTLIER_RATIONALE_COL_LABEL, OUTLIER_COUNT_COL_LABEL, 'Actions']);

      Object.keys(uniqueValues).map((key: string) => {
        groupsListGrid.dataFrame?.rows.addNew([key, uniqueValues[key], '', true]);
      });
    };

    inputData.onDataChanged.subscribe(updateGroupsDf);

    inputData.onSelectionChanged.subscribe(() => {
      if (!groupsListGrid.dataFrame || !inputData.selection.trueCount) return;

      const balloon = getBalloonContainer();
      if (!balloon || $(balloon).is(':empty')) {
        grok.shell.info(ui.div([
          ui.divText('Confirm selected points as outliers?'),
          ui.divH([
            ui.bigButton('Confirm', confirmOutliers),
            ui.button('Cancel', cancelOutliers),
          ]),
        ]));
      }

      const isConfirmedColumn = (groupsListGrid.dataFrame.columns as DG.ColumnList)
        .byName(IS_GROUP_CONFIRMED_LABEL).toList();

      if (isConfirmedColumn.some((value) => value === false)) {
        groupsListGrid.dataFrame.set(
          OUTLIER_COUNT_COL_LABEL,
          groupsListGrid.dataFrame.rowCount-1,
          inputData.selection.trueCount,
        );
      } else {
        const newRow = groupsListGrid.dataFrame.rows.addNew(['Manual', inputData.selection.trueCount, '', false]);
        groupsListGrid.dataFrame.currentCell = groupsListGrid.dataFrame.cell(newRow.idx, OUTLIER_RATIONALE_COL_LABEL);
      }
    });

    inputData.selection.setAll(false);

    updateGroupsDf();
  }
}
