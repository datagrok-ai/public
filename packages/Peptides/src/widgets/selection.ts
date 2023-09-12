import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PeptidesModel} from '../model';

import wu from 'wu';

export function getSelectionWidget(table: DG.DataFrame, model: PeptidesModel): DG.Widget {
  const compBitset = model.getCompoundBitset();
  if (compBitset.trueCount === 0)
    return new DG.Widget(ui.divText('No compounds selected'));
  const newTable = DG.DataFrame.create(table.rowCount);
  newTable.filter.copyFrom(compBitset);
  const sourceGrid = model.analysisView.grid;
  const numericalCols = wu(table.columns.numerical);
  for (let gridColIdx = 1; gridColIdx < sourceGrid.columns.length; gridColIdx++) {
    const gridCol = sourceGrid.columns.byIndex(gridColIdx)!;
    if (!gridCol.visible)
      continue;
    const sourceCol = gridCol.column!;
    const sourceColRawData = sourceCol.getRawData();
    const sourceColCategories = sourceCol.categories;
    const getValue = numericalCols.some((col) => col.name === sourceCol.name) ? (i: number): number => sourceColRawData[i] :
      (i: number): string => sourceColCategories[sourceColRawData[i]];
    const col = newTable.columns.addNewVirtual(gridCol.name, (i) => getValue(i), sourceCol.type as DG.TYPE);
    for (const [tag, value] of sourceCol.tags)
      col.setTag(tag, value);
  }
  const newGrid = newTable.plot.grid();
  return new DG.Widget(ui.box(newGrid.root, {style: {width: '100%'}}));
}
