import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PeptidesModel} from '../model';

import wu from 'wu';

export function getDataSliceWidget(table: DG.DataFrame, model: PeptidesModel): DG.Widget {
  const newTable = DG.DataFrame.create(table.rowCount);
  newTable.filter.copyFrom(model.getCompoundBitset());
  const sourceGrid = model.analysisView.grid;
  const numericalCols = wu(table.columns.numerical);
  for (let gridColIdx = 1; gridColIdx < sourceGrid.columns.length; gridColIdx++) {
    const gridCol = sourceGrid.columns.byIndex(gridColIdx)!;
    if (!gridCol.visible)
      continue;
    const sourceCol = gridCol.column!;
    const sourceColRawData = sourceCol.getRawData();
    const sourceColCategories = sourceCol.categories;
    const getValue = numericalCols.has(sourceCol) ? (i: number): number => sourceColRawData[i] :
      (i: number): string => sourceColCategories[sourceColRawData[i]];
    const col = newTable.columns.addNewVirtual(gridCol.name, (i) => getValue(i), sourceCol.type as DG.TYPE);
    for (const [tag, value] of sourceCol.tags)
      col.setTag(tag, value);
  }
  const newGrid = newTable.plot.grid();
  return new DG.Widget(newGrid.root);
}
