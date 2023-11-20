import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {PeptidesModel} from '../model';

import wu from 'wu';
import {COLUMNS_NAMES} from '../utils/constants';
import {addExpandIcon} from '../utils/misc';
import {CellRendererOptions, setWebLogoRenderer, WebLogoBounds} from '../utils/cell-renderer';
import {CachedWebLogoTooltip, SelectionItem} from '../utils/types';
import {TooltipOptions} from '../utils/tooltips';
import {calculateMonomerPositionStatistics} from '../utils/algorithms';

export function getSelectionWidget(table: DG.DataFrame, model: PeptidesModel): HTMLElement {
  const compBitset = model.getVisibleSelection();
  if (compBitset.trueCount === 0)
    return ui.divText('No compounds selected');
  const newTable = DG.DataFrame.create(table.rowCount);
  newTable.name = 'Selected compounds';
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
    const col = sourceCol.name === COLUMNS_NAMES.ACTIVITY ?
      newTable.columns.addNewFloat(gridCol.name).init((i) => getValue(i)) :
      newTable.columns.addNewVirtual(gridCol.name, (i) => getValue(i), sourceCol.type as DG.TYPE);
    for (const [tag, value] of sourceCol.tags)
      col.setTag(tag, value);
  }
  const grid = newTable.plot.grid();
  grid.props.showRowHeader = false;
  grid.root.style.maxWidth = '100%';

  DG.debounce(ui.onSizeChanged(grid.root), 50).subscribe((_) => {
    const panel = grid.root.parentElement;
    if (panel?.parentElement?.classList.contains('panel-content')) {
      grid.root.style.height = 'calc(100% - 20px)';
      grid.root.style.width = 'calc(100% - 20px)';
      grid.root.style.position = 'absolute';
      grid.root.style.right = '10px';
      grid.root.style.top = '10px';
    }
  });

  addExpandIcon(grid);

  const gridHost = ui.box(grid.root);
  gridHost.style.marginLeft = '0px';
  setTimeout(() => {
    for (let gridColIdx = 1; gridColIdx < sourceGrid.columns.length; gridColIdx++) {
      const gridCol = sourceGrid.columns.byIndex(gridColIdx)!;
      if (!gridCol.visible)
        continue;
      grid.col(gridCol.name)!.width = gridCol.width;
    }
  }, 500);

  const mpStats = calculateMonomerPositionStatistics(grid.dataFrame, model.positionColumns.toArray());

  const cachedWebLogoTooltip: CachedWebLogoTooltip = {bar: '', tooltip: null};
  const webLogoBounds: WebLogoBounds = {};
  const cellRendererOptions: CellRendererOptions = {isSelectionTable: true, cachedWebLogoTooltip, webLogoBounds};
  const tooltipOptions: TooltipOptions = {x: 0, y: 0, monomerPosition: {} as SelectionItem, mpStats};

  setWebLogoRenderer(grid, model, cellRendererOptions, tooltipOptions);

  return gridHost;
}
