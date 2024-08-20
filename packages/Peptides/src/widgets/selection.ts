import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import wu from 'wu';
import {addExpandIcon} from '../utils/misc';
import {setWebLogoRenderer, WebLogoBounds, WebLogoCellRendererOptions} from '../utils/cell-renderer';
import {CachedWebLogoTooltip, SelectionItem} from '../utils/types';
import {TooltipOptions} from '../utils/tooltips';
import {calculateMonomerPositionStatistics} from '../utils/algorithms';
import {AggregationColumns} from '../utils/statistics';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';

export type SelectionWidgetOptions = {
  tableSelection: DG.BitSet, gridColumns: DG.GridColumnList, positionColumns: DG.Column<string>[],
  activityColumn: DG.Column<number>, columns: AggregationColumns, colorPalette: SeqPalette, isAnalysis: boolean,
};

/**
 * Creates selection grid with WebLogo rendered in its header.
 * @param table - table with selected sequences.
 * @param options - options for selection widget.
 * @return - selection grid.
 */
export function getSelectionWidget(table: DG.DataFrame, options: SelectionWidgetOptions): HTMLElement {
  if (options.tableSelection.trueCount === 0)
    return ui.divText('No compounds selected');


  const newTable = DG.DataFrame.create(table.rowCount);
  newTable.name = 'Selected compounds';
  newTable.filter.copyFrom(options.tableSelection);
  const numericalCols = wu(table.columns.numerical);
  let gridSortOrder: {cols: DG.Column[], types: boolean[]} | null = null;
  for (let gridColIdx = 1; gridColIdx < options.gridColumns.length; gridColIdx++) {
    const gridCol = options.gridColumns.byIndex(gridColIdx)!;
    if (!gridCol.visible)
      continue;


    const sourceCol = gridCol.column!;
    if (sourceCol.type === DG.COLUMN_TYPE.BOOL)
      continue;

    // restore sorting
    if (gridSortOrder === null && gridCol.grid && gridCol.grid.sortByColumns?.length > 0 &&
       gridCol.grid.sortTypes?.length === gridCol.grid.sortByColumns.length)
      gridSortOrder = {cols: gridCol.grid.sortByColumns, types: gridCol.grid.sortTypes};


    // bigint or qnum columns can throw an error when getting raw data (might be fixed in future versions)
    let sourceColRawData: ArrayLike<number> | null = null;
    let sourceColCategories: string[] | null = null;

    try {
      sourceColRawData = sourceCol.getRawData();
      sourceColCategories = sourceCol.categories;
    } catch (_e) {
    }
    const getValue = !sourceColRawData || !sourceColCategories ?
      (i: number): any => sourceCol.get(i) :
      numericalCols
        .some((col) => col.name === sourceCol.name) ? (i: number): number => sourceColRawData[i] :
        (i: number): string => sourceColCategories[sourceColRawData[i]];
    const col = sourceCol.name === options.activityColumn.name ?
      newTable.columns.addNewFloat(gridCol.name).init((i) => getValue(i) as number) :
      newTable.columns.addNewVirtual(gridCol.name, (i) => getValue(i) as unknown, sourceCol.type as DG.TYPE);
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
    if (gridSortOrder !== null)
      grid.sort(gridSortOrder.cols, gridSortOrder.types);
    for (let gridColIdx = 1; gridColIdx < options.gridColumns.length; gridColIdx++) {
      const originalGridCol = options.gridColumns.byIndex(gridColIdx)!;
      if (!originalGridCol.visible)
        continue;


      if (!grid.col(originalGridCol.name))
        continue;


      grid.col(originalGridCol.name)!.width = originalGridCol.width;
    }
  }, 500);

  const mpStats = calculateMonomerPositionStatistics(options.activityColumn, newTable.filter,
    options.positionColumns, {isFiltered: newTable.filter.anyTrue || newTable.filter.anyFalse});

  const cachedWebLogoTooltip: () => CachedWebLogoTooltip = () => {
    return {bar: '', tooltip: null};
  };
  const webLogoBounds: () => WebLogoBounds = () => {
    return {};
  };
  const cellRendererOptions: WebLogoCellRendererOptions = {
    isSelectionTable: true, cachedWebLogoTooltip, webLogoBounds,
    colorPalette: () => options.colorPalette,
  };
  const tooltipOptions: TooltipOptions = {x: 0, y: 0, monomerPosition: {} as SelectionItem, mpStats};

  if (options.isAnalysis) {
    setWebLogoRenderer(grid, mpStats, options.positionColumns, options.activityColumn, cellRendererOptions,
      tooltipOptions);
  }

  return gridHost;
}
