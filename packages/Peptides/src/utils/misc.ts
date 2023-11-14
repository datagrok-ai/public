import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as C from './constants';
import * as type from './types';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';

export function getTypedArrayConstructor(
  maxNum: number): Uint8ArrayConstructor | Uint16ArrayConstructor | Uint32ArrayConstructor {
  return maxNum < 256 ? Uint8Array :
    maxNum < 65536 ? Uint16Array :
      Uint32Array;
}

export function getSeparator(col: DG.Column<string>): string {
  return col.getTag(C.TAGS.SEPARATOR) ?? '';
}

export function scaleActivity(activityCol: DG.Column<number>, scaling: C.SCALING_METHODS = C.SCALING_METHODS.NONE,
): DG.Column<number> {
  let formula = (x: number): number => x;
  switch (scaling) {
  case C.SCALING_METHODS.NONE:
    break;
  case C.SCALING_METHODS.LG:
    formula = (x: number): number => Math.log10(x);
    break;
  case C.SCALING_METHODS.MINUS_LG:
    formula = (x: number): number => -Math.log10(x);
    break;
  default:
    throw new Error(`ScalingError: method \`${scaling}\` is not available.`);
  }
  const activityColData = activityCol.getRawData();
  const scaledCol: DG.Column<number> = DG.Column.float(C.COLUMNS_NAMES.ACTIVITY, activityCol.length)
    .init((i) => {
      const val = activityColData[i];
      return val === DG.FLOAT_NULL || val === DG.INT_NULL ? val : formula(val);
    });
  scaledCol.setTag(C.TAGS.ANALYSIS_COL, `${true}`);
  scaledCol.setTag(DG.TAGS.FORMULA, scaling);
  return scaledCol;
}

//TODO: optimize
export function calculateSelected(df: DG.DataFrame): type.SelectionStats {
  const monomerColumns: DG.Column<string>[] = df.columns.bySemTypeAll(C.SEM_TYPES.MONOMER);
  const selectedObj: type.SelectionStats = {};
  const selectedIndexes = df.filter.clone().and(df.selection).getSelectedIndexes();
  for (const idx of selectedIndexes) {
    for (const col of monomerColumns) {
      const monomer = col.get(idx);
      if (!monomer)
        continue;

      selectedObj[col.name] ??= {};
      selectedObj[col.name][monomer] ??= 0;
      selectedObj[col.name][monomer] += 1;
    }
  }

  return selectedObj;
}

export function extractColInfo(col: DG.Column<string>): type.RawColumn {
  return {
    name: col.name,
    cat: col.categories,
    rawData: col.getRawData(),
  };
}

export function getStatsSummary(legend: HTMLDivElement, hist: DG.Viewer<DG.IHistogramLookSettings>,
  statsMap: StringDictionary): HTMLDivElement {
  const result = ui.divV([legend, hist.root, ui.tableFromMap(statsMap)]);
  hist.root.style.maxHeight = '75px';
  return result;
}

/* Creates a table to plot activity distribution. */
export function getDistributionTable(activityCol: DG.Column<number>, selection: DG.BitSet, peptideSelection?: DG.BitSet): DG.DataFrame {
  // const activityCol: DG.Column<number> = table.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
  // const splitCol: DG.Column<boolean> = table.getCol(C.COLUMNS_NAMES.SPLIT_COL);
  const isCustomSelection = peptideSelection?.clone().xor(selection).anyTrue ?? false;

  const rowCount = activityCol.length;
  const activityColData = activityCol.getRawData();
  const expandedData = new Float32Array(rowCount + selection.trueCount + (peptideSelection?.trueCount ?? 0));
  const expandedMasks = new Uint8Array(expandedData.length);

  for (let i = 0, j = 0, k = 0; i < rowCount; ++i) {
    const isSplit = selection.get(i);
    expandedData[i] = activityColData[i];
    expandedMasks[i] = isSplit ? 1 : 0;
    if (isSplit) {
      expandedData[rowCount + j] = activityColData[i];
      expandedMasks[rowCount + j] = 0;
      ++j;
    }
    if (isCustomSelection) {
      expandedData[rowCount + selection.trueCount + k] = activityColData[i];
      expandedMasks[rowCount + selection.trueCount + k] = 2;
      ++k;
    }
  }

  return DG.DataFrame.fromColumns([
    DG.Column.fromFloat32Array(activityCol.name, expandedData, rowCount),
    DG.Column.string(C.COLUMNS_NAMES.SPLIT_COL, rowCount).init((i) => expandedMasks[i].toString()),
  ]);
}

export function addExpandIcon(grid: DG.Grid): void {
  const fullscreenIcon = ui.iconFA('expand-alt', () => {
    const fullscreenGrid = grid.dataFrame.plot.grid();
    setGridProps(fullscreenGrid);
    fullscreenGrid.root.style.height = '100%';
    const pairsFullscreenDialog = ui.dialog(grid.dataFrame.name);
    pairsFullscreenDialog.add(fullscreenGrid.root);
    pairsFullscreenDialog.showModal(true);
    fullscreenGrid.invalidate();
  });
  grid.root.appendChild(fullscreenIcon);
  fullscreenIcon.style.position = 'absolute';
  fullscreenIcon.style.right = '0px';
  fullscreenIcon.style.top = '0px';
  fullscreenIcon.style.visibility = 'hidden';
  grid.root.addEventListener('mouseenter', (_) => {
    fullscreenIcon.style.visibility = 'visible';
  });
  grid.root.addEventListener('mouseleave', (_) => {
    fullscreenIcon.style.visibility = 'hidden';
  });
}

export function setGridProps(grid: DG.Grid): void {
  grid.props.allowEdit = false;
  grid.props.allowRowSelection = false;
  grid.props.allowBlockSelection = false;
  grid.props.allowColSelection = false;
  grid.props.showRowHeader = false;
  grid.props.showCurrentRowIndicator = false;
  grid.root.style.width = '100%';
  grid.root.style.maxWidth = '100%';
  grid.autoSize(1000, 1000, 0, 0, true);
}
