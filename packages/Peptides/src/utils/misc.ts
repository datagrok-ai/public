import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as C from './constants';
import * as type from './types';
import {PeptidesSettings} from './types';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import BitArray from '@datagrok-libraries/utils/src/bit-array';
import {AggregationColumns, MasksInfo, MonomerPositionStats, PositionStats} from './statistics';
import {PeptideViewer} from '../widgets/distribution';
import wu from 'wu';

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

export enum SPLIT_CATEGORY {
  SELECTION = 'Selection',
  ALL = 'All',
  PEPTIDES_SELECTION = 'Peptides selection',
}

export type DistributionLabelMap = { [key in SPLIT_CATEGORY]?: string };

export function getDistributionPanel(hist: DG.Viewer<DG.IHistogramLookSettings>, statsMap: StringDictionary,
  labelMap: DistributionLabelMap = {}): HTMLDivElement {
  const splitCol = hist.dataFrame.getCol(C.COLUMNS_NAMES.SPLIT_COL);
  const labels = [];
  const categories = splitCol.categories as SPLIT_CATEGORY[];
  const rawData = splitCol.getRawData();
  for (let categoryIdx = 0; categoryIdx < categories.length; ++categoryIdx) {
    if (!Object.values(SPLIT_CATEGORY).includes(categories[categoryIdx]))
      continue;
    const color = DG.Color.toHtml(splitCol.colors.getColor(rawData.indexOf(categoryIdx)));
    const label = ui.label(labelMap[categories[categoryIdx]] ?? categories[categoryIdx], {style: {color}});
    labels.push(label);
  }

  const result = ui.divV([ui.divV(labels), hist.root, ui.tableFromMap(statsMap)]);
  hist.root.style.maxHeight = '75px';
  return result;
}

/* Creates a table to plot activity distribution. */
export function getDistributionTable(activityCol: DG.Column<number>, selection: DG.BitSet, peptideSelection?: DG.BitSet,
): DG.DataFrame {
  const selectionMismatch = peptideSelection?.clone().xor(selection).anyTrue ?? false;
  const rowCount = activityCol.length;
  const activityColData = activityCol.getRawData();
  const activityData = new Float32Array(rowCount + selection.trueCount +
    (selectionMismatch ? (peptideSelection?.trueCount ?? 0) : 0));
  const categories: string[] = new Array(activityData.length);

  for (let i = 0, j = 0, k = 0; i < rowCount; ++i) {
    const isSelected = selection.get(i);
    activityData[i] = activityColData[i];
    categories[i] = isSelected ? SPLIT_CATEGORY.SELECTION : SPLIT_CATEGORY.ALL;
    if (isSelected) {
      activityData[rowCount + j] = activityColData[i];
      categories[rowCount + j] = SPLIT_CATEGORY.ALL;
      ++j;
    }
    if (selectionMismatch && peptideSelection?.get(i)) {
      activityData[rowCount + selection.trueCount + k] = activityColData[i];
      categories[rowCount + selection.trueCount + k] = SPLIT_CATEGORY.PEPTIDES_SELECTION;
      ++k;
    }
  }

  const splitCol = DG.Column.fromStrings(C.COLUMNS_NAMES.SPLIT_COL, categories);
  const categoryOrder = [SPLIT_CATEGORY.ALL, SPLIT_CATEGORY.SELECTION];
  if (selectionMismatch)
    categoryOrder.push(SPLIT_CATEGORY.PEPTIDES_SELECTION);
  splitCol.setCategoryOrder(categoryOrder);

  splitCol.colors.setCategorical();
  return DG.DataFrame.fromColumns([DG.Column.fromFloat32Array(C.COLUMNS_NAMES.ACTIVITY, activityData), splitCol]);
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
  grid.autoSize(1000, 400, 0, 0, true);
}

export function isSelectionEmpty(selection: type.Selection): boolean {
  for (const selectionList of Object.values(selection)) {
    if (selectionList.length !== 0)
      return false;
  }
  return true;
}

export function modifySelection(selection: type.Selection, clusterOrMonomerPosition: type.SelectionItem,
  options: type.SelectionOptions): type.Selection {
  const monomerList = selection[clusterOrMonomerPosition.positionOrClusterType];
  const monomerIndex = monomerList.indexOf(clusterOrMonomerPosition.monomerOrCluster);
  if (options.shiftPressed && options.ctrlPressed) {
    if (monomerIndex !== -1)
      monomerList.splice(monomerIndex, 1);
  } else if (options.ctrlPressed) {
    if (monomerIndex === -1)
      monomerList.push(clusterOrMonomerPosition.monomerOrCluster);
    else
      monomerList.splice(monomerIndex, 1);
  } else if (options.shiftPressed) {
    if (monomerIndex === -1)
      monomerList.push(clusterOrMonomerPosition.monomerOrCluster);
  } else {
    const selectionKeys = Object.keys(selection);
    selection = {};
    for (const posOrClustType of selectionKeys) {
      selection[posOrClustType] = [];
      if (posOrClustType === clusterOrMonomerPosition.positionOrClusterType)
        selection[posOrClustType].push(clusterOrMonomerPosition.monomerOrCluster);
    }
  }
  return selection;
}

export function highlightMonomerPosition(monomerPosition: type.SelectionItem, dataFrame: DG.DataFrame,
  monomerPositionStats: MonomerPositionStats): void {
  const bitArray = new BitArray(dataFrame.rowCount);
  if (monomerPosition.positionOrClusterType === C.COLUMNS_NAMES.MONOMER) {
    const positionStats = Object.values(monomerPositionStats);
    for (const posStat of positionStats) {
      const monomerPositionStats = (posStat as PositionStats)[monomerPosition.monomerOrCluster];
      if (typeof monomerPositionStats !== 'undefined')
        bitArray.or(monomerPositionStats.mask);
    }
  } else {
    const positionStats = monomerPositionStats[monomerPosition.positionOrClusterType];
    if (typeof positionStats !== 'undefined') {
      const monomerPositionStats = positionStats[monomerPosition.monomerOrCluster];
      if (typeof monomerPositionStats !== 'undefined')
        bitArray.or(monomerPositionStats.mask);
    }
  }

  dataFrame.rows.highlight((i) => bitArray.getBit(i));
}

export function initSelection(positionColumns: DG.Column<string>[]): type.Selection {
  const tempSelection: type.Selection = {};
  for (const posCol of positionColumns)
    tempSelection[posCol.name] = [];

  return tempSelection;
}

export function getSelectionBitset(selection: type.Selection, stats: MasksInfo): DG.BitSet | null {
  let combinedBitset: BitArray | null = null;
  const selectionEntries = Object.entries(selection);
  for (const [positionOrClusterType, selected] of selectionEntries) {
    const statsType = stats[positionOrClusterType];
    for (const monomerOrCluster of selected) {
      const statsItem = statsType[monomerOrCluster];
      combinedBitset ??= new BitArray(statsItem.mask.length, false);
      combinedBitset!.or(statsType[monomerOrCluster].mask);
    }
  }

  return (combinedBitset != null) ? DG.BitSet.fromBytes(combinedBitset!.buffer.buffer, combinedBitset!.length) : null;
}

export function areParametersEqual(o1: PeptideViewer | PeptidesSettings, o2: PeptideViewer | PeptidesSettings,
): boolean {
  return o1.sequenceColumnName === o2.sequenceColumnName && o1.activityColumnName === o2.activityColumnName &&
    o1.activityScaling === o2.activityScaling;
}

export function mutationCliffsToMaskInfo(mutationCliffs: type.MutationCliffs, rowCount: number): MasksInfo {
  const result: MasksInfo = {};
  for (const [position, monomerMap] of mutationCliffs) {
    for (const [monomer, indexMap] of monomerMap) {
      result[monomer] ??= {};
      const bitArray = new BitArray(rowCount, false);
      for (const [index, indexList] of indexMap) {
        bitArray.setTrue(index);
        for (const i of indexList)
          bitArray.setTrue(i);
      }
      result[monomer][position] = {mask: bitArray};
    }
  }
  return result;
}

export function getTotalAggColumns(viewerSelectedColNames: string[], aggColsViewer: AggregationColumns,
  aggColsModel?: AggregationColumns): [string, DG.AggregationType][] {
  const aggColsEntries = Object.entries(aggColsViewer);
  const aggColsEntriesFromSettings = aggColsModel ?
    Object.entries(aggColsModel).filter((it) => !viewerSelectedColNames.includes(it[0]) || aggColsViewer[it[0]] !== it[1]) : [];
  return aggColsEntries.concat(aggColsEntriesFromSettings);
}

export function isApplicableDataframe(table: DG.DataFrame, minRows: number = 2): boolean {
  return table.columns.bySemTypeAll(DG.SEMTYPE.MACROMOLECULE).length > 0 &&
    wu(table.columns.numerical).toArray().length > 0 && table.rowCount >= minRows;
}
