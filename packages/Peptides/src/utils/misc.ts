import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as C from './constants';
import * as type from './types';
import {getSplitterForColumn} from '@datagrok-libraries/bio/src/utils/macromolecule';
import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import {ISeqSplitted} from '@datagrok-libraries/bio/src/utils/macromolecule/types';

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
  const scaledCol: DG.Column<number> = DG.Column.float(C.COLUMNS_NAMES.ACTIVITY_SCALED, activityCol.length)
    .init((i) => {
      const val = activityColData[i];
      return val === DG.FLOAT_NULL || val === DG.INT_NULL ? val : formula(val);
    });

  return scaledCol;
}

//TODO: optimize
export function calculateSelected(df: DG.DataFrame): type.MonomerSelectionStats {
  const monomerColumns: DG.Column<string>[] = df.columns.bySemTypeAll(C.SEM_TYPES.MONOMER);
  const selectedObj: type.MonomerSelectionStats = {};
  for (const idx of df.selection.getSelectedIndexes()) {
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
  statsMap: StringDictionary, isTooltip: boolean = false): HTMLDivElement {
  const result = ui.divV([legend, hist.root, ui.tableFromMap(statsMap)]);
  result.style.minWidth = '200px';
  if (isTooltip)
    hist.root.style.maxHeight = '150px';
  return result;
}

export function prepareTableForHistogram(table: DG.DataFrame): DG.DataFrame {
  const activityCol: DG.Column<number> = table.getCol(C.COLUMNS_NAMES.ACTIVITY_SCALED);
  const splitCol: DG.Column<boolean> = table.getCol(C.COLUMNS_NAMES.SPLIT_COL);

  const rowCount = activityCol.length;
  const activityColData = activityCol.getRawData();
  const expandedData: number[] = new Array(rowCount + splitCol.stats.sum);
  const expandedMasks: boolean[] = new Array(expandedData.length);
  for (let i = 0, j = 0; i < rowCount; ++i) {
    const isSplit = splitCol.get(i)!;
    expandedData[i] = activityColData[i];
    expandedMasks[i] = isSplit;
    if (isSplit) {
      expandedData[rowCount + j] = activityColData[i];
      expandedMasks[rowCount + j] = false;
      ++j;
    }
  }

  return DG.DataFrame.fromColumns([
    DG.Column.fromList(DG.TYPE.FLOAT, activityCol.name, expandedData),
    DG.Column.fromList(DG.TYPE.BOOL, C.COLUMNS_NAMES.SPLIT_COL, expandedMasks),
  ]);
}

export async function getTemplate(sequence: string): Promise<ISeqSplitted> {
  const tempDf = DG.DataFrame.fromCsv(`sequence\n${new Array(10).fill(sequence).join('\n')}`);
  await grok.data.detectSemanticTypes(tempDf);
  const splitter = getSplitterForColumn(tempDf.getCol('sequence'));
  return splitter(sequence);
}
