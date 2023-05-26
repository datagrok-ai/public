import {StringDictionary} from '@datagrok-libraries/utils/src/type-declarations';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as C from './constants';
import * as type from './types';

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
  const scaledCol: DG.Column<number> = DG.Column.float(C.COLUMNS_NAMES.ACTIVITY_SCALED, activityCol.length)
    .init((i) => {
      const val = activityCol.get(i);
      return val ? formula(val) : val;
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

// export function isGridCellInvalid(gc: DG.GridCell | null): boolean {
//   return !gc || !gc.cell.value || !gc.tableColumn || gc.tableRowIndex === null || gc.tableRowIndex === -1 ||
//     gc.cell.value === DG.INT_NULL || gc.cell.value === DG.FLOAT_NULL;
// }

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
