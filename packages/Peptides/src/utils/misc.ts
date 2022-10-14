import * as DG from 'datagrok-api/dg';
import * as C from './constants';
import * as type from './types';

import {AminoacidsPalettes} from '@datagrok-libraries/bio/src/aminoacids';
import {NucleotidesPalettes} from '@datagrok-libraries/bio/src/nucleotides';
import {UnknownSeqPalettes} from '@datagrok-libraries/bio/src/unknown';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';

export function getTypedArrayConstructor(
  maxNum: number): Uint8ArrayConstructor | Uint16ArrayConstructor | Uint32ArrayConstructor {
  return maxNum < 256 ? Uint8Array :
    maxNum < 65536 ? Uint16Array :
      Uint32Array;
}

export function getSeparator(col: DG.Column<string>): string {
  return col.getTag(C.TAGS.SEPARATOR) ?? '';
}

export function scaleActivity(activityScaling: string, activityCol: DG.Column<number>, indexes?: number[],
): [DG.DataFrame, (x: number) => number, string] {
  const tempDf = DG.DataFrame.create(activityCol.length);

  let formula = (x: number): number => x;
  let newColName = 'activity';
  switch (activityScaling) {
  case 'none':
    break;
  case 'lg':
    formula = (x: number): number => Math.log10(x);
    newColName = `Log10(${newColName})`;
    break;
  case '-lg':
    formula = (x: number): number => -Math.log10(x);
    newColName = `-Log10(${newColName})`;
    break;
  default:
    throw new Error(`ScalingError: method \`${activityScaling}\` is not available.`);
  }
  tempDf.columns.addNewVirtual(
    C.COLUMNS_NAMES.ACTIVITY_SCALED, (i) => {
      const val = activityCol.get(indexes ? indexes[i] : i);
      return val ? formula(val) : val;
    }, DG.TYPE.FLOAT);

  return [tempDf, formula, newColName];
}

export function calculateBarsData(columns: DG.Column<string>[], selection: DG.BitSet): type.MonomerDfStats {
  const dfStats: type.MonomerDfStats = {};
  const columnsLen = columns.length;

  for (let colIndex = 0; colIndex < columnsLen; colIndex++) {
    const col = columns[colIndex];
    dfStats[col.name] = calculateSingleBarData(col, selection);
  }

  return dfStats;
}

export function calculateSingleBarData(col: DG.Column<string>, selection: DG.BitSet): type.MonomerColStats {
  const colLen = col.length;
  const colStats: type.MonomerColStats = {};
  col.categories.forEach((monomer) => colStats[monomer] = {count: 0, selected: 0});

  for (let rowIndex = 0; rowIndex < colLen; rowIndex++) {
    const monomerStats = colStats[col.get(rowIndex)!];
    monomerStats.count += 1;
    monomerStats.selected += +selection.get(rowIndex);
  }

  return colStats;
}

export function isGridCellInvalid(gc: DG.GridCell | null): boolean {
  return !gc || !gc.cell.value || !gc.tableColumn || gc.tableRowIndex == null || gc.tableRowIndex == -1 ||
    gc.cell.value == DG.INT_NULL || gc.cell.value == DG.FLOAT_NULL;
}
