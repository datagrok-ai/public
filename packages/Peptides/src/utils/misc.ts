import * as DG from 'datagrok-api/dg';
import * as C from './constants';
import * as type from './types';

import {AminoacidsPalettes} from '@datagrok-libraries/bio/src/aminoacids';
import {NucleotidesPalettes} from '@datagrok-libraries/bio/src/nucleotides';
import {UnknownSeqPalettes} from '@datagrok-libraries/bio/src/unknown';
import {SeqPalette} from '@datagrok-libraries/bio/src/seq-palettes';
import {WebLogo} from '@datagrok-libraries/bio/src/viewers/web-logo';

export function getPalleteByType(paletteType: string): SeqPalette {
  switch (paletteType) {
  case 'PT':
    return AminoacidsPalettes.GrokGroups;
  case 'NT':
  case 'DNA':
  case 'RNA':
    return NucleotidesPalettes.Chromatogram;
    // other
  default:
    return UnknownSeqPalettes.Color;
  }
}

export function getTypedArrayConstructor(
  maxNum: number): Uint8ArrayConstructor | Uint16ArrayConstructor | Uint32ArrayConstructor {
  return maxNum < 256 ? Uint8Array :
    maxNum < 65536 ? Uint16Array :
      Uint32Array;
}

export function getSeparator(col: DG.Column<string>): string {
  return col.getTag(C.TAGS.SEPARATOR) ?? '';
}

export function scaleActivity(
  activityScaling: string, df: DG.DataFrame, originalActivityName?: string, cloneBitset = false,
): [DG.DataFrame, string] {
  let currentActivityColName = originalActivityName ?? C.COLUMNS_NAMES.ACTIVITY;
  const flag = df.columns.names().includes(currentActivityColName) &&
    currentActivityColName === originalActivityName;
  currentActivityColName = flag ? currentActivityColName : C.COLUMNS_NAMES.ACTIVITY;
  const tempDf = df.clone(cloneBitset ? df.filter : null, [currentActivityColName]);

  let formula = (v: number) => v;
  let newColName = 'activity';
  switch (activityScaling) {
  case 'none':
    break;
  case 'lg':
    formula = (v: number) => Math.log10(v);
    newColName = `Log10(${newColName})`;
    break;
  case '-lg':
    formula = (v: number) => -Math.log10(v);
    newColName = `-Log10(${newColName})`;
    break;
  default:
    throw new Error(`ScalingError: method \`${activityScaling}\` is not available.`);
  }

  const asCol = tempDf.columns.addNewFloat(C.COLUMNS_NAMES.ACTIVITY_SCALED);
  const activityCol = df.getCol(currentActivityColName);
  asCol.init((i) => formula(activityCol.get(i)));
  df.tags['scaling'] = activityScaling;

  return [tempDf, newColName];
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
