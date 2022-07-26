import * as DG from 'datagrok-api/dg';
import * as C from './constants';

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

export function splitAlignedPeptides(peptideColumn: DG.Column<string>): DG.DataFrame {
  const splitter = WebLogo.getSplitterForColumn(peptideColumn);
  const colLen = peptideColumn.length;
  const resultDf = DG.DataFrame.create(colLen);
  let monomerList = splitter(peptideColumn.get(0)!);
  const columnList: DG.Column<string>[] = [];

  // create columns and fill the first row for faster values filling in the next loop
  for (let i = 0; i < monomerList.length; i++) {
    const col = resultDf.columns.addNewString((i + 1).toString());
    col.set(0, monomerList[i] || '-', false);
    columnList.push(col);
  }

  for (let rowIndex = 1; rowIndex < colLen; rowIndex++) {
    monomerList = splitter(peptideColumn.get(rowIndex)!);
    monomerList.forEach((monomer, colIndex) => columnList[colIndex].set(rowIndex, monomer || '-', false));
  }

  return resultDf;
} 
