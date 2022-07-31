import * as DG from 'datagrok-api/dg';
import {map, normalizedObj} from './map';
import {COL_NAMES} from './additional-modifications';

export const UNITS = {
  OPTICAL_DENSITY: 'OD',
  MICRO_GRAM: 'µg',
  MILLI_GRAM: 'mg',
  MICRO_MOLE: 'µmole',
  NANO_MOLE: 'nmole',
};

export function normalizeSequence(sequence: string, synthesizer: string | null, technology: string | null,
  additionalModsDf: DG.DataFrame): string {
  const additionalCodesCol = additionalModsDf.getCol(COL_NAMES.ABBREVIATION);
  const baseModifsCol = additionalModsDf.getCol(COL_NAMES.BASE_MODIFICATION);

  const codes = (technology == null) ?
    getAllCodesOfSynthesizer(synthesizer!).concat(additionalCodesCol.toList()) :
    Object.keys(map[synthesizer!][technology]);

  for (let i = 0; i < additionalModsDf.rowCount; i++)
    normalizedObj[additionalCodesCol.getString(i)] = (baseModifsCol.get(i) != 'NO') ? baseModifsCol.get(i) : '';

  for (let i = 0; i < codes.length; i++) {
    codes[i] = codes[i].replace('(', '\\(');
    codes[i] = codes[i].replace(')', '\\)');
  }

  const sortedCodes = sortByStringLengthInDescOrder(codes);
  const regExp = new RegExp('(' + sortedCodes.join('|') + ')', 'g');
  return sequence.replace(regExp, function(code) {return normalizedObj[code];});
}

export function getAllCodesOfSynthesizer(synthesizer: string): string[] {
  let codes: string[] = [];
  for (const technology of Object.keys(map[synthesizer]))
    codes = codes.concat(Object.keys(map[synthesizer][technology]));
  return codes;
}

export function deleteWord(sequence: string, searchTerm: string): string {
  let n = sequence.search(searchTerm);
  while (sequence.search(searchTerm) > -1) {
    n = sequence.search(searchTerm);
    sequence = sequence.substring(0, n) + sequence.substring(n + searchTerm.length, sequence.length);
  }
  return sequence;
}

export function saveAsCsv(table: DG.DataFrame): void {
  const link = document.createElement('a');
  link.setAttribute('href', 'data:text/csv;charset=utf-8,\uFEFF' + encodeURI(table.toCsv()));
  link.setAttribute('download', 'Oligo Properties.csv');
  link.click();
}

export function sortByStringLengthInDescOrder(array: string[]): string[] {
  return array.sort(function(a, b) {return b.length - a.length;});
}

export function mergeOptions(obj1: {[ind: string]: number}, obj2: {[ind: string]: number}): {[ind: string]: number} {
  const obj3: {[index: string]: number} = {};
  for (const attrname in obj1) {
    if (Object.prototype.hasOwnProperty.call(obj1, attrname))
      obj3[attrname] = obj1[attrname];
  }
  for (const attrname in obj2) {
    if (Object.prototype.hasOwnProperty.call(obj2, attrname))
      obj3[attrname] = obj2[attrname];
  }
  return obj3;
}
