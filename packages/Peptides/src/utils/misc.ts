import * as DG from 'datagrok-api/dg';

import * as C from './constants';

export function stringToBool(str: string): boolean {
  return str === 'true' ? true : false;
}

export function getSeparator(col: DG.Column): string {
  const separator = col.tags[C.TAGS.SEPARATOR];
  if (separator)
    return separator as string;

  const defaultSeparators = ['-', ' '];
  const categories = col.categories;
  for (const potentialSeparator of defaultSeparators) {
    if (categories.filter((sequence) => sequence.includes(potentialSeparator)).length)
      return potentialSeparator;
  }
  return separator as string ?? '';
}

export function getTypedArrayConstructor(
  maxNum: number): Uint8ArrayConstructor | Uint16ArrayConstructor | Uint32ArrayConstructor {
  return maxNum < 256 ? Uint8Array :
    maxNum < 65536 ? Uint16Array :
      Uint32Array;
}
