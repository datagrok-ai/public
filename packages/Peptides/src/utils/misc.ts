import * as DG from 'datagrok-api/dg';

import * as C from './constants';

export function stringToBool(str: string) {
  return str === 'true' ? true : false;
}

export function getSeparator(col: DG.Column): string {
  const separator = col.tags[C.TAGS.SEPARATOR];
  if (separator)
    return separator;

  const defaultSeparators = ['-', ' '];
  const categories = col.categories;
  for (const potentialSeparator of defaultSeparators) {
    if (categories.filter((sequence) => sequence.includes(potentialSeparator)).length)
      return potentialSeparator;
  }
  return separator ?? '';
}
