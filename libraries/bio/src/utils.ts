import * as DG from 'datagrok-api/dg';

export const SEPARATOR = 'sequence-separator';

export function getColumnSeparator(col: DG.Column<string>): string {
  const separator: string | undefined = col.tags[SEPARATOR];
  if (separator)
    return separator;

  const defaultSeparators = ['.', '-', ' '];
  const categories = col.categories;
  const catLen = categories.length;
  for (const potentialSeparator of defaultSeparators) {
    if (categories.filter(sequence => sequence.includes(potentialSeparator) || sequence == '').length == catLen)
      return potentialSeparator;
  }
  return separator ?? '';
}
