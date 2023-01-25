import * as DG from 'datagrok-api/dg';

export function sortByStringLengthInDescendingOrder(array: string[]): string[] {
  return array.sort(function(a, b) { return b.length - a.length; });
}

export function stringify(items: string[]): string {
  return '["' + items.join('", "') + '"]';
}

export function download(name: string, href: string): void {
  const element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + href);
  element.setAttribute('download', name);
  element.click();
}

export function removeEmptyRows(t: DG.DataFrame, colToCheck: DG.Column): DG.DataFrame {
  for (let i = t.rowCount - 1; i > -1; i--) {
    if (colToCheck.getString(i) === '')
      t.rows.removeAt(i, 1, false);
  }
  return t;
}

export function differenceOfTwoArrays(a: string[], b: string[]): string[] {
  return a.filter((x) => !b.includes(x));
}
