import * as DG from 'datagrok-api/dg';

export function sortByReverseLength(array: string[]): string[] {
  return array.sort((a, b) => b.length - a.length);
}

export function download(name: string, href: string): void {
  const element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + href);
  element.setAttribute('download', name);
  element.click();
}
