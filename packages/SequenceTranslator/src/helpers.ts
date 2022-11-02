export function sortByStringLengthInDescendingOrder(array: string[]): string[] {
  return array.sort(function(a, b) {return b.length - a.length;});
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
