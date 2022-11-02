export function sortByStringLengthInDescendingOrder(array: string[]): string[] {
  return array.sort(function(a, b) {return b.length - a.length;});
}
  
export function stringify(items: string[]): string {
  return '["' + items.join('", "') + '"]';
}
