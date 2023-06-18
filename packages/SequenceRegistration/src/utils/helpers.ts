/** Builds JSON-like string from string array and doubles back-slash specifically for DG.Column '.choices' tag */
export function stringify(items: string[]): string {
  return '["' + items.map((v) => v.replace('\\', '\\\\')).join('", "') + '"]';
}

export function download(name: string, href: string): void {
  const element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + href);
  element.setAttribute('download', name);
  element.click();
}
