/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
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

export async function tryCatch<T>(func: () => Promise<T>, finallyFunc?: () => any, callbackName: string = 'Oligo app'): Promise<T> {
  try {
    return await func();
  } catch (err: any) {
    const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
    grok.shell.error(`${callbackName} error: ` + errMsg);
    throw err;
  } finally {
    if (finallyFunc)
      finallyFunc();
  }
};

