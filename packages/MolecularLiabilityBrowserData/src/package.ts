/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();
export let _properties: any;

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//tags: init
export async function initMlbData(): Promise<void> {
  _properties = await _package.getProperties();
}

export async function initMlbDataHelper(): Promise<void> {

}

//name: getPackageProperty
//input: string propertyName
//output: string result
export function getPackageProperty(propertyName): string {
  const value: string = _properties[propertyName];
  return value;
}
