import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}
