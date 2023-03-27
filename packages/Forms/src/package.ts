/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {FormsViewer} from './forms';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: Forms
//description: Forms viewer
//tags: viewer
//output: viewer result
export function formsViewer() {
  return new FormsViewer();
}
