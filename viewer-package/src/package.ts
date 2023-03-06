import {SecondViewer} from './second-viewer';
import {FirstViewer} from './first-viewer';
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

export const _package = new DG.Package();
// export const _functions: string[] = [];

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

//name: FirstViewer
//description: Creates FirstViewer viewer
//tags: viewer
//output: viewer result
export function _FirstViewer() {
  return new FirstViewer();
}

//name: SecondViewer
//description: Creates SecondViewer viewer
//tags: viewer
//output: viewer result
export function _SecondViewer() {
  return new SecondViewer();
}
