/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {GisViewer} from '../src/gis-viewer';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info('GIS Package info: ' +_package.webRoot);
}

//name: GISViewer
//description: GIS map viewer
//tags: viewer
//output: viewer result
export function gisViewer(): GisViewer {
  return new GisViewer();
}
