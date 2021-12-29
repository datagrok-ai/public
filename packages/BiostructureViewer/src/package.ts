/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {BioStructureViewer } from './biostructure-viewer'
import { PdbEntry } from './pdb-entry';

export const _package = new DG.Package();

//name: BioStructure Viewer
//tags: app
export async function biostructureApp() {

  let pi = DG.TaskBarProgressIndicator.create('Opening BioStructure Viewer');
  let app = new BioStructureViewer();
  await app.init();
  pi.close();
}
