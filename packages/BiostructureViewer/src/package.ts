/* Do not change these import lines to match external modules in webpack configuration */
//import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {BioStructureViewer} from './biostructure-viewer';
import {createViewer} from './viewers/molstar-viewer';

export const _package = new DG.Package();

//name: BioStructure Viewer
//tags: app
export async function biostructureApp() {
  const pi = DG.TaskBarProgressIndicator.create('Opening BioStructure Viewer');
  const app = new BioStructureViewer();
  await app.init();
  pi.close();
}

//name: Mol* BioStructure Viewer.
//tags: viewer
//input: string pdbID
export async function molstarView(pdbID: string) {
  const pi = DG.TaskBarProgressIndicator.create('Opening BioStructure* Viewer');
  await createViewer(pdbID);
  pi.close();
}
