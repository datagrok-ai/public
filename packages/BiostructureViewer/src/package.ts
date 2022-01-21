/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
//import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

//@ts-ignore
import {Viewer} from 'molstar/build/viewer/molstar';

import {BioStructureViewer} from './biostructure-viewer';

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
//tags: app
export async function molstarView() {
  const view = grok.shell.newView('Mol*');
  const viewerContainer = view.root;
  const v = new Viewer(viewerContainer);//Viewer();
  //await v.loadPdb('2v0a');
  v.handleResize();
}
