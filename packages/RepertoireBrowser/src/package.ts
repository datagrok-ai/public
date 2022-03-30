import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {MolecularLiabilityBrowser} from './molecular-liability-browser';
import {NucleotidesWebLogo, AminoacidsWebLogo} from './viewers/web-logo';

export const _package = new DG.Package();

function getPathSegments(path: string) {
  const parser = document.createElement('a');
  parser.href = path;
  const pathSegments = parser.pathname.split('/');
  if (pathSegments.length > 4)
    return pathSegments[4];
  else
    return null;
}

//name: Repertoire Browser
//tags: app
export async function RepertoireBrowserApp() {
  grok.shell.windows.showToolbox = false;
  const vid = getPathSegments(<string><unknown>window.location);

  const pi = DG.TaskBarProgressIndicator.create('Opening Molecular Liability Browser');
  const app = new MolecularLiabilityBrowser();
  await app.init(vid);
  pi.close();
}

//name: NucleotidesWebLogo
//tags: viewer,panel
//output: viewer result
export function nucleotidesWebLogoViewer() {
  return new NucleotidesWebLogo();
}

//name: AminoacidsWebLogo
//tags: viewer,panel
//output: viewer result
export function aminoacidsWebLogoViewer() {
  return new AminoacidsWebLogo();
}
