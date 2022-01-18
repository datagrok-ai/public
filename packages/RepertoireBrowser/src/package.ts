import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import { MolecularLiabilityBrowser } from './molecular-liability-browser';

export let _package = new DG.Package();

function getPathSegments(path: string) {
  let parser = document.createElement('a');
  parser.href = path;
  let pathSegments = parser.pathname.split('/');
  if (pathSegments.length > 4)
    return pathSegments[4];
  else
    return null;
}

//name: Repertoire Browser
//tags: app
export async function RepertoireBrowserApp() {
  grok.shell.windows.showToolbox = false;
  let vid = getPathSegments(<string><unknown>window.location);

  let pi = DG.TaskBarProgressIndicator.create('Opening Molecular Liability Browser');
  let app = new MolecularLiabilityBrowser();
  await app.init(vid);
  pi.close();
}
