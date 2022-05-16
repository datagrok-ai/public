import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {MolecularLiabilityBrowser} from './molecular-liability-browser';
import {PtmFilter} from './custom-filters';
import {DataLoader} from './utils/data-loader';
import {DataLoaderFiles} from './utils/data-loader-files';

import {WebLogo} from '@datagrok-libraries/bio';

// import {DataLoaderDb} from './utils/data-loader-jnj';

export const _package = new DG.Package();

/** DataLoader instance
 */
let dl: DataLoader;

function getPathSegments(path: string) {
  const parser = document.createElement('a');
  parser.href = path;
  const pathSegments = parser.pathname.split('/');
  if (pathSegments.length > 4)
    return pathSegments[4];
  else
    return null;
}

//tags: init
export async function init() {
  const pi = DG.TaskBarProgressIndicator.create('Loading filters data...');

  dl = new DataLoaderFiles();
  // dl = new DataLoaderDb();
  await dl.init();
  pi.close();
}

//name: PTM filter
//description: PTM filter
//tags: filter
//output: filter result
export function ptmFilter() {
  if (!(dl.ptmMap && dl.cdrMap && dl.refDf))
    throw new Error(`Filter data is not initialized!`);

  return new PtmFilter(dl.ptmMap, dl.cdrMap, dl.refDf);
}

//name: Molecular Liability Browser
//tags: app
export async function MolecularLiabilityBrowserApp() {
  grok.shell.windows.showToolbox = false;
  const vid = getPathSegments(<string><unknown>window.location);

  const app = new MolecularLiabilityBrowser(dl);
  await app.init(vid);
}

//name: WebLogo
//tags: viewer, panel
//output: viewer result
export function webLogoViewer() {
  return new WebLogo();
}
