import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {MolecularLiabilityBrowser} from './molecular-liability-browser';
import {MLBFilter} from './custom-filters';

export const _package = new DG.Package();
let ptmMap: {[key: string]: string};
let cdrMap: {[key: string]: string};
let referenceDf: DG.DataFrame;

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
  ptmMap = JSON.parse(await _package.files.readAsText('ptm_map.json'));
  cdrMap = JSON.parse(await _package.files.readAsText('cdr_map.json'));
  referenceDf = (await _package.files.readBinaryDataFrames(`ptm_in_cdr.d42`))[0];
  pi.close();
}

//name: MLB Filter
//description: MLB Filter
//tags: filter
//output: filter result
export function mlbFilter() {
  if (!(ptmMap && cdrMap && referenceDf))
    throw new Error(`Filter data is not initialized!`);

  return new MLBFilter(ptmMap, cdrMap, referenceDf);
}

//name: Molecular Liability Browser
//tags: app
export async function MolecularLiabilityBrowserApp() {
  grok.shell.windows.showToolbox = false;
  const vid = getPathSegments(<string><unknown>window.location);

  const app = new MolecularLiabilityBrowser();
  await app.init(vid);
}
