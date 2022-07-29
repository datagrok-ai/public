import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {PtmFilter} from './custom-filters';
import {catchToLog, DataLoader, DataLoaderType} from './utils/data-loader';
import {DataLoaderFiles} from './utils/data-loader-files';
import {DataLoaderDb} from './utils/data-loader-db';
import {MolecularLiabilityBrowser} from './molecular-liability-browser';
import {TreeBrowser} from './mlb-tree';

export let _startInit: number;
export const _package = new DG.Package();
const dataPackageName: string = 'MolecularLiabilityBrowserData';
const packageName: string = 'MolecularLiabilityBrowser';

/** DataLoader instance
 */
let dl: DataLoader;

// function getPathSegments(path: string) {
//   const parser = document.createElement('a');
//   parser.href = path;
//   const pathSegments = parser.pathname.split('/');
//   if (pathSegments.length > 4)
//     return pathSegments[4];
//   else
//     return null;
// }

function fromStartInit(): string {
  return ((Date.now() - _startInit) / 1000).toString();
}

//tags: init
export async function initMlb() {
  _startInit = Date.now();
  const pi = DG.TaskBarProgressIndicator.create('MLB: initMlb() Loading filters data...');
  try {
    // const dataSourceSettings: string = await grok.functions.call(`${dataPackageName}:getPackageProperty`,
    //   {propertyName: 'DataSource'});

    console.debug('MLB: initMlb() start, ' + `${fromStartInit()} s`);
    // Package property is a long call as well as getting serverVersionDf from database
    const [dataSourceSettings, serverListVersionDf]: [string, DG.DataFrame] = await Promise.all([
      grok.functions.call(`${dataPackageName}:getPackageProperty`, {propertyName: 'DataSource'}),
      catchToLog( // this call checks database connection also
        'MLB database error \'getListVersion\': ',
        () => { return grok.functions.call(`${packageName}:getListVersion`); }),
    ]);
    console.debug('MLB: initMlb() MLB-Data property + serverListVersion, ' + `${fromStartInit()}`);

    switch (dataSourceSettings) {
    case DataLoaderType.Files:
      dl = new DataLoaderFiles();
      break;
    case DataLoaderType.Database:
      dl = new DataLoaderDb();
      break;
    default:
      throw new Error(`MLB: Unexpected data package property 'DataSource' value '${dataSourceSettings}'.`);
    }
    console.debug(`MLB: initMLB() data loaded before init ${((Date.now() - _startInit) / 1000).toString()} s`);

    await dl.init(_startInit, serverListVersionDf);
    console.debug(`MLB: initMLB() after init ${((Date.now() - _startInit) / 1000).toString()} s`);
  } catch (err: unknown) {
    dl = null;
    const msg: string = 'MolecularLiabilityBrowser package init error: ' +
      `${err instanceof Error ? err.message : (err as Object).toString()}`;
    grok.shell.error(msg);
    console.error(err);
  } finally {
    pi.close();
  }
}

//name: PTM filter
//description: PTM filter
//tags: filter
//output: filter result
export function ptmFilter() {
  if (!(dl.ptmMap && dl.cdrMap && dl.refDf))
    throw new Error(`MLB: Filter data is not initialized!`);

  const flt: PtmFilter = new PtmFilter(dl.ptmMap, dl.cdrMap, dl.refDf);
  return flt;
}

//name: Molecular Liability Browser
//tags: app
export async function MolecularLiabilityBrowserApp() {
  try {
    if (!dl)
      throw new Error('Data loader is not initialized.');

    console.debug('MLB.package.MolecularLiabilityBrowserApp()');
    grok.shell.windows.showToolbox = false;
    const urlParams: URLSearchParams = new URLSearchParams(window.location.search);

    const app = new MolecularLiabilityBrowser(dl);
    console.debug(`MLB.package.MolecularLiabilityBrowserApp() before ` +
      `app init ${((Date.now() - _startInit) / 1000).toString()} s`);
    await app.init(urlParams);
    console.debug(`MLB.package.MolecularLiabilityBrowserApp() after ` +
      `app.init ${((Date.now() - _startInit) / 1000).toString()} s`);
  } catch (err: unknown) {
    const msg: string = 'MolecularLiabilityBrowser app error: ' +
      `${err instanceof Error ? err.message : (err as Object).toString()}`;
    grok.shell.error(msg);
    console.error(msg);
  }
}

/* WebLogo viewer is registered in Bio package */

/* VdRegions viewer is registered in Bio package */

//name: MlbTree
//description: Molecular Liability Browser clone tree viewer
//tags: viewer, panel
//output: viewer result
export function MlbTreeViewer() {
  return new TreeBrowser();
}
