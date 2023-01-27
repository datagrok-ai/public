import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {PtmFilter} from './custom-filters';
import {catchToLog, DataLoader, DataLoaderType, initPackageQueries, DataQueryDict} from './utils/data-loader';
import {DataLoaderFiles} from './utils/data-loader-files';
import {DataLoaderDb} from './utils/data-loader-db';
import {DataLoaderTest} from './tests/utils/data-loader-test';
import {MolecularLiabilityBrowser} from './molecular-liability-browser';
import {TreeBrowser} from './mlb-tree';
import {PtmFilterApp} from './apps/ptm-filter-app';
import {errorToConsole} from '@datagrok-libraries/utils/src/to-console';
import {MlbPtmFilter} from './mlb-ptm-filter';
import {checkVrForTreeUI} from './utils/check-vr-for-tree';
import {TaskBarProgressIndicator} from 'datagrok-api/dg';

export type PackagePropertiesType = { [name: string]: any };

export let _startInit: number;
export const _package = new DG.Package();
export let _properties: PackagePropertiesType;
export const dataPackageName: string = 'MolecularLiabilityBrowserData';
export const packageName: string = 'MolecularLiabilityBrowser';

let _mlbQueries: DataQueryDict;

/** DataLoader instance */
let dl: DataLoader | null;

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

async function getDataSourceSettings(): Promise<string> {
  return await catchToLog<Promise<string>>(
    'MLB: getDataSourceSettings',
    async () => {
      const funcList: DG.Func[] = DG.Func.find({package: dataPackageName, name: 'getPackageProperty'});
      if (funcList.length === 0) {
        return DataLoaderType.Test as string;
      } else {
        const call: DG.FuncCall = await funcList[0].prepare({propertyName: 'DataSource'}).call();
        const res = call.getOutputParamValue();
        return res;
      }
      // return grok.functions
      //   .call(`${dataPackageName}:getPackageProperty`, {propertyName: 'DataSource'});
    });
}

async function getMlbQueries(): Promise<DataQueryDict> {
  return await catchToLog(
    'MLB: package getMlbQueries()',
    async () => {
      return initPackageQueries(packageName);
    }
  );
}

//tags: init
export async function initMlb() {
  _startInit = Date.now();
  const pi = DG.TaskBarProgressIndicator.create('init MLB package');
  try {
    _mlbQueries = await getMlbQueries();
  } catch (err: any) {
    const msg: string = 'MolecularLiabilityBrowser package init error: ' +
      `${err instanceof Error ? err.message : (err as Object).toString()}`;
    grok.shell.error(msg);
    console.error(err);
  } finally {
    pi.close();
  }
}

/** Ensure dl is created and initialized */
async function initDataLoader(): Promise<void> {
  return await catchToLog<Promise<void>>(
    'MLB: initDataLoader',
    async () => {
      try {
        _properties = await catchToLog<Promise<PackagePropertiesType>>(
          'MLB: package getProperties',
          async () => { return _package.getProperties(); });

        if (!dl) {
          const serverListVersionDf: DG.DataFrame | undefined = _properties['IndexedDB'] ?
            await catchToLog( // this call checks database connection also
              'MLB: database query \'getListVersion\': ',
              async () => {
                const funcCall: DG.FuncCall = await _mlbQueries['getListVersion'].prepare().call();
                const res: DG.DataFrame = funcCall.getOutputParamValue() as DG.DataFrame;
                return res;
              }) : undefined;

          const dataSourceSettings: string = _properties['DataSource'];
          let resDl: DataLoader;
          switch (dataSourceSettings) {
          case DataLoaderType.Files:
            resDl = new DataLoaderFiles(_mlbQueries, serverListVersionDf);
            break;

          case DataLoaderType.Database:
            resDl = new DataLoaderDb(_mlbQueries, serverListVersionDf);
            break;

          case DataLoaderType.Test:
            resDl = new DataLoaderTest();
            break;

          default:
            throw new Error(`MLB: Unexpected data package property 'DataSource' value '${dataSourceSettings}'.`);
          }
          dl = resDl;

          // ptmFilter() function cannot be async and cannot await for any promise from init()
          await dl.init(_startInit);
        }
      } catch (err: any) {
        dl = null;
        throw err;
      }
    });
}

// -- Filters --

//name: PTM filter
//description: PTM filter
//tags: filter
//output: filter result
export function ptmFilter() {
  try {
    if (!dl)
      throw new Error(`MLB: PTM Filter: data loader is exist`);

    if (!(dl.predictedPtmMap && dl.predictedCdrMap && dl.observedPtmMap && dl.observedCdrMap))
      throw new Error(`MLB: PTM Filter: data loader is not initialized!`);

    const flt: PtmFilter = new PtmFilter(dl.predictedPtmMap, dl.predictedCdrMap, dl.observedPtmMap, dl.observedCdrMap);
    return flt;
  } catch (err: any) {
    const errMsg: string = 'MLB: ptmFilter() error:\n' + errorToConsole(err);
    console.error(errMsg);
    throw err;
  }
}

//name: mlbPtmFilter
//description: MLB PTM Filter
//tags: filter
//output: filter result
export function mlbPtmFilter(): DG.Filter {
  try {
    if (!dl)
      throw new Error(`MLB: mlbPtmFilter(), data loader is exist`);

    if (!(dl.predictedPtmMap && dl.predictedCdrMap && dl.observedPtmMap && dl.observedCdrMap))
      throw new Error(`MLB: mlbPtmFilter(), data loader is not initialized!`);

    const flt: MlbPtmFilter = new MlbPtmFilter();
    return flt;
  } catch (err: any) {
    const errMsg: string = 'MLB: mlbPtmFilter() error:\n' + errorToConsole(err);
    console.error(errMsg);
    throw err;
  }
}

// -- Apps --

//name: Molecular Liability Browser
//tags: app
export async function MolecularLiabilityBrowserApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('init MLB application');
  try {
    // Try to create data loader on every app start, if it is not created earlier
    await initDataLoader(); // at MolecularLiabilityBrowserApp()
    if (!dl)
      throw new Error(`MLB: MolecularLiabilityBrowser app: data loader does not exist`);

    const t1 = Date.now();

    console.debug('MLB: package.MolecularLiabilityBrowserApp()');
    grok.shell.windows.showToolbox = false;
    const urlParams: URLSearchParams = new URLSearchParams(window.location.search);

    const app = new MolecularLiabilityBrowser(dl);
    console.debug('MLB: package.MolecularLiabilityBrowserApp() before app init ' + `${fromStartInit()} s`);
    await app.init(urlParams);
    console.debug('MLB: package.MolecularLiabilityBrowserApp() after app init ' + `${fromStartInit()} s`);

    const t2 = Date.now();
    console.debug(`MLB: package.MolecularLiabilityBrowserApp(), ${((t2 - t1) / 1000).toString()} s`);
  } catch (err: unknown) {
    const msg: string = 'MolecularLiabilityBrowser app error: ' +
      `${err instanceof Error ? err.message : (err as Object).toString()}`;
    grok.shell.error(msg);
    console.error(msg);
  } finally {
    pi.close();
  }
}

// -- Apps for tests --

//name: ptmFilterApp
//description: Test/demo app for ptmFilterApp (ptmFilter)
export async function ptmFilterApp(): Promise<void> {
  const pi = DG.TaskBarProgressIndicator.create('open ptmFilter app');
  try {
    // Try to create data loader on every app start, if it is not created earlier
    await initDataLoader(); // at ptmFilterApp()
    if (!dl)
      throw new Error(`MLB: MolecularLiabilityBrowser app: data loader is exist`);

    const app = new PtmFilterApp(dl);
    await app.init();
  } catch (err: any) {
    const errMsg: string = 'MLB: ptmFilterApp() error:\n' + errorToConsole(err);
    grok.shell.error(errMsg);
    console.error(errMsg);
  } finally {
    pi.close();
  }
}

//name: dbGetListVersion
//description: Test query database
//output: dataframe result
export async function dbGetListVersion(): Promise<DG.DataFrame> {
  const t1: number = Date.now();
  const df: DG.DataFrame = await catchToLog( // this call checks database connection also
    'MLB: database query \'getListVersion\': ',
    async () => {
      const funcCall: DG.FuncCall = await _mlbQueries['getListVersion'].prepare().call();
      const res: DG.DataFrame = funcCall.getOutputParamValue() as DG.DataFrame;
      return res;
    });
  const t2: number = Date.now();
  console.debug(`MLB: dbGetListVersion() ET: ${((t2 - t1) / 1000)} s`);

  // grok.shell.addTableView(df);
  return df;
}

// --

/* WebLogo viewer is registered in Bio package */

/* VdRegions viewer is registered in Bio package */

//name: MlbTree
//description: Molecular Liability Browser clone tree viewer
//tags: viewer, panel
//output: viewer result
export function MlbTreeViewer() {
  return new TreeBrowser();
}

// -- Database checks --
//name: checkVrForTree
//description: Gets all trees and check for all VR ids in 'mlb' table
export async function checkVrForTree(): Promise<void> {
  const pi = TaskBarProgressIndicator.create('check VR for tree');
  try {
    await initDataLoader(); // at checkVrForTrees()
    if (!dl)
      throw new Error(`MLB: checkVrForTree() data loader does not exist`);

    await checkVrForTreeUI(dl!, pi);
  } catch (err: any) {
    const errMsg: string = errorToConsole(err);
    grok.shell.error(errMsg);
    console.error(errMsg);
  } finally {
    pi.close();
  }
}