import * as ui from 'datagrok-api/ui';
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {PtmFilter} from './custom-filters';
import {catchToLog, DataLoader, DataLoaderType, initPackageQueries, DataQueryDict} from './utils/data-loader';
import {DataLoaderFiles} from './utils/data-loader-files';
import {DataLoaderDb} from './utils/data-loader-db';
import {DataLoaderTest} from './utils/data-loader-test';
import {MolecularLiabilityBrowser} from './molecular-liability-browser';
import {TreeBrowser} from './mlb-tree';

export let _startInit: number;
export const _package = new DG.Package();
export const dataPackageName: string = 'MolecularLiabilityBrowserData';
export const packageName: string = 'MolecularLiabilityBrowser';

let dataSourceSettings: string;
let mlbQueries: DataQueryDict;

/** DataLoader instance */
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
    'MLB: getMlbQueries',
    async () => {
      return initPackageQueries(packageName);
    }
  );
}

async function initDataLoader(): Promise<void> {
  return await catchToLog<Promise<void>>(
    'MLB: getDataLoader',
    async () => {
      try {
        if (!dl) {
          const [serverListVersionDf]: [DG.DataFrame] = await Promise.all([
            catchToLog( // this call checks database connection also
              'MLB: database query \'getListVersion\': ',
              async () => {
                const funcCall: DG.FuncCall = await mlbQueries['getListVersion'].prepare().call();
                const res: DG.DataFrame = funcCall.getOutputParamValue() as DG.DataFrame;
                return res;
              }),
          ]);
          // const serverListVersionDf: DG.DataFrame = null;

          switch (dataSourceSettings) {
          case DataLoaderType.Files:
            dl = new DataLoaderFiles(mlbQueries, serverListVersionDf);
            break;

          case DataLoaderType.Database:
            dl = new DataLoaderDb(mlbQueries, serverListVersionDf);
            break;

          case DataLoaderType.Test:
            dl = new DataLoaderTest();
            break;

          default:
            throw new Error(`MLB: Unexpected data package property 'DataSource' value '${dataSourceSettings}'.`);
          }

          // ptmFilter() function cannot be async and cannot await for any promise
          await dl.init(_startInit);
        }
      } catch (err: unknown) {
        dl = null;
        throw err;
      }
    });
}

//tags: init
export async function initMlb() {
  _startInit = Date.now();
  const pi = DG.TaskBarProgressIndicator.create('init MLB package');
  try {
    console.debug('MLB: initMlb() start, ' + `${fromStartInit()} s`);

    [dataSourceSettings, mlbQueries] = await Promise.all([
      // Package property is a long call as well as getting serverVersionDf from database
      getDataSourceSettings(),
      getMlbQueries()
    ]);

    let k = 11;

    // console.debug('MLB: initMlb() MLB-Data property + serverListVersion, ' + `${fromStartInit()}`);
    //
    //
    // console.debug(`MLB: initMLB() data loaded before init ${((Date.now() - _startInit) / 1000).toString()} s`);
    //
    // await dl.init(_startInit, serverListVersionDf);
    // console.debug(`MLB: initMLB() after init ${((Date.now() - _startInit) / 1000).toString()} s`);
  } catch (err: unknown) {
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
  if (!dl)
    throw new Error(`MLB: PTM Filter: data loader is exist`);

  if (!(dl.predictedPtmMap && dl.predictedCdrMap && dl.observedPtmMap && dl.observedCdrMap))
    throw new Error(`MLB: PTM Filter: data loader is not initialized!`);

  const flt: PtmFilter = new PtmFilter(dl.predictedPtmMap, dl.predictedCdrMap, dl.observedPtmMap, dl.observedCdrMap);
  return flt;
}

//name: Molecular Liability Browser
//tags: app
export async function MolecularLiabilityBrowserApp() {
  const pi = DG.TaskBarProgressIndicator.create('init MLB application');
  try {
    await initDataLoader();
    if (!dl)
      throw new Error(`MLB: MolecularLiabilityBrowser app: data loader is exist`);

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

/* WebLogo viewer is registered in Bio package */

/* VdRegions viewer is registered in Bio package */

//name: MlbTree
//description: Molecular Liability Browser clone tree viewer
//tags: viewer, panel
//output: viewer result
export function MlbTreeViewer() {
  return new TreeBrowser();
}
