import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DataLoaderBase, DataLoaderDB} from './utils/data-loader';

import {autostartOligoSdFileSubscription} from './autostart/registration';
import {OligoSdFileApp} from './apps/oligo-sd-file-app';
import {SequenceTranslatorUI} from './view/view';
import {LIB_PATH, DEFAULT_LIB_FILENAME} from './utils/const';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {MonomerWorks} from '@datagrok-libraries/bio/src/monomer-works/monomer-works';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

export class StPackage extends DG.Package {
  private _dataLoader?: DataLoaderBase;
  get dataLoader(): DataLoaderBase {
    if (!this._dataLoader)
      throw new Error('dataLoader is not initialized');
    return this._dataLoader!;
  };

  public async initDataLoader(): Promise<void> {
    if (this._dataLoader == undefined) {
      const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(
        'Initializing Sequence Translator data loader ...');
      try {
        const dl = new DataLoaderDB();
        await dl.init();
        this._dataLoader = dl;
      } catch (err: any) {
        const errMsg = err.hasOwnProperty('message') ? err.message : err.toString();
        grok.shell.error(errMsg);
        throw new Error('Initializing Sequence Translator data loader error: ' + errMsg);
      } finally {
        pi.close();
      }
    }
  }
}

export const _package: StPackage = new StPackage();

let monomerLib: IMonomerLib | null = null;
export let monomerWorks: MonomerWorks | null = null;

export function getMonomerWorks() {
  return monomerWorks;
}

export function getMonomerLib() {
  return monomerLib;
}


//name: Sequence Translator
//tags: app
export async function sequenceTranslator(): Promise<void> {
  const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create('Loading Sequence Translator app ...');
  try {
    try {
      const libHelper: IMonomerLibHelper = await getMonomerLibHelper();
      monomerLib = await libHelper.readLibrary(LIB_PATH, DEFAULT_LIB_FILENAME);
    } catch (err: any) {
      const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
      throw new Error('Loading monomer library error: ' + errMsg);
    }

    if (monomerWorks === null)
      monomerWorks = new MonomerWorks(monomerLib);

    const v = new SequenceTranslatorUI();
    await v.createLayout();
  } catch (err: any) {
    const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
    grok.shell.error(`Loading Sequence Translator application error: ` + errMsg);
    const k = 11;
    throw err;
  } finally {
    pi.close();
  }
}

//tags: autostart
export async function autostartST() {
  autostartOligoSdFileSubscription();
}

//name: oligoSdFileApp
//description: Test/demo app for oligoSdFile
export async function oligoSdFileApp() {
  const pi = DG.TaskBarProgressIndicator.create('open oligoSdFile app');
  try {
    grok.shell.windows.showProperties = false;
    grok.shell.windows.showHelp = false;

    const app = new OligoSdFileApp();
    await app.init();
  } finally {
    pi.close();
  }
}
