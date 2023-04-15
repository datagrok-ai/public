import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DataLoaderBase, DataLoaderDB} from './model/data-loader/data-loader';

import {engageViewForOligoSdFileUI} from './model/registration/registration';
import {SequenceTranslatorUI} from './view/view';
import {LIB_PATH, DEFAULT_LIB_FILENAME} from './model/data-loader/const';
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
    if (this._dataLoader === undefined) {
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
export async function sequenceTranslatorApp(): Promise<void> {
  _package.logger.debug('ST: sequenceTranslatorApp() ');

  const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create('Loading Sequence Translator app ...');
  try {
    _package.initDataLoader().then(() => {}); // Do not wait for lists loaded from the database
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
    throw err;
  } finally {
    pi.close();
  }
}

//name: engageViewForOligoSdFile()
//description: Function to modify the view for SequenceTranslator registration
export async function engageViewForOligoSdFile(view: DG.TableView): Promise<void> {
  await engageViewForOligoSdFileUI(view);
}
