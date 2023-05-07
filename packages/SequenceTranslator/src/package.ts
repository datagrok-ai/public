import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {DBLoaderBase, DataLoaderDB} from './model/database-utils/loader';

import {engageViewForOligoSdFileUI} from './model/registration/registration';
import {SequenceTranslatorUI} from './view/view';
import {LIB_PATH, DEFAULT_LIB_FILENAME} from './model/database-utils/const';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {FormatDictionaryLoader} from './model/translation-utils/dictionary-loader';

class StPackage extends DG.Package {
  private _dbLoader?: DBLoaderBase;
  private _monomerLib?: IMonomerLib;

  get dbLoader(): DBLoaderBase {
    if (!this._dbLoader)
      throw new Error('ST: dataLoader is not initialized');
    return this._dbLoader!;
  };

  get monomerLib(): IMonomerLib {
    if (!this._monomerLib)
      throw new Error ('ST: monomer lib not loaded')
    return this._monomerLib!;
  }

  public async initMonomerLib(): Promise<void> {
    if (this._monomerLib === undefined) {
      const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(
        'Initializing Sequence Translator monomer library ...');
      try {
        const libHelper: IMonomerLibHelper = await getMonomerLibHelper();
        this._monomerLib = await libHelper.readLibrary(LIB_PATH, DEFAULT_LIB_FILENAME);
      } catch (err: any) {
        const errMsg: string = err.hasOwnProperty('message') ? err.message : err.toString();
        throw new Error('ST: Loading monomer library error: ' + errMsg);
      } finally {
        pi.close();
      }
    }
  }

  public async initDBLoader(): Promise<void> {
    if (this._dbLoader === undefined) {
      const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create(
        'Initializing Sequence Translator data loader ...');
      try {
        const dl = new DataLoaderDB();
        await dl.init();
        this._dbLoader = dl;
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

//tags: init
export async function initSequenceTranslator(): Promise<void> {
  await FormatDictionaryLoader.getInstance().init();
  await _package.initMonomerLib();
  _package.initDBLoader().then(() => {}); // Do not wait for lists loaded from the database
}

//name: Sequence Translator
//tags: app
export async function sequenceTranslatorApp(): Promise<void> {
  const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create('Loading Sequence Translator app ...');
  try {
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
