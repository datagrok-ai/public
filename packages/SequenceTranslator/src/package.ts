import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SequenceTranslatorUI} from './view/view';
import {LIB_PATH, DEFAULT_LIB_FILENAME} from './model/data-loading-utils/const';
import {IMonomerLib} from '@datagrok-libraries/bio/src/types';
import {getMonomerLibHelper, IMonomerLibHelper} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';
import {getJsonData} from './model/data-loading-utils/json-loader';

class StPackage extends DG.Package {
  private _monomerLib?: IMonomerLib;

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
}

export const _package: StPackage = new StPackage();

//name: Sequence Translator
//tags: app
export async function sequenceTranslatorApp(): Promise<void> {
  const pi: DG.TaskBarProgressIndicator = DG.TaskBarProgressIndicator.create('Loading Sequence Translator app ...');
  try {
    await getJsonData();
    await _package.initMonomerLib();
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

