import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// datagrok libraries dependencies
import {IMonomerLib, MonomerWorks, readLibrary} from '@datagrok-libraries/bio';
import {DataLoaderBase, DataLoaderDB} from './utils/data-loader';

import {autostartOligoSdFileSubscription} from './autostart/registration';
import {OligoSdFileApp} from './apps/oligo-sd-file-app';
import {SequenceTranslatorUI} from './view/view';
import {LIB_PATH, DEFAULT_LIB_FILENAME} from './utils/const';

export const _package = new class extends DG.Package {
  private _dataLoader?: DataLoaderBase;
  get dataLoader(): DataLoaderBase {
    if (!this._dataLoader)
      throw new Error('dataLoader is not initialized');
    return this._dataLoader!;
  };

  public async initDataLoader(): Promise<void> {
    if (this._dataLoader == undefined) {
      const dl = new DataLoaderDB();
      await dl.init();
      this._dataLoader = dl;
    }
  }
}();

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
  monomerLib = await readLibrary(LIB_PATH, DEFAULT_LIB_FILENAME);

  if (monomerWorks === null)
    monomerWorks = new MonomerWorks(monomerLib);

  const v = new SequenceTranslatorUI();
  await v.createLayout();
}

//tags: autostart
export async function autostartST() {
  autostartOligoSdFileSubscription();
};

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
