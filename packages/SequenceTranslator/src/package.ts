import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {autostartOligoSdFileSubscription} from './autostart/registration';
import {defineAxolabsPattern} from './axolabs/define-pattern';
import {saveSenseAntiSense} from './structures-works/save-sense-antisense';
import {mainView} from './main/main-view';
import {IMonomerLib, MonomerWorks, readLibrary} from '@datagrok-libraries/bio';

export const _package = new DG.Package();

const LIB_PATH = 'System:AppData/SequenceTranslator';

let monomerLib: IMonomerLib | null = null;
export let monomerWorks: MonomerWorks | null = null;

export function getMonomerWorks() {
  return monomerWorks;
};

//name: Sequence Translator
//tags: app
export async function sequenceTranslator(): Promise<void> {
  monomerLib = await readLibrary(LIB_PATH, 'helmLib.json');

  if (monomerWorks == null)
    monomerWorks = new MonomerWorks(monomerLib);

  const windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  const v = grok.shell.newView('Sequence Translator', []);
  v.box = true;
  v.append(ui.tabControl({
    'MAIN': await mainView(),
    'AXOLABS': defineAxolabsPattern(),
    'SDF': saveSenseAntiSense(),
  }));
}

//tags: autostart
autostartOligoSdFileSubscription();
