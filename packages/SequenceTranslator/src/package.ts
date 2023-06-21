import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {autostartOligoSdFileSubscription} from './autostart/registration';
import {OligoSdFileApp} from './apps/oligo-sd-file-app';

// three tabs of the app's view
import {getMainTab} from './main-tab/main-tab';
import {getAxolabsTab} from './axolabs-tab/axolabs-tab';
import {getSdfTab} from './sdf-tab/sdf-tab';

const MAIN = 'MAIN';
const AXOLABS = 'AXOLABS';
const SDF = 'SDF';

const SEQUENCE_TRANSLATOR = 'Sequence Translator';
const DEFAULT_SEQUENCE = 'fAmCmGmAmCpsmU';
const DEFAULT_LIB_FILENAME = 'helmLib.json';

import {IMonomerLib, MonomerWorks, readLibrary} from '@datagrok-libraries/bio';

export const _package = new DG.Package();

const LIB_PATH = 'System:AppData/SequenceTranslator';

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

  const windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  let urlParams = new URLSearchParams(window.location.search);

  let mainSeq: string = DEFAULT_SEQUENCE;
  const view = grok.shell.newView(SEQUENCE_TRANSLATOR, []);
  view.box = true;

  const tabControl = ui.tabControl({
    [MAIN]: await getMainTab((seq) => {
      mainSeq = seq;
      urlParams = new URLSearchParams();
      urlParams.set('seq', mainSeq);
      updatePath();
    }),
    [AXOLABS]: getAxolabsTab(),
    [SDF]: getSdfTab(),
  });

  tabControl.onTabChanged.subscribe(() => {
    if (tabControl.currentPane.name !== MAIN)
      urlParams.delete('seq');
    else
      urlParams.set('seq', mainSeq);
    updatePath();
  });

  function updatePath() {
    const urlParamsTxt: string = Object.entries(urlParams)
      .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');
    view.path = '/apps/SequenceTranslator/SequenceTranslator' + `/${tabControl.currentPane.name}/?${urlParamsTxt}`;
  }

  const pathParts: string[] = window.location.pathname.split('/');
  if (pathParts.length >= 5) {
    const tabName: string = pathParts[5];
    tabControl.currentPane = tabControl.getPane(tabName);
  }

  view.append(tabControl);
  // console.debug('SequenceTranslator: app sequenceTranslator() ' + `view.path='${view.path}', view.basePath='${view.basePath}'.`);
}

//tags: autostart
export async function autostartST() {
  autostartOligoSdFileSubscription();
};

//name: oligoSdFileApp
//description: Test/demo app for oligoSdFile
export async function oligoSdFileApp() {
  // console.debug('SequenceTranslator: package.ts oligoSdFileApp()');

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
