import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {autostartOligoSdFileSubscription} from './autostart/registration';
import {defineAxolabsPattern} from './axolabs/define-pattern';
import {saveSenseAntiSense} from './structures-works/save-sense-antisense';
import {mainView} from './main/main-view';
import {IMonomerLib, MonomerWorks, readLibrary} from '@datagrok-libraries/bio';
import {OligoSdFileApp} from './apps/oligo-sd-file-app';

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
  monomerLib = await readLibrary(LIB_PATH, 'helmLib.json');

  if (monomerWorks == null)
    monomerWorks = new MonomerWorks(monomerLib);

  const windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  let urlParams: URLSearchParams = new URLSearchParams(window.location.search);
  let mainSeq: string = 'fAmCmGmAmCpsmU';
  const v = grok.shell.newView('Sequence Translator', []);
  v.box = true;
  const tc = ui.tabControl({
    'MAIN': await mainView((seq) => {
      mainSeq = seq;
      urlParams = new URLSearchParams();
      urlParams.set('seq', mainSeq);
      updatePath();
    }),
    'AXOLABS': defineAxolabsPattern(),
    'SDF': saveSenseAntiSense(),
  });
  tc.onTabChanged.subscribe((value) => {
    if (tc.currentPane.name != 'MAIN')
      urlParams.delete('seq');
    else
      urlParams.set('seq', mainSeq);
    updatePath();
  });

  function updatePath() {
    const urlParamsTxt: string = Object.entries(urlParams)
      .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');
    v.path = '/apps/SequenceTranslator/SequenceTranslator' + `/${tc.currentPane.name}/?${urlParamsTxt}`;
  }

  const pathParts: string[] = window.location.pathname.split('/');
  if (pathParts.length >= 5) {
    const tabName: string = pathParts[5];
    tc.currentPane = tc.getPane(tabName);
  }

  v.append(tc);
  console.debug('SequenceTranslator: app sequenceTranslator() ' + `v.path='${v.path}', v.basePath='${v.basePath}'.`);
}

//tags: autostart
autostartOligoSdFileSubscription();

//name: oligoSdFileApp
//description: Test/demo app for oligoSdFile
export async function oligoSdFileApp() {
  console.debug('SequenceTranslator: package.ts oligoSdFileApp()');

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