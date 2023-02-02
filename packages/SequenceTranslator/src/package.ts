import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

// datagrok libraries dependencies
import {IMonomerLib, MonomerWorks, readLibrary} from '@datagrok-libraries/bio';

// inner dependencies
import {autostartOligoSdFileSubscription} from './autostart/registration';
import {OligoSdFileApp} from './apps/oligo-sd-file-app';
import {SequenceTranslatorUI} from './view/view';
import {LIB_PATH, DEFAULT_LIB_FILENAME} from './utils/const';

export const _package = new DG.Package();

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
  const pathParts: string[] = window.location.pathname.split('/');
  // here the value is '' and ''

  if (monomerWorks === null)
    monomerWorks = new MonomerWorks(monomerLib);

  // let urlParams = new URLSearchParams(window.location.search);

  // let mainSeq: string = DEFAULT_SEQUENCE;
  // const view = grok.shell.newView(SEQUENCE_TRANSLATOR, []);
  // view.box = true;

  // const tabControl = ui.tabControl({
  //   [MAIN]: await getMainTab((seq) => {
  //     mainSeq = seq;
  //     urlParams = new URLSearchParams();
  //     urlParams.set('seq', mainSeq);
  //     updatePath();
  //   }),
  //   [AXOLABS]: getAxolabsTab(),
  //   [SDF]: getSdfTab(),
  // });

  // const sdfPane = tabControl.getPane(SDF);
  // ui.tooltip.bind(sdfPane.header, 'Get atomic-level structure for SS + AS/AS2 and save SDF');

  // tabControl.onTabChanged.subscribe(() => {
  //   if (tabControl.currentPane.name !== MAIN)
  //     urlParams.delete('seq');
  //   else
  //     urlParams.set('seq', mainSeq);
  //   updatePath();
  // });

  // function updatePath() {
  //   const urlParamsTxt: string = Object.entries(urlParams)
  //     .map(([key, value]) => `${key}=${encodeURIComponent(value)}`).join('&');
  //   view.path = '/apps/SequenceTranslator' + `/${tabControl.currentPane.name}/?${urlParamsTxt}`;
  // }

  // if (pathParts.length >= 4) {
  //   const tabName: string = pathParts[3];
  //   tabControl.currentPane = tabControl.getPane(tabName);
  // }

  // view.append(tabControl);
  
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
