import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {autostartOligoSdFileSubscription} from './autostart/registration';
import {defineAxolabsPattern} from './axolabs/define-pattern';
import {saveSenseAntiSense} from './structures-works/save-sense-antisense';
import {mainView} from './main/main-view';

export const _package = new DG.Package();


//name: Sequence Translator
//tags: app
export function sequenceTranslator(): void {
  const windows = grok.shell.windows;
  windows.showProperties = false;
  windows.showToolbox = false;
  windows.showHelp = false;

  const v = grok.shell.newView('Sequence Translator', [
    ui.tabControl({
      'MAIN': mainView(),
      'AXOLABS': defineAxolabsPattern(),
      'SDF': saveSenseAntiSense(),
    }),
  ]);
  v.box = true;
}

//tags: autostart
autostartOligoSdFileSubscription();
