/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {MonomerLibHelper} from './monomer-lib-helper';

import {manageFiles, getLibFileNameList} from './helpers';

export async function getLibraryPanelUI(): Promise<DG.Widget> {
  const filesButton: HTMLButtonElement = ui.button('Manage', manageFiles);
  const inputsForm: HTMLDivElement = ui.inputs([]);
  const libFileNameList: string[] = await getLibFileNameList();

  const settings = await getUserLibSettings();

  for (const libFileName of libFileNameList) {
    const libInput: DG.InputBase<boolean | null> = ui.boolInput(libFileName, !settings.exclude.includes(libFileName),
      () => {
        if (libInput.value == true) {
          // Checked library remove from excluded list
          settings.exclude = settings.exclude.filter((l) => l != libFileName);
        } else {
          // Unchecked library add to excluded list
          if (!settings.exclude.includes(libFileName)) settings.exclude.push(libFileName);
        }
        setUserLibSettings(settings).then(async () => {
          await MonomerLibHelper.instance.loadLibraries(true); // from libraryPanel()
          grok.shell.info('Monomer library user settings saved.');
        });
      });
    inputsForm.append(libInput.root);
  }
  return new DG.Widget(ui.divV([inputsForm, ui.div(filesButton)]));
}
