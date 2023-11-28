/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {MonomerLibHelper} from './monomer-lib-helper';

import {getLibFileNameList} from './helpers';

export async function getLibraryPanelUI(): Promise<DG.Widget> {
  return new Widget().getWidget();
}

class Widget {
  async getWidget() {
    const libChoiceInputsForm = await this.getInputs();
    const addLibrariesBtn: HTMLButtonElement = ui.button('Add', this.addFiles);
    return new DG.Widget(ui.divV([libChoiceInputsForm, ui.div([addLibrariesBtn])]));
  }

  private async getInputs(): Promise<HTMLDivElement> {
    const inputsForm = ui.form([]);
    const settings = await getUserLibSettings();
    const libFileNameList: string[] = await getLibFileNameList();

    for (const libFileName of libFileNameList) {
      const libInput: DG.InputBase<boolean | null> = ui.boolInput(libFileName, !settings.exclude.includes(libFileName),
        () => {
          if (libInput.value == true) {
            // Checked library remove from excluded list
            settings.exclude = settings.exclude.filter((l) => l != libFileName);
          } else {
            // Unchecked library add to excluded list
            if (!settings.exclude.includes(libFileName))
              settings.exclude.push(libFileName);
          }
          setUserLibSettings(settings).then(async () => {
            await MonomerLibHelper.instance.loadLibraries(true); // from libraryPanel()
            grok.shell.info('Monomer library user settings saved.');
          });
        });
      const deleteIcon = ui.iconFA('trash-alt', () => this.getDeleteDialog(libFileName));
      ui.tooltip.bind(deleteIcon, `Delete ${libFileName}`);
      libInput.addOptions(deleteIcon);
      inputsForm.append(libInput.root);
    }
    return inputsForm;
  }

  private async addFiles(): Promise<void> {
    const popup = DG.Menu.popup();
    popup.item('Add standard', () => {
      const libFile = DG.Utils.openFile({
        accept: '.json',
        open: async (libFile) => { },
      });
    });
    popup.separator();
    popup.item('Add custom', () => {
      const libFile = DG.Utils.openFile({
        accept: '.csv',
        open: async (libFile) => { },
      });
    });
    popup.show();
  }

  private getDeleteDialog(fileName: string): void {
    const dialog = ui.dialog('Warning');
    dialog.add(ui.divText(`Are you sure you want to delete ${fileName}`))
      .onOK(() => {})
      .showModal(false);
  }
}
