/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {MonomerLibManager} from './lib-manager';

import {MonomerLibFileManager} from './lib-file-manager';

import * as rxjs from 'rxjs';

export async function getLibraryPanelUI(): Promise<DG.Widget> {
  return new Widget().getWidget();
}

class Widget {
  constructor() { }

  private onFileListChange = new rxjs.Subject<void>();

  private monomerLibFileManager: MonomerLibFileManager;

  async getWidget() {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance();
    const libChoiceInputsForm = await this.getInputs();
    const addLibrariesBtn: HTMLButtonElement = ui.button('Add', async () => await this.addFiles());
    return new DG.Widget(ui.divV([libChoiceInputsForm, ui.div([addLibrariesBtn])]));
  }

  private async getInputs(): Promise<HTMLDivElement> {
    const inputsForm = ui.form([]);
    const settings = await getUserLibSettings();
    const fileManager = await MonomerLibFileManager.getInstance();
    const libFileNameList: string[] = fileManager.getValidFiles();

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
            await MonomerLibManager.instance.loadLibraries(true); // from libraryPanel()
            grok.shell.info('Monomer library user settings saved.');
          });
        });
      const deleteIcon = ui.iconFA('trash-alt', async () => this.getDeleteDialog(libFileName));
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
        open: async (libFile) => {
          const content = await libFile.text();
          const name = libFile.name;
          try {
            await this.monomerLibFileManager.addLibFile(content, name);
          } catch (e) {
            grok.shell.error(`File ${name} is not a valid monomer library, verify it is aligned to HELM standard`);
          }
        },
      });
    });
    popup.separator();
    popup.item('Add custom', () => {
      const libFile = DG.Utils.openFile({
        accept: '.csv',
        open: async (libFile) => {
          // const content = await libFile.text();
          // const name = libFile.name;
          // await this.monomerLibFileManager.addCustomLibFile(content, name);
        },
      });
    });
    popup.show();
  }

  private getDeleteDialog(fileName: string): void {
    const dialog = ui.dialog('Warning');
    dialog.root.style.zIndex = '9999';
    dialog.add(ui.divText(`Are you sure you want to delete ${fileName}`))
      .onOK(async () => await this.monomerLibFileManager.deleteLibFile(fileName))
      .showModal(false);
  }
}
