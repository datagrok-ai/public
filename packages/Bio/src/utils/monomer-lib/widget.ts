/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {
  getUserLibSettings, LIB_PATH, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {MonomerLibManager} from './lib-manager';

import {MonomerLibFileManager} from './lib-file-manager';

export async function getLibraryPanelUI(): Promise<DG.Widget> {
  return new Widget().getWidget();
}

class Widget {
  constructor() { }

  private monomerLibFileManager: MonomerLibFileManager;

  async getWidget() {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance();
    const content = await this.getWidgetContent();
    return new DG.Widget(content);
  }

  private async getWidgetContent(): Promise<HTMLElement> {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance();
    const formHandler = new ControlsFormHandler();
    const libControlsForm = await formHandler.getInputsForm();
    $(libControlsForm).addClass('monomer-lib-controls-form');
    const addLibrariesBtn: HTMLButtonElement = ui.button('Add', async () => await this.addFiles());
    ui.tooltip.bind(addLibrariesBtn, 'Load new monomer libraries');
    const refreshIcon = ui.iconFA('sync-alt', async () => {
      await formHandler.refreshInputsForm();
    });
    const refreshBtn = ui.button(refreshIcon, () => {} );
    ui.tooltip.bind(refreshIcon, 'Refresh libraries list');
    const result = ui.divV([libControlsForm, ui.div([addLibrariesBtn, refreshBtn])]);
    return result;
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
}

class ControlsFormHandler {
  private monomerLibFileManager: MonomerLibFileManager;
  private inputsForm: HTMLDivElement;

  async getInputsForm(): Promise<HTMLDivElement> {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance();
    const controlList = await this.getControlList();
    const inputsForm = ui.form(controlList);

    const onFileListChange = this.monomerLibFileManager.onFileListChange;
    DG.debounce<void>(onFileListChange, 1000).subscribe(
      async () => await this.refreshInputsForm()
    );
    return inputsForm;
  }

  async refreshInputsForm(): Promise<void> {
    const pi = DG.TaskBarProgressIndicator.create('Updating monomer library list');
    // WARNING: this is necessary to prevent sync issues with the file system
    await this.monomerLibFileManager.refreshValidFilePaths();
    const updatedForm = await this.getInputsForm();
    $(this.inputsForm).replaceWith(updatedForm);
    grok.shell.info('Updated of monomer libraries');
    pi.close();
  }

  private async getControlList(): Promise<DG.InputBase<boolean | null>[]> {
    const settings = await getUserLibSettings();
    const fileManager = await MonomerLibFileManager.getInstance();
    await fileManager.refreshValidFilePaths();
    const libFileNameList: string[] = fileManager.getRelativePathsOfValidFiles();
    const libInputList: DG.InputBase<boolean | null>[] = [];

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
      libInputList.push(libInput);
    }
    return libInputList;
  }

  private getDeleteDialog(fileName: string): void {
    const dialog = ui.dialog('Warning');
    dialog.add(ui.divText(`Are you sure you want to delete ${fileName}?\nThis will delete the file from ${LIB_PATH}`))
      .onOK(async () => await this.monomerLibFileManager.deleteLibFile(fileName))
      .showModal(false);
  }
}
