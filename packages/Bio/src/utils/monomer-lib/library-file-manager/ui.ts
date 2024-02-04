/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {
  getUserLibSettings, LIB_PATH, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {MonomerLibManager} from '../lib-manager';

import {MonomerLibFileManager} from './file-manager';
import {MonomerLibFileEventManager} from './event-manager';

export async function getLibraryPanelUI(): Promise<DG.Widget> {
  const eventManager = MonomerLibFileEventManager.getInstance();
  return new MonomerLibraryManagerWidget(eventManager).createWidget();
}

class MonomerLibraryManagerWidget {
  constructor(
    private eventManager: MonomerLibFileEventManager
  ) { }

  private monomerLibFileManager: MonomerLibFileManager;

  async createWidget() {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    const content = await this.getWidgetContent();
    return new DG.Widget(content);
  }

  private async getWidgetContent(): Promise<HTMLElement> {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    const formHandler = new ControlsFormManager(this.eventManager);
    const libControlsForm = await formHandler.createControlsForm();
    $(libControlsForm).addClass('monomer-lib-controls-form');
    const addLibraryFilesButton: HTMLButtonElement = ui.button('Add', async () => await this.addLibraryFiles());
    ui.tooltip.bind(addLibraryFilesButton, 'Load new monomer libraries');
    const refreshLibraryListIcon = ui.iconFA('sync-alt', async () => {
      const progressIndicator = DG.TaskBarProgressIndicator.create('Updating monomer library list');
      await formHandler.updateControlsForm();
      grok.shell.info('List of monomer libraries updated');
      progressIndicator.close();
    });
    const refreshLibraryListButton = ui.button(refreshLibraryListIcon, () => {} );
    ui.tooltip.bind(refreshLibraryListIcon, 'Refresh libraries list');
    const widgetContent = ui.divV([libControlsForm, ui.div([addLibraryFilesButton, refreshLibraryListButton])]);
    return widgetContent;
  }

  private async addLibraryFiles(): Promise<void> {
    const addFilesMenuPopup = DG.Menu.popup();
    addFilesMenuPopup.item('Add standard', () => {
      DG.Utils.openFile({
        accept: '.json',
        open: async (selectedFile) => {
          const content = await selectedFile.text();
          const name = selectedFile.name;
          const progressIndicator = DG.TaskBarProgressIndicator.create(`Adding ${name} as a monomer library`);
          try {
            await this.monomerLibFileManager.addLibraryFile(content, name);
          } catch (e) {
            grok.shell.error(`File ${name} is not a valid monomer library, verify it is aligned to HELM standard`);
          } finally {
            progressIndicator.close();
          }
        },
      });
    });
    addFilesMenuPopup.separator();
    addFilesMenuPopup.item('Add custom', () => {
      const libFile = DG.Utils.openFile({
        accept: '.csv',
        open: async (libFile) => {
          // const content = await libFile.text();
          // const name = libFile.name;
          // await this.monomerLibFileManager.addCustomLibFile(content, name);
        },
      });
    });
    addFilesMenuPopup.show();
  }
}

class ControlsFormManager {
  constructor(private eventManager: MonomerLibFileEventManager) { }
  private monomerLibFileManager: MonomerLibFileManager;
  private inputsForm: HTMLDivElement;

  async createControlsForm(): Promise<HTMLElement> {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    const controlList = await this.getControlList();
    const inputsForm = ui.form(controlList);

    this.eventManager.updateUIControlsRequested$.subscribe(
      async () => await this.updateControlsForm()
    );
    return inputsForm;
  }

  async updateControlsForm(): Promise<void> {
    // WARNING: this is necessary to prevent sync issues with the file system
    // await this.monomerLibFileManager.refreshValidFilePaths();
    const updatedForm = await this.createControlsForm();
    $(this.inputsForm).replaceWith(updatedForm);
  }

  private async getControlList(): Promise<DG.InputBase<boolean | null>[]> {
    const settings = await getUserLibSettings();
    const fileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    await fileManager.refreshLibraryFilePaths();
    const libFileNameList: string[] = fileManager.getRelativePathsOfValidLibraryFiles();
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
      const deleteIcon = ui.iconFA('trash-alt', async () => this.createLibraryFileDeletionDialog(libFileName));
      ui.tooltip.bind(deleteIcon, `Delete ${libFileName}`);
      libInput.addOptions(deleteIcon);
      libInputList.push(libInput);
    }
    return libInputList;
  }

  private createLibraryFileDeletionDialog(fileName: string): void {
    const dialog = ui.dialog('Warning');
    dialog.add(ui.divText(`Are you sure you want to delete ${fileName}?\nThis will delete the file from ${LIB_PATH}`))
      .onOK(async () => {
        const progressIndicator = DG.TaskBarProgressIndicator.create(`Deleting ${fileName} library`);
        await this.monomerLibFileManager.deleteLibraryFile(fileName);
        progressIndicator.close();
      })
      .showModal(false);
  }
}
