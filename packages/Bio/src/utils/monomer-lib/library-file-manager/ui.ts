/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';

import {
  getUserLibSettings, LIB_PATH, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {MonomerLibManager} from '../lib-manager';

import {MonomerLibFileManager} from './file-manager';
import {MonomerLibFileEventManager} from './event-manager';

export async function getLibraryPanelUI(): Promise<DG.Widget> {
  const eventManager = MonomerLibFileEventManager.getInstance();
  return new MonomerLibraryManagerWidget(eventManager).createWidget();
}

async function showManageLibrariesDialog() {
  const eventManager = MonomerLibFileEventManager.getInstance();
  const widget = await new MonomerLibraryManagerWidget(eventManager).createWidget();
  const dialog = ui.dialog('Manage monomer libraries');
  dialog.clear();
  dialog.add(widget);
  dialog.show();
}

export async function getMonomerLibraryManagerLink(): Promise<DG.Widget> {
  const link = ui.label('Monomer Library Manager');
  $(link).addClass('d4-link-action');
  link.onclick = async () => await showManageLibrariesDialog();
  return new DG.Widget(
    link
  );
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
    const formHandler = new LibraryControlsManager(this.eventManager);
    const libControlsForm = await formHandler.createControlsForm();
    $(libControlsForm).addClass('monomer-lib-controls-form');
    const addLibraryFilesButton: HTMLButtonElement = ui.button('Add', async () => await this.promptToAddLibraryFiles());
    ui.tooltip.bind(addLibraryFilesButton, 'Load new monomer libraries');
    const widgetContent = ui.divV([libControlsForm,
      ui.div([
        addLibraryFilesButton,
      ])
    ]);
    return widgetContent;
  }

  private async promptToAddLibraryFiles(): Promise<void> {
    DG.Utils.openFile({
      accept: '.json',
      open: async (selectedFile) => {
        const content = await selectedFile.text();
        const name = selectedFile.name;
        const progressIndicator = DG.TaskBarProgressIndicator.create(`Adding ${name} as a monomer library`);
        try {
          await this.monomerLibFileManager.addLibraryFile(content, name);
        } catch (e) {
          grok.shell.error(`File ${name} is not a valid monomer library, verify it is aligned to HELM JSON schema.`);
        } finally {
          progressIndicator.close();
        }
      },
    });
  }
}

class LibraryControlsManager {
  constructor(private eventManager: MonomerLibFileEventManager) {
    this.eventManager.updateUIControlsRequested$.subscribe(
      async () => await this.updateControlsForm()
    );
  }
  private monomerLibFileManager: MonomerLibFileManager;
  private inputsForm: HTMLDivElement;

  async createControlsForm(): Promise<HTMLElement> {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    const libraryControls = await this.createLibraryControls();
    const inputsForm = ui.form(libraryControls);

    return inputsForm;
  }

  async updateControlsForm(): Promise<void> {
    const updatedForm = await this.createControlsForm();
    $(this.inputsForm).replaceWith(updatedForm);
  }

  private async createLibraryControls(): Promise<DG.InputBase<boolean | null>[]> {
    const settings = await getUserLibSettings();
    const fileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    const libFileNameList: string[] = fileManager.getValidLibraryPaths();
    return libFileNameList.map((libFileName) => this.createLibInput(libFileName, settings));
  }

  private createLibInput(libFileName: string, settings: UserLibSettings): DG.InputBase<boolean | null> {
    const isMonomerLibrarySelected = !settings.exclude.includes(libFileName);
    const libInput = ui.boolInput(
      libFileName,
      isMonomerLibrarySelected,
      () => this.updateLibrarySelectionStatus(libInput, libFileName, settings)
    );
    const deleteIcon = ui.iconFA('trash-alt', () => this.promptForLibraryDeletion(libFileName));
    ui.tooltip.bind(deleteIcon, `Delete ${libFileName}`);
    libInput.addOptions(deleteIcon);
    return libInput;
  }

  private async updateLibrarySelectionStatus(
    libInput: DG.InputBase<boolean | null>,
    libFileName: string,
    settings: UserLibSettings
  ): Promise<void> {
    this.updateLibrarySettings(libInput.value, libFileName, settings);
    await setUserLibSettings(settings);
    await MonomerLibManager.instance.loadLibraries(true);
    grok.shell.info('Monomer library user settings saved.');
  }

  private updateLibrarySettings(
    isLibrarySelected: boolean | null,
    libFileName: string,
    settings: UserLibSettings
  ): void {
    if (isLibrarySelected) {
      // Remove selected library from exclusion list
      settings.exclude = settings.exclude.filter((libName) => libName !== libFileName);
    } else if (!settings.exclude.includes(libFileName)) {
      // Add unselected library to exclusion list
      settings.exclude.push(libFileName);
    }
  }

  private promptForLibraryDeletion(fileName: string): void {
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
