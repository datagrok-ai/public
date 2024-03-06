/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import * as rxjs from 'rxjs';
import './style.css';

import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {MonomerLibManager} from '../lib-manager';

import {MonomerLibFileManager} from './file-manager';
import {MonomerLibFileEventManager} from './event-manager';

export async function showManageLibrariesDialog(): Promise<void> {
  await DialogWrapper.showDialog();
}

export async function getMonomerLibraryManagerLink(): Promise<DG.Widget> {
  const link = ui.label('Manage monomer libraries');
  $(link).addClass('d4-link-action');
  link.onclick = async () => await showManageLibrariesDialog();
  return new DG.Widget(
    link
  );
}

class MonomerLibraryManagerWidget {
  private constructor(
    private eventManager: MonomerLibFileEventManager
  ) { }

  private static _instance: MonomerLibraryManagerWidget;

  static async getContent(eventManager: MonomerLibFileEventManager): Promise<DG.Widget> {
    if (!MonomerLibraryManagerWidget._instance)
      MonomerLibraryManagerWidget._instance = new MonomerLibraryManagerWidget(eventManager);

    if (!MonomerLibraryManagerWidget._instance.widget)
      MonomerLibraryManagerWidget._instance.widget = await MonomerLibraryManagerWidget._instance.createWidget();

    return MonomerLibraryManagerWidget._instance.widget;
  }

  private monomerLibFileManager: MonomerLibFileManager;
  private widget: DG.Widget | undefined;

  private async createWidget() {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    const content = await this.getWidgetContent();
    this.eventManager.addLibraryFileRequested$.subscribe(
      async () => await this.promptToAddLibraryFiles()
    );
    return new DG.Widget(content);
  }

  private async getWidgetContent(): Promise<HTMLElement> {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    const libControlsForm = await LibraryControlsManager.createControlsForm(this.eventManager);
    $(libControlsForm).addClass('monomer-lib-controls-form');
    const widgetContent = ui.divV([libControlsForm]);
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
          // this.eventManager.updateLibrarySelectionStatus(name, true);
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
  private constructor(private eventManager: MonomerLibFileEventManager) {
    this.eventManager.updateUIControlsRequested$.subscribe(
      async () => await this.updateControlsForm()
    );
    this.eventManager.librarySelectionRequested$.subscribe(
      async ([fileName, isSelected]) => await this.updateLibrarySelectionStatus(isSelected, fileName)
    );
  }
  private monomerLibFileManager: MonomerLibFileManager;
  private userLibSettings: UserLibSettings;

  static async createControlsForm(eventManager: MonomerLibFileEventManager): Promise<HTMLElement> {
    const manager = new LibraryControlsManager(eventManager);
    await manager.initialize();

    return await manager._createControlsForm();
  }

  private async _createControlsForm(): Promise<HTMLElement> {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    const libraryControls = await this.createLibraryControls();
    const inputsForm = ui.form(libraryControls);
    $(inputsForm).addClass('monomer-lib-controls-form');

    return inputsForm;
  }

  private async initialize(): Promise<void> {
    this.userLibSettings = await getUserLibSettings();
  };

  private async updateControlsForm(): Promise<void> {
    const updatedForm = await this._createControlsForm();
    $('.monomer-lib-controls-form').replaceWith(updatedForm);
  }

  private async createLibraryControls(): Promise<DG.InputBase<boolean | null>[]> {
    const fileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    const libFileNameList: string[] = fileManager.getValidLibraryPaths();
    return libFileNameList.map((libFileName) => this.createLibInput(libFileName));
  }

  private createLibInput(libFileName: string): DG.InputBase<boolean | null> {
    const isMonomerLibrarySelected = !this.userLibSettings.exclude.includes(libFileName);
    const libInput = ui.boolInput(
      libFileName,
      isMonomerLibrarySelected,
      (isSelected: boolean) => this.eventManager.updateLibrarySelectionStatus(libFileName, isSelected)
    );
    ui.tooltip.bind(libInput.root, `Include monomers from ${libFileName}`);
    const deleteIcon = ui.iconFA('trash-alt', () => this.promptForLibraryDeletion(libFileName));
    ui.tooltip.bind(deleteIcon, `Delete ${libFileName}`);
    libInput.addOptions(deleteIcon);
    return libInput;
  }

  private async updateLibrarySelectionStatus(
    isMonomerLibrarySelected: boolean,
    libFileName: string
  ): Promise<void> {
    this.updateLibrarySettings(isMonomerLibrarySelected, libFileName);
    await setUserLibSettings(this.userLibSettings);
    await MonomerLibManager.instance.loadLibraries(true);
    grok.shell.info('Monomer library user settings saved');
  }

  private updateLibrarySettings(
    isLibrarySelected: boolean | null,
    libFileName: string,
  ): void {
    if (isLibrarySelected) {
      // Remove selected library from exclusion list
      this.userLibSettings.exclude = this.userLibSettings.exclude.filter((libName) => libName !== libFileName);
    } else if (!this.userLibSettings.exclude.includes(libFileName)) {
      // Add unselected library to exclusion list
      this.userLibSettings.exclude.push(libFileName);
    }
  }

  private promptForLibraryDeletion(fileName: string): void {
    const dialog = ui.dialog('Warning');
    dialog.add(ui.divText(`Delete file ${fileName}?`))
      .onOK(async () => {
        try {
          const progressIndicator = DG.TaskBarProgressIndicator.create(`Deleting ${fileName} library`);
          this.updateLibrarySelectionStatus(false, fileName);
          await this.monomerLibFileManager.deleteLibraryFile(fileName);
          progressIndicator.close();
        } catch (e) {
          console.error(e);
          grok.shell.error(`Failed to delete ${fileName} library`);
        }
      })
      .showModal(false);
  }
}

class DialogWrapper {
  private constructor() { }

  private static _instance: DialogWrapper;
  private dialog?: DG.Dialog;
  private closeDialogSubject$ = new rxjs.Subject<void>();

  static async showDialog(): Promise<void> {
    if (!DialogWrapper._instance) {
      DialogWrapper._instance = new DialogWrapper();
      DialogWrapper._instance.closeDialogSubject$.subscribe(
        () => { DialogWrapper._instance.dialog = undefined; }
      );
    }

    if (!DialogWrapper._instance.dialog)
      DialogWrapper._instance.dialog = await DialogWrapper._instance.getDialog();

    DialogWrapper._instance.dialog.show();
  }

  private async getDialog(): Promise<DG.Dialog> {
    const eventManager = MonomerLibFileEventManager.getInstance();
    const widget = await MonomerLibraryManagerWidget.getContent(eventManager);
    const dialog = ui.dialog(
      {
        title: 'Manage monomer libraries',
        helpUrl: '/help/datagrok/solutions/domains/bio/bio.md#manage-monomer-libraries'
      }
    );
    $(dialog.root).css('width', '350px');
    dialog.clear();
    dialog.addButton(
      'Add',
      () => eventManager.addLibraryFile(),
      undefined,
      'Upload new HELM monomer library'
    );
    dialog.add(widget);
    dialog.onClose.subscribe(() => this.closeDialogSubject$.next());
    return dialog;
  }
}
