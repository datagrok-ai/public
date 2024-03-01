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
    this.eventManager.saveLibrarySettingsRequested$.subscribe(
      async () => await this.saveLibrarySettings()
    );

    this.eventManager.resetLibrarySettingsRequested$.subscribe(
      () => this.resetLibrarySettings()
    );
  }
  private monomerLibFileManager: MonomerLibFileManager;
  private originalLibSettings: UserLibSettings;
  private updatedLibSettings: UserLibSettings;

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
    this.originalLibSettings = await getUserLibSettings();
    this.updatedLibSettings = {...this.originalLibSettings};
  };

  private async updateControlsForm(): Promise<void> {
    const updatedForm = await this._createControlsForm();
    $('.monomer-lib-controls-form').replaceWith(updatedForm);
  }

  private async getFileNameList(): Promise<string[]> {
    const fileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    return fileManager.getValidLibraryPaths();
  }

  private async createLibraryControls(): Promise<DG.InputBase<boolean | null>[]> {
    const libFileNameList = await this.getFileNameList();
    return libFileNameList.map((libFileName) => this.createLibInput(libFileName));
  }

  private createLibInput(libFileName: string): DG.InputBase<boolean | null> {
    const isMonomerLibrarySelected = !this.updatedLibSettings.exclude.includes(libFileName);
    const libInput = ui.switchInput(
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
    const libFileNameList = await this.getFileNameList();
    await MonomerLibManager.instance.applyUserLibSettings(libFileNameList, this.updatedLibSettings, true);
  }

  private async saveLibrarySettings(): Promise<void> {
    console.log('updatedLibSettings', this.updatedLibSettings);
    console.log('originalLibSettings', this.originalLibSettings);
    if (JSON.stringify(this.updatedLibSettings) === JSON.stringify(this.originalLibSettings))
      return;

    await setUserLibSettings(this.updatedLibSettings);
    await MonomerLibManager.instance.loadLibraries(true);
    this.originalLibSettings = {...this.updatedLibSettings};
    grok.shell.info('Monomer library settings saved');
  }

  private resetLibrarySettings():void {
    this.updatedLibSettings = {...this.originalLibSettings};
  }

  private updateLibrarySettings(
    isLibrarySelected: boolean | null,
    libFileName: string,
  ): void {
    if (isLibrarySelected) {
      // Remove selected library from exclusion list
      this.updatedLibSettings.exclude = this.updatedLibSettings.exclude.filter((libName) => libName !== libFileName);
    } else if (!this.updatedLibSettings.exclude.includes(libFileName)) {
      // Add unselected library to exclusion list
      this.updatedLibSettings.exclude.push(libFileName);
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
  private eventManager = MonomerLibFileEventManager.getInstance();

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

  private async getWidgetContent(): Promise<DG.Widget> {
    const widget = await MonomerLibraryManagerWidget.getContent(this.eventManager);
    return widget;
  }

  private getAddFileControl(): HTMLElement {
    const addFileControl = ui.link(
      'Add file',
      () => this.eventManager.addLibraryFile(),
      'Upload new HELM monomer library',
    );
    $(addFileControl).addClass('d4-link-action');
    $(addFileControl).removeClass(' d4-link-external');
    $(addFileControl).css('padding', '10px 0 0 20px');

    return addFileControl;
  }


  private async getDialog(): Promise<DG.Dialog> {
    const dialog = ui.dialog(
      {
        title: 'Manage monomer libraries',
        helpUrl: '/help/datagrok/solutions/domains/bio/bio.md#manage-monomer-libraries'
      }
    );
    $(dialog.root).css('width', '350px');
    dialog.clear();

    const widget = await this.getWidgetContent();
    dialog.add(widget);

    const addFileControl = this.getAddFileControl();
    dialog.add(addFileControl);

    dialog.addButton(
      'APPLY',
      () => this.eventManager.saveLibrarySettings(),
      undefined,
      'Apply selection'
    );

    dialog.onOK(() => this.eventManager.saveLibrarySettings());

    const closeDialogHandler = () => {
      this.closeDialogSubject$.next();
      this.eventManager.resetLibrarySettings();
    };
    dialog.onClose.subscribe(() => closeDialogHandler());
    dialog.onCancel(() => closeDialogHandler());
    return dialog;
  }
}
