/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import * as rxjs from 'rxjs';
import './style.css';

import {
  getUserLibSettings, LIB_PATH, setUserLibSettings
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
    const formHandler = new LibraryControlsManager(this.eventManager);
    const libControlsForm = await formHandler.createControlsForm();
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

  async createControlsForm(): Promise<HTMLElement> {
    this.monomerLibFileManager = await MonomerLibFileManager.getInstance(this.eventManager);
    const libraryControls = await this.createLibraryControls();
    const inputsForm = ui.form(libraryControls);
    $(inputsForm).addClass('monomer-lib-controls-form');

    return inputsForm;
  }

  private async updateControlsForm(): Promise<void> {
    const updatedForm = await this.createControlsForm();
    $('.monomer-lib-controls-form').replaceWith(updatedForm);
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
    dialog.add(ui.divText(`Delete file ${fileName}?`))
      .onOK(async () => {
        const progressIndicator = DG.TaskBarProgressIndicator.create(`Deleting ${fileName} library`);
        await this.monomerLibFileManager.deleteLibraryFile(fileName);
        progressIndicator.close();
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
    const dialog = ui.dialog('Manage monomer libraries');
    $(dialog.root).css('width', '350px');
    dialog.clear();
    dialog.addButton(
      'Add',
      () => eventManager.addLibraryFile(),
      undefined,
      'Load new monomer libraries'
    );
    dialog.add(widget);
    dialog.onClose.subscribe(() => this.closeDialogSubject$.next());
    return dialog;
  }
}
