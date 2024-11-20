/* eslint-disable max-lines */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {Subject} from 'rxjs';
import './style.css';

import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {getMonomerLibHelper, IMonomerLibFileManager} from '@datagrok-libraries/bio/src/monomer-works/monomer-utils';

import {MonomerLibFileEventManager} from './event-manager';
import {_package} from '../../../package';
import {MonomerManager} from '../monomer-manager/monomer-manager';
import {DuplicateMonomerManager} from '../monomer-manager/duplicate-monomer-manager';
import {MonomerLibManager} from '../lib-manager';

export async function showManageLibrariesDialog(): Promise<void> {
  await DialogWrapper.showDialog();
}

export async function showManageLibrariesView() {
  await LibManagerView.showView();
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
  private _fileManager: IMonomerLibFileManager;

  private _widget: DG.Widget;
  public get widget(): DG.Widget { return this._widget; }

  private constructor() {}

  private static instancePromise?: Promise<MonomerLibraryManagerWidget>;

  static async getInstance(): Promise<MonomerLibraryManagerWidget> {
    if (MonomerLibraryManagerWidget.instancePromise === undefined) {
      MonomerLibraryManagerWidget.instancePromise = (async () => {
        const instance = new MonomerLibraryManagerWidget();
        const libHelper = await getMonomerLibHelper();
        instance._fileManager = await libHelper.getFileManager();
        instance._widget = await instance.createWidget();
        return instance;
      })();
    }
    return MonomerLibraryManagerWidget.instancePromise;
  }

  static async reloadWidget(): Promise<void> {
    const instance = await MonomerLibraryManagerWidget.getInstance();
    instance._widget = await instance.createWidget();
  }

  private async createWidget() {
    const content = await this.getWidgetContent();
    const monomerLibHelper = await getMonomerLibHelper();
    monomerLibHelper.eventManager.addLibraryFileRequested$.subscribe(
      async () => await this.promptToAddLibraryFiles()
    );
    return new DG.Widget(content);
  }

  private async getWidgetContent(): Promise<HTMLElement> {
    const libControlsForm = await LibraryControlsManager.createControlsForm();
    $(libControlsForm).addClass('monomer-lib-controls-form');
    setTimeout(() => {
      libControlsForm && $(libControlsForm) && $(libControlsForm).removeClass('ui-form-condensed');
    }, 200);
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
          await this._fileManager.addLibraryFile(content, name);
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
  private constructor(
    private fileManager: IMonomerLibFileManager,
    private readonly userLibSettings: UserLibSettings,
  ) {
    this.fileManager.eventManager.updateUIControlsRequested$.subscribe(() => {
      this.updateControlsForm();
    });
    this.fileManager.eventManager.librarySelectionRequested$.subscribe(async ([fileName, isSelected]) => {
      await this.updateLibrarySelectionStatus(isSelected, fileName);
    });
  }

  private toLog(): string {
    return `LibraryControlsManager<#>`;
  }

  static async createControlsForm(): Promise<HTMLElement> {
    const logPrefix = 'LibraryControlsForm.createControlsForm()';
    _package.logger.debug(`${logPrefix}, start`);
    const [fileManager, userLibSettings] = await Promise.all([
      getMonomerLibHelper().then((libHelper) => libHelper.getFileManager()),
      await getUserLibSettings(),
    ]);
    const manager = new LibraryControlsManager(fileManager, userLibSettings);

    return manager._createControlsForm();
  }

  private _createControlsForm(): HTMLElement {
    const libraryControls = this.createLibraryControls();
    const inputsForm = ui.wideForm(libraryControls, undefined);
    $(inputsForm).addClass('monomer-lib-controls-form');

    return inputsForm;
  }

  public updateControlsForm(): void {
    const updatedForm = this._createControlsForm();
    $('.monomer-lib-controls-form').replaceWith(updatedForm);
  }

  private createLibraryControls(): DG.InputBase<boolean | null>[] {
    const libFileNameList: string[] = this.fileManager.getValidLibraryPaths();
    return libFileNameList.map((libFileName) => this.createLibInput(libFileName));
  }

  private createLibInput(libFileName: string): DG.InputBase<boolean | null> {
    const logPrefix = `${this.toLog()}.createLibInput()`;
    _package.logger.debug(`${logPrefix}, libFileName = '${libFileName}', start`);
    const isMonomerLibrarySelected = !this.userLibSettings.exclude.includes(libFileName);
    const libInput = ui.input.bool(libFileName, {value: isMonomerLibrarySelected, onValueChanged: () => {
      this.fileManager.eventManager.updateLibrarySelectionStatus(libFileName, libInput.value);
    }});
    ui.tooltip.bind(libInput.root, `Include monomers from ${libFileName}`);
    const deleteIcon = ui.iconFA('trash-alt', () => this.promptForLibraryDeletion(libFileName));
    const editIcon = ui.icons.edit(async () => {
      grok.shell.v = await (await MonomerManager.getInstance()).getViewRoot(libFileName);
    }, 'Edit monomer library');
    ui.tooltip.bind(deleteIcon, `Delete ${libFileName}`);
    libInput.addOptions(editIcon);
    libInput.addOptions(deleteIcon);
    _package.logger.debug(`${logPrefix}, libFileName = '${libFileName}', end`);
    return libInput;
  }

  private async updateLibrarySelectionStatus(
    isMonomerLibrarySelected: boolean,
    libFileName: string
  ): Promise<void> {
    this.updateLibrarySettings(isMonomerLibrarySelected, libFileName);
    await setUserLibSettings(this.userLibSettings);
    const monomerLibHelper = await getMonomerLibHelper();
    await monomerLibHelper.loadMonomerLib(true);
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
          await this.updateLibrarySelectionStatus(false, fileName);
          await this.fileManager.deleteLibraryFile(fileName);
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
  private closeDialogSubject$ = new Subject<void>();

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
    const widget = (await MonomerLibraryManagerWidget.getInstance()).widget;
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

class LibManagerView {
  private constructor() {};
  private static _instance: LibManagerView;
  private _view: DG.View;
  private _duplicateManager: DuplicateMonomerManager;
  private libManager: MonomerLibManager;
  private async getView() {
    const eventManager = MonomerLibFileEventManager.getInstance();
    const widget = (await MonomerLibraryManagerWidget.getInstance()).widget;
    const addButton = ui.bigButton('Add',
      () => eventManager.addLibraryFile(), 'Upload new HELM monomer library');
    const mergeButton =
      ui.bigButton('Merge', () => { this.mergeSelectedLibs(); }, 'Merge selected libraries into one');

    const v = ui.splitH(
      [ui.divV([widget.root, ui.buttonsInput([addButton, mergeButton])], {classes: 'ui-form'}),
        this._duplicateManager.root],
      {style: {width: '100%', height: '100%'}},
      true);
    this._view = grok.shell.newView('Manage Monomer Libraries', [v]);

    ui.tools.waitForElementInDom(v).then(() => {
      setTimeout(() => {
        const children = Array.from(v.children as HTMLCollectionOf<HTMLElement>)
          .filter((el) => el.classList.contains('ui-box'));
        if (children.length !== 2)
          return;
        const [left, right] = children;
        const combinedWidth = left.getBoundingClientRect().width + right.getBoundingClientRect().width;
        const leftWidth = combinedWidth * 0.3;
        left.style.width = `${leftWidth}px`;
        const rightWidth = combinedWidth - leftWidth;
        right.style.width = `${rightWidth}px`;
      }, 100);
      this._view.subs.push(grok.events.onCurrentViewChanged.subscribe(() => {
        try {
          const inst = LibManagerView._instance;
          if (inst && inst._view && 'id' in grok.shell.v && grok.shell.v.id === inst._view.id)
            inst._duplicateManager?.refresh();
        } catch (e) {
          console.error(e);
        }
      }));
    });
    //grok.shell.dockManager.dock(this._duplicateManager.root, DG.DOCK_TYPE.RIGHT, null, '', 0.4);
  }

  static async showView() {
    if (!LibManagerView._instance)
      LibManagerView._instance = new LibManagerView();
    if (!LibManagerView._instance._duplicateManager)
      LibManagerView._instance._duplicateManager = await DuplicateMonomerManager.getInstance();
    if (!LibManagerView._instance.libManager)
      LibManagerView._instance.libManager = await MonomerLibManager.getInstance();
    if (LibManagerView._instance._view &&
        Array.from(grok.shell.views).find((v) => v.id && v.id === LibManagerView._instance._view.id)) {
      grok.shell.v = LibManagerView._instance._view;
      await LibManagerView._instance._duplicateManager.refresh();
      return;
    }
    LibManagerView._instance.getView();
  }
  async mergeSelectedLibs() {
    const libraryExistsError = 'Library with this name already exists';
    const libManager = await MonomerLibManager.getInstance();
    await libManager.awaitLoaded();
    await libManager.loadLibrariesPromise;
    if (!libManager.duplicatesHandled) {
      grok.shell.warning(`Selected libraries contain repeating symbols with different monomers.
        Please choose the correct monomer for each symbol using duplicate monomomer manager.`);
      return;
    }
    const libJSON = libManager.getBioLib().toJSON();
    const dialog = ui.dialog('Merge selected libraries');
    const newFileNameInput = ui.input.string('Library Name', {
      placeholder: 'Enter new library name',
      nullable: false,
      onValueChanged: () => {
        const res = validateInput(newFileNameInput.value);
        dialog.getButton('Download')?.classList?.toggle('d4-disabled', !!res && res !== libraryExistsError);
        dialog.getButton('Save')?.classList?.toggle('d4-disabled', !!res);
      }
    });
    const validLibPaths = (await this.libManager.getFileManager()).getValidLibraryPaths();
    newFileNameInput.addValidator(validateInput);
    function getFileNameInputValue() {
      let fileName = newFileNameInput.value;
      if (!fileName.endsWith('.json'))
        fileName += '.json';
      return fileName;
    };

    function validateInput(v: string) {
      if (!v || !v.trim()) return 'Library name cannot be empty';
      if ((v.endsWith('.json') && validLibPaths.includes(v)) || validLibPaths.includes(v + '.json'))
        return libraryExistsError;
      return null;
    }
    dialog
      .add(newFileNameInput)
      .add(ui.divText(`Total monomers: ${libJSON.length}`))
      .addButton('Download', () => { DG.Utils.download(getFileNameInputValue(), JSON.stringify(libJSON)); })
      .addButton('Save', async () => {
        dialog.close();
        const fileName = getFileNameInputValue();
        const content = JSON.stringify(libJSON);
        const fileManager = await this.libManager.getFileManager();
        this._view && ui.setUpdateIndicator(this._view.root, true);
        try {
          await fileManager.addLibraryFile(content, fileName, false); // we will reload after updating settings
          const settings = await getUserLibSettings();
          settings.exclude = validLibPaths; // exclude all previous libraries
          await setUserLibSettings(settings);
          await this.libManager.loadLibraries(true);
          await MonomerLibraryManagerWidget.reloadWidget();
        } catch (e) {
          grok.shell.error(`Failed to save library ${fileName}. see console for details.`);
          console.error(e);
        } finally {
          this._view && ui.setUpdateIndicator(this._view.root, false);
        }
      })
      .show();
    dialog.getButton('Download')?.classList?.add('d4-disabled');
    dialog.getButton('Save')?.classList?.add('d4-disabled');
  }
}


