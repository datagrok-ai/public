/* eslint-disable rxjs/no-async-subscribe */
/* eslint-disable rxjs/no-ignored-subscription */
/* eslint-disable max-lines */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import $ from 'cash-dom';
import {Subject, Subscription} from 'rxjs';

import {
  getUserLibSettings, setUserLibSettings
} from '@datagrok-libraries/bio/src/monomer-works/lib-settings';
import {UserLibSettings} from '@datagrok-libraries/bio/src/monomer-works/types';
import {findProviderWithLibraryName,
  getMonomerLibHelper, IMonomerLibHelper,
  IMonomerLibProvider} from '@datagrok-libraries/bio/src/types/monomer-library';

import {_package} from '../../../package';
import {MonomerManager} from '../monomer-manager/monomer-manager';
import {DuplicateMonomerManager} from '../monomer-manager/duplicate-monomer-manager';
import {MonomerLibManager} from '../lib-manager';
// @ts-ignore
import './style.css';

export async function showManageLibrariesDialog(): Promise<void> {
  await DialogWrapper.showDialog();
}

export async function showManageLibrariesView(addView = true) {
  return await LibManagerView.showView(addView);
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
  private _widget: DG.Widget;
  public get widget(): DG.Widget { return this._widget; }

  private constructor() {}

  private static instancePromise?: Promise<MonomerLibraryManagerWidget>;

  private libHelper: IMonomerLibHelper;

  static async getInstance(): Promise<MonomerLibraryManagerWidget> {
    if (MonomerLibraryManagerWidget.instancePromise == undefined) {
      MonomerLibraryManagerWidget.instancePromise = (async () => {
        const instance = new MonomerLibraryManagerWidget();
        const libHelper = await getMonomerLibHelper();
        instance.libHelper = libHelper;
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

  private _fileUploadSubscription: Subscription | null = null;
  private async createWidget() {
    const content = await this.getWidgetContent();
    const monomerLibHelper = await getMonomerLibHelper();
    this._fileUploadSubscription?.unsubscribe();
    this._fileUploadSubscription =
    monomerLibHelper.fileUploadRequested.subscribe(
      () => this.promptToAddLibraryFiles()
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
        const doAdd = async (provider: IMonomerLibProvider) => {
          const name = selectedFile.name;
          const existingLibs = await provider.listLibraries();
          // chech if library already exists
          if (existingLibs.includes(name)) {
            const confirm = await new Promise<boolean>((resolve) => {
              ui.dialog('Confirm Library Update')
                .add(ui.divText(`Library '${name}' already exists. Do you want to overwrite it?`))
                .onOK(() => resolve(true))
                .onCancel(() => resolve(false))
                .show();
            });
            if (!confirm)
              return;
          }

          const content = await selectedFile.text();
          const progressIndicator = DG.TaskBarProgressIndicator.create(`Adding ${name} as a monomer library`);
          try {
            await provider.addOrUpdateLibraryString(name, content);
          // this.eventManager.updateLibrarySelectionStatus(name, true);
          } catch (e) {
            grok.shell.error(`File ${name} is not a valid monomer library, verify it is aligned to HELM JSON schema.`);
            console.error(e);
          } finally {
            progressIndicator.close();
          }
        };
        const providers = await this.libHelper.getProviders();
        if (providers.length === 0) {
          grok.shell.error('No monomer library providers available to add the library.');
          return;
        }
        if (providers.length === 1) {
          await doAdd(providers[0]);
          return;
        }
        const dialog = ui.dialog('Select storage for new monomer library');
        const providersInput =
        ui.input.choice('Storage', {items: providers.map((p) => p.name), value: providers[0].name,
          nullable: false, tooltipText: 'Storage provider for new monomer library'});
        dialog
          .add(providersInput)
          .onOK(async () => {
            const provider = providers.find((p) => p.name === providersInput.value)!; // should not be null
            await doAdd(provider);
          }).show();
      },
    });
  }
}

class LibraryControlsManager {
  private constructor(
    private readonly libHelper: IMonomerLibHelper,
    private readonly userLibSettings: UserLibSettings,
  ) {
    this.libHelper.providersDataChanged.subscribe(async () => {
      await this.updateControlsForm();
    });
  }

  private toLog(): string {
    return `LibraryControlsManager<#>`;
  }

  static async createControlsForm(): Promise<HTMLElement> {
    const logPrefix = 'LibraryControlsForm.createControlsForm()';
    _package.logger.debug(`${logPrefix}, start`);
    const userLibSettings = await getUserLibSettings();
    const libHelper = await getMonomerLibHelper();
    const manager = new LibraryControlsManager(libHelper, userLibSettings);

    return manager._createControlsForm();
  }

  private async _createControlsForm(): Promise<HTMLElement> {
    const libraryControls = await this.createLibraryControls();
    const inputsForm = ui.wideForm(libraryControls, undefined);
    $(inputsForm).addClass('monomer-lib-controls-form');

    return inputsForm;
  }

  public async updateControlsForm(): Promise<void> {
    const updatedForm = await this._createControlsForm();
    $('.monomer-lib-controls-form').replaceWith(updatedForm);
  }

  private async createLibraryControls(): Promise<DG.InputBase<boolean | null>[]> {
    const libFileNameList: string[] = await this.libHelper.getAvaliableLibraryNames();
    return libFileNameList.map((libFileName) => this.createLibInput(libFileName));
  }

  private createLibInput(libFileName: string): DG.InputBase<boolean | null> {
    const logPrefix = `${this.toLog()}.createLibInput()`;
    _package.logger.debug(`${logPrefix}, libFileName = '${libFileName}', start`);
    const isMonomerLibrarySelected = !this.userLibSettings.exclude.includes(libFileName);
    const libInput = ui.input.bool(libFileName, {value: isMonomerLibrarySelected, onValueChanged: () => {
      updateLibrarySelectionStatus(libInput.value, libFileName);
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


  private promptForLibraryDeletion(fileName: string): void {
    const dialog = ui.dialog('Warning');
    dialog.add(ui.divText(`Are you sure you want to delete library "${fileName}"?`))
      .onOK(async () => {
        try {
          const progressIndicator = DG.TaskBarProgressIndicator.create(`Deleting ${fileName} library`);
          await updateLibrarySelectionStatus(false, fileName);
          const provider = await findProviderWithLibraryName(await this.libHelper.getProviders(), fileName);
          if (!provider)
            throw new Error(`Cannot find provider for library ${fileName}`);
          await provider.deleteLibrary(fileName);
          // await this.fileManager.deleteLibraryFile(fileName);
          progressIndicator.close();
        } catch (e) {
          console.error(e);
          grok.shell.error(`Failed to delete ${fileName} library`);
        }
      })
      .showModal(false);
  }
}

async function updateLibrarySelectionStatus(
  isMonomerLibrarySelected: boolean,
  libFileName: string
): Promise<void> {
  const userLibSettings = await getUserLibSettings();
  updateLibrarySettings(userLibSettings, isMonomerLibrarySelected, libFileName);
  await setUserLibSettings(userLibSettings);
  const monomerLibHelper = await getMonomerLibHelper();
  await monomerLibHelper.loadMonomerLib(true);
  grok.shell.info('Monomer library user settings saved');
  monomerLibHelper.notifyLibrarySelectionChanged();
}

function updateLibrarySettings(
  userLibSettings: UserLibSettings,
  isLibrarySelected: boolean | null,
  libFileName: string,
): void {
  if (isLibrarySelected) {
    // Remove selected library from exclusion list
    userLibSettings.exclude = userLibSettings.exclude.filter((libName) => libName !== libFileName);
  } else if (!userLibSettings.exclude.includes(libFileName)) {
    // Add unselected library to exclusion list
    userLibSettings.exclude.push(libFileName);
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
      // eslint-disable-next-line rxjs/no-ignored-subscription
      DialogWrapper._instance.closeDialogSubject$.subscribe(
        () => { DialogWrapper._instance.dialog = undefined; }
      );
    }

    if (!DialogWrapper._instance.dialog)
      DialogWrapper._instance.dialog = await DialogWrapper._instance.getDialog();

    DialogWrapper._instance.dialog.show();
  }

  private async getDialog(): Promise<DG.Dialog> {
    const widget = (await MonomerLibraryManagerWidget.getInstance()).widget;
    const libHelper = await getMonomerLibHelper();
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
      () => libHelper.requestFileUpload(),
      undefined,
      'Upload new HELM monomer library'
    );
    dialog.add(widget);
    const sub = dialog.onClose.subscribe(() => {
      this.closeDialogSubject$.next();
      sub.unsubscribe();
    });
    return dialog;
  }
}

class LibManagerView {
  private constructor(
  ) {};
  private static _instance: LibManagerView;
  static viewName = 'Manage Monomer Libraries';
  private _view: DG.View;
  private _duplicateManager: DuplicateMonomerManager;
  private libManager: MonomerLibManager;
  private async getView(addView = true) {
    const widget = (await MonomerLibraryManagerWidget.getInstance()).widget;
    const addButton = ui.bigButton('Add',
      () => this.libManager.requestFileUpload(), 'Upload new HELM monomer library');
    const mergeButton =
      ui.bigButton('Merge', () => { this.mergeSelectedLibs(); }, 'Merge selected libraries into one');

    const v = ui.splitH(
      [ui.divV([widget.root, ui.buttonsInput([addButton, mergeButton])], {classes: 'ui-form'}),
        this._duplicateManager.root],
      {style: {width: '100%', height: '100%'}},
      true);
    if (this._view) {
      try {
        this._view.subs.forEach((s) => s.unsubscribe());
        this._view.detach();
        this._view.close();
      } catch (_e) {
      }
    }
    this._view = DG.View.fromRoot(v);
    this._view.name = LibManagerView.viewName;
    if (addView)
      grok.shell.addView(this._view);

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
          if (inst && inst._view && grok.shell.v && 'id' in grok.shell.v && grok.shell.v.id === inst._view.id)
            inst._duplicateManager?.refresh();
        } catch (e) {
          console.error(e);
        }
      }));
    });
    return this._view;
    //grok.shell.dockManager.dock(this._duplicateManager.root, DG.DOCK_TYPE.RIGHT, null, '', 0.4);
  }

  static async showView(addView = true) {
    if (!LibManagerView._instance)
      LibManagerView._instance = new LibManagerView();
    if (!LibManagerView._instance._duplicateManager)
      LibManagerView._instance._duplicateManager = await DuplicateMonomerManager.getInstance();
    if (!LibManagerView._instance.libManager)
      LibManagerView._instance.libManager = await MonomerLibManager.getInstance();
    if (addView && LibManagerView._instance._view &&
        Array.from(grok.shell.views).find((v) => v.id && v.id === LibManagerView._instance._view.id)) {
      grok.shell.v = LibManagerView._instance._view;
      await LibManagerView._instance._duplicateManager.refresh();
      return LibManagerView._instance._view;
    }
    // something can conflict with browse view, so need to make sure that we close all existing views
    LibManagerView.closeExistingViews();
    return LibManagerView._instance.getView(addView);
  }

  private static closeExistingViews() {
    Array.from(grok.shell.views).filter((v) => v.name === LibManagerView.viewName).forEach((v) => v.close());
  }

  async mergeSelectedLibs() {
    const libraryExistsError = 'Library with this name already exists';
    const libManager = await MonomerLibManager.getInstance();
    await libManager.awaitLoaded();
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
    const validLibPaths = await this.libManager.getAvaliableLibraryNames();
    newFileNameInput.addValidator(validateInput);

    function validateInput(v: string) {
      if (!v || !v.trim()) return 'Library name cannot be empty';
      if (validLibPaths.includes(v) || validLibPaths.includes(v + '.json'))
        return libraryExistsError;
      return null;
    }
    const providers = await this.libManager.getProviders();
    if (providers.length === 0) {
      grok.shell.error('No monomer library providers available to save the merged library.');
      return; // I mean, this should not happen, but ....
    }
    const providersInput = ui.input.choice('Storage', {items: providers.map((p) => p.name), value: providers[0].name,
      nullable: false, tooltipText: 'Storage provider for saving new monomer library'});

    const getFileNameInputValue = () => {
      let fileName = newFileNameInput.value.trim();
      if (!fileName)
        fileName = 'New';
      if (!fileName.toLowerCase().endsWith('.json'))
        fileName += '.json';
      return fileName;
    };
    dialog
      .add(providersInput)
      .add(newFileNameInput)
      .add(ui.divText(`Total monomers: ${libJSON.length}`))
      .addButton('Download', () => { DG.Utils.download(getFileNameInputValue(), JSON.stringify(libJSON)); })
      .addButton('Save', async () => {
        if (!newFileNameInput.value || !newFileNameInput.value.trim() || !providersInput.value) {
          providersInput.validate();
          newFileNameInput.validate(); // this will force showing validation error
          return;
        }
        dialog.close();
        const fileName = newFileNameInput.value!;
        const content = JSON.stringify(libJSON);
        this._view && ui.setUpdateIndicator(this._view.root, true);
        try {
          const provider = providers.find((p) => p.name === providersInput.value)!; // should not be null
          await provider.addOrUpdateLibraryString(fileName, content); // we will reload after updating settings
          const settings = await getUserLibSettings();
          settings.exclude = validLibPaths; // exclude all previous libraries
          await setUserLibSettings(settings);
          await this.libManager.loadMonomerLib(true);
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


