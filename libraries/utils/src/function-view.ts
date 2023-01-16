/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import wu from 'wu';
import {historyUtils} from './history-utils';
import {UiUtils} from './shared-components/ui-utils';

/**
   * Decorator to pass all thrown errors to grok.shell.error
   * @returns The actual funccall associated with the view
   * @stability Experimental
 */
export const passErrorToShell = () => {
  return (target: any, memberName: string, descriptor: PropertyDescriptor) => {
    const original = descriptor.value;

    descriptor.value = async function(...args: any[]) {
      try {
        return await original.call(this, ...args);
      } catch (err: any) {
        grok.shell.error((err as Error).message);
        throw err;
      }
    };
  };
};

export const INTERACTIVE_CSS_CLASS = 'cv-interactive';

export const defaultUsersIds = {
  'Test': 'ca1e672e-e3be-40e0-b79b-d2c68e68d380',
  'Admin': '878c42b0-9a50-11e6-c537-6bf8e9ab02ee',
  'System': '3e32c5fa-ac9c-4d39-8b4b-4db3e576b3c3',
};

export const defaultGroupsIds = {
  'All users': 'a4b45840-9a50-11e6-9cc9-8546b8bf62e6',
  'Developers': 'ba9cd191-9a50-11e6-9cc9-910bf827f0ab',
  'Need to create': '00000000-0000-0000-0000-000000000000',
  'Test': 'ca1e672e-e3be-40e0-b79b-8546b8bf62e6',
  'Admin': 'a4b45840-9a50-11e6-c537-6bf8e9ab02ee',
  'System': 'a4b45840-ac9c-4d39-8b4b-4db3e576b3c3',
  'Administrators': '1ab8b38d-9c4e-4b1e-81c3-ae2bde3e12c5',
};

export abstract class FunctionView extends DG.ViewBase {
  protected _funcCall?: DG.FuncCall;
  protected _lastCall?: DG.FuncCall;
  protected _type: string = 'function';

  constructor(funcCall?: DG.FuncCall) {
    super();
    this.box = true;

    if (!funcCall) return;
    this.name = funcCall.func.friendlyName;
  }

  /**
   * Get current function call of the view
   * @returns The actual funccall associated with the view
   * @stability Stable
 */
  public get funcCall(): DG.FuncCall | undefined {
    return this._funcCall;
  }

  /**
   * Get Func of the view
   * @returns The actual func associated with the view
   * @stability Stable
 */
  get func() {
    return this.funcCall?.func;
  }

  /**
   * Get data of last call of associated function
   * @returns The actual func associated with the view
   * @stability Stable
 */
  get lastCall() {
    return this._lastCall;
  }

  /**
   * Set data of last call of associated function
   * @stability Stable
 */
  set lastCall(lastCall: DG.FuncCall | undefined) {
    this._lastCall = lastCall;
  }

  /**
   * View type
   * @stability Stable
 */
  public get type(): string {
    return this._type;
  }

  /** Export options. Could be overriden partially, using default implementation of each option.
    * @stability Stable
  */
  exportConfig: {
    /** Override to provide custom export logic.
      *
      *  Default implementation {@link defaultExport} heavily relies on the default implementation of {@link buildIO}.
      * @returns Blob with data to be exported into the file.
      * @stability Stable
    */
    export: ((format: string) => Promise<Blob>);


    /** Filename for exported files. Override for custom filenames.
      * Default implementation is {@link defaultExportFilename}
      * @param format Format name to be exported
      * @returns The actual filename to be used for the generated file.
      * @stability Stable
    */
    filename: ((format: string) => string);

    /** Override to provide custom list of supported export formats.
     * Default implementation is {@link defaultSupportedExportFormats}
     * These formats are available under the "Export" popup on the ribbon panel.
     * @returns The array of formats available for the export.
     * @stability Stable
    */
    supportedFormats: string[];

    /** Override to provide custom file extensions for exported formats.
       * Default implementation is {@link defaultSupportedExportExtensions}
       * These extensions are used in filenames {@link exportFilename}.
       * @returns The mapping between supported export formats and their extensions.
       * @stability Stable
     */
    supportedExtensions: Record<string, string>;
  } | null = null;

  /**
   * Link FuncCall to the view
   * @param funcCall The actual funccall to be associated with the view
   * @stability Stable
 */
  public linkFunccall(funcCall: DG.FuncCall) {
    const isPreviousHistorical = this._funcCall?.options['isHistorical'];
    this._funcCall = funcCall;

    if (funcCall.options['isHistorical']) {
      if (!isPreviousHistorical)
        this.name = `${this.name} — ${funcCall.options['title'] ?? new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})}`;
      else
        this.name = `${this.name.substring(0, this.name.indexOf(' — '))} — ${funcCall.options['title'] ?? new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})}`;


      // FIX ME: view name does not change in models
      document.querySelector('div.d4-ribbon-name')?.replaceChildren(ui.span([this.name]));
      this.path = `?id=${this._funcCall.id}`;
    } else {
      this.path = ``;

      this.name = `${this.name.substring(0, (this.name.indexOf(' — ') > 0) ? this.name.indexOf(' — ') : undefined)}`;
    }
    this.buildRibbonPanels();
  }

  /**
   * Method for custom logic that could not be placed in the constructor.
   * Any async methods and most of the logic should be placed here.
   * @stability Stable
 */
  public async init() {}

  /**
   * Override to create a fully custom UI including ribbon menus and panels
   * @stability Stable
 */
  public build(): void {
    ui.empty(this.root);
    this.root.appendChild(this.buildIO());

    this.buildHistoryBlock();
    this.buildRibbonMenu();
  }

  /**
   * Override to create a custom input-output block
   * @returns The HTMLElement with whole UI excluding ribbon menus and panels
   * @stability Stable
 */
  public abstract buildIO(): HTMLElement;

  /**
   * Override to create a custom historical runs control.
   * @returns The HTMLElement with history block UI
   * @stability Stable
 */
  public buildHistoryBlock(): HTMLElement {
    const newHistoryBlock = UiUtils.historyPanel(this.func!);

    newHistoryBlock.onRunChosen.subscribe(async (id) => this.linkFunccall(await this.loadRun(id)));

    newHistoryBlock.beforeRunAddToFavorites.subscribe(async (funcCall) => {
      ui.setUpdateIndicator(newHistoryBlock.historyTab, true);
      ui.setUpdateIndicator(newHistoryBlock.favTab, true);

      funcCall = await this.addRunToFavorites(funcCall);

      newHistoryBlock.afterRunAddToFavorites.next(funcCall);

      ui.setUpdateIndicator(newHistoryBlock.historyTab, false);
      ui.setUpdateIndicator(newHistoryBlock.favTab, false);
    });
    newHistoryBlock.beforeRunAddToShared.subscribe(async (funcCall) => {
      ui.setUpdateIndicator(newHistoryBlock.historyTab, true);
      ui.setUpdateIndicator(newHistoryBlock.sharedTab, true);

      funcCall = await this.addRunToShared(funcCall);

      newHistoryBlock.afterRunAddToShared.next(funcCall);

      ui.setUpdateIndicator(newHistoryBlock.sharedTab, false);
      ui.setUpdateIndicator(newHistoryBlock.historyTab, false);
    });

    newHistoryBlock.beforeRunDeleted.subscribe(async (id) => {
      ui.setUpdateIndicator(newHistoryBlock.tabs.root, true);
      await this.deleteRun(await historyUtils.loadRun(id, true));

      newHistoryBlock.afterRunDeleted.next(id);
      ui.setUpdateIndicator(newHistoryBlock.tabs.root, false);
    });

    newHistoryBlock.beforeRunRemoveFromFavorites.subscribe(async (id) => {
      ui.setUpdateIndicator(newHistoryBlock.historyTab, true);
      ui.setUpdateIndicator(newHistoryBlock.favTab, true);
      await this.removeRunFromFavorites(await historyUtils.loadRun(id, true));

      newHistoryBlock.afterRunRemoveFromFavorites.next(id);

      ui.setUpdateIndicator(newHistoryBlock.favTab, false);
      ui.setUpdateIndicator(newHistoryBlock.historyTab, false);
    });

    newHistoryBlock.beforeRunRemoveFromShared.subscribe(async (id) => {
      ui.setUpdateIndicator(newHistoryBlock.historyTab, true);
      ui.setUpdateIndicator(newHistoryBlock.sharedTab, true);

      await this.removeRunFromShared(await historyUtils.loadRun(id, true));

      newHistoryBlock.afterRunRemoveFromShared.next(id);

      ui.setUpdateIndicator(newHistoryBlock.historyTab, false);
      ui.setUpdateIndicator(newHistoryBlock.sharedTab, false);
    });

    ui.empty(this.historyRoot);
    this.historyRoot.style.removeProperty('justify-content');
    this.historyRoot.style.width = '100%';
    this.historyRoot.append(newHistoryBlock.root);
    return newHistoryBlock.root;
  }

  /**
   * Looks for {@link supportedExportFormats} members and creates ribbon panel
   * @returns The HTMLElements of ribbonPanels
   * @stability Stable
 */
  buildRibbonPanels(): HTMLElement[][] {
    const newRibbonPanels: HTMLElement[][] = [
      [...(this.exportConfig && this.exportConfig.supportedFormats.length > 0) ? [ui.divH([
        ui.comboPopup(
          ui.iconFA('arrow-to-bottom'),
          this.exportConfig.supportedFormats,
          async (format: string) => DG.Utils.download(this.exportConfig!.filename(format), await this.exportConfig!.export(format))),
      ])]: []
      ]];

    if (this.func?.id) {
      const historyButton = ui.iconFA('history', () => {
        grok.shell.windows.showProperties = !grok.shell.windows.showProperties;
        historyButton.classList.toggle('d4-current');
        grok.shell.o = this.historyRoot;
      });

      historyButton.classList.add('d4-toggle-button');
      if (grok.shell.windows.showProperties) historyButton.classList.add('d4-current');

      newRibbonPanels.push([
        historyButton
      ]);
    }

    this.setRibbonPanels(newRibbonPanels);
    return newRibbonPanels;
  }

  /**
   * Override to create a custom ribbon menu on the top.
   * @stability Stable
 */
  public buildRibbonMenu() {

  }

  public async onBeforeRemoveRunFromFavorites(callToFavorite: DG.FuncCall) { }

  public async onAfterRemoveRunFromFavorites(favoriteCall: DG.FuncCall) { }

  /**
   * Saves the run as usual run
   * @param callToUnfavorite FuncCall object to remove from favorites
   * @returns Saved FuncCall
   * @stability Experimental
 */
  public async removeRunFromFavorites(callToUnfavorite: DG.FuncCall): Promise<DG.FuncCall> {
    callToUnfavorite.options['title'] = null;
    callToUnfavorite.options['annotation'] = null;
    callToUnfavorite.options['isFavorite'] = false;
    await this.onBeforeRemoveRunFromFavorites(callToUnfavorite);
    const favoriteSave = await grok.dapi.functions.calls.allPackageVersions().save(callToUnfavorite);
    await this.onAfterRemoveRunFromFavorites(favoriteSave);
    return favoriteSave;
  }

  public async onBeforeAddingToFavorites(callToAddToFavorites: DG.FuncCall) { }

  public async onAfterAddingToFavorites(favoriteCall: DG.FuncCall) { }

  /**
   * Saves the run as favorite
   * @param callToFavorite FuncCall object to add to favorites
   * @returns Saved FuncCall
   * @stability Experimental
 */

  public async addRunToFavorites(callToFavorite: DG.FuncCall): Promise<DG.FuncCall> {
    callToFavorite.options['isFavorite'] = true;
    await this.onBeforeAddingToFavorites(callToFavorite);
    const savedFavorite = await grok.dapi.functions.calls.allPackageVersions().save(callToFavorite);
    await this.onAfterAddingToFavorites(savedFavorite);
    return savedFavorite;
  }

  public async onBeforeRemoveRunFromShared(callToShare: DG.FuncCall) { }

  public async onAfterRemoveRunFromSahred(sharedCall: DG.FuncCall) { }

  /**
   * Removes run from shared
   * @param callToUnshare FuncCall object to remove from shared
   * @returns Saved FuncCall
   * @stability Experimental
 */

  public async removeRunFromShared(callToUnshare: DG.FuncCall): Promise<DG.FuncCall> {
    callToUnshare.options['title'] = null;
    callToUnshare.options['annotation'] = null;
    callToUnshare.options['isShared'] = false;
    await this.onBeforeRemoveRunFromFavorites(callToUnshare);
    const savedShared = await grok.dapi.functions.calls.allPackageVersions().save(callToUnshare);
    await this.onAfterRemoveRunFromFavorites(savedShared);
    return savedShared;
  }

  public async onBeforeAddingToShared(callToAddToShared: DG.FuncCall) { }

  public async onAfterAddingToShared(sharedCall: DG.FuncCall) { }

  /**
   * Saves the run as shared
   * @param callToShare FuncCall object to add to shared
   * @returns Saved FuncCall
   * @stability Experimental
 */

  public async addRunToShared(callToShare: DG.FuncCall): Promise<DG.FuncCall> {
    callToShare.options['isShared'] = true;
    await this.onBeforeAddingToShared(callToShare);

    const allGroup = await grok.dapi.groups.find(defaultGroupsIds['All users']);

    const dfOutputs = wu(callToShare.outputParams.values() as DG.FuncCallParam[])
      .filter((output) => output.property.propertyType === DG.TYPE.DATA_FRAME);

    for (const output of dfOutputs) {
      const df = callToShare.outputs[output.name] as DG.DataFrame;
      await grok.dapi.permissions.grant(df.getTableInfo(), allGroup, false);
    }

    const dfInputs = wu(callToShare.inputParams.values() as DG.FuncCallParam[])
      .filter((input) => input.property.propertyType === DG.TYPE.DATA_FRAME);
    for (const input of dfInputs) {
      const df = callToShare.inputs[input.name] as DG.DataFrame;
      await grok.dapi.permissions.grant(df.getTableInfo(), allGroup, false);
    }

    const savedShared = await grok.dapi.functions.calls.allPackageVersions().save(callToShare);
    await this.onAfterAddingToShared(savedShared);
    return savedShared;
  }

  /**
   * Called before saving the FUncCall results to the historical results, returns the saved call. See also {@link saveRun}.
   * @param callToSave FuncCall object to save
   * @returns Saved FuncCall
   * @stability Stable
 */
  public async onBeforeSaveRun(callToSave: DG.FuncCall) { }

  /**
   * Saves the computation results to the historical results, returns the saved call. See also {@link saveRun}.
   * @param savedCall FuncCall object to save
   * @returns Saved FuncCall
   * @stability Stable
 */
  public async onAfterSaveRun(savedCall: DG.FuncCall) { }

  /**
   * Saves the computation results to the historical results, returns the saved call. See also {@link loadRun}.
   * @param callToSave FuncCall object to save
   * @returns Saved FuncCall
   * @stability Stable
 */
  public async saveRun(callToSave: DG.FuncCall): Promise<DG.FuncCall> {
    await this.onBeforeSaveRun(callToSave);
    const savedCall = await historyUtils.saveRun(callToSave);
    savedCall.options['isHistorical'] = false;
    console.log(savedCall);
    this.linkFunccall(savedCall);
    this.buildHistoryBlock();
    this.path = `?id=${savedCall.id}`;
    await this.onAfterSaveRun(savedCall);
    return savedCall;
  }

  /**
   * Called before deleting the computation results from history, returns its id. See also {@link loadRun}.
   * @param callToDelete FuncCall object to be deleted
   * @stability Stable
 */
  public async onBeforeDeleteRun(callToDelete: DG.FuncCall) { }

  /**
   * Called after deleting the computation results from history, returns its id. See also {@link loadRun}.
   * @param deletedCall deleted FuncCall value
   * @stability Stable
 */
  public async onAfterDeleteRun(deletedCall: DG.FuncCall) { }

  /**
   * Deletes the computation results from history, returns its id. See also {@link loadRun}.
   * @param callToDelete FuncCall object to delete
   * @returns ID of deleted historical run
   * @stability Stable
 */

  public async deleteRun(callToDelete: DG.FuncCall): Promise<string> {
    await this.onBeforeDeleteRun(callToDelete);
    await historyUtils.deleteRun(callToDelete);
    await this.onAfterDeleteRun(callToDelete);
    return callToDelete.id;
  }

  /**
   * Called before fetching the historical run data in {@link loadRun}.
   * @stability Stable
 */
  public async onBeforeLoadRun() {}

  /**
   * Called after fetching the historical run data in {@link loadRun}.
   * @param funcCall FuncCall fetched from server during {@link loadRun}
   * @stability Stable
 */
  public async onAfterLoadRun(funcCall: DG.FuncCall) {}

  /**
   * Loads the specified historical run. See also {@link saveRun}.
   * @param funcCallId ID of FuncCall to look for. Get it using {@see funcCall.id} field
   * @returns FuncCall augemented with inputs' and outputs' values
   * @stability Stable
 */

  public async loadRun(funcCallId: string): Promise<DG.FuncCall> {
    await this.onBeforeLoadRun();
    const pulledRun = await historyUtils.loadRun(funcCallId);
    await this.onAfterLoadRun(pulledRun);
    return pulledRun;
  }

  /**
   * Called before actual computations are made {@link run}.
   * @param funcToCall FuncCall object to be called {@see DG.FuncCall.call()}
   * @stability Experimental
  */
  public async onBeforeRun(funcToCall: DG.FuncCall) {}

  /**
    * Called after actual computations are made {@link run}.
    * @param runFunc FuncCall object after call method {@see DG.FuncCall.call()}
    * @stability Experimental
   */
  public async onAfterRun(runFunc: DG.FuncCall) {}

  public async run(): Promise<void> {
    if (!this.funcCall) throw new Error('The correspoding function is not specified');

    await this.onBeforeRun(this.funcCall);
    const pi = DG.TaskBarProgressIndicator.create('Calculating...');
    this.funcCall.newId();
    await this.funcCall.call(); // mutates the funcCall field
    pi.close();
    await this.onAfterRun(this.funcCall);

    this.lastCall = await this.saveRun(this.funcCall);
  }

  protected historyRoot: HTMLDivElement = ui.divV([], {style: {'justify-content': 'center'}});

  protected defaultExportFilename = (format: string) => {
    return `${this.name} - ${new Date().toLocaleString()}.${this.exportConfig!.supportedExtensions[format]}`;
  };

  protected defaultSupportedExportExtensions = () => {
    return {
      'Excel': 'xlsx'
    };
  };

  protected defaultSupportedExportFormats = () => {
    return ['Excel'];
  };
}
