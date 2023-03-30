/* eslint-disable valid-jsdoc */
/* eslint-disable max-len */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {Subject, BehaviorSubject} from 'rxjs';
import {historyUtils} from './history-utils';
import {UiUtils} from './shared-components/ui-utils';

// Getting inital URL user entered with
const startUrl = new URL(grok.shell.startUri);

export abstract class FunctionView extends DG.ViewBase {
  protected _funcCall?: DG.FuncCall;
  protected _lastCall?: DG.FuncCall;
  protected _type: string = 'function';

  // emitted when after a new FuncCall is linked
  protected funcCallReplaced = new Subject<true>();

  // emitted when after an initial FuncCall is linked
  public onFuncCallReady = new BehaviorSubject<false>(false);

  /**
   * Constructs a new view using function with the given {@link funcName}. An fully-specified name is expected.
   * Search of the function is async, so async {@link init} function is used.
   * All other functions are called only when initialization is over and {@link this.onFuncCallReady} is emitted.
   * @param funcName Name of DG.Func (either script or package function) to use as view foundation
   * @param options Configuration object for the view.
   */
  constructor(
    protected funcName: string,
    public options: {historyEnabled: boolean, isTabbed: boolean} = {historyEnabled: true, isTabbed: false}
  ) {
    super();
    this.box = true;

    // Changing view and building IO are reasonable only after FuncCall is linked
    this.onFuncCallReady.subscribe({
      complete: async () => {
        this.changeViewName(this.funcCall.func.friendlyName);
        this.build();

        if (this.getStartId()) {
          await this.onBeforeLoadRun();
          this.lastCall = this.funcCall;
          await this.onAfterLoadRun(this.funcCall);

          this.setAsLoaded();
        }
      }
    });

    this.init();
  }

  /**
   * Changes the name of the view. This method also deals with rare bug when view name is not updated after change.
   * @param newName New name for the view
   */
  protected changeViewName(newName: string) {
    // TODO: Find a reproducible sample of the bug
    this.name = newName;
    document.querySelector('div.d4-ribbon-name')?.replaceChildren(ui.span([newName]));
  }

  private getStartId(): string | undefined {
    // To prevent loading same ID on opening different package,
    // we should check if we have already loaded run by this ID

    //@ts-ignore
    return (!grok.shell.getVar('isLoaded')) ? startUrl.searchParams.get('id'): undefined;
  }

  private setAsLoaded(): string | undefined {
    // @ts-ignore
    return grok.shell.setVar('isLoaded', true);
  }

  /**
   * Get current function call of the view
   * @returns The actual funccall associated with the view
   * @stability Stable
 */
  public get funcCall(): DG.FuncCall {
    return this._funcCall!;
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
      * There is no default implementation, since, in general, export is dependent on the UI.
      *
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
   * Links FuncCall to the view. In addition, sets "path" and "name" properties to corresponding ones.
   * After linking, emits {@link this.funcCallReplaced} event.
   * @param funcCall The actual funccall to be associated with the view
   * @stability Stable
 */
  public linkFunccall(funcCall: DG.FuncCall) {
    const isPreviousHistorical = this._funcCall?.options['isHistorical'];
    this._funcCall = funcCall;

    if (!this.options.isTabbed) {
      if (funcCall.options['isHistorical']) {
        if (!isPreviousHistorical)
          this.changeViewName(`${this.name} — ${funcCall.options['title'] ?? new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})}`);
        else
          this.changeViewName(`${this.name.substring(0, (this.name.indexOf(' — ') > 0) ? this.name.indexOf(' — ') : undefined)} — ${funcCall.options['title'] ?? new Date(funcCall.started.toString()).toLocaleString('en-us', {month: 'short', day: 'numeric', hour: 'numeric', minute: 'numeric'})}`);
        this.path = `?id=${this._funcCall!.id}`;
      } else {
        this.path = ``;

        this.changeViewName(`${this.name.substring(0, (this.name.indexOf(' — ') > 0) ? this.name.indexOf(' — ') : undefined)}`);
      }
    }

    this.funcCallReplaced.next(true);
  }

  /**
   * Method loads corresponding FuncCall from DB if "id" param is provided in URL.
   * @stability Stable
  */
  protected async loadFuncCallById() {
    ui.setUpdateIndicator(this.root, true);

    const runId = this.getStartId();
    if (runId && !this.options.isTabbed) {
      this.linkFunccall(await historyUtils.loadRun(runId));
    } else {
      const func: DG.Func = await grok.functions.eval(this.funcName);
      this.linkFunccall(func.prepare({}));
    }

    ui.setUpdateIndicator(this.root, false);
  }

  /**
   * Method for any async logic that could not be placed in the constructor directly.
   * It is only called in the constructor, but not awaited.
   * A soon as {@link this.funcCall} is set, {@link this.onFuncCallReady} is emitted.
   * @stability Stable
 */
  public async init() {
    await this.loadFuncCallById();

    this.onFuncCallReady.complete();
  }

  /**
   * Override to create a fully custom UI including ribbon menus and panels
   * @stability Stable
 */
  public build(): void {
    ui.empty(this.root);
    this.root.appendChild(this.buildIO());

    if (this.options.historyEnabled) this.buildHistoryBlock();
    this.buildRibbonMenu();
    this.buildRibbonPanels();
  }

  /**
   * Override to create a custom input-output block.
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
    this.linkFunccall(savedCall);

    if (this.options.historyEnabled) this.buildHistoryBlock();
    if (!this.options.isTabbed) this.path = `?id=${savedCall.id}`;

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
    this.lastCall = pulledRun;
    await this.onAfterLoadRun(pulledRun);
    return pulledRun;
  }

  /**
   * Called before actual computations are made {@link run}.
   * @param funcToCall FuncCall object to be called {@see DG.FuncCall.call()}
   * @stability Stable
  */
  public async onBeforeRun(funcToCall: DG.FuncCall) {}

  /**
    * Called after actual computations are made {@link run}.
    * @param runFunc FuncCall object after call method {@see DG.FuncCall.call()}
    * @stability Stable
   */
  public async onAfterRun(runFunc: DG.FuncCall) {}

  /**
    * Called to perform actual computations.
    * @stability Stable
   */
  public async run(): Promise<void> {
    if (!this.funcCall) throw new Error('The correspoding function is not specified');

    await this.onBeforeRun(this.funcCall);
    const pi = DG.TaskBarProgressIndicator.create('Calculating...');
    this.funcCall.newId();
    await this.funcCall.call(); // CAUTION: mutates the funcCall field
    pi.close();
    await this.onAfterRun(this.funcCall);

    // If a view is incapuslated into a tab (e.g. in PipelineView),
    // there is no need to save run till an entire pipeline is over.
    this.lastCall = this.options.isTabbed ? this.funcCall.clone() : await this.saveRun(this.funcCall);
  }

  protected historyRoot: HTMLDivElement = ui.divV([], {style: {'justify-content': 'center'}});

  /**
    * Default export filename generation method.
    * It automatically replaces all symbols unsupported by Windows filesystem.
    * @param format A format listed in {@link defaultSupportedExportFormats}.
    * @stability Stable
   */
  protected defaultExportFilename = (format: string) => {
    return `${this.name} - ${new Date().toLocaleString('en-US').replaceAll(/:|\//g, '-')}.${this.exportConfig!.supportedExtensions[format]}`;
  };

  protected defaultSupportedExportExtensions: () => Record<string, string> = () => {
    return {
      'Excel': 'xlsx'
    };
  };

  protected defaultSupportedExportFormats = () => {
    return ['Excel'];
  };
}
