import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject} from 'rxjs';
import {CUSTOM_VIEW_OPTION, deepCopy} from '../../shared-utils/utils';
import {historyUtils} from '../../history-utils';
import type {RunMetadata, SaveTarget} from '../../history-utils/src/history-utils';
import {initNewCall} from '../../reactive-tree-driver/src/runtime/funccall-utils';

/** The Save and Open dialogs a host gives the view; each resolves null when the user cancels. */
export interface RunStorageHost {
  // where and under which name to save a run
  saveDialog: (call: DG.FuncCall) => Promise<{metadata: RunMetadata, target: SaveTarget} | null>;
  // a saved run of the function, loaded with its inputs and outputs
  openDialog: (func: DG.Func) => Promise<DG.FuncCall | null>;
  // the function that opens this view, recorded on saved runs so that they open here again
  customView?: string;
}

/** An item of the view's Export menu; the file it returns is downloaded under its own name. */
export interface ViewExport {
  name: string;
  export: () => Promise<DG.FileInfo>;
}

export abstract class CustomFunctionView extends DG.ViewBase {
  public isReady = new BehaviorSubject(false);

  public funcCall?: DG.FuncCall;

  private _runStorage?: RunStorageHost;

  private _exportEnabled = true;

  /** Set by the host that provides the Save and Open dialogs; without it the view has no History and Save buttons. */
  get runStorage(): RunStorageHost | undefined {
    return this._runStorage;
  }

  set runStorage(host: RunStorageHost | undefined) {
    this._runStorage = host;
    if (this.funcCall)
      this.buildRibbonPanels();
  }

  /** The Export menu; the ribbon's Export icon runs the first one. */
  exports: ViewExport[] = [];

  /** Off greys out the Export icon and menu, for runs the view can't export yet. */
  get exportEnabled(): boolean {
    return this._exportEnabled;
  }

  set exportEnabled(enabled: boolean) {
    this._exportEnabled = enabled;
    if (this.funcCall)
      this.buildRibbonPanels();
  }

  constructor(public funcNqName: string) {
    super();
    this.box = true;
    this.init();
  }

  // mandatory
  public abstract buildIO(): HTMLElement; // returns view content
  public abstract onAfterLoadRun(funcCall: DG.FuncCall): Promise<void>; // replacing funccall

  // optional hooks
  public async onBeforeSaveRun(callToSave: DG.FuncCall) { }
  public async onAfterSaveRun(savedCall: DG.FuncCall) { }

  getMocks: ({mockName: string, action: () => Promise<void>}[]) | null = null;

  getTemplates: ({name: string, action: () => Promise<void>}[]) | null = null;

  getHelp: (() => Promise<void>) | null = null;

  reportBug: (() => Promise<void>) | null = null;

  requestFeature: (() => Promise<void>) | null = null;

  public async init() {
    const call = DG.Func.byName(this.funcNqName).prepare({});
    await initNewCall(call, {});
    this.linkFunccall(call);
    await this.onFuncCallReady();
    this.isReady.next(true);
  }

  public linkFunccall(funcCall: DG.FuncCall) {
    this.funcCall = funcCall;
    this.path = funcCall.author ? `?id=${funcCall.id}` : '?';
    const modelName = funcCall?.func?.friendlyName ?? funcCall?.func?.name;
    this.name = modelName;
  }

  public async onFuncCallReady() {
    await this.getPackageData();
    this.build();
  }

  public async saveRun(callToSave: DG.FuncCall, metadata?: RunMetadata, target?: SaveTarget): Promise<DG.FuncCall> {
    const callCopy = deepCopy(callToSave);
    await this.onBeforeSaveRun(callCopy);

    if (callCopy.id) callCopy.newId();
    if (this.runStorage?.customView)
      callCopy.options[CUSTOM_VIEW_OPTION] = this.runStorage.customView;

    const savedCall = await historyUtils.saveRun(callCopy, {metadata, target});
    const loadedCall = await historyUtils.loadRun(savedCall.id);

    this.linkFunccall(loadedCall);

    await this.onAfterSaveRun(loadedCall);
    return loadedCall;
  }

  public async saveRunWithDialog(currentFuncCall: DG.FuncCall): Promise<DG.FuncCall | null> {
    try {
      const result = await this.runStorage?.saveDialog(currentFuncCall);
      return result ? await this.saveRun(currentFuncCall, result.metadata, result.target) : null;
    } catch (e: any) {
      grok.shell.error(`Couldn't save the run: ${e?.message ?? e}`);
      return null;
    }
  }

  /** Opens a saved run picked in the host's Open dialog. */
  public async openRunWithDialog(): Promise<DG.FuncCall | null> {
    try {
      const run = await this.runStorage?.openDialog(DG.Func.byName(this.funcNqName));
      if (!run)
        return null;
      this.linkFunccall(run);
      await this.onAfterLoadRun(run);
      return run;
    } catch (e: any) {
      grok.shell.error(`Couldn't open the run: ${e?.message ?? e}`);
      return null;
    }
  }

  public async loadRun(funcCallId: string): Promise<DG.FuncCall> {
    const pulledRun = await historyUtils.loadRun(funcCallId);

    this.linkFunccall(pulledRun);

    await this.onAfterLoadRun(pulledRun);
    return pulledRun;
  }

  public build(): void {
    ui.empty(this.root);
    const rootItem = ui.div([
      this.buildIO(),
    ], {style: {width: '100%', height: '100%', display: 'flex'}});
    this.root.appendChild(rootItem);

    this.buildRibbonMenu();
    this.buildRibbonPanels();
  }

  public buildRibbonPanels(): HTMLElement[][] {
    const historyButton = ui.iconFA('history', () => this.openRunWithDialog(), 'Open a saved run');
    const saveButton = ui.iconFA('save', () => this.funcCall ? this.saveRunWithDialog(this.funcCall) : undefined, 'Save');
    const exportButton = ui.iconFA('arrow-to-bottom', () => {
      if (this.exportEnabled)
        this.exportRun(this.exports[0]);
    }, this.exports[0]?.name);
    exportButton.classList.toggle('d4-disabled', !this.exportEnabled);
    const newRibbonPanels: HTMLElement[][] =
      [[...(this.runStorage ? [historyButton, saveButton] : []), ...(this.exports.length > 0 ? [exportButton] : [])]];
    this.setRibbonPanels(newRibbonPanels);
    return newRibbonPanels;
  }

  public buildRibbonMenu() {
    this.ribbonMenu.clear();

    if (this.exports.length === 0 && !this.reportBug && !this.requestFeature && !this.getHelp && !this.getMocks &&
      !this.getTemplates) return;

    const ribbonMenu = this.ribbonMenu.group('Model');

    if (this.getMocks && this.getMocks.length > 0) {
      if (this.getMocks.length === 1)
        ribbonMenu.item('Input data mock', this.getMocks[0].action);
      else {
        const dataGroup = ribbonMenu.group('Input data mocks');
        this.getMocks.forEach((val) => {
          dataGroup.item(val.mockName, val.action);
        });
        ribbonMenu.endGroup();
      }
    }

    if (this.getTemplates && this.getTemplates.length > 0) {
      if (this.getTemplates.length === 1)
        ribbonMenu.item('Input data template', this.getTemplates[0].action);
      else {
        const dataGroup = ribbonMenu.group('Input data templates');
        this.getTemplates.forEach((val) => {
          dataGroup.item(val.name, val.action);
        });
        ribbonMenu.endGroup();
      }
    }

    if (this.exports.length > 0) {
      ribbonMenu
        .group('Export')
        .items(this.exports, (item) => this.exportRun(item), {
          toString: (item) => item.name,
          isValid: () => this.exportEnabled ? null : 'Export is turned off for this run',
        })
        .endGroup();
    }

    if (this.reportBug)
      ribbonMenu.item('Report a bug', () => this.reportBug!());

    if (this.requestFeature)
      ribbonMenu.item('Request a feature', () => this.requestFeature!());

    if (this.getHelp)
      ribbonMenu.item('Help', () => this.getHelp!());
  }

  public async exportRun(item: ViewExport) {
    try {
      const file = await item.export();
      DG.Utils.download(file.fileName, new Blob([file.data as any]));
    } catch (e: any) {
      grok.shell.error(`Couldn't export the run: ${e?.message ?? e}`);
    }
  }

  private async getPackageData() {
    const pack = this.func?.package;
    if (!pack)
      return;

    const props = await pack.getProperties() as any;

    const reportBugUrl = props?.REPORT_BUG_URL;
    if (reportBugUrl && !this.reportBug)
      this.reportBug = async () => {window.open(reportBugUrl, '_blank');};

    const reqFeatureUrl = props?.REQUEST_FEATURE_URL;
    if (reqFeatureUrl && !this.requestFeature)
      this.requestFeature = async () => {window.open(reqFeatureUrl, '_blank');};
  }

  get func() {
    return this.funcCall?.func;
  }
}
