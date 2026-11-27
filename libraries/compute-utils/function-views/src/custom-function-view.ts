import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {BehaviorSubject} from 'rxjs';
import {CUSTOM_VIEW_OPTION, deepCopy} from '../../shared-utils/utils';
import {historyUtils} from '../../history-utils';
import type {RunMetadata, SaveTarget} from '../../history-utils/src/history-utils';
import {initNewCall, OUTPUT_OUTDATED_PATH, readOutputOutdated, RUN_ERROR_PATH} from '../../reactive-tree-driver/src/runtime/funccall-utils';

/** Whether Compute2, which has the Save and Open dialogs, is installed. */
export const hasRunDialogs = () => DG.Func.find({package: 'Compute2', name: 'OpenRunDialog'}).length > 0;

/** Compute2's Open dialog: a saved run of the function, loaded with its inputs and outputs, or null. */
export const openRunDialog = (func: DG.Func): Promise<DG.FuncCall | null> =>
  grok.functions.call('Compute2:OpenRunDialog', {func});

/** Compute2's Save dialog: where and under which name to save the run, or null. */
const saveRunDialog = (call: DG.FuncCall): Promise<{metadata: RunMetadata, target: SaveTarget} | null> =>
  grok.functions.call('Compute2:SaveRunDialog', {call});

/** An item of the view's Export menu; the file it returns is downloaded under its own name. */
export interface ViewExport {
  name: string;
  export: () => Promise<DG.FileInfo>;
}

const OUTDATED_MSG = 'Output outdated';
const EXPORT_DISABLED_MSG = 'Export is turned off for this run';

export abstract class CustomFunctionView extends DG.ViewBase {
  public isReady = new BehaviorSubject(false);

  public funcCall?: DG.FuncCall;

  private _exportEnabled = true;
  private _isOutputOutdated: boolean | undefined;
  private _runError: string | undefined;

  /** The function that opens this view, recorded on saved runs so that they open here again; set by the editor. */
  public customView?: string;

  /** The Export menu; the ribbon's Export icon runs the first one. */
  exports: ViewExport[] = [];

  get exportEnabled(): boolean {
    return this._exportEnabled;
  }

  set exportEnabled(enabled: boolean) {
    this._exportEnabled = enabled;
  }

  get outputOutdated(): boolean | undefined {
    return this._isOutputOutdated;
  }

  set outputOutdated(isOutdated: boolean) {
    this._isOutputOutdated = isOutdated;
  }

  get runError(): string | undefined {
    return this._runError;
  }

  set runError(error: string) {
    this._runError = error;
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
    const outdatedRaw = funcCall.options[OUTPUT_OUTDATED_PATH];
    this._isOutputOutdated = outdatedRaw ? JSON.parse(outdatedRaw) : undefined;
    this._runError = funcCall.options[RUN_ERROR_PATH] || undefined;
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
    if (this.customView)
      callCopy.options[CUSTOM_VIEW_OPTION] = this.customView;
    if (this._isOutputOutdated != null)
      callCopy.options[OUTPUT_OUTDATED_PATH] = this._isOutputOutdated;
    if (this._runError != null)
      callCopy.options[RUN_ERROR_PATH] = this._runError;

    const savedCall = await historyUtils.saveRun(callCopy, {metadata, target});
    const loadedCall = await historyUtils.loadRun(savedCall.id);

    this.linkFunccall(loadedCall);

    await this.onAfterSaveRun(loadedCall);
    return loadedCall;
  }

  public async saveRunWithDialog(currentFuncCall: DG.FuncCall): Promise<DG.FuncCall | null> {
    try {
      const result = await saveRunDialog(currentFuncCall);
      return result ? await this.saveRun(currentFuncCall, result.metadata, result.target) : null;
    } catch (e: any) {
      grok.shell.error(`Couldn't save the run: ${e?.message ?? e}`);
      return null;
    }
  }

  /** Opens a saved run picked in Compute2's Open dialog. */
  public async openRunWithDialog(): Promise<DG.FuncCall | null> {
    try {
      const run = await openRunDialog(DG.Func.byName(this.funcNqName));
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
    this.root.appendChild(this.buildIO());

    this.buildRibbonMenu();
    this.buildRibbonPanels();
  }

  public buildRibbonPanels(): HTMLElement[][] {
    const historyButton = ui.iconFA('history', () => this.openRunWithDialog(), 'Open a saved run');
    const saveButton = ui.iconFA('save', () => this.funcCall ? this.saveRunWithDialog(this.funcCall) : undefined, 'Save');
    const exportButton = ui.iconFA('arrow-to-bottom', () => {
      if (this.outputOutdated)
        grok.shell.error(OUTDATED_MSG)
      else if (!this.exportEnabled)
        grok.shell.error(EXPORT_DISABLED_MSG)
      else
        this.exportRun(this.exports[0]);
    }, this.exports[0]?.name);
    const newRibbonPanels: HTMLElement[][] =
      [[...(hasRunDialogs() ? [historyButton, saveButton] : []), ...(this.exports.length > 0 ? [exportButton] : [])]];
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
          isValid: () => this.outputOutdated ? OUTDATED_MSG :
            (this.exportEnabled ? null : EXPORT_DISABLED_MSG),
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
