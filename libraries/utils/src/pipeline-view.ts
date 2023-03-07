/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {RichFunctionView} from './rich-function-view';
import {UiUtils} from './shared-components/ui-utils';
import {historyUtils} from './history-utils';
import {FunctionView} from './function-view';
import JSZip from 'jszip';

export class PipelineView extends RichFunctionView {
  public stepViews: { [nqName: string]: FunctionView } = {};
  public stepTabs: DG.TabControl | null = null;

  protected defaultExportFilename = (format: string) => {
    return `${this.name} - ${new Date()
      .toLocaleString('en-US')
      .replaceAll(/:|\//g, '-')}.${this.exportConfig!.supportedExtensions[format]}`;
  };

  protected defaultSupportedExportExtensions: () => Record<string, string> = () => {
    return {
      'Archive': 'zip'
    };
  };

  protected defaultSupportedExportFormats = () => {
    return ['Archive'];
  };

  protected defaultExport = async (format: string) => {
    if (format !== 'Archive')
      throw new Error('This export format is not supported');

    if (!this.stepTabs)
      throw new Error('Set step tabs please for export');

    const zip = new JSZip();

    for (const [nqName, stepView] of Object.entries(this.stepViews)) {
      this.stepTabs.currentPane = this.stepTabs?.getPane(nqName);
      await new Promise((r) => setTimeout(r, 100));
      const stepBlob = await stepView.exportConfig!.export('Excel');

      zip.file(stepView.exportConfig!.filename('Excel'), stepBlob, {binary: true, createFolders: false});
    };

    return await zip.generateAsync({type: 'blob'});
  };

  constructor(
    funcCall: DG.FuncCall
  ) {
    super(funcCall, {exportEnabled: true, historyEnabled: false, isTabbed: true});

    this.exportConfig = {
      supportedExtensions: this.defaultSupportedExportExtensions(),
      supportedFormats: this.defaultSupportedExportFormats(),
      export: this.defaultExport,
      filename: this.defaultExportFilename,
    };
  }

  // REVIEW INHERITANCE HIERARCHY
  public override buildIO() {
    return ui.div();
  }

  // Should override to prevent saving right after model loading
  // Actually, it does nothing
  public async run(): Promise<void> {
  }

  /**
   * Override to create a custom historical runs control.
   * @returns The HTMLElement with history block UI
  */
  // Overrrided since we should load child funcCalls
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
   * Loads the specified historical run. See also {@link saveRun}.
   * @param funcCallId ID of FuncCall to look for. Get it using {@see funcCall.id} field
   * @returns FuncCall augemented with inputs' and outputs' values
   * @stability Stable
 */
  public async loadRun(funcCallId: string): Promise<DG.FuncCall> {
    await this.onBeforeLoadRun();
    const {parentRun: pulledParentRun, childRuns: pulledChildRuns} = await historyUtils.loadChildRuns(funcCallId);

    pulledChildRuns.forEach((pulledChildRun) => {
      this.stepViews[pulledChildRun.nqName].linkFunccall(pulledChildRun);
    });

    await this.onAfterLoadRun(pulledParentRun);
    return pulledParentRun;
  }
}
