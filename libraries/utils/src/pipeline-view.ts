/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import JSZip from 'jszip';
import {Subject} from 'rxjs';
import {historyUtils} from './history-utils';
import {FunctionView} from './function-view';
import {ComputationView} from './computation-view';

export class PipelineView extends ComputationView {
  public steps = {} as {[scriptNqName: string]: { funcCall: DG.FuncCall, editor: string, view: FunctionView }};
  public get stepFuncNames() {
    return this.stepsConfig.map((step) => step.funcName);
  }
  public onStepCompleted = new Subject<DG.FuncCall>();

  private stepTabs: DG.TabControl | null = null;

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

    for (const {nqName, stepView} of Object.entries(this.steps)
      .map(([nqName, step]) => ({nqName, stepView: step.view}))) {
      this.stepTabs.currentPane = this.stepTabs?.getPane(nqName);
      await new Promise((r) => setTimeout(r, 100));
      const stepBlob = await stepView.exportConfig!.export('Excel');

      zip.file(stepView.exportConfig!.filename('Excel'), stepBlob, {binary: true, createFolders: false});
    };

    return await zip.generateAsync({type: 'blob'});
  };

  constructor(
    funcName: string,
    private stepsConfig: {funcName: string}[]
  ) {
    super(funcName);

    this.exportConfig = {
      supportedExtensions: this.defaultSupportedExportExtensions(),
      supportedFormats: this.defaultSupportedExportFormats(),
      export: this.defaultExport,
      filename: this.defaultExportFilename,
    };

    stepsConfig.forEach((stepConfig) => {
      //@ts-ignore
      this.steps[stepConfig.funcName] = {};
    });
  }

  public override async init() {
    const stepScripts = Object.keys(this.steps).map((stepNqName) => {
      const stepScript = (grok.functions.eval(stepNqName) as Promise<DG.Func>);
      return stepScript;
    });
    const allStepsLoading = Promise.all(stepScripts) as Promise<DG.Script[]>;
    this.root.classList.remove('ui-panel');

    const loadedScripts = await allStepsLoading;

    const editorFuncs = {} as {[editor: string]: DG.Func};

    const EDITOR_TAG = 'editor:' as const;
    const NEWLINE = '\n' as const;
    const DEFAULT_EDITOR = 'Compute:PipelineStepEditor';
    const extractEditor = (script: DG.Script) => {
      const scriptCode = script.script;
      const editorTagIndex = scriptCode.indexOf(EDITOR_TAG);
      if (editorTagIndex < 0)
        return DEFAULT_EDITOR;

      const newlineIndex = scriptCode.indexOf(NEWLINE, editorTagIndex);
      const editorFuncName = scriptCode.substring(editorTagIndex + EDITOR_TAG.length, newlineIndex).trim();

      return editorFuncName;
    };

    const editorsLoading = loadedScripts.map(async (loadedScript) => {
      // TO DO: replace for type guard
      const editorName = (loadedScript.script) ? extractEditor(loadedScript): DEFAULT_EDITOR;
      if (!editorFuncs[editorName])
        editorFuncs[editorName] = await(grok.functions.eval(editorName.split(' ').join('')) as Promise<DG.Func>);
      this.steps[loadedScript.nqName].editor = editorName;

      return Promise.resolve();
    });

    await Promise.all(editorsLoading);

    const viewsLoading = loadedScripts.map(async (loadedScript) => {
      const scriptCall: DG.FuncCall = loadedScript.prepare();

      this.steps[loadedScript.nqName].funcCall = scriptCall;
      this.steps[loadedScript.nqName].view =
        await editorFuncs[this.steps[loadedScript.nqName].editor].apply({'call': scriptCall}) as FunctionView;

      return Promise.resolve();
    });

    await Promise.all(viewsLoading);
  }

  public override buildIO() {
    const tabs = Object.entries(this.steps)
      .reduce((prev, [funcName, step]) => ({
        ...prev,
        [funcName]: step.view.root
      }), {} as Record<string, HTMLElement>);

    const pipelineTabs = ui.tabControl(tabs);

    for (let i = 0; i < pipelineTabs.panes.length - 1; i++) {
      pipelineTabs.panes[i].header.classList.add('arrow-tab');
      pipelineTabs.panes[i].header.insertAdjacentElement('afterend', ui.div(undefined, 'empty-box'));
    }
    pipelineTabs.root.style.height = '100%';
    pipelineTabs.root.style.width = '100%';

    this.stepTabs = pipelineTabs;

    return pipelineTabs.root;
  }

  public override async run(): Promise<void> {
    if (!this.funcCall) throw new Error('The correspoding function is not specified');

    await this.onBeforeRun(this.funcCall);
    const pi = DG.TaskBarProgressIndicator.create('Calculating...');
    this.funcCall.newId();
    await this.funcCall.call(); // mutates the funcCall field
    pi.close();

    Object.values(this.steps)
      .map((step) => step.funcCall)
      .forEach(async (scriptCall) => {
        scriptCall.options['parentCallId'] = this.funcCall!.id;
        await this.steps[scriptCall.func.nqName].view.saveRun(scriptCall);
      });

    await this.onAfterRun(this.funcCall);

    this.lastCall = await this.saveRun(this.funcCall);
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

    pulledChildRuns.forEach(async (pulledChildRun) => {
      await this.steps[pulledChildRun.func.nqName].view.onAfterLoadRun(pulledChildRun);
      this.steps[pulledChildRun.func.nqName].view.linkFunccall(pulledChildRun);
    });

    await this.onAfterLoadRun(pulledParentRun);
    return pulledParentRun;
  }
}
