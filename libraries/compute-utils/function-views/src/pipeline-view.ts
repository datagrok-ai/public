/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import JSZip from 'jszip';
import {Subject} from 'rxjs';
import {filter} from 'rxjs/operators';
import {FunctionView} from './function-view';
import {ComputationView} from './computation-view';
import {historyUtils} from '../../history-utils';
import '../css/pipeline-view.css';

export class PipelineView extends ComputationView {
  public steps = {} as {[scriptNqName: string]: { editor: string, view: FunctionView }};
  public onStepCompleted = new Subject<DG.FuncCall>();

  private stepTabs: DG.TabControl | null = null;

  // PipelineView unites several export files into single ZIP file
  protected pipelineViewExportExtensions: () => Record<string, string> = () => {
    return {
      'Archive': 'zip',
    };
  };

  protected pipelineViewExportFormats = () => {
    return ['Archive'];
  };

  protected pipelineViewExport = async (format: string) => {
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

  exportConfig = {
    supportedExtensions: this.pipelineViewExportExtensions(),
    supportedFormats: this.pipelineViewExportFormats(),
    export: this.pipelineViewExport,
    filename: this.defaultExportFilename,
  };

  constructor(
    funcName: string,
    private stepsConfig: {funcName: string}[],
  ) {
    super(
      funcName,
      {historyEnabled: true, isTabbed: false},
    );
  }

  public override async init() {
    await this.loadFuncCallById();

    this.stepsConfig.forEach((stepConfig) => {
      //@ts-ignore
      this.steps[stepConfig.funcName] = {};
    });

    this.subs.push(
      grok.functions.onAfterRunAction.pipe(
        filter((run) => Object.keys(this.steps).includes(run.func.nqName)),
      ).subscribe((run) => {
        this.onStepCompleted.next(run);

        if (run.func.nqName === this.stepsConfig[this.stepsConfig.length-1].funcName) this.run();
      }),
    );

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

      this.steps[loadedScript.nqName].view =
        await editorFuncs[this.steps[loadedScript.nqName].editor].apply({'call': scriptCall}) as FunctionView;

      if (!this.steps[loadedScript.nqName].view.onFuncCallReady.value) {
        const prom = this.steps[loadedScript.nqName].view.onFuncCallReady.toPromise();
        return prom;
      } else return Promise.resolve();
    });

    await Promise.all(viewsLoading);

    this.onFuncCallReady.complete();
  }

  public override buildIO() {
    const tabs = Object.entries(this.steps)
      .reduce((prev, [funcName, step]) => ({
        ...prev,
        [funcName]: step.view.root,
      }), {} as Record<string, HTMLElement>);

    const pipelineTabs = ui.tabControl(tabs);

    const tabsLine = pipelineTabs.panes[0].header.parentElement!;
    tabsLine.classList.add('d4-ribbon', 'pipeline-view');
    tabsLine.classList.remove('d4-tab-header-stripe');
    tabsLine.firstChild!.remove();
    for (let i = 0; i < pipelineTabs.panes.length; i++) {
      pipelineTabs.panes[i].header.classList.add('d4-ribbon-name');
      pipelineTabs.panes[i].header.classList.remove('d4-tab-header');
    }
    pipelineTabs.panes[0].header.style.marginLeft = '12px';

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
      .forEach(async (step) => {
        const scriptCall = step.view.funcCall;

        scriptCall.options['parentCallId'] = this.funcCall!.id;

        this.steps[scriptCall.func.nqName].view.lastCall =
          await this.steps[scriptCall.func.nqName].view.saveRun(scriptCall);
      });

    await this.onAfterRun(this.funcCall);

    this.lastCall = await this.saveRun(this.funcCall);
  }

  /**
   * Overrided to use {@link loadChildRuns} during run load.
   * This implementation takes "parentCallId" and looks for the funcCalls with options.parentCallId = parentCallId.
   * Each child run is related to the particular pipeline step.
   * @param funcCallId ID of the parent FuncCall
   * @returns Parent FuncCall
   */
  public async loadRun(funcCallId: string): Promise<DG.FuncCall> {
    const {parentRun: pulledParentRun, childRuns: pulledChildRuns} = await historyUtils.loadChildRuns(funcCallId);

    this.subs.push(
      this.onFuncCallReady.subscribe({
        complete: async () => {
          await this.onBeforeLoadRun();

          pulledChildRuns.forEach(async (pulledChildRun) => {
            this.steps[pulledChildRun.func.nqName].view.loadRun(pulledChildRun.id);
          });

          await this.onAfterLoadRun(pulledParentRun);
        },
      }),
    );
    return pulledParentRun;
  }
}
