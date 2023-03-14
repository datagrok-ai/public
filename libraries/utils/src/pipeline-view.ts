/* eslint-disable valid-jsdoc */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {from} from 'rxjs';
import {filter} from 'rxjs/operators';
import JSZip from 'jszip';
import {historyUtils} from './history-utils';
import {FunctionView} from './function-view';
import {ComputationView} from './computation-view';

export class PipelineView extends ComputationView {
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

  scriptCalls = {} as {[scriptNqName: string]: DG.FuncCall};
  editorFuncs = {} as {[editor: string]: DG.Func};
  scriptEditors = {} as {[scriptName: string]: string};
  views = {} as {[script: string]: HTMLElement};
  scripts: DG.Script[] = [];

  constructor(
    funcName: string,
    private scriptStepNames: string[]
  ) {
    super(funcName);

    this.exportConfig = {
      supportedExtensions: this.defaultSupportedExportExtensions(),
      supportedFormats: this.defaultSupportedExportFormats(),
      export: this.defaultExport,
      filename: this.defaultExportFilename,
    };
  }

  public override async init() {
    const scriptStepNames = this.scriptStepNames;

    const stepScripts = scriptStepNames.map((stepScriptName) => {
      const stepScript = (grok.functions.eval(stepScriptName) as Promise<DG.Func>);
      return stepScript;
    });
    const allStepsLoaded = Promise.all(stepScripts) as Promise<DG.Script[]>;

    this.root.classList.remove('ui-panel');

    this.scripts = await allStepsLoaded;

    const editorTag = 'editor:' as const;
    const newlineSymbol = '\n' as const;
    const defaultEditor = 'Compute:PipelineStepEditor';
    const extractEditor = (script: DG.Script) => {
      const scriptCode = script.script;
      const editorTagIndex = scriptCode.indexOf(editorTag);
      if (editorTagIndex < 0)
        return defaultEditor;

      const newlineIndex = scriptCode.indexOf(newlineSymbol, editorTagIndex);
      const editorFuncName = scriptCode.substring(editorTagIndex + editorTag.length, newlineIndex).trim();

      return editorFuncName;
    };

    for (const script of this.scripts) {
      // TO DO: replace for type guard
      const editorName = (script.script) ? extractEditor(script): defaultEditor;

      if (!this.editorFuncs[editorName])
        this.editorFuncs[editorName] = await(grok.functions.eval(editorName.split(' ').join('')) as Promise<DG.Func>);

      this.scriptEditors[script.nqName] = editorName;
    };

    for (const script of this.scripts) {
      const scriptCall = script.prepare();
      const view: FunctionView = await this.editorFuncs[this.scriptEditors[script.nqName]].apply({'call': scriptCall});

      // scriptCall.options['parentCallId'] = wrapperCall.id;
      this.scriptCalls[script.nqName] = scriptCall;
      this.views[script.nqName] = view.root;

      // TEMP WAY TO DEAL WITH IT
      this.stepViews[script.nqName] = view;
    }
  }

  public override buildIO() {
    const pipelineTabs = ui.tabControl(this.views);
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

    Object.values(this.scriptCalls).forEach(async (scriptCall) => {
      scriptCall.options['parentCallId'] = this.funcCall!.id;
      await this.stepViews[scriptCall.func.nqName].saveRun(scriptCall);
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
      await this.stepViews[pulledChildRun.func.nqName].onAfterLoadRun(pulledChildRun);
      this.stepViews[pulledChildRun.func.nqName].linkFunccall(pulledChildRun);
    });

    await this.onAfterLoadRun(pulledParentRun);
    return pulledParentRun;
  }
}
