/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

import {FuncFlowView} from './funcflow-view';
import {FlowEntityHandler} from './entity/flow-entity-handler';
import {parseFlowBody} from './serialization/flow-script-format';
import { getFilesBrowser } from './utils/files-browser-tree';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

/* The 'flow' script-language handler: core wraps this into a
 * PackageScriptHandler at func sync, so flow scripts run like any script
 * (grok.functions.call, Run pane, funccall dialog). Registered via plain
 * annotations (not decorators) because the templateScript meta value —
 * which itself contains '//' header lines — does not survive the decorator
 * code generator. */
//name: flowScriptHandler
//input: funccall scriptCall
//meta.role: scriptHandler
//meta.scriptHandler.language: flow
//meta.scriptHandler.extensions: flow
//meta.scriptHandler.commentStart: //
//meta.scriptHandler.codeEditorMode: javascript
//meta.scriptHandler.editorFunction: Flow:flowScriptEditor
//meta.scriptHandler.templateScript: //name: New Flow\n//language: flow\n//tags: flow\n{"version":"2.0","name":"New Flow","description":"","author":"","created":"","modified":"","nodes":[],"connections":[],"metadata":{"settings":{"scriptName":"New Flow","scriptDescription":"","tags":["flow"]}}}
//meta.icon: package.png
//meta.includeInFlow: false
export async function flowScriptHandler(scriptCall: DG.FuncCall): Promise<void> {
  await FlowEntityHandler.instance.run(scriptCall);
}

export class PackageFunctions {
  @grok.decorators.app({
    name: 'Flow',
    description: 'Interactive function chain designer',
    tags: ['app'],
  })
  static funcflowApp(@grok.decorators.param({options: {metaUrl: true, optional: true}}) path?: string): DG.ViewBase {
    const url = new URL(window.location.href);
    const params = url.searchParams;
    console.log(params);
    setTimeout(() => {
      grok.shell.windows.showToolbox = false;
      grok.shell.windows.showBrowse = false;
      grok.shell.windows.showContextPanel = true;
      grok.shell.windows.showHelp = false;
      grok.shell.windows.showBrowse = true;
      grok.shell.windows.showToolbox = true;
    }, 200);
    return new FuncFlowView();
  }

  @grok.decorators.fileViewer({fileViewer: 'ffjson'})
  static viewFuncFlow(file: DG.FileInfo): DG.ViewBase {
    const view = new FuncFlowView();
    file.readAsString().then((json) => view.loadFromJson(json));
    setTimeout(() => {
      grok.shell.windows.showToolbox = false;
      grok.shell.windows.showBrowse = false;
      grok.shell.windows.showContextPanel = true;
      grok.shell.windows.showHelp = false;
      grok.shell.windows.showBrowse = true;
      grok.shell.windows.showToolbox = true;
    }, 200);
    return view;
  }

  /** Builds a flow from a table-creation script (the function-call cascade
   *  Datagrok records for reproducibly-created tables, used by data sync)
   *  and opens it in the Flow editor. */
  @grok.decorators.func({
    name: 'flowFromCreationScript',
    description: 'Builds a flow diagram from a table creation script and opens it in the Flow editor',
  })
  static async flowFromCreationScript(script: string): Promise<DG.ViewBase> {
    const view = new FuncFlowView();
    await view.loadFromCreationScript(script);
    return view;
  }

  @grok.decorators.func({
    name: 'openCreationScriptFlowDialog',
    // includeInFlow: Flow-internal dialog opener — hide it from Flow's own toolbox.
    meta: {role: 'creationScriptEditor', includeInFlow: 'false'},
  })
  static async openCreationScriptFlowDialog(script: string, tableIds: string[], show: boolean = true): Promise<DG.Dialog> {
    // Load the tables being edited so the view can split the flow back into a
    // creation script per table and save each via TableInfo.saveCreationScript.
    const loaded = await Promise.all((tableIds ?? []).map((id) => grok.dapi.tables.find(id)));
    const tableInfos = loaded.filter((t): t is DG.TableInfo => t != null);
    const view = new FuncFlowView(tableInfos);
    view.name = `Creation Script`;
    // Inside the cramped dialog the overview adds clutter — start it minimized;
    // expand it once the flow is opened in the full editor.
    view.setMinimapCollapsed(true);
    try {
      await view.loadFromCreationScript(script);
    } catch (e) {
      grok.shell.error(`Failed to load flow from creation script`);
      console.error(e);
    }
    const d = ui.dialog({title: 'Creation Script Flow'})
      .add(view.root)
      .addButton('Open In Editor', () => {
        view.setMinimapCollapsed(false);
        grok.shell.addView(view);
        d.close();
      });
    if (show)
      d.show({resizable: true, width: 800, height: 600});
    return d;
  }

  @grok.decorators.func()
  static testDialog() {
    ui.dialog().add(getFilesBrowser((n) => {console.log(n.name)}, (n) => {console.log('dblclick', n.name)}, 'test-dialog-files').root).show();
  }

  // ---------- first-class Flow entity (Script with language 'flow') ----------
  // (the scriptHandler function itself is annotation-registered above,
  //  next to `info` — see the note there)

  /** The visual editor for a flow script entity — consumed by core through the
   *  `scriptHandler.editorFunction` seam (double-click, Edit, /script/<id>). */
  @grok.decorators.func({
    name: 'flowScriptEditor',
    description: 'Opens the visual Flow editor for a flow script entity',
    meta: {includeInFlow: 'false'},
  })
  static flowScriptEditor(
    @grok.decorators.param({type: 'script'}) script: DG.Script): DG.ViewBase {
    return FlowEntityHandler.instance.editorView(script);
  }

  /** Browse-preview view for a flow script entity (FlowScriptMeta.renderPreview). */
  @grok.decorators.func({
    name: 'flowScriptPreview',
    meta: {includeInFlow: 'false'},
  })
  static flowScriptPreview(
    @grok.decorators.param({type: 'script'}) script: DG.Script): DG.ViewBase {
    return FlowEntityHandler.instance.previewView(script);
  }

  /** Context-panel pane content for a flow script entity (FlowScriptMeta). */
  @grok.decorators.func({
    name: 'flowScriptWidget',
    meta: {includeInFlow: 'false'},
  })
  static flowScriptWidget(
    @grok.decorators.param({type: 'script'}) script: DG.Script): DG.Widget {
    return FlowEntityHandler.instance.widget(script);
  }

  /** `.flow` exports sitting in file shares open in the editor too. */
  @grok.decorators.fileViewer({fileViewer: 'flow'})
  static viewFlowFile(file: DG.FileInfo): DG.ViewBase {
    const view = new FuncFlowView();
    file.readAsString()
      .then((text) => view.loadFromDoc(parseFlowBody(text).doc))
      .catch((e) => grok.shell.error(`Cannot open flow file: ${e?.message ?? e}`));
    return view;
  }
}
