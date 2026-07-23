/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

import {FuncFlowView} from './funcflow-view';
import {FlowEntityHandler} from './entity/flow-entity-handler';
import {parseFlowBody, FLOW_LANGUAGE} from './serialization/flow-script-format';
import { getFilesBrowser } from './utils/files-browser-tree';
import {readUploadedFileBytes, parseFileToDataFrame, syncFlowFilePermissions} from './utils/uploaded-files';
import * as aiTools from './ai-tools';

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

/* The node behind local-file uploads: dropping a file onto the canvas stores
 * its bytes (in memory until the flow is saved, then in the server's
 * GUID-addressed file store) and adds this function as a node, so the flow
 * replays and shares like any other. See utils/uploaded-files.ts. */
//name: readUploadedFile
//friendlyName: Uploaded File
//description: Reads a file uploaded into a flow and parses it into a table
//input: string fileId
//input: string fileName
//output: dataframe result
//meta.includeInFlow: true
export async function readUploadedFile(fileId: string, fileName: string): Promise<DG.DataFrame> {
  const bytes = await readUploadedFileBytes(fileId, fileName);
  return parseFileToDataFrame(fileName, bytes);
}

/* Runs at platform startup (not just when a Flow view is open): sharing a flow
 * script from anywhere — Browse, a link, the context panel — must extend read
 * access to the uploaded-file blobs its nodes reference. */
//name: flowShareSync
//tags: autostart
//description: Keeps uploaded-file permissions in sync when a flow script is shared
//meta.includeInFlow: false
export function flowShareSync(): void {
  grok.events.onEntityShared.subscribe((e) => {
    if (e instanceof DG.Script && (e.language as string) === FLOW_LANGUAGE)
      void syncSharedFlow(e.id);
  });
}

/** The shared entity from the dialog can be shallow — re-fetch for the body. */
async function syncSharedFlow(id: string): Promise<void> {
  const script = await grok.dapi.scripts.find(id).catch(() => null);
  if (script?.script)
    await syncFlowFilePermissions(script);
}

export class PackageFunctions {
  @grok.decorators.app({
    name: 'Flow',
    description: 'Interactive function chain designer',
    tags: ['app'],
  })
  static funcflowApp(@grok.decorators.param({options: {metaUrl: true, optional: true}}) path?: string): DG.ViewBase {
    return new FuncFlowView();
  }

  @grok.decorators.fileViewer({fileViewer: 'ffjson'})
  static viewFuncFlow(file: DG.FileInfo): DG.ViewBase {
    const view = new FuncFlowView();
    file.readAsString().then((json) => view.loadFromJson(json));
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
    // No output panel inside the dialog — run results belong to the real
    // editor view only; it is re-enabled below when promoted via Open In Editor.
    const view = new FuncFlowView(tableInfos, {outputPanel: false});
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
        view.enableOutputPanel();
        grok.shell.addView(view);
        setTimeout(() => view.fitToScreen(), 100);
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

  // ---------- Flow view functions (AI) ----------
  // Returned by FuncFlowView.getFunctions() (found by the 'flowViewFunction' tag) so the
  // AI assistant can act on the open editor. Each takes the generic current view and
  // reaches the FuncFlowView instance through `view.jsView`; `includeInFlow: false`
  // keeps them out of Flow's own node catalog.

  @grok.decorators.func({tags: ['flowViewFunction'], meta: {includeInFlow: 'false'},
    description: 'List the current flow graph: all nodes (id, label, type, status, set input values) and connections. Call this first to understand what is on the canvas'})
  static listFlowNodes(@grok.decorators.param({type: 'view'}) view: any): any {
    return aiTools.listFlowNodes(view);
  }

  @grok.decorators.func({tags: ['flowViewFunction'], meta: {includeInFlow: 'false'},
    description: 'Ports (with DG types), editable input values, unmet requirements, and last-run outputs of one node'})
  static getFlowNodeDetails(@grok.decorators.param({type: 'view'}) view: any, nodeId: string): any {
    return aiTools.getFlowNodeDetails(view, nodeId);
  }

  @grok.decorators.func({tags: ['flowViewFunction'], meta: {includeInFlow: 'false'},
    description: 'Search the flow node catalog (a curated subset of platform functions plus input/output/utility nodes). ALWAYS filter: pass a query with what the node should do (e.g. "join tables", "open file"), and/or a DG type it must accept or produce (dataframe, column, string, ...). Returns at most limit (default 15) matches with their input/output types'})
  static findFlowNodeTypes(
    @grok.decorators.param({type: 'view'}) view: any,
    @grok.decorators.param({type: 'string', options: {optional: true, description: 'Words describing what the node does'}}) query?: string,
    @grok.decorators.param({type: 'string', options: {optional: true, description: 'DG type one of its inputs must accept'}}) acceptsInputType?: string,
    @grok.decorators.param({type: 'string', options: {optional: true, description: 'DG type one of its outputs must produce'}}) producesOutputType?: string,
    @grok.decorators.param({type: 'int', options: {optional: true}}) limit?: number): any {
    return aiTools.findFlowNodeTypes(view, query, acceptsInputType, producesOutputType, limit);
  }

  @grok.decorators.func({tags: ['flowViewFunction'], meta: {includeInFlow: 'false'},
    description: 'Add a node to the canvas by its registered typeName (from findFlowNodeTypes). Optionally set editable input values right away. Returns the new node id and its ports'})
  static async addFlowNode(
    @grok.decorators.param({type: 'view'}) view: any,
    typeName: string,
    @grok.decorators.param({type: 'string', options: {optional: true, description: 'Optional custom title'}}) label?: string,
    @grok.decorators.param({type: 'map', options: {optional: true, description: 'Editable primitive inputs, key to value'}}) inputValues?: object): Promise<any> {
    return aiTools.addFlowNode(view, typeName, label, inputValues);
  }

  @grok.decorators.func({tags: ['flowViewFunction'], meta: {includeInFlow: 'false'},
    description: 'Connect a source node output to a target node input (port keys from getFlowNodeDetails / addFlowNode). Types must be compatible'})
  static async connectFlowNodes(
    @grok.decorators.param({type: 'view'}) view: any,
    sourceNodeId: string, sourceOutput: string, targetNodeId: string, targetInput: string): Promise<any> {
    return aiTools.connectFlowNodes(view, sourceNodeId, sourceOutput, targetNodeId, targetInput);
  }

  @grok.decorators.func({tags: ['flowViewFunction'], meta: {includeInFlow: 'false'},
    description: 'Set editable input values of a node (key to value; keys from getFlowNodeDetails). Marks the node and its downstream stale'})
  static async setFlowNodeInputs(
    @grok.decorators.param({type: 'view'}) view: any,
    nodeId: string,
    @grok.decorators.param({type: 'map'}) values: object): Promise<any> {
    return aiTools.setFlowNodeInputs(view, nodeId, values);
  }

  @grok.decorators.func({tags: ['flowViewFunction'], meta: {includeInFlow: 'false'},
    description: 'Select a node on the canvas so the user sees it (opens its properties panel)'})
  static async selectFlowNode(@grok.decorators.param({type: 'view'}) view: any, nodeId: string): Promise<any> {
    return aiTools.selectFlowNode(view, nodeId);
  }

  @grok.decorators.func({tags: ['flowViewFunction'], meta: {includeInFlow: 'false'},
    description: 'List Flow built-in interactive guides: step-by-step tutorials and short "how do I" walkthroughs that highlight the actual UI. When the user asks how to do something in Flow, check here first — a matching guide beats a textual explanation'})
  static listFlowGuides(
    @grok.decorators.param({type: 'view'}) view: any,
    @grok.decorators.param({type: 'string', options: {optional: true, description: 'Words to filter by'}}) query?: string): any {
    return aiTools.listFlowGuides(view, query);
  }

  @grok.decorators.func({tags: ['flowViewFunction'], meta: {includeInFlow: 'false'},
    description: 'Start an interactive guide (id from listFlowGuides) — it highlights the real UI step by step and waits for the user to act. ALWAYS confirm with the user first before starting one; never launch it unasked'})
  static startFlowGuide(@grok.decorators.param({type: 'view'}) view: any, guideId: string): any {
    return aiTools.startFlowGuide(view, guideId);
  }

  @grok.decorators.func({tags: ['flowViewFunction'], meta: {includeInFlow: 'false'},
    description: 'Validate and execute the whole flow. Returns validation problems instead of running if the graph is invalid; otherwise waits for the run and reports per-node failures'})
  static async runFlow(@grok.decorators.param({type: 'view'}) view: any): Promise<any> {
    return aiTools.runFlow(view);
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
