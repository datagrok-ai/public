/* eslint-disable max-len */
/* Do not change these import lines to match external modules in webpack configuration */
import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
export * from './package.g';

import {FuncFlowView} from './funcflow-view';
import {FlowEntityHandler} from './entity/flow-entity-handler';
import {parseFlowBody, FLOW_LANGUAGE} from './serialization/flow-script-format';
import {readUploadedFileBytes, parseFileToDataFrame, syncFlowFilePermissions} from './utils/uploaded-files';
import * as aiTools from './ai-tools';
import * as dataOps from './ops/data-ops';

export const _package = new DG.Package();

//name: info
export function info() {
  grok.shell.info(_package.webRoot);
}

/* Plain annotations, not decorators: the templateScript meta value contains
 * '//' header lines that don't survive the decorator code generator. */
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

//name: readUploadedFile
//friendlyName: Uploaded File
//description: Reads a file uploaded into a flow and parses it into a table
//input: string fileId
//input: string fileName
//output: dataframe result
//meta.includeInFlow: true
//meta.autorun: true
export async function readUploadedFile(fileId: string, fileName: string): Promise<DG.DataFrame> {
  const bytes = await readUploadedFileBytes(fileId, fileName);
  return parseFileToDataFrame(fileName, bytes);
}

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

  @grok.decorators.fileViewer({fileViewer: 'flow'})
  static viewFuncFlow(file: DG.FileInfo): DG.ViewBase {
    const view = new FuncFlowView();
    file.readAsString().then((json) => view.loadFromJson(json))
      .catch((e) => grok.shell.error(`Cannot open ${file.name}: ${e instanceof Error ? e.message : e}`));
    return view;
  }

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
    meta: {role: 'creationScriptEditor', includeInFlow: 'false'},
  })
  static async openCreationScriptFlowDialog(script: string, tableIds: string[], show: boolean = true): Promise<DG.Dialog> {
    const loaded = await Promise.all((tableIds ?? []).map((id) => grok.dapi.tables.find(id)));
    const tableInfos = loaded.filter((t): t is DG.TableInfo => t != null);
    const view = new FuncFlowView(tableInfos, {outputPanel: false});
    view.name = `Creation Script`;
    view.setMinimapCollapsed(true);
    try {
      await view.loadFromCreationScript(script);
    } catch (e) {
      grok.shell.error(`Failed to load flow from creation script`);
      console.error(e);
    }
    let promoted = false;
    const d = ui.dialog({title: 'Creation Script Flow'})
      .add(view.root)
      .addButton('Open In Editor', () => {
        view.setMinimapCollapsed(false);
        view.enableOutputPanel();
        promoted = true;
        grok.shell.addView(view);
        setTimeout(() => view.fitToScreen(), 100);
        d.close();
      });
    // A dialog-hosted view is never detached by the shell — release it on close.
    d.onClose.subscribe(() => {
      if (!promoted) view.detach();
    });
    if (show)
      d.show({resizable: true, width: 800, height: 600});
    return d;
  }

  /** Consumed by core through the `scriptHandler.editorFunction` seam. */
  @grok.decorators.func({
    name: 'flowScriptEditor',
    description: 'Opens the visual Flow editor for a flow script entity',
    meta: {includeInFlow: 'false'},
  })
  static flowScriptEditor(
    @grok.decorators.param({type: 'script'}) script: DG.Script): DG.ViewBase {
    return FlowEntityHandler.instance.editorView(script);
  }

  @grok.decorators.func({
    name: 'flowScriptPreview',
    meta: {includeInFlow: 'false'},
  })
  static flowScriptPreview(
    @grok.decorators.param({type: 'script'}) script: DG.Script): DG.ViewBase {
    return FlowEntityHandler.instance.previewView(script);
  }

  @grok.decorators.func({
    name: 'flowScriptWidget',
    meta: {includeInFlow: 'false'},
  })
  static flowScriptWidget(
    @grok.decorators.param({type: 'script'}) script: DG.Script): DG.Widget {
    return FlowEntityHandler.instance.widget(script);
  }

  // 'flowViewFunction'-tagged: the AI assistant acts on the open editor through these.

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

  // Flow-native data operations (bodies in ops/data-ops.ts). Declared HERE because the
  // server scans only a fixed file list for annotations — a separate module registers
  // nothing. No friendlyName: it doesn't survive publishing, so names are verbs the
  // humanizer titles well (`filterRows` → `Filter Rows`).

  @grok.decorators.func({
    name: 'filterRows',
    description: 'Keeps the rows matching a condition, as a new table',
    meta: {includeInFlow: 'true'},
  })
  static filterRows(
    @grok.decorators.param({options: {nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({options: {nullable: false, description: 'Boolean expression over the table columns'}}) condition: string,
  ): Promise<DG.DataFrame> {
    return dataOps.filterRows(table, condition);
  }

  @grok.decorators.func({
    name: 'deleteRows',
    description: 'Removes the rows matching a condition, as a new table',
    meta: {includeInFlow: 'true'},
  })
  static deleteRows(
    @grok.decorators.param({options: {nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({options: {nullable: false, description: 'Boolean expression selecting the rows to drop'}}) condition: string,
  ): Promise<DG.DataFrame> {
    return dataOps.deleteRows(table, condition);
  }

  @grok.decorators.func({
    name: 'extractRows',
    description: 'Rows matching a condition, keeping only the chosen columns',
    meta: {includeInFlow: 'true'},
  })
  static extractRows(
    @grok.decorators.param({options: {nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({options: {nullable: false, description: 'Boolean expression over the table columns'}}) condition: string,
    @grok.decorators.param({type: 'column_list', options: {nullable: true, description: 'Columns to keep. Leave empty to keep all of them'}}) columns?: DG.Column[],
  ): Promise<DG.DataFrame> {
    return dataOps.extractRows(table, condition, columns);
  }

  @grok.decorators.func({
    name: 'selectRows',
    description: 'Selects the rows matching a condition and passes the table on',
    meta: {includeInFlow: 'true'},
  })
  static selectRows(
    @grok.decorators.param({options: {nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({options: {nullable: false, description: 'Boolean expression over the table columns'}}) condition: string,
    @grok.decorators.param({options: {initialValue: 'true', description: 'Drop the previous selection instead of adding to it'}}) clearSelection: boolean,
  ): Promise<DG.DataFrame> {
    return dataOps.selectRows(table, condition, clearSelection);
  }

  @grok.decorators.func({
    name: 'filterRandomRows',
    description: 'A reproducible random sample of rows, as a new table',
    meta: {includeInFlow: 'true'},
  })
  static filterRandomRows(
    @grok.decorators.param({options: {nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'int', options: {nullable: false, min: '1', description: 'How many rows to keep'}}) count: number,
    @grok.decorators.param({type: 'int', options: {initialValue: '42', description: 'Random seed. The same seed always draws the same rows'}}) seed: number,
  ): DG.DataFrame {
    return dataOps.filterRandomRows(table, count, seed);
  }

  @grok.decorators.func({
    name: 'selectRandomRows',
    description: 'Selects a reproducible random sample of rows and passes the table on',
    meta: {includeInFlow: 'true'},
  })
  static selectRandomRows(
    @grok.decorators.param({options: {nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'int', options: {nullable: false, min: '1', description: 'How many rows to select'}}) count: number,
    @grok.decorators.param({type: 'int', options: {initialValue: '42', description: 'Random seed. The same seed always draws the same rows'}}) seed: number,
    @grok.decorators.param({options: {initialValue: 'true', description: 'Drop the previous selection instead of adding to it'}}) clearSelection: boolean,
  ): DG.DataFrame {
    return dataOps.selectRandomRows(table, count, seed, clearSelection);
  }

  @grok.decorators.func({
    name: 'deleteColumns',
    description: 'A copy of the table without the chosen columns. Removes selected columns',
    meta: {includeInFlow: 'true'},
  })
  static deleteColumns(
    @grok.decorators.param({options: {nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column_list', options: {nullable: false, description: 'Columns to remove'}}) columns: DG.Column[],
  ): DG.DataFrame {
    return dataOps.deleteColumns(table, columns);
  }

  @grok.decorators.func({
    name: 'tagColumns',
    description: 'Sets a tag on the chosen columns and passes the table on',
    meta: {includeInFlow: 'true'},
  })
  static tagColumns(
    @grok.decorators.param({options: {nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column_list', options: {nullable: false, description: 'Columns to tag'}}) columns: DG.Column[],
    @grok.decorators.param({options: {nullable: false, description: 'Tag name, for example units or .formula'}}) tag: string,
    @grok.decorators.param({options: {nullable: true, description: 'Tag value'}}) value: string,
  ): DG.DataFrame {
    return dataOps.tagColumns(table, columns, tag, value);
  }

  @grok.decorators.func({
    name: 'expressionToColumn',
    description: 'Computes an expression into a new column of the table',
    meta: {includeInFlow: 'true'},
  })
  static expressionToColumn(
    @grok.decorators.param({options: {nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({options: {nullable: false, description: 'Formula over the table columns'}}) expression: string,
    @grok.decorators.param({options: {nullable: false, description: 'Name of the resulting column'}}) name: string,
    @grok.decorators.param({options: {choices: ['auto', 'string', 'int', 'double', 'bool', 'datetime', 'qnum'], initialValue: 'auto', description: 'Column type. auto infers it from the expression'}}) type: string,
  ): Promise<DG.Column> {
    return dataOps.expressionToColumn(table, expression, name, type);
  }

  @grok.decorators.func({
    name: 'aggregate',
    description: 'Groups rows and aggregates columns. Add a pivot column to build a pivot table',
    meta: {includeInFlow: 'true'},
  })
  static aggregate(
    @grok.decorators.param({options: {nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column_list', options: {nullable: true, description: 'Columns to group by. Leave empty to aggregate the whole table into one row'}}) groupByColumns: DG.Column[],
    @grok.decorators.param({options: {nullable: false, description: 'Aggregations to compute, as a list of column and function pairs'}}) aggregations: string,
    @grok.decorators.param({type: 'column_list', options: {nullable: true, description: 'Columns whose values become result columns'}}) pivotColumns?: DG.Column[],
  ): DG.DataFrame {
    return dataOps.aggregate(table, groupByColumns, aggregations, pivotColumns);
  }

  @grok.decorators.func({
    name: 'unpivot',
    description: 'Wide to long. Each merged column becomes a category and value row pair',
    meta: {includeInFlow: 'true'},
  })
  static unpivot(
    @grok.decorators.param({options: {nullable: false}}) table: DG.DataFrame,
    @grok.decorators.param({type: 'column_list', options: {nullable: true, description: 'Columns repeated alongside every produced row'}}) copyColumns: DG.Column[],
    @grok.decorators.param({type: 'column_list', options: {nullable: false, description: 'Columns folded into the category and value pair'}}) mergeColumns: DG.Column[],
    @grok.decorators.param({options: {initialValue: 'Category', description: 'Name of the column holding the source column names'}}) categoryColumnName: string,
    @grok.decorators.param({options: {initialValue: 'Value', description: 'Name of the column holding the values'}}) valueColumnName: string,
  ): Promise<DG.DataFrame> {
    return dataOps.unpivot(table, copyColumns, mergeColumns, categoryColumnName, valueColumnName);
  }

  @grok.decorators.fileViewer({fileViewer: 'flow'})
  static viewFlowFile(file: DG.FileInfo): DG.ViewBase {
    const view = new FuncFlowView();
    file.readAsString()
      .then((text) => view.loadFromDoc(parseFlowBody(text).doc))
      .catch((e) => grok.shell.error(`Cannot open flow file: ${e?.message ?? e}`));
    return view;
  }
}
