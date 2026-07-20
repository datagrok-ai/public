import {PackageFunctions} from './package';
import * as DG from 'datagrok-api/dg';

//name: Flow
//description: Interactive function chain designer
//tags: app
//input: string path { meta.url: true; optional: true }
//output: view result
//meta.role: app
export function funcflowApp(path?: string) : any {
  return PackageFunctions.funcflowApp(path);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: ffjson
export function viewFuncFlow(file: DG.FileInfo) : any {
  return PackageFunctions.viewFuncFlow(file);
}

//description: Builds a flow diagram from a table creation script and opens it in the Flow editor
//input: string script 
//output: view result
export async function flowFromCreationScript(script: string) : Promise<any> {
  return await PackageFunctions.flowFromCreationScript(script);
}

//input: string script 
//input: list<string> tableIds 
//input: bool show 
//output: dynamic result
//meta.role: creationScriptEditor
//meta.includeInFlow: false
export async function openCreationScriptFlowDialog(script: string, tableIds: string[], show: boolean) : Promise<any> {
  return await PackageFunctions.openCreationScriptFlowDialog(script, tableIds, show);
}

//name: testDialog
export function testDialog() : void {
  PackageFunctions.testDialog();
}

//description: Opens the visual Flow editor for a flow script entity
//input: script script 
//output: view result
//meta.includeInFlow: false
export function flowScriptEditor(script: any) : any {
  return PackageFunctions.flowScriptEditor(script);
}

//input: script script 
//output: view result
//meta.includeInFlow: false
export function flowScriptPreview(script: any) : any {
  return PackageFunctions.flowScriptPreview(script);
}

//input: script script 
//output: widget result
//meta.includeInFlow: false
export function flowScriptWidget(script: any) : any {
  return PackageFunctions.flowScriptWidget(script);
}

//description: List the current flow graph: all nodes (id, label, type, status, set input values) and connections. Call this first to understand what is on the canvas
//tags: flowViewFunction
//input: view view 
//output: dynamic result
//meta.includeInFlow: false
export function listFlowNodes(view: any) : any {
  return PackageFunctions.listFlowNodes(view);
}

//description: Ports (with DG types), editable input values, unmet requirements, and last-run outputs of one node
//tags: flowViewFunction
//input: view view 
//input: string nodeId 
//output: dynamic result
//meta.includeInFlow: false
export function getFlowNodeDetails(view: any, nodeId: string) : any {
  return PackageFunctions.getFlowNodeDetails(view, nodeId);
}

//description: Search the flow node catalog (a curated subset of platform functions plus input/output/utility nodes). ALWAYS filter: pass a query with what the node should do (e.g. "join tables", "open file"), and/or a DG type it must accept or produce (dataframe, column, string, ...). Returns at most limit (default 15) matches with their input/output types
//tags: flowViewFunction
//input: view view 
//input: string query { optional: true; description: Words describing what the node does }
//input: string acceptsInputType { optional: true; description: DG type one of its inputs must accept }
//input: string producesOutputType { optional: true; description: DG type one of its outputs must produce }
//input: int limit { optional: true }
//output: dynamic result
//meta.includeInFlow: false
export function findFlowNodeTypes(view: any, query?: string, acceptsInputType?: string, producesOutputType?: string, limit?: number) : any {
  return PackageFunctions.findFlowNodeTypes(view, query, acceptsInputType, producesOutputType, limit);
}

//description: Add a node to the canvas by its registered typeName (from findFlowNodeTypes). Optionally set editable input values right away. Returns the new node id and its ports
//tags: flowViewFunction
//input: view view 
//input: string typeName 
//input: string label { optional: true; description: Optional custom title }
//input: map inputValues { optional: true; description: Editable primitive inputs, key to value }
//output: dynamic result
//meta.includeInFlow: false
export async function addFlowNode(view: any, typeName: string, label?: string, inputValues?: any) : Promise<any> {
  return await PackageFunctions.addFlowNode(view, typeName, label, inputValues);
}

//description: Connect a source node output to a target node input (port keys from getFlowNodeDetails / addFlowNode). Types must be compatible
//tags: flowViewFunction
//input: view view 
//input: string sourceNodeId 
//input: string sourceOutput 
//input: string targetNodeId 
//input: string targetInput 
//output: dynamic result
//meta.includeInFlow: false
export async function connectFlowNodes(view: any, sourceNodeId: string, sourceOutput: string, targetNodeId: string, targetInput: string) : Promise<any> {
  return await PackageFunctions.connectFlowNodes(view, sourceNodeId, sourceOutput, targetNodeId, targetInput);
}

//description: Set editable input values of a node (key to value; keys from getFlowNodeDetails). Marks the node and its downstream stale
//tags: flowViewFunction
//input: view view 
//input: string nodeId 
//input: map values 
//output: dynamic result
//meta.includeInFlow: false
export async function setFlowNodeInputs(view: any, nodeId: string, values: any) : Promise<any> {
  return await PackageFunctions.setFlowNodeInputs(view, nodeId, values);
}

//description: Select a node on the canvas so the user sees it (opens its properties panel)
//tags: flowViewFunction
//input: view view 
//input: string nodeId 
//output: dynamic result
//meta.includeInFlow: false
export async function selectFlowNode(view: any, nodeId: string) : Promise<any> {
  return await PackageFunctions.selectFlowNode(view, nodeId);
}

//description: List Flow built-in interactive guides: step-by-step tutorials and short "how do I" walkthroughs that highlight the actual UI. When the user asks how to do something in Flow, check here first — a matching guide beats a textual explanation
//tags: flowViewFunction
//input: view view 
//input: string query { optional: true; description: Words to filter by }
//output: dynamic result
//meta.includeInFlow: false
export function listFlowGuides(view: any, query?: string) : any {
  return PackageFunctions.listFlowGuides(view, query);
}

//description: Start an interactive guide (id from listFlowGuides) — it highlights the real UI step by step and waits for the user to act. ALWAYS confirm with the user first before starting one; never launch it unasked
//tags: flowViewFunction
//input: view view 
//input: string guideId 
//output: dynamic result
//meta.includeInFlow: false
export function startFlowGuide(view: any, guideId: string) : any {
  return PackageFunctions.startFlowGuide(view, guideId);
}

//description: Validate and execute the whole flow. Returns validation problems instead of running if the graph is invalid; otherwise waits for the run and reports per-node failures
//tags: flowViewFunction
//input: view view 
//output: dynamic result
//meta.includeInFlow: false
export async function runFlow(view: any) : Promise<any> {
  return await PackageFunctions.runFlow(view);
}

//input: file file 
//output: view result
//meta.role: fileViewer
//meta.fileViewer: flow
export function viewFlowFile(file: DG.FileInfo) : any {
  return PackageFunctions.viewFlowFile(file);
}
