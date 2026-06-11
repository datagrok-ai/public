/** Builds a flow graph from a Datagrok table-creation script.
 *
 * Creation scripts are the linear cascades of function calls Datagrok records
 * when a table is produced in a reproducible way (file share, query, …) and
 * then transformed (calculated columns, chem properties, …):
 *
 *   Mol1K = OpenFile("System:AppData/Chem/mol1K.csv") //{"timestamp": …}
 *   Chem:addChemPropertiesColumns(Mol1K, "molecule", true, …)
 *   AddNewColumn(Mol1K, "${HBA}+${HBD}+${LogP}", "sumOfSome", subscribeOnChanges = true)
 *
 * Each line parses (`grok.functions.parse`) into a `DG.FuncCall`. Confirmed
 * shapes (from the live parser):
 *   - An assignment is a `SetVar` call: `variableName:string` + `value:dynamic`
 *     (the assigned call/literal).
 *   - A variable read is a `GetVar` call (`variableName:string`).
 *   - A column argument parses to `ResolveColumn(value:dynamic, parentTable:dataframe)`
 *     where `value` is the column name string and `parentTable` is often **null**
 *     (resolved at runtime from the enclosing call's table). The `ResolveColumn`
 *     function misbehaves on the Datagrok side, so we substitute the built-in
 *     **Select Column** utility (`table.col('name')`): its `columnName` property
 *     gets the column name and its `table` input is wired to the enclosing call's
 *     table (or `parentTable` when the script provides one). `ResolveColumnList`
 *     maps to **Select Columns** the same way.
 *
 * Build rules:
 *   1. Each call becomes a `FuncNode`. Primitive inputs on *editable* slots
 *      (string/int/double/num/bool) are stored in `node.inputValues` (shown and
 *      editable in the property panel). Primitive inputs on other slots
 *      (e.g. `dynamic` column names) get a **Constant node** wired in — those
 *      slots aren't editable in the panel and need an explicit producer.
 *   2. `FuncCall` inputs are resolved recursively; `GetVar` resolves through the
 *      variable table; a resolver's missing dataframe input is wired to the
 *      enclosing call's table.
 *   3. **Bare (mutating) calls** — `addChemPropertiesColumns(Mol1K, …)` mutate
 *      their table in place, so after the call the consumed variable advances to
 *      that node's `<input>__pt` pass-through. The next consumer connects there,
 *      so the topological sort reproduces the script's line order while the
 *      compiler resolves the pass-through to the same expression (no extra
 *      variable). **Assignments** (`y = f(x)`) produce a new variable and never
 *      advance their inputs.
 *   4. The first assigned variable is the script's result, wired to an output
 *      node.
 *
 * The graph is built in memory (DOM-free, synchronous — see
 * `buildCreationScriptGraph`) so it is directly testable; `applyGraphToEditor`
 * pushes it into a live `FlowEditor`. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {TypedSocket} from '../rete/sockets';
import {constLabel} from '../rete/nodes/utility-nodes';
import {
  registerBuiltinNodes, registerAllFunctions, createNode, ensureFuncNodeType,
} from '../rete/node-factory';

const PASSTHROUGH_SUFFIX = '__pt';

/** DG slot types whose primitive value is editable in the property panel and
 *  emitted as a literal by the compiler — kept as `inputValues` rather than a
 *  wired-in Constant node (see `func-node.ts` PRIMITIVE_DEFAULTS). */
const EDITABLE_PRIMITIVE_SLOTS = new Set(['string', 'int', 'double', 'num', 'bool']);

/** Platform resolvers replaced by the Select Column / Select Columns /
 *  Select Table utilities (the `ResolveColumn` / `ResolveTable` functions
 *  misbehave at runtime). */
const SINGLE_COLUMN_RESOLVER = 'resolvecolumn';
const LIST_COLUMN_RESOLVER = 'resolvecolumnlist';
const TABLE_RESOLVER = 'resolvetable';
const SELECT_COLUMN_TYPE = 'Utilities/Select Column';
const SELECT_COLUMNS_TYPE = 'Utilities/Select Columns';
const SELECT_TABLE_TYPE = 'Utilities/Select Table';

/** A value's location in the graph: a specific output slot of a node. */
interface OutputRef {
  node: FlowNode;
  outputKey: string;
}

interface BuiltConnection {
  source: FlowNode;
  sourceKey: string;
  target: FlowNode;
  targetKey: string;
}

export interface BuiltGraph {
  nodes: FlowNode[];
  connections: BuiltConnection[];
  /** Every script variable, in first-assignment order — each is wired to its
   *  own output node. */
  outputVariables: string[];
  warnings: string[];
}

export interface ImportResult {
  nodesAdded: number;
  connectionsAdded: number;
  outputVariables: string[];
  warnings: string[];
}

/** Build the in-memory graph for a creation script. Synchronous and DOM-free:
 *  it constructs `FlowNode` instances and connection records but touches no
 *  editor, so it is the unit-test entry point. */
export function buildCreationScriptGraph(script: string): BuiltGraph {
  registerBuiltinNodes();
  registerAllFunctions();
  return new CreationScriptBuilder().build(script);
}

/** Push a built graph into a live editor (clears nothing — caller decides). */
export async function applyGraphToEditor(graph: BuiltGraph, flow: FlowEditor): Promise<number> {
  for (const node of graph.nodes)
    await flow.addNodeAt(node, node.pos.x, node.pos.y);
  let connectionsAdded = 0;
  for (const c of graph.connections) {
    if (await flow.addConnectionByKeys(c.source.id, c.sourceKey, c.target.id, c.targetKey))
      connectionsAdded++;
  }
  return connectionsAdded;
}

/** Build the graph and add it to the editor. */
export async function buildFlowFromCreationScript(flow: FlowEditor, script: string): Promise<ImportResult> {
  const graph = buildCreationScriptGraph(script);
  const connectionsAdded = await applyGraphToEditor(graph, flow);
  return {
    nodesAdded: graph.nodes.length,
    connectionsAdded,
    outputVariables: graph.outputVariables,
    warnings: graph.warnings,
  };
}

/** What kind of resolution an input value produced. */
type Resolution =
  | {kind: 'ref'; ref: OutputRef; varKey: string | null}
  | {kind: 'literal'; value: string | number | boolean}
  | {kind: 'skip'};

interface CallContext {
  /** True for bare top-level statements (in-place mutators): consumed variables
   *  advance to the node's pass-through. False for assignment RHS and nested
   *  calls. */
  advanceConsumedVars: boolean;
  /** Table the enclosing call operates on, used to wire a resolver's missing
   *  `parentTable`. */
  contextTable: OutputRef | null;
}

class CreationScriptBuilder {
  /** Lowercased variable name → output slot currently holding its value. */
  private readonly variables = new Map<string, OutputRef>();
  /** Variable names (original casing) in first-assignment order. */
  private readonly variableOrder: string[] = [];
  private readonly nodes: FlowNode[] = [];
  private readonly connections: BuiltConnection[] = [];
  private readonly layer = new Map<FlowNode, number>();
  private readonly warnings: string[] = [];

  build(script: string): BuiltGraph {
    const calls = this.parseLines(script);
    if (calls.length === 0)
      throw new Error('No function calls could be parsed from the creation script');

    for (const fc of calls) {
      const assignment = this.asAssignment(fc);
      if (assignment)
        this.processAssignment(assignment.name, assignment.value);
      else
        this.addCall(fc, {advanceConsumedVars: true, contextTable: null});
    }

    // Every variable becomes a script result: each gets its own output node,
    // connected to the variable's final ref (the last pass-through in its
    // chain, so the output carries the fully transformed value).
    const outputVariables: string[] = [];
    const usedParamNames = new Set<string>();
    for (const name of this.variableOrder) {
      if (!this.variables.has(name.toLowerCase())) continue;
      this.wireOutput(name, usedParamNames);
      outputVariables.push(name);
    }
    if (outputVariables.length === 0)
      this.warnings.push('No variable assignments found — the flow has no output node');

    this.layout();
    return {
      nodes: this.nodes,
      connections: this.connections,
      outputVariables,
      warnings: this.warnings,
    };
  }

  // ---------- parsing ----------

  private parseLines(script: string): DG.FuncCall[] {
    const calls: DG.FuncCall[] = [];
    for (const rawLine of script.split('\n')) {
      const line = rawLine.trim();
      if (line === '' || line.startsWith('//')) continue;
      const fc = this.parseLine(line) ?? this.parseLine(stripTrailingComment(line));
      if (fc) calls.push(fc);
      else this.warnings.push(`Skipped unparsable line: ${ellipsis(line, 80)}`);
    }
    return calls;
  }

  private parseLine(line: string): DG.FuncCall | null {
    if (line === '') return null;
    try {
      const parsed = grok.functions.parse(line, false);
      return parsed instanceof DG.FuncCall ? parsed : null;
    } catch {
      return null;
    }
  }

  // ---------- call classification ----------

  /** A variable assignment parses into a SetVar-shaped call: exactly two
   *  inputs, a non-empty string `variableName` and a non-null `value`. */
  private asAssignment(fc: DG.FuncCall): {name: string; value: unknown} | null {
    try {
      if ((fc.func?.name?.toLowerCase() ?? '') !== 'setvar') return null;
      if (Object.entries(fc.inputs).length !== 2) return null;
      const name: unknown = fc.inputs['variableName'];
      const value: unknown = fc.inputs['value'];
      if (typeof name !== 'string' || name === '' || value === null || value === undefined) return null;
      return {name, value};
    } catch {
      return null;
    }
  }

  /** For GetVar calls, the referenced variable name; null for anything else. */
  private variableNameOf(fc: DG.FuncCall): string | null {
    try {
      if ((fc.func?.name?.toLowerCase() ?? '') !== 'getvar') return null;
      const name: unknown = fc.inputs['variableName'] ?? fc.inputs['name'];
      return typeof name === 'string' && name !== '' ? name : null;
    } catch {
      return null;
    }
  }

  /** Whether the call is a column resolver we substitute with Select Column(s).
   *  Returns the matching utility node type, or null. */
  private columnResolverNodeType(fc: DG.FuncCall): string | null {
    const name = fc.func?.name?.toLowerCase() ?? '';
    if (name === SINGLE_COLUMN_RESOLVER) return SELECT_COLUMN_TYPE;
    if (name === LIST_COLUMN_RESOLVER) return SELECT_COLUMNS_TYPE;
    return null;
  }

  private isSingleColumnResolver(v: unknown): v is DG.FuncCall {
    if (!(v instanceof DG.FuncCall)) return false;
    try {
      return (v.func?.name?.toLowerCase() ?? '') === SINGLE_COLUMN_RESOLVER;
    } catch {
      return false;
    }
  }

  // ---------- graph building ----------

  private processAssignment(name: string, value: unknown): void {
    // Assignment RHS: resolve the producing value, do NOT advance inputs.
    const res = this.resolveValue(value, {advanceConsumedVars: false, contextTable: null}, `variable "${name}"`);
    if (res.kind !== 'ref') {
      if (res.kind === 'literal') {
        const ref = this.addConstant(res.value);
        if (ref) this.registerVariable(name, ref);
      }
      return;
    }
    this.registerVariable(name, res.ref);
  }

  private registerVariable(name: string, ref: OutputRef): void {
    if (!this.variables.has(name.toLowerCase())) this.variableOrder.push(name);
    this.variables.set(name.toLowerCase(), ref);
  }

  /** Resolve an input value to a graph location or a literal. */
  private resolveValue(value: unknown, ctx: CallContext, context: string): Resolution {
    if (value instanceof DG.FuncCall) {
      const varName = this.variableNameOf(value);
      if (varName !== null) {
        const ref = this.variables.get(varName.toLowerCase());
        if (!ref) {
          this.warnings.push(`Unknown variable "${varName}" used by ${context}`);
          return {kind: 'skip'};
        }
        return {kind: 'ref', ref, varKey: varName.toLowerCase()};
      }
      const assignment = this.asAssignment(value);
      if (assignment) {
        this.processAssignment(assignment.name, assignment.value);
        const ref = this.variables.get(assignment.name.toLowerCase());
        return ref ? {kind: 'ref', ref, varKey: null} : {kind: 'skip'};
      }
      const colResolverType = this.columnResolverNodeType(value);
      if (colResolverType !== null) {
        const ref = this.addColumnSelection(value, colResolverType, ctx);
        return ref ? {kind: 'ref', ref, varKey: null} : {kind: 'skip'};
      }
      if ((value.func?.name?.toLowerCase() ?? '') === TABLE_RESOLVER) {
        const ref = this.addTableSelection(value);
        return ref ? {kind: 'ref', ref, varKey: null} : {kind: 'skip'};
      }
      const ref = this.addCall(value, {advanceConsumedVars: false, contextTable: ctx.contextTable});
      return ref ? {kind: 'ref', ref, varKey: null} : {kind: 'skip'};
    }
    if (Array.isArray(value)) {
      // A column_list argument parses to an array of ResolveColumn calls —
      // map the whole array to a single Select Columns utility.
      if (value.length > 0 && value.every((v) => this.isSingleColumnResolver(v))) {
        const items = value as DG.FuncCall[];
        const names = items.flatMap((item) => this.columnNames(item));
        const tableRef = this.explicitParentTable(items, ctx) ?? ctx.contextTable;
        const ref = this.addSelect(SELECT_COLUMNS_TYPE, names, tableRef);
        return ref ? {kind: 'ref', ref, varKey: null} : {kind: 'skip'};
      }
      if (value.every(isPrimitive)) {
        const ref = this.addListConstant(value as Array<string | number | boolean>);
        return ref ? {kind: 'ref', ref, varKey: null} : {kind: 'skip'};
      }
      this.warnings.push(`Unsupported list value for ${context}`);
      return {kind: 'skip'};
    }
    if (isPrimitive(value)) return {kind: 'literal', value};
    if (value === null || value === undefined) return {kind: 'skip'};
    this.warnings.push(`Unsupported value for ${context}: ${typeof value}`);
    return {kind: 'skip'};
  }

  /** Create a FuncNode for the call, wire inputs (constants for non-editable
   *  primitive slots, connections for FuncCall inputs), infer a resolver's
   *  missing table, and advance consumed variables for mutating calls. */
  private addCall(fc: DG.FuncCall, ctx: CallContext): OutputRef | null {
    const func = fc.func;
    const node = createNode(ensureFuncNodeType(func));
    if (!node) {
      this.warnings.push(`Could not create a node for function "${func?.name}"`);
      return null;
    }
    this.addNode(node);

    const params = func.inputs;
    // Resolve dataframe inputs first to establish this call's table context for
    // sibling column/resolver inputs.
    const ordered = [...params].sort((a, b) =>
      (this.isDataframeParam(b) ? 1 : 0) - (this.isDataframeParam(a) ? 1 : 0));

    let ownTable: OutputRef | null = null;
    /** Dataframe param name → its resolved source, for sibling pairing. */
    const dfSources = new Map<string, OutputRef>();
    /** Table context for a sibling param. Numbered params pair with the
     *  dataframe param sharing the numeric suffix (JoinTables: keys2/values2
     *  resolve against table2, not the first table); otherwise the call's
     *  first table, else the outer context. */
    const tableCtxFor = (paramName: string): OutputRef | null => {
      const suffix = /(\d+)$/.exec(paramName)?.[1];
      if (suffix !== undefined) {
        for (const [dfName, ref] of dfSources)
          if (dfName.endsWith(suffix)) return ref;
      }
      return ownTable ?? ctx.contextTable;
    };
    let maxSourceLayer = -1;

    for (const param of ordered) {
      if (!(param.name in node.inputs)) continue;
      const value = this.safeInput(fc, param.name);
      const slotType = this.slotType(node, param.name);
      const res = this.resolveValue(
        value, {advanceConsumedVars: false, contextTable: tableCtxFor(param.name)},
        `input "${param.name}" of ${func.name}`);

      if (res.kind === 'literal') {
        if (!EDITABLE_PRIMITIVE_SLOTS.has(slotType)) {
          const constRef = this.addConstant(res.value);
          if (constRef) {
            this.connect(constRef, node, param.name);
            maxSourceLayer = Math.max(maxSourceLayer, this.layer.get(constRef.node) ?? 0);
          }
        } else
          node.inputValues[param.name] = res.value;
        continue;
      }
      if (res.kind !== 'ref') continue;

      this.connect(res.ref, node, param.name);
      maxSourceLayer = Math.max(maxSourceLayer, this.layer.get(res.ref.node) ?? 0);
      if (this.isDataframeParam(param)) {
        dfSources.set(param.name, res.ref);
        if (!ownTable) ownTable = res.ref;
      }
      if (ctx.advanceConsumedVars && res.varKey !== null) {
        const ptKey = param.name + PASSTHROUGH_SUFFIX;
        if (ptKey in node.outputs) this.variables.set(res.varKey, {node, outputKey: ptKey});
      }
    }

    this.layer.set(node, maxSourceLayer + 1);
    return this.primaryOutput(node);
  }

  /** Substitute a `ResolveColumn` / `ResolveColumnList` call with a Select
   *  Column / Select Columns utility (`table.col('name')`) — the platform
   *  resolver misbehaves at runtime. The column name(s) come from the
   *  resolver's `value` input; the `table` input is the explicit `parentTable`
   *  if present, else the enclosing call's table. */
  private addColumnSelection(fc: DG.FuncCall, typeName: string, ctx: CallContext): OutputRef | null {
    const names = this.columnNames(fc);
    const tableRef = this.explicitParentTable([fc], ctx) ?? ctx.contextTable;
    return this.addSelect(typeName, names, tableRef);
  }

  /** Substitute a `ResolveTable` call with the Select Table utility
   *  (`grok.shell.tableByName(name)`) — the platform resolver misbehaves at
   *  runtime. Titled `table: <name>` so collapsed nodes stay readable. */
  private addTableSelection(fc: DG.FuncCall): OutputRef | null {
    const node = createNode(SELECT_TABLE_TYPE);
    if (!node) return null;
    const value = this.safeInput(fc, 'value');
    const name = typeof value === 'string' ? value : String(value ?? '');
    node.properties['tableName'] = name;
    if (name !== '') node.label = `table: ${name}`;
    this.addNode(node);
    this.layer.set(node, 0);
    return {node, outputKey: 'table'};
  }

  /** Create a Select Column / Select Columns utility for the given column
   *  names, wired to `tableRef`. */
  private addSelect(typeName: string, names: string[], tableRef: OutputRef | null): OutputRef | null {
    const node = createNode(typeName);
    if (!node) return null;
    this.addNode(node);

    if (typeName === SELECT_COLUMNS_TYPE)
      node.properties['columnNames'] = names.join(', ');
    else
      node.properties['columnName'] = names[0] ?? '';

    let layer = 0;
    if (tableRef && 'table' in node.inputs) {
      this.connect(tableRef, node, 'table');
      layer = (this.layer.get(tableRef.node) ?? 0) + 1;
    } else if (!tableRef)
      this.warnings.push(`Could not resolve the table for column "${names.join(', ')}"`);
    this.layer.set(node, layer);
    return this.primaryOutput(node);
  }

  /** The first explicit (non-null) `parentTable` among the resolver calls,
   *  resolved to its graph location. */
  private explicitParentTable(items: DG.FuncCall[], ctx: CallContext): OutputRef | null {
    for (const item of items) {
      const parentTable = this.safeInput(item, 'parentTable');
      if (parentTable === null || parentTable === undefined) continue;
      const res = this.resolveValue(
        parentTable, {advanceConsumedVars: false, contextTable: ctx.contextTable},
        'parentTable of a column reference');
      if (res.kind === 'ref') return res.ref;
    }
    return null;
  }

  /** Column name(s) from a resolver's `value` input (string or list). */
  private columnNames(fc: DG.FuncCall): string[] {
    const value = this.safeInput(fc, 'value');
    if (typeof value === 'string') return [value];
    if (Array.isArray(value)) return value.map((v) => String(v));
    if (value != null) {
      this.warnings.push(`Unexpected column reference of type ${typeof value}`);
      return [String(value)];
    }
    return [];
  }

  private primaryOutput(node: FlowNode): OutputRef | null {
    const keys = Object.keys(node.outputs);
    const key = keys.find((k) => !k.endsWith(PASSTHROUGH_SUFFIX)) ?? keys[0];
    return key ? {node, outputKey: key} : null;
  }

  /** Literal value → constant utility node, titled after its value. */
  private addConstant(value: string | number | boolean): OutputRef | null {
    let typeName: string;
    if (typeof value === 'boolean')
      typeName = 'Constants/Boolean';
    else if (typeof value === 'number')
      typeName = Number.isInteger(value) ? 'Constants/Int' : 'Constants/Double';
    else
      typeName = 'Constants/String';
    const node = createNode(typeName);
    if (!node) return null;
    node.properties['value'] = value;
    node.label = constLabel(typeName.split('/')[1], value);
    // ConstStringNode renders an inline widget — keep it in sync.
    const control = node.controls['value'] as {value?: string} | undefined;
    if (control) control.value = String(value);
    this.addNode(node);
    this.layer.set(node, 0);
    return {node, outputKey: 'value'};
  }

  /** Array of primitive literals → List constant (comma-separated). */
  private addListConstant(values: Array<string | number | boolean>): OutputRef | null {
    const node = createNode('Constants/List');
    if (!node) return null;
    const joined = values.map((v) => String(v)).join(', ');
    node.properties['value'] = joined;
    node.label = constLabel('List', joined);
    this.addNode(node);
    this.layer.set(node, 0);
    return {node, outputKey: 'value'};
  }

  /** Wire one script variable to its own Table/Value output node. The
   *  validator requires unique output param names, so sanitized names are
   *  deduplicated across outputs. */
  private wireOutput(varName: string, usedParamNames: Set<string>): void {
    const ref = this.variables.get(varName.toLowerCase());
    if (!ref) return;
    const dgType = this.outputType(ref);
    const isTable = dgType === 'dataframe';
    const node = createNode(isTable ? 'Outputs/Table Output' : 'Outputs/Value Output');
    if (!node) return;
    let paramName = toParamName(varName);
    for (let i = 2; usedParamNames.has(paramName.toLowerCase()); i++)
      paramName = `${toParamName(varName)}_${i}`;
    usedParamNames.add(paramName.toLowerCase());
    node.properties['paramName'] = paramName;
    node.label = isTable ? `Table Output: ${paramName}` : `Output: ${paramName}`;
    if (!isTable && dgType && dgType !== 'dynamic' && dgType !== 'object')
      node.properties['outputType'] = dgType;
    this.addNode(node);
    const inputKey = Object.keys(node.inputs)[0];
    this.connect(ref, node, inputKey);
    this.layer.set(node, (this.layer.get(ref.node) ?? 0) + 1);
  }

  // ---------- low-level graph ops ----------

  private addNode(node: FlowNode): void {
    // Imported nodes start collapsed (title bar only) so the whole flow fits
    // a view; the user expands individual nodes as needed. Connections keep
    // their endpoints — collapsed nodes still expose socket DOM at the title
    // bar edges (see node-component.tsx).
    node.collapsed = true;
    this.nodes.push(node);
  }

  private connect(source: OutputRef, target: FlowNode, targetKey: string): void {
    this.connections.push({source: source.node, sourceKey: source.outputKey, target, targetKey});
  }

  private safeInput(fc: DG.FuncCall, name: string): unknown {
    try {
      return fc.inputs[name];
    } catch {
      return undefined;
    }
  }

  private slotType(node: FlowNode, inputKey: string): string {
    const slot = node.inputs[inputKey] as {socket: TypedSocket} | undefined;
    return slot?.socket.dgType ?? 'dynamic';
  }

  private isDataframeParam(param: DG.Property): boolean {
    return param.propertyType === 'dataframe' || String(param.propertyType) === 'dataframe';
  }

  private outputType(ref: OutputRef): string | null {
    const slot = ref.node.outputs[ref.outputKey] as {socket: TypedSocket} | undefined;
    return slot?.socket.dgType ?? null;
  }

  // ---------- layout ----------

  /** Layered left-to-right layout. Build layers (longest path from the
   *  sources — every connection's source is created before its target) become
   *  columns, so every edge points rightward. Within a column, nodes are
   *  ordered by the average vertical center of their predecessors (barycenter)
   *  and then *greedily aligned* to it: each node lands at its predecessors'
   *  center unless that would overlap the node above. Chains read as straight
   *  horizontal lanes, branches fan out below, nothing overlaps. Column
   *  positions accumulate from the widest (estimated) node in each column —
   *  collapsed titles vary a lot (`const: …`, `table: …`). */
  private layout(): void {
    const marginX = 40;
    const marginY = 40;
    const columnGap = 60;
    const rowGap = 24;

    const byLayer = new Map<number, FlowNode[]>();
    for (const node of this.nodes) {
      const layer = this.layer.get(node) ?? 0;
      let bucket = byLayer.get(layer);
      if (!bucket) byLayer.set(layer, bucket = []);
      bucket.push(node);
    }

    const predecessors = new Map<FlowNode, FlowNode[]>();
    for (const c of this.connections) {
      let list = predecessors.get(c.target);
      if (!list) predecessors.set(c.target, list = []);
      list.push(c.source);
    }

    const centerY = new Map<FlowNode, number>();
    let x = marginX;
    for (const layer of Array.from(byLayer.keys()).sort((a, b) => a - b)) {
      const keyed = byLayer.get(layer)!.map((node, index) => {
        const centers = (predecessors.get(node) ?? [])
          .map((p) => centerY.get(p))
          .filter((y): y is number => y !== undefined);
        const barycenter = centers.length > 0 ?
          centers.reduce((a, b) => a + b, 0) / centers.length :
          Number.POSITIVE_INFINITY;
        return {node, index, barycenter};
      });
      keyed.sort((a, b) => a.barycenter === b.barycenter ? a.index - b.index : a.barycenter - b.barycenter);

      let nextFreeY = marginY;
      let columnWidth = 0;
      for (const {node, barycenter} of keyed) {
        const height = estimateNodeHeight(node);
        const desiredTop = Number.isFinite(barycenter) ? barycenter - height / 2 : nextFreeY;
        const top = Math.max(nextFreeY, desiredTop, marginY);
        node.pos = {x, y: top};
        centerY.set(node, top + height / 2);
        nextFreeY = top + height + rowGap;
        columnWidth = Math.max(columnWidth, estimateNodeWidth(node));
      }
      x += columnWidth + columnGap;
    }
  }
}

// ---------- helpers ----------

function isPrimitive(v: unknown): v is string | number | boolean {
  return typeof v === 'string' || typeof v === 'number' || typeof v === 'boolean';
}

/** Estimated rendered height in canvas units, used for stacking before the
 *  DOM exists. Collapsed nodes render as a bare title bar. Expanded: title
 *  ≈ 28px, description ≈ 22px, each socket row ≈ 20px, body padding 12px
 *  (see funcflow.css: .ff-node-title / .ff-socket-row). Input and output
 *  columns sit side by side, so rows = max of the two. Exported for the
 *  layout-invariant tests. */
export function estimateNodeHeight(node: FlowNode): number {
  if (node.collapsed) return 30;
  const rows = Math.max(Object.keys(node.inputs).length, Object.keys(node.outputs).length, 1);
  return 28 + (node.description ? 22 : 0) + 12 + rows * 20;
}

/** Estimated rendered width. Collapsed nodes are title-driven (CSS min-width
 *  160px, ≈6.5px/char at the 12px title font plus status dot and paddings);
 *  expanded nodes are dominated by their socket-label rows. */
export function estimateNodeWidth(node: FlowNode): number {
  const labelWidth = 44 + String(node.label ?? '').length * 6.5;
  return Math.max(node.collapsed ? 160 : 220, labelWidth);
}

/** The graph validator requires output param names to be JS identifiers. */
function toParamName(name: string): string {
  const cleaned = name.replace(/[^A-Za-z0-9_$]/g, '_');
  return /^[A-Za-z_$]/.test(cleaned) ? cleaned : `_${cleaned}`;
}

/** Remove a trailing `// …` comment (creation scripts carry `//{"timestamp"}`
 *  metadata), ignoring `//` inside string literals such as URLs. */
function stripTrailingComment(line: string): string {
  let quote: string | null = null;
  for (let i = 0; i < line.length - 1; i++) {
    const ch = line[i];
    if (quote !== null) {
      if (ch === '\\') i++;
      else if (ch === quote) quote = null;
    } else if (ch === '"' || ch === '\'' || ch === '`')
      quote = ch;
    else if (ch === '/' && line[i + 1] === '/')
      return line.slice(0, i).trim();
  }
  return line;
}

function ellipsis(s: string, max: number): string {
  return s.length <= max ? s : `${s.slice(0, max - 1)}…`;
}
