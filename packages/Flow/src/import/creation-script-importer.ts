/** Builds a flow graph from a Datagrok table-creation script — the reverse of
 *  the creation-script emitter. */

import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode, EXEC_IN_KEY, EXEC_OUT_KEY} from '../rete/scheme';
import {TypedSocket} from '../rete/sockets';
import {layoutGraph} from '../rete/graph-layout';
import {constLabel} from '../rete/nodes/utility-nodes';

export {estimateNodeWidth, estimateNodeHeight} from '../rete/graph-layout';
import {
  registerBuiltinNodes, registerAllFunctions, createNode, ensureFuncNodeType,
} from '../rete/node-factory';

const PASSTHROUGH_SUFFIX = '__pt';

/** Slot types kept as editable `inputValues` rather than a wired-in Constant node. */
const EDITABLE_PRIMITIVE_SLOTS = new Set(['string', 'int', 'double', 'num', 'bool']);

// Platform resolvers replaced by the Select utilities — the originals misbehave at runtime.
const SINGLE_COLUMN_RESOLVER = 'resolvecolumn';
const LIST_COLUMN_RESOLVER = 'resolvecolumnlist';
const TABLE_RESOLVER = 'resolvetable';
const SELECT_COLUMN_TYPE = 'Utilities/Select Column';
const SELECT_COLUMNS_TYPE = 'Utilities/Select Columns';
const SELECT_TABLE_TYPE = 'Utilities/Select Table';

interface OutputRef {
  node: FlowNode;
  outputKey: string;
}

interface BuiltConnection {
  source: FlowNode;
  sourceKey: string;
  target: FlowNode;
  targetKey: string;
  /** Pure run-order edge (exec-out → exec-in), excluded from the layout. */
  order?: boolean;
}

export interface BuiltGraph {
  nodes: FlowNode[];
  connections: BuiltConnection[];
  /** Script variables in first-assignment order, each wired to a SetVar terminal. */
  outputVariables: string[];
  warnings: string[];
}

export interface ImportResult {
  nodesAdded: number;
  connectionsAdded: number;
  outputVariables: string[];
  warnings: string[];
}

/** Synchronous and DOM-free — touches no editor, so it is the unit-test entry point. */
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

type Resolution =
  | {kind: 'ref'; ref: OutputRef; varKey: string | null}
  | {kind: 'literal'; value: string | number | boolean}
  | {kind: 'skip'};

interface CallContext {
  /** True for bare top-level statements (in-place mutators): consumed variables advance to the pass-through. */
  advanceConsumedVars: boolean;
  /** Table the enclosing call operates on — wires a resolver's missing `parentTable`. */
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
  /** Lowercased variable name → its SetVar node, for inferring order edges. */
  private readonly setVarNodes = new Map<string, FlowNode>();

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

    // Each variable terminates in a SetVar node so running the flow registers
    // the value in the context under its original name.
    const outputVariables: string[] = [];
    for (const name of this.variableOrder) {
      if (!this.variables.has(name.toLowerCase())) continue;
      this.wireSetVar(name);
      outputVariables.push(name);
    }
    if (outputVariables.length === 0)
      this.warnings.push('No variable assignments found in the creation script');

    this.inferOrderEdges();
    this.layout();
    return {
      nodes: this.nodes,
      connections: this.connections,
      outputVariables,
      warnings: this.warnings,
    };
  }

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

  /** Don't gate on the input *count*: SetVarFunc has optional `outputName`/`outputIndex`,
   *  so an assignment surfaces up to four inputs — a strict `=== 2` check drops every assignment. */
  private asAssignment(fc: DG.FuncCall): {name: string; value: unknown} | null {
    try {
      if ((fc.func?.name?.toLowerCase() ?? '') !== 'setvar') return null;
      const name: unknown = fc.inputs['variableName'];
      const value: unknown = fc.inputs['value'];
      if (typeof name !== 'string' || name === '' || value === null || value === undefined) return null;
      return {name, value};
    } catch {
      return null;
    }
  }

  /** For GetVar calls, the referenced variable name; null otherwise. */
  private variableNameOf(fc: DG.FuncCall): string | null {
    try {
      if ((fc.func?.name?.toLowerCase() ?? '') !== 'getvar') return null;
      const name: unknown = fc.inputs['variableName'] ?? fc.inputs['name'];
      return typeof name === 'string' && name !== '' ? name : null;
    } catch {
      return null;
    }
  }

  /** The matching Select utility type for a column-resolver call, or null. */
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

  private processAssignment(name: string, value: unknown): void {
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
      // A column_list arg parses to an array of ResolveColumn calls — map it to one Select Columns.
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

  private addCall(fc: DG.FuncCall, ctx: CallContext): OutputRef | null {
    const func = fc.func;
    const node = createNode(ensureFuncNodeType(func));
    if (!node) {
      this.warnings.push(`Could not create a node for function "${func?.name}"`);
      return null;
    }
    this.addNode(node);

    const params = func.inputs;
    // Dataframe inputs first — they establish the table context for sibling column inputs.
    const ordered = [...params].sort((a, b) =>
      (this.isDataframeParam(b) ? 1 : 0) - (this.isDataframeParam(a) ? 1 : 0));

    let ownTable: OutputRef | null = null;
    const dfSources = new Map<string, OutputRef>();
    // Numbered params pair with the dataframe param sharing the numeric suffix
    // (JoinTables: keys2 resolves against table2, not the first table).
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

      // `param.name in node.inputValues` holds only when FuncNode seeded the slot,
      // i.e. the func has a dataframe input; the table-less case falls through to a Select Column node.
      if ((slotType === 'column' || slotType === 'column_list') && param.name in node.inputValues) {
        const names = this.extractColumnNames(value);
        node.inputValues[param.name] = slotType === 'column' ? (names[0] ?? '') : names.join(', ');
        continue;
      }

      // Inlined as a comma-separated editable value so emit → import → emit round-trips.
      if (slotType === 'string_list' && param.name in node.inputValues) {
        const items = Array.isArray(value) ?
          (value as unknown[]).map((v) => String(v).trim()).filter(Boolean) :
          String(value ?? '').split(',').map((s) => s.trim()).filter(Boolean);
        node.inputValues[param.name] = items.join(', ');
        continue;
      }

      // An array of primitives inlines as the editable array value; arrays of
      // calls/objects (e.g. column resolvers) fall through to generic resolution.
      if (slotType === 'list' && param.name in node.inputValues && Array.isArray(value) &&
          (value as unknown[]).every((v) => v === null || ['string', 'number', 'boolean'].includes(typeof v))) {
        node.inputValues[param.name] = (value as unknown[]).filter((v) => v !== null);
        continue;
      }

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

    // Title Open File nodes with their file — otherwise several are indistinguishable until they run.
    const path = String(node.inputValues['fullPath'] ?? '');
    if ((node.dgFuncName ?? '').toLowerCase() === 'openfile' && path !== '')
      node.label = `${node.label}: ${path.split('/').pop() || path}`;

    this.layer.set(node, maxSourceLayer + 1);
    return this.primaryOutput(node);
  }

  private addColumnSelection(fc: DG.FuncCall, typeName: string, ctx: CallContext): OutputRef | null {
    const names = this.columnNames(fc);
    const tableRef = this.explicitParentTable([fc], ctx) ?? ctx.contextTable;
    return this.addSelect(typeName, names, tableRef);
  }

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

  private extractColumnNames(value: unknown): string[] {
    if (value === null || value === undefined) return [];
    if (typeof value === 'string') return [value];
    if (Array.isArray(value)) return value.flatMap((v) => this.extractColumnNames(v));
    if (value instanceof DG.FuncCall) return this.columnNames(value);
    return [String(value)];
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

  private addNode(node: FlowNode): void {
    // Imported nodes start collapsed so the whole flow fits a view.
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

  /** Wires the variable's final ref into a real SetVar call — the only terminal node per variable. */
  private wireSetVar(varName: string): void {
    const ref = this.variables.get(varName.toLowerCase());
    if (!ref) return;
    const setVarFunc = this.findSetVarFunc();
    if (!setVarFunc) {
      this.warnings.push('SetVar function not found — variables will not be registered at run time');
      return;
    }
    const node = createNode(ensureFuncNodeType(setVarFunc));
    if (!node) return;
    node.label = `set: ${varName}`;
    node.inputValues['variableName'] = varName;
    this.addNode(node);
    this.connect(ref, node, 'value');
    this.layer.set(node, (this.layer.get(ref.node) ?? 0) + 1);
    this.setVarNodes.set(varName.toLowerCase(), node);
  }

  /** Order edge from a table-producing SetVar to each Select Table reading it by (normalized)
   *  name — nothing else forces the producer to run first. Creation scripts are acyclic. */
  private inferOrderEdges(): void {
    if (this.setVarNodes.size === 0) return;
    const byNorm = new Map<string, FlowNode>();
    for (const [name, node] of this.setVarNodes) {
      const key = normalizeName(name);
      if (key !== '' && !byNorm.has(key)) byNorm.set(key, node);
    }
    for (const node of this.nodes) {
      if (node.dgTypeName !== SELECT_TABLE_TYPE) continue;
      const producer = byNorm.get(normalizeName(node.properties['tableName']));
      if (!producer || producer === node) continue;
      this.connections.push({
        source: producer, sourceKey: EXEC_OUT_KEY,
        target: node, targetKey: EXEC_IN_KEY, order: true,
      });
    }
  }

  private _setVarFunc: DG.Func | null | undefined = undefined;
  private findSetVarFunc(): DG.Func | null {
    if (this._setVarFunc === undefined) {
      try {
        this._setVarFunc = DG.Func.find({name: 'SetVar'})[0] ?? null;
      } catch {
        this._setVarFunc = null;
      }
    }
    return this._setVarFunc;
  }

  private layout(): void {
    const dataEdges = this.connections
      .filter((c) => !c.order)
      .map((c) => ({source: c.source, target: c.target}));
    layoutGraph(this.nodes, dataEdges, this.layer);
  }
}

function isPrimitive(v: unknown): v is string | number | boolean {
  return typeof v === 'string' || typeof v === 'number' || typeof v === 'boolean';
}

/** Matches names across the name↔friendlyName convention: "Mol1KLocal" ≡ "mol1K local". */
function normalizeName(s: unknown): string {
  return String(s ?? '').toLowerCase().replace(/[^a-z0-9]/g, '');
}


/** Removes a trailing `// …` comment (the `//{"timestamp"}` metadata), quote-aware so URLs survive. */
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
