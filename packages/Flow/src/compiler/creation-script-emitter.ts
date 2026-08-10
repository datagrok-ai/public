/** Compiles the FuncFlow graph into a Datagrok creation script, using the
 *  platform's own serializer (`DG.Func.prepare(params).toString()`) — Dart-backed,
 *  so a live backend is required. */

import * as DG from 'datagrok-api/dg';

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode, FlowConnection, isExecKey} from '../rete/scheme';
import {FuncNode} from '../rete/nodes/func-node';
import {topologicalSort} from './topological-sort';

const PASSTHROUGH_SUFFIX = '__pt';

export interface CreationScriptResult {
  script: string;
  warnings: string[];
}

/** A resolved argument value; `var` renders as a bare `GetVar` reference — the
 *  only way to render a variable (a string or a real DataFrame both render quoted). */
type ArgValue =
  | {kind: 'var'; name: string}
  | {kind: 'literal'; value: unknown}
  | {kind: 'funcCall'; call: DG.FuncCall}
  | {kind: 'skip'};

/** An optional emitted line plus the table variable it builds/mutates (`owner`). */
interface NodeEmission {
  line: string | null;
  owner: string | null;
}
const NONE: NodeEmission = {line: null, owner: null};

interface FuncCache {
  setVar: DG.Func | null;
  getVar: DG.Func | null;
}

/** Utility node type names that have no creation-script (function) form. */
const UNSUPPORTED_UTILITY_TYPES = new Set([
  'Utilities/Log', 'Utilities/Info', 'Utilities/Warning', 'Utilities/Add Table View',
  'Utilities/ToString', 'Utilities/FromJSON', 'Utilities/ToJSON',
]);

export function emitCreationScript(flow: FlowEditor): CreationScriptResult {
  return new CreationScriptEmitter(flow).emit();
}

/** A single table's creation script, keyed by the variable name that builds it. */
export interface TableCreationScript {
  variableName: string;
  /** Each line suffixed with a strictly increasing `//{"timestamp": …}` comment. */
  script: string;
}

export interface PerTableResult {
  /** One entry per requested variable name, in the same order. */
  tables: TableCreationScript[];
  warnings: string[];
  /** Emitted lines the per-table split couldn't attribute to a requested table. */
  unassigned: string[];
}

/** Compile a separate creation script per table, splitting the combined cascade
 *  by the variable each line builds or mutates. */
export function emitCreationScriptsForTables(
  flow: FlowEditor, tableVarNames: string[],
): PerTableResult {
  return new CreationScriptEmitter(flow).emitPerTable(tableVarNames);
}

class CreationScriptEmitter {
  private readonly flow: FlowEditor;
  private readonly funcs: FuncCache;
  private readonly warnings: string[] = [];

  /** "<nodeId>:<outputKey>" → the value its slot yields. */
  private readonly outExpr = new Map<string, ArgValue>();
  /** target nodeId → inputKey → incoming connection. */
  private readonly incoming = new Map<string, Map<string, FlowConnection>>();
  /** "<sourceId>:<sourceOutput>" → consumer (nodeId,inputKey) list. */
  private readonly consumers = new Map<string, Array<{nodeId: string; inputKey: string}>>();
  /** producer nodeId → the variable name an anchoring SetVar gives it. */
  private readonly anchorName = new Map<string, string>();
  private readonly usedNames = new Set<string>();
  /** Emitted lines in topological order, with the owning table variable. */
  private readonly records: Array<{owner: string | null; line: string}> = [];
  private walked = false;

  constructor(flow: FlowEditor) {
    this.flow = flow;
    this.funcs = {
      setVar: findFunc('SetVar'),
      getVar: findFunc('GetVar'),
    };
  }

  private run(): void {
    if (this.walked) return;
    this.walked = true;
    const sorted = this.safeTopoSort();
    if (sorted === null) return;
    this.indexConnections();
    this.computeAnchors();
    for (const nodeId of sorted) {
      const node = this.flow.getNodeById(nodeId);
      if (!node) continue;
      const {line, owner} = this.emitNode(node);
      if (line) this.records.push({owner, line});
    }
  }

  emit(): CreationScriptResult {
    this.run();
    return {script: this.records.map((r) => r.line).join('\n'), warnings: this.warnings};
  }

  /** Split the emitted lines by owning table variable. Timestamps come from a
   *  single global counter across ALL tables — never per-table — because the
   *  value tells Datagrok which table to build first. */
  emitPerTable(tableVarNames: string[]): PerTableResult {
    this.run();
    const targets = new Set(tableVarNames);
    const byOwner = new Map<string, string[]>();
    const unassigned: string[] = [];
    let ts = Date.now();
    for (const {owner, line} of this.records) {
      const stamped = `${line} //{"timestamp": ${ts++}}`;
      if (owner !== null && targets.has(owner)) {
        let arr = byOwner.get(owner);
        if (!arr) byOwner.set(owner, arr = []);
        arr.push(stamped);
      } else
        unassigned.push(stamped);
    }
    const tables = tableVarNames.map((variableName) => ({
      variableName,
      script: (byOwner.get(variableName) ?? []).join('\n'),
    }));
    return {tables, warnings: this.warnings, unassigned};
  }

  private safeTopoSort(): string[] | null {
    try {
      return topologicalSort(this.flow);
    } catch (e) {
      this.warnings.push((e as Error).message);
      return null;
    }
  }

  private indexConnections(): void {
    for (const c of this.flow.getConnections()) {
      let perNode = this.incoming.get(c.target);
      if (!perNode) this.incoming.set(c.target, perNode = new Map());
      perNode.set(String(c.targetInput), c);

      const key = `${c.source}:${String(c.sourceOutput)}`;
      let list = this.consumers.get(key);
      if (!list) this.consumers.set(key, list = []);
      list.push({nodeId: c.target, inputKey: String(c.targetInput)});
    }
  }

  /** Propagate each anchor's variable name back through the pass-through chain
   *  to the producer at its head. Anchors are SetVar AND Output nodes — an
   *  Output's paramName IS its variable name. */
  private computeAnchors(): void {
    for (const node of this.flow.getNodes()) {
      let varName: string;
      let inputKey: string;
      if (isSetVar(node)) {
        varName = String(node.inputValues['variableName'] ?? '').trim();
        inputKey = 'value';
      } else if (node.dgNodeType === 'output') {
        varName = String(node.properties['paramName'] ?? '').trim();
        inputKey = dataInputKeys(node)[0] ?? 'value';
      } else
        continue;
      if (varName === '') continue;
      this.usedNames.add(varName);
      const conn = this.incoming.get(node.id)?.get(inputKey);
      if (!conn) continue;
      const producer = this.walkToProducer(conn.source, String(conn.sourceOutput), new Set());
      if (producer && !this.anchorName.has(producer))
        this.anchorName.set(producer, varName);
    }
  }

  /** Follow a value back through pass-through outputs to the func node that produced it. */
  private walkToProducer(nodeId: string, outKey: string, seen: Set<string>): string | null {
    const guard = `${nodeId}:${outKey}`;
    if (seen.has(guard)) return null;
    seen.add(guard);

    const node = this.flow.getNodeById(nodeId);
    if (!node) return null;
    if (!outKey.endsWith(PASSTHROUGH_SUFFIX)) {
      // Only func producers get renamed.
      return node.dgNodeType === 'func' && !isSetVar(node) && !isGetVar(node) ? nodeId : null;
    }
    const inputName = FuncNode.passthroughInputName(outKey);
    if (!inputName) return null;
    const conn = this.incoming.get(nodeId)?.get(inputName);
    if (!conn) return null;
    return this.walkToProducer(conn.source, String(conn.sourceOutput), seen);
  }

  private emitNode(node: FlowNode): NodeEmission {
    switch (node.dgNodeType) {
    case 'input':
      return this.emitInput(node);
    case 'output':
      return this.emitOutput(node);
    case 'utility':
      return this.emitUtility(node);
    case 'func':
      return this.emitFunc(node);
    default:
      return NONE;
    }
  }

  /** Input param node: every output slot resolves to a bare GetVar reference. */
  private emitInput(node: FlowNode): NodeEmission {
    const paramName = String(node.properties['paramName'] ?? 'input');
    for (const key of dataOutputKeys(node))
      this.outExpr.set(`${node.id}:${key}`, {kind: 'var', name: paramName});
    return NONE;
  }

  private emitOutput(node: FlowNode): NodeEmission {
    const paramName = String(node.properties['paramName'] ?? 'result');
    const firstInKey = dataInputKeys(node)[0] ?? 'value';
    const arg = this.resolveInput(node, firstInKey);
    if (arg.kind === 'skip') return NONE;
    if (arg.kind === 'var' && arg.name === paramName) return NONE;
    return {line: this.assignmentLine(paramName, arg), owner: paramName};
  }

  private emitUtility(node: FlowNode): NodeEmission {
    const type = node.dgTypeName ?? '';
    const firstOut = dataOutputKeys(node)[0];

    if (type === 'Debug/Breakpoint' || node.label === 'Breakpoint') {
      const arg = this.resolveInput(node, dataInputKeys(node)[0] ?? 'in');
      for (const key of dataOutputKeys(node)) this.outExpr.set(`${node.id}:${key}`, arg);
      return NONE;
    }

    if (type === 'Utilities/Select Table') {
      this.setOut(node, firstOut, {kind: 'literal', value: String(node.properties['tableName'] ?? '')});
      return NONE;
    }
    if (type === 'Utilities/Select Column') {
      this.setOut(node, firstOut, {kind: 'literal', value: String(node.properties['columnName'] ?? '')});
      return NONE;
    }
    if (type === 'Utilities/Select Columns') {
      this.setOut(node, firstOut, {kind: 'literal', value: splitList(node.properties['columnNames'])});
      return NONE;
    }
    if (type.startsWith('Constants/')) {
      this.setOut(node, firstOut, {kind: 'literal', value: constantValue(node)});
      return NONE;
    }

    if (UNSUPPORTED_UTILITY_TYPES.has(type) || firstOut) {
      this.warnings.push(`Node "${node.label}" (${type || 'utility'}) has no creation-script equivalent — skipped`);
      for (const key of dataOutputKeys(node)) this.outExpr.set(`${node.id}:${key}`, {kind: 'skip'});
    }
    return NONE;
  }

  private emitFunc(node: FlowNode): NodeEmission {
    if (isSetVar(node)) return this.emitSetVarNode(node);
    if (isGetVar(node)) {
      const name = String(node.inputValues['variableName'] ?? '').trim();
      for (const key of dataOutputKeys(node)) this.setOut(node, key, {kind: 'var', name});
      return NONE;
    }

    const func = node.dgFunc;
    if (!func) {
      this.warnings.push(`Func node "${node.label}" has no underlying function — skipped`);
      return NONE;
    }

    // A wrapped func's node inputs don't match the real signature — there is no
    // faithful creation-script call.
    if (node.funcWrapper) {
      this.warnings.push(`Node "${node.label}" wraps ${func.name} with reshaped inputs — ` +
        'no creation-script equivalent, skipped');
      for (const key of dataOutputKeys(node)) this.outExpr.set(`${node.id}:${key}`, {kind: 'skip'});
      return NONE;
    }

    const params = this.buildParams(node);
    let call: DG.FuncCall;
    try {
      call = func.prepare(params);
    } catch (e) {
      this.warnings.push(`Could not prepare "${node.label}": ${(e as Error).message}`);
      return NONE;
    }
    this.forwardPassthroughs(node);

    // Assignment (`X = f(...)`) iff a *real* (non-pass-through) output is
    // consumed/anchored — not merely declared: AddNewColumn rarely uses its
    // real output and must stay a bare statement preserving the in-place mutation.
    const realKeys = realOutputKeys(node);
    const realConsumed = this.anchorName.has(node.id) ||
      realKeys.some((k) => (this.consumers.get(`${node.id}:${k}`)?.length ?? 0) > 0);

    // Bare mutator: the line belongs to the table threaded through its primary
    // dataframe input.
    if (!realConsumed)
      return {line: safeToString(call, node, this.warnings), owner: this.primaryTableOwner(node)};

    const name = this.anchorName.get(node.id) ?? uniqueName(toCamelCase(node.label) || 'result', this.usedNames);
    this.usedNames.add(name);
    for (let i = 0; i < realKeys.length; i++) {
      const varName = i === 0 ? name : uniqueName(`${name}_${realKeys[i]}`, this.usedNames);
      this.usedNames.add(varName);
      this.setOut(node, realKeys[i], {kind: 'var', name: varName});
    }
    if (realKeys.length > 1)
      this.warnings.push(`Func "${node.label}" has multiple outputs — only the first is named in the assignment`);

    return {line: this.assignmentLine(name, {kind: 'funcCall', call}), owner: name};
  }

  private emitSetVarNode(node: FlowNode): NodeEmission {
    const varName = String(node.inputValues['variableName'] ?? '').trim();
    if (varName === '') return NONE;
    const arg = this.resolveInput(node, 'value');
    if (arg.kind === 'skip') return NONE;
    // Already produced under this exact name — no redundant `Mol1K = Mol1K`.
    if (arg.kind === 'var' && arg.name === varName) return NONE;
    return {line: this.assignmentLine(varName, arg), owner: varName};
  }

  /** The table variable a bare call mutates — the variable reference resolved
   *  on its first dataframe input, or null. */
  private primaryTableOwner(node: FlowNode): string | null {
    const func = node.dgFunc;
    if (!func) return null;
    for (const param of func.inputs) {
      if (slotTypeOf(node, param.name) !== 'dataframe') continue;
      const arg = this.resolveInput(node, param.name);
      return arg.kind === 'var' ? arg.name : null;
    }
    return null;
  }

  /** Forward each pass-through output to its input's resolved value, so a
   *  mutator chain keeps referencing the same variable. */
  private forwardPassthroughs(node: FlowNode): void {
    for (const key of Object.keys(node.outputs)) {
      if (isExecKey(key) || !key.endsWith(PASSTHROUGH_SUFFIX)) continue;
      const inputName = FuncNode.passthroughInputName(key);
      if (!inputName) continue;
      this.outExpr.set(`${node.id}:${key}`, this.resolveInput(node, inputName));
    }
  }

  /** Build the `prepare()` params map for a func node from its resolved inputs. */
  private buildParams(node: FlowNode): Record<string, unknown> {
    const params: Record<string, unknown> = {};
    const func = node.dgFunc!;
    for (const param of func.inputs) {
      const key = param.name;
      if (!(key in node.inputs)) continue;
      const arg = this.resolveInput(node, key);
      const value = this.toPrepareValue(arg);
      if (value === OMIT) continue;
      // The importer seeds every primitive with `''`-style defaults — passing
      // one would leak `sheetName = ""` etc.; drop empty literals so the
      // serializer falls back to the param's true default.
      if (arg.kind === 'literal' && isEmptyLiteral(value)) continue;
      params[key] = value;
    }
    return params;
  }

  private resolveInput(node: FlowNode, key: string): ArgValue {
    const conn = this.incoming.get(node.id)?.get(key);
    if (conn)
      return this.outExpr.get(`${conn.source}:${String(conn.sourceOutput)}`) ?? {kind: 'skip'};

    if (key in node.inputValues) {
      const slotType = slotTypeOf(node, key);
      const raw = node.inputValues[key];
      if (slotType === 'column') {
        const name = String(raw ?? '').trim();
        return name === '' ? {kind: 'skip'} : {kind: 'literal', value: name};
      }
      if (slotType === 'column_list') {
        const names = splitList(raw);
        return names.length === 0 ? {kind: 'skip'} : {kind: 'literal', value: names};
      }
      if (slotType === 'string_list') {
        const items = splitList(raw);
        return items.length === 0 ? {kind: 'skip'} : {kind: 'literal', value: items};
      }
      if (slotType === 'list') {
        const items = Array.isArray(raw) ? raw : splitList(raw);
        return items.length === 0 ? {kind: 'skip'} : {kind: 'literal', value: items};
      }
      if (raw === undefined || raw === null) return {kind: 'skip'};
      return {kind: 'literal', value: raw};
    }
    return {kind: 'skip'};
  }

  /** Convert an `ArgValue` into the JS value handed to `prepare()`, or `OMIT`. */
  private toPrepareValue(arg: ArgValue): unknown {
    switch (arg.kind) {
    case 'literal':
      return arg.value;
    case 'funcCall':
      return arg.call;
    case 'var': {
      if (this.funcs.getVar) return this.funcs.getVar.prepare({variableName: arg.name});
      this.warnOnce('GetVar function not found — variable references emitted as quoted names');
      return arg.name;
    }
    default:
      return OMIT;
    }
  }

  /** `Name = <value>` using SetVar's own serializer; falls back to concatenation. */
  private assignmentLine(name: string, arg: ArgValue): string {
    const value = this.toPrepareValue(arg);
    if (this.funcs.setVar && value !== OMIT) {
      try {
        return this.funcs.setVar.prepare({variableName: name, value}).toString();
      } catch {
        // fall through to manual concatenation
      }
    }
    const rhs = arg.kind === 'funcCall' ? safeToString(arg.call, null, this.warnings) : String(value);
    return `${name} = ${rhs}`;
  }

  private setOut(node: FlowNode, key: string | undefined, value: ArgValue): void {
    if (key) this.outExpr.set(`${node.id}:${key}`, value);
  }

  private readonly warnedOnce = new Set<string>();
  private warnOnce(msg: string): void {
    if (this.warnedOnce.has(msg)) return;
    this.warnedOnce.add(msg);
    this.warnings.push(msg);
  }
}

/** Sentinel: "do not pass this param to prepare()". */
const OMIT = Symbol('omit');

function findFunc(name: string): DG.Func | null {
  try {
    return DG.Func.find({name})[0] ?? null;
  } catch {
    return null;
  }
}

function safeToString(call: DG.FuncCall, node: FlowNode | null, warnings: string[]): string {
  try {
    return call.toString();
  } catch (e) {
    warnings.push(`Could not serialize ${node ? `"${node.label}"` : 'a call'}: ${(e as Error).message}`);
    return node ? `// ${node.label}` : '';
  }
}

function isSetVar(node: FlowNode): boolean {
  return (node.dgFunc?.name?.toLowerCase() ?? '') === 'setvar';
}
function isGetVar(node: FlowNode): boolean {
  return (node.dgFunc?.name?.toLowerCase() ?? '') === 'getvar';
}

/** A seeded-empty literal to drop so the param's true default applies. */
function isEmptyLiteral(value: unknown): boolean {
  if (value === '' || value === null || value === undefined) return true;
  return Array.isArray(value) && value.length === 0;
}

function dataOutputKeys(node: FlowNode): string[] {
  return Object.keys(node.outputs).filter((k) => !isExecKey(k));
}
function dataInputKeys(node: FlowNode): string[] {
  return Object.keys(node.inputs).filter((k) => !isExecKey(k));
}
function realOutputKeys(node: FlowNode): string[] {
  return dataOutputKeys(node).filter((k) => !k.endsWith(PASSTHROUGH_SUFFIX));
}

function slotTypeOf(node: FlowNode, key: string): string | undefined {
  return (node.inputs as Record<string, {socket: {dgType: string}} | undefined>)[key]?.socket.dgType;
}

function constantValue(node: FlowNode): unknown {
  const type = node.dgTypeName ?? '';
  const raw = node.properties['value'];
  if (type === 'Constants/Int') return Math.round(Number(raw ?? 0));
  if (type === 'Constants/Double') return Number(raw ?? 0);
  if (type === 'Constants/Boolean') return raw === true || raw === 'true';
  if (type === 'Constants/List') return splitList(raw);
  return String(raw ?? '');
}

function splitList(raw: unknown): string[] {
  if (Array.isArray(raw)) return raw.map((s) => String(s).trim()).filter(Boolean);
  return String(raw ?? '').split(',').map((s) => s.trim()).filter(Boolean);
}

function toCamelCase(s: string): string {
  return s
    .replace(/[^a-zA-Z0-9]+/g, ' ')
    .trim()
    .split(' ')
    .map((word, i) => (i === 0 ? word.toLowerCase() : word.charAt(0).toUpperCase() + word.slice(1).toLowerCase()))
    .join('');
}

function uniqueName(base: string, used: Set<string>): string {
  if (!base) base = 'v';
  if (!used.has(base)) return base;
  let i = 2;
  while (used.has(`${base}${i}`)) i++;
  return `${base}${i}`;
}
