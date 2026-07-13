/** Compiles the FuncFlow graph into a flat ordered list of `CompiledStep`s.
 *
 * Rete-aware: input/output keys are strings; pass-through outputs use the
 * `<inputName>__pt` convention (see `func-node.ts`). The compiler resolves
 * input expressions by walking incoming connections and looking up the
 * source node's output variable in `outputVarMap`.
 *
 * The output expressions are JavaScript snippets that the script-emitter
 * inserts into the generated body. */

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode, FlowConnection, isExecKey} from '../rete/scheme';
import {FuncNode, defaultTableParam} from '../rete/nodes/func-node';
import {topologicalSort} from './topological-sort';
import {stringListToArrayLiteral} from '../types/type-map';

export type StepKind = 'input' | 'output' | 'utility' | 'func';

export interface CompiledStep {
  nodeId: string;
  nodeType: StepKind;
  /** Qualified DG function name for `func` steps; node title for utilities. */
  funcName: string;
  /** Variable name we declared in the body (or empty for non-producing steps). */
  variableName: string;
  /** input slot key → source expression to feed in. */
  inputs: Map<string, string>;
  /** real output slot key → expression that yields its value. For a func step
   *  this is the call variable (single-output funcs return the value directly)
   *  or a property read off it (`<var>.<outputName>` for a multi-output func —
   *  `grok.functions.call` resolves to an object keyed by the declared output
   *  names). Utilities/inputs map to their declared variable. */
  outputs: Map<string, string>;
  /** With `cloneDataframeInputs` (instrumented runs): snapshot variable →
   *  original source expression, one per connected dataframe input. The call,
   *  the pass-through outputs, and the live-stash all use the snapshot, so an
   *  in-place function mutates its own copy — never the upstream node's
   *  captured value (which must stay "the state at that node"). */
  cloneInputs?: Map<string, string>;
  /** Snapshot of node.properties at compile time. */
  properties: Record<string, unknown>;
  /** Snapshot of node.inputValues (hardcoded primitive defaults). */
  inputValues: Record<string, unknown>;
}

const PASSTHROUGH_SUFFIX = '__pt';

/** Every node that must run to produce `targetId`'s output: the target plus all
 *  its transitive predecessors (walking connections backward — data,
 *  pass-through, and order edges alike). The set is closed under "ancestors",
 *  so emitting only these steps yields a self-contained sub-script. */
export function sliceUpTo(flow: FlowEditor, targetId: string): Set<string> {
  const connections = flow.getConnections();
  const result = new Set<string>();
  const stack = [targetId];
  while (stack.length > 0) {
    const id = stack.pop()!;
    if (result.has(id)) continue;
    result.add(id);
    for (const c of connections)
      if (c.target === id) stack.push(c.source);
  }
  return result;
}

/** Every node whose result can be affected by a change at `rootId`: the node
 *  itself plus all its transitive successors (walking connections forward —
 *  data, pass-through, and order edges alike). The mirror of {@link sliceUpTo};
 *  used to invalidate run results precisely after a graph edit. */
export function sliceDownFrom(flow: FlowEditor, rootId: string): Set<string> {
  const connections = flow.getConnections();
  const result = new Set<string>();
  const stack = [rootId];
  while (stack.length > 0) {
    const id = stack.pop()!;
    if (result.has(id)) continue;
    result.add(id);
    for (const c of connections)
      if (c.source === id) stack.push(c.target);
  }
  return result;
}

/** When `liveBoundary` is set, a connection whose *source* node is not in the
 *  set resolves to a `_ffLive(nodeId, outputKey)` call — reading the value that
 *  node produced in a prior instrumented run from the live-value registry —
 *  instead of an in-script variable. This is what lets a single node be re-run
 *  in isolation (its upstream inputs come from the registry, not a re-run). */
export interface CompileOptions {
  /** Snapshot-clone every connected dataframe input of a func step (see
   *  `CompiledStep.cloneInputs`). On for instrumented emission: many platform
   *  functions transform tables **in place**, and without per-step clones they
   *  would mutate the very instance stashed as the upstream node's live value —
   *  corrupting node previews (a node's preview must show the state AT that
   *  node) and making autorun slice re-runs non-idempotent (e.g. re-running a
   *  descriptor calc appends its columns a second time). Off for clean scripts,
   *  which run once from scratch and keep the platform's in-place idiom. */
  cloneDataframeInputs?: boolean;
}

export function compileGraph(
  flow: FlowEditor, liveBoundary?: Set<string>, opts?: CompileOptions,
): CompiledStep[] {
  const sortedIds = topologicalSort(flow);
  const connections = flow.getConnections();

  // Index incoming connections per (targetId → targetInput → connection)
  const incoming = new Map<string, Map<string, FlowConnection>>();
  for (const c of connections) {
    let perNode = incoming.get(c.target);
    if (!perNode) incoming.set(c.target, perNode = new Map());
    perNode.set(String(c.targetInput), c);
  }

  const steps: CompiledStep[] = [];
  /** "<nodeId>:<outputKey>" → expression that yields its value. */
  const outputVarMap = new Map<string, string>();
  const usedVarNames = new Set<string>();

  for (const nodeId of sortedIds) {
    const node = flow.getNodeById(nodeId);
    if (!node) continue;

    const kind = node.dgNodeType;

    if (kind === 'input') {
      const paramName = String(node.properties['paramName'] ?? 'input');
      usedVarNames.add(paramName);
      // Map every data output slot of an input node to the param name (the
      // exec-out ordering port carries no data and is skipped).
      const outKeys = Object.keys(node.outputs).filter((k) => !isExecKey(k));
      for (const key of outKeys)
        outputVarMap.set(`${nodeId}:${key}`, paramName);
      const firstOutKey = outKeys[0] ?? 'value';
      steps.push({
        nodeId, nodeType: 'input', funcName: '',
        variableName: paramName,
        inputs: new Map(),
        outputs: new Map([[firstOutKey, paramName]]),
        properties: {...node.properties},
        inputValues: {...node.inputValues},
      });
      continue;
    }

    if (kind === 'output') {
      const paramName = String(node.properties['paramName'] ?? 'result');
      const firstInKey = Object.keys(node.inputs).find((k) => !isExecKey(k)) ?? 'value';
      const inputExpr = resolveInputExpr(nodeId, firstInKey, incoming, outputVarMap, liveBoundary);
      steps.push({
        nodeId, nodeType: 'output', funcName: '',
        variableName: paramName,
        inputs: new Map([[firstInKey, inputExpr]]),
        outputs: new Map(),
        properties: {...node.properties},
        inputValues: {...node.inputValues},
      });
      continue;
    }

    if (kind === 'utility') {
      // Breakpoint: pure pass-through — output slot resolves to its input expr.
      if (utilityKind(node) === 'Breakpoint') {
        const inputExpr = resolveInputExpr(nodeId, 'in', incoming, outputVarMap, liveBoundary);
        for (const key of Object.keys(node.outputs))
          if (!isExecKey(key)) outputVarMap.set(`${nodeId}:${key}`, inputExpr);
        steps.push({
          nodeId, nodeType: 'utility', funcName: 'Breakpoint',
          variableName: '',
          inputs: new Map([['in', inputExpr]]),
          outputs: new Map(),
          properties: {...node.properties},
          inputValues: {...node.inputValues},
        });
        continue;
      }

      const step = compileUtilityNode(node, incoming, outputVarMap, usedVarNames, liveBoundary);
      steps.push(step);
      for (const key of Object.keys(node.outputs))
        if (!isExecKey(key)) outputVarMap.set(`${nodeId}:${key}`, step.variableName);
      continue;
    }

    // -------- func step --------
    // Exec ordering ports carry no data — exclude them from data resolution.
    const funcName = node.dgFuncName || node.label;
    const ptCount = node.passthroughCount;
    const inputKeys = Object.keys(node.inputs).filter((k) => !isExecKey(k));
    const outputKeys = Object.keys(node.outputs).filter((k) => !isExecKey(k));

    // Variable name based on title + first real output's name.
    const realOutputKeys = outputKeys.filter((k) => !k.endsWith(PASSTHROUGH_SUFFIX));
    const firstRealKey = realOutputKeys[0];
    let varBase = toCamelCase(node.label);
    if (firstRealKey) varBase = `${toCamelCase(node.label)}_${firstRealKey}`;
    const varName = uniqueVarName(varBase, usedVarNames);
    usedVarNames.add(varName);

    // Resolve inputs: connected → expression from outputVarMap, else hardcoded.
    // Column / column-list values are deferred to a second pass: they compile to
    // `table.col(...)` against an associated dataframe input, which must be
    // resolved first.
    const inputMap = new Map<string, string>();
    for (const key of inputKeys) {
      const conn = incoming.get(nodeId)?.get(key);
      if (conn) {
        inputMap.set(key, resolveConnExpr(conn, outputVarMap, liveBoundary));
        continue;
      }
      if (!(key in node.inputValues)) continue;
      const slotType = slotTypeOf(node, key);
      if (slotType === 'column' || slotType === 'column_list') continue; // pass 2
      const val = node.inputValues[key];
      if (slotType === 'string_list') {
        // Comma-separated → JS array of trimmed, non-empty strings. Empty → omit
        // so the function falls back to its own default (don't force `[]`).
        const lit = stringListToArrayLiteral(val);
        if (lit !== '[]') inputMap.set(key, lit);
        continue;
      }
      if (slotType === 'list') {
        // Native DG List input (incl. list<string> params) → a JS array; emit
        // as a JSON array literal. Empty → omit so the function default applies.
        if (Array.isArray(val) && val.length > 0) inputMap.set(key, JSON.stringify(val));
        continue;
      }
      if (val !== undefined)
        inputMap.set(key, formatLiteral(val, slotType ?? 'dynamic'));
    }
    // Snapshot-clone connected dataframe inputs (instrumented runs): rewrite
    // the input expression to a fresh snapshot variable BEFORE pass 2, so
    // inlined `table.col(...)` args select from the very instance the call
    // receives. The pass-through mapping below then also picks up the snapshot
    // — the (possibly mutated) copy flows on, the upstream value stays intact.
    const cloneMap = new Map<string, string>();
    if (opts?.cloneDataframeInputs) {
      for (const key of inputKeys) {
        if (!incoming.get(nodeId)?.get(key)) continue; // literals aren't shared
        if (slotTypeOf(node, key) !== 'dataframe') continue;
        const srcExpr = inputMap.get(key);
        if (!srcExpr || srcExpr === 'undefined') continue;
        const snapVar = uniqueVarName(`${varName}_${key}_in`, usedVarNames);
        usedVarNames.add(snapVar);
        cloneMap.set(snapVar, srcExpr);
        inputMap.set(key, snapVar);
      }
    }

    // Pass 2: inline unconnected column / column-list values as `table.col(...)`
    // expressions, resolving the table from the node's `columnTables` association.
    for (const key of inputKeys) {
      if (inputMap.has(key) || !(key in node.inputValues)) continue;
      const slotType = slotTypeOf(node, key);
      if (slotType !== 'column' && slotType !== 'column_list') continue;
      const raw = String(node.inputValues[key] ?? '').trim();
      if (raw === '') continue;
      const tableExpr = tableExprForColumnParam(node, key, inputMap);
      if (tableExpr) inputMap.set(key, columnSelectionExpr(slotType, raw, tableExpr));
    }

    // Map output slots: pass-throughs → corresponding input expr; real → varName.
    const outputMap = new Map<string, string>();
    for (const key of outputKeys) {
      const ptInput = key.endsWith(PASSTHROUGH_SUFFIX)
        ? FuncNode.passthroughInputName(key)
        : null;
      if (ptInput) {
        const inExpr = inputMap.get(ptInput) ?? 'undefined';
        outputVarMap.set(`${nodeId}:${key}`, inExpr);
      } else {
        // `grok.functions.call` returns the value directly for a single-output
        // func, but an object keyed by the declared output names for a
        // multi-output one — so a lone output reads straight off the call
        // variable, while several read a property (`<var>.<outputName>`). The
        // previous `<var>_<key>` scheme referenced variables never declared.
        const outExpr = realOutputKeys.length === 1 ? varName : `${varName}${propertyAccessor(key)}`;
        outputMap.set(key, outExpr);
        outputVarMap.set(`${nodeId}:${key}`, outExpr);
      }
    }

    steps.push({
      nodeId, nodeType: 'func', funcName,
      variableName: varName,
      inputs: inputMap, outputs: outputMap,
      cloneInputs: cloneMap.size > 0 ? cloneMap : undefined,
      properties: {...node.properties, _passthroughCount: ptCount},
      inputValues: {...node.inputValues},
    });
  }

  return steps;
}

function compileUtilityNode(
  node: FlowNode,
  incoming: Map<string, Map<string, FlowConnection>>,
  outputVarMap: Map<string, string>,
  usedVarNames: Set<string>,
  liveBoundary?: Set<string>,
): CompiledStep {
  const varName = uniqueVarName(toCamelCase(node.label), usedVarNames);
  usedVarNames.add(varName);

  const inputMap = new Map<string, string>();
  for (const key of Object.keys(node.inputs)) {
    if (isExecKey(key)) continue;
    const conn = incoming.get(node.id)?.get(key);
    if (conn) inputMap.set(key, resolveConnExpr(conn, outputVarMap, liveBoundary));
  }

  // Side-effect-only utilities (Log / Info / Warning) have no data output
  // socket — their step must declare NO outputs, or downstream instrumentation
  // (`__ff_stash`, output summaries) would reference an undeclared variable.
  const firstOutKey = Object.keys(node.outputs).find((k) => !isExecKey(k));
  return {
    nodeId: node.id, nodeType: 'utility', funcName: utilityKind(node),
    variableName: varName,
    inputs: inputMap,
    outputs: new Map(firstOutKey ? [[firstOutKey, varName]] : []),
    properties: {...node.properties},
    inputValues: {...node.inputValues},
  };
}

/** Stable utility node kind — the trailing segment of the registered type name
 *  (`Constants/String` → `String`). Labels are user-editable (constant nodes
 *  title themselves `const: <value>`), so emission dispatch must not read
 *  `node.label`; it stays only a fallback for nodes created outside the
 *  factory (tests). */
function utilityKind(node: FlowNode): string {
  return node.dgTypeName?.split('/').pop() ?? node.label;
}

function slotTypeOf(node: FlowNode, key: string): string | undefined {
  return (node.inputs as Record<string, {socket: {dgType: string}} | undefined>)[key]?.socket.dgType;
}

/** Resolve the table expression a column/column-list input selects from: the
 *  dataframe input recorded in the node's `columnTables` association, falling
 *  back (older graphs / hand-built nodes) to the numeric-suffix pairing then
 *  the first connected dataframe input. Returns undefined when no dataframe
 *  input is resolvable, in which case the column value is dropped. */
function tableExprForColumnParam(
  node: FlowNode, paramName: string, inputMap: Map<string, string>,
): string | undefined {
  const associations = node.properties['columnTables'] as Record<string, string> | undefined;
  const explicit = associations?.[paramName];
  if (explicit && inputMap.has(explicit)) return inputMap.get(explicit);

  const dataframeKeys = (Object.entries(node.inputs) as Array<[string, {socket: {dgType: string}} | undefined]>)
    .filter(([, inp]) => inp?.socket.dgType === 'dataframe')
    .map(([k]) => k);
  if (dataframeKeys.length === 0) return undefined;
  const fallback = defaultTableParam(paramName, dataframeKeys);
  if (fallback && inputMap.has(fallback)) return inputMap.get(fallback);
  for (const k of dataframeKeys) if (inputMap.has(k)) return inputMap.get(k);
  return undefined;
}

/** `table.col('name')` for a column, `[table.col('a'), table.col('b')]` for a
 *  column-list (comma-separated, trimmed) — identical to the Select Column(s)
 *  utility emission, so inlined columns and explicit nodes generate the same code. */
function columnSelectionExpr(slotType: string, raw: string, tableExpr: string): string {
  if (slotType === 'column') return `${tableExpr}.col('${escapeColumnName(raw)}')`;
  const names = raw.split(',').map((s) => s.trim()).filter(Boolean);
  return `[${names.map((n) => `${tableExpr}.col('${escapeColumnName(n)}')`).join(', ')}]`;
}

function escapeColumnName(name: string): string {
  return name.replace(/\\/g, '\\\\').replace(/'/g, '\\\'');
}

function resolveInputExpr(
  nodeId: string, inputKey: string,
  incoming: Map<string, Map<string, FlowConnection>>,
  outputVarMap: Map<string, string>,
  liveBoundary?: Set<string>,
): string {
  const conn = incoming.get(nodeId)?.get(inputKey);
  if (!conn) return 'undefined';
  return resolveConnExpr(conn, outputVarMap, liveBoundary);
}

function resolveConnExpr(
  c: FlowConnection, outputVarMap: Map<string, string>, liveBoundary?: Set<string>,
): string {
  // Source outside the compiled slice → read its captured value from the
  // live-value registry rather than an (absent) in-script variable.
  if (liveBoundary && !liveBoundary.has(c.source))
    return `_ffLive(${JSON.stringify(c.source)}, ${JSON.stringify(String(c.sourceOutput))})`;
  return outputVarMap.get(`${c.source}:${String(c.sourceOutput)}`) ?? 'undefined';
}

function toCamelCase(s: string): string {
  return s
    .replace(/[^a-zA-Z0-9]+/g, ' ')
    .trim()
    .split(' ')
    .map((word, i) => i === 0
      ? word.toLowerCase()
      : word.charAt(0).toUpperCase() + word.slice(1).toLowerCase(),
    )
    .join('');
}

const JS_IDENTIFIER_RE = /^[A-Za-z_$][A-Za-z0-9_$]*$/;

/** `.name` for a valid identifier, else `['name']` — reads a named output off a
 *  multi-output call result (`grok.functions.call` returns an object keyed by
 *  the declared output names). */
function propertyAccessor(key: string): string {
  return JS_IDENTIFIER_RE.test(key) ? `.${key}` : `[${JSON.stringify(key)}]`;
}

function uniqueVarName(base: string, used: Set<string>): string {
  if (!base) base = 'v';
  // `toCamelCase` keeps digits, so a func/node named "2 Inputs Flow" yields
  // "2InputsFlow" — illegal as a variable name (can't start with a digit).
  if (/^[0-9]/.test(base)) base = `_${base}`;
  if (!used.has(base)) return base;
  let i = 2;
  while (used.has(`${base}${i}`)) i++;
  return `${base}${i}`;
}

function formatLiteral(value: unknown, type: string): string {
  if (value === null || value === undefined) return 'null';
  switch (type) {
  case 'string':
    return JSON.stringify(String(value));
  case 'bool':
    return value ? 'true' : 'false';
  case 'int':
    return String(Math.round(Number(value)));
  case 'double':
  case 'num':
    return String(Number(value));
  default:
    if (typeof value === 'string') return JSON.stringify(value);
    return String(value);
  }
}
