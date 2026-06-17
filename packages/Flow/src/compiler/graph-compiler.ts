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
import {FlowNode, FlowConnection} from '../rete/scheme';
import {FuncNode} from '../rete/nodes/func-node';
import {topologicalSort} from './topological-sort';

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
  /** real output slot key → variable name we declared. */
  outputs: Map<string, string>;
  /** Snapshot of node.properties at compile time. */
  properties: Record<string, unknown>;
  /** Snapshot of node.inputValues (hardcoded primitive defaults). */
  inputValues: Record<string, unknown>;
}

const PASSTHROUGH_SUFFIX = '__pt';

export function compileGraph(flow: FlowEditor): CompiledStep[] {
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
      // Map every output slot of an input node to the param name.
      for (const key of Object.keys(node.outputs))
        outputVarMap.set(`${nodeId}:${key}`, paramName);
      const firstOutKey = Object.keys(node.outputs)[0] ?? 'value';
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
      const firstInKey = Object.keys(node.inputs)[0] ?? 'value';
      const inputExpr = resolveInputExpr(nodeId, firstInKey, incoming, outputVarMap);
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
        const inputExpr = resolveInputExpr(nodeId, 'in', incoming, outputVarMap);
        for (const key of Object.keys(node.outputs))
          outputVarMap.set(`${nodeId}:${key}`, inputExpr);
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

      const step = compileUtilityNode(node, incoming, outputVarMap, usedVarNames);
      steps.push(step);
      for (const key of Object.keys(node.outputs))
        outputVarMap.set(`${nodeId}:${key}`, step.variableName);
      continue;
    }

    // -------- func step --------
    const funcName = node.dgFuncName || node.label;
    const ptCount = node.passthroughCount;
    const inputKeys = Object.keys(node.inputs);
    const outputKeys = Object.keys(node.outputs);

    // Variable name based on title + first real output's name.
    const realOutputKeys = outputKeys.filter((k) => !k.endsWith(PASSTHROUGH_SUFFIX));
    const firstRealKey = realOutputKeys[0];
    let varBase = toCamelCase(node.label);
    if (firstRealKey) varBase = `${toCamelCase(node.label)}_${firstRealKey}`;
    const varName = uniqueVarName(varBase, usedVarNames);
    usedVarNames.add(varName);

    // Resolve inputs: connected → expression from outputVarMap, else hardcoded.
    const inputMap = new Map<string, string>();
    for (const key of inputKeys) {
      const conn = incoming.get(nodeId)?.get(key);
      if (conn) {
        inputMap.set(key, resolveConnExpr(conn, outputVarMap));
      } else if (key in node.inputValues) {
        const val = node.inputValues[key];
        const slotType = (node.inputs as Record<string, {socket: {dgType: string}} | undefined>)[key]?.socket.dgType;
        if (val !== undefined)
          inputMap.set(key, formatLiteral(val, slotType ?? 'dynamic'));
      }
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
        const realCount = realOutputKeys.length;
        const outVarName = realCount === 1 ? varName : `${varName}_${key}`;
        outputMap.set(key, outVarName);
        outputVarMap.set(`${nodeId}:${key}`, outVarName);
      }
    }

    steps.push({
      nodeId, nodeType: 'func', funcName,
      variableName: varName,
      inputs: inputMap, outputs: outputMap,
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
): CompiledStep {
  const varName = uniqueVarName(toCamelCase(node.label), usedVarNames);
  usedVarNames.add(varName);

  const inputMap = new Map<string, string>();
  for (const key of Object.keys(node.inputs)) {
    const conn = incoming.get(node.id)?.get(key);
    if (conn) inputMap.set(key, resolveConnExpr(conn, outputVarMap));
  }

  const firstOutKey = Object.keys(node.outputs)[0] ?? 'value';
  return {
    nodeId: node.id, nodeType: 'utility', funcName: utilityKind(node),
    variableName: varName,
    inputs: inputMap,
    outputs: new Map([[firstOutKey, varName]]),
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

function resolveInputExpr(
  nodeId: string, inputKey: string,
  incoming: Map<string, Map<string, FlowConnection>>,
  outputVarMap: Map<string, string>,
): string {
  const conn = incoming.get(nodeId)?.get(inputKey);
  if (!conn) return 'undefined';
  return resolveConnExpr(conn, outputVarMap);
}

function resolveConnExpr(c: FlowConnection, outputVarMap: Map<string, string>): string {
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

function uniqueVarName(base: string, used: Set<string>): string {
  if (!base) base = 'v';
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
