import {LGraph, LGraphNode} from 'litegraph.js';
import {topologicalSort} from './topological-sort';
import {FuncFlowNode} from '../types/funcflow-node';

export interface CompiledStep {
  nodeId: number;
  nodeType: 'func' | 'input' | 'output' | 'utility';
  funcName: string;
  variableName: string;
  inputs: Map<string, string>; // paramName → source expression
  outputs: Map<string, string>; // paramName → variable name
  properties: Record<string, any>; // node.properties snapshot for widget values
}

/** Compiles the graph into an ordered list of execution steps */
export function compileGraph(graph: LGraph): CompiledStep[] {
  const sortedIds = topologicalSort(graph);
  const links = graph.links || {};

  const steps: CompiledStep[] = [];
  // Map: nodeId:outputSlotIndex → variable name expression
  const outputVarMap = new Map<string, string>();
  const usedVarNames = new Set<string>();

  for (const nodeId of sortedIds) {
    const node = graph.getNodeById(nodeId);
    if (!node) continue;

    const n = node as FuncFlowNode;
    const nodeType = n.dgNodeType || 'func';

    if (nodeType === 'input') {
      // Input nodes produce a variable from their param name
      const paramName = node.properties['paramName'] || 'input';
      usedVarNames.add(paramName);
      // Map all output slots to this param name
      if (node.outputs) {
        for (let i = 0; i < node.outputs.length; i++)
          outputVarMap.set(`${nodeId}:${i}`, paramName);
      }
      steps.push({
        nodeId,
        nodeType: 'input',
        funcName: '',
        variableName: paramName,
        inputs: new Map(),
        outputs: new Map([[node.outputs?.[0]?.name || 'value', paramName]]),
        properties: {...node.properties},
      });
      continue;
    }

    if (nodeType === 'output') {
      // Output nodes just capture their input as the output variable
      const paramName = node.properties['paramName'] || 'result';
      const inputExpr = resolveInputExpression(node, 0, links, outputVarMap);
      steps.push({
        nodeId,
        nodeType: 'output',
        funcName: '',
        variableName: paramName,
        inputs: new Map([['value', inputExpr]]),
        outputs: new Map(),
        properties: {...node.properties},
      });
      continue;
    }

    if (nodeType === 'utility') {
      const step = compileUtilityNode(node, links, outputVarMap, usedVarNames);
      if (step) {
        steps.push(step);
        // Map outputs
        if (node.outputs) {
          for (let i = 0; i < node.outputs.length; i++)
            outputVarMap.set(`${nodeId}:${i}`, step.variableName);
        }
      }
      continue;
    }

    // Function node
    const funcName = n.dgFuncName || node.title;

    // Generate variable name from func name + first output
    let varBase = toCamelCase(node.title);
    if (node.outputs && node.outputs.length > 0)
      varBase = `${toCamelCase(node.title)}_${node.outputs[0].name}`;

    const varName = uniqueVarName(varBase, usedVarNames);
    usedVarNames.add(varName);

    // Resolve inputs
    const inputMap = new Map<string, string>();
    if (node.inputs) {
      for (let i = 0; i < node.inputs.length; i++) {
        const inp = node.inputs[i];
        if (node.isInputConnected(i)) {
          const expr = resolveInputExpression(node, i, links, outputVarMap);
          inputMap.set(inp.name, expr);
        } else {
          // Use hardcoded value from widget
          const val = node.properties[`_input_${inp.name}`];
          if (val !== undefined)
            inputMap.set(inp.name, formatLiteral(val, inp.type as string));
        }
      }
    }

    // Map output slots to variable names
    const outputMap = new Map<string, string>();
    if (node.outputs) {
      for (let i = 0; i < node.outputs.length; i++) {
        const outVarName = node.outputs.length === 1 ? varName : `${varName}_${node.outputs[i].name}`;
        outputMap.set(node.outputs[i].name, outVarName);
        outputVarMap.set(`${nodeId}:${i}`, outVarName);
      }
    }

    steps.push({
      nodeId,
      nodeType: 'func',
      funcName,
      variableName: varName,
      inputs: inputMap,
      outputs: outputMap,
      properties: {...node.properties},
    });
  }

  return steps;
}

function resolveInputExpression(
  node: LGraphNode,
  slotIndex: number,
  links: Record<number, any>,
  outputVarMap: Map<string, string>,
): string {
  const inp = node.inputs[slotIndex];
  if (!inp || inp.link === null || inp.link === undefined) return 'undefined';

  const link = links[inp.link];
  if (!link) return 'undefined';

  const key = `${link.origin_id}:${link.origin_slot}`;
  return outputVarMap.get(key) || 'undefined';
}

function compileUtilityNode(
  node: LGraphNode,
  links: Record<number, any>,
  outputVarMap: Map<string, string>,
  usedVarNames: Set<string>,
): CompiledStep | null {
  const title = node.title;
  const varName = uniqueVarName(toCamelCase(title), usedVarNames);
  usedVarNames.add(varName);

  const inputMap = new Map<string, string>();
  if (node.inputs) {
    for (let i = 0; i < node.inputs.length; i++) {
      const inp = node.inputs[i];
      if (node.isInputConnected(i))
        inputMap.set(inp.name, resolveInputExpression(node, i, links, outputVarMap));
    }
  }

  // Map outputs
  if (node.outputs) {
    for (let i = 0; i < node.outputs.length; i++)
      outputVarMap.set(`${node.id}:${i}`, varName);
  }

  return {
    nodeId: node.id,
    nodeType: 'utility',
    funcName: title,
    variableName: varName,
    inputs: inputMap,
    outputs: new Map([[node.outputs?.[0]?.name || 'value', varName]]),
    properties: {...node.properties},
  };
}

function toCamelCase(s: string): string {
  return s
    .replace(/[^a-zA-Z0-9]+/g, ' ')
    .trim()
    .split(' ')
    .map((word, i) => i === 0 ? word.toLowerCase() : word.charAt(0).toUpperCase() + word.slice(1).toLowerCase())
    .join('');
}

function uniqueVarName(base: string, used: Set<string>): string {
  if (!base) base = 'v';
  if (!used.has(base)) return base;
  let i = 2;
  while (used.has(`${base}${i}`)) i++;
  return `${base}${i}`;
}

function formatLiteral(value: any, type: string): string {
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
