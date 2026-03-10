/* eslint-disable max-len */
import {LGraph} from 'litegraph.js';
import {topologicalSort} from './topological-sort';
import {getGraphNodes} from './graph-utils';
import {FuncFlowNode} from '../types/funcflow-node';

export type ValidationSeverity = 'error' | 'warning';

export interface ValidationResult {
  severity: ValidationSeverity;
  message: string;
  nodeId?: number;
}

/** Validates the graph before compilation */
export function validateGraph(graph: LGraph): ValidationResult[] {
  const results: ValidationResult[] = [];
  const nodes = getGraphNodes(graph) as FuncFlowNode[];

  if (nodes.length === 0) {
    results.push({severity: 'warning', message: 'Graph is empty \u2014 nothing to generate'});
    return results;
  }

  try {
    topologicalSort(graph);
  } catch (e: any) {
    results.push({severity: 'error', message: e.message});
    return results;
  }

  const hasTableInput = nodes.some((n) => n.dgNodeType === 'input' && n.dgOutputType === 'dataframe');

  for (const node of nodes) {
    if (node.dgNodeType === 'input' && node.dgOutputType === 'column') {
      if (!hasTableInput) {
        results.push({
          severity: 'error',
          message: `Column input '${node.properties['paramName']}' requires a Table input in the graph`,
          nodeId: node.id,
        });
      }
    }

    if (node.dgFunc && node.inputs) {
      const func = node.dgFunc;
      for (let i = 0; i < node.inputs.length; i++) {
        const inp = node.inputs[i];
        const connected = node.isInputConnected(i);
        if (!connected) {
          const propKey = `_input_${inp.name}`;
          const storedVal = node.properties[propKey];
          // 0 and false are valid values; only undefined/null/empty-string mean "no value"
          const hasValue = storedVal !== undefined && storedVal !== null &&
            (typeof storedVal !== 'string' || storedVal !== '');

          // Check if the function parameter is nullable
          const funcParam = func.inputs.find((p: any) => p.name === inp.name);
          const isNullable = funcParam?.nullable === true || funcParam?.options?.optional === true || funcParam?.options?.nullable === true;

          if (!hasValue && !isNullable) {
            results.push({
              severity: 'error',
              message: `Required input '${inp.name}' on node '${node.title}' is not connected and has no value`,
              nodeId: node.id,
            });
          } else if (!hasValue)
            node.properties[propKey] = node.properties[propKey] ?? null;
        }
      }
    }

    if (node.inputs && node.outputs) {
      const anyInputConnected = node.inputs.some((_, idx) => node.isInputConnected(idx));
      const anyOutputConnected = node.outputs.some((_, idx) => node.isOutputConnected(idx));
      if (!anyInputConnected && !anyOutputConnected && node.dgNodeType !== 'input') {
        results.push({
          severity: 'warning',
          message: `Node '${node.title}' is disconnected from the graph`,
          nodeId: node.id,
        });
      }
    }

    if (node.dgNodeType === 'output') {
      const anyConnected = node.inputs && node.inputs.some((_, idx) => node.isInputConnected(idx));
      if (!anyConnected) {
        results.push({
          severity: 'warning',
          message: `Output node '${node.title}' has no incoming connection`,
          nodeId: node.id,
        });
      }
    }

    if (node.dgNodeType === 'input' || node.dgNodeType === 'output') {
      const paramName = node.properties['paramName'];
      if (!paramName || paramName.trim() === '') {
        results.push({
          severity: 'error',
          message: `Node '${node.title}' has an empty parameter name`,
          nodeId: node.id,
        });
      }
      if (paramName && !/^[a-zA-Z_][a-zA-Z0-9_]*$/.test(paramName)) {
        results.push({
          severity: 'error',
          message: `Node '${node.title}' has an invalid parameter name '${paramName}' (must be a valid JS identifier)`,
          nodeId: node.id,
        });
      }
    }
  }

  const paramNames = new Map<string, string>();
  for (const node of nodes) {
    if (node.dgNodeType === 'input' || node.dgNodeType === 'output') {
      const name = node.properties['paramName'];
      if (name && paramNames.has(name)) {
        results.push({
          severity: 'error',
          message: `Duplicate parameter name '${name}' used by '${node.title}' and '${paramNames.get(name)}'`,
          nodeId: node.id,
        });
      } else if (name)
        paramNames.set(name, node.title);
    }
  }

  return results;
}
