/* eslint-disable max-len */
/** Pre-compilation validation. Runs against the FlowEditor data layer. */

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {topologicalSort} from './topological-sort';

export type ValidationSeverity = 'error' | 'warning';

export interface ValidationResult {
  severity: ValidationSeverity;
  message: string;
  nodeId?: string;
}

export function validateGraph(flow: FlowEditor): ValidationResult[] {
  const results: ValidationResult[] = [];
  const nodes = flow.getNodes();

  if (nodes.length === 0) {
    results.push({severity: 'warning', message: 'Graph is empty — nothing to generate'});
    return results;
  }

  try {
    topologicalSort(flow);
  } catch (e: unknown) {
    results.push({severity: 'error', message: (e as Error).message});
    return results;
  }

  const hasTableInput = nodes.some(
    (n) => n.dgNodeType === 'input' && n.dgOutputType === 'dataframe',
  );

  for (const node of nodes) {
    if (node.dgNodeType === 'input' && node.dgOutputType === 'column' && !hasTableInput) {
      results.push({
        severity: 'error',
        message: `Column input '${node.properties['paramName']}' requires a Table input in the graph`,
        nodeId: node.id,
      });
    }

    // Disconnected-node warning (skip pure inputs).
    if (node.dgNodeType !== 'input' && isDisconnected(node, flow))
      results.push({
        severity: 'warning',
        message: `Node '${node.label}' is disconnected from the graph`,
        nodeId: node.id,
      });

    if (node.dgNodeType === 'output' && !hasIncoming(node, flow))
      results.push({
        severity: 'warning',
        message: `Output node '${node.label}' has no incoming connection`,
        nodeId: node.id,
      });

    if (node.dgNodeType === 'input' || node.dgNodeType === 'output') {
      const paramName = String(node.properties['paramName'] ?? '');
      if (!paramName.trim())
        results.push({severity: 'error', message: `Node '${node.label}' has an empty parameter name`, nodeId: node.id});
      else if (!/^[a-zA-Z_][a-zA-Z0-9_]*$/.test(paramName))
        results.push({severity: 'error', message: `Node '${node.label}' has an invalid parameter name '${paramName}' (must be a valid JS identifier)`, nodeId: node.id});
    }
  }

  // Duplicate input/output param names.
  const seen = new Map<string, string>();
  for (const node of nodes) {
    if (node.dgNodeType !== 'input' && node.dgNodeType !== 'output') continue;
    const name = String(node.properties['paramName'] ?? '');
    if (!name) continue;
    if (seen.has(name))
      results.push({
        severity: 'error',
        message: `Duplicate parameter name '${name}' used by '${node.label}' and '${seen.get(name)}'`,
        nodeId: node.id,
      });
    else seen.set(name, node.label);
  }

  return results;
}

function isDisconnected(node: FlowNode, flow: FlowEditor): boolean {
  for (const c of flow.getConnections())
    if (c.source === node.id || c.target === node.id) return false;
  return true;
}

function hasIncoming(node: FlowNode, flow: FlowEditor): boolean {
  for (const c of flow.getConnections())
    if (c.target === node.id) return true;
  return false;
}
