/** Shared helpers for the Flow test suites. */
import * as ui from 'datagrok-api/ui';

import {FlowEditor} from '../rete/flow-editor';
import {FlowNode} from '../rete/scheme';
import {createNode} from '../rete/node-factory';
import {BuiltGraph} from '../import/creation-script-importer';

export interface TestEditor {
  flow: FlowEditor;
  container: HTMLElement;
}

/** Live FlowEditor on a detached container; the data layer is populated synchronously,
 *  so tests can read it without waiting for the React render. Pair with `destroyEditor`. */
export function makeEditor(): TestEditor {
  const container = ui.div([], {style: {width: '1000px', height: '700px', position: 'absolute', left: '-10000px'}});
  document.body.appendChild(container);
  const flow = new FlowEditor(container);
  return {flow, container};
}

export function destroyEditor(e: TestEditor): void {
  try {
    e.flow.destroy();
  } finally {
    e.container.remove();
  }
}

export async function addNode(flow: FlowEditor, typeName: string, x = 0, y = 0): Promise<FlowNode> {
  const node = createNode(typeName);
  if (!node) throw new Error(`Unknown node type: ${typeName}`);
  await flow.addNodeAt(node, x, y);
  return node;
}

/** Resolves false on timeout — assert the returned boolean. */
export async function until(cond: () => boolean, timeoutMs = 3000, stepMs = 50): Promise<boolean> {
  const deadline = Date.now() + timeoutMs;
  for (;;) {
    if (cond()) return true;
    if (Date.now() > deadline) return false;
    await new Promise((r) => setTimeout(r, stepMs));
  }
}

export function nodesByFunc(graph: BuiltGraph, funcName: string): FlowNode[] {
  const lc = funcName.toLowerCase();
  return graph.nodes.filter((n) => (n.dgFunc?.name ?? '').toLowerCase() === lc);
}

export function oneNodeByFunc(graph: BuiltGraph, funcName: string): FlowNode {
  const matches = nodesByFunc(graph, funcName);
  if (matches.length !== 1)
    throw new Error(`Expected exactly one "${funcName}" node, found ${matches.length}`);
  return matches[0];
}

export function nodesByLabel(graph: BuiltGraph, label: string): FlowNode[] {
  return graph.nodes.filter((n) => n.label === label);
}

export function sourceOf(graph: BuiltGraph, target: FlowNode, targetKey: string):
  {node: FlowNode; key: string} | null {
  const c = graph.connections.find((x) => x.target === target && x.targetKey === targetKey);
  return c ? {node: c.source, key: c.sourceKey} : null;
}

export function isConnected(
  graph: BuiltGraph, source: FlowNode, sourceKey: string, target: FlowNode, targetKey: string): boolean {
  return graph.connections.some((c) =>
    c.source === source && c.sourceKey === sourceKey && c.target === target && c.targetKey === targetKey);
}
