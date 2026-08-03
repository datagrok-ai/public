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

/** Create a live FlowEditor mounted on a detached container. The data layer
 *  (nodes / connections) is populated synchronously by `addNode` /
 *  `addConnection`, so compiler / serializer tests can read it without waiting
 *  for the React render. Always pair with `destroyEditor` in a `finally`. */
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

/** Add a registered node type to the editor at a position. */
export async function addNode(flow: FlowEditor, typeName: string, x = 0, y = 0): Promise<FlowNode> {
  const node = createNode(typeName);
  if (!node) throw new Error(`Unknown node type: ${typeName}`);
  await flow.addNodeAt(node, x, y);
  return node;
}

/** Poll until `cond()` is true; resolves false on timeout. For assertions on
 *  asynchronously rendered DOM (React mounts, socket-position propagation). */
export async function until(cond: () => boolean, timeoutMs = 3000, stepMs = 50): Promise<boolean> {
  const deadline = Date.now() + timeoutMs;
  for (;;) {
    if (cond()) return true;
    if (Date.now() > deadline) return false;
    await new Promise((r) => setTimeout(r, stepMs));
  }
}

// ---------- BuiltGraph query helpers (for importer tests) ----------

/** All nodes whose underlying DG function name (case-insensitive) matches. */
export function nodesByFunc(graph: BuiltGraph, funcName: string): FlowNode[] {
  const lc = funcName.toLowerCase();
  return graph.nodes.filter((n) => (n.dgFunc?.name ?? '').toLowerCase() === lc);
}

/** Exactly one node for the given DG function — throws otherwise. */
export function oneNodeByFunc(graph: BuiltGraph, funcName: string): FlowNode {
  const matches = nodesByFunc(graph, funcName);
  if (matches.length !== 1)
    throw new Error(`Expected exactly one "${funcName}" node, found ${matches.length}`);
  return matches[0];
}

/** Nodes by their display label (e.g. constant nodes 'String', 'Boolean'). */
export function nodesByLabel(graph: BuiltGraph, label: string): FlowNode[] {
  return graph.nodes.filter((n) => n.label === label);
}

/** The source ref feeding a node's input, or null if unconnected. */
export function sourceOf(graph: BuiltGraph, target: FlowNode, targetKey: string):
  {node: FlowNode; key: string} | null {
  const c = graph.connections.find((x) => x.target === target && x.targetKey === targetKey);
  return c ? {node: c.source, key: c.sourceKey} : null;
}

/** Whether a connection exists from a node's output to a target's input. */
export function isConnected(
  graph: BuiltGraph, source: FlowNode, sourceKey: string, target: FlowNode, targetKey: string): boolean {
  return graph.connections.some((c) =>
    c.source === source && c.sourceKey === sourceKey && c.target === target && c.targetKey === targetKey);
}
