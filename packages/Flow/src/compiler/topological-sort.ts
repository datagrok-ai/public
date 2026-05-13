/** Kahn's algorithm topological sort over a Rete `NodeEditor`'s connections. */

import {FlowEditor} from '../rete/flow-editor';

export function topologicalSort(flow: FlowEditor): string[] {
  const nodes = flow.getNodes();
  if (nodes.length === 0) return [];

  const inDegree = new Map<string, number>();
  const outEdges = new Map<string, string[]>();
  for (const node of nodes) {
    inDegree.set(node.id, 0);
    outEdges.set(node.id, []);
  }

  for (const c of flow.getConnections()) {
    if (inDegree.has(c.target))
      inDegree.set(c.target, (inDegree.get(c.target) ?? 0) + 1);
    if (outEdges.has(c.source))
      outEdges.get(c.source)!.push(c.target);
  }

  const queue: string[] = [];
  for (const [id, deg] of inDegree)
    if (deg === 0) queue.push(id);

  const sorted: string[] = [];
  while (queue.length > 0) {
    const id = queue.shift()!;
    sorted.push(id);
    for (const t of outEdges.get(id) ?? []) {
      const next = (inDegree.get(t) ?? 1) - 1;
      inDegree.set(t, next);
      if (next === 0) queue.push(t);
    }
  }

  if (sorted.length !== nodes.length) {
    const remaining = nodes.filter((n) => !sorted.includes(n.id)).map((n) => n.label || n.id);
    throw new Error(`Cycle detected in graph involving nodes: ${remaining.join(', ')}`);
  }
  return sorted;
}
