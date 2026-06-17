/** Kahn's algorithm topological sort over a Rete `NodeEditor`'s connections.
 *
 * Deterministic and layout-aware:
 *  - **Disjoint subgraphs** (weakly connected components) execute one after
 *    another, ordered by their topmost node on the canvas — a path placed
 *    above another finishes completely before the lower one starts. Data can
 *    flow between components implicitly (e.g. a Select Table node reading a
 *    table that an upper path opened), so this order is part of the contract,
 *    not just cosmetics.
 *  - **Within a component**, ready nodes are processed top-to-bottom
 *    (y, then x, then insertion order), so the emitted script is stable
 *    across runs and follows the visual reading order. */

import {FlowEditor} from '../rete/flow-editor';

export function topologicalSort(flow: FlowEditor): string[] {
  const nodes = flow.getNodes();
  if (nodes.length === 0) return [];

  const insertionIndex = new Map<string, number>();
  nodes.forEach((n, i) => insertionIndex.set(n.id, i));

  // ---- weakly connected components, ranked by their topmost node ----
  const adjacency = new Map<string, string[]>();
  for (const n of nodes) adjacency.set(n.id, []);
  for (const c of flow.getConnections()) {
    if (adjacency.has(c.source) && adjacency.has(c.target)) {
      adjacency.get(c.source)!.push(c.target);
      adjacency.get(c.target)!.push(c.source);
    }
  }

  const componentOf = new Map<string, number>();
  let componentCount = 0;
  for (const n of nodes) {
    if (componentOf.has(n.id)) continue;
    const component = componentCount++;
    const queue = [n.id];
    componentOf.set(n.id, component);
    while (queue.length > 0) {
      const current = queue.pop()!;
      for (const next of adjacency.get(current) ?? []) {
        if (!componentOf.has(next)) {
          componentOf.set(next, component);
          queue.push(next);
        }
      }
    }
  }

  const tops: Array<{y: number; x: number}> =
    Array.from({length: componentCount}, () => ({y: Number.POSITIVE_INFINITY, x: Number.POSITIVE_INFINITY}));
  for (const n of nodes) {
    const top = tops[componentOf.get(n.id)!];
    if (n.pos.y < top.y || (n.pos.y === top.y && n.pos.x < top.x)) {
      top.y = n.pos.y;
      top.x = n.pos.x;
    }
  }
  const componentRank = new Array<number>(componentCount);
  Array.from({length: componentCount}, (_, i) => i)
    .sort((a, b) => (tops[a].y - tops[b].y) || (tops[a].x - tops[b].x))
    .forEach((component, rank) => {
      componentRank[component] = rank;
    });

  // ---- Kahn's with a (componentRank, y, x, insertion) priority pick ----
  // A DAG component always has a ready node while it has unprocessed nodes,
  // so preferring the lowest component rank drains each component fully
  // before the next one starts.
  const inDegree = new Map<string, number>();
  const outEdges = new Map<string, string[]>();
  const byId = new Map(nodes.map((n) => [n.id, n]));
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

  const precedes = (a: string, b: string): boolean => {
    const rankA = componentRank[componentOf.get(a)!];
    const rankB = componentRank[componentOf.get(b)!];
    if (rankA !== rankB) return rankA < rankB;
    const nodeA = byId.get(a)!;
    const nodeB = byId.get(b)!;
    if (nodeA.pos.y !== nodeB.pos.y) return nodeA.pos.y < nodeB.pos.y;
    if (nodeA.pos.x !== nodeB.pos.x) return nodeA.pos.x < nodeB.pos.x;
    return (insertionIndex.get(a) ?? 0) < (insertionIndex.get(b) ?? 0);
  };

  const ready: string[] = [];
  for (const [id, degree] of inDegree)
    if (degree === 0) ready.push(id);

  const sorted: string[] = [];
  while (ready.length > 0) {
    let best = 0;
    for (let i = 1; i < ready.length; i++)
      if (precedes(ready[i], ready[best])) best = i;
    const id = ready.splice(best, 1)[0];
    sorted.push(id);
    for (const target of outEdges.get(id) ?? []) {
      const next = (inDegree.get(target) ?? 1) - 1;
      inDegree.set(target, next);
      if (next === 0) ready.push(target);
    }
  }

  if (sorted.length !== nodes.length) {
    const remaining = nodes.filter((n) => !sorted.includes(n.id)).map((n) => n.label || n.id);
    throw new Error(`Cycle detected in graph involving nodes: ${remaining.join(', ')}`);
  }
  return sorted;
}
