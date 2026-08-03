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

/** The minimal node/edge shape the sort needs — so it can run on the live editor
 *  *and* on plain lists (e.g. the flow summary), guaranteeing one canonical order. */
export interface SortNode {
  id: string;
  pos?: {x: number; y: number};
  label?: string;
}
export interface SortEdge {
  source: string;
  target: string;
}

/** The execution order, the editor way (strict: throws on a cycle). */
export function topologicalSort(flow: FlowEditor): string[] {
  const nodes = flow.getNodes();
  const sorted = topologicalSortNodes(nodes, flow.getConnections());
  if (sorted.length !== nodes.length) {
    const done = new Set(sorted);
    const remaining = nodes.filter((n) => !done.has(n.id)).map((n) => n.label || n.id);
    throw new Error(`Cycle detected in graph involving nodes: ${remaining.join(', ')}`);
  }
  return sorted;
}

/** Pure core: the same deterministic, layout-aware ordering on plain node/edge
 *  lists. On a cycle it returns the acyclic prefix (shorter than `nodes`) rather
 *  than throwing — callers choose strict vs lenient handling. */
export function topologicalSortNodes(nodes: SortNode[], connections: SortEdge[]): string[] {
  if (nodes.length === 0) return [];

  const posY = (n: SortNode): number => n.pos?.y ?? 0;
  const posX = (n: SortNode): number => n.pos?.x ?? 0;

  const insertionIndex = new Map<string, number>();
  nodes.forEach((n, i) => insertionIndex.set(n.id, i));

  // ---- weakly connected components, ranked by their topmost node ----
  const adjacency = new Map<string, string[]>();
  for (const n of nodes) adjacency.set(n.id, []);
  for (const c of connections) {
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
    if (posY(n) < top.y || (posY(n) === top.y && posX(n) < top.x)) {
      top.y = posY(n);
      top.x = posX(n);
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
  for (const c of connections) {
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
    if (posY(nodeA) !== posY(nodeB)) return posY(nodeA) < posY(nodeB);
    if (posX(nodeA) !== posX(nodeB)) return posX(nodeA) < posX(nodeB);
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

  // Lenient: on a cycle, `sorted` is the acyclic prefix. The strict wrapper
  // (topologicalSort) turns that into a thrown error; the summary appends the
  // remainder instead.
  return sorted;
}
