import {LGraph} from 'litegraph.js';
import {getGraphNodes} from './graph-utils';

/** Performs topological sort on the graph using Kahn's algorithm.
 *  Returns ordered node IDs or throws if a cycle is detected. */
export function topologicalSort(graph: LGraph): number[] {
  const nodes = getGraphNodes(graph);
  if (nodes.length === 0) return [];

  const links = graph.links || {};

  // Build adjacency: inDegree for each node, and outgoing edges
  const inDegree = new Map<number, number>();
  const outEdges = new Map<number, number[]>();

  for (const node of nodes) {
    inDegree.set(node.id, 0);
    outEdges.set(node.id, []);
  }

  // Count incoming links per node
  for (const [, link] of Object.entries(links)) {
    if (!link) continue;
    const targetId = link.target_id;
    const originId = link.origin_id;
    if (inDegree.has(targetId))
      inDegree.set(targetId, (inDegree.get(targetId) || 0) + 1);

    if (outEdges.has(originId))
      outEdges.get(originId)!.push(targetId);
  }

  // Start with all nodes that have no incoming edges
  const queue: number[] = [];
  for (const [nodeId, degree] of inDegree.entries())
    if (degree === 0) queue.push(nodeId);


  const sorted: number[] = [];
  while (queue.length > 0) {
    const nodeId = queue.shift()!;
    sorted.push(nodeId);

    const targets = outEdges.get(nodeId) || [];
    for (const targetId of targets) {
      const newDegree = (inDegree.get(targetId) || 1) - 1;
      inDegree.set(targetId, newDegree);
      if (newDegree === 0) queue.push(targetId);
    }
  }

  if (sorted.length !== nodes.length) {
    const remaining = nodes.filter((n) => !sorted.includes(n.id)).map((n) => n.title || n.id);
    throw new Error(`Cycle detected in graph involving nodes: ${remaining.join(', ')}`);
  }

  return sorted;
}
