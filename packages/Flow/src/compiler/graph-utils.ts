import {LGraph, LGraphNode} from 'litegraph.js';

/** Safely get all nodes from a graph. LGraph._nodes is private, so we use serialize(). */
export function getGraphNodes(graph: LGraph): LGraphNode[] {
  const data = graph.serialize();
  if (!data || !data.nodes) return [];
  return data.nodes
    .map((n: {id: number}) => graph.getNodeById(n.id))
    .filter((n: LGraphNode | undefined): n is LGraphNode => n != null);
}
