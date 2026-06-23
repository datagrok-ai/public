/** Layered, banded graph layout — extracted from the creation-script importer
 *  so the same arrangement can be applied to any live graph (the ribbon's
 *  "Clean Layout" action) as well as to freshly-imported scripts.
 *
 *  The algorithm:
 *   - **Layers** (columns): `layer(n) = 0` for sources, else `max(layer(pred))+1`
 *     (longest path). Every edge therefore points strictly rightward. Column x
 *     is global — shared per layer, width = widest node in that layer.
 *   - **Bands** (rows): one horizontal band per weakly connected component,
 *     stacked in **dependency order** (a path producing a table via `SetVar`
 *     sits above the path that reads it through a `Select Table` node — matched
 *     by normalized name), ties broken by node order. Because the execution
 *     topological sort ranks components by topmost `y`, this band order *is* the
 *     execution order.
 *   - Within a band/column, nodes order by predecessor barycenter and greedily
 *     stack so chains read as straight lanes, branches fan out, nothing overlaps.
 *
 *  Pure and DOM-free: it only reads node metadata and mutates `node.pos`. */

import {FlowNode} from './scheme';

/** A directed layout edge between two nodes (source → target). */
export interface LayoutEdge {
  source: FlowNode;
  target: FlowNode;
}

const SELECT_TABLE_TYPE = 'Utilities/Select Table';

/** Longest-path layering over a DAG: `layer(n) = 0` for nodes with no incoming
 *  edge, else `max(layer(pred)) + 1`. Nodes left unlayered (cycles) fall to 0. */
export function computeLayers(nodes: FlowNode[], edges: LayoutEdge[]): Map<FlowNode, number> {
  const outgoing = new Map<FlowNode, FlowNode[]>();
  const indegree = new Map<FlowNode, number>();
  for (const node of nodes) {
    outgoing.set(node, []);
    indegree.set(node, 0);
  }
  for (const edge of edges) {
    if (!outgoing.has(edge.source) || !indegree.has(edge.target)) continue;
    outgoing.get(edge.source)!.push(edge.target);
    indegree.set(edge.target, (indegree.get(edge.target) ?? 0) + 1);
  }

  const layer = new Map<FlowNode, number>();
  const remaining = new Map(indegree);
  const queue: FlowNode[] = [];
  for (const node of nodes) {
    if ((indegree.get(node) ?? 0) === 0) {
      layer.set(node, 0);
      queue.push(node);
    }
  }
  while (queue.length > 0) {
    const node = queue.shift()!;
    const here = layer.get(node) ?? 0;
    for (const next of outgoing.get(node) ?? []) {
      layer.set(next, Math.max(layer.get(next) ?? 0, here + 1));
      const left = (remaining.get(next) ?? 1) - 1;
      remaining.set(next, left);
      if (left === 0) queue.push(next);
    }
  }
  for (const node of nodes) if (!layer.has(node)) layer.set(node, 0);
  return layer;
}

/** Arrange `nodes` in place (mutating `node.pos`) using the layered/banded
 *  layout, given a precomputed `layer` per node (see `computeLayers`, or the
 *  importer's incremental layering). */
export function layoutGraph(nodes: FlowNode[], edges: LayoutEdge[], layer: Map<FlowNode, number>): void {
  const marginX = 40;
  const marginY = 40;
  const columnGap = 60;
  const rowGap = 24;
  const bandGap = 56;

  // Global column x positions from the widest node per layer.
  const columnWidth = new Map<number, number>();
  for (const node of nodes) {
    const l = layer.get(node) ?? 0;
    columnWidth.set(l, Math.max(columnWidth.get(l) ?? 0, estimateNodeWidth(node)));
  }
  const columnX = new Map<number, number>();
  let x = marginX;
  for (const l of Array.from(columnWidth.keys()).sort((a, b) => a - b)) {
    columnX.set(l, x);
    x += columnWidth.get(l)! + columnGap;
  }

  const predecessors = new Map<FlowNode, FlowNode[]>();
  for (const c of edges) {
    let list = predecessors.get(c.target);
    if (!list) predecessors.set(c.target, list = []);
    list.push(c.source);
  }

  const centerY = new Map<FlowNode, number>();
  let bandTop = marginY;
  for (const component of orderedComponents(nodes, edges)) {
    const byLayer = new Map<number, FlowNode[]>();
    for (const node of component) {
      const l = layer.get(node) ?? 0;
      let bucket = byLayer.get(l);
      if (!bucket) byLayer.set(l, bucket = []);
      bucket.push(node);
    }

    let bandBottom = bandTop;
    for (const l of Array.from(byLayer.keys()).sort((a, b) => a - b)) {
      const keyed = byLayer.get(l)!.map((node, index) => {
        const centers = (predecessors.get(node) ?? [])
          .map((p) => centerY.get(p))
          .filter((y): y is number => y !== undefined);
        const barycenter = centers.length > 0 ?
          centers.reduce((a, b) => a + b, 0) / centers.length :
          Number.POSITIVE_INFINITY;
        return {node, index, barycenter};
      });
      keyed.sort((a, b) => a.barycenter === b.barycenter ? a.index - b.index : a.barycenter - b.barycenter);

      let nextFreeY = bandTop;
      for (const {node, barycenter} of keyed) {
        const height = estimateNodeHeight(node);
        const desiredTop = Number.isFinite(barycenter) ? barycenter - height / 2 : nextFreeY;
        const top = Math.max(nextFreeY, desiredTop, bandTop);
        node.pos = {x: columnX.get(layer.get(node) ?? 0)!, y: top};
        centerY.set(node, top + height / 2);
        nextFreeY = top + height + rowGap;
        bandBottom = Math.max(bandBottom, top + height);
      }
    }
    bandTop = bandBottom + bandGap;
  }
}

/** Weakly connected components, each a node list in stable (input) order,
 *  ordered top-to-bottom for layout: a component that produces a table (a
 *  `SetVar` variable) precedes the component that reads it via a `Select Table`
 *  node; ties break by earliest node index. */
export function orderedComponents(nodes: FlowNode[], edges: LayoutEdge[]): FlowNode[][] {
  const indexOf = new Map<FlowNode, number>();
  nodes.forEach((n, i) => indexOf.set(n, i));

  const adjacency = new Map<FlowNode, FlowNode[]>();
  for (const node of nodes) adjacency.set(node, []);
  for (const c of edges) {
    adjacency.get(c.source)?.push(c.target);
    adjacency.get(c.target)?.push(c.source);
  }

  const componentOf = new Map<FlowNode, number>();
  const components: FlowNode[][] = [];
  for (const start of nodes) {
    if (componentOf.has(start)) continue;
    const id = components.length;
    const members: FlowNode[] = [];
    const stack = [start];
    componentOf.set(start, id);
    while (stack.length > 0) {
      const node = stack.pop()!;
      members.push(node);
      for (const next of adjacency.get(node) ?? []) {
        if (!componentOf.has(next)) {
          componentOf.set(next, id);
          stack.push(next);
        }
      }
    }
    members.sort((a, b) => indexOf.get(a)! - indexOf.get(b)!);
    components.push(members);
  }

  // Producer → consumer edges between components, matched by normalized name
  // (a Select Table's tableName ↔ a SetVar's variableName).
  const norm = (s: unknown): string => String(s ?? '').toLowerCase().replace(/[^a-z0-9]/g, '');
  const producerOf = new Map<string, number>();
  components.forEach((members, id) => {
    for (const node of members) {
      if (node.dgFunc?.name?.toLowerCase() === 'setvar') {
        const key = norm(node.inputValues['variableName']);
        if (key !== '' && !producerOf.has(key)) producerOf.set(key, id);
      }
    }
  });

  const inDegree = components.map(() => 0);
  const downstream = components.map(() => new Set<number>());
  components.forEach((members, id) => {
    for (const node of members) {
      if (node.dgTypeName !== SELECT_TABLE_TYPE) continue;
      const producer = producerOf.get(norm(node.properties['tableName']));
      if (producer !== undefined && producer !== id && !downstream[producer].has(id)) {
        downstream[producer].add(id);
        inDegree[id]++;
      }
    }
  });

  // Kahn over components; among ready ones, pick the earliest by node index.
  const earliest = components.map((m) => (m.length > 0 ? indexOf.get(m[0])! : 0));
  const ready: number[] = [];
  for (let id = 0; id < components.length; id++)
    if (inDegree[id] === 0) ready.push(id);
  const order: number[] = [];
  const placed = new Set<number>();
  while (ready.length > 0) {
    let best = 0;
    for (let i = 1; i < ready.length; i++)
      if (earliest[ready[i]] < earliest[ready[best]]) best = i;
    const id = ready.splice(best, 1)[0];
    order.push(id);
    placed.add(id);
    for (const next of downstream[id])
      if (--inDegree[next] === 0) ready.push(next);
  }
  for (let id = 0; id < components.length; id++)
    if (!placed.has(id)) order.push(id);

  return order.map((id) => components[id]);
}

/** Estimated rendered height in canvas units, used for stacking before the DOM
 *  exists. Collapsed nodes render as a bare title bar. Expanded: title ≈ 28px,
 *  description ≈ 22px, each socket row ≈ 20px, body padding 12px (see
 *  funcflow.css). Input and output columns sit side by side, so rows = max of
 *  the two. */
export function estimateNodeHeight(node: FlowNode): number {
  if (node.collapsed) return 30;
  const rows = Math.max(Object.keys(node.inputs).length, Object.keys(node.outputs).length, 1);
  return 28 + (node.description ? 22 : 0) + 12 + rows * 20;
}

/** Estimated rendered width. Collapsed nodes are title-driven (CSS min-width
 *  160px, ≈6.5px/char at the 12px title font plus status dot and paddings);
 *  expanded nodes are dominated by their socket-label rows. */
export function estimateNodeWidth(node: FlowNode): number {
  const labelWidth = 44 + String(node.label ?? '').length * 6.5;
  return Math.max(node.collapsed ? 160 : 220, labelWidth);
}
