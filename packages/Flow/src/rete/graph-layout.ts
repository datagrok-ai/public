/** Layered, banded graph layout shared by the creation-script importer and the ribbon's
 *  "Clean Layout". Pure and DOM-free: reads node metadata, mutates `node.pos`. */

import {FlowNode, isExecKey, hiddenSocketRow} from './scheme';

export interface LayoutEdge {
  source: FlowNode;
  target: FlowNode;
}

const SELECT_TABLE_TYPE = 'Utilities/Select Table';

/** Longest-path layering over a DAG; nodes left unlayered (cycles) fall to 0. */
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

/** Arranges nodes in place (mutating `node.pos`), given a precomputed layer per node. */
export function layoutGraph(nodes: FlowNode[], edges: LayoutEdge[], layer: Map<FlowNode, number>): void {
  const marginX = 40;
  const marginY = 40;
  const columnGap = 60;
  const rowGap = 24;
  const bandGap = 56;

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

/** Weakly connected components in layout order: a table-producing component precedes
 *  the one that reads it via a Select Table node; ties break by earliest node index. */
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

  // Producer → consumer edges: a Select Table's tableName ↔ a SetVar's variableName, normalized.
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

/** Estimated rendered height in canvas units (pre-DOM); constants mirror funcflow.css. */
export function estimateNodeHeight(node: FlowNode): number {
  if (node.collapsed) return 30;
  // Hidden rows take no space; import-time layout runs pre-bridge, so unwired.
  const conn = (side: 'input' | 'output', key: string): boolean =>
    node.editorBridge?.isSocketConnected(node.id, side, key) ?? false;
  const visible = (side: 'input' | 'output', keys: string[]): number =>
    keys.filter((k) => !isExecKey(k) && !hiddenSocketRow(node, side, k, conn)).length;
  const rows = Math.max(visible('input', Object.keys(node.inputs)),
    visible('output', Object.keys(node.outputs)), 1);
  // title bar + the always-present one-line info row + body padding + rows.
  return 28 + 19 + 12 + rows * 20;
}

/** Estimated rendered width; constants mirror the funcflow.css min-widths and title font. */
export function estimateNodeWidth(node: FlowNode): number {
  const labelWidth = 44 + String(node.label ?? '').length * 6.5;
  // Must cap with the CSS max-width 280 or long-titled nodes get phantom column width.
  return Math.min(280, Math.max(node.collapsed ? 160 : 220, labelWidth));
}
