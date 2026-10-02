// The filter pipeline: the left pane's state becomes two masks over the blob, then the compact arrays the
// renderer draws. Everything here is typed arrays over indices; ids only appear at the edges of the page.

export const PRESETS = {
  spine: {label: 'Spine', types: ['feature', 'concept', 'initiative', 'package', 'library', 'person', 'team', 'customer', 'release']},
  code: {label: 'Code', types: ['feature', 'package', 'library', 'source-file', 'declaration', 'function', 'app', 'cell-renderer', 'editor',
    'file-handler', 'file-viewer', 'filter', 'lifecycle-hook', 'panel', 'query', 'script', 'script-handler', 'sem-type-detector', 'viewer',
    'semantic-type', 'connection', 'container', 'script-environment', 'endpoint', 'db-table', 'migration']},
  docs: {label: 'Docs', types: ['feature', 'concept', 'doc-page', 'web-page', 'doc-anchor', 'media', 'sample', 'tutorial']},
  tests: {label: 'Tests', types: ['feature', 'scenario', 'test', 'test-suite', 'sample']},
  work: {label: 'Work', types: ['feature', 'initiative', 'ticket', 'pull-request', 'commit', 'release', 'report', 'person', 'customer']},
  everything: {label: 'Everything', types: null, confirm: 60000},
};

export class Filter {
  constructor(data) {
    this.data = data;
    this.nodeTypes = new Set();
    this.edgeKinds = new Set(data.edgeKinds);
    this.facets = {layer: new Set(data.enums.layer), visibility: new Set(data.enums.visibility), status: new Set(data.enums.status),
      provenance: new Set(data.enums.provenance), derivedBy: new Set(data.enums.derivedBy)};
    this.minConfidence = 0;
    this.scope = null;
    this.pinned = new Set();
    this.restrict = null;
    /** Explore mode: `open` nodes are shown as they are, `expanded` ones bring their neighbours along. */
    this.explore = null;
  }

  snapshot() {
    return {nodeTypes: new Set(this.nodeTypes), edgeKinds: new Set(this.edgeKinds),
      facets: Object.fromEntries(Object.entries(this.facets).map(([k, v]) => [k, new Set(v)])),
      minConfidence: this.minConfidence, scope: this.scope, pinned: new Set(this.pinned), restrict: this.restrict,
      explore: this.explore ? {open: new Set(this.explore.open), expanded: new Set(this.explore.expanded)} : null};
  }

  restore(s) {
    this.nodeTypes = s.nodeTypes;
    this.edgeKinds = s.edgeKinds;
    this.facets = s.facets;
    this.minConfidence = s.minConfidence;
    this.scope = s.scope;
    this.pinned = s.pinned;
    this.restrict = s.restrict;
    this.explore = s.explore;
  }

  /** Both masks and the counts the panes show. */
  apply() {
    const d = this.data;
    const n = d.nodes;
    const typeOn = mask(d.nodeTypes, this.nodeTypes);
    const layerOn = mask(d.enums.layer, this.facets.layer);
    const visOn = mask(d.enums.visibility, this.facets.visibility);
    const statusOn = mask(d.enums.status, this.facets.status);
    const provOn = mask(d.enums.provenance, this.facets.provenance);
    const kindOn = mask(d.edgeKinds, this.edgeKinds);
    const derivedOn = mask(d.enums.derivedBy, this.facets.derivedBy);
    const edgeOn = new Uint8Array(d.edges);
    for (let e = 0; e < d.edges; e++)
      edgeOn[e] = kindOn[d.kind[e]] && derivedOn[d.derivedBy[e]] && d.confidence[e] >= this.minConfidence ? 1 : 0;
    const nodeOn = new Uint8Array(n);
    for (let i = 0; i < n; i++)
      nodeOn[i] = typeOn[d.type[i]] && layerOn[d.layer[i]] && visOn[d.visibility[i]] && statusOn[d.status[i]] && provOn[d.provenance[i]] ? 1 : 0;
    if (this.explore) {
      const keep = this.exploreMask(nodeOn, edgeOn);
      for (let i = 0; i < n; i++) nodeOn[i] = keep[i];
    }
    if (this.scope) {
      const inScope = this.scopeMask(edgeOn);
      for (let i = 0; i < n; i++) nodeOn[i] &= inScope[i];
    }
    if (this.restrict) {
      for (let i = 0; i < n; i++)
        if (!this.restrict.has(i)) nodeOn[i] = 0;
    }
    for (const i of this.pinned) nodeOn[i] = 1;
    let nodeCount = 0;
    const typeCounts = new Uint32Array(d.nodeTypes.length);
    for (let i = 0; i < n; i++)
      if (nodeOn[i]) {
        nodeCount++;
        typeCounts[d.type[i]]++;
      }
    let edgeCount = 0;
    const kindCounts = new Uint32Array(d.edgeKinds.length);
    for (let e = 0; e < d.edges; e++) {
      if (edgeOn[e] && nodeOn[d.from[e]] && nodeOn[d.to[e]]) {
        edgeCount++;
        kindCounts[d.kind[e]]++;
      }
      else edgeOn[e] = 0;
    }
    return {nodeOn, edgeOn, nodeCount, edgeCount, typeCounts, kindCounts};
  }

  /** Open nodes, expanded nodes, and the neighbours of the expanded ones that pass the filters. */
  exploreMask(nodeOn, edgeOn) {
    const d = this.data;
    const keep = new Uint8Array(d.nodes);
    const {offsets, node, edge} = d.adjacency;
    for (const i of this.explore.open) keep[i] = 1;
    for (const i of this.explore.expanded) {
      keep[i] = 1;
      for (let k = offsets[i]; k < offsets[i + 1]; k++)
        if (edgeOn[edge[k]] && nodeOn[node[k]]) keep[node[k]] = 1;
    }
    return keep;
  }

  /** The scope node, everything under it by part-of, and whatever the enabled edges reach within `hops`. */
  scopeMask(edgeOn) {
    const d = this.data;
    const on = new Uint8Array(d.nodes);
    const root = this.scope.index;
    on[root] = 1;
    let frontier = [root];
    for (let i = 0; i < d.nodes; i++) {
      for (let p = d.parent[i]; p >= 0; p = d.parent[p])
        if (p === root) {
          on[i] = 1;
          frontier.push(i);
          break;
        }
    }
    const {offsets, node, edge} = d.adjacency;
    for (let hop = 0; hop < this.scope.hops; hop++) {
      const next = [];
      for (const i of frontier)
        for (let k = offsets[i]; k < offsets[i + 1]; k++) {
          const j = node[k];
          if (on[j] || !edgeOn[edge[k]]) continue;
          on[j] = 1;
          next.push(j);
        }
      frontier = next;
    }
    return on;
  }

  /** Indices of every neighbour of [index] over every edge kind, for Expand. */
  neighbours(index) {
    const {offsets, node} = this.data.adjacency;
    const out = new Set();
    for (let k = offsets[index]; k < offsets[index + 1]; k++) out.add(node[k]);
    return out;
  }
}

function mask(names, enabled) {
  const m = new Uint8Array(names.length);
  names.forEach((name, i) => m[i] = enabled.has(name) ? 1 : 0);
  return m;
}

/** The arrays the renderer takes: visible nodes renumbered densely, links as pairs of new indices. */
export function compact(data, masks, style, positions) {
  const {nodeOn, edgeOn, nodeCount, edgeCount} = masks;
  const map = new Int32Array(data.nodes).fill(-1);
  const back = new Uint32Array(nodeCount);
  const pos = new Float32Array(nodeCount * 2);
  const colors = new Float32Array(nodeCount * 4);
  const sizes = new Float32Array(nodeCount);
  const clusters = new Array(nodeCount);
  let k = 0;
  for (let i = 0; i < data.nodes; i++) {
    if (!nodeOn[i]) continue;
    map[i] = k;
    back[k] = i;
    pos[k * 2] = positions[i * 2];
    pos[k * 2 + 1] = positions[i * 2 + 1];
    const c = style.nodeColor[data.type[i]];
    colors[k * 4] = c[0]; colors[k * 4 + 1] = c[1]; colors[k * 4 + 2] = c[2]; colors[k * 4 + 3] = 1;
    sizes[k] = style.size(data.degree[i]);
    clusters[k] = style.cluster[data.type[i]];
    k++;
  }
  const links = new Float32Array(edgeCount * 2);
  const linkColors = new Float32Array(edgeCount * 4);
  const linkEdge = new Uint32Array(edgeCount);
  let m = 0;
  for (let e = 0; e < data.edges; e++) {
    if (!edgeOn[e]) continue;
    links[m * 2] = map[data.from[e]];
    links[m * 2 + 1] = map[data.to[e]];
    const c = style.edgeColor[data.kind[e]];
    linkColors[m * 4] = c[0]; linkColors[m * 4 + 1] = c[1]; linkColors[m * 4 + 2] = c[2];
    linkColors[m * 4 + 3] = 0.12 + 0.4 * data.confidence[e] / 100;
    linkEdge[m] = e;
    m++;
  }
  return {map, back, positions: pos, colors, sizes, clusters, links, linkColors, linkEdge, nodeCount, edgeCount};
}
