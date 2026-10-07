// The render tier as the page sees it: the blob's typed arrays, the ids and names, the type tree, and the
// adjacency built once so neighbours never need the server.

export async function loadAll() {
  const [manifest, schema, index, blob] = await Promise.all([
    getJson('api/manifest'), getJson('api/schema'), getJson('api/index'),
    fetch('api/graph.bin').then(ok).then((r) => r.arrayBuffer()),
  ]);
  const graph = parseBlob(blob);
  if (graph.nodes !== index.ids.length) throw new Error(`graph.bin has ${graph.nodes} nodes, index.json ${index.ids.length}: restart grok kg serve`);
  return {manifest, schema, ids: index.ids, names: index.names, ...graph, adjacency: buildAdjacency(graph.nodes, graph.from, graph.to)};
}

export function parseBlob(buffer) {
  const view = new DataView(buffer);
  const magic = String.fromCharCode(view.getUint8(0), view.getUint8(1), view.getUint8(2), view.getUint8(3));
  if (magic !== 'KGV1') throw new Error('graph.bin: not a graph blob');
  const length = view.getUint32(4, true);
  const header = JSON.parse(new TextDecoder().decode(new Uint8Array(buffer, 8, length)));
  const base = (8 + length + 3) & ~3;
  const section = (name) => {
    const s = header.sections[name];
    const ctor = s.dtype === 'u8' ? Uint8Array : s.dtype === 'u32' ? Uint32Array : Int32Array;
    return new ctor(buffer, base + s.offset, s.length);
  };
  return {
    batch: header.batch, nodes: header.nodes, edges: header.edges, dropped: header.dropped,
    nodeTypes: header.nodeTypes, edgeKinds: header.edgeKinds, enums: header.enums,
    type: section('type'), layer: section('layer'), visibility: section('visibility'), status: section('status'),
    provenance: section('provenance'), degree: section('degree'), parent: section('parent'),
    from: section('from'), to: section('to'), kind: section('kind'), confidence: section('confidence'), derivedBy: section('derivedBy'),
  };
}

/** CSR over both directions: for node i, entries offsets[i]..offsets[i+1] hold (neighbour, edge) pairs. */
export function buildAdjacency(n, from, to) {
  const offsets = new Uint32Array(n + 1);
  for (let e = 0; e < from.length; e++) {
    offsets[from[e] + 1]++;
    offsets[to[e] + 1]++;
  }
  for (let i = 0; i < n; i++) offsets[i + 1] += offsets[i];
  const fill = offsets.slice(0, n);
  const node = new Uint32Array(offsets[n]);
  const edge = new Uint32Array(offsets[n]);
  for (let e = 0; e < from.length; e++) {
    const a = from[e], b = to[e];
    node[fill[a]] = b; edge[fill[a]++] = e;
    node[fill[b]] = a; edge[fill[b]++] = e;
  }
  return {offsets, node, edge};
}

export const api = {
  node: (id) => getJson(`api/node?id=${encodeURIComponent(id)}`),
  edge: (from, to, kind) => getJson(`api/edge?from=${encodeURIComponent(from)}&to=${encodeURIComponent(to)}&kind=${encodeURIComponent(kind)}`),
  op: (name, target, limit = 50) => getJson(`api/op/${name}?target=${encodeURIComponent(target)}&limit=${limit}`),
  questions: () => getJson('api/questions'),
  ask: (id, params) => getJson(`api/ask/${encodeURIComponent(id)}?${new URLSearchParams(params).toString()}`),
  query: async (cypher, limit) => {
    const r = await fetch('api/query', {method: 'POST', headers: {'Content-Type': 'application/json'}, body: JSON.stringify({cypher, limit})});
    return (await ok(r)).json();
  },
};

async function getJson(url) {
  return (await ok(await fetch(url))).json();
}

async function ok(r) {
  if (r.ok) return r;
  let message = `${r.status} ${r.statusText}`;
  try {
    message = (await r.json()).error ?? message;
  }
  catch {
    // the body was not JSON; the status line is the message
  }
  throw new Error(message);
}
