"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.batchId = batchId;
exports.gitRevisions = gitRevisions;
exports.projectPublic = projectPublic;
exports.toolsVersion = toolsVersion;
exports.writeBuild = writeBuild;
var fs = _interopRequireWildcard(require("fs"));
var path = _interopRequireWildcard(require("path"));
var _crypto = require("crypto");
var _child_process = require("child_process");
var _types = require("../types");
function _interopRequireWildcard(e, t) { if ("function" == typeof WeakMap) var r = new WeakMap(), n = new WeakMap(); return (_interopRequireWildcard = function (e, t) { if (!t && e && e.__esModule) return e; var o, i, f = { __proto__: null, default: e }; if (null === e || "object" != typeof e && "function" != typeof e) return f; if (o = t ? n : r) { if (o.has(e)) return o.get(e); o.set(e, f); } for (const t in e) "default" !== t && {}.hasOwnProperty.call(e, t) && ((i = (o = Object.defineProperty) && Object.getOwnPropertyDescriptor(e, t)) && (i.get || i.set) ? o(f, t, i) : f[t] = e[t]); return f; })(e, t); }
/// The build outputs (build-plan.md WO-1, Decisions): `data/nodes/<type>.jsonl`, `data/edges/<type|refname>.jsonl`,
/// `manifest.json` and `reports/`, deterministic line by line; the public projection of a graph.

/** Node types the public snapshot carries (Decisions "Public mode"); a doc-page only when it is a help page. */
const PUBLIC_TYPES = ['feature', 'concept', 'package', 'library', 'doc-page', 'doc-anchor', 'sample', 'scenario'];
const HEAD_KEYS = ['id', 'type', 'name', 'from', 'to'];
const INVALID_CAP = 500;
function toolsVersion() {
  return JSON.parse(fs.readFileSync(path.join(__dirname, '..', '..', '..', '..', 'package.json'), 'utf8')).version;
}

/** HEAD of the monorepo, HEAD of the public checkout and the submodule pin; `unknown` where git cannot answer. */
function gitRevisions(repoRoot) {
  const sha = (cwd, ref) => {
    const r = (0, _child_process.spawnSync)('git', ['rev-parse', '--verify', ref], {
      cwd,
      encoding: 'utf8'
    });
    return r.status === 0 ? r.stdout.trim() : 'unknown';
  };
  return {
    reddata: sha(repoRoot, 'HEAD'),
    public: sha(path.join(repoRoot, 'public'), 'HEAD'),
    public_pin: sha(repoRoot, 'HEAD:public')
  };
}

/** Content-addressed: two builds of the same revisions with the same builder and schema share a batch. */
function batchId(revisions, schemaVersion, builder) {
  const hash = (0, _crypto.createHash)('sha1').update([revisions.reddata, revisions.public, schemaVersion, builder].join('\n')).digest('hex');
  return `b-${hash.slice(0, 12)}`;
}

/** Public logical nodes with their descriptions; the home, the owner, non-public paths and every edge with a non-public end stay behind. */
function projectPublic(graph, system) {
  const keep = row => {
    const type = String(row.type);
    if (row.visibility !== 'public' || !PUBLIC_TYPES.some(t => (0, _types.isSubtype)(system, type, t))) return false;
    return !(0, _types.isSubtype)(system, type, 'doc-page') || row.kind === 'help';
  };
  const nodes = graph.nodes.filter(keep);
  const ids = new Set(nodes.map(n => n.id));
  const isPublicPath = p => typeof p === 'string' && (p.startsWith('public/') || p.startsWith('landing/'));
  const projected = nodes.map(row => {
    const out = {};
    const members = system.nodes.get(row.type).members;
    for (const [k, v] of Object.entries(row)) {
      if (k === 'home' || k === 'owner') continue;
      const m = members[k];
      if (m?.kind === 'ref' && !(m.list ? v.every(id => ids.has(id)) : ids.has(v))) continue;
      if (m?.scalar === 'Path' && !(m.list ? v.every(isPublicPath) : isPublicPath(v))) continue;
      out[k] = v;
    }
    return out;
  });
  const edges = graph.edges.filter(e => ids.has(e.from) && ids.has(e.to) && !(e.type === 'ref' && e.name === 'owner')).map(e => {
    const evidence = e.evidence?.filter(isPublicPath);
    const out = {
      ...e
    };
    delete out.evidence;
    return evidence?.length ? {
      ...out,
      evidence
    } : out;
  });
  return {
    ...graph,
    nodes: projected,
    edges
  };
}

/** Writes everything under [outRoot] and returns the manifest. `data/` and `reports/` are replaced whole. */
function writeBuild(graph, outRoot, info) {
  const dataDir = path.join(outRoot, 'data');
  const reportsDir = path.join(outRoot, 'reports');
  for (const dir of [dataDir, reportsDir]) fs.rmSync(dir, {
    recursive: true,
    force: true
  });
  fs.mkdirSync(path.join(dataDir, 'nodes'), {
    recursive: true
  });
  fs.mkdirSync(path.join(dataDir, 'edges'), {
    recursive: true
  });
  fs.mkdirSync(reportsDir, {
    recursive: true
  });
  const nodes = groupBy(graph.nodes, n => String(n.type));
  const edges = groupBy(graph.edges, e => e.type === 'ref' ? String(e.name) : String(e.type));
  const counts = {
    nodes: {},
    edges: {}
  };
  for (const [type, rows] of nodes) {
    rows.sort((a, b) => compare(String(a.id), String(b.id)));
    writeJsonl(path.join(dataDir, 'nodes', `${type}.jsonl`), rows);
    counts.nodes[type] = rows.length;
  }
  for (const [name, rows] of edges) {
    rows.sort((a, b) => compare(String(a.type), String(b.type)) || compare(String(a.from), String(b.from)) || compare(String(a.to), String(b.to)) || compare(String(a.name ?? ''), String(b.name ?? '')));
    writeJsonl(path.join(dataDir, 'edges', `${name}.jsonl`), rows);
    counts.edges[name] = rows.length;
  }
  const claims = [...graph.claims].sort((a, b) => compare(a.file, b.file) || compare(a.feature, b.feature) || a.rung - b.rung);
  writeJsonl(path.join(reportsDir, 'claims.jsonl'), claims);
  if (graph.invalid.length) writeJsonl(path.join(reportsDir, 'invalid.jsonl'), graph.invalid.slice(0, INVALID_CAP));
  if (Object.keys(graph.details).length) fs.writeFileSync(path.join(reportsDir, 'problems.json'), `${JSON.stringify(graph.details, null, 2)}\n`);
  const manifest = {
    built_at: new Date().toISOString(),
    schema_version: info.schemaVersion,
    builder: info.builder,
    batch: info.batch,
    mode: info.mode,
    revisions: info.revisions,
    sources: sortKeys(graph.sources),
    counts: {
      nodes: sortKeys(counts.nodes),
      edges: sortKeys(counts.edges)
    },
    problems: sortKeys(graph.problems)
  };
  fs.writeFileSync(path.join(outRoot, 'manifest.json'), `${JSON.stringify(manifest, null, 2)}\n`);
  return manifest;
}
function groupBy(rows, key) {
  const out = new Map();
  for (const row of rows) {
    let list = out.get(key(row));
    if (!list) out.set(key(row), list = []);
    list.push(row);
  }
  return new Map([...out].sort(([a], [b]) => compare(a, b)));
}

/** One object per line, keys in a fixed order: id, type, name, from, to, then the rest alphabetically. */
function writeJsonl(file, rows) {
  const lines = rows.map(row => {
    const keys = [...HEAD_KEYS.filter(k => row[k] !== undefined), ...Object.keys(row).filter(k => !HEAD_KEYS.includes(k)).sort(compare)];
    return JSON.stringify(Object.fromEntries(keys.map(k => [k, row[k]])));
  });
  fs.writeFileSync(file, lines.length ? `${lines.join('\n')}\n` : '');
}
function sortKeys(record) {
  return Object.fromEntries(Object.entries(record).sort(([a], [b]) => compare(a, b)));
}

/** Code-point order, the same on every platform and locale. */
function compare(a, b) {
  return a < b ? -1 : a > b ? 1 : 0;
}