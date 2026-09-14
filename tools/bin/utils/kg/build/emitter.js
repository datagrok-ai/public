"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.PROVENANCE_RANK = exports.Emitter = void 0;
var _types = require("../types");
var _normalize = require("./normalize");
var _ids = require("./ids");
/// Where the extractors put what they find (build-plan.md WO-1). Rows are normalized and shape-checked
/// on the way in, merged by id (nodes) or by type/from/to/name (edges), and `finalize()` derives the
/// hierarchy, inherits down it, materializes reference properties, checks endpoints, enforces required
/// members once and settles visibility.

const PROVENANCE_RANK = exports.PROVENANCE_RANK = ['annotation', 'ast', 'registry', 'filesystem', 'external', 'git', 'manual', 'llm'];
const VISIBILITY_ORDER = ['public', 'dev', 'internal'];
const EVIDENCE_CAP = 20;
const INVALID_CAP = 500;
const PROBLEM_KINDS = ['invalid_rows', 'dangling_edges', 'unresolved_ids', 'ambiguous_owners', 'orphans', 'partial_stubs'];
class Emitter {
  nodes = new Map();
  edges = new Map();
  claims = [];
  problems = Object.fromEntries(PROBLEM_KINDS.map(k => [k, 0]));
  details = {};
  invalid = [];
  sources = {};
  constructor(system, batch) {
    this.system = system;
    this.batch = batch;
  }
  node(row) {
    const {
      row: r,
      problems,
      defaulted
    } = (0, _normalize.normalizeRow)(this.system, row);
    if (typeof r.id !== 'string' || !r.id) problems.push({
      key: 'id',
      code: 'missing-key',
      message: 'no id'
    });
    if (problems.length) {
      this.reject(row, problems);
      return;
    }
    this.merge(r, new Set(defaulted));
  }

  /** A node that exists only because something references it; exempt from required members, never dropped. A second stub fills the gaps of the first. */
  stub(id, type, name, provenance, extra = {}) {
    const {
      row,
      problems,
      defaulted
    } = (0, _normalize.normalizeRow)(this.system, {
      ...extra,
      id,
      type,
      name,
      status: 'proposed',
      provenance,
      source_layer: 'synthetic'
    });
    if (problems.some(p => p.key === 'type')) {
      this.reject({
        ...extra,
        id,
        type
      }, problems);
      return;
    }
    const existing = this.nodes.get(id);
    if (existing) {
      if (existing.partial) for (const [k, v] of Object.entries(row)) if (existing.row[k] === undefined) {
        existing.row[k] = v;
        existing.prov[k] = provenance;
      }
      return;
    }
    const prov = Object.fromEntries(Object.keys(row).map(k => [k, provenance]));
    this.nodes.set(id, {
      row,
      type: this.system.nodes.get(type),
      prov,
      weak: new Set([...defaulted, 'status', 'source_layer', 'name']),
      partial: true
    });
  }
  edge(row) {
    const type = String(row.type ?? '');
    const edgeType = this.system.edges.get(type);
    if (!edgeType || edgeType.abstract && !(type === 'ref' && typeof row.name === 'string')) {
      this.reject(row, [{
        key: 'type',
        code: 'bad-type',
        message: edgeType ? `edge type '${type}' is abstract` : `unknown edge type '${type}'`
      }]);
      return;
    }
    const {
      row: r,
      problems
    } = (0, _normalize.normalizeEdgeRow)(this.system, edgeType, row);
    for (const end of ['from', 'to']) if (typeof r[end] !== 'string' || !r[end]) problems.push({
      key: end,
      code: 'missing-key',
      message: `no ${end}`
    });
    if (problems.length) {
      this.reject(row, problems);
      return;
    }
    const key = [type, r.from, r.to, r.name ?? ''].join('\u0000');
    const existing = this.edges.get(key);
    if (!existing) {
      this.edges.set(key, r);
      return;
    }
    const [winner, loser] = Number(r.confidence) > Number(existing.confidence) ? [r, existing] : [existing, r];
    const merged = {
      ...loser,
      ...winner
    };
    const evidence = [...new Set([...(existing.evidence ?? []), ...(r.evidence ?? [])])];
    if (evidence.length) merged.evidence = evidence.slice(0, EVIDENCE_CAP);
    this.edges.set(key, merged);
  }

  /** A reference property as an edge line named by the property (conventions.md §7.6). */
  ref(from, name, to, provenance, confidence = 1) {
    this.edge({
      type: 'ref',
      name,
      from,
      to,
      derived_by: provenance,
      confidence
    });
  }
  claim(c) {
    this.claims.push(c);
  }
  source(name, status) {
    this.sources[name] = status;
  }
  problem(kind, detail) {
    this.problems[kind] = (this.problems[kind] ?? 0) + 1;
    if (detail) (this.details[kind] ??= []).push(detail);
  }
  finalize() {
    this.resolveEndpoints();
    const parents = this.derivePartOf();
    this.inherit(parents);
    this.refLines();
    this.resolveEndpoints();
    for (const entry of this.nodes.values()) entry.row.batch ??= this.batch;
    for (const row of this.edges.values()) row.batch ??= this.batch;
    this.requireMembers();
    for (const entry of this.nodes.values()) this.settleVisibility(entry);
    this.problems.partial_stubs = [...this.nodes.values()].filter(e => e.partial).length;
    return {
      nodes: [...this.nodes.values()].map(e => e.row),
      edges: [...this.edges.values()],
      claims: this.claims,
      sources: this.sources,
      problems: this.problems,
      details: this.details,
      invalid: this.invalid
    };
  }
  reject(row, problems) {
    this.problems.invalid_rows++;
    if (this.invalid.length < INVALID_CAP) this.invalid.push({
      ...row,
      problems: problems.map(p => p.message)
    });
  }
  merge(incoming, defaulted) {
    const id = incoming.id;
    const incomingType = this.system.nodes.get(incoming.type);
    const provenance = String(incoming.provenance ?? '');
    const existing = this.nodes.get(id);
    if (!existing) {
      this.nodes.set(id, {
        row: incoming,
        type: incomingType,
        prov: Object.fromEntries(Object.keys(incoming).map(k => [k, provenance])),
        weak: defaulted,
        partial: false
      });
      return;
    }
    let type = existing.type;
    if (incomingType !== existing.type) {
      if ((0, _types.isSubtype)(this.system, incomingType.name, existing.type.name)) type = incomingType;else if (!(0, _types.isSubtype)(this.system, existing.type.name, incomingType.name)) {
        this.reject(incoming, [{
          key: 'type',
          code: 'bad-type',
          message: `type '${incomingType.name}' conflicts with '${existing.type.name}' already asserted for ${id}`
        }]);
        return;
      }
    }
    if (existing.partial) {
      const prov = Object.fromEntries(Object.keys(incoming).map(k => [k, provenance]));
      this.nodes.set(id, {
        row: {
          ...existing.row,
          ...incoming,
          type: type.name
        },
        type,
        prov: {
          ...existing.prov,
          ...prov
        },
        weak: defaulted,
        partial: false
      });
      return;
    }
    const {
      row,
      prov,
      weak
    } = existing;
    for (const [k, v] of Object.entries(incoming)) {
      if (k === 'type') continue;
      const cur = row[k];
      const incomingWeak = defaulted.has(k);
      if (cur === undefined) {
        row[k] = v;
        prov[k] = provenance;
        if (incomingWeak) weak.add(k);
        continue;
      }
      if (Array.isArray(cur) && Array.isArray(v)) {
        row[k] = [...new Set([...cur, ...v].map(x => JSON.stringify(x)))].map(x => JSON.parse(x));
        if (rank(provenance) < rank(prov[k])) prov[k] = provenance;
        continue;
      }
      if (JSON.stringify(cur) === JSON.stringify(v)) {
        if (!incomingWeak) weak.delete(k);
        continue;
      }
      const incomingWins = weak.has(k) !== incomingWeak ? weak.has(k) : rank(provenance) < rank(prov[k]);
      if (!incomingWins) continue;
      row[k] = v;
      prov[k] = provenance;
      if (incomingWeak) weak.add(k);else weak.delete(k);
    }
    row.type = type.name;
    existing.type = type;
  }

  /** Part-of from the id path of every hierarchical node; missing parents become stubs. Returns child -> parent. */
  derivePartOf() {
    const parents = new Map();
    if (!this.system.edges.has('part-of')) return parents;
    for (const [id, entry] of [...this.nodes]) {
      if (!entry.type.hierarchical) continue;
      const prefixed = _ids.PREFIXED_ID.exec(id);
      const prefix = prefixed ? `${prefixed[1]}:` : '';
      const segments = (prefixed ? prefixed[2] : id).split('/');
      for (let i = segments.length - 1; i >= 1; i--) {
        const child = prefix + segments.slice(0, i + 1).join('/');
        const parent = prefix + segments.slice(0, i).join('/');
        if (!this.nodes.has(parent)) this.stub(parent, entry.type.name, (0, _ids.titleCase)(segments[i - 1]), 'filesystem');
        parents.set(child, parent);
        this.edge({
          type: 'part-of',
          from: child,
          to: parent,
          derived_by: 'filesystem',
          confidence: 1
        });
      }
    }
    return parents;
  }

  /** status, visibility and owner flow down part-of unless the child sets them (nodes/node.yaml `inherit`). */
  inherit(parents) {
    const depth = id => id.split('/').length;
    for (const child of [...parents.keys()].sort((a, b) => depth(a) - depth(b))) {
      const entry = this.nodes.get(child);
      const parent = this.nodes.get(parents.get(child));
      const inheritable = new Set(entry.type.chain.flatMap(t => this.system.nodes.get(t)?.inherit ?? []));
      for (const m of inheritable) {
        if (parent.row[m] === undefined || parent.weak.has(m)) continue;
        if (entry.row[m] !== undefined && !entry.weak.has(m)) continue;
        entry.row[m] = parent.row[m];
        entry.prov[m] = parent.prov[m];
        entry.weak.delete(m);
      }
    }
  }
  refLines() {
    for (const [id, entry] of this.nodes) for (const member of Object.values(entry.type.members)) {
      if (member.kind !== 'ref' || entry.row[member.name] === undefined) continue;
      const targets = member.list ? entry.row[member.name] : [entry.row[member.name]];
      for (const to of targets) this.ref(id, member.name, to, String(entry.row.provenance ?? entry.prov[member.name] ?? 'filesystem'));
    }
  }

  /** Endpoint types per edge type with subtypes; a missing target of an authored edge becomes a stub, anything else is dropped. */
  resolveEndpoints() {
    for (const [key, row] of [...this.edges]) {
      const edgeType = this.system.edges.get(row.type);
      const from = this.nodes.get(row.from);
      const expectedTo = row.type === 'ref' ? from?.type.members[row.name]?.refs ?? [] : edgeType.to;
      const ok = this.endpoint(row, 'from', edgeType.from) && this.endpoint(row, 'to', expectedTo);
      if (!ok) {
        this.edges.delete(key);
        this.problems.dangling_edges++;
      }
    }
  }
  endpoint(row, side, expected) {
    const id = row[side];
    const entry = this.nodes.get(id);
    if (entry) return expected.some(t => (0, _types.isSubtype)(this.system, entry.type.name, t));
    if (row.derived_by !== 'annotation') return false;
    const type = this.typeOfId(id, expected);
    if (!type || !expected.some(t => (0, _types.isSubtype)(this.system, type, t))) return false;
    this.stub(id, type, (0, _ids.stubName)(id), 'annotation');
    return this.nodes.has(id);
  }

  /** The concrete type an id's shape implies, within [expected]. */
  typeOfId(id, expected) {
    const prefixed = _ids.PREFIXED_ID.exec(id);
    if (prefixed) return this.system.prefixes.get(prefixed[1]);
    if (_ids.JIRA_KEY.test(id)) return 'ticket';
    const schemed = _ids.SCHEMED_ID.exec(id);
    if (schemed) return _ids.SCHEME_TYPES[schemed[1]]?.find(t => this.system.nodes.has(t));
    const type = (0, _types.concreteAuthored)(this.system, expected).find(t => !t.prefix);
    return type?.name;
  }

  /** A row that is neither complete nor a stub after merging is invalid; edges then left without an end follow it. */
  requireMembers() {
    for (const [id, entry] of [...this.nodes]) {
      if (entry.partial) continue;
      const required = [...Object.values(entry.type.members), ...this.system.buildFields.node].filter(m => !m.nullable && entry.row[m.name] === undefined);
      if (!required.length) continue;
      this.reject(entry.row, required.map(m => ({
        key: m.name,
        code: 'missing-key',
        message: `missing required member '${m.name}' for type ${entry.type.name}`
      })));
      this.nodes.delete(id);
    }
    for (const [key, row] of [...this.edges]) {
      const edgeType = this.system.edges.get(row.type);
      const required = [...Object.values(edgeType.properties), ...this.system.buildFields.edge].filter(m => !m.nullable && row[m.name] === undefined);
      if (required.length) {
        this.reject(row, required.map(m => ({
          key: m.name,
          code: 'missing-key',
          message: `missing required property '${m.name}' for edge ${edgeType.name}`
        })));
        this.edges.delete(key);
      } else if (!this.nodes.has(row.from) || !this.nodes.has(row.to)) {
        this.edges.delete(key);
        this.problems.dangling_edges++;
      }
    }
  }

  /** A logical node: its type default narrowed by the home; a node with a path: narrowed by its location too. */
  settleVisibility(entry) {
    const {
      row,
      type
    } = entry;
    const candidates = [type.visibility, row.visibility];
    if (typeof row.path === 'string' && type.members.path?.scalar === 'Path') candidates.push((0, _ids.locationVisibility)(row.path));
    const ranks = candidates.map(v => VISIBILITY_ORDER.indexOf(v ?? '')).filter(r => r >= 0);
    row.visibility = ranks.length ? VISIBILITY_ORDER[Math.max(...ranks)] : 'dev';
  }
}
exports.Emitter = Emitter;
function rank(provenance) {
  const i = PROVENANCE_RANK.indexOf(provenance ?? '');
  return i < 0 ? PROVENANCE_RANK.length : i;
}