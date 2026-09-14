/// Where the extractors put what they find (build-plan.md WO-1). Rows are normalized and shape-checked
/// on the way in, merged by id (nodes) or by type/from/to/name (edges), and `finalize()` derives the
/// hierarchy, inherits down it, materializes reference properties, checks endpoints, enforces required
/// members once and settles visibility.
import {TypeSystem, NodeType, isSubtype, concreteAuthored} from '../types';
import {normalizeRow, normalizeEdgeRow, Row, RowProblem} from './normalize';
import {PREFIXED_ID, SCHEMED_ID, SCHEME_TYPES, JIRA_KEY, parseId, stubName, titleCase, locationVisibility, sourceLayerOf} from './ids';

export interface Claim {
  /** Posix path relative to the monorepo root. */
  file: string;
  feature: string;
  rung: 1 | 2 | 3;
  source: 'home' | 'marker';
  props: Record<string, unknown>;
  line?: number;
  /** An inline `// ~id` marker (conventions.md §6): participates-in only, never ownership. */
  mode?: 'participates';
}

export interface Graph {
  nodes: Row[];
  edges: Row[];
  /** Ids of the partial stubs among the nodes. */
  stubs: string[];
  claims: Claim[];
  sources: Record<string, string>;
  problems: Record<string, number>;
  /** Free-text detail per problem kind, for the reports. */
  details: Record<string, string[]>;
  /** The first rows rejected, with the reasons. */
  invalid: Row[];
  /** Structured reports an extractor built, by file name under `reports/`. */
  reports: Record<string, unknown>;
}

export const PROVENANCE_RANK = ['annotation', 'ast', 'registry', 'filesystem', 'external', 'git', 'manual', 'llm'];
const VISIBILITY_ORDER = ['public', 'dev', 'internal'];
const EVIDENCE_CAP = 20;
const INVALID_CAP = 500;
const PROBLEM_KINDS = ['invalid_rows', 'dangling_edges', 'unresolved_ids', 'ambiguous_owners', 'orphans', 'partial_stubs', 'duplicate_ids'];

interface NodeEntry {
  row: Row;
  type: NodeType;
  /** Provenance of the row each field came from. */
  prov: Record<string, string>;
  /** Fields holding a declared default or a stub placeholder, so inheritance and a real row may replace them. */
  weak: Set<string>;
  partial: boolean;
}

export class Emitter {
  private nodes = new Map<string, NodeEntry>();
  private edges = new Map<string, Row>();
  private claims: Claim[] = [];
  private problems: Record<string, number> = Object.fromEntries(PROBLEM_KINDS.map((k) => [k, 0]));
  private details: Record<string, string[]> = {};
  private invalid: Row[] = [];
  private reports: Record<string, unknown> = {};
  readonly sources: Record<string, string> = {};

  constructor(private system: TypeSystem, readonly batch: string) {}

  node(row: Row): void {
    const {row: r, problems, defaulted} = normalizeRow(this.system, row);
    if (typeof r.id !== 'string' || !r.id) problems.push({key: 'id', code: 'missing-key', message: 'no id'});
    if (problems.length) {
      this.reject(row, problems);
      return;
    }
    this.merge(r, new Set(defaulted));
  }

  /** A node that exists only because something references it: carries only what created it (no defaults), is exempt from
   * required members and never dropped. A second stub fills the gaps of the first. */
  stub(id: string, type: string, name: string, provenance: string, extra: Row = {}): void {
    const layer = typeof extra.path === 'string' ? sourceLayerOf(extra.path) : 'synthetic';
    const {row, problems} = normalizeRow(this.system, {...extra, id, type, name, status: 'proposed', provenance, source_layer: layer}, {defaults: false});
    if (problems.some((p) => p.key === 'type')) {
      this.reject({...extra, id, type}, problems);
      return;
    }
    const existing = this.nodes.get(id);
    if (existing) {
      if (existing.partial)
        for (const [k, v] of Object.entries(row))
          if (existing.row[k] === undefined) {
            existing.row[k] = v;
            existing.prov[k] = provenance;
          }
      return;
    }
    const prov = Object.fromEntries(Object.keys(row).map((k) => [k, provenance]));
    this.nodes.set(id, {row, type: this.system.nodes.get(type)!, prov, weak: new Set(['status', 'source_layer', 'name']), partial: true});
  }

  edge(row: Row): void {
    const type = String(row.type ?? '');
    const edgeType = this.system.edges.get(type);
    if (!edgeType || (edgeType.abstract && !(type === 'ref' && typeof row.name === 'string'))) {
      this.reject(row, [{key: 'type', code: 'bad-type', message: edgeType ? `edge type '${type}' is abstract` : `unknown edge type '${type}'`}]);
      return;
    }
    const {row: r, problems} = normalizeEdgeRow(this.system, edgeType, row);
    for (const end of ['from', 'to'])
      if (typeof r[end] !== 'string' || !r[end]) problems.push({key: end, code: 'missing-key', message: `no ${end}`});
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
    const merged: Row = {...loser, ...winner};
    const evidence = [...new Set([...(existing.evidence as string[] ?? []), ...(r.evidence as string[] ?? [])])];
    if (evidence.length) merged.evidence = evidence.slice(0, EVIDENCE_CAP);
    this.edges.set(key, merged);
  }

  /** A reference property as an edge line named by the property (conventions.md §7.6). */
  ref(from: string, name: string, to: string, provenance: string, confidence = 1): void {
    this.edge({type: 'ref', name, from, to, derived_by: provenance, confidence});
  }

  claim(c: Claim): void {
    this.claims.push(c);
  }

  /** What membership resolution (WO-4) reads: every claim made so far. */
  get claimed(): Claim[] {
    return this.claims;
  }

  /** The rows emitted so far whose type is [type] or narrows it; membership resolution reads files and tests this way. */
  rowsOf(type: string): Row[] {
    return [...this.nodes.values()].filter((e) => isSubtype(this.system, e.type.name, type)).map((e) => e.row);
  }

  /** A structured report an extractor built; the writer puts it in `reports/<name>.json`. */
  report(name: string, data: unknown): void {
    this.reports[name] = data;
  }

  source(name: string, status: string): void {
    this.sources[name] = status;
  }

  problem(kind: string, detail?: string): void {
    this.problems[kind] = (this.problems[kind] ?? 0) + 1;
    if (detail) (this.details[kind] ??= []).push(detail);
  }

  finalize(): Graph {
    this.resolveEndpoints();
    const parents = this.derivePartOf();
    this.inherit(parents);
    this.refLines();
    this.resolveEndpoints();
    for (const entry of this.nodes.values()) entry.row.batch ??= this.batch;
    for (const row of this.edges.values()) row.batch ??= this.batch;
    this.requireMembers();
    for (const entry of this.nodes.values()) this.settleVisibility(entry);
    const stubs = [...this.nodes.values()].filter((e) => e.partial).map((e) => e.row.id as string);
    this.problems.partial_stubs = stubs.length;
    return {
      nodes: [...this.nodes.values()].map((e) => e.row),
      edges: [...this.edges.values()],
      stubs,
      claims: this.claims,
      sources: this.sources,
      problems: this.problems,
      details: this.details,
      invalid: this.invalid,
      reports: this.reports,
    };
  }

  private reject(row: Row, problems: RowProblem[]): void {
    this.problems.invalid_rows++;
    if (this.invalid.length < INVALID_CAP) this.invalid.push({...row, problems: problems.map((p) => p.message)});
  }

  private merge(incoming: Row, defaulted: Set<string>): void {
    const id = incoming.id as string;
    const incomingType = this.system.nodes.get(incoming.type as string)!;
    const provenance = String(incoming.provenance ?? '');
    const existing = this.nodes.get(id);
    if (!existing) {
      this.nodes.set(id, {row: incoming, type: incomingType, prov: Object.fromEntries(Object.keys(incoming).map((k) => [k, provenance])), weak: defaulted, partial: false});
      return;
    }
    let type = existing.type;
    if (incomingType !== existing.type) {
      if (isSubtype(this.system, incomingType.name, existing.type.name)) type = incomingType;
      else if (!isSubtype(this.system, existing.type.name, incomingType.name)) {
        this.duplicateId(existing, incoming, incomingType);
        return;
      }
    }
    if (!existing.partial && elsewhere(id, existing.row, incoming)) {
      this.duplicateId(existing, incoming, incomingType);
      return;
    }
    if (existing.partial) {
      const prov = Object.fromEntries(Object.keys(incoming).map((k) => [k, provenance]));
      const entry: NodeEntry = {row: {...existing.row, ...incoming, type: type.name}, type, prov: {...existing.prov, ...prov}, weak: defaulted, partial: false};
      this.nodes.set(id, entry);
      this.revalidate(id, entry);
      return;
    }
    const {row, prov, weak} = existing;
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
        row[k] = [...new Set([...cur, ...v].map((x) => JSON.stringify(x)))].map((x) => JSON.parse(x));
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
      if (incomingWeak) weak.add(k);
      else weak.delete(k);
    }
    row.type = type.name;
    existing.type = type;
    this.revalidate(id, existing);
  }

  /** Two rows for one id that are different declarations: the first one stands, the second is reported, never merged. */
  private duplicateId(existing: NodeEntry, incoming: Row, incomingType: NodeType): void {
    this.problem('duplicate_ids', `${incoming.id}: ${incomingType.name} at ${where(incoming)} ignored; ${existing.type.name} at ${where(existing.row)} kept`);
  }

  /** A merged row must still satisfy the winning type, which a member narrowed by a subtype (a query's language) can break. */
  private revalidate(id: string, entry: NodeEntry): void {
    const {problems} = normalizeRow(this.system, entry.row, {defaults: false});
    if (!problems.length) return;
    this.reject(entry.row, problems);
    this.nodes.delete(id);
  }

  /** Part-of from the id path of every hierarchical node; missing parents become stubs. Returns child -> parent. */
  private derivePartOf(): Map<string, string> {
    const parents = new Map<string, string>();
    if (!this.system.edges.has('part-of')) return parents;
    for (const [id, entry] of [...this.nodes]) {
      if (!entry.type.hierarchical) continue;
      const prefixed = PREFIXED_ID.exec(id);
      const prefix = prefixed ? `${prefixed[1]}:` : '';
      const segments = (prefixed ? prefixed[2] : id).split('/');
      for (let i = segments.length - 1; i >= 1; i--) {
        const child = prefix + segments.slice(0, i + 1).join('/');
        const parent = prefix + segments.slice(0, i).join('/');
        if (!this.nodes.has(parent)) this.stub(parent, entry.type.name, titleCase(segments[i - 1]), 'filesystem');
        parents.set(child, parent);
        this.edge({type: 'part-of', from: child, to: parent, derived_by: 'filesystem', confidence: 1});
      }
    }
    return parents;
  }

  /** status, visibility and owner flow down part-of unless the child sets them (nodes/node.yaml `inherit`). */
  private inherit(parents: Map<string, string>): void {
    const depth = (id: string) => id.split('/').length;
    for (const child of [...parents.keys()].sort((a, b) => depth(a) - depth(b))) {
      const entry = this.nodes.get(child)!;
      const parent = this.nodes.get(parents.get(child)!)!;
      const inheritable = new Set(entry.type.chain.flatMap((t) => this.system.nodes.get(t)?.inherit ?? []));
      for (const m of inheritable) {
        if (parent.row[m] === undefined || parent.weak.has(m)) continue;
        if (entry.row[m] !== undefined && !entry.weak.has(m)) continue;
        entry.row[m] = parent.row[m];
        entry.prov[m] = parent.prov[m];
        entry.weak.delete(m);
      }
    }
  }

  private refLines(): void {
    for (const [id, entry] of this.nodes)
      for (const member of Object.values(entry.type.members)) {
        if (member.kind !== 'ref' || entry.row[member.name] === undefined) continue;
        const targets = member.list ? entry.row[member.name] as string[] : [entry.row[member.name] as string];
        for (const to of targets) this.ref(id, member.name, to, String(entry.row.provenance ?? entry.prov[member.name] ?? 'filesystem'));
      }
  }

  /** Endpoint types per edge type with subtypes; a missing target of an authored edge becomes a stub, anything else is dropped. */
  private resolveEndpoints(): void {
    for (const [key, row] of [...this.edges]) {
      const edgeType = this.system.edges.get(row.type as string)!;
      const from = this.nodes.get(row.from as string);
      const expectedTo = row.type === 'ref' ? from?.type.members[row.name as string]?.refs ?? [] : edgeType.to;
      const ok = this.endpoint(row, 'from', edgeType.from) && this.endpoint(row, 'to', expectedTo);
      if (!ok) {
        this.edges.delete(key);
        this.problems.dangling_edges++;
      }
    }
  }

  private endpoint(row: Row, side: 'from' | 'to', expected: string[]): boolean {
    const id = row[side] as string;
    const entry = this.nodes.get(id);
    if (entry) return expected.some((t) => isSubtype(this.system, entry.type.name, t));
    if (row.derived_by !== 'annotation') return false;
    const type = this.typeOfId(id, expected);
    if (!type || !expected.some((t) => isSubtype(this.system, type, t))) return false;
    this.stub(id, type, stubName(id), 'annotation');
    return this.nodes.has(id);
  }

  /** The concrete type an id's shape implies, within [expected]. */
  private typeOfId(id: string, expected: string[]): string | undefined {
    const prefixed = PREFIXED_ID.exec(id);
    if (prefixed) return this.system.prefixes.get(prefixed[1]);
    if (JIRA_KEY.test(id)) return 'ticket';
    const schemed = SCHEMED_ID.exec(id);
    if (schemed) return SCHEME_TYPES[schemed[1]]?.find((t) => this.system.nodes.has(t));
    const type = concreteAuthored(this.system, expected).find((t) => !t.prefix);
    return type?.name;
  }

  /** A row that is neither complete nor a stub after merging is invalid; edges then left without an end follow it. */
  private requireMembers(): void {
    for (const [id, entry] of [...this.nodes]) {
      if (entry.partial) continue;
      const required = [...Object.values(entry.type.members), ...this.system.buildFields.node].filter((m) => !m.nullable && entry.row[m.name] === undefined);
      if (!required.length) continue;
      this.reject(entry.row, required.map((m) => ({key: m.name, code: 'missing-key', message: `missing required member '${m.name}' for type ${entry.type.name}`})));
      this.nodes.delete(id);
    }
    for (const [key, row] of [...this.edges]) {
      const edgeType = this.system.edges.get(row.type as string)!;
      const required = [...Object.values(edgeType.properties), ...this.system.buildFields.edge].filter((m) => !m.nullable && row[m.name] === undefined);
      if (required.length) {
        this.reject(row, required.map((m) => ({key: m.name, code: 'missing-key', message: `missing required property '${m.name}' for edge ${edgeType.name}`})));
        this.edges.delete(key);
      }
      else if (!this.nodes.has(row.from as string) || !this.nodes.has(row.to as string)) {
        this.edges.delete(key);
        this.problems.dangling_edges++;
      }
    }
  }

  /** A logical node: its type default narrowed by the home; a node with a path: narrowed by its location too. */
  private settleVisibility(entry: NodeEntry): void {
    const {row, type} = entry;
    const candidates = [type.visibility, row.visibility as string | undefined];
    if (typeof row.path === 'string' && type.members.path?.scalar === 'Path') candidates.push(locationVisibility(row.path));
    const ranks = candidates.map((v) => VISIBILITY_ORDER.indexOf(v ?? '')).filter((r) => r >= 0);
    row.visibility = ranks.length ? VISIBILITY_ORDER[Math.max(...ranks)] : 'dev';
  }
}

/** Whether [incoming] names another place than the row already held; an id that embeds its path (decl:, file:, doc:) never does. */
function elsewhere(id: string, existing: Row, incoming: Row): boolean {
  if (parseId(id).path !== undefined) return false;
  const differs = (key: string) => existing[key] !== undefined && incoming[key] !== undefined && existing[key] !== incoming[key];
  return differs('path') || differs('line');
}

function where(row: Row): string {
  return row.path === undefined ? String(row.provenance ?? 'unknown') : `${row.path}${row.line === undefined ? '' : `:${row.line}`}`;
}

function rank(provenance: string | undefined): number {
  const i = PROVENANCE_RANK.indexOf(provenance ?? '');
  return i < 0 ? PROVENANCE_RANK.length : i;
}
