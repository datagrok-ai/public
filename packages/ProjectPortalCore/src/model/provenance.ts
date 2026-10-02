/** Provenance (F-28, FR-75): the portal's own records as a graph, and the walk from any of them to
 * what it rests on and what rests on it. Data → model run → document → decision, and back.
 *
 * Nodes are the portal's records: a program, its decisions (gate-log entries), model runs and the
 * tables they read and wrote, the catalogued model a run used, documents, snapshots, gate
 * criteria, the program dataset, a compound, and the files and links a decision cites. An edge
 * runs from what something rests on to what rests on it:
 * - a compound is in the program dataset; the dataset is frozen in a snapshot;
 * - a run reads its input tables and the model, and writes its output tables;
 * - a run selected for a criterion feeds the documents answering that criterion (or the criterion
 *   when none does); a document answers its criterion;
 * - every piece of a decision's evidence points at the decision; a decision moves the program.
 * Evidence is the edge, not a node: starting from a piece of evidence starts from its target.
 * Pure: the caller reads the records. */

export type ProvenanceKind = 'program' | 'decision' | 'model_run' | 'snapshot' | 'document' | 'compound' | 'criterion' |
  'dataset' | 'table' | 'model' | 'file' | 'url';
/** What a caller may start from: a record of the portal, or a piece of evidence. */
export type ProvenanceStart = 'program' | 'decision' | 'evidence' | 'model_run' | 'snapshot' | 'document' | 'compound';
export const PROVENANCE_STARTS: ProvenanceStart[] = ['program', 'decision', 'evidence', 'model_run', 'snapshot', 'document',
  'compound'];

export const PROVENANCE_LABELS: Record<ProvenanceKind, string> = {
  program: 'Program', decision: 'Decision', model_run: 'Run', snapshot: 'Snapshot', document: 'Document', compound: 'Compound',
  criterion: 'Criterion', dataset: 'Data', table: 'Table', model: 'Model', file: 'File', url: 'Link',
};

export interface ProvenanceNode {
  key: string;
  kind: ProvenanceKind;
  /** The record's id: a model_run row id, a stage_transition id; a path or URL for a file or a link. */
  id: string;
  name: string;
  detail?: string;
}

/** [from] is what [to] rests on. */
export interface ProvenanceEdge {
  from: string;
  to: string;
  relation: string;
}

export const nodeKey = (kind: ProvenanceKind, id: string) => `${kind}:${id}`;

/** The records of one or more programs, as the caller read them. */
export interface ProvenanceRecords {
  programs: {id: string, code: string, name: string}[];
  decisions: {id: string, programId: string, decision: string, on: string, from?: string | null, to?: string | null}[];
  evidence: {id: string, decisionId: string | null, kind: string, target: string, name: string,
    context?: Record<string, any> | null}[];
  /** Model runs: the row id, the saved call's id, the tables it read and wrote. */
  runs: {id: string, programId: string, runId?: string | null, title: string, model?: {code: string, name: string} | null,
    inputs?: {id: string, name: string}[], outputs?: {id: string, name: string}[]}[];
  /** The run selected for a criterion: model_run row id and criterion id. */
  selections: {programId: string, criterionId: string, runId: string}[];
  documents: {id: string, programId: string, title: string, version?: string | null, status: string,
    criterionId?: string | null, path: string}[];
  snapshots: {id: string, programId: string, label: string, takenOn: string}[];
  criteria: {id: string, code: string, name: string, stage?: string}[];
  /** The program dataset of each program, when the portfolio declares one. */
  datasets: {programId: string, name: string}[];
  /** Compounds held by the programs, with the program and role; only those asked about. */
  compounds: {id: string, code: string, programId: string, role?: string | null}[];
}

export interface ProvenanceGraph {
  nodes: Map<string, ProvenanceNode>;
  edges: ProvenanceEdge[];
}

const datasetId = (programId: string) => `${programId}:data`;

/** The key of an evidence row's target, as the graph names it. */
export function evidenceTargetKey(e: ProvenanceRecords['evidence'][number], runs: ProvenanceRecords['runs']): string {
  switch (e.kind) {
  case 'model_run': {
    const run = runs.find((r) => r.runId === e.target || r.id === e.target);
    return nodeKey('model_run', run?.id ?? e.target);
  }
  case 'snapshot':
    return nodeKey('snapshot', e.target);
  case 'criterion':
    return nodeKey('criterion', e.target);
  case 'url':
    return nodeKey('url', e.target);
  default:
    return e.context?.document ? nodeKey('document', String(e.context.document)) : nodeKey('file', e.target);
  }
}

/** The graph over [r]. An evidence target that is not among the records (a run deleted since) is
 * still a node, named as the evidence froze it. */
export function provenanceGraph(r: ProvenanceRecords): ProvenanceGraph {
  const nodes = new Map<string, ProvenanceNode>();
  const edges: ProvenanceEdge[] = [];
  const seen = new Set<string>();
  const add = (n: Omit<ProvenanceNode, 'key'>): string => {
    const key = nodeKey(n.kind, n.id);
    if (!nodes.has(key))
      nodes.set(key, {key, ...n});
    return key;
  };
  const link = (from: string, to: string, relation: string) => {
    const k = `${from}>${to}`;
    if (from !== to && !seen.has(k) && seen.add(k))
      edges.push({from, to, relation});
  };
  const programCode = new Map(r.programs.map((p) => [p.id, p.code]));
  for (const p of r.programs)
    add({kind: 'program', id: p.id, name: `${p.code} · ${p.name}`});
  for (const d of r.datasets)
    add({kind: 'dataset', id: datasetId(d.programId), name: `${d.name} of ${programCode.get(d.programId) ?? 'the program'}`});
  for (const c of r.criteria)
    add({kind: 'criterion', id: c.id, name: c.name, detail: [c.stage, c.code].filter((x) => x).join(' · ')});
  for (const c of r.compounds) {
    const k = add({kind: 'compound', id: c.id, name: c.code, detail: c.role ?? undefined});
    if (r.datasets.some((d) => d.programId === c.programId))
      link(k, nodeKey('dataset', datasetId(c.programId)), 'is in');
  }
  for (const s of r.snapshots) {
    const k = add({kind: 'snapshot', id: s.id, name: s.label, detail: s.takenOn.slice(0, 16).replace('T', ' ')});
    if (r.datasets.some((d) => d.programId === s.programId))
      link(nodeKey('dataset', datasetId(s.programId)), k, 'frozen in');
  }
  for (const run of r.runs) {
    const k = add({kind: 'model_run', id: run.id, name: run.title, detail: run.model?.name});
    if (run.model)
      link(add({kind: 'model', id: run.model.code, name: run.model.name}), k, 'model of');
    for (const t of run.inputs ?? [])
      link(add({kind: 'table', id: t.id, name: t.name}), k, 'input of');
    for (const t of run.outputs ?? [])
      link(k, add({kind: 'table', id: t.id, name: t.name}), 'output');
  }
  for (const d of r.documents) {
    const k = add({kind: 'document', id: d.id, name: d.title,
      detail: [d.version ? `v${d.version}` : '', d.status].filter((x) => x).join(' ')});
    if (d.criterionId && nodes.has(nodeKey('criterion', d.criterionId)))
      link(k, nodeKey('criterion', d.criterionId), 'answers');
  }
  for (const s of r.selections) {
    const run = nodeKey('model_run', s.runId);
    if (!nodes.has(run))
      continue;
    const docs = r.documents.filter((d) => d.programId === s.programId && d.criterionId === s.criterionId);
    if (docs.length)
      docs.forEach((d) => link(run, nodeKey('document', d.id), 'selected for its criterion'));
    else if (nodes.has(nodeKey('criterion', s.criterionId)))
      link(run, nodeKey('criterion', s.criterionId), 'selected for');
  }
  for (const d of r.decisions) {
    const k = add({kind: 'decision', id: d.id, name: `${d.decision} ${d.on.slice(0, 10)}`,
      detail: d.from || d.to ? `${d.from ?? '—'} → ${d.to ?? '—'}` : undefined});
    if (nodes.has(nodeKey('program', d.programId)))
      link(k, nodeKey('program', d.programId), 'moved');
  }
  for (const e of r.evidence) {
    if (!e.decisionId || !nodes.has(nodeKey('decision', e.decisionId)))
      continue;
    const key = evidenceTargetKey(e, r.runs);
    if (!nodes.has(key)) {
      const [kind, ...id] = key.split(':');
      // a record deleted since reads as the evidence froze it; a file or link is only ever that
      add({kind: kind as ProvenanceKind, id: id.join(':'), name: e.name,
        detail: kind === 'file' || kind === 'url' ? undefined : 'as recorded'});
    }
    link(key, nodeKey('decision', e.decisionId), 'evidence for');
  }
  return {nodes, edges};
}

/** Where a walk starts: the record itself, or for a piece of evidence its target. Null when [id]
 * is not among the records. */
export function startKey(g: ProvenanceGraph, r: ProvenanceRecords, kind: ProvenanceStart, id: string): string | null {
  if (kind === 'evidence') {
    const e = r.evidence.find((x) => x.id === id);
    return e ? evidenceTargetKey(e, r.runs) : null;
  }
  if (kind === 'model_run') {
    const run = r.runs.find((x) => x.id === id || x.runId === id);
    return run ? nodeKey('model_run', run.id) : null;
  }
  if (kind === 'compound') {
    const c = r.compounds.find((x) => x.id === id || x.code === id);
    return c ? nodeKey('compound', c.id) : null;
  }
  const key = nodeKey(kind, id);
  return g.nodes.has(key) ? key : null;
}

/** A node reached from the start: how far, which way, and the edge it was reached by. */
export interface ReachedNode extends ProvenanceNode {
  depth: number;
  direction: 'start' | 'up' | 'down';
  /** The node it was reached from, and the edge's relation; empty for the start. */
  via?: string;
  relation?: string;
}

export interface Provenance {
  start: string;
  nodes: ReachedNode[];
  /** Every edge between two reached nodes. */
  edges: ProvenanceEdge[];
}

/** Walks [depth] steps up (what [start] rests on) and down (what rests on it), breadth first. A node
 * reached both ways is kept the way it was reached first, up before down. */
export function traverse(g: ProvenanceGraph, start: string, depth = 4): Provenance {
  const reached = new Map<string, ReachedNode>();
  const s = g.nodes.get(start);
  if (!s)
    return {start, nodes: [], edges: []};
  reached.set(start, {...s, depth: 0, direction: 'start'});
  for (const direction of ['up', 'down'] as const) {
    let frontier = [start];
    for (let d = 1; d <= depth && frontier.length; d++) {
      const next: string[] = [];
      for (const k of frontier) {
        const steps = direction === 'up' ? g.edges.filter((e) => e.to === k).map((e) => ({key: e.from, relation: e.relation})) :
          g.edges.filter((e) => e.from === k).map((e) => ({key: e.to, relation: e.relation}));
        for (const step of steps) {
          if (reached.has(step.key))
            continue;
          reached.set(step.key, {...g.nodes.get(step.key)!, depth: d, direction, via: k, relation: step.relation});
          next.push(step.key);
        }
      }
      frontier = next;
    }
  }
  return {start, nodes: [...reached.values()],
    edges: g.edges.filter((e) => reached.has(e.from) && reached.has(e.to))};
}

/** One line of the indented tree. */
export interface TreeLine {
  node: ReachedNode;
  level: number;
}

/** The walk as an indented tree for one direction, children in the order they were reached. */
export function provenanceTree(p: Provenance, direction: 'up' | 'down'): TreeLine[] {
  const children = new Map<string, ReachedNode[]>();
  for (const n of p.nodes.filter((x) => x.direction === direction))
    children.set(n.via!, [...children.get(n.via!) ?? [], n]);
  const out: TreeLine[] = [];
  const walk = (key: string, level: number) => {
    for (const c of children.get(key) ?? []) {
      out.push({node: c, level});
      walk(c.key, level + 1);
    }
  };
  walk(p.start, 0);
  return out;
}

/** `Run: Feed 4 (Fed-batch titer)`, a node as one line. */
export const describeNode = (n: Pick<ProvenanceNode, 'kind' | 'name' | 'detail'>): string =>
  `${PROVENANCE_LABELS[n.kind]}: ${n.name}${n.detail ? ` (${n.detail})` : ''}`;

/** The tree as text, for reports: `Rests on` and `Used by`, two spaces per level, each line with the
 * relation that links it to the line above it: `Run: Feed 4 · evidence for`. */
export function summarizeProvenance(p: Provenance): string {
  const start = p.nodes.find((n) => n.key === p.start);
  if (!start)
    return 'Nothing found';
  const lines = [describeNode(start)];
  for (const [title, dir] of [['Rests on', 'up'], ['Used by', 'down']] as const) {
    const tree = provenanceTree(p, dir);
    lines.push(`${title}:${tree.length ? '' : ' nothing'}`,
      ...tree.map((t) => `${'  '.repeat(t.level + 1)}${describeNode(t.node)} · ${t.node.relation}`));
  }
  return lines.join('\n');
}
