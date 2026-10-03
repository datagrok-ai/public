/** A program's data package manifest (F-43, FR-43, PORT-11), started as a view: one line per thing
 * the program's gate case rests on, read from the portal's records each time, nothing stored.
 * - criteria of every stage up to the next gate, with their state and whether it was recorded or
 *   evaluated from data;
 * - documents, with version, status, who linked and approved them;
 * - generated evidence tables, with what they were written from and whether that evidence changed;
 * - saved model runs, with model, unit operation, run type, who ran them and the criteria they were
 *   selected for;
 * - the program's datasets (where its data comes from) and its snapshots (the data as frozen).
 * The same lines export to CSV. An authored manifest (a table people edit) is later.
 * Pure: the caller reads the records. */

export type ManifestSection = 'criterion' | 'document' | 'evidence_table' | 'run' | 'dataset' | 'snapshot';

export const MANIFEST_SECTIONS: ManifestSection[] = ['criterion', 'document', 'evidence_table', 'run', 'dataset', 'snapshot'];

export const MANIFEST_LABELS: Record<ManifestSection, string> = {
  criterion: 'Criterion', document: 'Document', evidence_table: 'Evidence table', run: 'Run', dataset: 'Dataset',
  snapshot: 'Snapshot',
};

export interface ManifestLine {
  section: ManifestSection;
  name: string;
  /** The stage it belongs to: a criterion's stage, the stage a run was saved at. */
  stage: string | null;
  /** The criterion it answers or was selected for, by code; several '|'-separated. */
  criterion: string | null;
  status: string;
  version: string | null;
  /** Where it comes from, in a few words: `recorded`, `model X · unit op · run type`, a connection. */
  source: string;
  by: string | null;
  /** yyyy-mm-dd. */
  on: string | null;
  /** What opens it: a file path, a saved call's id, a snapshot id. */
  ref: string;
}

export interface ManifestInput {
  stages: {id: string, code: string, ordinal: number}[];
  /** The program's stage; null when unstaged. */
  stageId: string | null;
  criteria: {id: string, stageId: string, code: string, name: string, kind: string, required: boolean}[];
  /** Criterion states by criterion id, and which were evaluated from data now. */
  states: Record<string, string>;
  evaluated: Set<string>;
  documents: {id: string, title: string, kind: string, version?: string | null, status: string, path: string,
    criterionId?: string | null, linkedBy?: string | null, linkedOn?: string | null, approvedBy?: string | null,
    approvedOn?: string | null, writtenBy?: string | null,
    /** For a generated evidence table: the kinds of evidence it was written from, and its last check. */
    generated?: {sources: string[], state?: 'current' | 'stale' | null} | null}[];
  runs: {id: string, runId?: string | null, title: string, model: string, modelVersion?: string | null,
    stageId?: string | null, unitOperation?: string | null, runType?: string | null, purpose?: string | null,
    by?: string | null, on?: string | null, isReference?: boolean | null,
    /** Criterion ids the run is the selected run for. */
    selectedFor: string[]}[];
  datasets: {code: string, name: string, connection: string}[];
  snapshots: {id: string, label: string, by?: string | null, on?: string | null}[];
  /** User display names by id. */
  names?: Map<string, string>;
}

const day = (s: string | null | undefined) => s ? String(s).slice(0, 10) : null;

/** The criteria a package covers: every stage the program has reached, and the next one (its gate).
 * An unstaged program: the first stage. */
export function manifestStages(stages: ManifestInput['stages'], stageId: string | null): Set<string> {
  const sorted = [...stages].sort((a, b) => a.ordinal - b.ordinal);
  const at = sorted.findIndex((s) => s.id === stageId);
  return new Set(sorted.slice(0, at + 2).map((s) => s.id));
}

/** The manifest's lines, section by section in MANIFEST_SECTIONS order; criteria in stage order,
 * the rest newest first as given. */
export function buildManifest(x: ManifestInput): ManifestLine[] {
  const name = (id: string | null | undefined) => id ? x.names?.get(id) ?? id : null;
  const stageCode = new Map(x.stages.map((s) => [s.id, s.code]));
  const ordinal = new Map(x.stages.map((s) => [s.id, s.ordinal]));
  const criterionCode = new Map(x.criteria.map((c) => [c.id, c.code]));
  const covered = manifestStages(x.stages, x.stageId);
  const out: ManifestLine[] = [];

  const criteria = x.criteria.filter((c) => covered.has(c.stageId))
    .sort((a, b) => (ordinal.get(a.stageId) ?? 0) - (ordinal.get(b.stageId) ?? 0));
  for (const c of criteria) {
    out.push({section: 'criterion', name: c.name, stage: stageCode.get(c.stageId) ?? null, criterion: c.code,
      status: (x.states[c.id] ?? 'not_started').replace(/_/g, ' '), version: null,
      source: `${c.kind.replace(/_/g, ' ')}${c.required ? '' : ', optional'} · ${x.evaluated.has(c.id) ? 'evaluated from data now' : 'recorded'}`,
      by: null, on: null, ref: c.code});
  }
  for (const d of x.documents) {
    const g = d.generated;
    out.push({section: g ? 'evidence_table' : 'document', name: d.title,
      stage: null, criterion: d.criterionId ? criterionCode.get(d.criterionId) ?? null : null,
      status: g?.state === 'stale' ? `${d.status}, evidence changed` : d.status, version: d.version ?? null,
      source: g ? `generated from ${g.sources.length ? [...new Set(g.sources)].join(', ') : 'nothing recorded'}` :
        `${d.kind.replace(/_/g, ' ')}${d.writtenBy ? ` · written by ${name(d.writtenBy)}` : ''}`,
      by: name(d.approvedBy ?? d.linkedBy), on: day(d.approvedOn ?? d.linkedOn), ref: d.path});
  }
  for (const r of x.runs) {
    const selected = r.selectedFor.map((id) => criterionCode.get(id)).filter((c): c is string => !!c);
    out.push({section: 'run', name: r.title, stage: r.stageId ? stageCode.get(r.stageId) ?? null : null,
      criterion: selected.length ? selected.join('|') : null,
      status: selected.length ? 'selected' : r.isReference ? 'reference' : 'saved', version: r.modelVersion ?? null,
      source: [r.model, r.unitOperation, r.runType, r.purpose].filter((s) => s).join(' · '),
      by: name(r.by), on: day(r.on), ref: r.runId ?? r.id});
  }
  for (const d of x.datasets) {
    out.push({section: 'dataset', name: d.name, stage: null, criterion: null, status: 'live', version: null,
      source: d.connection, by: null, on: null, ref: d.code});
  }
  for (const s of x.snapshots) {
    out.push({section: 'snapshot', name: s.label, stage: null, criterion: null, status: 'frozen', version: null,
      source: 'program snapshot', by: name(s.by), on: day(s.on), ref: s.id});
  }
  return out;
}

/** `3 criteria · 2 documents · 1 run`: the sections that have lines. */
export function summarizeManifest(lines: ManifestLine[]): string {
  const plural: Record<ManifestSection, [string, string]> = {criterion: ['criterion', 'criteria'],
    document: ['document', 'documents'], evidence_table: ['evidence table', 'evidence tables'], run: ['run', 'runs'],
    dataset: ['dataset', 'datasets'], snapshot: ['snapshot', 'snapshots']};
  const parts = MANIFEST_SECTIONS.map((s) => [s, lines.filter((l) => l.section === s).length] as const)
    .filter(([, n]) => n > 0).map(([s, n]) => `${n} ${plural[s][n === 1 ? 0 : 1]}`);
  return parts.length ? parts.join(' · ') : 'Nothing recorded yet';
}

export const MANIFEST_COLUMNS: {key: keyof ManifestLine, name: string}[] = [
  {key: 'section', name: 'Section'}, {key: 'name', name: 'Item'}, {key: 'stage', name: 'Stage'},
  {key: 'criterion', name: 'Criterion'}, {key: 'status', name: 'Status'}, {key: 'version', name: 'Version'},
  {key: 'source', name: 'Source'}, {key: 'by', name: 'By'}, {key: 'on', name: 'On'}, {key: 'ref', name: 'Reference'},
];

/** A cell as the manifest shows it: the section by its label. */
export const manifestCell = (l: ManifestLine, key: keyof ManifestLine): string =>
  key === 'section' ? MANIFEST_LABELS[l.section] : String(l[key] ?? '');

const csvField = (v: string) => /[",\n\r]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v;

/** The manifest as CSV (RFC 4180): a header, then a line each. */
export function manifestCsv(lines: ManifestLine[]): string {
  return [MANIFEST_COLUMNS.map((c) => c.name), ...lines.map((l) => MANIFEST_COLUMNS.map((c) => manifestCell(l, c.key)))]
    .map((row) => row.map(csvField).join(',')).join('\n') + '\n';
}

/** The download's file name: `EXB-201-data-package-2026-10-03.csv`. */
export const manifestFileName = (programCode: string, on: Date): string =>
  `${programCode.replace(/[^A-Za-z0-9_.-]+/g, '_')}-data-package-${on.toISOString().slice(0, 10)}.csv`;
