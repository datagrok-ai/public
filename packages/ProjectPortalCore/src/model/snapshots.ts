/** A dated program snapshot (F-32, FR-15): the program's stage and status, the state of every
 * criterion that applies to it, the runs selected for them, and the program dataset as a file,
 * frozen at one time. A snapshot is written once; its table lives in the program's workspace
 * folder with its SHA-256 on the record, so a file changed afterwards is noticed. A decision may
 * take one as evidence. Pure: the caller reads the context and writes the file. */
import type {EvidenceItem, SelectedRun} from './evidence';
import type {Readiness} from './readiness';

/** One criterion as the snapshot froze it. */
export interface SnapshotCriterion {
  code: string;
  name: string;
  stage: string;
  kind: string;
  required: boolean;
  state: string;
  /** `evaluated` for a data condition of the next gate computed at the time; else `recorded`. */
  basis: 'evaluated' | 'recorded';
  summary?: string;
}

export interface SnapshotCriterionInput {
  id: string;
  stageId: string;
  code: string;
  name: string;
  kind: string;
  required: boolean;
}

/** Every criterion that applies to the program, in stage order then declaration order, with its
 * state now: recorded, or for the next gate's data conditions as evaluated from the data. */
export function snapshotCriteria(criteria: SnapshotCriterionInput[], stages: {id: string, code: string, ordinal: number}[],
  recorded: Record<string, string>, readiness: Readiness | null): SnapshotCriterion[] {
  const order = new Map(stages.map((s) => [s.id, s.ordinal]));
  const code = new Map(stages.map((s) => [s.id, s.code]));
  const evaluated = new Map((readiness?.criteria ?? []).map((c) => [c.id, c]));
  return criteria.map((c, i) => ({c, i}))
    .sort((a, b) => (order.get(a.c.stageId) ?? 0) - (order.get(b.c.stageId) ?? 0) || a.i - b.i)
    .map(({c}) => {
      const r = evaluated.get(c.id);
      return {code: c.code, name: c.name, stage: code.get(c.stageId) ?? '', kind: c.kind, required: c.required,
        state: r ? r.state : recorded[c.id] ?? 'not_started', basis: r?.basis ?? 'recorded',
        ...(r?.result ? {summary: r.result.summary} : {})};
    });
}

/** The file a snapshot's table is written to in the program's folder: `snapshot-EXS-KIN1-2026-10-02T103105Z.csv`. */
export function snapshotFileName(programCode: string, takenOn: string): string {
  const stamp = takenOn.slice(0, 19).replace(/[-:]/g, '').replace(/^(\d{8})T/, (_, d) =>
    `${d.slice(0, 4)}-${d.slice(4, 6)}-${d.slice(6, 8)}T`);
  return `snapshot-${programCode.replace(/[^A-Za-z0-9_-]+/g, '_')}-${stamp}Z.csv`;
}

/** Whether the table file still holds what was frozen. */
export type TableCheck = 'intact' | 'changed' | 'missing' | 'not kept';

export function checkTable(recordedSha: string | null | undefined, actualSha: string | null | undefined, kept: boolean): TableCheck {
  if (!kept)
    return 'not kept';
  if (!actualSha)
    return 'missing';
  return actualSha === recordedSha ? 'intact' : 'changed';
}

/** What the program page and the dialog say about one snapshot. */
export interface SnapshotSummary {
  takenOn: string;
  label: string;
  stageName: string | null;
  criteria: SnapshotCriterion[];
  tableRows: number | null;
}

export const snapshotDate = (iso: string) => iso.slice(0, 16).replace('T', ' ');

/** `2026-10-02 10:31 · Target Validation · 1 of 1 required met or waived · 12 rows` */
export function describeSnapshot(s: SnapshotSummary): string {
  const required = s.criteria.filter((c) => c.required);
  const done = required.filter((c) => c.state === 'met' || c.state === 'waived').length;
  return [snapshotDate(s.takenOn), s.stageName ?? 'unstaged',
    required.length ? `${done} of ${required.length} required met or waived` : 'no required criteria',
    s.tableRows === null ? 'no table' : `${s.tableRows} row${s.tableRows === 1 ? '' : 's'}`].join(' · ');
}

/** A snapshot as evidence for a decision. */
export function snapshotEvidence(id: string, s: SnapshotSummary): EvidenceItem {
  return {kind: 'snapshot', target: id, name: s.label, role: 'supporting',
    context: {takenOn: s.takenOn, stage: s.stageName ?? 'unstaged', rows: s.tableRows}};
}

/** The label a snapshot gets when none is given. */
export const defaultSnapshotLabel = (programCode: string, takenOn: string, reason?: string) =>
  `${programCode} ${snapshotDate(takenOn)}${reason ? ` (${reason})` : ''}`;

/** The selected runs as the snapshot keeps them. */
export const snapshotArtifacts = (runs: SelectedRun[]) => runs.map((r) => ({criterion: r.criterionCode, runId: r.runId,
  title: r.title, model: r.model, selectedOn: r.selectedOn ?? null, selectedBy: r.selectedBy ?? null, rationale: r.rationale ?? null}));
