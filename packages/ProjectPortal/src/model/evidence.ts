/** Decisions with evidence (F-25, FR-71/74). A gate decision is the gate-log entry
 * (`stage_transition`); evidence rows say what it rested on. Each evidence row points at its target
 * by a soft reference (an id or address as text, no foreign key) and freezes the target's name and
 * context when attached, so it reads the same after the target is renamed or deleted.
 *
 * What the gate panel suggests: the readiness of the gate (each required criterion as the gate saw
 * it) and the run selected for each criterion of the stage. The decider may untick any of them and
 * add a link. A decision recorded without a list (by API, or by data) gets the suggestion.
 * Pure: the caller reads readiness, selections and what still exists. */
import type {Readiness} from './readiness';

export type EvidenceKind = 'model_run' | 'snapshot' | 'file' | 'url' | 'criterion';
export const EVIDENCE_KINDS: EvidenceKind[] = ['model_run', 'snapshot', 'file', 'url', 'criterion'];
export const EVIDENCE_LABELS: Record<EvidenceKind, string> = {
  model_run: 'Run', snapshot: 'Snapshot', file: 'File', url: 'Link', criterion: 'Criterion',
};

/** `readiness`: a criterion as the gate saw it; `supporting`: attached by the decider. */
export type EvidenceRole = 'readiness' | 'supporting';

/** One piece of evidence before it is written. */
export interface EvidenceItem {
  kind: EvidenceKind;
  /** The soft reference: a run's run_id, a snapshot's id, a file path, a URL, a gate criterion's id. */
  target: string;
  /** The target's name now, frozen when written. */
  name: string;
  role: EvidenceRole;
  /** What the target says now, frozen when written. */
  context?: Record<string, unknown>;
}

/** Problems with one item, as sentences. */
export function validateEvidence(e: Partial<EvidenceItem>): string[] {
  const errors: string[] = [];
  if (!e.kind || !EVIDENCE_KINDS.includes(e.kind))
    errors.push(`evidence kind "${e.kind ?? ''}" is not one of ${EVIDENCE_KINDS.join(', ')}`);
  if (!e.target?.trim())
    errors.push('evidence needs a target');
  if (!e.name?.trim())
    errors.push('evidence needs a name');
  if (e.kind === 'url' && e.target && !/^https?:\/\/\S+$/i.test(e.target.trim()))
    errors.push(`"${e.target}" is not an http(s) address`);
  if (e.role && !['readiness', 'supporting'].includes(e.role))
    errors.push(`evidence role "${e.role}" is not readiness or supporting`);
  return errors;
}

/** A link typed by the decider: an http(s) address is a URL, anything else a file path
 * (`System:AppData/...`, a folder's file). The name is the last segment of the path, or the host. */
export function linkEvidence(text: string): EvidenceItem {
  const t = text.trim();
  if (/^https?:\/\//i.test(t)) {
    const m = /^https?:\/\/([^/?#]+)([^?#]*)/i.exec(t);
    const last = (m?.[2] ?? '').split('/').filter((s) => s).pop();
    return {kind: 'url', target: t, name: last ? `${m![1]}/…/${last}` : m?.[1] ?? t, role: 'supporting'};
  }
  const name = t.split(/[/\\:]/).filter((s) => s).pop() ?? t;
  return {kind: 'file', target: t, name, role: 'supporting', context: {path: t}};
}

/** The readiness of a gate as evidence: one item per required criterion, with its state, basis and
 * what the data said. Optional criteria are left out. */
export function readinessEvidence(r: Readiness | null, stageName: string): EvidenceItem[] {
  if (!r)
    return [];
  return r.criteria.filter((c) => c.required).map((c) => ({
    kind: 'criterion' as const, target: c.id, name: c.name, role: 'readiness' as const,
    context: {code: c.code, stage: stageName, state: c.state, basis: c.basis,
      ...(c.result ? {summary: c.result.summary} : {}), ...(c.gap ? {gap: c.gap} : {})},
  }));
}

/** A run selected for a criterion of the gate. */
export interface SelectedRun {
  criterionId: string;
  criterionCode: string;
  criterionName: string;
  /** The saved call's id: the soft reference. */
  runId: string;
  title: string;
  model: string;
  selectedOn?: string | null;
  selectedBy?: string | null;
  rationale?: string | null;
}

export function runEvidence(s: SelectedRun): EvidenceItem {
  return {kind: 'model_run', target: s.runId, name: s.title, role: 'supporting',
    context: {model: s.model, criterion: s.criterionCode, criterionName: s.criterionName,
      ...(s.selectedOn ? {selectedOn: s.selectedOn} : {}), ...(s.selectedBy ? {selectedBy: s.selectedBy} : {}),
      ...(s.rationale ? {rationale: s.rationale} : {})}};
}

/** What the gate panel offers ticked: the gate's readiness, then the run selected for each
 * criterion of the gate's stage. */
export function suggestEvidence(readiness: Readiness | null, stageName: string, selected: SelectedRun[]): EvidenceItem[] {
  const ids = new Set((readiness?.criteria ?? []).map((c) => c.id));
  return [...readinessEvidence(readiness, stageName),
    ...selected.filter((s) => ids.has(s.criterionId)).map(runEvidence)];
}

/** A stable key for an item, so a list holds each target once. */
export const evidenceKey = (e: Pick<EvidenceItem, 'kind' | 'target'>) => `${e.kind}:${e.target}`;

export function uniqueEvidence(items: EvidenceItem[]): EvidenceItem[] {
  const seen = new Set<string>();
  return items.filter((e) => !seen.has(evidenceKey(e)) && !!seen.add(evidenceKey(e)));
}

/** Whether a target still resolves, and under which name: `live` (same name), `renamed`, `gone`,
 * or `unchecked` (a URL is not fetched; a target that could not be looked up). */
export type TargetState = 'live' | 'renamed' | 'gone' | 'unchecked';

/** [current] is the target's name now: null when it is gone, undefined when it was not looked up. */
export function targetState(kind: EvidenceKind, frozenName: string, current: string | null | undefined): TargetState {
  if (kind === 'url' || current === undefined)
    return 'unchecked';
  if (current === null)
    return 'gone';
  return current === frozenName ? 'live' : 'renamed';
}

/** One line for the history: `Run: Titer at 3% feed (renamed: Titer v2)`, `Criterion: Assay ready, met`,
 * `File: Process description · v1.0 approved` for a document. */
export function describeEvidence(e: {kind: EvidenceKind, name: string, context?: Record<string, any> | null},
  state: TargetState = 'live', current?: string | null): string {
  const c = e.context ?? {};
  const detail = e.kind === 'criterion' ? [String(c.state ?? '').replace(/_/g, ' '), c.summary].filter((x) => x).join(': ') :
    e.kind === 'model_run' ? [c.criterion ? `for ${c.criterion}` : '', c.model].filter((x) => x).join(', ') :
      e.kind === 'snapshot' ? String(c.takenOn ?? '').slice(0, 16).replace('T', ' ') :
        e.kind === 'file' && c.document ? [c.version ? `v${c.version}` : '', c.status].filter((x) => x).join(' ') : '';
  const tail = state === 'gone' ? ' (gone; as recorded)' : state === 'renamed' ? ` (now ${current})` : '';
  return `${EVIDENCE_LABELS[e.kind]}: ${e.name}${detail ? ` · ${detail}` : ''}${tail}`;
}

/** The evidence of one decision in the order the history shows it: supporting first, then readiness. */
export function orderEvidence<T extends {role?: string | null, kind: string, name: string}>(items: T[]): T[] {
  const rank = (e: T) => (e.role === 'readiness' ? 1 : 0);
  return [...items].sort((a, b) => rank(a) - rank(b) || a.kind.localeCompare(b.kind) || a.name.localeCompare(b.name));
}
