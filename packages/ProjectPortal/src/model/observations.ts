/** Observations (F-26, FR-76/106): a curated conclusion on a program, one of its series, campaigns
 * or studies, or a compound it holds. "KIN-A: para-substitution buys potency, not solubility."
 * Written by a person, with the evidence it rests on; never edited, retracted instead, so what the
 * team concluded at the time stays readable. Pure: the caller reads the program's units and
 * compounds and says who the user is. */
import type {EvidenceItem} from './evidence';
import {validateEvidence} from './evidence';

export type ObservationSubject = 'program' | 'series' | 'compound' | 'campaign' | 'study';
export const OBSERVATION_SUBJECTS: ObservationSubject[] = ['program', 'series', 'compound', 'campaign', 'study'];

/** Evidence an observation may cite: what exists outside the gate (no criterion). */
export const OBSERVATION_EVIDENCE = ['model_run', 'snapshot', 'file', 'url'] as const;

export interface ObservationInput {
  subject: ObservationSubject;
  /** The unit's or compound's code; empty for the program. */
  code?: string | null;
  title: string;
  finding: string;
  evidence?: Pick<EvidenceItem, 'kind' | 'target' | 'name'>[];
}

/** What the program has that an observation can be about, by kind. */
export type ObservationTargets = Partial<Record<Exclude<ObservationSubject, 'program'>, string[]>>;

/** One subject the dialog offers: `series KIN-A`. */
export interface SubjectOption {
  subject: ObservationSubject;
  code: string | null;
  label: string;
}

export function subjectLabel(subject: ObservationSubject, code?: string | null): string {
  return subject === 'program' || !code ? 'the program' : `${subject} ${code}`;
}

/** The program first, then its series, compounds, campaigns and studies, each sorted. */
export function subjectOptions(targets: ObservationTargets): SubjectOption[] {
  const out: SubjectOption[] = [{subject: 'program', code: null, label: 'The program'}];
  for (const kind of ['series', 'compound', 'campaign', 'study'] as const) {
    for (const code of [...new Set(targets[kind] ?? [])].sort())
      out.push({subject: kind, code, label: subjectLabel(kind, code)});
  }
  return out;
}

/** Problems with an observation before it is written, as sentences. [targets] are the program's
 * own units and compounds: an observation is about something the program holds. */
export function validateObservation(o: Partial<ObservationInput>, targets: ObservationTargets): string[] {
  const errors: string[] = [];
  if (!o.subject || !OBSERVATION_SUBJECTS.includes(o.subject))
    errors.push(`an observation is about one of ${OBSERVATION_SUBJECTS.join(', ')}`);
  else if (o.subject !== 'program') {
    if (!o.code?.trim())
      errors.push(`an observation on a ${o.subject} names it`);
    else if (!(targets[o.subject] ?? []).includes(o.code.trim()))
      errors.push(`the program has no ${o.subject} ${o.code.trim()}`);
  }
  if (!o.title?.trim())
    errors.push('an observation needs a title');
  if (!o.finding?.trim())
    errors.push('an observation needs the finding itself');
  for (const e of o.evidence ?? []) {
    if (!(OBSERVATION_EVIDENCE as readonly string[]).includes(e.kind))
      errors.push(`evidence of an observation is one of ${OBSERVATION_EVIDENCE.join(', ')}, not ${e.kind}`);
    else
      errors.push(...validateEvidence({...e, role: 'supporting'}));
  }
  return errors;
}

/** Why the user may not add an observation (empty when they may): contributors and approvers do,
 * also where the registry owns the stage. An ended program still takes observations: a conclusion
 * is often drawn after the end. */
export function checkObservationAccess(a: {isApprover: boolean, isContributor: boolean}): string[] {
  return a.isApprover || a.isContributor ? [] : ['Only the program\'s contributors and approvers record observations'];
}

/** Why the user may not retract one: its author or an approver may. */
export function checkRetraction(o: {status: string, by: string | null}, user: string, isApprover: boolean, reason: string): string[] {
  const errors: string[] = [];
  if (o.status === 'retracted')
    errors.push('It is retracted already');
  if (o.by !== user && !isApprover)
    errors.push('Only its author or an approver retracts an observation');
  if (!reason.trim())
    errors.push('Retracting needs a reason');
  return errors;
}

/** An observation as the page and the feed show it. */
export interface ObservationView {
  id: string;
  subject: ObservationSubject;
  code: string | null;
  title: string;
  finding: string;
  evidence: Pick<EvidenceItem, 'kind' | 'target' | 'name'>[];
  status: 'current' | 'retracted';
  by: string | null;
  /** ISO time. */
  on: string | null;
  retractedBy?: string | null;
  retractedOn?: string | null;
  retraction?: string | null;
}

/** Current ones first, each group newest first. */
export function sortObservations(list: ObservationView[]): ObservationView[] {
  const t = (s: string | null | undefined) => s ? Date.parse(s) : 0;
  return [...list].sort((a, b) => (a.status === 'current' ? 0 : 1) - (b.status === 'current' ? 0 : 1) || t(b.on) - t(a.on));
}

/** The ones about [subject] [code], or about anything when no subject is given. */
export function observationsAbout(list: ObservationView[], subject?: ObservationSubject | null, code?: string | null): ObservationView[] {
  return list.filter((o) => !subject || (o.subject === subject && (subject === 'program' || o.code === code)));
}

/** "series KIN-A · Para-substitution buys potency". */
export function describeObservation(o: ObservationView): string {
  return `${subjectLabel(o.subject, o.code)} · ${o.title}${o.status === 'retracted' ? ' (retracted)' : ''}`;
}

/** The evidence as stored (`{items: [...]}`, a json column holds a map): kind, target and name. */
export function evidenceOfObservation(raw: unknown): ObservationView['evidence'] {
  const items = Array.isArray(raw) ? raw : Array.isArray((raw as any)?.items) ? (raw as any).items : [];
  return (items as any[]).filter((e) => e && typeof e === 'object' && e.kind && e.target)
    .map((e: any) => ({kind: e.kind, target: String(e.target), name: String(e.name ?? e.target)}));
}
