/** Decision records where the registry owns the stage, and searching decisions.
 *
 * In a portfolio whose stages come from its registry, the registry stays the authority: a gate
 * decision is taken there and a sync moves the program. What the app adds is the record: the
 * outcome, the rationale and the evidence the gate review rested on, kept with the program. A
 * record is a gate-log entry (`stage_transition`) with decision `record`: from the program's stage,
 * to the stage the gate leads into, and the program stays where it is. It never moves the program
 * and never changes its status.
 *
 * Pure: the caller reads the program, its stages and the user's access. */

export type DecisionOutcome = 'go' | 'hold' | 'no_go';
export const DECISION_OUTCOMES: DecisionOutcome[] = ['go', 'hold', 'no_go'];
export const OUTCOME_LABELS: Record<DecisionOutcome, string> = {go: 'Go', hold: 'Hold', no_go: 'No-go'};

/** What the page says beside a record, so nobody reads it as a stage change. */
export const RECORD_NOTE = 'decision recorded; stage set by registry';

export interface RecordCheck {
  stageAuthority: string | null | undefined;
  isApprover: boolean;
  programStatus: string;
  outcome: string;
  rationale: string;
  /** The stage the gate leads into; null when there is none. */
  gateStageId: string | null;
  stages: {id: string, ordinal: number, isTerminal?: boolean | null}[];
  currentStageId: string | null;
}

/** Why the decision record may not be written; empty when it may. */
export function checkRecord(c: RecordCheck): string[] {
  const reasons: string[] = [];
  if (c.stageAuthority !== 'feed')
    reasons.push('Decisions here move the stage: record them as gate decisions');
  if (!c.isApprover)
    reasons.push('Only the portfolio\'s approvers record gate decisions');
  if (!DECISION_OUTCOMES.includes(c.outcome as DecisionOutcome))
    reasons.push(`The outcome is one of ${DECISION_OUTCOMES.map((o) => OUTCOME_LABELS[o]).join(', ')}`);
  if (!c.rationale.trim())
    reasons.push('A decision record needs a rationale');
  if (['terminated', 'completed', 'transitioned'].includes(c.programStatus))
    reasons.push(`The program is ${c.programStatus}`);
  if (c.gateStageId === null)
    reasons.push('There is no gate after the program\'s stage');
  else if (!c.stages.some((s) => s.id === c.gateStageId))
    reasons.push('The gate belongs to a different stage model');
  return reasons;
}

/** One line for a record: "Go at the gate to Phase 3 — rationale (decision recorded; stage set by registry)". */
export function describeRecord(outcome: string | null | undefined, gateName: string, rationale?: string | null): string {
  const label = OUTCOME_LABELS[outcome as DecisionOutcome] ?? 'Decision';
  return `${label} at the gate to ${gateName}${rationale ? ` — ${rationale}` : ''} (${RECORD_NOTE})`;
}

// ---------------------------------------------------------------------------------------------
// Searching decisions

export interface DecisionView {
  id: string;
  portfolio: string;
  programId: string;
  programCode: string;
  programName: string;
  decision: string;
  outcome: string | null;
  /** Stage codes and names; for a record, `to` is the gate's stage. */
  from: {code: string, name: string} | null;
  to: {code: string, name: string} | null;
  on: string;
  by: string | null;
  source: string | null;
  rationale: string | null;
  evidence: {kind: string, name: string, target: string}[];
}

export interface DecisionQuery {
  /** Words that must all appear in the decision, its rationale, its stages or its evidence. */
  text?: string;
  portfolio?: string;
  /** A program's code, or part of its name. */
  program?: string;
  /** A gate: the code or name of the stage it leads into (or comes from). */
  gate?: string;
  /** Leave out the registry's own moves (sync rows). */
  withoutSync?: boolean;
  /** Only decisions with evidence. */
  withEvidence?: boolean;
}

const norm = (s: string | null | undefined): string => String(s ?? '').toLowerCase().replace(/[\s_\-–—·.]+/g, ' ').trim();

/** The text a decision is found by. */
export function decisionText(d: DecisionView): string {
  return norm([d.decision, d.outcome ? OUTCOME_LABELS[d.outcome as DecisionOutcome] ?? d.outcome : '',
    d.decision === 'record' ? RECORD_NOTE : '', d.rationale, d.programCode, d.programName,
    d.from?.code, d.from?.name, d.to?.code, d.to?.name, d.source,
    ...d.evidence.flatMap((e) => [e.kind, e.name, e.target])].join(' '));
}

/** The decisions [q] finds, newest first. */
export function searchDecisions(decisions: DecisionView[], q: DecisionQuery): DecisionView[] {
  const words = norm(q.text).split(' ').filter((w) => w.length > 0);
  const program = norm(q.program);
  const gate = norm(q.gate);
  return decisions.filter((d) => {
    if (q.portfolio && d.portfolio !== q.portfolio)
      return false;
    if (q.withoutSync && d.decision === 'sync')
      return false;
    if (q.withEvidence && !d.evidence.length)
      return false;
    if (program && norm(d.programCode) !== program && !norm(d.programName).includes(program))
      return false;
    if (gate && ![d.to, d.from].some((s) => s && (norm(s.code) === gate || norm(s.name).includes(gate))))
      return false;
    if (words.length) {
      const text = decisionText(d);
      if (!words.every((w) => text.includes(w)))
        return false;
    }
    return true;
  }).sort((a, b) => b.on.localeCompare(a.on));
}

/** One line for a decision in the search results. */
export function describeDecision(d: DecisionView): string {
  const when = d.on.slice(0, 10);
  if (d.decision === 'record')
    return `${when} · ${d.programCode}: ${describeRecord(d.outcome, d.to?.name ?? '?', d.rationale)}`;
  const how = d.decision === 'sync' ? ' (registry)' : d.source === 'data' ? ' (by data)' : '';
  return `${when} · ${d.programCode}: ${d.decision} ${d.from?.name ?? '—'} → ${d.to?.name ?? '—'}${how}` +
    `${d.rationale ? ` — ${d.rationale}` : ''}`;
}
