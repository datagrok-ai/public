/** Decision records with evidence where the registry owns the stage, and the decision search, against
 * the platform. The rules are in ./model/decision-records. A record is written like a gate decision
 * (gate-log entry plus evidence) but leaves the program's stage and status alone: the registry stays
 * the authority. Writes need what any gate-log write needs on the server (Edit on the program row). */
import * as DG from 'datagrok-api/dg';
import {portalDb} from './generated/db';
import {programContext} from './portal';
import {suggestedEvidence, writeEvidence} from './evidence';
import {EvidenceItem} from './model/evidence';
import {nextStage} from './model/gates';
import {summaryOf, takeSnapshot} from './snapshots';
import {snapshotEvidence} from './model/snapshots';
import {checkRecord, DecisionOutcome, DecisionView, describeRecord} from './model/decision-records';
import {queryIn, ROW_LIMIT} from './rows';

export interface RecordOptions {
  outcome: DecisionOutcome;
  rationale: string;
  /** The stage the gate leads into; the next stage by default. */
  gateStageId?: string | null;
  /** What the record rests on; the gate's suggestion (readiness, selected runs, documents) by default. */
  items?: EvidenceItem[];
  /** Also take a snapshot now and attach it. */
  snapshot?: boolean;
}

export interface RecordedRecord {
  transitionId: string;
  line: string;
  evidence: {written: number, errors: string[]};
}

/** Records a gate decision with its evidence in a portfolio whose registry owns the stage. The
 * program does not move. The portfolio's or the program's approvers only. */
export async function recordGateDecision(programId: string, o: RecordOptions): Promise<RecordedRecord> {
  const ctx = await programContext(programId);
  const next = nextStage(ctx.stages.map((s) => ({...s, isTerminal: s.is_terminal})), ctx.program.stage_id ?? null);
  const gateStageId = o.gateStageId ?? next?.id ?? null;
  const refused = checkRecord({stageAuthority: ctx.portfolio.stage_authority, isApprover: ctx.approver,
    programStatus: ctx.program.status ?? 'active', outcome: o.outcome, rationale: o.rationale ?? '', gateStageId,
    stages: ctx.stages, currentStageId: ctx.program.stage_id ?? null});
  if (refused.length)
    throw new Error(refused.join('; '));
  const gate = ctx.stages.find((s) => s.id === gateStageId)!;
  let items = o.items ?? (gateStageId === next?.id ? await suggestedEvidence(programId, ctx.readiness, gate.name).catch(() => []) : []);
  if (o.snapshot) {
    const snap = await takeSnapshot(ctx.program, {reason: `at the decision recorded for the gate to ${gate.name}`});
    items = [...items, snapshotEvidence(snap.snapshot.id, summaryOf(snap.snapshot))];
  }
  const [t] = await portalDb.stageTransitions.insert({program_id: programId, from_stage_id: ctx.program.stage_id ?? undefined,
    to_stage_id: gate.id, decision: 'record', outcome: o.outcome, decided_by: DG.User.current().id,
    decided_on: new Date().toISOString() as any, rationale: o.rationale.trim(), source: 'approver'});
  const evidence = items.length ? await writeEvidence(programId, t.id, items)
    .catch((e: any) => ({written: 0, errors: [`evidence not written: ${e?.message ?? e}`]})) : {written: 0, errors: []};
  return {transitionId: t.id, evidence,
    line: `${ctx.program.code}: ${describeRecord(o.outcome, gate.name)}` +
      (evidence.written ? `, with ${evidence.written} piece(s) of evidence` : '')};
}

const iso = (v: any): string => {
  if (!v)
    return '';
  const d = typeof v?.toDate === 'function' ? v.toDate() : new Date(v);
  return Number.isNaN(d.getTime()) ? '' : d.toISOString();
};

/** Every decision the user can read, with its evidence, in [portfolio] or every portfolio they can
 * see. What the server does not show them (another portfolio's programs) is never read. */
export async function decisionsOf(portfolio?: string | null): Promise<DecisionView[]> {
  const portfolios = (await portalDb.portfolios.query({limit: ROW_LIMIT})).filter((p) => !portfolio || p.code === portfolio);
  if (!portfolios.length)
    return [];
  const programs = await queryIn(portalDb.programs, 'portfolio_id', portfolios.map((p) => p.id));
  const ids = programs.map((p) => p.id);
  const [transitions, evidence, stages] = await Promise.all([queryIn(portalDb.stageTransitions, 'program_id', ids),
    queryIn(portalDb.evidences, 'program_id', ids),
    queryIn(portalDb.stages, 'stage_model_id', programs.map((p) => p.stage_model_id))]);
  const program = new Map(programs.map((p) => [p.id, p]));
  const code = new Map(portfolios.map((p) => [p.id, p.code]));
  const stage = (id?: string | null) => {
    const s = id ? stages.find((x) => x.id === id) : null;
    return s ? {code: s.code, name: s.name} : null;
  };
  return transitions.filter((t) => program.has(t.program_id)).map((t) => {
    const p = program.get(t.program_id)!;
    return {id: t.id, portfolio: code.get(p.portfolio_id) ?? '?', programId: p.id, programCode: p.code, programName: p.name,
      decision: t.decision, outcome: t.outcome ?? null, from: stage(t.from_stage_id), to: stage(t.to_stage_id),
      on: iso(t.decided_on ?? t.created_on), by: t.decided_by ?? null, source: t.source ?? null, rationale: t.rationale ?? null,
      evidence: evidence.filter((e) => e.stage_transition_id === t.id).map((e) => ({kind: e.kind, name: e.name, target: e.target}))};
  });
}
