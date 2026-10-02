/** Decisions with evidence against the platform: what the gate panel suggests (readiness and the
 * runs selected for the gate's criteria), writing evidence beside a gate-log entry, and reading it
 * back with each target looked up again, so the history can say what was renamed or is gone.
 * The rules live in ./model/evidence. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import dayjs from 'dayjs';
import {EvidenceRow, portalDb, ProgramRow, StageTransitionRow} from './generated/db';
import {programContext} from './portal';
import {
  describeEvidence, EvidenceItem, EvidenceKind, linkEvidence, orderEvidence, runEvidence, SelectedRun, suggestEvidence,
  targetState, TargetState, uniqueEvidence, validateEvidence,
} from './model/evidence';
import {Readiness} from './model/readiness';

const inList = (ids: string[]) => ids.map((id) => `"${id}"`).join(', ');

/** The run selected now for each criterion of the program, with the names the evidence freezes. */
export async function selectedRunsOf(programId: string): Promise<SelectedRun[]> {
  const picks = (await portalDb.criterionArtifacts.query({filter: `program_id = "${programId}"`}))
    .filter((a) => a.role === 'selected');
  if (!picks.length)
    return [];
  const [statuses, runs] = await Promise.all([
    portalDb.criterionStatuses.query({filter: `id in (${inList([...new Set(picks.map((p) => p.criterion_status_id))])})`}),
    portalDb.modelRuns.query({filter: `id in (${inList([...new Set(picks.map((p) => p.model_run_id))])})`}),
  ]);
  const criteria = await portalDb.gateCriterions.query(
    {filter: `id in (${inList([...new Set(statuses.map((s) => s.gate_criterion_id))])})`});
  const users = new Map<string, string>();
  for (const id of new Set(picks.map((p) => p.selected_by).filter((x): x is string => !!x)))
    users.set(id, (await grok.dapi.users.find(id).catch(() => null))?.friendlyName ?? '');
  const out: SelectedRun[] = [];
  for (const p of picks) {
    const status = statuses.find((s) => s.id === p.criterion_status_id);
    const c = criteria.find((x) => x.id === status?.gate_criterion_id);
    const run = runs.find((r) => r.id === p.model_run_id);
    if (!c || !run)
      continue;
    out.push({criterionId: c.id, criterionCode: c.code, criterionName: c.name, runId: run.run_id ?? run.id,
      title: run.title ?? run.run_id ?? run.model_name, model: run.model_name,
      selectedOn: p.selected_on ? dayjs(p.selected_on as any).toISOString() : null,
      selectedBy: p.selected_by ? users.get(p.selected_by) || null : null, rationale: p.rationale ?? null});
  }
  return out;
}

/** What the gate panel offers ticked for [readiness]'s stage: its readiness, and the runs selected
 * for its criteria. */
export async function suggestedEvidence(programId: string, readiness: Readiness | null, stageName: string):
  Promise<EvidenceItem[]> {
  return suggestEvidence(readiness, stageName, await selectedRunsOf(programId).catch(() => []));
}

/** Writes [items] beside the gate-log entry [transitionId]. Every item is checked first; one that
 * fails is reported, not written. Returns what was written and the refusals. */
export async function writeEvidence(programId: string, transitionId: string, items: EvidenceItem[]):
  Promise<{written: number, errors: string[]}> {
  const errors: string[] = [];
  const now = dayjs().toISOString();
  const good = uniqueEvidence(items).filter((e) => {
    const problems = validateEvidence(e);
    errors.push(...problems.map((p) => `${e.name || e.target}: ${p}`));
    return problems.length === 0;
  });
  if (good.length) {
    await portalDb.evidences.insert(good.map((e) => ({program_id: programId, stage_transition_id: transitionId,
      kind: e.kind, target: e.target.trim(), name: e.name.trim(), role: e.role, context: e.context ?? {},
      added_by: DG.User.current().id, added_on: now as any})));
  }
  return {written: good.length, errors};
}

/** The current name of each target, by evidence id: null when it is gone; absent when it cannot be
 * looked up (a URL, a file the caller cannot reach). One query per kind. */
export async function currentNames(rows: EvidenceRow[]): Promise<Map<string, string | null>> {
  const out = new Map<string, string | null>();
  const of = (kind: EvidenceKind) => rows.filter((r) => r.kind === kind);
  const runs = of('model_run');
  if (runs.length) {
    const found = await portalDb.modelRuns.query({filter: runs.map((r) => `run_id = "${r.target}"`).join(' or ')}).catch(() => null);
    if (found) {
      for (const r of runs) {
        const m = found.find((x) => x.run_id === r.target);
        out.set(r.id, m ? m.title ?? m.run_id ?? m.model_name : null);
      }
    }
  }
  const criteria = of('criterion');
  if (criteria.length) {
    const found = await portalDb.gateCriterions.query({filter: `id in (${inList([...new Set(criteria.map((r) => r.target))])})`})
      .catch(() => null);
    if (found)
      for (const r of criteria) out.set(r.id, found.find((x) => x.id === r.target)?.name ?? null);
  }
  for (const r of of('file')) {
    const exists = await grok.dapi.files.exists(r.target).catch(() => undefined);
    if (exists !== undefined)
      out.set(r.id, exists ? r.name : null);
  }
  return out;
}

/** One piece of evidence as the history shows it. */
export interface EvidenceLine {
  row: EvidenceRow;
  state: TargetState;
  current: string | null | undefined;
  text: string;
}

/** The program's evidence by gate-log entry, each target looked up again. */
export async function evidenceOf(programId: string): Promise<Map<string, EvidenceLine[]>> {
  const rows = await portalDb.evidences.query({filter: `program_id = "${programId}"`});
  const names = await currentNames(rows);
  const byDecision = new Map<string, EvidenceLine[]>();
  for (const row of orderEvidence(rows)) {
    const current = names.has(row.id) ? names.get(row.id) : undefined;
    const state = targetState(row.kind, row.name, current);
    const line = {row, state, current, text: describeEvidence(row as any, state, current)};
    const key = row.stage_transition_id ?? '';
    byDecision.set(key, [...byDecision.get(key) ?? [], line]);
  }
  return byDecision;
}

/** Adds [items] to a decision already recorded: by the program's approvers, as for the decision.
 * [transitionId] defaults to the program's latest decision. */
export async function addEvidenceTo(program: ProgramRow, items: EvidenceItem[], transitionId?: string | null):
  Promise<{transition: StageTransitionRow, written: number, errors: string[]}> {
  const ctx = await programContext(program.id);
  if (!ctx.approver)
    throw new Error('Only the portfolio\'s or the program\'s approvers add evidence to a decision');
  const transitions = (await portalDb.stageTransitions.query({filter: `program_id = "${program.id}"`}))
    .filter((t) => t.decision !== 'sync')
    .sort((a, b) => String(b.decided_on ?? b.created_on).localeCompare(String(a.decided_on ?? a.created_on)));
  const t = transitionId ? transitions.find((x) => x.id === transitionId) : transitions[0];
  if (!t)
    throw new Error(transitionId ? `${program.code} has no decision ${transitionId}` : `${program.code} has no decision yet`);
  return {transition: t, ...await writeEvidence(program.id, t.id, items)};
}

/** The program's decisions with their evidence, one row per piece (a decision without evidence has
 * one row with the evidence columns empty), newest decision first. */
export async function decisionEvidenceFrame(program: ProgramRow, transitionId?: string | null): Promise<DG.DataFrame> {
  const [transitions, stages, evidence] = await Promise.all([
    portalDb.stageTransitions.query({filter: `program_id = "${program.id}"`}),
    portalDb.stages.query({filter: `stage_model_id = "${program.stage_model_id}"`}),
    evidenceOf(program.id),
  ]);
  const code = (id?: string) => stages.find((s) => s.id === id)?.code ?? '';
  const rows: Record<string, any>[] = [];
  const shown = transitions.filter((t) => t.decision !== 'sync' && (!transitionId || t.id === transitionId))
    .sort((a, b) => String(b.decided_on ?? b.created_on).localeCompare(String(a.decided_on ?? a.created_on)));
  for (const t of shown) {
    const head = {decision_id: t.id, decided_on: (t.decided_on ?? t.created_on)?.toISOString?.() ?? '', decision: t.decision,
      from: code(t.from_stage_id), to: code(t.to_stage_id), source: t.source ?? '', rationale: t.rationale ?? ''};
    const lines = evidence.get(t.id) ?? [];
    if (!lines.length)
      rows.push({...head, kind: '', name: '', target: '', role: '', state: '', now: '', evidence: ''});
    for (const l of lines) {
      rows.push({...head, kind: l.row.kind, name: l.row.name, target: l.row.target, role: l.row.role ?? '', state: l.state,
        now: l.current ?? '', evidence: l.text});
    }
  }
  const df = rows.length ? DG.DataFrame.fromObjects(rows)! : DG.DataFrame.fromColumns(['decision_id', 'decided_on', 'decision',
    'from', 'to', 'source', 'rationale', 'kind', 'name', 'target', 'role', 'state', 'now', 'evidence']
    .map((n) => DG.Column.fromStrings(n, [])));
  df.name = `${program.code} · Decisions and evidence`;
  return df;
}

export function summarizeAdded(program: ProgramRow, r: {transition: StageTransitionRow, written: number, errors: string[]}): string {
  const when = (r.transition.decided_on ?? r.transition.created_on)?.format?.('YYYY-MM-DD') ?? '';
  return [`${program.code}: ${r.written} piece(s) of evidence added to the ${r.transition.decision} of ${when}`,
    ...r.errors.map((e) => `  not added: ${e}`)].join('\n');
}

/** An item from what a caller types: a run by its run_id (or title), a criterion by its code, a URL
 * or a file path. The name and context are read now, to be frozen. */
export async function evidenceItemFor(program: ProgramRow, kind: EvidenceKind, target: string, name?: string):
  Promise<EvidenceItem> {
  const t = target.trim();
  switch (kind) {
  case 'url':
  case 'file': {
    const item = linkEvidence(t);
    if (item.kind !== kind)
      throw new Error(kind === 'url' ? `"${t}" is not an http(s) address` : `"${t}" is an address; add it as url`);
    return {...item, name: name?.trim() || item.name};
  }
  case 'model_run': {
    const runs = await portalDb.modelRuns.query({filter: `program_id = "${program.id}"`});
    const run = runs.find((r) => r.run_id === t || r.id === t) ?? runs.find((r) => r.title === t);
    if (!run)
      throw new Error(`${program.code} has no saved run ${t}`);
    const selected = (await selectedRunsOf(program.id)).find((s) => s.runId === (run.run_id ?? run.id));
    return selected ? {...runEvidence(selected), name: name?.trim() || selected.title} :
      {kind, target: run.run_id ?? run.id, name: name?.trim() || run.title || run.model_name, role: 'supporting',
        context: {model: run.model_name}};
  }
  case 'criterion': {
    const ctx = await programContext(program.id);
    const c = ctx.criteria.find((x) => x.code === t || x.id === t);
    if (!c)
      throw new Error(`no criterion ${t} applies to ${program.code}`);
    const stage = ctx.stages.find((s) => s.id === c.stageId);
    return {kind, target: c.id, name: name?.trim() || c.name, role: 'supporting',
      context: {code: c.code, stage: stage?.name ?? '', state: ctx.states[c.id] ?? 'not_started'}};
  }
  default:
    return {kind, target: t, name: name?.trim() || t, role: 'supporting'};
  }
}
