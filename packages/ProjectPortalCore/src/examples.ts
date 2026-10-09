/** Seeds the example history (files/examples/history.json) into the EX-* portfolios: gate decisions
 * with rationale and evidence, criteria recorded met, observations and a follow, dated in the past,
 * plus the example documents. Admins only; EX-* only; idempotent. What is missing is worked out in
 * ./model/examples. Decisions go straight to the gate log with their own dates, so no approval
 * request, workspace growth or snapshot is made for them. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {portalDb, PortfolioRow, StageRow} from './generated/db';
import {portfolioByCode, programContext, programData, readinessFor, setCriterion, stagesOf} from './portal';
import {follow, myFollowsAll} from './new-data';
import {loadExampleDocuments} from './documents';
import {EvidenceItem, readinessEvidence, uniqueEvidence, validateEvidence} from './model/evidence';
import {describeSeedPlan, exampleKey, ExampleDecision, ExampleHistory, isExamplePortfolio, planExampleSeed,
  SeedPlan, SeedState} from './model/examples';
import {queryIn, ROW_LIMIT, filterValue} from './rows';

export const EXAMPLE_HISTORY_FILE = 'examples/history.json';

export async function readExampleHistory(pkg: DG.Package): Promise<ExampleHistory> {
  const h = JSON.parse(await pkg.files.readAsText(EXAMPLE_HISTORY_FILE));
  return {criteria: h.criteria ?? [], decisions: h.decisions ?? [], observations: h.observations ?? [], follows: h.follows ?? []};
}

async function isAdmin(): Promise<boolean> {
  const mine = await grok.dapi.groups.currentUserGroups();
  return mine.some((g) => g.name === 'Administrators' || g.friendlyName === 'Administrators');
}

/** A date of the history as a time on that day. */
const at = (day: string) => `${day}T10:00:00.000Z`;

interface Loaded {
  state: SeedState;
  portfolios: Map<string, PortfolioRow>;
  stages: Map<string, StageRow[]>;
}

/** What the platform holds for the programs [h] names. */
async function loadState(h: ExampleHistory): Promise<Loaded> {
  const codes = [...new Set([...h.criteria, ...h.decisions, ...h.observations, ...h.follows].map((x) => x.portfolio))]
    .filter(isExamplePortfolio);
  const portfolios = new Map<string, PortfolioRow>();
  const stages = new Map<string, StageRow[]>();
  const programs: SeedState['programs'] = new Map();
  for (const code of codes) {
    const pf = await portfolioByCode(code).catch(() => null);
    if (!pf)
      continue;
    portfolios.set(code, pf);
    const st = await stagesOf(pf.stage_model_id);
    stages.set(code, st);
    for (const p of await portalDb.programs.query({filter: `portfolio_id = ${filterValue(pf.id)}`, limit: ROW_LIMIT})) {
      programs.set(exampleKey({portfolio: code, program: p.code}),
        {id: p.id, stage: st.find((s) => s.id === p.stage_id)?.code ?? null, status: p.status ?? 'active'});
    }
  }
  const ids = [...programs.values()].map((p) => p.id);
  const [transitions, statuses, observations, follows, criteria] = await Promise.all([
    queryIn(portalDb.stageTransitions, 'program_id', ids),
    queryIn(portalDb.criterionStatuses, 'program_id', ids),
    queryIn(portalDb.observations, 'program_id', ids),
    myFollowsAll(),
    queryIn(portalDb.gateCriterions, 'stage_id', [...stages.values()].flat().map((s) => s.id)),
  ]);
  const stageOf = new Map([...stages.values()].flat().map((s) => [s.id, s.code]));
  const criterionKey = new Map(criteria.map((c) => [c.id, `${stageOf.get(c.stage_id)} ${c.code}`]));
  return {portfolios, stages, state: {
    programs,
    transitions: transitions.map((t) => ({programId: t.program_id, decision: t.decision, rationale: t.rationale ?? null})),
    recorded: new Set(statuses.map((s) => `${s.program_id} ${criterionKey.get(s.gate_criterion_id)}`)),
    observations: observations.map((o) => ({programId: o.program_id, title: o.title})),
    follows: new Set(follows.filter((f) => !f.work_unit_kind).map((f) => f.program_id)),
  }};
}

/** Plans the example history and, unless [dryRun], writes what is missing. */
export async function seedExamples(pkg: DG.Package, dryRun = true): Promise<{plan: SeedPlan, report: string}> {
  if (!await isAdmin())
    throw new Error('Only an administrator seeds the example history');
  const h = await readExampleHistory(pkg);
  const loaded = await loadState(h);
  const plan = planExampleSeed(h, loaded.state);
  if (dryRun)
    return {plan, report: describeSeedPlan(plan, true)};
  const errors = await writePlan(plan, loaded);
  const documents = await loadExampleDocuments(pkg);
  const linked = documents.split('\n').filter((l) => / linked/.test(l) && !/already linked|not linked/.test(l)).length;
  const report = [describeSeedPlan(plan, false), `Example documents: ${linked} linked`,
    ...errors.map((e) => `  error: ${e}`)].join('\n');
  return {plan, report};
}

async function writePlan(plan: SeedPlan, loaded: Loaded): Promise<string[]> {
  const errors: string[] = [];
  const me = DG.User.current().id;
  const program = (x: {portfolio: string, program: string}) => loaded.state.programs.get(exampleKey(x))!;
  const stageId = (portfolio: string, code?: string | null) =>
    code ? loaded.stages.get(portfolio)?.find((s) => s.code === code)?.id ?? null : null;

  // criteria first: the decisions' evidence reads them
  for (const c of plan.criteria) {
    const p = program(c);
    const criterion = (await portalDb.gateCriterions.query({filter: `stage_id = ${filterValue(stageId(c.portfolio, c.stage) ?? '')} ` +
      `and code = ${filterValue(c.criterion)}`, limit: ROW_LIMIT}))[0];
    if (!criterion) {
      errors.push(`${exampleKey(c)}: no criterion ${c.stage} ${c.criterion} (load the portfolios first)`);
      continue;
    }
    await setCriterion(p.id, criterion.id, c.status, `${c.evidence} (${c.on})`);
  }

  // decisions in date order; each program walks from where it stands
  const walk = new Map<string, string | null>();
  for (const d of plan.decisions) {
    const p = program(d);
    const from = d.decision === 'record' ? stageId(d.portfolio, d.from) :
      walk.has(exampleKey(d)) ? walk.get(exampleKey(d))! : stageId(d.portfolio, p.stage);
    const to = d.decision === 'record' || d.decision === 'advance' ? stageId(d.portfolio, d.to) : from;
    if (d.decision !== 'record')
      walk.set(exampleKey(d), to);
    try {
      const [t] = await portalDb.stageTransitions.insert({program_id: p.id, from_stage_id: from ?? undefined,
        to_stage_id: to ?? undefined, decision: d.decision, ...(d.outcome ? {outcome: d.outcome} : {}),
        decided_by: me, decided_on: at(d.on) as any, rationale: d.rationale, source: 'approver'} as any);
      const items = uniqueEvidence([...await gateEvidence(d, p.id, to, loaded),
        ...(d.links ?? []).map((l): EvidenceItem => ({kind: 'url', target: l.url, name: l.name, role: 'supporting'}))])
        .filter((e) => validateEvidence(e).length === 0);
      if (items.length) {
        await portalDb.evidences.insert(items.map((e) => ({program_id: p.id, stage_transition_id: t.id, kind: e.kind,
          target: e.target, name: e.name, role: e.role, context: e.context ?? {}, added_by: me, added_on: at(d.on) as any})));
      }
    } catch (e: any) {
      errors.push(`${exampleKey(d)} ${d.on} ${d.decision}: ${e?.message ?? e}`);
    }
  }
  for (const m of plan.moves) {
    const [portfolio] = m.key.split(' ');
    const p = loaded.state.programs.get(m.key)!;
    await portalDb.programs.update(p.id, {stage_id: stageId(portfolio, m.stage) ?? null,
      ...(m.status ? {status: m.status} : {})} as any).catch((e: any) => errors.push(`${m.key}: not moved: ${e?.message ?? e}`));
  }

  for (const o of plan.observations) {
    const p = program(o);
    await portalDb.observations.insert({program_id: p.id, subject_kind: o.subject,
      subject_code: o.subject === 'program' ? undefined : o.code, title: o.title, finding: o.finding, evidence: {items: []},
      status: 'current', observed_by: me, observed_on: at(o.on) as any} as any)
      .catch((e: any) => errors.push(`${exampleKey(o)}: observation not written: ${e?.message ?? e}`));
  }
  for (const f of plan.follows) {
    const p = program(f);
    await follow(loaded.portfolios.get(f.portfolio)!, await portalDb.programs.get(p.id))
      .catch((e: any) => errors.push(`${exampleKey(f)}: not followed: ${e?.message ?? e}`));
  }
  return errors;
}

/** The readiness of the stage an advance or a record leads into, as the gate panel would offer it. */
async function gateEvidence(d: ExampleDecision, programId: string, to: string | null, loaded: Loaded): Promise<EvidenceItem[]> {
  if (!to || (d.decision !== 'advance' && d.decision !== 'record'))
    return [];
  const ctx = await programContext(programId);
  const pf = loaded.portfolios.get(d.portfolio)!;
  const rows = ctx.criteria.some((c) => c.stageId === to && c.condition) ?
    await programData(pf, ctx.program.code).catch(() => null) : null;
  const stage = loaded.stages.get(d.portfolio)!.find((s) => s.id === to)!;
  return readinessEvidence(readinessFor(ctx, to, rows), stage.name);
}
