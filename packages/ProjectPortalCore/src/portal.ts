/** The portal against the platform: loading portfolio files into the `portal` schema, running
 * registry feeds, applying sync plans, recording gate decisions, and the routing audit.
 * The rules live in ./model; this file only reads, writes and shares. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {
  FeedConfig, PortfolioConfig, SourceConfig, StageModelConfig, validatePortfolio, validatePortfolioSet,
} from './model/config';
import {FeedRow, planSync, ProgramRef, StageRef, summarizePlan, SyncPlan} from './model/sync-plan';
import {checkDecision, Decision, GateVerdict, CriterionState, statusAfter} from './model/gates';
import {portalDb, PortfolioRow, ProgramRow, StageRow} from './generated/db';

type Groups = {viewers?: string, contributors?: string, approvers?: string};

// ---------------------------------------------------------------------------------------------
// Groups

const groupCache = new Map<string, DG.Group>();

async function findGroup(name: string): Promise<DG.Group | null> {
  const cached = groupCache.get(name);
  if (cached)
    return cached;
  const g = await grok.dapi.groups.filter(`friendlyName = "${name}"`).first().catch(() => null) ??
    await grok.dapi.groups.filter(`name = "${name}"`).first().catch(() => null);
  if (g)
    groupCache.set(name, g);
  return g ?? null;
}

/** The portfolio's three groups, created when missing and nested approvers ⊂ contributors ⊂
 * viewers, so one View share to viewers and one Edit share to contributors covers everyone. */
async function ensureGroups(names: Groups | undefined, report: string[]): Promise<Record<keyof Groups, string | null>> {
  const ids: Record<keyof Groups, string | null> = {viewers: null, contributors: null, approvers: null};
  const made: Partial<Record<keyof Groups, DG.Group>> = {};
  for (const role of ['viewers', 'contributors', 'approvers'] as (keyof Groups)[]) {
    const name = names?.[role];
    if (!name)
      continue;
    let g = await findGroup(name);
    if (!g) {
      g = await grok.dapi.groups.createNew(name);
      groupCache.set(name, g);
      report.push(`created group ${name}`);
    }
    made[role] = g;
    ids[role] = g.id;
  }
  const nest = async (child?: DG.Group, parent?: DG.Group) => {
    if (!child || !parent)
      return;
    const fresh = await grok.dapi.groups.find(child.id);
    if (!fresh.memberships.some((m) => m.id === parent.id))
      await grok.dapi.groups.includeTo(fresh, parent);
  };
  await nest(made.approvers, made.contributors);
  await nest(made.contributors, made.viewers);
  return ids;
}

/** Whether the current user may record gate decisions in [portfolio]. */
export async function isApprover(portfolio: PortfolioRow): Promise<boolean> {
  const mine = await grok.dapi.groups.currentUserGroups();
  const ids = new Set(mine.map((g) => g.id));
  if (mine.some((g) => g.name === 'Administrators' || g.friendlyName === 'Administrators'))
    return true;
  return portfolio.approvers !== undefined && portfolio.approvers !== null && ids.has(portfolio.approvers);
}

/** Shares a program row with its portfolio: View to viewers, Edit to contributors.
 * Today the server refuses `permissions.grant` on a promoted domain row even for its author
 * (privileges_service.savePermission checks Share against the projects repository, which never
 * holds domain rows), so callers treat a failure here as a reported gap, not a fatal error:
 * the program stays private to its author, which fails closed. */
export async function shareProgram(programId: string, portfolio: PortfolioRow): Promise<void> {
  const targets: [string | undefined, boolean][] = [[portfolio.viewers, false], [portfolio.contributors, true]];
  if (!targets.some(([g]) => g))
    return;
  await portalDb.programs.promote(programId);
  const [entity] = await grok.dapi.getEntities([programId]);
  if (!entity)
    throw new Error(`program ${programId} is not an entity after promotion; cannot share it`);
  for (const [groupId, edit] of targets) {
    if (!groupId)
      continue;
    const group = await grok.dapi.groups.find(groupId);
    await grok.dapi.permissions.grant(entity, group, edit);
  }
}

// ---------------------------------------------------------------------------------------------
// Loading configuration

async function upsertByKey<T extends {id: string}>(
  client: {getByKey(k: any): Promise<T | null>, insert(v: any): Promise<{id: string}[]>,
    update(id: string, v: any): Promise<any>},
  key: Record<string, unknown>, values: Record<string, unknown>): Promise<string> {
  const existing = await client.getByKey(key);
  if (existing === null) {
    const [r] = await client.insert({...key, ...values});
    return (r as any).existingId ?? r.id;
  }
  const changed: Record<string, unknown> = {};
  for (const [k, v] of Object.entries(values)) {
    const was = (existing as any)[k];
    if (JSON.stringify(was ?? null) !== JSON.stringify(v ?? null))
      changed[k] = v;
  }
  if (Object.keys(changed).length > 0)
    await client.update(existing.id, changed);
  return existing.id;
}

export async function loadStageModel(m: StageModelConfig, report: string[]): Promise<string> {
  const modelId = await upsertByKey(portalDb.stageModels, {code: m.code},
    {name: m.name, model_version: m.version ?? null, setting: m.setting, description: m.description ?? null});
  for (const s of m.stages) {
    const stageId = await upsertByKey(portalDb.stages, {stage_model_id: modelId, code: s.code}, {
      name: s.name, ordinal: s.ordinal, phase: s.phase ?? null, is_gate: s.isGate ?? true,
      is_terminal: s.isTerminal ?? false, governance: s.governance ?? null,
      spellings: (s.aliases ?? []).join('|') || null,
      description: s.description ?? null,
    });
    for (const c of s.criteria ?? []) {
      await upsertByKey(portalDb.gateCriterions, {stage_id: stageId, code: c.code},
        {name: c.name, kind: c.kind, rule: c.rule ?? null, required: c.required ?? true});
    }
  }
  report.push(`stage model ${m.code}: ${m.stages.length} stages`);
  return modelId;
}

/** Loads one portfolio file: groups, stage model, portfolio row, feeds. Idempotent. */
export async function loadPortfolio(c: PortfolioConfig): Promise<string> {
  const errors = validatePortfolio(c);
  if (errors.length > 0)
    throw new Error(`${c.code}: ${errors.join('; ')}`);
  const report: string[] = [];
  const groups = await ensureGroups(c.groups, report);
  const modelId = await loadStageModel(c.stageModel, report);
  const portfolioId = await upsertByKey(portalDb.portfolios, {code: c.code}, {
    name: c.name, setting: c.setting, modality_class: c.modalityClass ?? null, stage_model_id: modelId,
    stage_authority: c.stageAuthority, program_noun: c.programNoun ?? 'Program',
    workspace_root: c.workspaceRoot ?? null, description: c.description ?? null,
    viewers: groups.viewers, contributors: groups.contributors, approvers: groups.approvers,
  });
  for (const f of c.feeds) {
    await upsertByKey(portalDb.registryFeeds, {portfolio_id: portfolioId, code: f.code},
      {system: f.system, connection: f.connection, sql: f.sql, mapping: f.mapping ?? null, enabled: true});
  }
  report.push(`portfolio ${c.code}: ${c.feeds.length} feed(s)`);
  return report.join('\n');
}

export async function readPortfolioFiles(pkg: DG.Package): Promise<PortfolioConfig[]> {
  // the pattern is a substring match, not a glob
  const files = (await pkg.files.list('portfolios', false, '.json')).filter((f) => f.name.endsWith('.json'));
  const configs: PortfolioConfig[] = [];
  for (const f of files.sort((a, b) => a.name.localeCompare(b.name)))
    configs.push(JSON.parse(await pkg.files.readAsText(`portfolios/${f.name}`)));
  return configs;
}

export async function readSources(pkg: DG.Package): Promise<SourceConfig[]> {
  return JSON.parse(await pkg.files.readAsText('sources.json'));
}

/** Loads every portfolio file the package ships, after checking them as a set. */
export async function loadAllPortfolios(pkg: DG.Package): Promise<string> {
  const configs = await readPortfolioFiles(pkg);
  if (configs.length === 0)
    throw new Error('no portfolio files found under files/portfolios');
  return loadPortfolioSet(configs);
}

/** Loads a set of portfolios after checking them together. A customer package hands its own
 * portfolio files to the core this way (see the loadPortfolioConfigs function). */
export async function loadPortfolioSet(configs: PortfolioConfig[]): Promise<string> {
  const errors = validatePortfolioSet(configs);
  if (errors.length > 0)
    throw new Error(errors.join('\n'));
  const out: string[] = [];
  for (const c of configs)
    out.push(await loadPortfolio(c));
  return out.join('\n');
}

// ---------------------------------------------------------------------------------------------
// Feeds and sync

/** A connection by name, friendly name, or `Namespace:Name`. The dapi `name =` filter matches the
 * friendly name, which a connection may spell differently (MyData vs my_data). */
async function findConnection(ref: string): Promise<DG.DataConnection | null> {
  if (ref.includes(':'))
    return await grok.functions.eval(ref).catch(() => null) as DG.DataConnection | null;
  const byFriendly = await grok.dapi.connections.filter(`name = "${ref}"`).list().catch(() => []);
  const exact = (cs: DG.DataConnection[]) => cs.find((c) => c.name === ref) ?? cs.find((c) => c.friendlyName === ref);
  const hit = exact(byFriendly);
  if (hit)
    return hit;
  const all = await grok.dapi.connections.filter(ref).list({pageSize: 100}).catch(() => []);
  return exact(all) ?? all.find((c) => c.name.toLowerCase() === ref.toLowerCase()) ?? null;
}

async function runSql(connectionName: string, sql: string): Promise<FeedRow[]> {
  const conn = await findConnection(connectionName);
  if (!conn)
    throw new Error(`data connection "${connectionName}" is not found or not shared with you`);
  const df: DG.DataFrame = await conn.query('portal-feed', sql).apply();
  const rows: FeedRow[] = [];
  const cols = df.columns.toList();
  for (let i = 0; i < df.rowCount; i++) {
    const row: FeedRow = {};
    for (const c of cols)
      row[c.name.toLowerCase()] = c.isNone(i) ? null : c.get(i);
    rows.push(row);
  }
  return rows;
}

export async function runFeed(feed: Pick<FeedConfig, 'connection' | 'sql'>): Promise<FeedRow[]> {
  return runSql(feed.connection, feed.sql);
}

export async function portfolioByCode(code: string): Promise<PortfolioRow> {
  const p = await portalDb.portfolios.getByKey({code});
  if (!p)
    throw new Error(`portfolio ${code} is not loaded`);
  return p;
}

export async function stagesOf(stageModelId: string): Promise<StageRow[]> {
  const stages = await portalDb.stages.query({filter: `stage_model_id = "${stageModelId}"`});
  return stages.sort((a, b) => a.ordinal - b.ordinal);
}

const toRef = (p: ProgramRow, portfolioCode: string): ProgramRef => ({
  id: p.id, code: p.code, portfolioCode, stageId: p.stage_id ?? null, status: (p.status ?? 'active') as any,
  source: p.source ?? null, externalRef: p.external_ref ?? null, predecessorId: p.predecessor_id ?? null,
  name: p.name, modality: p.modality ?? null, target: p.target ?? null,
  therapeuticArea: p.therapeutic_area ?? null, sourceStage: p.source_stage ?? null,
});

/** Plans (and unless [dryRun], applies) every enabled feed of one portfolio. */
export async function syncPortfolio(code: string, dryRun = false): Promise<SyncPlan[]> {
  const portfolio = await portfolioByCode(code);
  const [stages, feeds, portfolios] = await Promise.all([
    stagesOf(portfolio.stage_model_id),
    portalDb.registryFeeds.query({filter: `portfolio_id = "${portfolio.id}"`}),
    portalDb.portfolios.query({}),
  ]);
  const codeOf = new Map(portfolios.map((p) => [p.id, p.code]));
  const stageRefs: StageRef[] = stages.map((s) => ({id: s.id, code: s.code, name: s.name, ordinal: s.ordinal,
    aliases: (s.spellings ?? '').split('|').filter((a) => a.length > 0)}));
  const plans: SyncPlan[] = [];
  for (const feed of feeds.filter((f) => f.enabled !== false)) {
    const rows = await runFeed({connection: feed.connection, sql: feed.sql});
    const all = await portalDb.programs.query({filter: `external_ref like "${feed.system}:%"`});
    const mine = (await portalDb.programs.query({filter: `portfolio_id = "${portfolio.id}"`}))
      .map((p) => toRef(p, code));
    const plan = planSync({
      portfolioCode: code, system: feed.system, stageAuthority: portfolio.stage_authority as any,
      stages: stageRefs, statusMap: (feed.mapping as any)?.status, rows,
      existing: mine.filter((p) => p.source !== 'feed' || (p.externalRef ?? '').startsWith(`${feed.system}:`)),
      elsewhere: all.filter((p) => p.portfolio_id !== portfolio.id).map((p) => toRef(p, codeOf.get(p.portfolio_id) ?? '?')),
    });
    plans.push(plan);
    if (!dryRun) {
      await applyPlan(plan, portfolio);
      await portalDb.registryFeeds.update(feed.id, {last_synced_on: new Date().toISOString() as any});
    }
  }
  return plans;
}

async function applyPlan(plan: SyncPlan, portfolio: PortfolioRow): Promise<void> {
  const now = new Date().toISOString();
  const me = DG.User.current().id;
  const idOf = new Map<string, string>();
  for (const v of plan.inserts) {
    const [r] = await portalDb.programs.insert({...v, portfolio_id: portfolio.id,
      stage_model_id: portfolio.stage_model_id, source: 'feed', last_synced_on: now} as any);
    const id = (r as any).existingId ?? r.id;
    idOf.set(v.code, id);
    try {
      await shareProgram(id, portfolio);
    } catch (e: any) {
      plan.warnings.push(`${v.code}: not shared with ${portfolio.code}'s groups (${e?.message ?? e})`);
    }
  }
  for (const u of plan.updates) {
    idOf.set(u.code, u.id);
    await portalDb.programs.update(u.id, {...u.changes, last_synced_on: now} as any);
  }
  for (const s of plan.stageChanges) {
    const programId = s.programId ?? idOf.get(s.code);
    if (programId) {
      await portalDb.stageTransitions.insert({program_id: programId, from_stage_id: s.fromStageId ?? undefined,
        to_stage_id: s.toStageId ?? undefined, decision: 'sync', decided_by: me, decided_on: now as any,
        rationale: `Registry ${plan.system}`});
    }
  }
  for (const m of plan.missing)
    await portalDb.programs.update(m.id, {status: m.status});
  // a predecessor whose successor just arrived has moved on
  for (const l of plan.lineage) {
    if (!plan.overlaps.some((o) => o.code === l.code))
      await portalDb.programs.update(l.predecessorId, {status: 'transitioned'}).catch(() => undefined);
  }
}

// ---------------------------------------------------------------------------------------------
// Gate decisions

export interface ProgramContext {
  program: ProgramRow;
  portfolio: PortfolioRow;
  stages: StageRow[];
  criteria: {id: string, stageId: string, code: string, name: string, kind: string, required: boolean}[];
  states: Record<string, CriterionState>;
  approver: boolean;
}

export async function programContext(programId: string): Promise<ProgramContext> {
  const program = await portalDb.programs.get(programId);
  const [portfolio, stages, states] = await Promise.all([
    portalDb.portfolios.get(program.portfolio_id),
    stagesOf(program.stage_model_id),
    portalDb.criterionStatuses.query({filter: `program_id = "${programId}"`}),
  ]);
  const stageIds = stages.map((s) => `"${s.id}"`).join(', ');
  const criteria = stageIds === '' ? [] :
    await portalDb.gateCriterions.query({filter: `stage_id in (${stageIds})`});
  return {
    program, portfolio, stages,
    criteria: criteria.map((c) => ({id: c.id, stageId: c.stage_id, code: c.code, name: c.name, kind: c.kind,
      required: c.required !== false})),
    states: Object.fromEntries(states.map((s) => [s.gate_criterion_id, s.status as CriterionState])),
    approver: await isApprover(portfolio),
  };
}

export function verdictFor(ctx: ProgramContext, decision: Decision, targetStageId?: string | null,
  rationale?: string): GateVerdict {
  return checkDecision({
    decision, stageAuthority: ctx.portfolio.stage_authority as any, programStatus: (ctx.program.status ?? 'active') as any,
    stages: ctx.stages.map((s) => ({id: s.id, code: s.code, name: s.name, ordinal: s.ordinal, isTerminal: s.is_terminal})),
    currentStageId: ctx.program.stage_id ?? null, targetStageId, criteria: ctx.criteria,
    criterionStates: ctx.states, rationale, isApprover: ctx.approver,
  });
}

/** Records a gate decision after checking it again server-side of the UI. */
export async function recordDecision(programId: string, decision: Decision, targetStageId: string | null,
  rationale: string): Promise<GateVerdict> {
  const ctx = await programContext(programId);
  const verdict = verdictFor(ctx, decision, targetStageId, rationale);
  if (!verdict.allowed)
    throw new Error(verdict.reasons.join('; '));
  const now = new Date().toISOString();
  await portalDb.stageTransitions.insert({program_id: programId, from_stage_id: ctx.program.stage_id,
    to_stage_id: verdict.toStageId ?? ctx.program.stage_id, decision, decided_by: DG.User.current().id,
    decided_on: now as any, rationale});
  const status = statusAfter(decision, (ctx.program.status ?? 'active') as any);
  await portalDb.programs.update(programId, {
    ...(verdict.toStageId ? {stage_id: verdict.toStageId} : {}),
    ...(status !== ctx.program.status ? {status} : {}),
  });
  return verdict;
}

export async function setCriterion(programId: string, criterionId: string, status: CriterionState,
  evidence?: string): Promise<void> {
  await portalDb.criterionStatuses.upsert({program_id: programId, gate_criterion_id: criterionId, status,
    evidence: evidence ?? undefined});
}

// ---------------------------------------------------------------------------------------------
// Routing audit

export interface RoutingReport {
  system: string;
  total: number;
  claims: Record<string, string[]>;
  unclaimed: string[];
  claimedTwice: {code: string, portfolios: string[]}[];
}

/** For each source system: which keys no portfolio's feed takes, and which two feeds both take. */
export async function routingAudit(sources: SourceConfig[]): Promise<RoutingReport[]> {
  const [portfolios, feeds] = await Promise.all([portalDb.portfolios.query({}), portalDb.registryFeeds.query({})]);
  const codeOf = new Map(portfolios.map((p) => [p.id, p.code]));
  const reports: RoutingReport[] = [];
  for (const source of sources) {
    const universe = (await runSql(source.connection, source.sql)).map((r) => String(r.code));
    const claimedBy = new Map<string, string[]>();
    const claims: Record<string, string[]> = {};
    for (const feed of feeds.filter((f) => f.system === source.system && f.enabled !== false)) {
      const portfolio = codeOf.get(feed.portfolio_id) ?? '?';
      const codes = (await runFeed(feed)).map((r) => String(r.code));
      claims[portfolio] = codes;
      for (const c of codes)
        claimedBy.set(c, [...(claimedBy.get(c) ?? []), portfolio]);
    }
    reports.push({
      system: source.system, total: universe.length, claims,
      unclaimed: universe.filter((c) => !claimedBy.has(c)).sort(),
      claimedTwice: [...claimedBy.entries()].filter(([, ps]) => ps.length > 1)
        .map(([code, ps]) => ({code, portfolios: ps})),
    });
  }
  return reports;
}

export function summarizeRouting(reports: RoutingReport[]): string {
  return reports.map((r) => {
    const claimed = Object.entries(r.claims).map(([p, c]) => `${p} ${c.length}`).join(', ');
    return `${r.system}: ${r.total} keys; ${claimed}; unclaimed ${r.unclaimed.length}` +
      (r.unclaimed.length ? ` (${r.unclaimed.join(', ')})` : '') +
      `; claimed twice ${r.claimedTwice.length}` +
      (r.claimedTwice.length ? ` (${r.claimedTwice.map((c) => `${c.code}: ${c.portfolios.join('+')}`).join(', ')})` : '');
  }).join('\n');
}

export {summarizePlan};
