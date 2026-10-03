/** Profile tooling against the platform (F-18): a loaded portfolio read back into its file format,
 * diffs between loaded portfolios and files, and bundles that move a set of portfolios between
 * instances. The rules live in ./model/profile. */
import * as grok from 'datagrok-api/grok';
import {portalDb} from './generated/db';
import {loadPortfolioSet, portfolioByCode, stagesOf, COMPOUND_FEED} from './portal';
import {CriterionConfig, FeedConfig, PortfolioConfig, StageConfig} from './model/config';
import {CardParam, ModelCardConfig, ModelCatalogueConfig, ModelConfig, readParams, splitCodes} from './model/catalogue';
import {WorkUnitFeedConfig, WorkUnitKind} from './model/work-units';
import {funcOf, funcParams} from './catalogue';
import {cardParams} from './model/catalogue';
import {
  diffPortfolios, ImportPlan, makeBundle, planImport, PortfolioDiff, ProfileBundle, readBundle, summarizeImport,
} from './model/profile';
import {filterValue, ROW_LIMIT} from './rows';

const split = (s: string | null | undefined) => (s ?? '').split('|').filter((x) => x.length > 0);
const orNull = <T>(xs: T[]) => xs.length ? xs : undefined;

/** Drops the keys whose value is undefined, so the JSON reads like a hand-written file. */
function compact<T extends object>(o: T): T {
  return Object.fromEntries(Object.entries(o).filter(([, v]) => v !== undefined && v !== null)) as T;
}

async function groupName(id: string | null | undefined): Promise<string | undefined> {
  if (!id)
    return undefined;
  const g = await grok.dapi.groups.find(id).catch(() => null);
  return g ? g.friendlyName || g.name : undefined;
}

/** The catalogue entries a portfolio declares: unit operations and run types that name it (or,
 * loaded before 0.2.13 recorded that, that its models use), and the models that name it, each with
 * its current card. */
async function catalogueOf(code: string): Promise<ModelCatalogueConfig | undefined> {
  const [ops, types, models] = await Promise.all([portalDb.unitOperations.query({limit: ROW_LIMIT}), portalDb.runTypes.query({limit: ROW_LIMIT}),
    portalDb.models.query({limit: ROW_LIMIT})]);
  const mine = models.filter((m) => split(m.portfolios).includes(code)).sort((a, b) => a.code.localeCompare(b.code));
  const typeCode = new Map(types.map((t) => [t.id, t.code]));
  const usedOps = new Set(mine.flatMap((m) => split(m.unit_operations)));
  const usedTypes = new Set(mine.map((m) => m.default_run_type_id ? typeCode.get(m.default_run_type_id) : undefined));
  const claims = (row: any, used: Set<string | undefined>) => row.portfolios ? split(row.portfolios).includes(code) : used.has(row.code);
  const byOrdinal = (a: any, b: any) => (a.ordinal ?? 1e9) - (b.ordinal ?? 1e9) || a.code.localeCompare(b.code);
  const out: ModelCatalogueConfig = {};
  const unitOperations = ops.filter((o) => claims(o, usedOps)).sort(byOrdinal)
    .map((o) => compact({code: o.code, name: o.name, ordinal: o.ordinal, description: o.description}));
  const runTypes = types.filter((t) => claims(t, usedTypes)).sort(byOrdinal)
    .map((t) => compact({code: t.code, name: t.name, definition: t.definition, ordinal: t.ordinal,
      isExecutionMode: t.is_execution_mode}));
  const configs: ModelConfig[] = [];
  for (const m of mine) {
    const cards = await portalDb.modelCards.query({filter: `model_id = ${filterValue(m.id)}`, limit: ROW_LIMIT});
    const card = cards.find((c) => c.model_version === m.current_version);
    configs.push(compact({code: m.code, name: m.name, package: m.package, category: m.category, description: m.description,
      status: m.status as any, unitOperations: orNull(split(m.unit_operations)),
      defaultRunType: m.default_run_type_id ? typeCode.get(m.default_run_type_id) : undefined,
      card: card ? compact({version: card.model_version, question: card.question, intendedUse: card.intended_use,
        inputs: readParams(card.inputs) ?? undefined, outputs: readParams(card.outputs) ?? undefined,
        applicability: card.applicability, performance: card.performance, lineage: card.lineage,
        supersedes: card.supersedes, limitations: card.limitations, owner: card.owner, approval: card.approval as any,
        approvedOn: card.approved_on ? new Date(card.approved_on.valueOf()).toISOString().slice(0, 10) : undefined,
        appliesTo: (card.applies_to as ModelCardConfig['appliesTo']) ?? undefined}) : undefined}));
  }
  if (unitOperations.length) out.unitOperations = unitOperations;
  if (runTypes.length) out.runTypes = runTypes;
  if (configs.length) out.models = configs;
  return Object.keys(out).length ? out : undefined;
}

/** A loaded portfolio in its file format: what loading the file wrote, read back. Rows a later
 * file no longer declares are still there (loading never deletes), so they show here too. */
export async function exportPortfolioConfig(code: string): Promise<PortfolioConfig> {
  const p = await portfolioByCode(code);
  const [model, stages, feeds] = await Promise.all([
    portalDb.stageModels.get(p.stage_model_id),
    stagesOf(p.stage_model_id),
    portalDb.registryFeeds.query({filter: `portfolio_id = ${filterValue(p.id)}`, limit: ROW_LIMIT}),
  ]);
  const ids = stages.map((s) => `"${s.id}"`).join(', ');
  const criteria = ids ? await portalDb.gateCriterions.query({filter: `stage_id in (${ids})`, limit: ROW_LIMIT}) : [];
  const stageConfigs: StageConfig[] = stages.map((s) => compact({
    code: s.code, name: s.name, ordinal: s.ordinal, phase: s.phase, isGate: s.is_gate, isTerminal: s.is_terminal,
    governance: s.governance, aliases: orNull(split(s.spellings)), description: s.description,
    criteria: orNull(criteria.filter((c) => c.stage_id === s.id)
      .sort((a, b) => a.created_on.valueOf() - b.created_on.valueOf() || a.code.localeCompare(b.code))
      .map((c): CriterionConfig => compact({code: c.code, name: c.name, kind: c.kind, rule: c.rule,
        condition: c.condition as any, required: c.required, modality: c.modality, artifact: c.artifact,
        runFilter: c.run_filter as any}))),
    folders: orNull(split(s.folders)),
  }));
  const byCode = (a: {code: string}, b: {code: string}) => a.code.localeCompare(b.code);
  const programFeeds: FeedConfig[] = feeds.filter((f) => (f.target ?? 'program') === 'program').sort(byCode)
    .map((f) => compact({code: f.code, system: f.system, connection: f.connection, sql: f.sql, mapping: f.mapping as any}));
  const compoundFeed = feeds.find((f) => f.target === 'compound' && f.code === COMPOUND_FEED);
  const unitFeeds: WorkUnitFeedConfig[] = feeds.filter((f) => f.target && f.target !== 'program' && f.target !== 'compound')
    .sort(byCode).map((f) => compact({code: f.code, kind: f.target as WorkUnitKind,
      system: f.system === `${p.code}-${f.target}` ? undefined : f.system, connection: f.connection, sql: f.sql}));
  const groups = compact({viewers: await groupName(p.viewers), contributors: await groupName(p.contributors),
    approvers: await groupName(p.approvers)});
  return compact({
    code: p.code, name: p.name, setting: p.setting, modalityClass: p.modality_class, stageAuthority: p.stage_authority,
    programNoun: p.program_noun, workspaceRoot: p.workspace_root, workstreams: orNull(split(p.workstreams)),
    flowFolders: p.flow_folders as any,
    description: p.description, schedule: p.schedule, owners: orNull(split(p.owners)),
    autoAdvance: p.auto_advance as any,
    personas: p.personas as any,
    documentKinds: (p.document_kinds as any)?.kinds,
    approval: p.approval as any,
    groups: Object.keys(groups).length ? groups : undefined,
    handoff: p.handoff as any,
    workUnits: orNull(split(p.work_units)) as WorkUnitKind[] | undefined,
    seriesMaturity: orNull(split(p.series_maturity)),
    workUnitFeeds: orNull(unitFeeds),
    compoundFeed: compoundFeed ? {system: compoundFeed.system, connection: compoundFeed.connection, sql: compoundFeed.sql} : undefined,
    programData: p.program_data as any,
    programDatasets: (p.program_datasets as any)?.datasets,
    attention: (p.attention as any)?.rules,
    catalogue: await catalogueOf(p.code),
    stageModel: compact({code: model.code, name: model.name, version: model.model_version, setting: model.setting,
      description: model.description, stages: stageConfigs}),
    feeds: programFeeds,
  }) as PortfolioConfig;
}

/** A file's model cards completed from the functions, as loading stores them: inputs and outputs
 * the card leaves out come from the signature. Compare a file with a loaded portfolio after this. */
export function completeCards(c: PortfolioConfig): PortfolioConfig {
  if (!c.catalogue?.models?.length)
    return c;
  const out: PortfolioConfig = JSON.parse(JSON.stringify(c));
  for (const m of out.catalogue!.models!) {
    if (!m.card)
      continue;
    const func = funcOf(m.code);
    if (!func)
      continue;
    const params = funcParams(func);
    m.card.inputs = cardParams(m.card.inputs, params.inputs) as CardParam[];
    m.card.outputs = cardParams(m.card.outputs, params.outputs) as CardParam[];
  }
  return out;
}

/** A portfolio named by its code (the loaded one) or given as JSON (a file). */
async function portfolioFrom(ref: string): Promise<{label: string, config: PortfolioConfig}> {
  const text = ref.trim();
  if (text.startsWith('{')) {
    const config = completeCards(JSON.parse(text));
    return {label: `${config.code} (file)`, config};
  }
  return {label: `${text} (loaded)`, config: await exportPortfolioConfig(text)};
}

export async function diffPortfolioRefs(a: string, b: string): Promise<{title: string, diff: PortfolioDiff}> {
  const [x, y] = await Promise.all([portfolioFrom(a), portfolioFrom(b)]);
  return {title: `${x.label} → ${y.label}`, diff: diffPortfolios(x.config, y.config)};
}

/** A bundle of the loaded portfolios [codes], or of every loaded portfolio. */
export async function exportProfileBundle(codes?: string[]): Promise<ProfileBundle> {
  const all = (await portalDb.portfolios.query({limit: ROW_LIMIT})).map((p) => p.code).sort();
  const wanted = codes?.length ? codes : all;
  const missing = wanted.filter((c) => !all.includes(c));
  if (missing.length)
    throw new Error(`not loaded: ${missing.join(', ')}`);
  const portfolios: PortfolioConfig[] = [];
  for (const c of wanted)
    portfolios.push(await exportPortfolioConfig(c));
  return makeBundle(portfolios, new Date().toISOString());
}

/** Checks a bundle as a set, diffs each portfolio against the loaded one, and unless [dryRun]
 * loads the set. A bundle that does not validate is refused whole. */
export async function importProfileBundle(input: string | unknown, dryRun: boolean): Promise<{plan: ImportPlan, report: string}> {
  const {bundle, errors} = readBundle(input);
  if (!bundle) {
    const plan: ImportPlan = {errors, portfolios: []};
    return {plan, report: summarizeImport(plan, dryRun)};
  }
  bundle.portfolios = bundle.portfolios.map(completeCards);
  const loadedCodes = (await portalDb.portfolios.query({limit: ROW_LIMIT})).map((p) => p.code);
  const loaded = new Map<string, PortfolioConfig>();
  for (const p of bundle.portfolios)
    if (loadedCodes.includes(p.code)) loaded.set(p.code, await exportPortfolioConfig(p.code));
  const plan = planImport(bundle, loaded, loadedCodes);
  let report = summarizeImport(plan, dryRun);
  if (!plan.errors.length && !dryRun)
    report += `\n${await loadPortfolioSet(bundle.portfolios, loadedCodes)}`;
  return {plan, report};
}
