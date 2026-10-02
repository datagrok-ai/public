/** The model catalogue against the platform: loading a merged catalogue into `unit_operation`,
 * `run_type`, `model` and `model_card`, and reading a model's card back with what the function
 * itself says. The rules live in ./model/catalogue. */
import * as DG from 'datagrok-api/dg';
import {ModelCardRow, ModelRow, portalDb, RunTypeRow, UnitOperationRow} from './generated/db';
import {
  CardParam, cardMismatches, cardParams, FuncParamInfo, MergedCatalogue, mergeCatalogues, ModelCatalogueConfig,
  readParams, splitCodes, storedParams, unionCodes,
} from './model/catalogue';
// portal.ts imports from here too; the cycle is safe because neither uses the other at load
import {upsertByKey} from './portal';

/** The function a model code names, when its package is installed. */
export function funcOf(code: string): DG.Func | null {
  const [pkg, name] = code.split(':');
  return (pkg && name ? DG.Func.find({package: pkg, name})[0] : DG.Func.byName(code)) ?? null;
}

/** A number from the property, or from its annotation options (`units`, `min`, `max` live there). */
const num = (...xs: unknown[]): number | null => {
  for (const x of xs) {
    const n = typeof x === 'number' ? x : typeof x === 'string' && x.trim() !== '' ? Number(x) : NaN;
    if (!Number.isNaN(n))
      return n;
  }
  return null;
};

const paramInfo = (p: DG.Property): FuncParamInfo => {
  const o = ((p as any).options ?? {}) as Record<string, any>;
  return {name: p.name, caption: p.caption || o.caption || null, units: p.units || o.units || null,
    description: p.description || o.description || null, min: num(p.min, o.min), max: num(p.max, o.max)};
};

/** The function's inputs and outputs as a card shows them. */
export function funcParams(func: DG.Func | null): {inputs: FuncParamInfo[], outputs: FuncParamInfo[]} {
  if (!func)
    return {inputs: [], outputs: []};
  return {inputs: func.inputs.map(paramInfo), outputs: func.outputs.map(paramInfo)};
}

/** Loads a merged catalogue. Idempotent: rows are matched by code, and a model's portfolios and
 * unit operations only grow, so loading one portfolio never drops what another declared. A card
 * that lists no inputs or outputs gets the function's own when its package is installed. */
export async function loadCatalogue(c: MergedCatalogue, report: string[]): Promise<void> {
  if (!c.unitOperations.length && !c.runTypes.length && !c.models.length)
    return;
  // like a model's, the portfolios that declare an entry only grow
  for (const u of c.unitOperations) {
    const existing = await portalDb.unitOperations.getByKey({code: u.code});
    await upsertByKey(portalDb.unitOperations, {code: u.code},
      {name: u.name, ordinal: u.ordinal ?? null, description: u.description ?? null,
        portfolios: unionCodes(existing?.portfolios, u.portfolios ?? []) || null});
  }
  for (const r of c.runTypes) {
    const existing = await portalDb.runTypes.getByKey({code: r.code});
    await upsertByKey(portalDb.runTypes, {code: r.code}, {name: r.name, definition: r.definition ?? null,
      ordinal: r.ordinal ?? null, is_execution_mode: r.isExecutionMode ?? true,
      portfolios: unionCodes(existing?.portfolios, r.portfolios ?? []) || null});
  }
  let cards = 0;
  for (const m of c.models) {
    const existing = await portalDb.models.getByKey({code: m.code});
    const runType = m.defaultRunType ? await portalDb.runTypes.getByKey({code: m.defaultRunType}) : null;
    const modelId = await upsertByKey(portalDb.models, {code: m.code}, {
      name: m.name, package: m.package ?? m.code.split(':')[0], category: m.category ?? null,
      description: m.description ?? null, status: m.status ?? 'active',
      current_version: m.card?.version ?? existing?.current_version ?? null,
      portfolios: unionCodes(existing?.portfolios, m.portfolios) || null,
      unit_operations: unionCodes(existing?.unit_operations, m.unitOperations ?? []) || null,
      default_run_type_id: runType?.id ?? existing?.default_run_type_id ?? null,
    });
    // runs saved before the model was catalogued name it only by model_name: link them now
    const unlinked = (await portalDb.modelRuns.query({filter: `model_name = "${m.code}"`})).filter((r) => !r.model_id);
    for (const r of unlinked)
      await portalDb.modelRuns.update(r.id, {model_id: modelId});
    if (unlinked.length)
      report.push(`model ${m.code}: ${unlinked.length} earlier run(s) linked to the catalogue`);
    if (!m.card)
      continue;
    const func = funcOf(m.code);
    const params = funcParams(func);
    const card = m.card;
    if (func) {
      for (const problem of cardMismatches(card, {inputs: params.inputs.map((p) => p.name),
        outputs: params.outputs.map((p) => p.name)}))
        report.push(`model ${m.code} ${card.version}: ${problem}`);
    }
    // a card the approval workflow (F-30) has acted on keeps its approval: the file does not undo it
    const stored = await portalDb.modelCards.getByKey({model_id: modelId, model_version: card.version});
    const reviewed = stored ? (await portalDb.approvalRequests.query({filter: `kind = "model_card" and target = "${stored.id}"`})
      .catch(() => [])).length > 0 : false;
    await upsertByKey(portalDb.modelCards, {model_id: modelId, model_version: card.version}, {
      question: card.question ?? null, intended_use: card.intendedUse ?? null,
      inputs: storedParams(func || card.inputs?.length ? cardParams(card.inputs, params.inputs) : null),
      outputs: storedParams(func || card.outputs?.length ? cardParams(card.outputs, params.outputs) : null),
      applicability: card.applicability ?? null, performance: card.performance ?? null,
      lineage: card.lineage ?? null, supersedes: card.supersedes ?? null, limitations: card.limitations ?? null,
      owner: card.owner ?? null,
      ...(reviewed ? {} : {approval: card.approval ?? 'draft',
        approved_on: card.approvedOn ? new Date(card.approvedOn).toISOString() : null}),
      applies_to: card.appliesTo ?? null,
    });
    if (reviewed && stored?.approval !== (card.approval ?? 'draft'))
      report.push(`model ${m.code} ${card.version}: approval kept as ${stored!.approval} (set by the approval workflow, not the file)`);
    cards++;
  }
  report.push(`model catalogue: ${c.unitOperations.length} unit operation(s), ${c.runTypes.length} run type(s), ` +
    `${c.models.length} model(s), ${cards} card(s)`);
}

/** Loads catalogues handed over without a portfolio file, each optionally for a portfolio. */
export async function loadCatalogues(sources: {portfolio?: string | null, catalogue: ModelCatalogueConfig}[]):
  Promise<string> {
  const merged = mergeCatalogues(sources.map((s) => ({portfolio: s.portfolio ?? null, catalogue: s.catalogue})));
  if (merged.errors.length > 0)
    throw new Error(merged.errors.join('\n'));
  const report: string[] = [];
  await loadCatalogue(merged.catalogue, report);
  return report.join('\n');
}

/** A model with everything its card shows. */
export interface ModelCardView {
  model: ModelRow;
  /** The current version's card, else the newest; null when the model has none. */
  card: ModelCardRow | null;
  /** Every card of the model, newest first: the version lineage. */
  versions: ModelCardRow[];
  inputs: CardParam[];
  outputs: CardParam[];
  unitOperations: UnitOperationRow[];
  defaultRunType: RunTypeRow | null;
  portfolios: string[];
  /** Whether the function is installed on this instance. */
  installed: boolean;
}

/** The model a code names (a function's nqName), with its card; null when it is not catalogued. */
export async function modelCardOf(code: string): Promise<ModelCardView | null> {
  const model = await portalDb.models.getByKey({code});
  if (!model)
    return null;
  const versions = (await portalDb.modelCards.query({filter: `model_id = "${model.id}"`}))
    .sort((a, b) => b.created_on.valueOf() - a.created_on.valueOf());
  const card = versions.find((v) => v.model_version === model.current_version) ?? versions[0] ?? null;
  const func = funcOf(code);
  const params = funcParams(func);
  const opCodes = splitCodes(model.unit_operations);
  const ops = (await Promise.all(opCodes.map((op) => portalDb.unitOperations.getByKey({code: op}))))
    .filter((x): x is UnitOperationRow => !!x);
  const runType = model.default_run_type_id ? await portalDb.runTypes.get(model.default_run_type_id).catch(() => null) : null;
  return {model, card, versions,
    inputs: cardParams(readParams(card?.inputs), params.inputs),
    outputs: cardParams(readParams(card?.outputs), params.outputs),
    unitOperations: ops, defaultRunType: runType, portfolios: splitCodes(model.portfolios), installed: !!func};
}
