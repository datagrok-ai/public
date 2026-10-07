/** Model recommendation against the platform (F-13): the catalogue and its cards from the portal's
 * tables, the program's place and readiness from its context, and Run, which opens a model's
 * function dialog filled in from the program. The rules live in ./model/recommend. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as ui from 'datagrok-api/ui';
import {portalDb, PortfolioRow, ProgramRow, StageRow} from './generated/db';
import {enabledWorkUnits, portfolioByCode, programContext, ProgramContext, programData, readinessFor} from './portal';
import {funcOf} from './catalogue';
import {splitCodes} from './model/catalogue';
import {
  CandidateModel, missingAssays, MissingAssay, ModelApplicability, predictorsOf, prefillInputs, recommend, Recommendation,
} from './model/recommend';
import {DataRow, DEFAULT_LEAD_ROLES, ProgramDataConfig} from './model/conditions';
import {Readiness} from './model/readiness';
import {programTable} from './program-data';
import {describePrediction, Prediction, withPredictions} from './model/predictions';
import {predictionsOf} from './predictions';
import {filterValue, ROW_LIMIT} from './rows';

/** The catalogued models of [portfolio] (every catalogued model when omitted) with their current
 * card's approval, applicability, question and unit operations. Another portfolio's models are not
 * read into the list: the model tables are readable by all users, so the scoping is the app's. */
export async function candidateModels(portfolio?: string): Promise<CandidateModel[]> {
  const [all, cards] = await Promise.all([portalDb.models.query({limit: ROW_LIMIT}), portalDb.modelCards.query({limit: ROW_LIMIT})]);
  const models = portfolio ? all.filter((m) => splitCodes(m.portfolios).includes(portfolio)) : all;
  return models.map((m) => {
    const mine = cards.filter((c) => c.model_id === m.id);
    const card = mine.find((c) => c.model_version === m.current_version) ?? mine[0];
    return {code: m.code, name: m.name, status: m.status ?? null, approval: card?.approval ?? null,
      portfolios: splitCodes(m.portfolios), appliesTo: (card?.applies_to as ModelApplicability | undefined) ?? null,
      installed: !!funcOf(m.code), unitOperations: splitCodes(m.unit_operations), question: card?.question ?? null};
  });
}

export interface ProgramRecommendation {
  portfolio: PortfolioRow;
  program: ProgramRow;
  /** The stage whose gate the missing assays are read from: the one asked for, else the next. */
  stage: StageRow | null;
  series: string | null;
  missing: MissingAssay[];
  models: Recommendation[];
  /** Per missing column, the suggested models that predict it. */
  predictors: Map<string, Recommendation[]>;
  /** The program's stored predictions, the latest per compound, column and model: shown first. */
  predictions?: Prediction[];
  dataError?: string;
}

/** The models that fit a program (or one of its series) and the assays its gate still lacks.
 * [readiness] and [rows] spare the program page a second read of what it already holds. */
export async function recommendFor(ctx: ProgramContext, opts: {stageId?: string | null, series?: string | null,
  task?: string | null, readiness?: Readiness | null, rows?: DataRow[] | null, dataError?: string} = {}): Promise<ProgramRecommendation> {
  const {program, portfolio, stages} = ctx;
  const stageId = opts.stageId ?? opts.readiness?.stageId ?? ctx.readiness?.stageId ?? null;
  const stage = stages.find((s) => s.id === stageId) ?? null;
  let rows = opts.rows;
  let dataError = opts.dataError;
  const open = ctx.criteria.filter((c) => c.stageId === stageId && c.condition);
  if (rows === undefined && open.length) {
    try {
      rows = await programData(portfolio, program.code);
    } catch (e: any) {
      rows = null;
      dataError = e?.message ?? String(e);
    }
  }
  const readiness = opts.readiness && opts.readiness.stageId === stageId ? opts.readiness :
    stageId ? readinessFor(ctx, stageId, rows ?? null, dataError) : null;
  const byId = new Map(ctx.criteria.map((c) => [c.id, c]));
  const leadRoles = (portfolio.program_data as ProgramDataConfig | undefined)?.leadRoles ?? DEFAULT_LEAD_ROLES;
  const missing = readiness ? missingAssays(readiness.criteria.map((c) => ({...c, condition: byId.get(c.id)?.condition,
    rule: byId.get(c.id)?.rule})), rows ?? null, {leadRoles, series: opts.series}) : [];
  const current = stages.find((s) => s.id === program.stage_id);
  const models = recommend(await candidateModels(portfolio.code), {
    portfolio: portfolio.code, modality: program.modality,
    stages: [current?.code, stage?.code].filter((s): s is string => !!s),
    workUnits: opts.series ? ['series'] : enabledWorkUnits(portfolio),
    task: opts.task, missingColumns: missing.map((m) => m.column).filter((c): c is string => !!c),
    unitOperations: await unitOperationsOf(program.id),
  });
  const predictions = await predictionsOf(program.id).then((p) => p.latest).catch(() => [] as Prediction[]);
  return {portfolio, program, stage, series: opts.series ?? null, missing, models,
    predictors: predictorsOf(missing, models), predictions, dataError};
}

/** The unit operations a program works in: those of its process stages. Empty when it has none. */
async function unitOperationsOf(programId: string): Promise<string[]> {
  const stages = await portalDb.processStages.query({filter: `program_id = ${filterValue(programId)}`, limit: ROW_LIMIT}).catch(() => []);
  const ids = [...new Set(stages.map((s) => s.unit_operation_id).filter((x): x is string => !!x))];
  const ops = await Promise.all(ids.map((id) => portalDb.unitOperations.get(id).catch(() => null)));
  return ops.filter((o) => o).map((o) => o!.code);
}

/** The same, by codes, for the package function. */
export async function recommendByCode(portfolioCode: string, programCode: string,
  opts: {stage?: string, series?: string, task?: string} = {}): Promise<ProgramRecommendation> {
  const portfolio = await portfolioByCode(portfolioCode);
  const p = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)} and code = ${filterValue(programCode)}`, limit: ROW_LIMIT}))[0];
  if (!p)
    throw new Error(`${portfolioCode} has no program ${programCode}, or it is not shared with you`);
  const ctx = await programContext(p.id);
  let stageId: string | null = null;
  if (opts.stage) {
    const s = ctx.stages.find((x) => x.code === opts.stage);
    if (!s)
      throw new Error(`${opts.stage} is not a stage of ${programCode}'s model`);
    stageId = s.id;
  }
  return recommendFor(ctx, {stageId, series: opts.series || null, task: opts.task || null});
}

/** One row per suggested model: the stable contract of recommendModels. */
export function recommendationFrame(r: ProgramRecommendation): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('model', r.models.map((m) => m.code)),
    DG.Column.fromStrings('name', r.models.map((m) => m.name)),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'rank', r.models.map((_, i) => i + 1)),
    DG.Column.fromStrings('predicts', r.models.map((m) => m.predicts.join(', '))),
    DG.Column.fromStrings('matched', r.models.map((m) => m.matched.join(', '))),
    DG.Column.fromStrings('reasons', r.models.map((m) => m.reasons.join('; '))),
    DG.Column.fromStrings('approval', r.models.map((m) => m.approval ?? '')),
    DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'installed', r.models.map((m) => m.installed)),
    DG.Column.fromStrings('question', r.models.map((m) => m.question ?? '')),
  ]);
  df.name = `${r.program.code} · Suggested models`;
  return df;
}

/** One row per missing column (or per criterion with only a rule in words). */
export function missingAssayFrame(r: ProgramRecommendation): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('stage', r.missing.map(() => r.stage?.code ?? '')),
    DG.Column.fromStrings('criterion', r.missing.map((m) => m.criterion)),
    DG.Column.fromList(DG.COLUMN_TYPE.BOOL, 'required', r.missing.map((m) => m.required)),
    DG.Column.fromStrings('assay', r.missing.map((m) => m.column ?? '')),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'missing', r.missing.map((m) => m.missing)),
    DG.Column.fromStrings('items', r.missing.map((m) => m.items.join(', '))),
    DG.Column.fromStrings('predicted by', r.missing.map((m) => m.column ?
      (r.predictors.get(m.column) ?? []).map((x) => x.code).join(', ') : '')),
    DG.Column.fromStrings('predicted', r.missing.map((m) => withPredictions([m], r.predictions ?? [])[0].predicted
      .map((p) => `${p.compound} ${describePrediction(p)}`).join('; '))),
  ]);
  df.name = `${r.program.code} · Missing assays`;
  return df;
}

/** Opens [code]'s function dialog filled in from the program: a dataframe input gets the program
 * dataset (opened as a table, so the dialog can offer it), inputs named program, portfolio or
 * series their codes. */
export async function runModelFor(code: string, ctx: {portfolio: string, program: string, series?: string | null}):
  Promise<DG.FuncCall> {
  const func = funcOf(code);
  if (!func)
    throw new Error(`${code} is not installed on this instance`);
  const {values, needsData} = prefillInputs(func.inputs.map((p) => ({name: p.name, type: p.propertyType})), ctx);
  if (needsData.length) {
    const t = await programTable(ctx.portfolio, ctx.program);
    grok.shell.addTableView(t.df);
    for (const name of needsData)
      values[name] = t.df;
  }
  const call = func.prepare(values);
  // the platform's own call.edit() runs a call whose inputs are all set without asking, so the
  // dialog is built here: the standard editor, then Run
  const editor = document.createElement('div');
  await call.buildEditor(editor, {showTableSelectors: true});
  const dlg = ui.dialog({title: func.friendlyName || func.name}).add(editor);
  dlg.root.dataset.u2Name = 'run-model-dialog';
  dlg.onOK(async () => {
    try {
      await call.call(true);
      for (const v of call.getResultViews())
        grok.shell.addView(v);
    } catch (e: any) {
      grok.shell.error(`${func.friendlyName || func.name}: ${e?.message ?? e}`);
    }
  });
  dlg.show();
  return call;
}
