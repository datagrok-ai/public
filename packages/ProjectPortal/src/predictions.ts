/** Predictions against the platform (F-22): batch scoring runs a catalogued predictive model once
 * over a program's dataset, saves the call as a model run (so provenance leads prediction → run →
 * model), and stores one `prediction` row per compound and predicted column, idempotently per model
 * version. Predicted meets measured on the column the card's `appliesTo.predicts` names. The rules
 * live in ./model/predictions. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {textColumn} from './columns';
import dayjs from 'dayjs';
import {ModelRow, portalDb, PortfolioRow, PredictionRow, ProgramRow} from './generated/db';
import {compoundsOf as programCompoundsOf, portfolioByCode, programData} from './portal';
import {funcOf} from './catalogue';
import {ModelApplicability, prefillInputs} from './model/recommend';
import {columnValue, rowId} from './model/conditions';
import {programTable} from './program-data';
import {saveModelRun} from './model-runs';
import {
  compoundColumn, compoundsOf, latestPredictions, planScoring, Prediction, predictionColumn, unitsOf, unscored,
} from './model/predictions';
import {filterValue, ROW_LIMIT} from './rows';

const iso = (d: any): string | null => d ? dayjs(d).toISOString() : null;

/** A stored row as the rules read it. */
export const asPrediction = (p: PredictionRow, names?: Map<string, string>): Prediction => ({
  compound: p.compound_code, predicts: p.predicts, model: p.model_code, modelName: names?.get(p.model_code) ?? p.model_code,
  version: p.model_version, value: p.predicted_value ?? null, units: p.units ?? null, modelRunId: p.model_run_id ?? null,
  on: iso(p.predicted_on ?? p.created_on)});

/** Every prediction stored for a program, the models' names alongside. */
export async function predictionsOf(programId: string): Promise<{rows: PredictionRow[], latest: Prediction[]}> {
  const rows = await portalDb.predictions.query({filter: `program_id = ${filterValue(programId)}`, limit: ROW_LIMIT});
  const codes = [...new Set(rows.map((r) => r.model_code))];
  const models = codes.length ? await portalDb.models.query({filter: `code in (${codes.map((c) => `"${c}"`).join(', ')})`, limit: ROW_LIMIT})
    .catch(() => [] as ModelRow[]) : [];
  const names = new Map(models.map((m) => [m.code, m.name]));
  return {rows, latest: latestPredictions(rows.map((r) => asPrediction(r, names)))};
}

/** The columns a model's current card says it predicts, and the version that card is for. */
async function cardOf(model: ModelRow): Promise<{version: string | null, predicts: string[]}> {
  const cards = await portalDb.modelCards.query({filter: `model_id = ${filterValue(model.id)}`, limit: ROW_LIMIT});
  const card = cards.find((c) => c.model_version === model.current_version) ?? cards[0];
  return {version: card?.model_version ?? model.current_version ?? null,
    predicts: ((card?.applies_to as ModelApplicability | undefined)?.predicts ?? []).filter((p) => !!p)};
}

export interface ScoringResult {
  portfolio: string;
  program: string;
  model: string;
  version: string;
  /** The saved run; null when nothing needed scoring. */
  runId: string | null;
  modelRunId: string | null;
  /** New predictions stored, per predicted column. */
  stored: {predicts: string, count: number}[];
  /** Compounds this version had scored already. */
  kept: number;
  notes: string[];
  summary: string;
}

/** Scores the program's compounds with a catalogued model whose card says what it predicts. The
 * model runs once over the whole program dataset; only compounds this model version has not scored
 * are stored, so scoring again with the same version runs nothing and changes nothing. */
export async function scoreCompounds(portfolioCode: string, programCode: string, modelCode: string): Promise<ScoringResult> {
  const portfolio = await portfolioByCode(portfolioCode);
  const program = await programOf(portfolio, programCode);
  const model = await portalDb.models.getByKey({code: modelCode});
  if (!model)
    throw new Error(`${modelCode} is not in the model catalogue`);
  const {version: cardVersion, predicts} = await cardOf(model);
  if (!predicts.length)
    throw new Error(`${model.name}: its card names no program-data column it predicts (appliesTo.predicts)`);
  const func = funcOf(model.code);
  if (!func)
    throw new Error(`${model.code} is not installed on this instance`);
  const version = cardVersion ?? (func as any).package?.version ?? '?';
  const rows = await programData(portfolio, program.code);
  if (!rows)
    throw new Error(`${portfolio.code} declares no program data to score`);
  const compounds = compoundsOf(rows);
  const stored = (await portalDb.predictions.query({filter: `program_id = ${filterValue(program.id)} and model_code = ${filterValue(model.code)}`, limit: ROW_LIMIT}))
    .map((p) => asPrediction(p));
  const base = {portfolio: portfolio.code, program: program.code, model: model.code, version};
  if (predicts.every((p) => unscored(compounds, stored, model.code, version, p).length === 0)) {
    return {...base, runId: null, modelRunId: null, stored: [], kept: compounds.length, notes: [],
      summary: `${program.code}: all ${compounds.length} compounds already scored by ${model.name} ${version}; nothing run`};
  }

  // one run over the whole dataset, the dialog's own prefill: the table, and the program's codes
  const t = await programTable(portfolio.code, program.code);
  const {values, needsData} = prefillInputs(func.inputs.map((p) => ({name: p.name, type: p.propertyType})),
    {portfolio: portfolio.code, program: program.code});
  for (const name of needsData)
    values[name] = t.df;
  const call = func.prepare(values);
  await call.call();
  const out = func.outputs.map((p) => call.getParamValue(p.name)).find((v) => v instanceof DG.DataFrame) as DG.DataFrame | undefined;
  if (!out)
    throw new Error(`${model.name} returned no table`);
  const names = out.columns.names();
  const idCol = compoundColumn(names);
  if (!idCol)
    throw new Error(`${model.name}'s output has no compound column (${names.join(', ')})`);

  const plans = predicts.map((p) => {
    const col = predictionColumn(names, p);
    const scored = col ? Array.from({length: out.rowCount}, (_, i) => ({compound: out.get(idCol, i),
      value: out.col(col)!.isNone(i) ? null : out.get(col, i)})) : [];
    return {predicts: p, column: col, plan: planScoring(scored, stored, model.code, version, p)};
  });
  const notes = plans.filter((x) => !x.column).map((x) => `no output column predicts ${x.predicts} (looked for predicted_${x.predicts}, ` +
    `estimated_${x.predicts}, ${x.predicts}_pred, ...)`);
  if (plans.every((x) => !x.column))
    throw new Error(`${model.name}: ${notes.join('; ')}`);

  // nothing new to store (estimates missing, or ids the model writes differently): no run is saved,
  // else every call would save another empty one
  const kept = Math.max(0, ...plans.map((x) => x.plan.kept));
  if (plans.every((x) => x.plan.insert.length === 0)) {
    const skipped = Math.max(0, ...plans.map((x) => x.plan.skipped));
    return {...base, runId: null, modelRunId: null, stored: [], kept, notes,
      summary: `${program.code}: ${model.name} ${version} gave nothing new to store` +
        (kept ? `; ${kept} already scored, kept` : '') + (skipped ? `; ${skipped} without an estimate or a compound` : '') +
        '; no run saved' + (notes.length ? `\n  ${notes.join('\n  ')}` : '')};
  }

  // the run first, so each prediction can point at it
  const saved = await saveModelRun(call, portfolio.code, program.code, {purpose: 'screening',
    title: `Scoring · ${model.name} ${version} · ${program.code}`,
    description: `Batch scoring of ${compounds.length} compounds; predictions stored for ${plans.filter((x) => x.column)
      .map((x) => x.predicts).join(', ')}`});
  notes.push(...saved.notes);
  // the program's own compounds, by code: the compound table is closed to most users
  const known = new Map(!compounds.length ? [] : (await programCompoundsOf(program.id).catch(() => []))
    .map((c) => [String(c.compound.code), c.compound.id as string]));
  const now = new Date().toISOString();
  const me = DG.User.current().id;
  const result: ScoringResult = {...base, runId: saved.runId, modelRunId: saved.modelRun.id, stored: [], kept: 0, notes, summary: ''};
  for (const x of plans.filter((y) => y.column)) {
    if (x.plan.insert.length) {
      await portalDb.predictions.insert(x.plan.insert.map((v) => ({program_id: program.id, compound_id: known.get(v.compound),
        compound_code: v.compound, model_id: model.id, model_code: model.code, model_version: version, predicts: x.predicts,
        predicted_value: v.value, units: unitsOf(x.predicts) ?? undefined, model_run_id: saved.modelRun.id, run_id: saved.runId,
        predicted_on: now as any, predicted_by: me})));
    }
    result.stored.push({predicts: x.predicts, count: x.plan.insert.length});
    result.kept = Math.max(result.kept, x.plan.kept);
  }
  result.summary = `${program.code}: ${model.name} ${version} scored ${compounds.length} compounds in one run; stored ` +
    result.stored.map((s) => `${s.count} ${s.predicts}`).join(', ') + (result.kept ? `; ${result.kept} already scored, kept` : '') +
    (notes.length ? `\n  ${notes.join('\n  ')}` : '');
  return result;
}

async function programOf(portfolio: PortfolioRow, code: string): Promise<ProgramRow> {
  const p = (await portalDb.programs.query({filter: `portfolio_id = ${filterValue(portfolio.id)} and code = ${filterValue(code)}`, limit: ROW_LIMIT}))[0];
  if (!p)
    throw new Error(`${portfolio.code} has no program ${code}, or it is not shared with you`);
  return p;
}

/** The program's predictions as a table: one row per compound, column and model version, the latest
 * first, with the measured value where the program dataset has one. */
export async function predictionsFrame(portfolioCode: string, programCode: string, modelCode?: string | null): Promise<DG.DataFrame> {
  const portfolio = await portfolioByCode(portfolioCode);
  const program = await programOf(portfolio, programCode);
  const {latest} = await predictionsOf(program.id);
  const shown = latest.filter((p) => !modelCode || p.model === modelCode);
  const rows = await programData(portfolio, program.code).catch(() => null) ?? [];
  const measured = (p: Prediction) => {
    const i = rows.findIndex((x, j) => rowId(x, j) === p.compound);
    const v = i < 0 ? null : columnValue(rows[i], p.predicts);
    return typeof v === 'number' ? v : null;
  };
  const df = DG.DataFrame.fromColumns([
    textColumn('compound', shown.map((p) => p.compound)),
    textColumn('predicts', shown.map((p) => p.predicts)),
    textColumn('model', shown.map((p) => p.model)),
    textColumn('version', shown.map((p) => p.version)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'predicted', shown.map((p) => p.value)),
    DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, 'measured', shown.map(measured)),
    textColumn('units', shown.map((p) => p.units ?? '')),
    textColumn('model run', shown.map((p) => p.modelRunId ?? '')),
  ]);
  df.name = `${program.code} · Predictions`;
  return df;
}

/** Opens predicted against measured for one model and column as a table view with a scatter plot. */
export async function openPredictedVsMeasured(portfolioCode: string, programCode: string, modelCode: string,
  predicts: string): Promise<DG.TableView> {
  const all = await predictionsFrame(portfolioCode, programCode, modelCode);
  const df = all.clone(DG.BitSet.create(all.rowCount, (i) => all.get('predicts', i) === predicts));
  df.name = `${programCode} · ${predicts}: predicted vs measured`;
  const v = grok.shell.addTableView(df);
  v.scatterPlot({x: 'measured', y: 'predicted', color: 'version'});
  return v;
}
