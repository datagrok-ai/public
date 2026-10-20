/** Predictions (F-22, FR-58/103): what a catalogued predictive model said about each compound of a
 * program, per model version, kept as `prediction` rows beside the run that produced them. A batch
 * scoring runs the model once over the program dataset, saves the run, and lands one value per
 * compound and predicted column. The card's `appliesTo.predicts` names the program-data column a
 * prediction stands in for, so predicted and measured meet on it. Pure: the caller runs the model
 * and reads and writes the rows. */
import {columnValue, DataRow, isBlank, rowId} from './conditions';
import type {MissingAssay} from './recommend';

/** One stored prediction, as the rules need it. */
export interface Prediction {
  compound: string;
  /** The program-data column it stands in for. */
  predicts: string;
  model: string;
  /** The model's name, for display; the code when unknown. */
  modelName?: string;
  version: string;
  value: number | null;
  units?: string | null;
  /** The saved run (model_run id) that produced it. */
  modelRunId?: string | null;
  /** ISO time. */
  on?: string | null;
}

const norm = (s: string) => s.trim().toLowerCase();

/** The output column holding the prediction for [predicts]: `predicted_x`, `estimated_x`, `x_predicted`,
 * `x_pred`, `x_estimate`, `pred_x`; null when the model's output has none. The measured column itself,
 * which a model may copy beside its estimate, is never taken. */
export function predictionColumn(columns: string[], predicts: string): string | null {
  const p = norm(predicts);
  const wanted = [`predicted_${p}`, `estimated_${p}`, `${p}_predicted`, `${p}_pred`, `${p}_estimate`, `pred_${p}`];
  for (const w of wanted) {
    const hit = columns.find((c) => norm(c) === w);
    if (hit)
      return hit;
  }
  return null;
}

/** The output column naming the compound: compound, compound_id, id, code or name. */
export function compoundColumn(columns: string[]): string | null {
  for (const w of ['compound', 'compound_id', 'compound_code', 'id', 'code', 'name']) {
    const hit = columns.find((c) => norm(c) === w);
    if (hit)
      return hit;
  }
  return null;
}

const UNIT_SUFFIXES: [RegExp, string][] = [[/_um$/i, 'uM'], [/_nm$/i, 'nM'], [/_mm$/i, 'mM'], [/_pct$/i, '%'],
  [/_mg_ml$/i, 'mg/mL'], [/_ug_ml$/i, 'ug/mL'], [/_h$/i, 'h']];

/** Units read off a column name by its suffix (`solubility_um` → uM); null when it says none. */
export function unitsOf(column: string): string | null {
  return UNIT_SUFFIXES.find(([re]) => re.test(column))?.[1] ?? null;
}

const num = (v: unknown): number | null => typeof v === 'number' && Number.isFinite(v) ? v :
  typeof v === 'string' && v.trim() !== '' && Number.isFinite(Number(v)) ? Number(v) : null;

/** The compounds of the program dataset, in order, once each. */
export function compoundsOf(rows: DataRow[]): string[] {
  return [...new Set(rows.map((r, i) => rowId(r, i)))];
}

/** The compounds a model version has not scored yet for [predicts]. */
export function unscored(compounds: string[], stored: Prediction[], model: string, version: string, predicts: string): string[] {
  const done = new Set(stored.filter((p) => p.model === model && p.version === version && norm(p.predicts) === norm(predicts))
    .map((p) => p.compound));
  return compounds.filter((c) => !done.has(c));
}

export interface ScoringPlan {
  /** New predictions, one per compound. */
  insert: {compound: string, value: number}[];
  /** Compounds this version had scored already: kept, not rewritten. */
  kept: number;
  /** Rows of the output without a compound or a number. */
  skipped: number;
}

/** Which of a run's [scored] values to store: one per compound this model version has not scored
 * yet. What is stored stays: scoring the same version again changes nothing. */
export function planScoring(scored: {compound: unknown, value: unknown}[], stored: Prediction[], model: string,
  version: string, predicts: string): ScoringPlan {
  const done = new Set(stored.filter((p) => p.model === model && p.version === version && norm(p.predicts) === norm(predicts))
    .map((p) => p.compound));
  const plan: ScoringPlan = {insert: [], kept: 0, skipped: 0};
  const seen = new Set<string>();
  for (const s of scored) {
    const compound = isBlank(s.compound) ? null : String(s.compound);
    const value = num(s.value);
    if (compound === null || value === null) {
      plan.skipped++;
      continue;
    }
    if (seen.has(compound))
      continue;
    seen.add(compound);
    if (done.has(compound))
      plan.kept++;
    else
      plan.insert.push({compound, value});
  }
  return plan;
}

/** Per compound and predicted column, the prediction to show: the newest version, then the newest. */
export function latestPredictions(stored: Prediction[]): Prediction[] {
  const by = new Map<string, Prediction>();
  const newer = (a: Prediction, b: Prediction) => compareVersions(a.version, b.version) > 0 ||
    (compareVersions(a.version, b.version) === 0 && String(a.on ?? '') > String(b.on ?? ''));
  for (const p of stored) {
    const k = `${p.model}|${norm(p.predicts)}|${p.compound}`;
    const cur = by.get(k);
    if (!cur || newer(p, cur))
      by.set(k, p);
  }
  return [...by.values()].sort((a, b) => a.compound.localeCompare(b.compound) || a.predicts.localeCompare(b.predicts) ||
    a.model.localeCompare(b.model));
}

/** Dotted versions by number where they are numbers (1.10 after 1.9), else as text. */
export function compareVersions(a: string, b: string): number {
  const x = String(a ?? '').split('.');
  const y = String(b ?? '').split('.');
  for (let i = 0; i < Math.max(x.length, y.length); i++) {
    const p = x[i] ?? '';
    const q = y[i] ?? '';
    const d = /^\d+$/.test(p) && /^\d+$/.test(q) ? Number(p) - Number(q) : p.localeCompare(q);
    if (d !== 0)
      return d;
  }
  return 0;
}

export interface PredictedVsMeasured {
  model: string;
  modelName: string;
  version: string;
  predicts: string;
  units: string | null;
  /** One per compound predicted, measured value or null. */
  pairs: {compound: string, predicted: number, measured: number | null}[];
  /** Compounds with both. */
  paired: number;
  /** Mean absolute error over the pairs; null without one. */
  mae: number | null;
  /** Pearson's r over the pairs; null with fewer than three or no spread. */
  r: number | null;
}

/** Joins one model's latest predictions for a column with the program dataset's measured values of it. */
export function predictedVsMeasured(latest: Prediction[], rows: DataRow[]): PredictedVsMeasured[] {
  const groups = new Map<string, Prediction[]>();
  for (const p of latest) {
    const k = `${p.model}|${p.predicts}`;
    groups.set(k, [...(groups.get(k) ?? []), p]);
  }
  const out: PredictedVsMeasured[] = [];
  for (const ps of groups.values()) {
    const first = ps[0];
    const measuredOf = new Map<string, number>();
    rows.forEach((r, i) => {
      const v = num(columnValue(r, first.predicts));
      if (v !== null && !measuredOf.has(rowId(r, i)))
        measuredOf.set(rowId(r, i), v);
    });
    const pairs = ps.filter((p) => p.value !== null).map((p) => ({compound: p.compound, predicted: p.value!,
      measured: measuredOf.get(p.compound) ?? null}));
    const both = pairs.filter((p) => p.measured !== null) as {compound: string, predicted: number, measured: number}[];
    const mae = both.length ? both.reduce((s, p) => s + Math.abs(p.predicted - p.measured), 0) / both.length : null;
    out.push({model: first.model, modelName: first.modelName ?? first.model, version: ps.map((p) => p.version)
      .sort(compareVersions).pop()!, predicts: first.predicts, units: first.units ?? unitsOf(first.predicts), pairs,
    paired: both.length, mae, r: pearson(both.map((p) => p.predicted), both.map((p) => p.measured))});
  }
  return out.sort((a, b) => a.predicts.localeCompare(b.predicts) || a.model.localeCompare(b.model));
}

function pearson(x: number[], y: number[]): number | null {
  const n = x.length;
  if (n < 3)
    return null;
  const mx = x.reduce((s, v) => s + v, 0) / n;
  const my = y.reduce((s, v) => s + v, 0) / n;
  let sxy = 0;
  let sxx = 0;
  let syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
  }
  return sxx === 0 || syy === 0 ? null : sxy / Math.sqrt(sxx * syy);
}

const round = (v: number) => Math.abs(v) >= 100 ? v.toFixed(0) : Math.abs(v) >= 1 ? v.toFixed(1) : v.toPrecision(2);

/** `48.4 uM (Solubility estimate (example) 1.0)` */
export function describePrediction(p: Prediction): string {
  const units = p.units ?? unitsOf(p.predicts);
  return `${p.value === null ? '—' : round(p.value)}${units ? ` ${units}` : ''} (${p.modelName ?? p.model} ${p.version})`;
}

/** `6 of 11 measured; mean absolute error 12.3 uM; r 0.81` */
export function describeFit(f: PredictedVsMeasured): string {
  return `${f.paired} of ${f.pairs.length} measured` +
    (f.mae !== null ? `; mean absolute error ${round(f.mae)}${f.units ? ` ${f.units}` : ''}` : '') +
    (f.r !== null ? `; r ${f.r.toFixed(2)}` : '');
}

/** A missing assay with the stored predictions that stand in for it, per compound lacking it. */
export type MissingWithPredictions = MissingAssay & {predicted: Prediction[]};

/** Stored predictions first: for each missing column, the latest prediction of each compound that
 * lacks it. A compound predicted by several models shows each. */
export function withPredictions(missing: MissingAssay[], latest: Prediction[]): MissingWithPredictions[] {
  return missing.map((m) => ({...m, predicted: m.column === null ? [] :
    latest.filter((p) => norm(p.predicts) === norm(m.column!) && m.items.includes(p.compound) && p.value !== null)}));
}
