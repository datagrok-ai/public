// exposure per group from the NCA package, when it is installed; nothing here is saved or enters the analysis
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {drugKey, parseArmName, therapyOf} from './arm-name';
import type {ArmPart} from './arm-name';
import type {Therapy} from './archetypes';
import {crux, rows as cruxRows} from './crux';
import {finite, median, SEED} from './engine';
import {FIRST_DOSE_DAY, readDayOf} from './study-effect';
import type {PreviousStudy} from './study-effect';
import {dosesPerDay, readSchedule} from './schedule';
import type {DosingRules} from './schedule';
import {armParts} from './study-meta';
import type {StudyMeta} from './study-meta';
import type {CtdStudy} from './ctd';

/** A group's exposure as NCA computed it: AUClast and Cmax (mean and CV%, as NCA labels and units them). */
export interface Exposure {
  auc: number, aucCv: number | null, aucUnit: string, cmax: number | null, cmaxCv: number | null, cmaxUnit: string,
  /** NCA's mean AUCinf and its unit, null when NCA gives none (too few points on the terminal phase) */
  aucInf: number | null, aucInfUnit: string,
  /** the animals (profiles) behind AUCinf and AUClast, when NCA counts them */
  nAucInf?: number | null, nAucLast?: number | null,
  /** the dose and route NCA read, and the arm's own schedule: for the tooltip, never to scale */
  dose: string, route: string, schedule: string | null, profiles: number,
  /** the drug NCA names (the group's name without a drug column); a pooled profile of animals sampled sparsely */
  drug: string, pooled: boolean,
}
/** An arm: its exposure, or why it has none, in a plain sentence. */
export type ArmExposure = {exposure: Exposure, why?: undefined} | {exposure?: undefined, why: string};

const nca = (name: string) => DG.Func.find({package: 'Nca', name})[0] ?? null;
/** Whether the NCA package is installed (its functions are looked up, never imported). */
export const ncaInstalled = () => nca('ncaCompute') !== null;

// a column of the concentration table that names the group as the study does
const GROUP = /^(group|groups|arm|treatment|trt|grp|cohort|dosegroup)$/;
const norm = (s: string) => s.toLowerCase().replace(/[\s_-]+/g, '');
// NCA's routes (IV-bolus, IV-infusion, PO, SC, IM, IP, other) and what an arm's name says
const ROUTE_OF: {[arm: string]: string} = {'i.v.': 'IV', 'p.o.': 'PO', 's.c.': 'SC', 'i.m.': 'IM', 'i.p.': 'IP',
  'i.t.': 'other'};
const sameRoute = (arm: string | null, theirs: string) => arm === null ||
  (ROUTE_OF[arm] ?? arm) === (theirs.startsWith('IV') ? 'IV' : theirs);
const unitOf = (u: string) => u.replace(/[µμ]/g, 'u').toLowerCase();

/** NCA's mapping of `table`, or null when it does not read as concentrations over time; `issues`, what NCA missed. */
export async function concentrationMapping(table: DG.DataFrame, issues: string[] = []): Promise<any | null> {
  const detect = nca('ncaDetectMapping');
  if (!detect)
    return null;
  try {
    const d = JSON.parse(await detect.apply({table}));
    issues.push(...Array.isArray(d?.issues) ? d.issues : []);
    return d?.mapping && d.confidence === 'high' ? d.mapping : null;
  } catch {
    return null;
  }
}

/** The group column of a concentration table NCA reads with `mapping` (null: none): a text column it does not use,
 *  named as a group. */
export function groupColumnOf(table: DG.DataFrame, mapping: any): string | null {
  const used = new Set([mapping.subject, mapping.time, mapping.concentration, mapping.dose, mapping.route,
    ...mapping.groupBy ?? []]);
  return table.columns.toList().find((c) => c.type === DG.COLUMN_TYPE.STRING && !used.has(c.name) &&
    GROUP.test(norm(c.name)))?.name ?? null;
}

/** NCA's per-profile results and per-group summary of `table`, its group column (if any) kept as a key. */
export async function runNca(table: DG.DataFrame, mapping: any):
Promise<{results: DG.DataFrame, summary: DG.DataFrame, groupCol: string | null}> {
  const groupCol = groupColumnOf(table, mapping);
  const m = {...mapping, groupBy: [...mapping.groupBy ?? [], ...groupCol ? [groupCol] : []]};
  const results: DG.DataFrame = await nca('ncaCompute')!.apply({table, mapping: JSON.stringify(m), rules: ''});
  const summary: DG.DataFrame = await nca('ncaSummarize')!.apply({results});
  return {results, summary, groupCol};
}

/** NCA's reasons for a group it could not compute, in plain words. */
function plain(reasons: string[]): string {
  const r = reasons.join('; ');
  if (/fewer than \d+ animals/.test(r))
    return 'too few animals sampled (a pooled profile needs at least 3)';
  if (/route column/.test(r))
    return 'NCA does not compute doses given by this route';
  if (/mixed sampling/.test(r))
    return 'some animals sampled fully and others only a few times';
  if (/time points?|λz|lambda/i.test(r))
    return 'too few time points per animal';
  return r || 'no reason given';
}

/** Each arm's exposure from NCA's results and summary: the group of the same name, else (one drug only) the drug
 *  with the same dose and route, read from the study's group details (`meta`), else from the name. Never a guess:
 *  several candidates, or none, is said so. A combination's exposure is its drugs', each apart: exposureByDrug. */
export function exposureByArm(arms: string[], results: DG.DataFrame, summary: DG.DataFrame,
  groupCol: string | null, meta: StudyMeta | null = null): Map<string, ArmExposure> {
  const one = ncaReader(results, summary, groupCol, meta);
  return new Map(arms.map((arm) => [arm, one(arm, null)]));
}

/** Each drug of combination `arm`, its exposure from NCA: the rows of the group of the same name for that drug (the
 *  concentrations name the drug), else the drug at the same dose and route (no group column). */
export function exposureByDrug(arm: string, results: DG.DataFrame, summary: DG.DataFrame, groupCol: string | null,
  meta: StudyMeta | null = null): {part: ArmPart, x: ArmExposure}[] {
  const one = ncaReader(results, summary, groupCol, meta);
  return armParts(arm, meta).map((part) => ({part, x: one(arm, part)}));
}

/** A value read from a float32 column, to the digits it holds (null stays null). */
const float32 = (x: number | null) => x === null ? null : Number(x.toPrecision(7));

/** An arm's exposure in NCA's results and summary (`part`: one drug of a combination arm, null: the arm alone). */
function ncaReader(results: DG.DataFrame, summary: DG.DataFrame, groupCol: string | null,
  meta: StudyMeta | null): (arm: string, part: ArmPart | null) => ArmExposure {
  const keys: string[] = JSON.parse(summary.getTag('nca.keyColumns') || '[]');
  const units = JSON.parse(results.getTag('nca.units') || '{}');
  const doseUnit = typeof units.dose === 'string' ? unitOf(units.dose) : null;
  const drugCol = keys.find((k) => k !== groupCol && !/^(study|period|visit|occasion)/i.test(k)) ?? null;
  const str = (df: DG.DataFrame, c: string, i: number) => df.col(c)?.getString(i) ?? '';
  const num = (df: DG.DataFrame, c: string, i: number): number | null => {
    const col = df.col(c);
    return col && !col.isNone(i) && Number.isFinite(col.get(i)) ? col.get(i) : null;
  };
  // the rows of `df` that are this arm's, or why it cannot be matched; a row NCA blocked has no dose or route, so by
  // drug it is the drug's alone
  const rowsOf = (df: DG.DataFrame, arm: string, part: ArmPart | null, blocked = false): number[] | string => {
    const all = Array.from({length: df.rowCount}, (_, i) => i);
    if (groupCol && part) {
      return drugCol ? all.filter((i) => str(df, groupCol, i) === arm && drugKey(str(df, drugCol, i)) === part.key) :
        'the concentrations do not say which drug each is';
    }
    if (groupCol)
      return all.filter((i) => str(df, groupCol, i) === arm);
    const parts = part ? [part] : armParts(arm, meta);
    if (parts.length !== 1)
      return 'a combination is matched by its exact name only (add a group column to the concentrations)';
    const [p] = parts;
    if (!drugCol)
      return 'the concentrations name no drug or group';
    if (!p.dose || doseUnit === null || unitOf(p.dose.unit) !== doseUnit)
      return 'its dose is not given in the same unit as NCA’s';
    return all.filter((i) => drugKey(str(df, drugCol, i)) === p.key && (blocked ?
      str(df, 'VERDICT', i) === 'Blocked' : sameRoute(p.route, str(df, 'ROUTE', i)) &&
      Math.abs((num(df, 'DOSE', i) ?? NaN) - p.dose!.value) <= 1e-9 * Math.max(1, p.dose!.value)));
  };
  return (arm: string, part: ArmPart | null): ArmExposure => {
    const no = (why: string): ArmExposure =>
      ({why: `NCA could not compute exposure for ${arm}${part ? ` (${part.drug})` : ''}: ${why}`});
    const inSummary = rowsOf(summary, arm, part);
    if (typeof inSummary === 'string')
      return no(inSummary);
    if (inSummary.length > 1) {
      const sexes = new Set(inSummary.map((i) => str(summary, 'SEX', i)));
      return no(sexes.size > 1 ? 'NCA gives its males and females apart' : groupCol && !part && drugCol &&
        armParts(arm, meta).length > 1 ? 'a combination: each drug\'s exposure is given apart' :
        'more than one NCA group fits it (by period or study)');
    }
    if (inSummary.length === 0) {
      const inResults = rowsOf(results, arm, part) as number[];
      const blocked = groupCol ? inResults : rowsOf(results, arm, part, true) as number[];
      const reasons = blocked.map((i) => str(results, 'REASONS', i)).filter((r) => r);
      return no(reasons.length ? `${groupCol ? '' : 'some doses of this drug were not computed, '}` +
        `${plain(reasons)}` : inResults.length === 0 ? 'the concentrations file has none for it' : 'no reason given');
    }
    const i = inSummary[0];
    const auc = num(summary, 'MEAN_AUC_LAST', i);
    if (auc === null)
      return no('no AUClast in NCA’s results');
    const unitTag = (c: string) => summary.col(c)?.getTag('units') ?? '';
    const parts = part ? [part] : armParts(arm, meta);
    return {exposure: {auc, aucCv: num(summary, 'CV_PCT_AUC_LAST', i), aucUnit: unitTag('MEAN_AUC_LAST'),
      cmax: num(summary, 'MEAN_CMAX', i), cmaxCv: num(summary, 'CV_PCT_CMAX', i), cmaxUnit: unitTag('MEAN_CMAX'),
      aucInf: num(summary, 'MEAN_AUC_INF', i), aucInfUnit: unitTag('MEAN_AUC_INF'),
      nAucInf: num(summary, 'N_AUC_INF', i), nAucLast: num(summary, 'N_AUC_LAST', i),
      // (to float32's precision: the column reads 0.3 as 0.30000001192092896)
      dose: `${float32(num(summary, 'DOSE', i)) ?? '?'}${doseUnit ? ` ${units.dose}` : ''}`,
      route: str(summary, 'ROUTE', i),
      schedule: parts.length === 1 ? parts[0].schedule : null, profiles: num(summary, 'N', i) ?? 0,
      drug: drugCol ? str(summary, drugCol, i) : arm, pooled: (rowsOf(results, arm, part) as number[])
        .some((j) => str(results, 'SUBJECT', j) === '(composite)')}};
  };
}

/** "1,230 (CV 18%)": three significant figures; no CV for a single profile. */
export const meanCv = (mean: number, cv: number | null) =>
  `${Number(mean.toPrecision(3)).toLocaleString('en-US')}${cv === null ? '' : ` (CV ${Math.round(cv)}%)`}`;

/** A drug's PK as NCA summarised it: its dose levels and, by dose, the mean AUClast as a line of text. */
/** A drug's PK as NCA read it: its doses, their unit, the AUC text per dose (for a tooltip) and whether exposure rises
 *  in step with the dose (pk_proportional of Crux's function). */
export interface DrugPk {doses: number[], unit: string, auc: Map<number, string>, proportional: Proportional,
  /** the mean AUC per dose as a number, AUCinf when NCA gives it, else AUClast (`inf` says which) */
  perDose?: Map<number, {auc: number, unit: string, inf: boolean}>}

export type Proportional = 'yes' | 'no' | 'unknown';

/** Whether exposure rises in step with the dose, as fase2 judges it: AUC/dose within 25% for every pair of doses;
 *  unknown with fewer than 2 doses with an AUC. */
export function proportionalOf(aucByDose: Map<number, number>): Proportional {
  const perDose = [...aucByDose].filter(([d, a]) => d > 0 && Number.isFinite(a) && a > 0).map(([d, a]) => a / d);
  if (perDose.length < 2)
    return 'unknown';
  return perDose.every((x) => perDose.every((y) => Math.abs(x / y - 1) <= 0.25)) ? 'yes' : 'no';
}

/** The PK of the drug of key `key` (arm-name's drugKey) in NCA's results and summary: the summary's groups of that
 *  drug (its drug column, else a group named by the drug alone) given by `route` (an arm name's, null = any), their
 *  distinct doses, mean AUC (last) per dose and whether it rises in step with the dose; null when none. NCA leaves
 *  its blocked profiles out of the summary; those still to review are in it, and the AUC text says how many. */
export function drugPkOf(key: string, results: DG.DataFrame, summary: DG.DataFrame,
  groupCol: string | null, route: string | null = null): DrugPk | null {
  const keys: string[] = JSON.parse(summary.getTag('nca.keyColumns') || '[]');
  const units = JSON.parse(results.getTag('nca.units') || '{}');
  const drugCol = keys.find((k) => k !== groupCol && !/^(study|period|visit|occasion)/i.test(k)) ?? null;
  const dose = summary.col('DOSE');
  if (!dose || typeof units.dose !== 'string' || (!drugCol && !groupCol))
    return null;
  const auc = new Map<number, string>();
  const aucs = new Map<number, number>();
  const perDose = new Map<number, {auc: number, unit: string, inf: boolean}>();
  const routes = summary.col('ROUTE');
  const n = summary.col('N');
  const review = summary.col('N_REVIEW');
  for (let i = 0; i < summary.rowCount; i++) {
    const parts = drugCol ? [] : parseArmName(summary.col(groupCol!)!.getString(i));
    const theirs = drugCol ? drugKey(summary.col(drugCol)!.getString(i)) : parts.length === 1 ? parts[0].key : null;
    const d = dose.isNone(i) ? NaN : dose.get(i);
    if (theirs !== key || !Number.isFinite(d) || auc.has(d) || (routes && !sameRoute(route, routes.getString(i))))
      continue;
    const toReview = review && !review.isNone(i) ? review.get(i) : 0;
    const of = n && !n.isNone(i) ? ` of ${n.get(i)}` : '';
    const pending = toReview > 0 ? `; ${toReview}${of} profile${toReview === 1 ? '' : 's'} still to review in NCA` : '';
    const m = summary.col('MEAN_AUC_LAST');
    const unit = m?.getTag('units') ?? '';
    if (m && !m.isNone(i))
      aucs.set(d, m.get(i));
    const inf = summary.col('MEAN_AUC_INF');
    if (inf && !inf.isNone(i) && Number.isFinite(inf.get(i)) && inf.get(i) > 0)
      perDose.set(d, {auc: inf.get(i), unit: inf.getTag('units') ?? '', inf: true});
    else if (m && !m.isNone(i) && Number.isFinite(m.get(i)) && m.get(i) > 0)
      perDose.set(d, {auc: m.get(i), unit, inf: false});
    auc.set(d, m && !m.isNone(i) && Number.isFinite(m.get(i)) ?
      `${meanCv(m.get(i), summary.col('CV_PCT_AUC_LAST')?.isNone(i) === false ?
        summary.col('CV_PCT_AUC_LAST')!.get(i) : null)}${unit ? ` ${unit}` : ''}${pending}` : '');
  }
  return auc.size ? {doses: [...auc.keys()].sort((a, b) => a - b), unit: units.dose, auc,
    proportional: proportionalOf(aucs), perDose} : null;
}

// each table's NCA run, once until its data changes (forgetNca)
const runs = new WeakMap<DG.DataFrame, Promise<Awaited<ReturnType<typeof runNca>> | null>>();

/** Drops `table`'s NCA run: the next look reads it again (its data changed). */
export const forgetNca = (table: DG.DataFrame) => runs.delete(table);

/** The PK of the drug of key `key` by route `route` (an arm name's, null = any) in the concentration tables open in
 *  the workspace, as Analyze reads them (NCA's mapping and run); null without NCA or without such a table. */
export async function drugPkInWorkspace(key: string, route: string | null = null): Promise<DrugPk | null> {
  if (!ncaInstalled())
    return null;
  for (const t of grok.shell.tables) {
    if (!runs.has(t))
      runs.set(t, concentrationMapping(t).then((m) => m ? runNca(t, m) : null).catch(() => null));
    const run = await runs.get(t)!;
    const pk = run && drugPkOf(key, run.results, run.summary, run.groupCol, route);
    if (pk)
      return pk;
  }
  return null;
}

// ---- the mouse's exposure as a multiple of the human's (a description: nothing saved, nothing analyzed) ----

/** Doses a day of a schedule (schedule.ts, the one reader of schedules; here for those that read exposure). */
export {dosesPerDay};

// powers of ten against mg and L
const MASS: {[u: string]: number} = {g: 3, mg: 0, ug: -3, ng: -6, pg: -9};
const VOLUME: {[u: string]: number} = {l: 0, dl: -1, ml: -3, ul: -6};
const TIME: {[u: string]: number} = {h: 1, hr: 1, hrs: 1, min: 1 / 60, d: 24, day: 24};
const unitText = (u: string) => unitOf(u).replace(/\s+/g, '');

/** A concentration unit (mass per volume: ng/mL, µg/mL, mg/L) as its factor to mg/L; null for any other (a molar
 *  unit needs the molecular weight, and is refused). */
export function concentrationFactor(unit: string): number | null {
  const m = /^(g|mg|ug|ng|pg)\/(l|dl|ml|ul)$/.exec(unitText(unit));
  return m ? 10 ** (MASS[m[1]] - VOLUME[m[2]]) : null;
}

/** An AUC unit (a concentration times a time, however written: mg/L·h, h*ng/mL, ng·h/mL, µg·h/mL) as its factor to
 *  mg·h/L; null when it does not read so. */
export function aucFactor(unit: string): number | null {
  const u = unitText(unit);
  const m = /(^|[*·.×x])(hrs|hr|h|min|day|d)(?=$|[*·.×x/])/.exec(u);
  if (!m)
    return null;
  const rest = (u.slice(0, m.index) + u.slice(m.index + m[0].length)).replace(/^[*·.×x]+/, '');
  const c = concentrationFactor(rest);
  return c === null ? null : c * TIME[m[2]];
}

/** The patients' exposure of one drug, as typed under Advanced: AUC0–24 and (optionally) Cmax at
 *  steady state in one concentration unit (the AUC in it times hours), the unbound fractions if known, and whether it
 *  was measured (a label, a clinical study) or projected (a DMPK projection), and the patients' CV % of the AUC if
 *  known. */
export interface HumanExposure {
  drug: string, auc: number | null, cmax: number | null, unit: string, fuMouse: number | null,
  fuHuman: number | null, source: 'measured' | 'projected', cv?: number | null,
  /** filled from PMX Studio's reference table (files/human_exposure_v1.json): its row, and whether a value was
   *  changed since */
  reference?: {row: ReferenceRow, edited: boolean},
}

/** A drug of PMX Studio's reference table of human exposure (files/human_exposure_v1.json, built by the preclinical
 *  bundle's human_exposure_source/table/scripts/export_app_json.py from reviewed label, FDA review and EMA values):
 *  ng/mL, the AUC in ng·h/mL; `cycleDays`: the days its AUC covers (per infusion, per cycle), null when it is the
 *  daily AUC0–24 at steady state; `source`: document, version and section. */
export interface ReferenceRow {
  drug: string, synonyms: string[], regimen: string, auc: number | null, aucWindow: string,
  cycleDays: number | null, aucCv: number | null, cmax: number | null, cmaxCv: number | null, fu: number | null,
  matrix: 'plasma' | 'whole_blood', stat: '' | 'arith_mean' | 'geo_mean' | 'median' | 'model', tier: string,
  source: string,
  /** what was measured, e.g. "oxaliplatin (intact, in blood ultrafiltrate; not platinum)"; absent in rows kept by
   *  studies saved before it was exported */
  analyte?: string,
}

/** The analyte when it is more than the parent drug itself ("doxorubicin (parent)"): the mouse value must match it. */
export const notableAnalyte = (r: ReferenceRow): string | null =>
  r.analyte && !/^[\w-]+ \(parent\)$/i.test(r.analyte.trim()) ? r.analyte.trim() : null;

/** The row of drug key `key` (arm-name's drugKey) in `rows`, by name or synonym; undefined: not in the table. */
export const referenceRow = (key: string, rows: ReferenceRow[]): ReferenceRow | undefined =>
  rows.find((r) => [r.drug, ...r.synonyms].some((n) => drugKey(n) === key));

/** A drug's human exposure as the reference table gives it: measured, in ng/mL. */
export const fromReference = (drug: string, r: ReferenceRow): HumanExposure => ({drug, auc: r.auc, cmax: r.cmax,
  unit: 'ng/mL', fuMouse: null, fuHuman: r.fu, source: 'measured', cv: cvOk(r.aucCv) ? r.aucCv : null,
  reference: {row: r, edited: false}});

/** Whether `h` still holds its reference row's values (the mouse's unbound fraction is not in the table). */
export function referenceKept(h: HumanExposure, r: ReferenceRow): boolean {
  const f = fromReference(h.drug, r);
  return h.auc === f.auc && h.cmax === f.cmax && h.unit === f.unit && h.fuHuman === f.fuHuman &&
    h.source === f.source && (h.cv ?? null) === f.cv;
}

/** The days the AUC of `h` covers: its reference row's cycle while the AUC is the table's; null (a daily AUC0–24)
 *  once typed or edited. */
export const cycleOf = (h: HumanExposure): number | null =>
  h.reference && tableAuc(h) ? h.reference.row.cycleDays : null;

/** The patients' average daily AUC (mg·h/L): an AUC per administration ÷ the days between them (a daily AUC0–24 as
 *  is), unbound (× fu in patients) when asked; null: no AUC, or a unit it cannot convert. */
export function patientsDaily(h: HumanExposure, unbound = false): number | null {
  const f = concentrationFactor(h.unit);
  return f === null || h.auc === null ? null : h.auc * f / (cycleOf(h) ?? 1) * (unbound ? h.fuHuman ?? 1 : 1);
}

/** Whether the AUC of `h` is its reference row's (ng·h/mL), in whatever unit it is now written. */
function tableAuc(h: HumanExposure): boolean {
  const f = concentrationFactor(h.unit);
  const r = h.reference?.row.auc;
  return f !== null && h.auc !== null && r != null && Math.abs(h.auc * f - r * 1e-3) <= 1e-6 * r * 1e-3;
}

/** The table's per-administration AUC as a daily one. */
const cycleText = (days: number) => `AUC per administration ÷ ${days} days between administrations (on average)`;

export const STAT_TEXT: Record<string, string> = {arith_mean: 'mean', geo_mean: 'geometric mean', median: 'median',
  model: 'model estimate'};

/** Where a drug's patients' exposure comes from, in the one wording every view and report uses: the reference
 *  table's document and the regimen its values are for, else projected, else the user's value. */
export const sourceText = (h: HumanExposure): string => h.reference ?
  `${h.reference.row.source}${h.reference.row.regimen ? `, ${h.reference.row.regimen}` : ''}` +
    `${h.reference.edited ? ' (edited)' : ''}` : h.source === 'projected' ? 'Projected' : 'Your value';

/** The line under the dialog's table for a drug filled from the reference table. */
export const referenceLine = (h: HumanExposure): string => !h.reference ? '' :
  `${h.drug}: ` + [sourceText(h), h.reference.row.aucWindow, STAT_TEXT[h.reference.row.stat],
    ...h.reference.row.cycleDays ? [`table value: ${cycleText(h.reference.row.cycleDays)}` +
      `${cycleOf(h) === null ? '; your AUC is used as the daily AUC0–24' : ''}`] : [],
    ...!tableAuc(h) && h.auc === h.reference.row.auc ? ['the table\'s values are in ng/mL: the unit was ' +
      'changed, not the numbers'] : []].filter((t) => t).join(' · ');

/** What the reference row's values are, when it matters to the multiple: whole blood, a median, a model estimate. */
function referenceCaveats(r: ReferenceRow): string[] {
  const analyte = notableAnalyte(r);
  return [...analyte ? [`measured as ${analyte}: compare with the same in the mouse`] : [],
    ...r.matrix === 'whole_blood' ? ['whole blood, not plasma: compare with mouse whole-blood exposure'] : [],
    ...r.stat === 'median' ? ['a median across patients, not a mean'] : [],
    ...r.stat === 'model' ? ['a population PK model estimate, not a measured mean'] : []];
}

/** The units the human exposure can be typed in: a concentration, the AUC in it times hours. */
export const HUMAN_UNITS = ['ng/mL', 'µg/mL', 'mg/L'];

/** A group's exposure against the human's: the ratio and its tooltip, or why there is none. */
/** `weekly`: the study's no-dosing-at-weekends rule lowered the doses a day, so the multiple is averaged over the
 *  week; its line ("Mouse daily AUC: … × 0.714 doses a day (QD, no dosing at weekends) = …"), for the report. */
export type HumanMultiple = {ratio: number, tip: string, range?: 'above' | 'below', why?: undefined,
  weekly?: string} |
  {ratio?: undefined, tip?: undefined, range?: undefined, why: string, weekly?: undefined};

/** A CV % of the human AUC that is used: from 5 to under 300 (below 5, likely typed as a fraction). */
export const cvOk = (cv: number | null | undefined): cv is number => cv != null && cv >= 5 && cv < 300;

/** Whether a multiple of the typed human exposure (a label's arithmetic mean) lies outside the patients' spread:
 *  with the AUC's CV %, outside the 5th–95th percentiles of a log-normal with that CV centred on its geometric mean
 *  (mean/√(1+CV²)), exp(±1.645·√ln(1+CV²)); without it, beyond 3-fold either way, by convention. `notes` qualify
 *  the text. null: within. */
/** The patients' range as multiples of their (mean) exposure: the 5th to 95th percentile of a log-normal around the
 *  geometric mean when the CV is known, else the 3-fold convention. */
export function patientsRange(cv: number | null | undefined): {lo: number, hi: number} {
  if (!cvOk(cv))
    return {lo: 1 / 3, hi: 3};
  const s2 = Math.log(1 + (cv / 100) ** 2);
  const m = Math.exp(-s2 / 2); // the geometric mean, as a multiple of the mean
  const edge = Math.exp(1.645 * Math.sqrt(s2));
  return {lo: m / edge, hi: m * edge};
}

export function exposureRange(ratio: number, cv: number | null | undefined, notes: string[] = []):
{dir: 'above' | 'below', text: string} | null {
  const known = cvOk(cv);
  const {lo, hi} = patientsRange(cv);
  const dir = ratio > hi ? 'above' : ratio < lo ? 'below' : null;
  if (!dir)
    return null;
  const basis = known ? `${dir === 'above' ? '95th' : '5th'} percentile, CV ${cv}%` :
    'a convention, not patients\' spread';
  const how = `(${[basis, ...notes].join('; ')})`;
  return {dir, text: known ? `${dir} the range seen in patients ${how}` :
    `more than 3-fold ${dir} the patients' exposure ${how}`};
}

const sig = (x: number) => Number(x.toPrecision(3)).toLocaleString('en-US');
export const fuOk = (x: number | null): x is number => x !== null && Number.isFinite(x) && x > 0 && x <= 1;

/** The AUC of a group's NCA exposure that the multiple uses: AUCinf when NCA gives one, else AUClast. */
export const aucUsed = (x: Exposure): {auc: number, unit: string, inf: boolean} =>
  x.aucInf !== null && x.aucInf > 0 ? {auc: x.aucInf, unit: x.aucInfUnit || x.aucUnit, inf: true} :
    {auc: x.auc, unit: x.aucUnit, inf: false};

/** The mouse AUC per dose of a group, its unit and which AUC it is; else why there is none. */
function mouseAuc(arm: string, p: ArmPart, x: Exposure | undefined, pk: DrugPk | null):
{auc: number, unit: string, which: string, scaled: string} | string {
  if (x) {
    const {auc, unit, inf} = aucUsed(x);
    if (!(auc > 0))
      return `${arm}: NCA's AUC is not above 0`;
    // AUCinf needs a terminal phase: fewer animals may give it than AUClast
    const fewer = inf && x.nAucInf != null && x.nAucLast != null && x.nAucInf < x.nAucLast ?
      ` (from ${x.nAucInf} of the ${x.nAucLast} animals with an AUClast)` : '';
    return {auc, unit, scaled: '',
      which: inf ? `AUCinf${fewer}` : 'AUClast (NCA gives no AUCinf: it reads low when sampling stops early)'};
  }
  const levels = pk?.perDose;
  if (!pk || !levels?.size || !p.dose || unitText(p.dose.unit) !== unitText(pk.unit))
    return `${arm}: NCA has no exposure at its dose`;
  const dose = p.dose.value;
  const whichOf = (inf: boolean) => inf ? 'AUCinf' : 'AUClast (NCA gives no AUCinf)';
  const all = [...levels.entries()].filter(([d, v]) => d > 0 && v.auc > 0);
  const at = all.find(([d]) => Math.abs(d - dose) <= 1e-9 * Math.max(1, dose))?.[1];
  if (at)
    return {auc: at.auc, unit: at.unit, which: whichOf(at.inf), scaled: ''};
  const doses = all.map(([d]) => d);
  if (!all.length || dose < Math.min(...doses) || dose > Math.max(...doses)) {
    return `${arm}: ${dose} ${p.dose.unit} is outside NCA's doses${doses.length ? ` (${Math.min(...doses)} to ` +
      `${Math.max(...doses)})` : ''}, and exposure is not extrapolated`;
  }
  if (new Set(all.map(([, v]) => `${v.inf}|${v.unit}`)).size !== 1)
    return `${arm}: NCA's doses give their AUC in different ways, so it is not scaled to ${dose}`;
  // proportionality judged on the AUC used
  const proportional = proportionalOf(new Map(all.map(([d, v]) => [d, v.auc])));
  if (proportional !== 'yes') {
    return `${arm}: NCA has no exposure at ${dose} ${p.dose.unit}, and the drug's exposure ` +
      `${proportional === 'no' ? 'does not rise in step with the dose' : 'is known at one dose only'}, so it is ` +
      'not scaled';
  }
  return {auc: all.reduce((t, [d, v]) => t + v.auc / d, 0) / all.length * dose, unit: all[0][1].unit,
    which: whichOf(all[0][1].inf), scaled: `, scaled to ${dose} ${p.dose.unit} from NCA's doses (the ` +
      'exposure rises in step with the dose)'};
}

/** How many times the human exposure the mouse's daily exposure is in group `arm`: its AUC per dose (AUCinf, which
 *  after one dose equals the AUC over a dosing interval at steady state under linear PK; else AUClast) times the
 *  doses a day of its schedule, against the human AUC0–24,ss; on unbound exposure when both fractions are given.
 *  `x`: NCA's exposure matched to the group; else `pk`, the drug's PK, scaled to the group's dose only when exposure
 *  rises in step with the dose. null: nothing to say (no single drug, no human value for it). Never a ratio across
 *  units it cannot convert. */
export function humanMultiple(arm: string, x: Exposure | undefined, pk: DrugPk | null,
  humans: Map<string, HumanExposure>, meta: StudyMeta | null = null): HumanMultiple | null {
  const parts = armParts(arm, meta);
  if (parts.length === 0)
    return null;
  if (parts.length > 1) {
    return parts.some((q) => humans.get(q.key)?.auc) ?
      {why: `${arm} is a combination: its drugs' exposures are not added up into one multiple`} : null;
  }
  return partMultiple(arm, parts[0], x, pk, humans, meta);
}

/** Each drug of combination `arm` against its patients' exposure, as humanMultiple for one drug (`drugs`: NCA's
 *  exposure and PK of each, NcaByArm.drugs); null for a drug without the patients' exposure. */
export function drugMultiples(arm: string, drugs: DrugExposure[], humans: Map<string, HumanExposure>,
  meta: StudyMeta | null = null): {drug: string, key: string, m: HumanMultiple | null}[] {
  return armParts(arm, meta).map((p) => {
    const d = drugs.find((x) => x.key === p.key);
    return {drug: p.drug, key: p.key, m: partMultiple(`${arm}, ${p.drug}`, p, d?.x.exposure, d?.pk ?? null, humans,
      meta)};
  });
}

/** One drug's multiple (humanMultiple); `arm` names it in what is said. */
function partMultiple(arm: string, p: ArmPart, x: Exposure | undefined, pk: DrugPk | null,
  humans: Map<string, HumanExposure>, meta: StudyMeta | null): HumanMultiple | null {
  const h = humans.get(p.key);
  if (!h || h.auc === null || !(h.auc > 0))
    return null;
  // an AUC per infusion or per cycle, as an average daily exposure over the cycle, like the mouse's
  const cycle = cycleOf(h);
  const schedule = p.schedule ?? x?.schedule ?? null;
  // the doses over the schedule's span ÷ its days, under the study's rules (no dosing at weekends)
  const perDay = dosesPerDay(schedule, meta);
  const weekendsOff = !!(meta?.weekendsOff && meta.start) && perDay !== null && perDay !== dosesPerDay(schedule);
  if (perDay === null) {
    return {why: readSchedule(schedule, p.route ?? '')?.kind === 'continuous' ?
      `${arm}: no multiple is given for continuous exposure (pump, drinking water, diet): it has no doses to count` :
      `${arm}: its schedule is not read from its details or its name (write it as QD, BID, Q3D or QW), so its ` +
      'daily exposure is unknown'};
  }
  if (p.dose && !(p.dose.value > 0))
    return {why: `${arm}: its dose is not above 0`};
  const m = mouseAuc(arm, p, x, pk);
  if (typeof m === 'string')
    return {why: m};
  const mouseF = aucFactor(m.unit);
  const humanF = concentrationFactor(h.unit);
  if (mouseF === null)
    return {why: `${arm}: NCA's AUC unit (${m.unit || 'none'}) cannot be converted to compare it with the patients'`};
  if (humanF === null)
    return {why: `${arm}: the patients' exposure's unit (${h.unit}) cannot be converted`};
  const mouseDaily = m.auc * mouseF * perDay;
  const humanDaily = patientsDaily(h)!;
  const unbound = fuOk(h.fuMouse) && fuOk(h.fuHuman);
  const ratio = (unbound ? mouseDaily * h.fuMouse! : mouseDaily) / patientsDaily(h, unbound)!;
  const lines: string[] = [`Mouse daily AUC: ${m.which} per dose${m.scaled} ` +
    `× ${sig(perDay)} dose${perDay === 1 ? '' : 's'} ` +
    `a day (${schedule}${weekendsOff ? ', no dosing at weekends' : ''}) = ${sig(mouseDaily)} mg/L·h` +
    `${perDay < 1 ? ', averaged over the days between doses' : ''}.`,
  cycle ? `Patients' average daily exposure over the cycle, ${h.source}: ${cycleText(cycle)} ` +
    `(${h.reference!.row.aucWindow}) = ${sig(humanDaily)} mg/L·h.` :
    // (as typed too, in its unit: the patients' dialog takes ng·h/mL where NCA gives mg/L·h)
    `Patients' AUC0–24 at steady state, ${h.source}: ${humanF === 1 ? '' : `${sig(h.auc!)} ${h.unit}·h = `}` +
    `${sig(humanDaily)} mg/L·h.`,
  // an intermittent regimen: also one administration against one dose, whatever the spacing of each
  ...cycle ? [`Per administration: the mouse's AUC per dose against the patients' AUC per administration is ` +
    `${sig((unbound ? m.auc * mouseF * h.fuMouse! : m.auc * mouseF) / (patientsDaily(h, unbound)! * cycle))}× ` +
    '(the multiple above compares daily averages).'] : [],
  unbound ? `On unbound exposure: unbound fraction ${h.fuMouse} in mouse and ${h.fuHuman} in patients' plasma.` :
    'On total exposure (give both unbound fractions to compare unbound exposure).'];
  const notes = [...h.source === 'projected' ? ['projected'] : [], ...unbound ? [] : ['total, not free']];
  const range = exposureRange(ratio, h.cv, notes);
  if (range)
    lines.unshift(`${cycle ? 'Average daily exposure over the cycle' : 'Daily AUC'} ${range.text}.`);
  const caveats = h.reference ? referenceCaveats(h.reference.row) : [];
  if (caveats.length)
    lines.splice(range ? 1 : 0, 0, `Patients' values: ${caveats.join('; ')}.`);
  if (x && x.cmax !== null && h.cmax !== null && h.cmax > 0) {
    const cm = concentrationFactor(x.cmaxUnit);
    if (cm !== null) {
      const c = x.cmax * cm / (h.cmax * humanF) * (unbound ? h.fuMouse! / h.fuHuman! : 1);
      lines.push(`Single-dose Cmax: ${sig(c)}× the patients' Cmax at steady state (the mouse's after one dose; ` +
        `higher ` +
        'at steady state if the drug accumulates).');
      const cr = exposureRange(c, null, notes);
      if (cr)
        lines.push(`Cmax ${cr.text}.`);
    }
  }
  const daily = lines.find((l) => l.startsWith('Mouse daily AUC'))!;
  lines.push('A description: nothing is fitted or tested.');
  return {ratio, tip: lines.join(' '), range: range?.dir, ...weekendsOff ? {weekly: daily} : {}};
}

/** "3.2×": two significant figures under 10, whole numbers from 10. */
export const multipleText = (r: number) =>
  `${r >= 10 ? Math.round(r).toLocaleString('en-US') : Number(r.toPrecision(2))}×`;

// ---- exposure against effect: a curve within the doses tested (a description: nothing saved, nothing analyzed) ----
// The composition validated in the preclinical bundle's analysis/out_translation/exposure_target (REPORT.md, M-ER):
// each group's kill rate κ from its TGI in closed form against the control's Gompertz curve (Crux's
// stats.tumor_effect_from_study), the line κ = k2·dose through 0 fitted on the groups, the TGI at k2·dose. Between
// the doses tested 4 doses left out in 5 were within 10 points of TGI (12 doses, 9 families; some missed by 20-30);
// beyond them it adds nothing, so it is drawn only between them, with a band of ±10 points.

/** A group on the exposure chart: its name and NCA's mean AUClast (null: none). */
export interface CurveGroup {arm: string, auc: number | null}
/** A measurement of the study the curve reads (the analysis's rows, excluded ones left out). */
export interface CurveRow {animal: string, group: string, day: number, volume: number}
/** The chart's X: NCA's AUClast, or the dose when every group names one (in one unit) and exposure does not rise in
 *  step with it; `x` places a dose on it (between two doses tested, the AUC between theirs). */
export interface ExposureAxis {axis: 'auc' | 'dose', unit: string, of: (g: CurveGroup) => number,
  x: (dose: number) => number}
/** The expected TGI (1 − T/C of the medians, %) and its band, from the lowest dose drawn to the highest. */
export interface ExposureCurve {curve: {x: number, tgi: number, lo: number, hi: number}[], note?: string}

/** The band around the curve, in points of TGI: about 4 doses left out in 5 were within it in the validation. */
export const CURVE_BAND = 10;
/** Medians below it are raised to it (mm³), as Crux's function does. */
const FLOOR = 3.4;
const IMMUNE: Therapy[] = ['immunotherapy', 'immune_agonist', 'cell_therapy', 'vaccine_peptide'];
const allSame = (xs: (string | null)[]) =>
  new Set(xs.map((x) => (x ?? '').replace(/\s+/g, '').toLowerCase())).size === 1;
// a group's dose: its details first, then its name
const doseOf = (arm: string, meta: StudyMeta | null = null) => {
  const p = armParts(arm, meta);
  return p.length === 1 ? p[0].dose : null;
};

/** Why no curve is drawn through these groups, read from their names (null: nothing seen here; Crux checks the
 *  animals). Only what was validated: one drug alone, not acting through the immune system, at 3 doses or more of one
 *  unit, schedule and route. */
export function curveSilence(arms: string[], meta: StudyMeta | null = null): string | null {
  const parts = arms.map((a) => armParts(a, meta));
  if (parts.some((p) => p.length > 1))
    return 'a group is a combination';
  if (parts.some((p) => p.length === 0 || !p[0].dose))
    return 'a group\'s name gives no dose';
  const ps = parts.map((p) => p[0]);
  if (!allSame(ps.map((p) => p.key)))
    return 'the groups are not all one drug';
  if (IMMUNE.includes(therapyOf(arms[0]).therapy))
    return 'the drug acts through the immune system';
  if (!allSame(ps.map((p) => unitOf(p.dose!.unit))))
    return 'the doses are not all in one unit';
  if (!allSame(ps.map((p) => p.schedule)) || !allSame(ps.map((p) => p.route ? ROUTE_OF[p.route] ?? p.route : null)))
    return 'the groups are not all dosed on the same schedule and route';
  if (new Set(ps.map((p) => p.dose!.value)).size < 3)
    return 'fewer than 3 doses';
  return null;
}

/** The chart's axis for these groups (see ExposureAxis). */
export function exposureAxis(groups: CurveGroup[], meta: StudyMeta | null = null): ExposureAxis {
  const doses = groups.map((g) => doseOf(g.arm, meta));
  const dosed = doses.every((d) => d && d.value > 0) && allSame(doses.map((d) => unitOf(d!.unit)));
  const aucs = new Map<number, number[]>();
  groups.forEach((g, i) => dosed && aucs.set(doses[i]!.value, [...aucs.get(doses[i]!.value) ?? [], g.auc ?? NaN]));
  const byDose = [...aucs].map(([d, a]) => [d, a.reduce((t, x) => t + x, 0) / a.length] as [number, number])
    .sort((a, b) => a[0] - b[0]);
  if (dosed && (!byDose.every(([, a]) => a > 0) || proportionalOf(new Map(byDose)) !== 'yes'))
    return {axis: 'dose', unit: doses[0]!.unit, of: (g) => doseOf(g.arm, meta)!.value, x: (d) => d};
  return {axis: 'auc', unit: '', of: (g) => g.auc ?? NaN, x: (d) => {
    const found = byDose.findIndex(([e]) => e >= d);
    const j = found < 0 ? byDose.length - 1 : Math.max(1, found);
    const [[d0, a0], [d1, a1]] = [byDose[j - 1], byDose[j]];
    return a0 + (a1 - a0) * (d - d0) / (d1 - d0);
  }};
}

/** The rows from the first dose to `lastDay` (null: any). */
const inWindow = (rows: CurveRow[], lastDay: number | null) => rows.filter((x) => x.day >= FIRST_DOSE_DAY &&
  (lastDay === null || x.day <= lastDay) && Number.isFinite(x.volume));

/** The day a group's 1 − T/C of the medians is read on: its last up to `lastDay` with 3 animals in it and in the
 *  control (null: no such day). */
export const medianDay = (rows: CurveRow[], arm: string, control: string, lastDay: number | null): number | null =>
  readDayOf({name: '', line: null, meta: null, rows: inWindow(rows, lastDay)}, arm, control);

/** A group's 1 − T/C of the median volumes (%) on its last day up to `lastDay` with 3 animals in it and in the
 *  control, as the curve reads it (null: no such day). */
export function medianTgi(rows: CurveRow[], arm: string, control: string, lastDay: number | null): number | null {
  const study = {name: '', line: null, meta: null, rows: inWindow(rows, lastDay)};
  const day = readDayOf(study, arm, control);
  const at = (g: string) => Math.max(FLOOR, median(study.rows.filter((x) => x.group === g && x.day === day)
    .map((x) => x.volume)));
  return day === null ? null : 100 * (1 - at(arm) / at(control));
}

/** Crux's reasons, in this chart's words. */
const CURVE_REASONS: {[code: string]: string} = {
  too_few_animals: 'a group or the control has fewer than 3 animals measured on any one day',
  control_fit_failed: 'the control group\'s growth could not be described',
  ratio_out_of_range: 'the tested doses are too far apart to draw a curve between them',
  bootstrap_failed: 'the effect could not be read from the animals',
};
const curveReason = (code: string) => CURVE_REASONS[code.trim()] ?? 'the groups give no basis for a curve';

/** Points of the curve from the lowest dose to the highest. */
const CURVE_POINTS = 13;

/** The curve of the expected TGI through `groups` and their `control`, read on the last day up to `lastDay` (the
 *  evaluation day; null: any) with 3 animals in every group and the control, with a band of ±CURVE_BAND points; or
 *  why there is none. Doses from the first that stopped its tumors almost completely (TGI over 95%) are left out when
 *  3 doses remain below it, and the note says so. */
export async function exposureCurve(groups: CurveGroup[], rows: CurveRow[], control: string,
  lastDay: number | null, meta: StudyMeta | null = null): Promise<ExposureCurve | {why: string}> {
  const silent = curveSilence(groups.map((g) => g.arm), meta);
  if (silent)
    return {why: silent};
  const arms = new Set(groups.map((g) => g.arm));
  const r = inWindow(rows, lastDay).filter((x) => x.group === control || arms.has(x.group));
  const table = {animal: r.map((x) => `${x.group}|${x.animal}`), group: r.map((x) => x.group),
    day: Float64Array.from(r, (x) => x.day), volume: Float64Array.from(r, (x) => x.volume)};
  // κ_B = κ_A·dose_new/dose_prev: with dose_prev = κ_A and dose_new = κ_B any κ is asked of any group (the ratio is
  // all Crux reads of the doses); exposure "in step" lifts its window of doses with PK, which the curve stays inside.
  // Only the point is used: a few bootstrap replicates (Crux gives no point when every replicate fails)
  const ask = (anchor: string, prev: number, next: number, more: object = {}) => crux('stats.tumor_effect_from_study',
    table, {anchor, control, first_dose_day: FIRST_DOSE_DAY, dose_prev: prev, dose_new: next, pk_proportional: 'yes',
      boot: 20, seed: SEED, ...more}).then((o) => cruxRows(o.main)[0] ?? {});
  const unit = doseOf(groups[0].arm, meta)!.unit;
  try {
    // each group's own κ (the point at its own dose) and its TGI as Crux reads it
    const own = (await Promise.all(groups.map((g) => ask(g.arm, 1, 1)))).map((o, i) => ({arm: groups[i].arm,
      d: doseOf(groups[i].arm, meta)!.value, tgi: finite(o.tgi_anchor), reason: String(o.reason ?? ''),
      k: finite(o.kill_rate_point) ?? 0, day: finite(o.read_day)}));
    // from the first dose that stopped the tumors almost completely, no κ can be read: drawn below it, or not at all
    const cut = Math.min(...own.filter((o) => o.tgi !== null && o.tgi > 95).map((o) => o.d));
    const fit = own.filter((o) => o.d < cut);
    if (Number.isFinite(cut) && new Set(fit.map((o) => o.d)).size < 3) {
      return {why: `from ${cut} ${unit} the tumors stopped almost completely, which leaves fewer than 3 doses to ` +
        'draw a curve between'};
    }
    const bad = fit.find((o) => o.tgi === null || (o.reason && o.reason !== 'anchor_effect_out_of_range'));
    if (bad)
      return {why: curveReason(bad.reason)};
    const day = fit[0].day;
    if (day === null || fit.some((o) => o.day !== day))
      return {why: 'the groups cannot be read on one day with 3 animals in each'};
    // ponytail: a TGI under 5% (Crux will not scale from it) is taken as no kill; the validated fit kept its small κ
    fit.forEach((o) => o.tgi! < 5 && (o.k = 0));
    const k2 = fit.reduce((t, x) => t + x.d * x.k, 0) / fit.reduce((t, x) => t + x.d * x.d, 0);
    const anchors = fit.filter((x) => x.k > 0);
    if (!(k2 > 0) || !anchors.length)
      return {why: 'no group slowed its tumors enough to draw a curve through'};
    // the curve starts from the groups' median start (each group's own otherwise)
    const starts = fit.map((o) => {
      const first = Math.min(...r.filter((x) => x.group === o.arm).map((x) => x.day));
      return median(r.filter((x) => x.group === o.arm && x.day === first).map((x) => x.volume));
    });
    const [lo, hi] = [Math.min(...fit.map((x) => x.d)), Math.max(...fit.map((x) => x.d))];
    const grid = Array.from({length: CURVE_POINTS}, (_, i) => i === CURVE_POINTS - 1 ? hi :
      lo + (hi - lo) * i / (CURVE_POINTS - 1));
    const start = median(starts) > 0 ? {start_volume_new: median(starts)} : {};
    const curve = await Promise.all(grid.map((d) => {
      const k = k2 * d;
      const a = anchors.reduce((b, x) => Math.abs(Math.log(k / x.k)) < Math.abs(Math.log(k / b.k)) ? x : b);
      return ask(a.arm, a.k, k, start);
    }));
    const failed = curve.find((o) => o.reason || finite(o.tgi_point) === null);
    if (failed)
      return {why: curveReason(String(failed.reason ?? ''))};
    if (curve.some((o) => finite(o.read_day) !== day))
      return {why: 'the groups cannot be read on one day with 3 animals in each'};
    const axis = exposureAxis(groups, meta);
    return {curve: curve.map((o, i) => {
      const tgi = finite(o.tgi_point)!;
      return {x: axis.x(grid[i]), tgi, lo: Math.max(0, tgi - CURVE_BAND), hi: Math.min(100, tgi + CURVE_BAND)};
    }), ...Number.isFinite(cut) ? {note: `Drawn up to ${hi} ${unit}: from ${cut} ${unit} the tumors stopped almost ` +
      'completely.'} : {}};
  } catch (e) {
    if (/unknown function/.test(String((e as Error)?.message ?? e)))
      return {why: 'this version of PMX Studio cannot draw it yet'};
    console.error('PMX Studio: the exposure curve failed', e);
    return {why: 'it could not be computed'};
  }
}

// ---- what the study keeps of exposure, a sentence on it, and what would add to it ----

/** One drug of a combination group: NCA's exposure for it in that group, else its PK (drugPkOf) to scale from. */
export interface DrugExposure {key: string, drug: string, x: ArmExposure, pk: DrugPk | null}

/** NCA's exposure for each group as the groups table reads it: matched by exposureByArm and, for a group with none
 *  at its dose, its drug's PK (drugPkOf), from which humanMultiple may scale; a combination's drug by drug in
 *  `drugs` (then not in `arms`, unless NCA has nothing for any of its drugs: why); `table`: the concentrations'
 *  name. */
export interface NcaByArm {table: string, arms: Map<string, ArmExposure>, pk: Map<string, DrugPk>,
  drugs?: Map<string, DrugExposure[]>}

export function ncaByArm(arms: string[], run: {results: DG.DataFrame, summary: DG.DataFrame, groupCol: string | null},
  table: string, meta: StudyMeta | null = null): NcaByArm {
  const exposure = exposureByArm(arms, run.results, run.summary, run.groupCol, meta);
  const pk = new Map<string, DrugPk>();
  const drugs = new Map<string, DrugExposure[]>();
  for (const arm of [...exposure.keys()]) {
    const parts = armParts(arm, meta);
    if (parts.length > 1) {
      const each = exposureByDrug(arm, run.results, run.summary, run.groupCol, meta).map(({part, x}) =>
        ({key: part.key, drug: part.drug, x, pk: x.exposure ? null :
          drugPkOf(part.key, run.results, run.summary, run.groupCol, part.route)}));
      // drug by drug when NCA gives a drug its own exposure in the group, or (the group none of its own) its PK
      if (each.some((d) => d.x.exposure) || (!exposure.get(arm)?.exposure && each.some((d) => d.pk))) {
        drugs.set(arm, each);
        exposure.delete(arm);
      }
      continue;
    }
    const p = !exposure.get(arm)?.exposure && parts.length === 1 ?
      drugPkOf(parts[0].key, run.results, run.summary, run.groupCol, parts[0].route) : null;
    if (p)
      pk.set(arm, p);
  }
  return {table, arms: exposure, pk, ...drugs.size ? {drugs} : {}};
}

type SavedPk = Omit<DrugPk, 'auc' | 'perDose'> & {auc: [number, string][],
  perDose?: [number, {auc: number, unit: string, inf: boolean}][]};
const savedPk = (p: DrugPk): SavedPk => ({...p, auc: [...p.auc], perDose: p.perDose ? [...p.perDose] : undefined});
const pkOfSaved = (p: any): DrugPk | null => p && Array.isArray(p.doses) && Array.isArray(p.auc) &&
  typeof p.unit === 'string' ? {...p, auc: new Map(p.auc), perDose: Array.isArray(p.perDose) ? new Map(p.perDose) :
    undefined} : null;

/** Exposure as the study keeps it (study-file.ts's EXPOSURE_TAG and SavedStudy.exposure, JSON): the patients'
 *  exposure by drug key, as typed or filled from the reference table, and NCA's result for each group with the name
 *  of the concentrations it read (never the concentrations). */
export interface SavedExposure {
  humans?: {[key: string]: HumanExposure},
  /** the analysis's effect, kept for Compare: its evaluation day and each treated group's TGI (%); `ctd` the study's
   *  rows of the regulatory table, for the Overview's summary of a compound */
  effect?: {day: number | null, tgi: {[arm: string]: number}, ctd?: CtdStudy},
  nca?: {table: string, arms: {[arm: string]: ArmExposure}, pk: {[arm: string]: SavedPk},
    /** a combination group's drugs (NcaByArm.drugs) */
    drugs?: {[arm: string]: (Omit<DrugExposure, 'pk'> & {pk: SavedPk | null})[]}},
}

/** What the study keeps; undefined when there is nothing to keep. */
export function exposureJson(humans: Map<string, HumanExposure>, nca: NcaByArm | null): SavedExposure | undefined {
  const kept = [...humans].filter(([, h]) => h.auc !== null);
  if (!kept.length && !nca)
    return undefined;
  return {...kept.length ? {humans: Object.fromEntries(kept)} : {},
    ...nca ? {nca: {table: nca.table, arms: Object.fromEntries(nca.arms), pk: Object.fromEntries([...nca.pk]
      .map(([arm, p]) => [arm, savedPk(p)])), ...nca.drugs?.size ? {drugs: Object.fromEntries([...nca.drugs]
      .map(([arm, ds]) => [arm, ds.map((d) => ({...d, pk: d.pk && savedPk(d.pk)}))]))} : {}}} : {}};
}

/** What a study kept, read from a tag or file someone may have edited: what does not read is left out (an older
 *  study: nothing). */
export function parseExposure(x: unknown): {humans: Map<string, HumanExposure>, nca: NcaByArm | null} {
  const s = (x && typeof x === 'object' ? x : {}) as SavedExposure;
  const isNum = (v: unknown) => v === null || (typeof v === 'number' && Number.isFinite(v));
  const humans = new Map(Object.entries(s.humans && typeof s.humans === 'object' ? s.humans : {})
    .filter(([, h]) => h && typeof h.drug === 'string' && typeof h.unit === 'string' && isNum(h.auc) &&
      isNum(h.cmax) && isNum(h.fuMouse) && isNum(h.fuHuman) && (h.source === 'measured' || h.source === 'projected'))
    // a reference without its row says nothing: the values stay, as typed
    .map(([k, h]): [string, HumanExposure] => {
      if (!h.reference || (h.reference.row && typeof h.reference.row === 'object'))
        return [k, h];
      const {reference: _, ...rest} = h;
      return [k, rest];
    }));
  const n = s.nca;
  if (!n || typeof n !== 'object' || typeof n.table !== 'string' || !n.arms || typeof n.arms !== 'object')
    return {humans, nca: null};
  // a group's exposure only while its name still says the dose NCA read (a dose swapped between names: left out);
  // to float32's precision: NCA's dose column reads 0.3 as 0.30000001192092896
  const same = (named: {value: number, unit: string}, read: {value: number, unit: string}) =>
    Math.abs(named.value - read.value) <= 1e-6 * Math.max(Math.abs(named.value), Math.abs(read.value)) &&
    unitText(named.unit) === unitText(read.unit);
  const sameDose = (arm: string, dose: string) => {
    const [named, read] = [parseArmName(arm), parseArmName(`x ${dose}`)].map((p) => p.length === 1 ? p[0].dose : null);
    return !named || !read || same(named, read);
  };
  const readsAs = (a: ArmExposure | undefined, dose: (d: string) => boolean) => !!a && (typeof a.why === 'string' ||
    (!!a.exposure && typeof a.exposure.auc === 'number' && typeof a.exposure.aucUnit === 'string' &&
      typeof a.exposure.dose === 'string' && typeof a.exposure.route === 'string' && dose(a.exposure.dose)));
  const arms = new Map(Object.entries(n.arms).filter(([arm, a]) => readsAs(a, (d) => sameDose(arm, d))));
  const pk = new Map(Object.entries(n.pk && typeof n.pk === 'object' ? n.pk : {})
    .map(([arm, p]): [string, DrugPk | null] => [arm, pkOfSaved(p)])
    .filter((e): e is [string, DrugPk] => e[1] !== null));
  // a combination's drugs, each while the name still says the dose NCA read for it
  const doseOf = (arm: string, key: string) => parseArmName(arm).find((q) => q.key === key)?.dose ?? null;
  const sameDrugDose = (arm: string, key: string, dose: string) => {
    const [named, read] = [doseOf(arm, key), parseArmName(`x ${dose}`)[0]?.dose ?? null];
    return !named || !read || same(named, read);
  };
  const drugs = new Map(Object.entries(n.drugs && typeof n.drugs === 'object' ? n.drugs : {})
    .filter(([, ds]) => Array.isArray(ds))
    .map(([arm, ds]): [string, DrugExposure[]] => [arm, ds.filter((d) => d && typeof d.key === 'string' &&
      typeof d.drug === 'string' && readsAs(d.x, (dose) => sameDrugDose(arm, d.key, dose)))
      .map((d) => ({key: d.key, drug: d.drug, x: d.x, pk: pkOfSaved(d.pk)}))])
    .filter(([, ds]) => ds.length));
  return {humans, nca: {table: n.table, arms, pk, ...drugs.size ? {drugs} : {}}};
}

/** What a study keeps, with the patients' exposure set to `humans` (none: removed) and all else as it is. */
export function humansInto(x: SavedExposure | undefined, humans: Map<string, HumanExposure>): SavedExposure {
  const kept = exposureJson(humans, null)?.humans;
  const {humans: _, ...rest} = x ?? {};
  return kept ? {...rest, humans: kept} : rest;
}

/** What a study keeps, with its groups renamed (`renames`, old name to new: their NCA exposure, their PK and the
 *  analysis's TGI). */
export function renamedExposure(x: SavedExposure | undefined, renames: Map<string, string>):
SavedExposure | undefined {
  if (!x || !renames.size)
    return x;
  const swap = <T>(o: {[k: string]: T} | undefined) => o && Object.fromEntries(Object.entries(o)
    .map(([k, v]) => [renames.get(k) ?? k, v]));
  const effect = x.effect && {...x.effect, tgi: swap(x.effect.tgi)!, ...x.effect.ctd ? {ctd: {...x.effect.ctd,
    groups: x.effect.ctd.groups.map((g) => ({...g, name: renames.get(g.name) ?? g.name}))}} : {}};
  return {...x, ...effect ? {effect} : {},
    ...x.nca ? {nca: {...x.nca, arms: swap(x.nca.arms)!, pk: swap(x.nca.pk)!,
      ...x.nca.drugs ? {drugs: swap(x.nca.drugs)} : {}}} : {}};
}

/** A treated group as the read-out reads it: its TGI (%), the AUC the multiple uses (aucUsed) and, with the
 *  patients' exposure, its
 *  multiple (cvKnown: the patients' CV was given, so "their range" is theirs, not the 3-fold convention). */
export interface ReadoutGroup {arm: string, tgi: number | null, auc: number | null, aucUnit: string,
  multiple: {ratio: number, range?: 'above' | 'below', cvKnown: boolean, weekly?: boolean} | null}

/** One sentence on what the groups show of exposure and effect, from numbers already computed: with the patients'
 *  exposure, the best group within their range and the group furthest above it; without it, the TGI at the lowest
 *  and highest dose with their AUC; null with neither (or too little to say). It states what was measured, nothing
 *  beyond: no direction, no cause. */
export function readout(groups: ReadoutGroup[], meta: StudyMeta | null = null): string | null {
  const usable = groups.map((g) => ({g, p: armParts(g.arm, meta)}))
    .filter(({g, p}) => p.length === 1 && p[0].dose && g.tgi !== null && Number.isFinite(g.tgi) &&
      (g.auc !== null || g.multiple))
    .map(({g, p}) => ({...g, tgi: g.tgi!, part: p[0], dose: `${sig(p[0].dose!.value)} ${p[0].dose!.unit}`}));
  if (!usable.length)
    return null;
  // the drug with the most groups (each drug its own sentence would be a paragraph)
  const byDrug = new Map<string, typeof usable>();
  usable.forEach((g) => byDrug.set(g.part.key, [...byDrug.get(g.part.key) ?? [], g]));
  const drug = [...byDrug.values()].reduce((a, b) => b.length > a.length ? b : a);
  const lead = byDrug.size > 1 ? `${drug[0].part.drug}: ` : '';
  const tgi = (x: number) => `${Math.round(x)}%`;
  const best = <T extends {tgi: number}>(xs: T[]) => xs.reduce((a, b) => b.tgi > a.tgi ? b : a);
  // the no-dosing-at-weekends rule lowered the doses a day: the multiple is a week's average, said where it is read
  const week = (g: {multiple: ReadoutGroup['multiple']}) => g.multiple?.weekly ? ', averaged over the week' : '';
  const placed = drug.filter((g) => g.multiple);
  if (placed.length) {
    const within = placed.filter((g) => !g.multiple!.range);
    if (within.length) {
      const w = best(within);
      // the group furthest above the range, beside it
      const above = placed.filter((g) => g.multiple!.range === 'above');
      const m = above.length ? above.reduce((a, b) => b.multiple!.ratio > a.multiple!.ratio ? b : a) : null;
      return `${lead}At an exposure within ${w.multiple!.cvKnown ? 'the patients\' range' :
        '3-fold of the patients\' exposure'} (${multipleText(w.multiple!.ratio)}${week(w)}), ${w.dose} gave a TGI of ` +
        `${tgi(w.tgi)}${m ? `; ${m.dose} gave ${tgi(m.tgi)} at ${multipleText(m.multiple!.ratio)} the patients' ` +
          `exposure${week(m)}` : ''}.`;
    }
    const b = best(placed);
    const where = b.multiple!.cvKnown ? `${b.multiple!.range} the patients' range` :
      `more than 3-fold ${b.multiple!.range} it`;
    return `${lead}${b.dose} gave a TGI of ${tgi(b.tgi)} at ${multipleText(b.multiple!.ratio)} the patients' ` +
      `exposure${week(b)}, ${where}.`;
  }
  const unit = unitText(drug[0].part.dose!.unit);
  const withAuc = drug.filter((g) => g.auc !== null && unitText(g.part.dose!.unit) === unit)
    .sort((a, b) => a.part.dose!.value - b.part.dose!.value);
  if (withAuc.length < 2 || withAuc[0].part.dose!.value === withAuc[withAuc.length - 1].part.dose!.value)
    return null;
  const [lo, hi] = [withAuc[0], withAuc[withAuc.length - 1]];
  const aucs = withAuc.map((g) => g.auc!);
  return `${lead}The groups gave a TGI of ${tgi(lo.tgi)} at ${sig(lo.part.dose!.value)} and ${tgi(hi.tgi)} at ` +
    `${hi.dose} (AUC ${Number(Math.min(...aucs).toPrecision(2))}–${Number(Math.max(...aucs).toPrecision(2))} ` +
    `${lo.aucUnit}).`;
}

/** The one thing that would add most to the groups card, or null when it has it all: the mouse exposure (NCA
 *  installed, none yet), then the patients' (`inTable`: drugs of the groups PMX Studio's reference table holds). */
export function exposureHint(o: {ncaInstalled: boolean, mousePk: boolean, patients: boolean, inTable: string[],
  placed?: boolean, why?: string}): {text: string, add: 'nca' | 'patients' | null} | null {
  if (!o.mousePk)
    return o.ncaInstalled ? {text: 'Add the mouse concentrations to see each group\'s exposure', add: 'nca'} : null;
  // the patients' exposure is in: nothing to ask, but why no group is placed against it when none is
  if (o.patients) {
    return o.placed !== false ? null : {add: null, text: 'No group is placed against the patients\' exposure' +
      `${o.why ? `: ${o.why}` : ''}`};
  }
  const t = o.inTable;
  return {add: 'patients', text: 'Add the patients\' exposure (a projection is fine) to place each group against ' +
    `it${t.length ? `: ${t.join(', ')} ${t.length === 1 ? 'is' : 'are'} in PMX Studio's reference table, one ` +
      'click' : ''}`};
}

/** The mouse dose of `arm`'s drug whose daily AUC matches the patients' daily AUC, from an earlier study's kept
 *  exposure (`kept`, study `study`): NCA's AUC of one dose of its groups of the drug alone (AUCinf at a dose that has
 *  it, else AUClast, never both) times the doses a day, against the patients' (over their cycle), unbound when both
 *  fractions are known. The route and schedule are the arm's, else the study's when its groups share one (named in
 *  the text); none, or several: nothing. Two doses or more: log-dose against log-AUC between the two around the
 *  patients' value, none beyond them; one dose: proportional. `rank` orders studies: bracketing first, then the
 *  doses with PK. null: nothing to say. It places the patients' exposure on the mouse's dose scale, nothing more. */
export function patientsDose(arm: string, kept: SavedExposure | undefined, study: string,
  meta: StudyMeta | null = null, rules: DosingRules | null = null): {text: string, tip: string, rank: number} | null {
  const parts = parseArmName(arm);
  if (parts.length !== 1)
    return null;
  const [p] = parts;
  const {humans, nca} = parseExposure(kept);
  const h = humans.get(p.key);
  const daily = h && h.auc !== null && h.auc > 0 ? patientsDaily(h) : null;
  if (!h || daily === null || !nca)
    return null;
  // the drug's groups alone, by the arm's route when it names one, in one dose unit
  const groups: {dose: number, unit: string, route: string, schedule: string | null, auc: number, inf: boolean}[] = [];
  for (const [name, a] of nca.arms) {
    const q = armParts(name, meta);
    const x = a.exposure;
    const dose = x ? parseArmName(`x ${x.dose}`)[0]?.dose : null;
    if (!x || q.length !== 1 || q[0].key !== p.key || !dose || !(dose.value > 0) || !sameRoute(p.route, x.route))
      continue;
    const inf = x.aucInf !== null && x.aucInf > 0;
    const f = aucFactor(inf ? x.aucInfUnit || x.aucUnit : x.aucUnit);
    const auc = inf ? x.aucInf! : x.auc;
    if (f !== null && auc > 0 && (!groups.length || unitText(dose.unit) === unitText(groups[0].unit))) {
      groups.push({dose: dose.value, unit: dose.unit, route: x.route.startsWith('IV') ? 'IV' : x.route,
        schedule: q[0].schedule ?? x.schedule, auc: auc * f, inf});
    }
  }
  const routes = [...new Set(groups.map((g) => g.route))];
  const schedules = [...new Set(groups.map((g) => g.schedule).filter((s): s is string => !!s))];
  const schedule = p.schedule ?? (schedules.length === 1 ? schedules[0] : null);
  // the doses a day of the plan the arm is in, under its rules (not the earlier study's)
  const perDay = dosesPerDay(schedule, rules);
  if (!groups.length || (p.route === null && routes.length !== 1) || perDay === null)
    return null;
  const unit = groups[0].unit;
  const route = p.route === null ? ` ${routes[0]}` : '';
  const on = p.schedule ? `${schedule}${route} in mice` : `${route.trim() ? `${route.trim()} ` : ''}in mice on ` +
    `${study}'s ${schedule} schedule`;
  const unbound = fuOk(h.fuMouse) && fuOk(h.fuHuman);
  const cycle = cycleOf(h);
  const patients = patientsDaily(h, unbound)!;
  // per dose, AUCinf where any group has it, else AUClast: one kind per dose
  const pts = [...new Set(groups.map((g) => g.dose))].sort((a, b) => a - b).map((d) => {
    const at = groups.filter((g) => g.dose === d);
    const inf = at.some((g) => g.inf);
    const used = at.filter((g) => g.inf === inf);
    const auc = used.reduce((s, g) => s + g.auc, 0) / used.length;
    return {d, inf, auc, daily: auc * perDay * (unbound ? h.fuMouse! : 1)};
  });
  const round = (x: number) => Number(x.toPrecision(2));
  const who = `${h.source === 'projected' ? 'projected ' : ''}patients' exposure`;
  const Who = who[0].toUpperCase() + who.slice(1);
  const tip = [`Mouse AUC of one dose from "${study}": ${pts.map((x) => `${x.d} ${unit}, ${x.inf ? 'AUCinf' :
    'AUClast'} ${sig(x.auc)} mg/L·h`).join('; ')}; × ${sig(perDay)} dose${perDay === 1 ? '' : 's'} a day ` +
    `(${schedule}).`, ...pts.some((x) => !x.inf) ? ['AUClast: may underestimate exposure, so the dose may be ' +
      'lower.'] : [], `Patients' daily AUC, ${h.source}: ${sig(daily)} mg/L·h` +
    `${cycle ? ` (AUC per administration ÷ ${cycle} days between administrations, on average)` : ''}.`,
  unbound ? `On unbound exposure: unbound fraction ${h.fuMouse} in mouse and ${h.fuHuman} in patients.` :
    'On total exposure (both unbound fractions are needed to compare unbound exposure).'].join(' ');
  const from = `from ${study}'s PK`;
  // the mouse dose that gives the patients' exposure: proportional from a single dose, log-log between two, none
  // beyond the doses with PK
  let rank = 0;
  let equal: number | null = null;
  let old: string;
  if (pts.length === 1) {
    rank = 1;
    equal = pts[0].d * patients / pts[0].daily;
    old = `${Who} ≈ ${round(equal)} ${unit} ${on} (assuming exposure proportional to dose, ${from} at ` +
      `${pts[0].d} ${unit})`;
  } else {
    const i = pts.findIndex((a, k) => k + 1 < pts.length && a.daily !== pts[k + 1].daily &&
      patients >= Math.min(a.daily, pts[k + 1].daily) && patients <= Math.max(a.daily, pts[k + 1].daily));
    if (i >= 0) {
      const [a, b] = [pts[i], pts[i + 1]];
      rank = 100 + pts.length;
      equal = Math.exp(Math.log(a.d) + (Math.log(patients) - Math.log(a.daily)) * (Math.log(b.d) - Math.log(a.d)) /
        (Math.log(b.daily) - Math.log(a.daily)));
      old = `${Who} ≈ ${round(equal)} ${unit} ${on} (between ${a.d} and ${b.d} ${unit}, ${from})`;
    } else {
      const above = patients > Math.max(...pts.map((x) => x.daily));
      const edge = above ? pts[pts.length - 1] : pts[0];
      old = `The ${who} is ${above ? 'above' : 'below'} what ${edge.d} ${unit} ${on} reached ` +
        `(${multipleText(edge.daily / patients)}, ${from})`;
    }
  }
  // an arm that says its dose: said of that dose, in plain words; the study, route and doses with PK on hover
  const dose = p.dose && p.dose.value > 0 && unitText(p.dose.unit) === unitText(unit) ? p.dose.value : null;
  if (dose === null)
    return {tip, rank, text: old};
  // the arm's daily exposure: log-log between the doses with PK, in step with the dose from the nearest beyond them
  const up = pts.findIndex((x) => x.d >= dose);
  const k = Math.max(0, Math.min(pts.length - 2, (up < 0 ? pts.length - 1 : up) - 1));
  const [a, b] = pts.length === 1 ? [pts[0], pts[0]] : [pts[k], pts[k + 1]];
  const beyond = pts.length === 1 || dose < a.d || dose > b.d;
  const near = dose < a.d ? a : b;
  const mine = beyond ? near.daily * dose / near.d : Math.exp(Math.log(a.daily) + (Math.log(dose) - Math.log(a.d)) *
    (Math.log(b.daily) - Math.log(a.daily)) / (Math.log(b.d) - Math.log(a.d)));
  return {rank, tip: `${tip} Mouse doses with PK: ${pts.map((x) => x.d).join(', ')} ${unit} (${on}, ${from})` +
    `${beyond ? '; this dose is beyond them, its exposure taken in step with the dose' : ''}.`,
  text: `At ${dose} ${unit}, mice get about ${multipleText(mine / patients)} the ${who}` +
    `${equal === null ? '' : `; ${round(equal)} ${unit} would match it`}`};
}

/** The patients' exposure as a mouse dose for `arm`, from the earlier study with the highest rank (on a tie the
 *  first, in the order given); with none, and an arm that gives a dose but no schedule when a schedule would place
 *  it, a hint to add one; null: nothing to say. */
export function patientsDoseOf(arm: string, studies: PreviousStudy[], rules: DosingRules | null = null):
  {text: string, tip: string} | null {
  let best: ReturnType<typeof patientsDose> = null;
  for (const s of studies) {
    const d = patientsDose(arm, s.exposure, s.name, s.meta, rules);
    if (d && (!best || d.rank > best.rank))
      best = d;
  }
  if (best)
    return {text: best.text, tip: best.tip};
  const p = parseArmName(arm);
  const placed = p.length === 1 && p[0].dose && !p[0].schedule ?
    studies.find((s) => patientsDose(`${arm} QD`, s.exposure, s.name, s.meta, rules)) : undefined;
  return placed ? {text: 'Add the schedule (QD, BID, Q3D, QW…) to the arm\'s details or name to place it against ' +
    'the patients\' exposure', tip: `${placed.name} kept the patients' exposure and the mouse PK of this drug: ` +
    'with a schedule, Plan gives the mouse dose that matches it'} : null;
}

/** The analysis's effect a study kept with its exposure (null: none, or damaged). */
export function savedEffect(x: SavedExposure | undefined): {day: number | null, tgi: Map<string, number>} | null {
  const e = x?.effect;
  if (!e || typeof e !== 'object' || !e.tgi || typeof e.tgi !== 'object')
    return null;
  return {day: typeof e.day === 'number' && Number.isFinite(e.day) ? e.day : null, tgi: new Map(Object.entries(e.tgi)
    .filter((t): t is [string, number] => typeof t[1] === 'number' && Number.isFinite(t[1])))};
}
