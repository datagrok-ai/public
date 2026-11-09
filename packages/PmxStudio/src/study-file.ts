/* A followed study: a .pmxstudy file (JSON) with its plan as it was when the study started, the reviews fixed then, its
   limit, its measurements and the log of what was decided and who revealed the groups. Shared from Files like a plan.
   It holds the allocation: hiding the groups in Monitor is a display control, not an access control. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import type {FuncCallForm} from '@datagrok-libraries/u2/dg';

import {PLAN_TAG, SavedPlan, parsePlan} from './plan-file';
import {StudyMeta, metaOf, parseMeta, setMetaOn} from './study-meta';

export const STUDY_EXT = 'pmxstudy';
const FORMAT = 'pmx-study';
/** The tags of a measurements table that carry the study's reviews and log. */
export const LOOKS_TAG = 'pmx.looks';
export const LOG_TAG = 'pmx.log';
/** The name of the table's body weight column, and the protocol's limits ({humane, weightLimit}). */
export const WEIGHT_TAG = 'pmx.weight';
export const LIMITS_TAG = 'pmx.limits';
/** The marks, visible columns of the study table (study-table.ts sets them): ticked in Excluded, out of PMX Studio's
 *  analysis (the row stays); ticked in Flagged, to look at; why, in Mark reason. */
export const EXCLUDED = 'Excluded';
export const FLAGGED = 'Flagged';
export const MARK_REASON = 'Mark reason';

/** Whether a mark column is ticked in a row (a bool column, or a file's own text: true, yes, y, 1, x). */
export function ticked(col: DG.Column | null | undefined, i: number): boolean {
  const v = col?.get(i);
  return v === true || (typeof v === 'string' && /^(true|yes|y|1|x)$/i.test(v.trim())) || v === 1;
}

export interface LogEntry {when: string, who: string, what: string}

export interface SavedStudy {
  format: typeof FORMAT;
  version: 1;
  name: string;
  plan: SavedPlan | null;
  looks: string;
  humane: number | null;
  /** the % of body weight lost at which an animal is removed (older files: none) */
  weightLimit?: number | null;
  /** the study's ID, owner, dates and what each group receives (older files: none) */
  meta?: StudyMeta;
  /** `weight`: only when the study has body weights; the marks only on marked rows */
  rows: {animal: string, group: string | null, day: number | null, volume: number | null, weight?: number | null,
    excluded?: true, flagged?: true, reason?: string}[];
  log: LogEntry[];
}

export const logOf = (df: DG.DataFrame): LogEntry[] => JSON.parse(df.getTag(LOG_TAG) || '[]');

/** What a table says of itself into a form on it: its plan's tumor model, its body weight column and the protocol's
 *  limits. */
export function restoreFields(form: FuncCallForm, df: DG.DataFrame): void {
  const limits = JSON.parse(df.getTag(LIMITS_TAG) || '{}');
  let line: string | null = null;
  try {
    line = df.getTag(PLAN_TAG) ? parsePlan(df.getTag(PLAN_TAG)!).predicted.line : null;
  } catch {
    // a damaged plan: Analyze says so; the field stays as it is
  }
  const fields = {line, bodyWeight: df.getTag(WEIGHT_TAG) || null, humane: limits.humane ?? null,
    weightLimit: limits.weightLimit ?? null};
  for (const [name, value] of Object.entries(fields)) {
    const input = form.getInput(name);
    if (input && value !== null)
      input.value = value;
  }
}

export const REVEALED = 'Groups revealed';

/** Whether the groups may be shown: always, unless the study is blinded and no one has revealed them yet. */
export const groupsShown = (df: DG.DataFrame): boolean =>
  !metaOf(df)?.blinded || logOf(df).some((e) => e.what === REVEALED);

/** A log line, on the table it is about: when, who, what. */
export function addLog(df: DG.DataFrame, what: string): LogEntry[] {
  const log = [...logOf(df), {when: new Date().toISOString(), who: grok.shell.user?.login ?? '', what}];
  df.setTag(LOG_TAG, JSON.stringify(log));
  return log;
}

/** Volumes and weights to one decimal in every cell and tooltip (Crux's and the caliper's floats carry a dozen). */
export function oneDecimal(df: DG.DataFrame): DG.DataFrame {
  for (const c of df.columns.numerical) {
    if (c.type === DG.TYPE.FLOAT && !c.meta.format && /volume|weight/i.test(c.name))
      c.meta.format = '0.0';
  }
  return df;
}

/** The measurements table of a plan just followed: empty, carrying its plan, ready to paste into. */
export function newStudyTable(plan: SavedPlan): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([DG.Column.string('animal', 0), DG.Column.string('group', 0),
    DG.Column.float('day', 0), DG.Column.float('volume', 0)]);
  oneDecimal(df).name = plan.name;
  df.setTag(PLAN_TAG, JSON.stringify(plan));
  if (plan.study)
    setMetaOn(df, plan.study);
  return df;
}

export function studyOf(df: DG.DataFrame, humane: number | null, weightLimit: number | null = null): SavedStudy {
  const planText = df.getTag(PLAN_TAG);
  const [a, g, d, v] = ['animal', 'group', 'day', 'volume'].map((c) => {
    const col = df.col(c);
    if (!col)
      throw new Error(`The study table has no column "${c}"`);
    return col;
  });
  const weightName = df.getTag(WEIGHT_TAG);
  const w = weightName ? df.col(weightName) : null;
  const [e, f, r] = [df.col(EXCLUDED), df.col(FLAGGED), df.col(MARK_REASON)];
  const meta = metaOf(df);
  return {format: FORMAT, version: 1, name: df.name, plan: planText ? parsePlan(planText) : null,
    looks: df.getTag(LOOKS_TAG) ?? '', humane, weightLimit, log: logOf(df), ...meta ? {meta} : {},
    rows: Array.from({length: df.rowCount}, (_, i) => ({animal: String(a.get(i)), group: g.isNone(i) ? null : g.get(i),
      day: d.isNone(i) ? null : d.get(i), volume: v.isNone(i) ? null : v.get(i),
      ...w ? {weight: w.isNone(i) ? null : w.get(i)} : {},
      ...ticked(e, i) ? {excluded: true as const} : {}, ...ticked(f, i) ? {flagged: true as const} : {},
      ...r && !r.isNone(i) && String(r.get(i)).trim() ? {reason: String(r.get(i))} : {}}))};
}

/** A study read from a file someone may have edited or sent: checked before use. */
export function parseStudy(text: string): SavedStudy {
  let s: any;
  try {
    s = JSON.parse(text);
  } catch {
    throw new Error('Not a PMX Studio study: the file is not JSON');
  }
  if (s?.format !== FORMAT)
    throw new Error('Not a PMX Studio study');
  if (s.version !== 1)
    throw new Error(`A study of version ${s.version}, which this PMX Studio cannot read: update the package`);
  if (!Array.isArray(s.rows) || !Array.isArray(s.log))
    throw new Error('The study file is damaged');
  if (s.plan)
    parsePlan(JSON.stringify(s.plan));
  return s as SavedStudy;
}

/** The measurements table of a saved study, with its plan, reviews, limits and log in tags. */
export function tableOf(s: SavedStudy): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('animal', s.rows.map((r) => r.animal)),
    DG.Column.fromStrings('group', s.rows.map((r) => r.group ?? '')),
    DG.Column.fromFloat64Array('day', Float64Array.from(s.rows, (r) => r.day ?? DG.FLOAT_NULL)),
    DG.Column.fromFloat64Array('volume', Float64Array.from(s.rows, (r) => r.volume ?? DG.FLOAT_NULL))]);
  if (s.rows.some((r) => 'weight' in r)) {
    df.columns.add(DG.Column.fromFloat64Array('body_weight',
      Float64Array.from(s.rows, (r) => r.weight ?? DG.FLOAT_NULL)));
    df.setTag(WEIGHT_TAG, 'body_weight');
  }
  oneDecimal(df);
  if (s.rows.some((r) => r.excluded || r.flagged || r.reason)) {
    df.columns.add(DG.Column.fromList(DG.TYPE.BOOL, EXCLUDED, s.rows.map((r) => r.excluded === true)));
    df.columns.add(DG.Column.fromList(DG.TYPE.BOOL, FLAGGED, s.rows.map((r) => r.flagged === true)));
    df.columns.add(DG.Column.fromStrings(MARK_REASON, s.rows.map((r) => r.reason ?? '')));
  }
  const meta = parseMeta(s.meta);
  if (meta)
    setMetaOn(df, meta);
  df.setTag(LIMITS_TAG, JSON.stringify({humane: s.humane, weightLimit: s.weightLimit ?? null}));
  df.name = s.name;
  if (s.plan)
    df.setTag(PLAN_TAG, JSON.stringify(s.plan));
  df.setTag(LOOKS_TAG, s.looks);
  df.setTag(LOG_TAG, JSON.stringify(s.log));
  return df;
}

export const studiesFolder = () => `${grok.shell.user.project.name}:Home/PMX Studio/studies/`;
const fileName = (name: string) => `${name.replace(/[\\/:*?"<>|]+/g, ' ').trim() || 'Study'}.${STUDY_EXT}`;

/** Saves the study to the user's studies folder (replacing the same study's file); returns its path. A different
 *  study saved under the same name (another ID) is not replaced: the user is asked to rename this one. */
export async function saveStudy(s: SavedStudy): Promise<string> {
  const path = studiesFolder() + fileName(s.name);
  const other = await grok.dapi.files.exists(path) ?
    await grok.dapi.files.readAsText(path).then((t) => JSON.parse(t)?.meta?.id ?? null, () => null) : null;
  if (other && s.meta?.id && other !== s.meta.id)
    throw new Error(`Another study (${other}) is saved as "${s.name}": rename this one to save it`);
  await grok.dapi.files.writeAsText(path, JSON.stringify(s, null, 2));
  return path;
}
