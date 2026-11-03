/* A followed study: a .pmxstudy file (JSON) with its plan as it was when the study started, the reviews fixed then, its
   limit, its measurements and the log of what was decided and who revealed the groups. Shared from Files like a plan.
   It holds the allocation: hiding the groups in Monitor is a display control, not an access control. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import type {FuncCallForm} from '@datagrok-libraries/u2/dg';

import {PLAN_TAG, SavedPlan, parsePlan} from './plan-file';

export const STUDY_EXT = 'pmxstudy';
const FORMAT = 'pmx-study';
/** The tags of a measurements table that carry the study's reviews and log. */
export const LOOKS_TAG = 'pmx.looks';
export const LOG_TAG = 'pmx.log';
/** The name of the table's body weight column, and the protocol's limits ({humane, weightLimit}). */
export const WEIGHT_TAG = 'pmx.weight';
export const LIMITS_TAG = 'pmx.limits';

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
  /** `weight`: only when the study has body weights */
  rows: {animal: string, group: string | null, day: number, volume: number | null, weight?: number | null}[];
  log: LogEntry[];
}

export const logOf = (df: DG.DataFrame): LogEntry[] => JSON.parse(df.getTag(LOG_TAG) || '[]');

/** What a table says of itself into a form on it: its body weight column and the protocol's limits. */
export function restoreFields(form: FuncCallForm, df: DG.DataFrame): void {
  const limits = JSON.parse(df.getTag(LIMITS_TAG) || '{}');
  const fields = {bodyWeight: df.getTag(WEIGHT_TAG) || null, humane: limits.humane ?? null,
    weightLimit: limits.weightLimit ?? null};
  for (const [name, value] of Object.entries(fields)) {
    const input = form.getInput(name);
    if (input && value !== null)
      input.value = value;
  }
}

/** A log line, on the table it is about: when, who, what. */
export function addLog(df: DG.DataFrame, what: string): LogEntry[] {
  const log = [...logOf(df), {when: new Date().toISOString(), who: grok.shell.user?.login ?? '', what}];
  df.setTag(LOG_TAG, JSON.stringify(log));
  return log;
}

/** The measurements table of a plan just followed: empty, carrying its plan, ready to paste into. */
export function newStudyTable(plan: SavedPlan): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([DG.Column.string('animal', 0), DG.Column.string('group', 0),
    DG.Column.float('day', 0), DG.Column.float('volume', 0)]);
  df.name = plan.name;
  df.setTag(PLAN_TAG, JSON.stringify(plan));
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
  return {format: FORMAT, version: 1, name: df.name, plan: planText ? parsePlan(planText) : null,
    looks: df.getTag(LOOKS_TAG) ?? '', humane, weightLimit, log: logOf(df),
    rows: Array.from({length: df.rowCount}, (_, i) => ({animal: String(a.get(i)), group: g.isNone(i) ? null : g.get(i),
      day: d.get(i), volume: v.isNone(i) ? null : v.get(i), ...w ? {weight: w.isNone(i) ? null : w.get(i)} : {}}))};
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
    DG.Column.fromFloat64Array('day', Float64Array.from(s.rows, (r) => r.day)),
    DG.Column.fromFloat64Array('volume', Float64Array.from(s.rows, (r) => r.volume ?? NaN))]);
  if (s.rows.some((r) => 'weight' in r)) {
    df.columns.add(DG.Column.fromFloat64Array('body_weight', Float64Array.from(s.rows, (r) => r.weight ?? NaN)));
    df.setTag(WEIGHT_TAG, 'body_weight');
  }
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

/** Saves the study to the user's studies folder (replacing one of the same name); returns its path. */
export async function saveStudy(s: SavedStudy): Promise<string> {
  const path = studiesFolder() + fileName(s.name);
  await grok.dapi.files.writeAsText(path, JSON.stringify(s, null, 2));
  return path;
}
