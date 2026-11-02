/* A followed study: a .pmxstudy file (JSON) with its plan as it was when the study started, the reviews fixed then, its
   limit, its measurements and the log of what was decided and who revealed the groups. Shared from Files like a plan.
   It holds the allocation: hiding the groups in Monitor is a display control, not an access control. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {PLAN_TAG, SavedPlan, parsePlan} from './plan-file';

export const STUDY_EXT = 'pmxstudy';
const FORMAT = 'pmx-study';
/** The tags of a measurements table that carry the study's reviews and log. */
export const LOOKS_TAG = 'pmx.looks';
export const LOG_TAG = 'pmx.log';

export interface LogEntry {when: string, who: string, what: string}

export interface SavedStudy {
  format: typeof FORMAT;
  version: 1;
  name: string;
  plan: SavedPlan | null;
  looks: string;
  humane: number | null;
  rows: {animal: string, group: string | null, day: number, volume: number | null}[];
  log: LogEntry[];
}

export const logOf = (df: DG.DataFrame): LogEntry[] => JSON.parse(df.getTag(LOG_TAG) || '[]');

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

export function studyOf(df: DG.DataFrame, humane: number | null): SavedStudy {
  const planText = df.getTag(PLAN_TAG);
  const [a, g, d, v] = ['animal', 'group', 'day', 'volume'].map((c) => {
    const col = df.col(c);
    if (!col)
      throw new Error(`The study table has no column "${c}"`);
    return col;
  });
  return {format: FORMAT, version: 1, name: df.name, plan: planText ? parsePlan(planText) : null,
    looks: df.getTag(LOOKS_TAG) ?? '', humane, log: logOf(df),
    rows: Array.from({length: df.rowCount}, (_, i) => ({animal: String(a.get(i)), group: g.isNone(i) ? null : g.get(i),
      day: d.get(i), volume: v.isNone(i) ? null : v.get(i)}))};
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

/** The measurements table of a saved study, with its plan, reviews and log in tags. */
export function tableOf(s: SavedStudy): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('animal', s.rows.map((r) => r.animal)),
    DG.Column.fromStrings('group', s.rows.map((r) => r.group ?? '')),
    DG.Column.fromFloat64Array('day', Float64Array.from(s.rows, (r) => r.day)),
    DG.Column.fromFloat64Array('volume', Float64Array.from(s.rows, (r) => r.volume ?? NaN))]);
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
