/* A followed study: a .pmxstudy file (JSON) with its plan as it was when the study started, the reviews fixed then, its
   limit, its measurements and the log of what was decided and who revealed the groups. Shared from Files like a plan.
   It holds the allocation: hiding the groups in Monitor is a display control, not an access control. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import type {FuncCallForm} from '@datagrok-libraries/u2/dg';

import {_package} from './package';
import {PLAN_TAG, Provenance, SavedPlan, fileNameOf, parsePlan, provenance, provenanceOf} from './plan-file';
import {StudyMeta, metaOf, parseMeta, setMetaOn} from './study-meta';
import {VolumeSource, volumeSourceOf} from './study-reader';

export const STUDY_EXT = 'pmxstudy';
const FORMAT = 'pmx-study';
export {FORMAT as STUDY_FORMAT};
/** The tag with the columns' roles of a study table, as JSON ({animal, group, day, volume, bodyWeight}): what makes a
 *  table a study (study-table.ts writes it on every run; detectors.js reads it too). */
export const STUDY_TAG = 'pmx.study';
/** The tags of a measurements table that carry the study's reviews and log. */
export const LOOKS_TAG = 'pmx.looks';
const LOG_TAG = 'pmx.log';
/** The name of the table's body weight column, and the protocol's limits ({humane, weightLimit}). */
export const WEIGHT_TAG = 'pmx.weight';
export const LIMITS_TAG = 'pmx.limits';
/** The tumor model chosen for the study, kept for a study with no plan (with a plan, the plan's model is read). */
export const LINE_TAG = 'pmx.line';
/** Where the table's volumes come from, as Open a study read them (JSON, study-reader.ts's VolumeSource). */
export const VOLUME_TAG = 'pmx.volume';
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
  /** the atlas and app the study was saved with (older files: none; none either when the atlas could not be read);
   *  its plan's predictions carry their own. `volume`: where its volumes come from, given, computed from length and
   *  width with a formula, or given and found to follow one (studies opened from a table since 1.x; older: none) */
  provenance?: Partial<Provenance> & {volume?: VolumeSource};
  /** the tumor model chosen in Monitor (older files: none); a plan's model comes first */
  line?: string;
  /** `weight`: only when the study has body weights; the marks only on marked rows */
  rows: {animal: string, group: string | null, day: number | null, volume: number | null, weight?: number | null,
    excluded?: true, flagged?: true, reason?: string}[];
  /** the table's other columns (caliper readings, cage, sex, comments), by name, a value per row (older files: none) */
  extra?: {[column: string]: (string | number | boolean | null)[]};
  log: LogEntry[];
}

/** Where a table's volumes come from (VOLUME_TAG), when Open a study read it; undefined otherwise. */
function volumeOfTag(df: DG.DataFrame): VolumeSource | undefined {
  try {
    return volumeSourceOf(JSON.parse(df.getTag(VOLUME_TAG) || 'null'));
  } catch {
    return undefined; // a damaged tag says nothing
  }
}

export const logOf = (df: DG.DataFrame): LogEntry[] => JSON.parse(df.getTag(LOG_TAG) || '[]');

/** The lab's humane limit (mm3), set once in the package's settings (humaneLimit): where a plan or a study with no
 *  limit of its own starts; null: none set, the form's default stands. */
export function labHumane(): number | null {
  const x = Number(_package.settings?.['humaneLimit']);
  return Number.isFinite(x) && x > 0 ? x : null;
}

/** A study table's tumor model: its plan's, else the one chosen for it in Monitor (LINE_TAG), else null. */
export function lineOf(df: DG.DataFrame): string | null {
  try {
    if (df.getTag(PLAN_TAG))
      return parsePlan(df.getTag(PLAN_TAG)!).predicted.line;
  } catch {
    // a damaged plan: Analyze says so; the study's own model, if any, stands
  }
  return df.getTag(LINE_TAG) || null;
}

/** What a table says of itself into a form on it: its tumor model, its body weight column and the protocol's
 *  limits. */
export function restoreFields(form: FuncCallForm, df: DG.DataFrame): void {
  const limits = JSON.parse(df.getTag(LIMITS_TAG) || '{}');
  const fields = {line: lineOf(df), bodyWeight: df.getTag(WEIGHT_TAG) || null, humane: limits.humane ?? labHumane(),
    weightLimit: limits.weightLimit ?? null};
  for (const [name, value] of Object.entries(fields)) {
    const input = form.getInput(name);
    if (input && value !== null)
      input.value = value;
  }
}

export const REVEALED = 'Groups revealed';

/** Whether the groups may be shown: always, unless the study is blinded and no one has revealed them yet. */
export const groupsShownIn = (meta: StudyMeta | null | undefined, log: LogEntry[]): boolean =>
  !meta?.blinded || log.some((e) => e.what === REVEALED);
export const groupsShown = (df: DG.DataFrame): boolean => groupsShownIn(metaOf(df), logOf(df));
/** The same for a saved study: its details, else its plan's (as metaOf reads a table), and its log. */
export const savedGroupsShown = (s: SavedStudy): boolean =>
  groupsShownIn(parseMeta(s.meta) ?? parseMeta(s.plan?.study), s.log);

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

/** A study table's columns by role, as Monitor or Analyze last ran on it (STUDY_TAG); a table neither has run on
 *  yet, by their usual names. */
export function roleColumns(df: DG.DataFrame): {animal: DG.Column, group: DG.Column, day: DG.Column,
  volume: DG.Column} {
  let roles: {[role: string]: unknown} = {};
  try {
    roles = JSON.parse(df.getTag(STUDY_TAG) || '{}') ?? {};
  } catch {
    // a damaged tag: the usual names
  }
  const col = (role: string) => {
    const c = (typeof roles[role] === 'string' ? df.col(roles[role] as string) : null) ?? df.col(role);
    if (!c)
      throw new Error(`The study table has no ${role} column: pick it in the form, then try again.`);
    return c;
  };
  return {animal: col('animal'), group: col('group'), day: col('day'), volume: col('volume')};
}

export function studyOf(df: DG.DataFrame, humane: number | null, weightLimit: number | null = null): SavedStudy {
  const planText = df.getTag(PLAN_TAG);
  const {animal: a, group: g, day: d, volume: v} = roleColumns(df);
  const weightName = df.getTag(WEIGHT_TAG);
  const w = weightName ? df.col(weightName) : null;
  const [e, f, r] = [df.col(EXCLUDED), df.col(FLAGGED), df.col(MARK_REASON)];
  const meta = metaOf(df);
  const line = df.getTag(LINE_TAG);
  const volume = volumeOfTag(df);
  // what else the table holds goes with the study: the caliper readings behind a volume, a cage, a comment
  const kept = new Set([a, g, d, v, w, e, f, r].filter((c) => c).map((c) => c!.name));
  const plain: string[] = [DG.TYPE.FLOAT, DG.TYPE.INT, DG.TYPE.STRING, DG.TYPE.BOOL];
  const extra = Object.fromEntries(df.columns.toList()
    .filter((c) => !kept.has(c.name) && !c.name.startsWith('~') && plain.includes(c.type))
    .map((c) => [c.name, Array.from({length: df.rowCount}, (_, i) => c.isNone(i) ? null : c.get(i))]));
  return {format: FORMAT, version: 1, name: df.name, plan: planText ? parsePlan(planText) : null,
    looks: df.getTag(LOOKS_TAG) ?? '', humane, weightLimit, log: logOf(df), ...meta ? {meta} : {},
    ...line ? {line} : {}, ...volume ? {provenance: {volume}} : {},
    ...Object.keys(extra).length ? {extra} : {},
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
    throw new Error('Not a PMX Studio study: the file cannot be read (cut short, or another kind of file). ' +
      'Open the .pmxstudy file saved from Monitor.');
  }
  if (s?.format !== FORMAT)
    throw new Error('Not a PMX Studio study. Open the .pmxstudy file saved from Monitor.');
  if (s.version !== 1) {
    throw new Error('This study was saved by another version of PMX Studio, which this one cannot read. Ask your ' +
      'administrator to update the PMX Studio package.');
  }
  // what the views read without asking: a row's animal, a log line's time and text
  if (!Array.isArray(s.rows) || !Array.isArray(s.log) ||
    !s.rows.every((r: any) => r && typeof r.animal === 'string') ||
    !s.log.every((e: any) => e && typeof e.when === 'string' && typeof e.what === 'string')) {
    throw new Error('The study file is damaged: some of its rows or log lines are incomplete. Open another copy of ' +
      'it.');
  }
  if (s.plan)
    parsePlan(JSON.stringify(s.plan));
  // a field this version does not know is kept as it came; the model only as a name
  const {line, provenance: prov, ...rest} = s;
  const volume = volumeSourceOf(prov?.volume);
  const kept = {...provenanceOf(prov), ...volume ? {volume} : {}};
  return {...rest, ...Object.keys(kept).length ? {provenance: kept} : {},
    ...typeof line === 'string' && line ? {line} : {},
    name: typeof s.name === 'string' ? s.name : '', looks: typeof s.looks === 'string' ? s.looks : ''} as SavedStudy;
}

/** The atlas a plan (its file's text) was made with: undefined with no plan, null when it records none. */
export function planAtlas(text: string | null | undefined): string | null | undefined {
  try {
    const p = JSON.parse(text?.trim() || 'null');
    return p ? provenanceOf(p.provenance)?.atlas ?? null : undefined;
  } catch {
    return undefined; // a damaged plan: Analyze says so
  }
}

/** The study with what it is saved with (the atlas and app now); unchanged when the atlas cannot be read. */
export const stamped = async (s: SavedStudy): Promise<SavedStudy> => {
  const provenanceNow = await provenance().catch(() => undefined);
  const volume = s.provenance?.volume;
  return provenanceNow ? {...s, provenance: {...provenanceNow, ...volume ? {volume} : {}}} : s;
};

/** The measurements table of a saved study, with its plan, reviews, limits and log in tags. */
export function tableOf(s: SavedStudy): DG.DataFrame {
  // a value someone typed into the file as text is an empty cell, not NaN
  const floats = (name: string, values: unknown[]) => DG.Column.fromFloat64Array(name,
    Float64Array.from(values, (x) => typeof x === 'number' && Number.isFinite(x) ? x : DG.FLOAT_NULL));
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('animal', s.rows.map((r) => r.animal)),
    DG.Column.fromStrings('group', s.rows.map((r) => r.group ?? '')),
    floats('day', s.rows.map((r) => r.day)), floats('volume', s.rows.map((r) => r.volume))]);
  if (s.rows.some((r) => 'weight' in r)) {
    df.columns.add(floats('body_weight', s.rows.map((r) => r.weight)));
    df.setTag(WEIGHT_TAG, 'body_weight');
  }
  for (const [name, values] of Object.entries(s.extra ?? {})) {
    if (df.col(name) || !Array.isArray(values) || values.length !== s.rows.length)
      continue;
    // ponytail: a calculated column comes back as its values (its formula is not kept); numbers as floats
    df.columns.add(values.every((x) => x === null || typeof x === 'number') ? floats(name, values) :
      DG.Column.fromStrings(name, values.map((x) => x === null ? '' : String(x))));
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
  if (s.line)
    df.setTag(LINE_TAG, s.line);
  if (s.provenance?.volume)
    df.setTag(VOLUME_TAG, JSON.stringify(s.provenance.volume));
  df.setTag(LOOKS_TAG, s.looks);
  df.setTag(LOG_TAG, JSON.stringify(s.log));
  return df;
}

export const studiesFolder = () => `${grok.shell.user.project.name}:Home/PMX Studio/studies/`;
const fileName = (name: string) => fileNameOf(name, 'Study', STUDY_EXT);

/** Saves the study to the user's studies folder (replacing the same study's file); returns its path. A different
 *  study saved under the same name (another ID) is not replaced: the user is asked to rename this one. */
export async function saveStudy(s: SavedStudy): Promise<string> {
  const path = studiesFolder() + fileName(s.name);
  // a file there that does not read (cut short, edited by hand) is no one's study: saving replaces it
  const other = await grok.dapi.files.exists(path) ?
    await grok.dapi.files.readAsText(path).then((t) => JSON.parse(t)?.meta?.id ?? null).catch(() => null) : null;
  if (other && s.meta?.id && other !== s.meta.id)
    throw new Error(`Another study (${other}) is saved as "${s.name}": rename this one to save it`);
  await grok.dapi.files.writeAsText(path, JSON.stringify(await stamped(s), null, 2));
  // (its index row, now older than the file, is read again)
  return path;
}
