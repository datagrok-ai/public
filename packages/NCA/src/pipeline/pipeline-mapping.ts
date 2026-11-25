/**
 * Headless mapping for the pipeline functions (ncaDetectMapping / ncaCompute).
 *
 * The interactive app maps columns in the mapping pane; a pipeline has no person to
 * ask, so the mapping is a serialisable JSON object and detection reports a
 * confidence. Low confidence blocks the study instead of guessing.
 *
 * Column detection reuses `detectMapping` (the same synonym tables and unit-token
 * parsing as the pane); this module only adds what a pane does not need: the
 * grouping columns that make a profile key (study, compound, period) and the
 * ambiguity checks.
 *
 * **Never assume a unit or a route.** A person sees and picks units and route in
 * the app's form; a pipeline cannot, so a mapping that leaves either open blocks
 * every profile (`mappingGaps`) instead of falling back to a default. A malformed
 * mapping (bad JSON, unknown key, unknown unit) is a configuration error and throws.
 */

import * as DG from 'datagrok-api/dg';

import {detectMapping, isCompoundIdColumnName} from '../adapter/column-autodetect';
import {
  ROUTE_STRING_TO_CODE, CONC_TO_MG_L, TIME_TO_HOUR, DOSE_TO_MG, BW_TO_KG, SEX_SYNONYMS,
  defaultSexAssignment, uniqueSexValues, isOccasionName, splitsProfiles, numAt, normOccasionName,
} from '../adapter';
import type {BodyWeightUnit, ConcentrationUnit, DoseUnit, SexCode, TimeUnit} from '../adapter';

/** How samples were taken: `auto` decides per partition (see `pipeline-compute.ts`). */
export type SamplingDesign = 'auto' | 'serial' | 'sparse-destructive' | 'sparse-batch';
const DESIGNS: ReadonlyArray<SamplingDesign> = ['auto', 'serial', 'sparse-destructive', 'sparse-batch'];

export interface PipelineUnits {
  concentration: ConcentrationUnit;
  /** Applies to `time`, `nominalTime` and `infusionDuration`. */
  time: TimeUnit;
  dose: DoseUnit;
  bodyWeight: BodyWeightUnit;
}

/** Serialisable column mapping consumed by `ncaCompute`. */
export interface PipelineMapping {
  subject: string;
  time: string;
  concentration: string;
  dose?: string;
  route?: string;
  bodyWeight?: string;
  blq?: string;
  lloq?: string;
  /** Free-text CRO comment column; see `rules.reviewCommentPattern` for which comments route to review. */
  comment?: string;
  nominalTime?: string;
  /** Per-row infusion duration (in `units.time`); a bare `IV` route with a positive duration is an infusion. */
  infusionDuration?: string;
  /** Sex column. F% is stratified by sex, so without a recognised sex column there is no F%. */
  sex?: string;
  /** Columns that, with the subject, make the profile key (study, compound, period). */
  groupBy: string[];
  /** Units of the mapped quantities — never defaulted (see `mappingGaps`). `units.time` applies to
   *  `time`, `nominalTime` and `infusionDuration`. */
  units: Partial<PipelineUnits>;
  /** Route applied when no route column is mapped (one of `ROUTE_STRING_TO_CODE`'s keys). */
  defaultRoute?: string;
  /** Sampling design; default `auto`. */
  design?: SamplingDesign;
}

export type MappingConfidence = 'high' | 'low';

export interface MappingDetection {
  mapping: PipelineMapping | null;
  confidence: MappingConfidence;
  issues: string[];
}

/** `F_PCT_REASON` when the mapping has no sex column. */
export const NO_SEX_COLUMN = 'no sex column';

const norm = normOccasionName;

const STUDY_NAME = /^(study|studyid|studyno|studynumber|studycode)$/;
/** Occasion names a profile key was taken from unconditionally before the shared matcher
 *  (`isOccasionName`) existed. They stay unconditional so that no existing pipeline's profile
 *  key changes (AC-1.9 golden); the names the shared list adds (`StudyDay`, `Day`, `NOMDY`,
 *  `Week`) must also pass `splitsProfiles`, because a `Day` that is a per-row nominal day of
 *  one daily-sampled profile would otherwise split it. */
const LEGACY_PERIOD_KEY = /^(period|visitdy|visit|visitday|occasion|occ)$/;
const COMMENT_NAME = /^(comment|comments|remark|remarks|note|notes|crocomment)$/;
const NOMINAL_NAME = /^(nominal|nom|nominaltime|nomtime|nominaltime(h|hr|hrs|hours|min)?|pctptnum)$/;
const CONC_LIKE = /(conc|concentration|^dv$|^aval$)/;

/**
 * Grouping columns found by name, in profile-key order: study, compound, occasion. An occasion
 * column (the ONE matcher `isOccasionName`, shared with the app) becomes a profile key when it
 * is not the time column and, for a name outside the legacy period list, when it splits each
 * subject into genuine profiles (`splitsProfiles`: ≥ 3 times per block, time restarting per
 * block) — AC-1.9.
 */
export function detectGroupByColumns(df: DG.DataFrame, exclude: ReadonlySet<string>,
  roles: {subject?: string; time?: string} = {}): string[] {
  const cols = df.columns.toList().filter((c) => !exclude.has(c.name) && c.name !== roles.time);
  const pick = (test: (n: string) => boolean): string[] =>
    cols.filter((c) => test(norm(c.name))).map((c) => c.name);
  const subj = roles.subject !== undefined ? df.col(roles.subject) : null;
  const time = roles.time !== undefined ? df.col(roles.time) : null;
  const splits = (name: string): boolean => {
    if (subj === null || time === null) return false;
    const c = df.col(name)!;
    const n = df.rowCount;
    const t = new Float64Array(n);
    for (let i = 0; i < n; i++) t[i] = numAt(time, i);
    return splitsProfiles(
      Array.from({length: n}, (_, i) => subj.getString(i)),
      Array.from({length: n}, (_, i) => c.isNone(i) ? '' : c.getString(i)), t);
  };
  return [
    ...pick((n) => STUDY_NAME.test(n)),
    ...pick((n) => isCompoundIdColumnName(n)),
    ...cols.filter((c) => isOccasionName(c.name) &&
      (LEGACY_PERIOD_KEY.test(norm(c.name)) || splits(c.name))).map((c) => c.name),
  ];
}

function isNumeric(c: DG.Column): boolean {
  return c.type === DG.COLUMN_TYPE.FLOAT || c.type === DG.COLUMN_TYPE.INT ||
    c.type === DG.COLUMN_TYPE.BIG_INT || c.type === DG.COLUMN_TYPE.QNUM;
}

/**
 * Detect a pipeline mapping. `confidence` is `'low'` when a required role is
 * missing, a unit could not be read from a column name, there is no route column,
 * or more than one numeric column could be the concentration (e.g. a reassay
 * column). A low-confidence mapping is still returned when the required roles were
 * found, so a person can review it, but `ncaCompute` will not use it unattended.
 * Units that could not be read are left out of `units`, never defaulted.
 */
export function detectPipelineMapping(df: DG.DataFrame,
  opts: {perPartitionAmbiguity?: boolean} = {}): MappingDetection {
  const d = detectMapping(df);
  const issues: string[] = [];

  if (d.subject === undefined) issues.push('no subject column found');
  if (d.time === undefined) issues.push('no time column found');
  if (d.concentration === undefined) issues.push('no concentration column found');

  const concCandidates = df.columns.toList()
    .filter((c) => isNumeric(c) && CONC_LIKE.test(norm(c.name)) && !/lloq|lod/.test(norm(c.name)));
  // `ncaCompute` checks this per study instead (a reassay column in one study's rows
  // must not block the others), so it asks detection to skip the table-wide check.
  if (concCandidates.length > 1 && opts.perPartitionAmbiguity !== true) {
    issues.push(`ambiguous concentration: ${concCandidates.length} candidate columns ` +
      `(${concCandidates.map((c) => `'${c.name}'`).join(', ')})`);
  }

  if (d.concentration !== undefined && d.concentrationUnit === undefined)
    issues.push(`concentration unit not declared in '${d.concentration.name}'`);
  if (d.time !== undefined && d.timeUnit === undefined)
    issues.push(`time unit not declared in '${d.time.name}'`);
  if (d.dose !== undefined && d.doseUnit === undefined)
    issues.push(`dose unit not declared in '${d.dose.name}'`);
  if (d.dose === undefined) issues.push('no dose column found');
  if (d.bodyWeight !== undefined && d.bodyWeightUnit === undefined)
    issues.push(`body-weight unit not declared in '${d.bodyWeight.name}'`);
  if (d.route === undefined) issues.push('no route column found');

  if (d.subject === undefined || d.time === undefined || d.concentration === undefined)
    return {mapping: null, confidence: 'low', issues};

  const mapped = new Set<string>([d.subject, d.time, d.concentration, d.dose, d.route,
    d.bodyWeight, d.blqFlag, d.lloq].filter((c): c is DG.Column => c !== undefined).map((c) => c.name));
  const findBy = (re: RegExp, numeric?: boolean): string | undefined => df.columns.toList()
    .find((c) => !mapped.has(c.name) && re.test(norm(c.name)) &&
      (numeric === undefined || isNumeric(c) === numeric))?.name;
  const comment = findBy(COMMENT_NAME, false);
  const nominalTime = findBy(NOMINAL_NAME, true);
  if (comment !== undefined) mapped.add(comment);
  if (nominalTime !== undefined) mapped.add(nominalTime);
  // Sex only gates F% (no sex → no F%, with a reason), so it never lowers confidence.
  const sexCandidates = df.columns.toList().filter((c) => !mapped.has(c.name) && SEX_SYNONYMS.test(c.name) &&
    (c.type === DG.COLUMN_TYPE.STRING || c.type === DG.COLUMN_TYPE.INT));
  const sex = sexCandidates.length === 1 ? sexCandidates[0].name : undefined;
  if (sex !== undefined) mapped.add(sex);

  const units: Partial<PipelineUnits> = {
    ...(d.concentrationUnit !== undefined ? {concentration: d.concentrationUnit as ConcentrationUnit} : {}),
    ...(d.timeUnit !== undefined ? {time: d.timeUnit as TimeUnit} : {}),
    ...(d.doseUnit !== undefined ? {dose: d.doseUnit as DoseUnit} : {}),
    ...(d.bodyWeightUnit !== undefined ? {bodyWeight: d.bodyWeightUnit as BodyWeightUnit} : {}),
  };
  const mapping: PipelineMapping = {
    subject: d.subject.name,
    time: d.time.name,
    concentration: d.concentration.name,
    ...(d.dose !== undefined ? {dose: d.dose.name} : {}),
    ...(d.route !== undefined ? {route: d.route.name} : {}),
    ...(d.bodyWeight !== undefined ? {bodyWeight: d.bodyWeight.name} : {}),
    ...(d.blqFlag !== undefined ? {blq: d.blqFlag.name} : {}),
    ...(d.lloq !== undefined ? {lloq: d.lloq.name} : {}),
    ...(comment !== undefined ? {comment} : {}),
    ...(nominalTime !== undefined ? {nominalTime} : {}),
    ...(sex !== undefined ? {sex} : {}),
    groupBy: detectGroupByColumns(df, mapped, {subject: d.subject.name, time: d.time.name}),
    units,
  };
  return {mapping, confidence: issues.length === 0 ? 'high' : 'low', issues};
}

/** Parse a JSON configuration argument; malformed JSON is a configuration error. */
export function parseConfigJson(what: 'mapping' | 'rules', json: string): unknown {
  try {
    return JSON.parse(json);
  } catch (e) {
    throw new Error(`${what}: invalid JSON — ${e instanceof Error ? e.message : String(e)}`);
  }
}

const COLUMN_KEYS = ['subject', 'time', 'concentration', 'dose', 'route', 'bodyWeight', 'blq', 'lloq',
  'comment', 'nominalTime', 'infusionDuration', 'sex'] as const;
const MAPPING_KEYS: ReadonlySet<string> =
  new Set([...COLUMN_KEYS, 'groupBy', 'units', 'defaultRoute', 'design']);
/** Accepted unit strings per quantity: the keys of the adapter's own conversion tables. */
const UNIT_TABLES: Readonly<Record<keyof PipelineUnits, Readonly<Record<string, number>>>> = {
  concentration: CONC_TO_MG_L, time: TIME_TO_HOUR, dose: DOSE_TO_MG, bodyWeight: BW_TO_KG,
};

export const isPlainObject = (v: unknown): v is Record<string, unknown> =>
  typeof v === 'object' && v !== null && !Array.isArray(v);

/**
 * Validate a parsed mapping object. Throws (configuration error) on an unknown
 * key, a wrong type, an unknown unit / route / design, or a missing required
 * column. Missing units and route are NOT errors here — see `mappingGaps`.
 */
export function validateMapping(o: unknown): PipelineMapping {
  if (!isPlainObject(o)) throw new Error('mapping: expected a JSON object');
  for (const k of Object.keys(o))
    if (!MAPPING_KEYS.has(k)) throw new Error(`mapping: unknown key '${k}'`);
  for (const k of ['subject', 'time', 'concentration'] as const) {
    if (typeof o[k] !== 'string' || o[k] === '')
      throw new Error(`mapping: '${k}' is required`);
  }
  for (const k of COLUMN_KEYS) {
    if (o[k] !== undefined && typeof o[k] !== 'string')
      throw new Error(`mapping: '${k}' must be a column name`);
  }
  if (o.groupBy !== undefined &&
      (!Array.isArray(o.groupBy) || o.groupBy.some((g) => typeof g !== 'string')))
    throw new Error(`mapping: 'groupBy' must be an array of column names`);
  const units = o.units ?? {};
  if (!isPlainObject(units)) throw new Error(`mapping: 'units' must be an object`);
  for (const [k, v] of Object.entries(units)) {
    const table = UNIT_TABLES[k as keyof PipelineUnits];
    if (table === undefined) throw new Error(`mapping: unknown key 'units.${k}'`);
    if (typeof v !== 'string' || !(v in table)) {
      throw new Error(`mapping: unknown ${k} unit '${String(v)}' — expected one of ` +
        `${Object.keys(table).join(', ')}`);
    }
  }
  if (o.defaultRoute !== undefined &&
      (typeof o.defaultRoute !== 'string' || ROUTE_STRING_TO_CODE[o.defaultRoute] === undefined)) {
    throw new Error(`mapping: unknown defaultRoute '${String(o.defaultRoute)}' — expected one of ` +
      `${Object.keys(ROUTE_STRING_TO_CODE).join(', ')}`);
  }
  if (o.design !== undefined && !DESIGNS.includes(o.design as SamplingDesign))
    throw new Error(`mapping: unknown design '${String(o.design)}' — expected one of ${DESIGNS.join(', ')}`);
  return {...o, groupBy: (o.groupBy as string[] | undefined) ?? [], units} as unknown as PipelineMapping;
}

/** Parse a mapping JSON string; throws with a readable message on bad input. */
export function parseMapping(json: string): PipelineMapping {
  return validateMapping(parseConfigJson('mapping', json));
}

/**
 * What the mapping leaves open that the pipeline will not assume: the
 * concentration and time units always, the dose / body-weight unit when that
 * column is mapped, and the route (a route column or an explicit `defaultRoute`).
 * Non-empty ⇒ every partition is `Blocked`, listing all of them.
 */
export function mappingGaps(m: PipelineMapping): string[] {
  const gaps: string[] = [];
  if (m.units.concentration === undefined) gaps.push('concentration unit');
  if (m.units.time === undefined) gaps.push('time unit');
  if (m.dose !== undefined && m.units.dose === undefined) gaps.push('dose unit');
  if (m.bodyWeight !== undefined && m.units.bodyWeight === undefined) gaps.push('body-weight unit');
  if (m.route === undefined && m.defaultRoute === undefined)
    gaps.push('route (no route column and no defaultRoute)');
  return gaps;
}

/**
 * The mapping's units, for a mapping with no `mappingGaps`. The dose and
 * body-weight units of an UNMAPPED column are inert placeholders: no value is
 * ever read in them (the dose is then 0 and the body weight absent).
 */
export function unitsOf(m: PipelineMapping): PipelineUnits {
  const u = m.units;
  if (u.concentration === undefined || u.time === undefined)
    throw new Error(`unitsOf: mapping leaves ${mappingGaps(m).join(', ')} open`);
  return {concentration: u.concentration, time: u.time, dose: u.dose ?? 'mg', bodyWeight: u.bodyWeight ?? 'kg'};
}

/** Every column the mapping names (roles and groupBy). */
export function mappedColumns(m: PipelineMapping): Set<string> {
  return new Set([...COLUMN_KEYS.map((k) => m[k]), ...m.groupBy].filter((n): n is string => n !== undefined));
}

/** Columns of `df` the mapping names that do not exist. */
export function missingColumns(df: DG.DataFrame, m: PipelineMapping): string[] {
  return [...mappedColumns(m)].filter((n) => df.col(n) === null);
}

export type SexResolution =
  | {readonly ok: true; readonly col: string; readonly map: Record<string, SexCode>}
  | {readonly ok: false; readonly reason: string};

/**
 * The mapped sex column resolved to M / F / U with the app's default assignment
 * (`defaultSexAssignment`), or why it cannot be: no sex column, or values the
 * assignment does not recognise (e.g. coded `1` / `2`) — those are never guessed.
 */
export function sexOf(df: DG.DataFrame, m: PipelineMapping): SexResolution {
  const col = m.sex !== undefined ? df.col(m.sex) : null;
  if (col === null) return {ok: false, reason: NO_SEX_COLUMN};
  const map: Record<string, SexCode> = {};
  const unknown: string[] = [];
  for (const v of uniqueSexValues(col)) {
    const code = defaultSexAssignment(v);
    if (code === null) unknown.push(v); else map[v] = code;
  }
  if (unknown.length > 0) return {ok: false, reason: `sex values not recognised: ${unknown.join(', ')}`};
  return {ok: true, col: col.name, map};
}
