/**
 * Headless NCA for pipelines (`ncaCompute`): one long observations table in, one
 * row per profile out, each with a review verdict and the reasons for it.
 *
 * The math is the app's: every dense profile goes through `runNcaPipeline` and
 * every sparse arm through `computeSparseNca`, and the λz acceptance verdict is
 * the shared `evaluateAcceptance` predicate. What this module adds is the part a
 * pipeline needs and the app gets from a person:
 *
 * 1. **Profile key.** Rows are partitioned by the mapping's `groupBy` columns
 *    (study, compound, period) before the subject-keyed engine sees them, so a
 *    cassette or a crossover is never pooled.
 * 2. **Never pool silently.** A partition in which a subject has two samples at
 *    the same time is split on the one unmapped column that explains the
 *    duplicates (flagged for confirmation), or blocked when no column does.
 * 3. **Sampling design.** Serial (per-animal NCA) or sparse (one composite per
 *    route × dose), declared in the mapping or decided from the samples per animal.
 * 4. **Review by exception.** A profile is routed to review when its λz fails
 *    acceptance, the engine could not finish, a screening heuristic fires (terminal
 *    rise, exposure outlier — both site-tunable through `rules`), the CRO left a
 *    comment matching `rules.reviewCommentPattern`, or the engine raised a unit-error
 *    heuristic. Other comments and engine warnings are information.
 *
 * **Output contract (`RESULTS_BASIS = per-kg-1`, nca-cl-bw-normalised-units AD-6).** `CL` / `VZ` / `VSS`
 * are per kg (mL/min/kg, L/kg) on every row, read from the app's results frame (the pipeline has no
 * normalisation of its own); `CL_ABS` / `VZ_ABS` / `VSS_ABS` are absolute; `BW_KG` is kg; `DOSE_NORM`
 * is the ENGINE dose basis. `CL_ML_MIN_KG` / `VZ_L_KG` are one-release deprecated aliases. The basis
 * rides in-band as the `RESULTS_BASIS` column and the `nca.resultsBasis` frame tag.
 */

import * as DG from 'datagrok-api/dg';
import type {nca} from '@datagrok-libraries/sci-comp';

import {runNcaPipeline} from '../controller/nca-controller';
import {
  DEFAULT_RULES, DEFAULT_ACCEPTANCE_RSQUARED, evaluateAcceptance, normalizeRouteColumn,
  DERIVED_ROUTE_CODE_COL, ROUTE_STRING_TO_CODE, ROUTE_CODE_TO_STRING, computeSparseStrata,
  INSTALLED_ENGINE_VERSION, AuditLog, numAt, ACCEPT_NA, CANONICAL_UNIT_TAGS, isIntravascular,
  scaleProfileInputs, impliedClPerKgWarning, ensureBoolBlq, RESULTS_BASIS, CL_PER_KG_UNIT, VOLUME_PER_KG_UNIT,
} from '../adapter';
import type {
  AcceptanceFlagCode, AcceptanceReading, AcceptanceThresholds, ColumnContract, ProfileScalars,
} from '../adapter';
import {
  detectPipelineMapping, validateMapping, missingColumns, mappingGaps, mappedColumns, unitsOf, sexOf,
  parseConfigJson, isPlainObject, type PipelineMapping, type SexResolution,
} from './pipeline-mapping';

export const VERDICT_ACCEPTED = 'Accepted';
export const VERDICT_REVIEW = 'Review';
export const VERDICT_BLOCKED = 'Blocked';

export const COMPOSITE_SUBJECT = '(composite)';
/** Separator of the human-readable `REASONS` / `INFO_FLAGS` text (the `_JSON` columns hold the arrays). */
export const FLAG_SEPARATOR = ' | ';
/** `SEX` of a composite that pools animals of more than one sex. */
export const MIXED_SEX = 'mixed';
/** Review reason of a profile the engine could not compute (the app shows it as Not estimated). */
export const REASON_NCA_FAILED = 'NCA failed for this profile';
/** Review reason of a profile without a λz fit (the app shows it as Not estimated). */
export const REASON_LZ_NOT_ESTIMABLE = 'λz not estimable';

/** Screening heuristic default: AUClast/dose this many times above or below the arm median. */
const DEFAULT_EXPOSURE_OUTLIER_FOLD = 4;
/** Screening heuristic default: a measurable point after Tmax this many times its predecessor. */
const DEFAULT_TERMINAL_RISE_FOLD = 2;
/**
 * CRO comments that route a profile to review (case-insensitive regex): dosing and
 * sample incidents. Any other comment is information. `rules.reviewCommentPattern`
 * REPLACES this pattern (it does not extend it).
 */
const DEFAULT_REVIEW_COMMENT_PATTERN =
  'regurgitat|vomit|emesis|missed|partial dose|hemoly|haemoly|mislabel|deviation|repeat';

/** Pipeline rule overrides: the engine's `NcaRules`, the acceptance R² and the pipeline's own policy. */
export interface PipelineRules {
  readonly rules: nca.NcaRules;
  readonly acceptance: AcceptanceThresholds;
  /**
   * Acceptance criteria reported as information instead of routing to review — a
   * site policy, e.g. discovery screening accepts a 3-point λz (`NPTS_LE`) and a
   * short span (`SPAN_LT`). Violations still appear in INFO_FLAGS. Default: none,
   * i.e. the app's verdict.
   */
  readonly acceptanceIgnore: ReadonlySet<AcceptanceFlagCode>;
  /** Screening heuristic, site-tunable (`exposureOutlierFold`, > 1). */
  readonly exposureOutlierFold: number;
  /** Screening heuristic, site-tunable (`terminalRiseFold`, > 1). */
  readonly terminalRiseFold: number;
  /** CRO comments matching this route the profile to review (`reviewCommentPattern`). */
  readonly reviewComment: RegExp;
}

/** Numeric output columns copied from the app's results frame, in output order. `CL` / `VZ` / `VSS`
 *  are per kg and `_ABS` absolute, as the app forms them (nca-cl-bw-normalised-units AD-6): the
 *  pipeline has no normalisation of its own. */
const PARAM_COLS = [
  'CMAX', 'TMAX', 'C0_EXTRAP', 'AUC_LAST', 'AUC_INF', 'PCT_EXTRAP', 'LAMBDA_Z', 'HALF_LIFE',
  'CL', 'VZ', 'MRT', 'VSS', 'CL_ABS', 'VZ_ABS', 'VSS_ABS', 'LAMBDAZ_ADJ_R2', 'LAMBDAZ_NPOINTS', 'SPAN_RATIO',
] as const;
type ParamCol = typeof PARAM_COLS[number];

/** Every acceptance criterion code — a record so the compiler flags a code added upstream. */
const ACCEPTANCE_CODES: Readonly<Record<AcceptanceFlagCode, true>> =
  {SPAN_LT: true, ADJR2_LT: true, NPTS_LE: true, EXTRAP_GT: true};

interface OutRow {
  keys: Map<string, string>;
  subject: string;
  /** M / F / U, `mixed` for a composite pooling sexes, '' when unknown or varying. */
  sex: string;
  route: string;
  dose: number;
  /** Body weight in kg — the app's `BW_KG` (canonical); a composite's group mean. */
  bodyWeight: number;
  nSamples: number;
  params: Partial<Record<ParamCol, number>>;
  clBasis: string;
  /** How the ENGINE dose was formed: 'absolute' (mg) or 'per_kg'. Not the CL reporting basis. */
  doseNorm: string;
  status: string;
  verdict: string;
  reasons: string[];
  info: string[];
}

/** Keys of `o` must be keys of `defaults`, with values of the same JSON type. */
function checkKnownKeys(o: Record<string, unknown>, defaults: object, prefix: string): void {
  for (const [k, v] of Object.entries(o)) {
    if (!(k in defaults)) throw new Error(`rules: unknown key '${prefix}${k}'`);
    const expected = typeof (defaults as Record<string, unknown>)[k];
    if (typeof v !== expected) throw new Error(`rules: '${prefix}${k}' must be a ${expected}`);
  }
}

function foldRule(v: unknown, key: string, dflt: number): number {
  if (v === undefined) return dflt;
  if (typeof v !== 'number' || !(v > 1)) throw new Error(`rules: '${key}' must be a number > 1`);
  return v;
}

/**
 * Parse the rules JSON: a partial `NcaRules` plus the pipeline's own keys
 * (`acceptanceRSquared`, `acceptanceIgnore`, `exposureOutlierFold`,
 * `terminalRiseFold`, `reviewCommentPattern`). A configuration error — malformed
 * JSON, an unknown key, a wrong type, an invalid pattern — throws.
 */
export function parseRules(json: string): PipelineRules {
  const o = json.trim() === '' ? {} : parseConfigJson('rules', json);
  if (!isPlainObject(o)) throw new Error('rules: expected a JSON object');
  const {acceptanceRSquared, acceptanceIgnore, exposureOutlierFold, terminalRiseFold, reviewCommentPattern,
    lambdaZ, blq, ...rest} = o;
  checkKnownKeys(rest, DEFAULT_RULES, '');
  const nested = [['lambdaZ', lambdaZ, DEFAULT_RULES.lambdaZ], ['blq', blq, DEFAULT_RULES.blq]] as const;
  for (const [name, part, dflt] of nested) {
    if (part === undefined) continue;
    if (!isPlainObject(part)) throw new Error(`rules: '${name}' must be an object`);
    checkKnownKeys(part, dflt, `${name}.`);
  }
  if (acceptanceRSquared !== undefined &&
      (typeof acceptanceRSquared !== 'number' || acceptanceRSquared < 0 || acceptanceRSquared > 1))
    throw new Error(`rules: 'acceptanceRSquared' must be a number in [0, 1]`);
  if (acceptanceIgnore !== undefined && (!Array.isArray(acceptanceIgnore) ||
      acceptanceIgnore.some((c) => typeof c !== 'string' || !(c in ACCEPTANCE_CODES)))) {
    throw new Error(`rules: 'acceptanceIgnore' must be an array of ` +
      `${Object.keys(ACCEPTANCE_CODES).join(', ')}`);
  }
  let reviewComment: RegExp;
  if (reviewCommentPattern !== undefined && typeof reviewCommentPattern !== 'string')
    throw new Error(`rules: 'reviewCommentPattern' must be a string`);
  try {
    reviewComment = new RegExp((reviewCommentPattern as string | undefined) ?? DEFAULT_REVIEW_COMMENT_PATTERN, 'i');
  } catch (e) {
    throw new Error(`rules: invalid reviewCommentPattern — ${e instanceof Error ? e.message : String(e)}`);
  }

  const rules: nca.NcaRules = {
    ...DEFAULT_RULES,
    ...(rest as Partial<nca.NcaRules>),
    lambdaZ: {...DEFAULT_RULES.lambdaZ, ...(lambdaZ as Partial<nca.LambdaZStrategy> | undefined)},
    blq: {...DEFAULT_RULES.blq, ...(blq as Partial<nca.BlqStrategy> | undefined)},
  };
  return {
    rules,
    acceptance: {
      lambdaAcceptanceRSquared: (acceptanceRSquared as number | undefined) ?? DEFAULT_ACCEPTANCE_RSQUARED,
      lambdaMinPoints: rules.lambdaZ.minPoints,
      extrapWarnPct: rules.extrapWarnPct,
    },
    acceptanceIgnore: new Set((acceptanceIgnore ?? []) as AcceptanceFlagCode[]),
    exposureOutlierFold: foldRule(exposureOutlierFold, 'exposureOutlierFold', DEFAULT_EXPOSURE_OUTLIER_FOLD),
    terminalRiseFold: foldRule(terminalRiseFold, 'terminalRiseFold', DEFAULT_TERMINAL_RISE_FOLD),
    reviewComment,
  };
}

/** Split an acceptance reading into review reasons and informational flags per the policy. */
function splitAcceptance(acc: AcceptanceReading, pr: PipelineRules): {reasons: string[]; info: string[]} {
  const reasons: string[] = [];
  const info: string[] = [];
  acc.flagCodes.forEach((code, i) => (pr.acceptanceIgnore.has(code) ? info : reasons).push(acc.reasons[i]));
  return {reasons, info};
}

/**
 * Engine warnings that route a profile to review: each is a probable unit error
 * (concentration magnitude, body weight read as kg, implied CL outside the
 * physiological band, an infusion duration in the wrong time unit), not a PK
 * finding. User decisions 2026-10-14 (the infusion code at the review gate).
 */
const REVIEW_WARNING_CODES: ReadonlySet<string> = new Set([
  'CONC_MAGNITUDE_HEURISTIC_WARNING', 'BW_UNIT_HEURISTIC_WARNING', 'CL_IMPLAUSIBLE_HEURISTIC_WARNING',
  'INFUSION_DURATION_HEURISTIC_WARNING',
]);
/**
 * Engine warnings left out of the results: the first two restate an acceptance
 * criterion the row already reports (EXTRAP_GT / NPTS_LE); the bootstrap notes
 * describe confidence intervals the pipeline does not output; a topology mismatch
 * is already a `sampling design:` review reason. User decision 2026-10-14.
 */
const SKIPPED_WARNING_CODES: ReadonlySet<string> = new Set([
  'AUC_EXTRAP_HIGH', 'LAMBDAZ_FEW_POINTS', 'SPARSE_BOOTSTRAP_AUC_EXPLORATORY', 'SPARSE_BOOTSTRAP_N_TOO_SMALL',
  'SPARSE_TOPOLOGY_MISMATCH',
]);

/** Engine warnings as `CODE: message`, once per code: unit-error heuristics → reasons, the rest → info. */
function splitEngineWarnings(warnings: ReadonlyArray<{code: string; message: string}>):
  {reasons: string[]; info: string[]} {
  const reasons: string[] = [];
  const info: string[] = [];
  const seen = new Set<string>();
  for (const w of warnings) {
    if (seen.has(w.code) || SKIPPED_WARNING_CODES.has(w.code)) continue;
    seen.add(w.code);
    (REVIEW_WARNING_CODES.has(w.code) ? reasons : info).push(`${w.code}: ${w.message}`);
  }
  return {reasons, info};
}

/**
 * The adapter's unit-error guards on a composite's RAW inputs (the dense path gets
 * them from the extractor per profile). `scaleProfileInputs` owns the concentration
 * magnitude and body-weight checks and `impliedClPerKgWarning` the implied-CL check,
 * so the composite reuses them, scaled values discarded. Only the unit-error codes
 * are kept: the composite reports a missing body weight itself (`PER_KG_NOT_FORMED`).
 */
function compositeUnitGuards(part: DG.DataFrame, m: PipelineMapping, g: {meanBodyWeightKg: number | null;
  engineCl: number; engineDoseNorm: 'absolute' | 'per_kg'; enginePerKgBodyWeightKg: number | null}):
  nca.ParameterWarning[] {
  const u = unitsOf(m);
  const raw = (name: string): Float64Array => {
    const c = part.col(name)!;
    return Float64Array.from({length: part.rowCount}, (_, i) => numAt(c, i));
  };
  let bw: number | null = null;
  if (m.bodyWeight !== undefined) {
    const v = [...raw(m.bodyWeight)].filter((x) => Number.isFinite(x));
    if (v.length > 0) bw = v.reduce((a, b) => a + b, 0) / v.length;
  }
  const scaled = scaleProfileInputs(
    {timeUnits: u.time, concentrationUnits: u.concentration, doseUnits: u.dose},
    raw(m.nominalTime ?? m.time), raw(m.concentration), 0, 0, bw, null, u.bodyWeight);
  // The ENGINE-basis CL with its engine dose basis (architect N1): a `/kg` composite runs per kg
  // (read as is, never ÷ BW again), an absolute one is ÷ its per-kg divisor. Never the reported
  // mL/min/kg value (a ~16.7× band error).
  const cl = impliedClPerKgWarning(g.engineCl, g.enginePerKgBodyWeightKg, g.engineDoseNorm);
  return [...scaled.warnings, ...(cl !== null ? [cl] : [])].filter((w) => REVIEW_WARNING_CODES.has(w.code));
}

/** A dense result row's `WARNINGS_JSON` (the writer's serialized `ParameterWarning[]`). */
function warningsAt(res: DG.DataFrame, i: number): Array<{code: string; message: string}> {
  const json = strAt(res.col('WARNINGS_JSON'), i);
  return json === '' ? [] : JSON.parse(json) as Array<{code: string; message: string}>;
}

/** CRO comments as review reasons (matching the review pattern) or information (all others). */
function splitComments(comments: ReadonlyArray<string>, pr: PipelineRules): {reasons: string[]; info: string[]} {
  const reasons: string[] = [];
  const info: string[] = [];
  for (const c of comments) (pr.reviewComment.test(c) ? reasons : info).push(`CRO comment: "${c}"`);
  return {reasons, info};
}

export const strAt = (col: DG.Column | null, row: number): string =>
  col === null || col.isNone(row) ? '' : String(col.get(row));

/** Group row indices by a key function, preserving first-seen order. */
export function groupRows(rows: ReadonlyArray<number>, keyOf: (r: number) => string): Map<string, number[]> {
  const out = new Map<string, number[]>();
  for (const r of rows) {
    const k = keyOf(r);
    let list = out.get(k);
    if (list === undefined) out.set(k, list = []);
    list.push(r);
  }
  return out;
}

const allRowsOf = (df: DG.DataFrame): number[] => Array.from({length: df.rowCount}, (_, i) => i);

/** (subject, time) pairs that occur more than once among `rows`. */
export function duplicateSamples(df: DG.DataFrame, m: PipelineMapping, rows: ReadonlyArray<number>): number {
  const subj = df.col(m.subject)!;
  const time = df.col(m.time)!;
  const seen = new Set<string>();
  let dups = 0;
  for (const r of rows) {
    const k = `${strAt(subj, r)}\u0001${numAt(time, r)}`;
    if (seen.has(k)) dups++; else seen.add(k);
  }
  return dups;
}

/**
 * The unmapped categorical column that, used as an extra key, removes every
 * duplicate sample in `rows` — the cassette's compound or the crossover's period
 * that the mapping missed. `null` when no single column does.
 */
function findSplitColumn(df: DG.DataFrame, m: PipelineMapping, rows: ReadonlyArray<number>,
  exclude: ReadonlySet<string>): string | null {
  for (const col of df.columns.toList()) {
    if (exclude.has(col.name)) continue;
    if (col.type !== DG.COLUMN_TYPE.STRING && col.type !== DG.COLUMN_TYPE.INT) continue;
    const parts = groupRows(rows, (r) => strAt(col, r));
    if (parts.size < 2 || parts.size > 50) continue;
    let clean = true;
    for (const part of parts.values())
      if (duplicateSamples(df, m, part) > 0) {clean = false; break;}
    if (clean) return col.name;
  }
  return null;
}

/** Numeric, concentration-like columns other than the mapped one that hold data in `rows`. */
function extraConcColumns(df: DG.DataFrame, m: PipelineMapping, rows: ReadonlyArray<number>): string[] {
  const out: string[] = [];
  for (const col of df.columns.toList()) {
    if (col.name === m.concentration || col.name === m.lloq) continue;
    if (col.type !== DG.COLUMN_TYPE.FLOAT && col.type !== DG.COLUMN_TYPE.INT) continue;
    const n = col.name.toLowerCase().replace(/[\s_-]+/g, '');
    if (!/(conc|concentration)/.test(n) || /lloq|lod/.test(n)) continue;
    if (rows.some((r) => !col.isNone(r))) out.push(col.name);
  }
  return out;
}

/** Copy `rows` of `df` into a new frame (row order preserved). */
export function subFrame(df: DG.DataFrame, rows: ReadonlyArray<number>): DG.DataFrame {
  const mask = DG.BitSet.create(df.rowCount);
  for (const r of rows) mask.set(r, true, false);
  return df.clone(mask);
}

// `ensureBoolBlq` moved to `adapter/adapter-utils.ts` (BUG-19 AD-8): the app's sparse path reads
// the BLQ flag the same way. `addScaled` moved into `adapter/sparse-nca.ts` `computeSparseStrata`.

type RouteRead =
  | {readonly ok: false; readonly reason: string}
  | {readonly ok: true; readonly codes: Int32Array | null; readonly warnings?: Map<string, nca.ParameterWarning[]>};

/**
 * Per-row route codes of `sub` through the app's tolerant normaliser (`null`
 * without a route column; `defaultRoute` then applies). An empty route cell takes
 * `defaultRoute`; with no `defaultRoute` it blocks — the pipeline never assumes a route.
 */
function readRoutes(sub: DG.DataFrame, m: PipelineMapping): RouteRead {
  if (m.route === undefined) return {ok: true, codes: null};
  const subj = sub.col(m.subject)!;
  const nr = normalizeRouteColumn(sub.col(m.route)!, subj,
    {infusionDurationCol: m.infusionDuration !== undefined ? sub.col(m.infusionDuration) : null});
  if (nr.errors.length > 0)
    return {ok: false, reason: `route column: ${nr.errors.map((e) => e.message).join('; ')}`};
  if (m.defaultRoute === undefined) {
    const empty = new Set(allRowsOf(sub).filter((i) => nr.codes[i] === DG.INT_NULL).map((i) => strAt(subj, i)));
    if (empty.size > 0) {
      return {ok: false, reason: `route column: '${m.route}' is empty for ${[...empty].join(', ')} ` +
        `and the mapping sets no defaultRoute`};
    }
  }
  return {ok: true, codes: nr.codes, warnings: nr.warningsBySubject};
}

/**
 * Per-row route string of `sub` (the app's normalised code as `ROUTE_CODE_TO_STRING`;
 * `defaultRoute` where no route is read) — the same reading the sparse path splits its composites
 * by, so `ncaSaveToStudy` can save one sparse analysis per route (nca-pipeline-save-rework F6).
 */
export function rowRoutes(sub: DG.DataFrame, m: PipelineMapping):
  {readonly ok: true; readonly routes: string[]} | {readonly ok: false; readonly reason: string} {
  const r = readRoutes(sub, m);
  if (!r.ok) return r;
  return {ok: true, routes: allRowsOf(sub).map((i) => ROUTE_CODE_TO_STRING[
    r.codes !== null && r.codes[i] !== DG.INT_NULL ? r.codes[i] : ROUTE_STRING_TO_CODE[m.defaultRoute!]])};
}

function contractFor(sub: DG.DataFrame, m: PipelineMapping, blqName: string | null,
  routeName: string | null): ColumnContract {
  const u = unitsOf(m);
  const subjectCol = sub.col(m.subject)!;
  const opt = <T>(name: string | null | undefined, spec: T): T | undefined =>
    name != null ? spec : undefined;
  const obs: Record<string, unknown> = {
    subjectId: {name: m.subject, required: true,
      rawType: subjectCol.type === DG.COLUMN_TYPE.STRING ? 'string' : 'int32'},
    time: {name: m.time, required: true, rawType: 'float64', units: u.time},
    concentration: {name: m.concentration, required: true, rawType: 'float64', units: u.concentration},
    dose: opt(m.dose, {name: m.dose, required: false, rawType: 'float64', units: u.dose}),
    route: opt(routeName, {name: routeName, required: false, rawType: 'int32'}),
    bodyWeight: opt(m.bodyWeight, {name: m.bodyWeight, required: false, rawType: 'float64'}),
    blqFlag: opt(blqName, {name: blqName, required: false, rawType: 'uint8'}),
    lloq: opt(m.lloq, {name: m.lloq, required: false, rawType: 'float64'}),
    // Read in `units.time`: the extractor rescales the duration with the time factor.
    infusionDuration: opt(m.infusionDuration,
      {name: m.infusionDuration, required: false, rawType: 'float64'}),
  };
  for (const k of Object.keys(obs)) if (obs[k] === undefined) delete obs[k];
  return {observations: obs} as unknown as ColumnContract;
}

/** Concentration–time points of one subject in `sub`, measurable only, sorted by time. */
function measurablePoints(sub: DG.DataFrame, m: PipelineMapping, blqName: string | null,
  subject: string): Array<{t: number; c: number}> {
  const s = sub.col(m.subject)!;
  const t = sub.col(m.time)!;
  const c = sub.col(m.concentration)!;
  const b = blqName !== null ? sub.col(blqName) : null;
  const pts: Array<{t: number; c: number}> = [];
  for (let i = 0; i < sub.rowCount; i++) {
    if (strAt(s, i) !== subject) continue;
    if (b !== null && b.get(i) === true) continue;
    const ci = numAt(c, i);
    const ti = numAt(t, i);
    if (Number.isFinite(ci) && Number.isFinite(ti) && ci > 0) pts.push({t: ti, c: ci});
  }
  return pts.sort((a, b2) => a.t - b2.t);
}

/** Screening heuristic: a measurable point after Tmax more than `fold` × its predecessor. */
function terminalRise(pts: ReadonlyArray<{t: number; c: number}>, unitC: string, unitT: string,
  fold: number): string | null {
  if (pts.length < 3) return null;
  let iMax = 0;
  for (let i = 1; i < pts.length; i++) if (pts[i].c > pts[iMax].c) iMax = i;
  for (let i = iMax + 2; i < pts.length; i++) {
    if (pts[i].c > fold * pts[i - 1].c) {
      return `screening heuristic: terminal rise — at ${fmt(pts[i].t)} ${unitT} ` +
        `(${fmt(pts[i].c)} vs ${fmt(pts[i - 1].c)} ${unitC} at ${fmt(pts[i - 1].t)} ${unitT})`;
    }
  }
  return null;
}

function fmt(v: number): string {
  if (!Number.isFinite(v)) return 'NA';
  return Math.abs(v) >= 100 ? v.toFixed(0) : Number(v.toPrecision(3)).toString();
}

/** Distinct non-empty comment texts in `rows` of `sub`. */
function commentsIn(sub: DG.DataFrame, m: PipelineMapping, rows: ReadonlyArray<number>): string[] {
  if (m.comment === undefined) return [];
  const c = sub.col(m.comment)!;
  const out = new Set<string>();
  for (const i of rows) if (strAt(c, i).trim() !== '') out.add(strAt(c, i).trim());
  return [...out];
}

/** Distinct resolved sex codes (M / F / U) in `rows` of `sub`; empty when sex is not resolved. */
function sexCodesIn(sub: DG.DataFrame, sex: SexResolution, rows: ReadonlyArray<number>): Set<string> {
  const out = new Set<string>();
  if (!sex.ok) return out;
  const c = sub.col(sex.col)!;
  for (const i of rows) out.add(sex.map[c.isNone(i) ? '' : c.getString(i)]);
  return out;
}

/** First non-null value of `colName` in `rows` of `sub` (NaN if none). */
function firstNum(sub: DG.DataFrame, colName: string | undefined, rows: ReadonlyArray<number>): number {
  if (colName === undefined) return NaN;
  const c = sub.col(colName)!;
  for (const i of rows) if (!c.isNone(i)) return numAt(c, i);
  return NaN;
}

/** Row indices of `sub` per subject id. */
const rowsBySubject = (sub: DG.DataFrame, m: PipelineMapping): Map<string, number[]> =>
  groupRows(allRowsOf(sub), (r) => strAt(sub.col(m.subject), r));

/** One partition run through the app's own pipeline (contract + RouteCode column + runNcaPipeline). */
export type AppRun =
  | {readonly ok: false; readonly reason: string}
  | {readonly ok: true; readonly contract: ColumnContract;
     readonly outcome: Extract<ReturnType<typeof runNcaPipeline>, {status: 'completed'}>};

/** The inputs of the app's dense run for one partition (everything `runNcaPipeline` takes but the
 *  frame and the rules). Split out of {@link runAppPipeline} (nca-pipeline-save-rework F1) so the
 *  save can drive an `NcaController` on the SAME inputs and read the app's own review state. */
export type AppRunInputs =
  | {readonly ok: false; readonly reason: string}
  | {readonly ok: true; readonly contract: ColumnContract; readonly scalars: ProfileScalars;
     readonly routeWarnings?: Map<string, nca.ParameterWarning[]>};

/**
 * Run the app's dense pipeline on `sub` (adds the derived RouteCode column to it).
 * Used by `ncaCompute`; the save drives an `NcaController` on {@link prepareAppRun}'s inputs, so
 * a saved analysis's summary is the same compute the pipeline judged. `m` must have no
 * `mappingGaps`.
 */
export function runAppPipeline(sub: DG.DataFrame, m: PipelineMapping, pr: PipelineRules,
  blqName: string | null): AppRun {
  const prep = prepareAppRun(sub, m, blqName);
  if (!prep.ok) return prep;
  const outcome = runNcaPipeline(sub, prep.contract, prep.scalars, pr.rules,
    {log: new AuditLog(), ...(prep.routeWarnings !== undefined ? {routeWarningsBySubject: prep.routeWarnings} : {})});
  if (outcome.status === 'validation-failed')
    return {ok: false, reason: `validation: ${outcome.errors.map((e) => e.message).join('; ')}`};
  return {ok: true, contract: prep.contract, outcome};
}

/**
 * The app's dense-run inputs for `sub`: the column contract, the profile scalars and the
 * per-subject route warnings. Adds the derived RouteCode column to `sub` when the mapping has a
 * route column. Blocks (`ok: false`) when the route column cannot be read. `m` must have no
 * `mappingGaps`.
 */
export function prepareAppRun(sub: DG.DataFrame, m: PipelineMapping, blqName: string | null): AppRunInputs {
  const u = unitsOf(m);
  const routes = readRoutes(sub, m);
  if (!routes.ok) return routes;
  let routeName: string | null = null;
  if (routes.codes !== null) {
    if (sub.col(DERIVED_ROUTE_CODE_COL) !== null) sub.columns.remove(DERIVED_ROUTE_CODE_COL);
    sub.columns.add(DG.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL, routes.codes));
    routeName = DERIVED_ROUTE_CODE_COL;
  }
  const scalars: ProfileScalars = {
    dose: 0,
    // The scalar route is read only for a row without a column route: `defaultRoute`
    // by construction (`readRoutes` blocks an empty cell without one), so the
    // 'other' placeholder is never read.
    route: ROUTE_STRING_TO_CODE[m.defaultRoute ?? 'other'],
    doseUnits: u.dose, concentrationUnits: u.concentration, timeUnits: u.time,
    bodyWeight: null, bodyWeightUnit: u.bodyWeight, infusionDuration: null, lloq: 0,
  };
  const contract = contractFor(sub, m, blqName, routeName);
  return {ok: true, contract, scalars, ...(routes.warnings !== undefined ? {routeWarnings: routes.warnings} : {})};
}

/** Dense NCA on one partition. */
function computeDense(sub: DG.DataFrame, m: PipelineMapping, pr: PipelineRules, sex: SexResolution,
  keys: Map<string, string>, blqName: string | null, extraReasons: string[]): OutRow[] {
  const u = unitsOf(m);
  const run = runAppPipeline(sub, m, pr, blqName);
  if (!run.ok) return [blockedRow(keys, run.reason)];
  const res = run.outcome.results;
  const bySubject = rowsBySubject(sub, m);
  const rows: OutRow[] = [];
  for (let i = 0; i < res.rowCount; i++) {
    const subject = strAt(res.col('SUBJECT'), i);
    const subjectRows = bySubject.get(subject) ?? [];
    const params: Partial<Record<ParamCol, number>> = {};
    for (const p of PARAM_COLS) {
      const c = res.col(p);
      params[p] = c === null ? NaN : numAt(c, i);
    }
    const reasons = [...extraReasons];
    const status = strAt(res.col('STATUS'), i);
    if (status === 'failed') reasons.push(REASON_NCA_FAILED);
    const acc = evaluateAcceptance({
      spanRatio: params.SPAN_RATIO!, adjR2: params.LAMBDAZ_ADJ_R2!,
      nPoints: params.LAMBDAZ_NPOINTS!, pctExtrap: params.PCT_EXTRAP!,
    }, pr.acceptance);
    const split = splitAcceptance(acc, pr);
    reasons.push(...split.reasons);
    if (acc.verdict === ACCEPT_NA && status !== 'failed') reasons.push(REASON_LZ_NOT_ESTIMABLE);
    const pts = measurablePoints(sub, m, blqName, subject);
    const rise = terminalRise(pts, u.concentration, u.time, pr.terminalRiseFold);
    if (rise !== null) reasons.push(rise);
    const comments = splitComments(commentsIn(sub, m, subjectRows), pr);
    reasons.push(...comments.reasons);
    const warnings = splitEngineWarnings(warningsAt(res, i));
    reasons.push(...warnings.reasons);
    const sexCodes = sexCodesIn(sub, sex, subjectRows);
    rows.push({
      keys, subject, sex: sexCodes.size === 1 ? [...sexCodes][0] : '', route: strAt(res.col('ROUTE'), i),
      // The app's canonical-kg weight (E4: the raw column read gave grams on a `g` mapping).
      dose: firstNum(sub, m.dose, subjectRows), bodyWeight: numAt(res.col('BW_KG')!, i),
      nSamples: pts.length, params, clBasis: strAt(res.col('CL_BASIS'), i),
      doseNorm: strAt(res.col('DOSE_NORM'), i), status,
      verdict: reasons.length > 0 ? VERDICT_REVIEW : VERDICT_ACCEPTED, reasons,
      info: [...split.info, ...comments.info, ...warnings.info],
    });
  }
  return rows;
}

/**
 * Composite (sparse) NCA on one partition, through the adapter's `computeSparseStrata` — the SAME
 * canonicalization, route × dose strata and infusion refusal the app runs (BUG-19 AD-2): one
 * composite per route × dose cell. Sexes are pooled within a cell, as in the app and the PKNCA
 * mouse reference; the composite's `SEX` says so. A cell's comments, sexes, dose and unit guards
 * are read on `part`, the cell's rows of `sub`.
 */
function computeSparse(sub: DG.DataFrame, m: PipelineMapping, pr: PipelineRules, sex: SexResolution,
  keys: Map<string, string>, blqName: string | null, extraReasons: string[],
  declaredTopology: nca.SamplingTopology | undefined): OutRow[] {
  const u = unitsOf(m);
  const routes = readRoutes(sub, m);
  if (!routes.ok) return [blockedRow(keys, routes.reason)];
  if (routes.codes !== null) {
    if (sub.col(DERIVED_ROUTE_CODE_COL) !== null) sub.columns.remove(DERIVED_ROUTE_CODE_COL);
    sub.columns.add(DG.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL, routes.codes));
  }
  const strata = computeSparseStrata(sub, {
    concCol: m.concentration, nominalTimeCol: m.nominalTime ?? m.time, animalIdCol: m.subject,
    blqCol: blqName, doseCol: m.dose ?? null, bodyWeightCol: m.bodyWeight ?? null,
    lloq: {colName: m.lloq ?? null, scalar: 0},
    // A blank route cell takes `defaultRoute`; `readRoutes` blocks one when there is none, so the
    // 'other' placeholder is never read (as in `prepareAppRun`).
    route: {codesCol: routes.codes !== null ? DERIVED_ROUTE_CODE_COL : null, scalar: m.defaultRoute ?? 'other'},
    units: {time: u.time, concentration: u.concentration, dose: u.dose, bodyWeight: u.bodyWeight},
    ...(declaredTopology !== undefined ? {declaredTopology} : {}),
    acceptance: pr.acceptance,
    bootstrapIterations: 200,
  });
  if (!strata.ok) return [blockedRow(keys, strata.reason)];

  const out: OutRow[] = [];
  for (const b of strata.blocked) out.push({...blockedRow(keys, b.reason), route: b.route});
  for (const {group: g, route, rows} of strata.cells) {
    const part = subFrame(sub, rows);
    const allRows = allRowsOf(part);
    const split = splitAcceptance(g.acceptance, pr);
    const reasons = [...extraReasons, ...split.reasons];
    if (g.acceptance.verdict === ACCEPT_NA) reasons.push(REASON_LZ_NOT_ESTIMABLE);
    for (const w of g.warnings)
      if (w.code === 'SPARSE_TOPOLOGY_MISMATCH') reasons.push(`sampling design: ${w.message}`);
    const comments = splitComments(commentsIn(part, m, allRows), pr);
    reasons.push(...comments.reasons);
    const warnings = splitEngineWarnings([...g.warnings, ...compositeUnitGuards(part, m, g)]);
    reasons.push(...warnings.reasons);
    const sexCodes = sexCodesIn(part, sex, allRows);
    const params: Partial<Record<ParamCol, number>> = {
      CMAX: g.cmax, TMAX: g.tmax, C0_EXTRAP: NaN, AUC_LAST: g.aucLast, AUC_INF: g.aucInf,
      PCT_EXTRAP: g.pctExtrap, LAMBDA_Z: g.lambdaZ, HALF_LIFE: g.halfLife, CL: g.cl, VZ: g.vz,
      MRT: g.mrt, VSS: g.vss, CL_ABS: g.clAbs, VZ_ABS: g.vzAbs, VSS_ABS: g.vssAbs,
      LAMBDAZ_ADJ_R2: g.lambdaZAdjR2, LAMBDAZ_NPOINTS: g.lambdaZNPoints, SPAN_RATIO: g.spanRatio,
    };
    out.push({
      keys, subject: COMPOSITE_SUBJECT,
      sex: sexCodes.size === 1 ? [...sexCodes][0] : (sexCodes.size > 1 ? MIXED_SEX : ''),
      route, dose: firstNum(part, m.dose, allRows), bodyWeight: g.meanBodyWeightKg ?? NaN,
      nSamples: g.rawPoints.length, params, clBasis: isIntravascular(route) ? 'CL' : 'CL/F',
      // The ENGINE dose basis (architect N2/S4): `per_kg` for a `/kg` composite (the engine runs
      // on D_kg since AD-5), `absolute` for an absolute-mg one. CL is per kg either way.
      doseNorm: g.engineDoseNorm, status: 'composite',
      verdict: reasons.length > 0 ? VERDICT_REVIEW : VERDICT_ACCEPTED, reasons,
      info: [...split.info, ...comments.info, ...warnings.info],
    });
  }
  return out;
}

/** Fewest single-sample animals an auto-detected destructive design is pooled from (user decision,
 *  review gate 2026-10-14): fewer gives a 1–2 point curve, so the partition is blocked instead. */
const MIN_COMPOSITE_ANIMALS = 3;

type Design =
  | {readonly kind: 'serial'}
  | {readonly kind: 'sparse'; readonly topology?: nca.SamplingTopology}
  | {readonly kind: 'blocked'; readonly reason: string};

/**
 * The partition's sampling design. Declared in the mapping, or (`auto`) decided
 * from the measurable-or-BLQ samples per animal:
 * - every animal 1 sample → sparse (the engine derives destructive), blocked below
 *   `MIN_COMPOSITE_ANIMALS` (a 1–2 point "composite" is not a profile);
 * - more than half of the animals have `minPoints + 1` samples or more → serial: the
 *   short ones are dropouts, computed per animal (λz "not estimable" → review);
 * - every animal fewer than `minPoints + 1` samples: all on one nominal-time set →
 *   serial (a short common schedule stays per-animal rather than a composite that
 *   hides inter-animal variability), else sparse (batch);
 * - anything else → blocked.
 */
function samplingDesign(sub: DG.DataFrame, m: PipelineMapping, blqName: string | null,
  minPoints: number): Design {
  switch (m.design ?? 'auto') {
  case 'serial': return {kind: 'serial'};
  case 'sparse-destructive': return {kind: 'sparse', topology: 'destructive'};
  case 'sparse-batch': return {kind: 'sparse', topology: 'batch'};
  }
  const conc = sub.col(m.concentration)!;
  const blq = blqName !== null ? sub.col(blqName) : null;
  const nom = sub.col(m.nominalTime ?? m.time)!;
  const schedules: string[] = [];
  const counts: number[] = [];
  for (const rows of rowsBySubject(sub, m).values()) {
    const sampled = rows.filter((i) => !conc.isNone(i) || (blq !== null && blq.get(i) === true));
    if (sampled.length === 0) continue;
    counts.push(sampled.length);
    schedules.push([...new Set(sampled.map((i) => numAt(nom, i)))].sort((a, b) => a - b).join(','));
  }
  const dense = minPoints + 1;
  // No sample anywhere: the dense path reports each profile's failure.
  if (counts.length === 0) return {kind: 'serial'};
  if (counts.every((n) => n === 1)) {
    return counts.length >= MIN_COMPOSITE_ANIMALS ? {kind: 'sparse'} : {kind: 'blocked',
      reason: `sparse design with fewer than ${MIN_COMPOSITE_ANIMALS} animals — a composite needs ` +
        `≥ ${MIN_COMPOSITE_ANIMALS} (${counts.length} animal(s), one sample each)`};
  }
  // Strictly more than half: a 50/50 split of complete and short animals is a mixed design.
  if (2 * counts.filter((n) => n >= dense).length > counts.length) return {kind: 'serial'};
  if (counts.every((n) => n < dense))
    return new Set(schedules).size === 1 ? {kind: 'serial'} : {kind: 'sparse'};
  return {kind: 'blocked', reason: `mixed sampling design — declare 'design' in the mapping ` +
    `(samples per animal range ${Math.min(...counts)}–${Math.max(...counts)})`};
}

function blockedRow(keys: Map<string, string>, reason: string): OutRow {
  return {
    keys, subject: '', sex: '', route: '', dose: NaN, bodyWeight: NaN, nSamples: 0, params: {},
    clBasis: '', doseNorm: '', status: '', verdict: VERDICT_BLOCKED, reasons: [reason], info: [],
  };
}

/**
 * Screening heuristic: flag profiles whose AUClast/dose is `fold`× off their arm's
 * median. The arm is the summary's (keys × SEX × route × dose), so a sex difference
 * in exposure is not read as an outlier; without a sex column `sex` is '' throughout.
 */
function flagExposureOutliers(rows: OutRow[], fold: number): void {
  const arms = new Map<string, OutRow[]>();
  for (const r of rows) {
    if (r.verdict === VERDICT_BLOCKED || r.subject === COMPOSITE_SUBJECT) continue;
    const k = `${[...r.keys.values()].join('\u0001')}\u0001${r.sex}\u0001${r.route}\u0001${r.dose}`;
    let list = arms.get(k);
    if (list === undefined) arms.set(k, list = []);
    list.push(r);
  }
  for (const arm of arms.values()) {
    const vals = arm.map((r) => (r.params.AUC_LAST ?? NaN) / r.dose).filter((v) => Number.isFinite(v) && v > 0);
    if (vals.length < 3) continue;
    const sorted = [...vals].sort((a, b) => a - b);
    const mid = sorted.length >> 1;
    const median = sorted.length % 2 === 1 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
    for (const r of arm) {
      const v = (r.params.AUC_LAST ?? NaN) / r.dose;
      if (!Number.isFinite(v) || v <= 0) continue;
      const f = v > median ? v / median : median / v;
      if (f >= fold) {
        r.reasons.push(`screening heuristic: exposure outlier — AUClast/dose ${fmt(f)}x ` +
          `${v > median ? 'above' : 'below'} the arm median (n=${vals.length})`);
        r.verdict = VERDICT_REVIEW;
      }
    }
  }
}

/** Output column order for the key columns: mapping groupBy first, then auto-split columns. */
function keyColumnNames(rows: ReadonlyArray<OutRow>, groupBy: ReadonlyArray<string>): string[] {
  const names = [...groupBy];
  for (const r of rows) for (const k of r.keys.keys()) if (!names.includes(k)) names.push(k);
  return names;
}

/**
 * Results frame. `SEX` is present only when the mapping's sex column resolved
 * (every value recognised); otherwise the frame carries why in the
 * `nca.sexReason` tag, which `ncaSummarize` reports as the F% reason.
 */
function toFrame(rows: ReadonlyArray<OutRow>, m: PipelineMapping | null, sex: SexResolution | null,
  rulesJson?: string): DG.DataFrame {
  const keyNames = keyColumnNames(rows, m?.groupBy ?? []);
  const n = rows.length;
  const cols: DG.Column[] = keyNames.map((k) =>
    DG.Column.fromStrings(k, rows.map((r) => r.keys.get(k) ?? '')));
  const f64 = (name: string, get: (r: OutRow) => number): DG.Column => {
    const a = new Float64Array(n);
    rows.forEach((r, i) => {const v = get(r); a[i] = Number.isFinite(v) ? v : DG.FLOAT_NULL;});
    return DG.Column.fromFloat64Array(name, a);
  };
  cols.push(DG.Column.fromStrings('SUBJECT', rows.map((r) => r.subject)));
  if (sex?.ok === true) cols.push(DG.Column.fromStrings('SEX', rows.map((r) => r.sex)));
  cols.push(
    DG.Column.fromStrings('ROUTE', rows.map((r) => r.route)),
    f64('DOSE', (r) => r.dose),
    f64('BW_KG', (r) => r.bodyWeight),
    DG.Column.fromInt32Array('N_SAMPLES', Int32Array.from(rows.map((r) => r.nSamples))),
    ...PARAM_COLS.map((p) => f64(p, (r) => r.params[p] ?? NaN)),
    // Deprecated aliases, kept ONE release: `CL` / `VZ` are per kg
    // now, and these explicit-unit names are what a unit-blind downstream join can still trust
    // during the transition. Same values; removed in the next release.
    f64('CL_ML_MIN_KG', (r) => r.params.CL ?? NaN),
    f64('VZ_L_KG', (r) => r.params.VZ ?? NaN),
    DG.Column.fromStrings('DOSE_NORM', rows.map((r) => r.doseNorm)),
    DG.Column.fromStrings('CL_BASIS', rows.map((r) => r.clBasis)),
    DG.Column.fromStrings('STATUS', rows.map((r) => r.status)),
    DG.Column.fromStrings('VERDICT', rows.map((r) => r.verdict)),
    // Text for people; the `_JSON` twins (string arrays) for machines — an engine
    // message can contain any separator, so the text is never the parse source.
    DG.Column.fromStrings('REASONS', rows.map((r) => r.reasons.join(FLAG_SEPARATOR))),
    DG.Column.fromStrings('INFO_FLAGS', rows.map((r) => r.info.join(FLAG_SEPARATOR))),
    DG.Column.fromStrings('REASONS_JSON', rows.map((r) => JSON.stringify(r.reasons))),
    DG.Column.fromStrings('INFO_FLAGS_JSON', rows.map((r) => JSON.stringify(r.info))),
    DG.Column.fromStrings('ENGINE_VERSION', rows.map(() => INSTALLED_ENGINE_VERSION)),
    // In-band reporting-basis marker (plan-review F-2): constant per row so it survives a CSV
    // export or a join that drops the frame tag.
    DG.Column.fromStrings('RESULTS_BASIS', rows.map(() => RESULTS_BASIS)),
  );
  const df = DG.DataFrame.fromColumns(cols);
  df.name = 'NCA results';
  // Outputs are in canonical units whatever the input units (mg/L, h), and CL / VZ / VSS are
  // per kg on every row — one map, uniform tags (E5 fix).
  for (const [name, unit] of Object.entries(CANONICAL_UNIT_TAGS))
    df.col(name)?.setTag('units', unit);
  df.col('CL_ML_MIN_KG')!.setTag('units', CL_PER_KG_UNIT).setTag('nca.deprecated', 'use CL');
  df.col('VZ_L_KG')!.setTag('units', VOLUME_PER_KG_UNIT).setTag('nca.deprecated', 'use VZ');
  df.setTag('nca.resultsBasis', RESULTS_BASIS);
  if (m !== null) {
    df.setTag('nca.mapping', JSON.stringify(m));
    df.setTag('nca.keyColumns', JSON.stringify(keyNames));
    df.setTag('nca.units', JSON.stringify(m.units));
  }
  if (sex !== null && !sex.ok) df.setTag('nca.sexReason', sex.reason);
  if (rulesJson !== undefined) df.setTag('nca.rules', rulesJson);
  return df;
}

/**
 * Run NCA on an observations table. `mappingJson` empty ⇒ auto-detect;
 * `rulesJson` empty ⇒ app defaults. Never throws for **data** problems: a
 * partition that cannot be computed safely comes back as one `Blocked` row with
 * the reason. Configuration errors (malformed JSON, unknown keys, invalid values
 * in `mapping` or `rules`) throw.
 */
export function computePipeline(df: DG.DataFrame, mappingJson: string, rulesJson: string): DG.DataFrame {
  const pr = parseRules(rulesJson);
  let m: PipelineMapping | null;
  let tableIssues: string[] = [];
  if (mappingJson.trim() === '') {
    const det = detectPipelineMapping(df, {perPartitionAmbiguity: true});
    m = det.mapping;
    tableIssues = det.issues;
  } else {
    // Accept either a bare mapping or ncaDetectMapping's {mapping, confidence, issues}.
    const parsed = parseConfigJson('mapping', mappingJson);
    if (isPlainObject(parsed) && 'confidence' in parsed) {
      m = parsed.mapping != null ? validateMapping(parsed.mapping) : null;
      // A table-wide concentration ambiguity is re-checked per study below.
      const issues = Array.isArray(parsed.issues) ? parsed.issues as string[] : ['low-confidence mapping'];
      tableIssues = parsed.confidence === 'low' ? issues.filter((i) => !i.startsWith('ambiguous concentration')) : [];
    } else
      m = validateMapping(parsed);
  }

  if (m === null)
    return toFrame([blockedRow(new Map(), `mapping: ${tableIssues.join('; ')}`)], null, null);
  const missing = missingColumns(df, m);
  if (missing.length > 0) {
    return toFrame([blockedRow(new Map(),
      `mapping names columns not in the table: ${missing.map((c) => `'${c}'`).join(', ')}`)], m, null);
  }
  const gaps = mappingGaps(m);
  if (gaps.length > 0) tableIssues = [`missing ${gaps.join(', ')}`, ...tableIssues];
  const sex = sexOf(df, m);

  const groupCols = m.groupBy.map((g) => df.col(g)!);
  const partitions = groupRows(allRowsOf(df), (r) => groupCols.map((c) => strAt(c, r)).join('\u0001'));
  const mapped = mappedColumns(m);

  const out: OutRow[] = [];
  for (const rows of partitions.values()) {
    const keys = new Map(m.groupBy.map((g, i) => [g, strAt(groupCols[i], rows[0])]));
    if (tableIssues.length > 0) {out.push(blockedRow(keys, `mapping: ${tableIssues.join('; ')}`)); continue;}
    const extra = extraConcColumns(df, m, rows);
    if (extra.length > 0) {
      out.push(blockedRow(keys, `ambiguous concentration: '${m.concentration}' and ` +
        `${extra.map((c) => `'${c}'`).join(', ')} both hold values (reassay?) — confirm which is reportable`));
      continue;
    }

    let units: Array<{rows: number[]; keys: Map<string, string>; reasons: string[]}> =
      [{rows, keys, reasons: []}];
    if (duplicateSamples(df, m, rows) > 0) {
      const split = findSplitColumn(df, m, rows, mapped);
      if (split === null) {
        out.push(blockedRow(keys, 'duplicate samples: a subject has more than one sample at the same time ' +
          'and no column tells them apart — refusing to pool'));
        continue;
      }
      const splitCol = df.col(split)!;
      units = [...groupRows(rows, (r) => strAt(splitCol, r)).entries()].map(([v, part]) => ({
        rows: part,
        keys: new Map([...keys, [split, v]]),
        reasons: [`split by '${split}' (duplicate samples per subject; not in the mapping's groupBy) — confirm`],
      }));
    }

    for (const u of units) {
      const sub = subFrame(df, u.rows);
      const blqName = m.blq !== undefined ? ensureBoolBlq(sub, m.blq) : null;
      const design = samplingDesign(sub, m, blqName, pr.rules.lambdaZ.minPoints);
      try {
        if (design.kind === 'blocked') out.push(blockedRow(u.keys, design.reason));
        else if (design.kind === 'sparse')
          out.push(...computeSparse(sub, m, pr, sex, u.keys, blqName, u.reasons, design.topology));
        else
          out.push(...computeDense(sub, m, pr, sex, u.keys, blqName, u.reasons));
      } catch (e) {
        out.push(blockedRow(u.keys, `NCA error: ${e instanceof Error ? e.message : String(e)}`));
      }
    }
  }
  flagExposureOutliers(out, pr.exposureOutlierFold);
  return toFrame(out, m, sex, rulesJson);
}
