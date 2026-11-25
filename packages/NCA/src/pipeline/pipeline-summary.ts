/**
 * `ncaSummarize`: per-arm summary of `ncaCompute` results, with bioavailability.
 *
 * The mean / SD / %CV come from `aggregateNcaResults` (the study-grid aggregation)
 * and F% from `computeBioavailability` (the app's F% table), so the pipeline and
 * the app summarise identically. Added here: the review count per arm, and the
 * reason when an oral arm has no F%.
 *
 * An arm is the result keys (study, compound, period) × `SEX` (when the results
 * carry it) × route × dose. F% is stratified by sex, as in the app (AD-6): without
 * a recognised sex column there is no F%, never a pooled one.
 */

import * as DG from 'datagrok-api/dg';

import {
  aggregateNcaResults, computeBioavailability, isIntravascular, numAt,
  detectOccasionColumn, isOccasionName, profileRowsOf, resolveOccasion,
} from '../adapter';
import {MIXED_SEX, VERDICT_BLOCKED, VERDICT_REVIEW} from './pipeline-compute';
import {NO_SEX_COLUMN} from './pipeline-mapping';

/** `CL` / `VZ` / `VSS` are per kg on every row and the `_ABS` companions absolute
 *  so pooling across DOSE_NORM is one basis per column.
 *  `CL_ML_MIN_KG` / `VZ_L_KG` are the one-release deprecated aliases. */
const SUMMARY_PARAMS = ['CMAX', 'TMAX', 'AUC_LAST', 'AUC_INF', 'HALF_LIFE', 'CL', 'VZ', 'MRT', 'VSS',
  'CL_ABS', 'VZ_ABS', 'VSS_ABS', 'CL_ML_MIN_KG', 'VZ_L_KG'];
/** Sex strata a profile can enter F% with: one known sex per profile. */
const SINGLE_SEX = new Set(['M', 'F', 'U']);

/** Key columns of a results frame: the `nca.keyColumns` tag, else the string columns before SUBJECT. */
export function resultKeyColumns(results: DG.DataFrame): string[] {
  const tag = results.getTag('nca.keyColumns');
  if (tag !== null && tag !== undefined && tag !== '')
    return (JSON.parse(tag) as string[]).filter((c) => results.col(c) !== null);
  const out: string[] = [];
  for (const c of results.columns.toList()) {
    if (c.name === 'SUBJECT') break;
    if (c.type === DG.COLUMN_TYPE.STRING) out.push(c.name);
  }
  return out;
}

const cellKey = (vals: ReadonlyArray<string>, route: string, dose: number): string =>
  [...vals, route, String(dose)].join('\u0001');

interface ArmF {f: number; reason: string; ivN: number}

/**
 * F% per arm of `agg`, from `computeBioavailability` on the profile rows `src`,
 * stratified by the non-period keys and `SEX`. Only single-sex profiles enter
 * it: a composite pooling sexes or a subject whose sex varies must never become
 * a cross-sex reference or test arm.
 */
function bioavailabilityByArm(src: DG.DataFrame, agg: DG.DataFrame, keys: ReadonlyArray<string>,
  noSexReason: string | null, notes: string[]): ArmF[] {
  // F5 (SF-PIPELINE-FPCT): the occasion key (one shared matcher, `detectOccasionColumn`) stays a
  // match key only for a repeat-dose design (same-occasion pairing). For a crossover and for an
  // unclassifiable design it is dropped as before — the IV and oral arms of one crossover sit in
  // different periods. The Pipeline has no animal mapping: the design is read through rule (ii)
  // (the subject id repeats across occasions) and adopted without confirmation (nobody to ask);
  // a period-like name still proposes `unknown`, so crossover pairing is kept there.
  const occasionKey = detectOccasionColumn(keys.filter((k) => isOccasionName(k)));
  let sameOccasion = false;
  if (occasionKey !== null) {
    const occ = resolveOccasion(occasionKey, profileRowsOf(src, {subject: 'SUBJECT', occasion: occasionKey,
      sex: src.col('SEX') !== null ? 'SEX' : null, dose: 'DOSE', route: 'ROUTE'}), {}, /* autoConfirm */ true);
    sameOccasion = occ.fPctPairing === 'same-occasion';
    if (occ.levels.length >= 2 && (occ.confirmed ?? occ.proposal.design) === 'unknown') {
      notes.push(`F% pairs IV and extravascular profiles across '${occasionKey}' (${occ.levels.join(', ')}): ` +
        `the design could not be classified (${occ.proposal.evidence}).`);
    }
  }
  const matchKeys = keys.filter((k) => !isOccasionName(k) || (sameOccasion && k === occasionKey));
  const route = agg.col('ROUTE')!;
  const dose = agg.col('DOSE')!;
  const sexCol = agg.col('SEX');
  const out: ArmF[] = [];
  let cells: Map<string, {f: number; ivN: number}> | null = null;
  if (sexCol !== null) {
    const srcSex = src.col('SEX')!;
    const single = src.clone(DG.BitSet.create(src.rowCount, (i) => SINGLE_SEX.has(srcSex.getString(i))));
    const strat = [...matchKeys, 'SEX'];
    const {table, warnings} = computeBioavailability(single, {doseCol: 'DOSE', routeCol: 'ROUTE', stratifyCols: strat,
      ...(sameOccasion && occasionKey !== null ? {occasionCol: occasionKey} : {})});
    for (const w of warnings) if (w.code === 'CROSS_OCCASION_IV_REFERENCE') notes.push(w.message);
    cells = new Map();
    for (let i = 0; i < table.rowCount; i++) {
      cells.set(cellKey(strat.map((s) => table.col(s)!.getString(i)), table.col('ROUTE')!.getString(i),
        numAt(table.col('DOSE')!, i)), {f: numAt(table.col('F_PCT')!, i), ivN: table.col('IV_N')!.get(i)});
    }
  }
  for (let i = 0; i < agg.rowCount; i++) {
    const none = (reason: string): ArmF => ({f: NaN, reason, ivN: NaN});
    if (isIntravascular(route.getString(i))) {out.push(none('')); continue;}
    if (cells === null || sexCol === null) {out.push(none(noSexReason ?? NO_SEX_COLUMN)); continue;}
    const sex = sexCol.getString(i);
    if (sex === MIXED_SEX) {out.push(none('composite pools sexes')); continue;}
    if (!SINGLE_SEX.has(sex)) {out.push(none('sex varies within subject')); continue;}
    const cell = cells.get(cellKey([...matchKeys.map((k) => agg.col(k)!.getString(i)), sex],
      route.getString(i), numAt(dose, i)));
    if (cell === undefined) out.push(none('no AUCinf/dose in this arm'));
    else if (!Number.isFinite(cell.f)) out.push({f: NaN, reason: 'no IV reference', ivN: cell.ivN});
    else out.push({f: cell.f, reason: '', ivN: cell.ivN});
  }
  return out;
}

export function summarizePipeline(results: DG.DataFrame): DG.DataFrame {
  const keys = resultKeyColumns(results);
  const hasSex = results.col('SEX') !== null;
  const verdict = results.col('VERDICT')!;
  const usable = DG.BitSet.create(results.rowCount, (i) => verdict.get(i) !== VERDICT_BLOCKED);
  const src = results.clone(usable);

  const isReview = src.col('VERDICT')!;
  src.columns.add(DG.Column.fromInt32Array('IS_REVIEW',
    Int32Array.from({length: src.rowCount}, (_, i) => isReview.get(i) === VERDICT_REVIEW ? 1 : 0)));

  const groupBy = [...keys, ...(hasSex ? ['SEX'] : []), 'ROUTE', 'DOSE'];
  const agg = aggregateNcaResults(src, groupBy, {params: SUMMARY_PARAMS});
  const reviewAgg = src.groupBy(groupBy).sum('IS_REVIEW', 'N_REVIEW').aggregate();
  const reviewCol = reviewAgg.col('N_REVIEW')!;
  const keyOf = (df: DG.DataFrame, cols: string[], i: number): string =>
    cols.map((c) => String(df.col(c)!.get(i))).join('\u0001');
  const reviewByKey = new Map<string, number>();
  for (let i = 0; i < reviewAgg.rowCount; i++) reviewByKey.set(keyOf(reviewAgg, groupBy, i), reviewCol.get(i));
  agg.columns.add(DG.Column.fromInt32Array('N_REVIEW',
    Int32Array.from({length: agg.rowCount}, (_, i) => reviewByKey.get(keyOf(agg, groupBy, i)) ?? 0)));

  const sexReason = results.getTag('nca.sexReason');
  const fNotes: string[] = [];
  const f = bioavailabilityByArm(src, agg, keys, sexReason !== null && sexReason !== undefined && sexReason !== '' ?
    sexReason : null, fNotes);
  // Design notes (unknown occasion design, cross-occasion IV reference) ride a frame tag: the
  // summary schema is unchanged.
  if (fNotes.length > 0) agg.setTag('nca.fPctNotes', JSON.stringify(fNotes));
  const fCol = DG.Column.fromFloat64Array('F_PCT',
    Float64Array.from(f, (a) => Number.isFinite(a.f) ? a.f : DG.FLOAT_NULL));
  fCol.setTag('units', '%');
  agg.columns.add(fCol);
  agg.columns.add(DG.Column.fromStrings('F_PCT_REASON', f.map((a) => a.reason)));
  // The same notes as a visible column on the extravascular arms (a headless caller never reads
  // frame tags). Written only when there is a note, so a summary without one keeps its schema.
  if (fNotes.length > 0) {
    const routeCol = agg.col('ROUTE')!;
    const note = fNotes.join(' ');
    agg.columns.add(DG.Column.fromStrings('F_PCT_NOTE', Array.from({length: agg.rowCount}, (_, i) =>
      isIntravascular(routeCol.getString(i)) ? '' : note)));
  }
  // The IV-reference N behind each F% (the app's "N (IV ref.)").
  agg.columns.add(DG.Column.fromInt32Array('F_PCT_IV_N',
    Int32Array.from(f, (a) => Number.isFinite(a.ivN) ? a.ivN : DG.INT_NULL)));

  agg.name = 'NCA summary';
  const units = results.getTag('nca.units');
  // `getTag` on an unset tag returns undefined, not null.
  if (units !== null && units !== undefined) agg.setTag('nca.units', units);
  agg.setTag('nca.keyColumns', JSON.stringify(keys));
  const basis = results.getTag('nca.resultsBasis');
  if (basis !== null && basis !== undefined) agg.setTag('nca.resultsBasis', basis);
  return agg;
}
