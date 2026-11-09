/**
 * "Not reportable" — the Exclude λz review decision applied to the per-subject results.
 *
 * **One feed, not a mask per consumer.** {@link reportableResults} returns a CLONE of the
 * results frame in which every λz-dependent cell of an excluded profile is the platform
 * null. Every downstream consumer (grouped summary, dose proportionality, F%, TK, the
 * study-grid summary cache) receives that clone and handles the profile through the null
 * path it ALREADY has for a profile whose λz was not estimated — which is what "does not
 * enter the statistics" means in PKNCA's `exclude()` + `summary()` (v07). `N_<P>` of the
 * aggregator then becomes the per-parameter n for free.
 *
 * sci-comp math is untouched (thesis #1): no value is recomputed; only which values reach
 * the report changes. The source results frame is never mutated — the per-subject grid
 * keeps the value and marks it (F-S4.5).
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import {compositeGroupKey} from './composite-group-key';
import {sparseUnitKey} from './review-units';
import type {ReviewUnit} from './review-units';
import type {ReviewDecisionRecord} from './study-schema';
import {REVIEW_STATE_LABEL} from './review-reasons';

/** The results-frame `Review` column (D15) — one neutral word per profile. Written
 *  presentation-side by the results view; read here and by the study summary cache
 *  (`summary-cache.ts`), which counts the excluded profiles per group. Moved from
 *  `views/results-view.ts` (nca-pipeline-save-rework F1) so the cache builder needs no view. */
export const REVIEW_COL = 'REVIEW';
/** The `Review` state of a profile whose λz was excluded (`Exclude λz`). */
export const REVIEW_EXCLUDED = REVIEW_STATE_LABEL.exclude;

/**
 * Result columns derived from the terminal-phase λz fit — the parameters an Exclude λz
 * decision withholds. A superset of D11's named list (t½, AUC∞, CL, Vz, MRT, Vss):
 * %AUCextrap, AUMC∞, %AUMCextrap and λz itself are λz-derived too.
 *
 * **Not in the set, deliberately:** CMAX, TMAX, AUC_LAST, AUMC_LAST (observed-data
 * parameters), C0_EXTRAP (back-extrapolated from the first points, not the terminal
 * phase), TLAG, and every `LAMBDAZ_*` / `SPAN_RATIO` diagnostic — those describe the fit,
 * they are not reported parameters (AC-S1.2.2).
 *
 * Route applicability needs no second discriminator: VSS is already null on a non-IV row,
 * and blanking a null is harmless (C-10).
 */
export const LZ_DEPENDENT_PARAMS: ReadonlyArray<string> = [
  'LAMBDA_Z', 'HALF_LIFE', 'AUC_INF', 'PCT_EXTRAP', 'CL', 'VZ',
  'AUMC_INF', 'MRT', 'VSS', 'PCT_AUMC_EXTRAP',
  // The absolute companions rest on the same λz as their per-kg siblings (nca-cl-bw-normalised-units
  // architect B2), so an Excluded profile blanks them too. BW_KG is not λz-dependent.
  'CL_ABS', 'VZ_ABS', 'VSS_ABS',
];

const LZ_DEPENDENT_SET: ReadonlySet<string> = new Set(LZ_DEPENDENT_PARAMS);

/** True when `colName` is a λz-dependent reported parameter. */
export function isLzDependentParam(colName: string): boolean {
  return LZ_DEPENDENT_SET.has(colName);
}

/**
 * The λz-dependent columns of the UC-04 sparse **composite** parameter table (`sparse-nca.ts`
 * `buildParameterTable`) — what an Exclude λz decision on a composite withholds (a composite
 * exclusion must be as real as a dense one). The composite analogue of {@link LZ_DEPENDENT_PARAMS}:
 * AUC∞, t½, CL, Vz, MRT, Vss, %AUCextrap. **Not in the set:** AUClast (+ its SE / df / CI), Cmax,
 * Tmax (observed composite data) and the fit diagnostics (R², adj. R², ACCEPTANCE, FRAGILITY) — as on the dense grid.
 */
export const SPARSE_LZ_DEPENDENT_COLS: ReadonlyArray<string> = [
  'AUC_INF', 'HALF_LIFE', 'CL', 'VZ', 'MRT', 'VSS', 'PCT_AUCEXTRAP', 'CL_ABS', 'VZ_ABS', 'VSS_ABS',
];

const SPARSE_LZ_DEPENDENT_SET: ReadonlySet<string> = new Set(SPARSE_LZ_DEPENDENT_COLS);

/** True when `colName` is a λz-dependent column of the sparse composite parameter table. */
export function isSparseLzDependentCol(colName: string): boolean {
  return SPARSE_LZ_DEPENDENT_SET.has(colName);
}

/**
 * Which rows of the sparse parameter table are EXCLUDED composites: row `i` is the composite
 * `groupKeys[i]` (the table is built one row per `SparseNcaResult.groups[i]`), excluded when its
 * review key `sparse:<groupKey>` is in `excluded` (the controller's `excludedKeys()`, AD-1).
 */
export function excludedCompositeRows(
  groupKeys: ReadonlyArray<string>, excluded: ReadonlySet<string>,
): boolean[] {
  return groupKeys.map((g) => excluded.has(sparseUnitKey(g)));
}

export interface ReportableResults {
  /** The blanked clone — the ONLY frame the aggregate consumers may read. When nothing is
   *  excluded this is the input frame itself (no clone; consumers are read-only). */
  readonly df: DG.DataFrame;
  /** Per-row: true when that results row's SUBJECT is excluded. */
  readonly excludedMask: ReadonlyArray<boolean>;
  /** Number of excluded rows. */
  readonly excludedCount: number;
  /** Excluded-row count per composite group key over `groupKeys` (the key the aggregator
   *  groups by — {@link compositeGroupKey}). Returns an empty map when a key column is
   *  missing. */
  excludedCountBy(groupKeys: ReadonlyArray<string>): Map<string, number>;
}

/**
 * Blank the λz-dependent cells of every excluded profile, on a clone.
 *
 * @param results   the per-subject results frame (`SUBJECT` column required).
 * @param excluded  the results-row `SUBJECT` keys whose review decision is Exclude λz.
 */
export function reportableResults(
  results: DG.DataFrame, excluded: ReadonlySet<string>,
): ReportableResults {
  const subj = results.col('SUBJECT');
  const n = results.rowCount;
  const mask: boolean[] = new Array<boolean>(n).fill(false);
  let count = 0;
  if (subj !== null && excluded.size > 0) {
    for (let i = 0; i < n; i++)
      if (excluded.has(subj.getString(i))) {mask[i] = true; count++;}
  }
  let df = results;
  if (count > 0) {
    df = results.clone();
    for (const p of LZ_DEPENDENT_PARAMS) {
      const c = df.col(p);
      if (c === null) continue;
      for (let i = 0; i < n; i++) if (mask[i]) c.set(i, DGNS.FLOAT_NULL);
    }
  }
  const excludedCountBy = (groupKeys: ReadonlyArray<string>): Map<string, number> => {
    const out = new Map<string, number>();
    if (count === 0) return out;
    const keyCols: DG.Column[] = [];
    for (const k of groupKeys) {
      const c = results.col(k);
      if (c === null) return out;
      keyCols.push(c);
    }
    for (let i = 0; i < n; i++) {
      if (!mask[i]) continue;
      const key = keyCols.length === 0 ? '' : compositeGroupKey(keyCols, i);
      out.set(key, (out.get(key) ?? 0) + 1);
    }
    return out;
  };
  return {df, excludedMask: mask, excludedCount: count, excludedCountBy};
}

// ─── The Review column (moved from `views/results-view.ts`, nca-hub-demo-studies F4) ──────────
// The study summary cache counts excluded profiles per group from this column, so a headless
// save (the demo seed) must write it exactly as the view does — one writer, in the adapter.

/** Hidden `_lz*` service column carrying each row's review detail as JSON (state, category,
 *  reasons, who/when, BLQ-only disclosure) for the Review / excluded-cell tooltips. An `_lz*`
 *  column, so `isLzServiceColumn` hides it from both grid profiles — the same projection
 *  pattern as `_lz_overridden`. */
export const LZ_REVIEW_INFO_COL = '_lz_review';

/** The six Review states (D15 + SF-3). */
export const REVIEW_TO_REVIEW = 'To review';
const REVIEW_STATE_BY_ACTION: Readonly<Record<string, string>> = REVIEW_STATE_LABEL;

/** One row's review detail (the JSON in {@link LZ_REVIEW_INFO_COL}). */
export interface ReviewInfo {
  readonly state: string;
  readonly category: string;
  readonly decision: ReviewDecisionRecord | null;
  readonly blqOnly: boolean;
  readonly blqTimes: ReadonlyArray<string>;
}

/** The Review word for a classified unit. */
export function reviewStateOf(u: ReviewUnit): string {
  if (u.need === 'not-applicable') return 'Not applicable';
  if (u.need === 'not-needed') return 'Not needed';
  return u.decision !== null ? (REVIEW_STATE_BY_ACTION[u.decision.action] ?? REVIEW_TO_REVIEW) : REVIEW_TO_REVIEW;
}

/** Add-or-replace a STRING column on `df` (idempotent — `configureResultsDataFrame`
 *  claims idempotency, and re-Compute may re-present the same frame). Mirrors the
 *  `columns.add(Column.fromList(...))` idiom of the writer / grouping-columns. */
export function upsertStringColumn(
  df: DG.DataFrame, name: string, values: Array<string | null>,
): void {
  if (df.col(name) !== null) df.columns.remove(name);
  df.columns.add(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, name, values));
}

/**
 * Write (or refresh in place) the Review column + its hidden detail column from the classified
 * units (the controller's `reviewSnapshot`, never re-derived here). Keyed by `SUBJECT`; a row with
 * no unit (should not happen) reads blank. Returns true when the columns were newly ADDED (the
 * caller then re-applies the column profile).
 */
export function writeReviewColumns(df: DG.DataFrame, units: ReadonlyArray<ReviewUnit>): boolean {
  const subj = df.col('SUBJECT');
  if (subj === null) return false;
  const byKey = new Map(units.map((u) => [u.key, u]));
  const n = df.rowCount;
  const states: string[] = new Array(n);
  const infos: string[] = new Array(n);
  for (let i = 0; i < n; i++) {
    const u = byKey.get(subj.getString(i));
    if (u === undefined) {states[i] = ''; infos[i] = ''; continue;}
    const info: ReviewInfo = {
      state: reviewStateOf(u), category: u.category, decision: u.decision,
      blqOnly: u.blqOnly, blqTimes: u.blqTimes,
    };
    states[i] = info.state;
    infos[i] = JSON.stringify(info);
  }
  const existing = df.col(REVIEW_COL);
  if (existing !== null && existing.length === n) {
    const info = df.col(LZ_REVIEW_INFO_COL)!;
    for (let i = 0; i < n; i++) {existing.set(i, states[i], false); info.set(i, infos[i], false);}
    return false;
  }
  upsertStringColumn(df, REVIEW_COL, states);
  upsertStringColumn(df, LZ_REVIEW_INFO_COL, infos);
  df.col(REVIEW_COL)!.setTag('friendlyName', 'Review');
  return true;
}
