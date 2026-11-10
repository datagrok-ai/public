/**
 * The study-grid summary cache (decision h) — the two-tier projection of one compute that a
 * study entry persists, so the study grid renders without recomputing.
 *
 * Moved here from `views/nca-analysis-open.ts` (nca-pipeline-save-rework F1, review 7a) so the
 * headless pipeline (`ncaSaveToStudy`) builds the SAME cache the app writes, without importing
 * the view layer. The function bodies are unchanged; only their home moved.
 */

import * as DG from 'datagrok-api/dg';

import type {PerProfileResult} from '../controller/nca-controller';
import {ROUTE_STRING_TO_CODE} from './profile-inputs';
import {aggregateNcaResults, DEFAULT_AGG_PARAMS} from './results-aggregator';
import {computeBioavailability, isIntravascular} from './bioavailability';
import {buildMeanSdProfile} from './mean-sd-profile';
import {compositeGroupKey, KEY_SEP} from './composite-group-key';
import {dedupGroupValuesByCanonical} from './study-grid';
import {isLzDependentParam, REVIEW_COL, REVIEW_EXCLUDED} from './reportable';
import type {
  AnalysisRecipe, SummaryCache, CompoundTierRow, MeanCv, CompoundCurve, CompoundWarning,
  GroupKeyValue, KeyInventoryEntry, KeyValueCount, PassthroughEntry,
} from './study-schema';

/** The normalized observations frame + its role column names — the substrate for
 *  the v5 inline curve thumbnail (chemist-sar-grid S4). The SAME frame/keys the
 *  analysis view's own mean±SD chart pools from (`MeanSdCtx`), so the cached
 *  thumbnail IS the analysis's own geomean curve, just small. */
export interface ObservationsCtx {
  /** The normalized dense observations frame (TAD-timed, BLQ-resolved). */
  readonly df: DG.DataFrame;
  /** Time column (nominal when mapped, else observed/TAD). */
  readonly timeName: string;
  /** Concentration column. */
  readonly concName: string;
  /** Dose-group column — one curve series per value (BUG-08 no cross-dose
   *  pooling). `null` when no dose is mapped → a single pooled series. */
  readonly doseName: string | null;
  /** The per-row derived RouteCode column ({@link DERIVED_ROUTE_CODE_COL}) — the
   *  canonical `nca.RouteCode` the pipeline itself encoded (from a mapped Route
   *  column in standard mode, or the NONMEM route-resolution in event-driven mode),
   *  the SAME code that determines each tier row's ROUTE. When present, the curve is
   *  keyed per (compound × ROUTE × dose) so an IV decline curve is never pooled with
   *  a PO absorption curve of the same compound (the tier is route-split; the curve
   *  must be too) — for BOTH ingestion modes. `null` for a single-route (scalar)
   *  analysis → one route per compound, so the compound-level key is already
   *  route-correct. */
  readonly routeCodeName: string | null;
}

/**
 * Log-space sorted-gap dose-level clusterer (F2, multi-column-group-by). Counts the
 * distinct dose CLUSTERS among the positive dose values: sort distinct `log(dose)`,
 * split where an adjacent gap exceeds `GAP_FACTOR ×` the median adjacent gap. Runs on
 * `log(dose)` because per-kg jitter is PROPORTIONAL to dose and geometric spacing
 * (1/3/10/30/100/300) is the modal tox dose-ranging convention — a linear threshold
 * would over-split the high arm. theoph (10 per-kg values in a narrow band, 1 nominal
 * level) → 1 cluster (suppress); a jittered ≥2-arm escalation → ≥2 clusters (fire).
 * Heuristic (very-closely-spaced arms remain hard) — the notice is advisory, off the
 * correctness path (Research Gap 1: a mapping-time covariate-vs-factor designation is
 * the exact fix). `GAP_FACTOR` pinned via DATA-GAP-31.
 */
export function countDoseClusters(col: DG.Column): number {
  const set = new Set<number>();
  for (let i = 0; i < col.length; i++) {
    if (col.isNone(i)) continue;
    const v = col.get(i) as number;
    if (Number.isFinite(v) && v > 0) set.add(v);
  }
  const vals = Array.from(set).sort((a, b) => a - b);
  if (vals.length <= 1) return vals.length;
  const logs = vals.map((v) => Math.log(v));
  const gaps: number[] = [];
  for (let i = 1; i < logs.length; i++) gaps.push(logs[i] - logs[i - 1]);
  const sortedGaps = gaps.slice().sort((a, b) => a - b);
  const medGap = sortedGaps[Math.floor((sortedGaps.length - 1) / 2)];
  // GAP_FACTOR pinned against real theoph (DATA-GAP-31): its distinct per-kg doses
  // [3.1 … 5.86] have a largest log-gap ≈ 4.05× the median, so the factor must exceed
  // that to resolve theoph to 1 cluster; a genuine ≥2-arm escalation (linear
  // 9-11/28-32/95-105 or log-spaced 1/3/10/…) has inter-arm log-gaps ≫ 5× the median,
  // so it splits into ≥2. 5 sits between (theoph is "near-boundary" per R2 — accepted;
  // the notice is advisory, off the correctness path).
  const GAP_FACTOR = 5;
  const threshold = Math.max(GAP_FACTOR * medGap, 1e-9);
  let clusters = 1;
  for (const g of gaps) if (g > threshold) clusters++;
  return clusters;
}

/** Inverse of `ROUTE_STRING_TO_CODE` (first canonical label per code) for the
 *  summary's display route string. Built once. */
const ROUTE_CODE_TO_LABEL: ReadonlyMap<number, string> = (() => {
  const m = new Map<number, string>();
  for (const [label, code] of Object.entries(ROUTE_STRING_TO_CODE))
    if (!m.has(code)) m.set(code, label);
  return m;
})();

/** Results-frame columns for the always-present route + dose-normalization
 *  sub-keys. Every case is route-homogeneous (→ CL basis-homogeneous) and
 *  unit-homogeneous by construction (both are implicit group sub-keys). */
const ROUTE_COL = 'ROUTE';
const DOSE_NORM_COL = 'DOSE_NORM';

/** Synthetic composite-key column name (v8, multi-column-group-by) — a control-
 *  char-joined `compositeGroupKey` over the ≥2 real compound keys, added to a
 *  working results/observations clone so the single-column-keyed %F / curve helpers
 *  key by the full composite tuple (SF-1, no cross-dimension pooling). */
const COMPOSITE_COMPOUND_COL = '__nca_composite_group';

/**
 * The observations frame carrying the curve key column that {@link buildCurvesByCompound} reads:
 * the single compound key's own column, or {@link COMPOSITE_COMPOUND_COL} for ≥ 2 keys, so the
 * thumbnail keys byte-identically to the tier's compoundId. A key the observations frame lacks
 * (`CL_BASIS` and any other results-only column) takes, on every observation row of a profile,
 * that profile's value in `resultsDf` (row `i` ↔ `perProfile[i]`, whose `observationRows`
 * index this frame). Such a column is constant within a profile, so the projection is exact. A
 * row of no profile gets `''` and keys no tier row. Before this, a missing key left every row
 * of the analysis without a curve (the lookup fell back to the analysis label).
 */
function keyedObservations(observations: ObservationsCtx, compoundKeyCols: readonly string[],
  resultsDf: DG.DataFrame, perProfile: ReadonlyArray<PerProfileResult>): ObservationsCtx {
  const odf = observations.df;
  const missing = compoundKeyCols.filter((k) => odf.col(k) === null);
  if (compoundKeyCols.length === 0 || (compoundKeyCols.length === 1 && missing.length === 0))
    return observations;
  const profileOf = new Int32Array(odf.rowCount).fill(-1);
  if (missing.length > 0) {
    perProfile.forEach((p, i) => {
      if (i >= resultsDf.rowCount) return;
      for (const r of p.index.observationRows) if (r >= 0 && r < odf.rowCount) profileOf[r] = i;
    });
  }
  const valueOf = compoundKeyCols.map((k): (r: number) => string => {
    const own = odf.col(k);
    if (own !== null) return (r) => own.isNone(r) ? '' : own.getString(r);
    const res = resultsDf.col(k)!;
    return (r) => {
      const i = profileOf[r];
      return i < 0 || res.isNone(i) ? '' : res.getString(i);
    };
  });
  const name = compoundKeyCols.length === 1 ? compoundKeyCols[0] : COMPOSITE_COMPOUND_COL;
  const work = odf.clone();
  // The same join as compositeGroupKey (a missing cell reads '').
  work.columns.add(DG.Column.fromStrings(name, Array.from({length: odf.rowCount}, (_v, r) =>
    valueOf.map((v) => v(r)).join(KEY_SEP))));
  return {...observations, df: work};
}

/** Read a per-parameter aggregate triple `{mean, cvPct, n}` at aggregate row
 *  `i`. An absent column or a null cell → NaN mean/CV and n 0 (honest: no
 *  estimate, not a fake 0). `N_P` is the na.rm per-param contributor count. */
function meanCvAt(agg: DG.DataFrame, param: string, i: number, nExcluded = 0): MeanCv {
  const meanCol = agg.col(`MEAN_${param}`);
  const cvCol = agg.col(`CV_PCT_${param}`);
  const nCol = agg.col(`N_${param}`);
  const mean = (meanCol !== null && !meanCol.isNone(i)) ? (meanCol.get(i) as number) : NaN;
  const cvPct = (cvCol !== null && !cvCol.isNone(i)) ? (cvCol.get(i) as number) : NaN;
  const n = (nCol !== null && !nCol.isNone(i)) ? (nCol.get(i) as number) : 0;
  // v15 (F-S5.2): profiles an Exclude λz decision withheld from this λz-dependent statistic.
  return nExcluded > 0 && isLzDependentParam(param) ? {mean, cvPct, n, nExcluded} : {mean, cvPct, n};
}

/** The codebase's "no estimate" MeanCv sentinel (matches {@link meanCvAt}'s
 *  absent-column result) — a NaN mean renders blank / NE, honest (never a fake 0). */
const NA_MEANCV: MeanCv = {mean: NaN, cvPct: NaN, n: 0};

/** The bioavailability warning code surfaced on a compound-tier row whose EV
 *  %F could not be fully estimated because a stratum had no IV reference. Merged
 *  into the row's `warnings` (→ study-grid Flags tooltip) so a blank %F cell is
 *  honestly explained, not silently dropped (rule 14 / thesis #4). */
const NO_IV_REFERENCE_CODE = 'NO_IV_REFERENCE_FOR_STRATUM';
/** Human-readable message for {@link NO_IV_REFERENCE_CODE} (v6). This warning is
 *  added by the adapter, not sci-comp, so its message is defined here — shown on the
 *  study-grid Flags tooltip alongside sci-comp's own messages. */
const NO_IV_REFERENCE_MSG =
  'No IV reference arm for this compound — absolute %F is not estimable.';

/** Warning code for a study-grid thumbnail that pools across ≥2 dose clusters (E-warn,
 *  study-grid-grouping-columns). E3 pools all of a group's observations into ONE geomean
 *  thumbnail curve (no auto dose-split); when the group spans ≥2 dose levels that pool
 *  mixes dose-dependent concentrations, so it is DISCLOSED in the Warnings cell — BUG-08
 *  honesty via DISCLOSURE (not via splitting), the rule-14-ratified pool+disclose decision
 *  (multi-column-group-by 2026-07-09). Merged into the row's `warnings`. */
export const CROSS_DOSE_POOL_CODE = 'CROSS_DOSE_POOL';
/** Human-readable message for {@link CROSS_DOSE_POOL_CODE}. Reuses the analysis-view
 *  `MEAN_SD_CROSS_DOSE_POOL` vocabulary ("pools across … dose levels; concentration scales
 *  with dose — … Dose to stratify"). */
const CROSS_DOSE_POOL_MSG =
  'Concentration–time thumbnail pools across ≥2 dose levels; concentration scales with ' +
  'dose — group by Dose to stratify.';
/** Absolute dose-ratio (max/min) that marks two distinct dose values as genuinely
 *  different dose LEVELS — the degenerate-case fallback for the {@link CROSS_DOSE_POOL}
 *  disclosure (user decision 2026-07-13, threshold = 2×). `countDoseClusters` is a
 *  scale-invariant log-space clusterer that works from ≥3 distinct values, but is
 *  mathematically degenerate at EXACTLY 2 distinct values (one gap equals its own median →
 *  `gap > 5×median` is never true → it returns 1 cluster regardless of spacing), silently
 *  missing a clean 2-arm dose pool. This ratio test is applied ONLY when there are exactly
 *  2 distinct doses (where the clusterer cannot decide); ≥3 distinct values still route
 *  through the clusterer, so a wide-spread per-kg column (theoph: 10 distinct, ratio ≈1.89)
 *  never reaches this test and never false-fires. Two dose levels ≥2× apart are different
 *  levels; within-arm per-kg jitter from body weight is well under 2×. (Residual: a rare
 *  2-subject absolute-dose group at ≈2× body-weight spread could advisory-false-fire — an
 *  extra disclosure, never a hidden number; logged in TODO.) */
const CROSS_DOSE_RATIO_FACTOR = 2;

/** Result of {@link computeFpctByCompound}: the per-(compound × EV route) %F, plus
 *  the set of compound keys that had ≥1 EV stratum with no IV reference (so the
 *  tier can surface the honest NE warning on those compounds' EV rows). */
interface FpctByCompound {
  /** Keyed `${compoundId}${KEY_SEP}${route}` → the compound's EV-strata %F. */
  readonly fPct: ReadonlyMap<string, MeanCv>;
  /** Compound keys (the compound VALUE, or '' when ungrouped) with a no-IV
   *  stratum — their EV rows carry {@link NO_IV_REFERENCE_CODE}. */
  readonly noIvRefCompounds: ReadonlySet<string>;
}

/**
 * Cross-profile %F per (compound × extravascular route) for the compound tier
 * (v4, chemist-sar-grid S2). Reuses the validated {@link computeBioavailability}
 * (the sex-stratified IV reference is MANDATORY — pooling sex biases %F ~36%)
 * for the ratio, then the SAME {@link aggregateNcaResults} engine
 * every other tier stat uses to reduce the stratum-level F% to the compound grain
 * (mean/%CV/n over the compound's EV strata, FR-403). No ratio and no mean/CV is
 * re-derived here (thesis #1 / rule 5).
 *
 * `fPct` is keyed by `${compoundId}${KEY_SEP}${route}` → the compound's EV-strata
 * %F {@link MeanCv}; a present key with a NaN mean = EV compound with no IV
 * reference (NE, never faked — rule 14). `noIvRefCompounds` carries the honest
 * explanation: a stratum with no IV arm produces a `NO_IV_REFERENCE_FOR_STRATUM`
 * warning from `computeBioavailability`, and its compound's EV tier rows surface
 * that code (blank %F cell → explained on the Flags tooltip, not silently
 * dropped). An EMPTY result = no computable %F (no designated per-subject sex
 * column, no per-row dose, no route, or no EV strata) — every EV tier row then
 * gets {@link NA_MEANCV} and the study grid drops the whole %F column (Decision 2:
 * no paired arm anywhere → no column). `compoundCol` is the grouping column NAME,
 * or null when the analysis is ungrouped (%F matched by route within the whole
 * analysis).
 *
 * Sex resolution mirrors the analysis-view F% decorator: the recipe's designated
 * `sexCol` counts ONLY when it was carried into the results frame (constant-per-
 * subject); a non-carried designation → no column → no %F (safe degradation).
 */
function computeFpctByCompound(
  resultsDf: DG.DataFrame,
  recipe: AnalysisRecipe,
  compoundCol: string | null,
): FpctByCompound {
  const empty: FpctByCompound = {fPct: new Map(), noIvRefCompounds: new Set()};
  const sexCol = recipe.columns.sexCol;
  const doseCol = recipe.columns.doseCol;
  if (sexCol === null || resultsDf.col(sexCol) === null) return empty;
  if (doseCol === null || resultsDf.col(doseCol) === null) return empty;
  if (resultsDf.col(ROUTE_COL) === null) return empty;

  const stratifyCols = compoundCol !== null ? [compoundCol, sexCol] : [sexCol];
  let table: DG.DataFrame;
  try {
    ({table} = computeBioavailability(resultsDf, {doseCol, routeCol: ROUTE_COL, stratifyCols}));
  } catch {
    return empty; // preconditions held but compute failed — never break the summary
  }
  if (table.rowCount === 0) return empty; // no EV strata to divide by an IV reference

  // An EMPTY F_PCT cell in the stratum table = an EV stratum with NO IV reference
  // (bioavailability.ts emits NO_IV_REFERENCE_FOR_STRATUM for exactly these). Collect
  // the affected compounds directly off the table (no warning-string parsing) so
  // their EV tier rows can surface the honest explanation.
  //
  // **`isNone`, NOT `Number.isNaN`** (nca-dg-native-missing-values F2). `F_PCT` is
  // built by `f64Column`, so an absent %F is `DG.FLOAT_NULL` — and
  // `Number.isNaN(2.68e-34)` is FALSE. The old predicate would `continue` on every
  // row, leave this set EMPTY, and make the honest "no IV reference" disclosure
  // silently disappear from the study grid. A bare `get()` is no help either: it
  // returns the raw sentinel even on a null cell (probe PR-2).
  //
  // Deliberately `isNone` ALONE, not `isNone || !isFinite`: a `+Infinity` F_PCT
  // means the IV arm EXISTS but had zero exposure, which is a different fact.
  // Reporting it as "no IV reference" would fabricate a disclosure. That defect is
  // filed as its own gap.
  const fCol = table.col('F_PCT')!;
  const tblCompoundCol = compoundCol !== null ? table.col(compoundCol) : null;
  const noIvRefCompounds = new Set<string>();
  for (let i = 0; i < table.rowCount; i++) {
    if (!fCol.isNone(i)) continue;
    noIvRefCompounds.add(
      (tblCompoundCol !== null && !tblCompoundCol.isNone(i)) ? tblCompoundCol.getString(i) : '');
  }

  // Aggregate the stratum-level F% to the compound grain with the validated
  // aggregator (F_PCT as the only param). An unestimable %F cell is already the
  // platform null when `computeBioavailability` builds the column (`f64Column`),
  // so it is skipped by the na.rm mean: `N_F_PCT` counts only estimable strata and
  // a compound with no estimable stratum yields a null mean → NE. (The aggregator
  // used to perform that NaN→null conversion itself on a shadow clone; that
  // workaround was deleted once the writers became correct — F3.)
  const aggKeys = compoundCol !== null ? [compoundCol, ROUTE_COL] : [ROUTE_COL];
  const fAgg = aggregateNcaResults(table, aggKeys, {params: ['F_PCT']});
  const cCol = compoundCol !== null ? fAgg.col(compoundCol) : null;
  const rCol = fAgg.col(ROUTE_COL);
  const fPct = new Map<string, MeanCv>();
  for (let i = 0; i < fAgg.rowCount; i++) {
    const compoundKey = (cCol !== null && !cCol.isNone(i)) ? cCol.getString(i) : '';
    const routeKey = (rCol !== null && !rCol.isNone(i)) ? rCol.getString(i) : '';
    fPct.set(`${compoundKey}${KEY_SEP}${routeKey}`, meanCvAt(fAgg, 'F_PCT', i));
  }
  return {fPct, noIvRefCompounds};
}

/** Result of {@link buildCurvesByCompound}: the curve map plus whether it is keyed
 *  per (compound × route) (`routeKeyed`, when the analysis carries a source route
 *  column) or per compound alone (single-route analysis). The caller builds the
 *  matching lookup key. */
export interface CurvesByCompound {
  /** `true` → `byKey` is keyed `${compoundId}${KEY_SEP}${route}`; `false` → keyed
   *  by `compoundId`. */
  readonly routeKeyed: boolean;
  readonly byKey: Map<string, CompoundCurve>;
  /** Outer keys (SAME keying as `byKey`) whose single pooled series spans ≥2 dose
   *  clusters (`countDoseClusters`). The pooled thumbnail mixes dose-dependent
   *  concentrations for these — the caller emits a `CROSS_DOSE_POOL` disclosure on the
   *  matching tier rows (BUG-08 honesty via disclosure, not via splitting). */
  readonly crossDoseOuters: Set<string>;
}

/**
 * Build the inline per-compound curve thumbnails (v5, chemist-sar-grid S4; regrouped
 * study-grid-grouping-columns E3). Reuses the validated {@link buildMeanSdProfile}
 * (thesis #1 / rule 5 — NOT a new curve builder): the SAME geometric-mean semi-log curve
 * the analysis's own mean±SD chart plots, so the cached thumbnail is that curve, just small.
 *
 * **One pooled {@link CurveSeries} per (compound × ROUTE) identity (E3).** The vestigial
 * inner dose-split is DROPPED: each outer pools ALL its observations into ONE geomean
 * series, matching the analysis view's Mean±SD chart (the pool-when-ungrouped decision,
 * multi-column-group-by 2026-07-09). This finally propagates that fix to the grid — a
 * whole-analysis theoph thumbnail is ONE pooled curve, not ~10 dose lines. Route is STILL a
 * key dimension (an IV-bolus decline curve must never pool with a PO absorption curve): when
 * the observations carry the pipeline's per-row derived RouteCode column
 * ({@link ObservationsCtx.routeCodeName} = `DERIVED_ROUTE_CODE_COL`), each observation's
 * canonical `nca.RouteCode` — the SAME code that set the tier row's ROUTE — is read directly
 * and the curve key gains a route dimension, for BOTH standard and NONMEM ingestion.
 * When the analyst grouped BY dose, dose is a real compound key so the outer is already
 * per-dose-group and the pool is within one dose (no cross-dose mixing). Otherwise a pool
 * that spans ≥2 dose CLUSTERS is DISCLOSED via {@link CurvesByCompound.crossDoseOuters}
 * (BUG-08 honesty). Points are the geometric mean at each observed/nominal timepoint (NO
 * lossy downsample), sorted ascending by time, measurable (> 0) only.
 *
 * The single series' `group` label is the outer's compound identity; `buildSummaryCache`
 * relabels it to the tier row's group (its `groupValues` joined by ` · `, or `''` for a
 * whole-analysis row) so a composite key reads `M · 10`, not the internal ␁-join — and so a
 * FRESH curve matches what the v8→v9 migration writes on a preserved curve.
 *
 * Efficiency: ONE `buildMeanSdProfile` pass over the outer key (not a per-compound clone),
 * then partitioned back. `compoundColName` is the group-by column NAME on the observations
 * frame, or null when ungrouped (every observation → the single `fallbackLabel` compound).
 * Returns an empty map when no observations / no time or conc column — every tier row's
 * `curve` is then null.
 */
export function buildCurvesByCompound(
  observations: ObservationsCtx | null,
  compoundColName: string | null,
  fallbackLabel: string,
): CurvesByCompound {
  const empty: CurvesByCompound = {routeKeyed: false, byKey: new Map(), crossDoseOuters: new Set()};
  if (observations === null) return empty;
  const {df, timeName, concName, doseName, routeCodeName} = observations;
  if (df.col(timeName) === null || df.col(concName) === null) return empty;

  const compoundCol = compoundColName !== null ? df.col(compoundColName) : null;
  const doseCol = (doseName !== null && df.col(doseName) !== null) ? df.col(doseName) : null;

  // Per-observation canonical route label (matching the tier's ROUTE), read directly
  // from the pipeline's own per-row derived RouteCode column (`DERIVED_ROUTE_CODE_COL`)
  // — the SAME code that set each tier row's ROUTE, so the labels join exactly, for
  // BOTH standard and NONMEM ingestion. No re-normalization (the codes are already
  // encoded). Absent → compound-level keying (single-route analysis).
  let routeLabels: Array<string> | null = null;
  const routeCodeCol = (routeCodeName !== null && df.col(routeCodeName) !== null) ?
    df.col(routeCodeName) : null;
  if (routeCodeCol !== null) {
    routeLabels = new Array<string>(df.rowCount);
    for (let i = 0; i < df.rowCount; i++) {
      routeLabels[i] = routeCodeCol.isNone(i) ? '' :
        (ROUTE_CODE_TO_LABEL.get(routeCodeCol.get(i) as number) ?? '');
    }
  }
  const routeKeyed = routeLabels !== null;

  // E3 curve key = the OUTER identity ONLY: compound (or compound␁route when route-keyed).
  // No inner dose segment — each outer pools into one series.
  const CURVE_KEY = '__curveKey';
  const outerOf = (i: number): string => {
    const compound = (compoundCol !== null && !compoundCol.isNone(i)) ?
      compoundCol.getString(i) : fallbackLabel;
    return routeLabels !== null ? `${compound}${KEY_SEP}${routeLabels[i]}` : compound;
  };
  const keys = new Array<string>(df.rowCount);
  // Per-outer positive dose values for the cross-dose-pool disclosure (the curve no longer
  // splits by dose, so a ≥2-cluster pool must be DISCLOSED, not silently pooled).
  const dosesByOuter = new Map<string, number[]>();
  for (let i = 0; i < df.rowCount; i++) {
    const outer = outerOf(i);
    keys[i] = outer;
    if (doseCol !== null && !doseCol.isNone(i)) {
      const d = Number(doseCol.getString(i));
      if (Number.isFinite(d) && d > 0) {
        let arr = dosesByOuter.get(outer);
        if (arr === undefined) {arr = []; dosesByOuter.set(outer, arr);}
        arr.push(d);
      }
    }
  }
  const work = df.clone();
  work.columns.add(DG.Column.fromStrings(CURVE_KEY, keys));

  // ONE geomean profile over the outer key → rows of (Group=outer, Time, geoMean, band).
  // Partition back into per-outer point lists (one pooled series each). v11: carry the
  // geoLower/geoUpper ×/÷ gSD band alongside geoMean so the context panel can draw whiskers
  // cache-only (already computed here by buildMeanSdProfile, previously DISCARDED). Each band
  // element is FLOAT_NULL where the spread is inestimable (n=1 / degenerate) — the honest
  // no-whisker point (mean-sd-profile.ts null idiom); the FRESH curve still carries the arrays
  // (present, distinguishing it from a legacy band-absent curve).
  const profile = buildMeanSdProfile(work, {timeCol: timeName, concCol: concName, groupCol: CURVE_KEY});
  const gCol = profile.col('Group')!;
  const tCol = profile.col('Time')!;
  const geoCol = profile.col('geoMean')!;
  const loCol = profile.col('geoLower')!;
  const upCol = profile.col('geoUpper')!;
  const ptsByOuter = new Map<string, Array<{time: number; geoMean: number; geoLower: number; geoUpper: number}>>();
  for (let i = 0; i < profile.rowCount; i++) {
    if (geoCol.isNone(i)) continue; // all-BLQ / no positive sample → not plottable on a log axis
    const geo = geoCol.get(i) as number;
    if (!Number.isFinite(geo) || geo <= 0) continue;
    const outer = gCol.getString(i);
    let pts = ptsByOuter.get(outer);
    if (pts === undefined) {pts = []; ptsByOuter.set(outer, pts);}
    // FLOAT_NULL where the band element is inestimable (n=1) — kept as the honest no-whisker
    // sentinel so the panel omits the whisker there (never a fake collapsed-to-point band).
    pts.push({
      time: tCol.get(i) as number, geoMean: geo,
      geoLower: loCol.isNone(i) ? DG.FLOAT_NULL : (loCol.get(i) as number),
      geoUpper: upCol.isNone(i) ? DG.FLOAT_NULL : (upCol.get(i) as number),
    });
  }

  // The single-series group label = the outer's compound identity (strip the route
  // segment when route-keyed). buildSummaryCache relabels it to the tier's group.
  const compoundOf = (outer: string): string => {
    if (!routeKeyed) return outer;
    const sep = outer.lastIndexOf(KEY_SEP);
    return sep >= 0 ? outer.slice(0, sep) : outer;
  };
  const byKey = new Map<string, CompoundCurve>();
  const crossDoseOuters = new Set<string>();
  for (const [outer, pts] of ptsByOuter) {
    if (pts.length === 0) continue;
    pts.sort((a, b) => a.time - b.time); // polyline order = ascending time
    byKey.set(outer, {series: [{
      group: compoundOf(outer),
      time: pts.map((p) => p.time),
      geoMean: pts.map((p) => p.geoMean),
      // v11: persist the ×/÷ gSD band (present on a fresh curve — distinguishes it from a
      // legacy band-absent one); FLOAT_NULL elements are the honest no-whisker points.
      geoLower: pts.map((p) => p.geoLower),
      geoUpper: pts.map((p) => p.geoUpper),
    }]});
    // Disclose a pool spanning ≥2 dose levels. ≥3 distinct doses → the scale-invariant
    // log-space clusterer (handles per-kg jitter; theoph's 10 per-kg values → 1 cluster, no
    // false fire). EXACTLY 2 distinct doses → the clusterer is degenerate (returns 1
    // always), so fall back to an absolute max/min ratio ≥ CROSS_DOSE_RATIO_FACTOR (2×) —
    // the fix for the silent 2-arm-pool honesty gap (user decision 2026-07-13).
    const doses = dosesByOuter.get(outer);
    if (doses !== undefined && doses.length > 0) {
      const distinct = new Set(doses);
      let cross = false;
      if (distinct.size >= 3)
        cross = countDoseClusters(DG.Column.fromFloat64Array('d', Float64Array.from(doses))) >= 2;
      else if (distinct.size === 2) {
        let mn = Infinity; let mx = -Infinity;
        for (const d of distinct) {if (d < mn) mn = d; if (d > mx) mx = d;}
        cross = mx / mn >= CROSS_DOSE_RATIO_FACTOR;
      }
      if (cross) crossDoseOuters.add(outer);
    }
  }
  return {routeKeyed, byKey, crossDoseOuters};
}

/** Cardinality cap for an inventory candidate column (E4). A column with more than this
 *  many distinct non-null values is treated as continuous (subject-id / time / conc / a
 *  per-observation measure) and excluded. 20 gives ~2× headroom over the largest genuine
 *  grouping factor observed across the fixtures (Dose ≤ ~6–10 levels; Sex 2; Route ≤ 3)
 *  while still excluding continuous columns (distinct ≈ rowCount ≫ 20). Provisional /
 *  fixture-grounded (RG-3) — raise it, never silently truncate a real factor, if a fixture's
 *  genuine grouping factor ever exceeds it. */
const INVENTORY_CARDINALITY_CAP = 20;

/**
 * Build the per-tier-row candidate-column inventory (E4, study-grid-grouping-columns):
 * a map from a tier row's composite key (over `keys`, the SAME `compositeGroupKey` the
 * aggregator groups on) to the {@link KeyInventoryEntry}[] for that subpopulation. For each
 * candidate grouping column it records the distinct non-null values + SUBJECT counts among
 * the resultsDf rows of that partition — so the study grid can disclose a key present in the
 * data but NOT grouped by as its counted composition (`F (15); M (15)`), scoped to each
 * row's OWN subpopulation (a route/dose-norm-split whole-analysis case then shows per-route
 * counts, not one analysis-wide count — R2 finding 3).
 *
 * **Candidate-column gate (both must hold).** A column enters the inventory iff it is:
 *  (1) NOT a grouping key of THIS analysis (`keys` — already carried in `groupValues`, so
 *      the grid's precedence step 1 handles it), NOT a non-grouping recipe role
 *      (subject / time / concentration / body-weight / LLOQ / infusion-duration / BLQ-flag),
 *      NOT an internal grouping/basis column (ROUTE / DOSE_NORM / CL_BASIS / BW_KG), and NOT an NCA
 *      result parameter ({@link DEFAULT_AGG_PARAMS}); AND
 *  (2) low-cardinality: distinct non-null count ≤ {@link INVENTORY_CARDINALITY_CAP}.
 * Dose is NOT excluded (a grouping factor with few levels — it rides the cardinality cap
 * exactly like Sex/Route). The role gate is the belt to the cardinality suspenders (a tiny
 * study where a continuous column happens to have ≤ cap distinct values).
 */
function buildKeyInventoryByTier(
  resultsDf: DG.DataFrame,
  keys: readonly string[],
  recipe: AnalysisRecipe,
): Map<string, KeyInventoryEntry[]> {
  const c = recipe.columns;
  const nonGrouping = new Set<string>(
    [c.subject, c.time, c.concentration, c.bodyWeightCol, c.lloqCol,
      c.infusionDurationCol, c.blqFlagCol].filter((x): x is string => x !== null));
  // `BW_KG` is the mapped body-weight role in kg (nca-cl-bw-normalised-units AD-3): excluded like
  // the role itself (`nonGrouping`), never a grouping key.
  const internal = new Set<string>([ROUTE_COL, DOSE_NORM_COL, 'CL_BASIS', 'BW_KG']);
  const aggParams = new Set<string>(DEFAULT_AGG_PARAMS);
  const keySet = new Set<string>(keys);

  // Candidate columns (present, gate-passing, low-cardinality), in frame order. The
  // cardinality CAP is the E4-only half of the gate (a per-analysis pre-partition filter) —
  // it is deliberately NOT shared with the passthrough gate (RG-D: a user-selected SMILES
  // column ~1 distinct/compound ≫ 20 must survive there).
  const candidates: DG.Column[] = [];
  for (const name of resultsDf.columns.names()) {
    if (keySet.has(name) || nonGrouping.has(name) || internal.has(name) || aggParams.has(name))
      continue;
    const col = resultsDf.col(name);
    if (col === null) continue;
    const distinct = new Set<string>();
    let overCap = false;
    for (let i = 0; i < col.length; i++) {
      if (col.isNone(i)) continue;
      distinct.add(col.getString(i));
      if (distinct.size > INVENTORY_CARDINALITY_CAP) {overCap = true; break;}
    }
    if (!overCap && distinct.size > 0) candidates.push(col);
  }
  if (candidates.length === 0) return new Map();

  // The gate-agnostic partition-and-count is shared with buildPassthroughByTier via
  // accumulateByTier (only the candidate-list construction — capped here, cap-free there —
  // differs; RG-D preserved). The KeyInventoryEntry shape IS the accumulator's {key, values}.
  return accumulateByTier(resultsDf, keys, candidates, keyLabelOf);
}

/** A candidate column's friendly display label for the study grid (which normalizes it via
 *  `normalizeGroupKey` at merge time). `getTag` returns `undefined` for an unset tag — guard
 *  against null / undefined / '' (the "undefined undefined" header footgun). A passthrough
 *  column carried onto the results frame has NO friendlyName (`buildPerSubjectColumn` strips
 *  tags), so this returns its source name — which is exactly the recipe spec's `name`. */
function keyLabelOf(col: DG.Column): string {
  const fn = col.getTag('friendlyName');
  return (fn !== null && fn !== undefined && fn !== '') ? fn : col.name;
}

/**
 * The gate-agnostic partition-and-count shared by {@link buildKeyInventoryByTier} (E4) and
 * {@link buildPassthroughByTier} (v10) — architect finding 3. Partitions `resultsDf`'s subject
 * rows by the composite tier key (over `keys`, the SAME `compositeGroupKey` the aggregator groups
 * on) and, per partition per `candidate` column, accumulates each distinct non-null value's
 * SUBJECT count. Returns a map from tier key to a per-candidate `{key, values}` list (candidate
 * order preserved; a candidate with no value in a tier is omitted for that tier). The two callers
 * differ ONLY in how they build `candidates` (capped vs cap-free) — the count mechanics are
 * identical, so this extraction is behavior-preserving for the E4 path (AC-F3e).
 */
function accumulateByTier(
  resultsDf: DG.DataFrame,
  keys: readonly string[],
  candidates: readonly DG.Column[],
  labelOf: (col: DG.Column) => string,
): Map<string, {key: string; values: KeyValueCount[]}[]> {
  if (candidates.length === 0) return new Map();
  const keyCols = keys.map((k) => resultsDf.col(k)!);
  const byTier = new Map<string, Map<string, Map<string, number>>>();
  for (let i = 0; i < resultsDf.rowCount; i++) {
    const tierKey = compositeGroupKey(keyCols, i);
    let perCol = byTier.get(tierKey);
    if (perCol === undefined) {perCol = new Map(); byTier.set(tierKey, perCol);}
    for (const col of candidates) {
      if (col.isNone(i)) continue;
      const v = col.getString(i);
      let counts = perCol.get(col.name);
      if (counts === undefined) {counts = new Map(); perCol.set(col.name, counts);}
      counts.set(v, (counts.get(v) ?? 0) + 1);
    }
  }

  const result = new Map<string, {key: string; values: KeyValueCount[]}[]>();
  for (const [tierKey, perCol] of byTier) {
    const entries: {key: string; values: KeyValueCount[]}[] = [];
    for (const col of candidates) {
      const counts = perCol.get(col.name);
      if (counts === undefined || counts.size === 0) continue;
      const values: KeyValueCount[] = Array.from(counts, ([value, count]) => ({value, count}));
      entries.push({key: labelOf(col), values});
    }
    result.set(tierKey, entries);
  }
  return result;
}

/**
 * Build the per-tier-row PASSTHROUGH inventory (v10, study-grid-passthrough-columns) — the
 * cache-only render substrate for the user-selected annotation columns (SMILES → structures,
 * name, series, …). SEPARATE from {@link buildKeyInventoryByTier} in ONE load-bearing way (RG-D,
 * BLOCKING): the candidate list is the EXPLICIT user pick (`recipe.passthroughCols`) with NO
 * cardinality cap — a 25-compound cassette's SMILES column (~1 distinct/compound ≫ 20) MUST
 * survive, whereas E4's cap would drop it. The dialog gate already vetted the columns
 * (constant-within-subject, not-a-role, no name-collision), so the build trusts the recipe. The
 * gate-agnostic partition-and-count is shared via {@link accumulateByTier}; only this cap-free
 * candidate construction is bespoke.
 *
 * `isStructure` is copied verbatim from each {@link PassthroughColumnSpec} (single source of
 * truth, LBC-2) — NEVER re-derived from `resultsDf` (its columns lost their `semType` in
 * `buildPerSubjectColumn`) nor re-detected at recompute (would flicker). A passthrough column
 * the recipe names but the results frame does not carry is simply absent (renders `—`).
 */
function buildPassthroughByTier(
  resultsDf: DG.DataFrame,
  keys: readonly string[],
  recipe: AnalysisRecipe,
): Map<string, PassthroughEntry[]> {
  const specs = recipe.passthroughCols ?? [];
  if (specs.length === 0) return new Map();
  // Cap-free candidate list: the selected columns present on the results frame, in pick order.
  // The structure designation rides the spec (keyed by name = the results-frame column name).
  const isStructureByName = new Map<string, boolean>(specs.map((s) => [s.name, s.isStructure]));
  const candidates: DG.Column[] = [];
  // One candidate per DISTINCT spec name — the function's own contract, not a guard against a
  // specific caller (BUG-11 stress). `accumulateByTier` keys its counts by `col.name`, so two
  // same-named candidates would increment the SAME counter and DOUBLE every subject count in the
  // disclosure (`F (30); M (30)` for 15 + 15) — a wrong number that reads as real data. The
  // dialog (add-to-study) and `MIGRATIONS[11]` both dedupe upstream, so a duplicate should never
  // arrive; this makes the count correct by construction rather than by the caller's good
  // behaviour, because the failure mode here is silent rather than loud.
  const seen = new Set<string>();
  for (const s of specs) {
    if (seen.has(s.name)) continue;
    const col = resultsDf.col(s.name);
    if (col === null) continue;
    seen.add(s.name);
    candidates.push(col);
  }
  if (candidates.length === 0) return new Map();

  const raw = accumulateByTier(resultsDf, keys, candidates, keyLabelOf);
  // Decorate each accumulated entry with its persisted structure flag. `entry.key` is the
  // column's source name (keyLabelOf → col.name, friendlyName stripped), so it keys back into
  // the spec map by name.
  const result = new Map<string, PassthroughEntry[]>();
  for (const [tierKey, entries] of raw) {
    result.set(tierKey, entries.map((e) => ({
      key: e.key, values: e.values, isStructure: isStructureByName.get(e.key) ?? false,
    })));
  }
  return result;
}

/**
 * Derive the two-tier {@link SummaryCache} from a compute (decision h).
 *
 * **The compound tier is per-(`groupByCol` × route × dose-norm) via the
 * VAL-01-validated {@link aggregateNcaResults}** (v2, study-unified-grid AD-1) —
 * the single grouping/statistics engine; this function computes NO mean/CV/n of
 * its own (thesis #1 / rule 5). `route` 1:1-determines the CL basis (`CL` vs
 * `CL/F`) and `DOSE_NORM` gives unit-homogeneity, so every case row is basis-
 * AND unit-uniform (honest Mean|%CV|n, thesis #4). Keys are deduped + filtered
 * to present columns, so a stale `groupByCol` (mapping changed) safely falls
 * back to per-(route × dose-norm), and grouping by ROUTE/DOSE_NORM itself never
 * double-adds a key.
 *
 * The `analysisTier` (Name/status the overview shows; the compute time is the stored
 * `results_computed_at`, ems-modification) is derived
 * from `perProfile` as before. `resultsDf` is the post-`appendGroupingColumns`
 * frame (carries the grouping metadata); `perProfile[i]` aligns with
 * `resultsDf` row `i` (both in profiles order), so per-group warning codes are
 * collected by the same group key the aggregate carries through.
 */
export function buildSummaryCache(
  resultsDf: DG.DataFrame,
  perProfile: ReadonlyArray<PerProfileResult>,
  recipe: AnalysisRecipe,
  label: string,
  observations: ObservationsCtx | null = null,
): SummaryCache {
  const nSubjects = perProfile.length;
  const failures = perProfile.filter((p) => p.result.status === 'failed').length;
  const warningCount = perProfile.reduce(
    (a, p) => a + p.result.provenance.warnings.length, 0);
  const status: SummaryCache['analysisTier']['status'] =
    (nSubjects > 0 && failures === nSubjects) ? 'error' :
      (failures > 0 || warningCount > 0) ? 'warn' : 'ok';
  const flagSummary = warningCount === 0 && failures === 0 ?
    'No flags' :
    `${failures} failed, ${warningCount} warning${warningCount === 1 ? '' : 's'}`;

  // Modal route across profiles → analysis-tier display label.
  const routeCounts = new Map<number, number>();
  for (const p of perProfile) routeCounts.set(p.route, (routeCounts.get(p.route) ?? 0) + 1);
  let modalRoute = -1; let modalN = -1;
  for (const [code, n] of routeCounts) if (n > modalN) {modalN = n; modalRoute = code;}
  const route = ROUTE_CODE_TO_LABEL.get(modalRoute) ?? '';

  const dose = recipe.columns.doseSource === 'value' ?
    `${recipe.columns.doseScalar} ${recipe.units.dose}` : '(per-row)';

  // Group key = [...groupByCols, ROUTE, DOSE_NORM], deduped + present-only (v8:
  // groupByCols is a set, so a composite (Sex, Dose, …) selection contributes each
  // user key ahead of the implicit route / dose-norm sub-keys).
  const keys = Array.from(new Set(
    [...recipe.groupByCols, ROUTE_COL, DOSE_NORM_COL]
      .filter((c): c is string => resultsDf.col(c) !== null)));

  // Real compound key columns = the user's selected keys that are NOT the implicit
  // ROUTE / DOSE_NORM sub-keys (those ride their own Route / unit columns), present
  // on the results frame, in pick order. Empty ⇒ whole-analysis.
  const compoundKeyCols = recipe.groupByCols.filter(
    (k) => k !== ROUTE_COL && k !== DOSE_NORM_COL && resultsDf.col(k) !== null);

  const compoundTier: CompoundTierRow[] = [];
  if (resultsDf.rowCount > 0 && keys.length > 0) {
    // Per-group union of warning codes. perProfile[i] ↔ resultsDf row i
    // (appendGroupingColumns aligns by index); the aggregate carries the same
    // group columns through, so the join is by the identical key string.
    // The join delimiter is U+0001 (a control char absent from any route /
    // dose-norm / group value), so ['IV','bolus'] can never collide with
    // ['IVbolus','']. (It renders invisibly in diffs — it is NOT an empty join.)
    const keyStr = (col: (k: string) => DG.Column | null, i: number): string =>
      // [multi-column-group-by: unified onto the shared compositeGroupKey helper so
      // every composite-key consumer agrees byte-for-byte. Already U+0001 (NOT the
      // empty join the architect gate flagged — that was a misread of the invisible
      // control char); this refactor is behaviour-preserving.]
      compositeGroupKey(keys.map((k) => col(k)!), i);
    // Per-group warnings as code→message (v6): dedupe by code, keep sci-comp's own
    // human-readable message for the study grid's Flags tooltip.
    const perGroupWarnings = new Map<string, Map<string, string>>();
    for (let i = 0; i < perProfile.length; i++) {
      const gk = keyStr((k) => resultsDf.col(k), i);
      let m = perGroupWarnings.get(gk);
      if (m === undefined) {m = new Map<string, string>(); perGroupWarnings.set(gk, m);}
      for (const w of perProfile[i].result.provenance.warnings)
        if (!m.has(w.code)) m.set(w.code, w.message);
    }

    const agg = aggregateNcaResults(resultsDf, keys);
    // v15 (F-S5.2): excluded profiles per group, from the Review column the reportable frame
    // carries (the per-subject grid's own state — never re-derived), keyed like the aggregate.
    const excludedByGroup = new Map<string, number>();
    {
      const rc = resultsDf.col(REVIEW_COL);
      if (rc !== null) {
        for (let r = 0; r < resultsDf.rowCount; r++) {
          if (rc.isNone(r) || rc.getString(r) !== REVIEW_EXCLUDED) continue;
          const gk = keyStr((k) => resultsDf.col(k), r);
          excludedByGroup.set(gk, (excludedByGroup.get(gk) ?? 0) + 1);
        }
      }
    }
    // Real compound key columns AS PRESENT ON the aggregate (the `keys` group-by
    // includes them, so the aggregate carries them through). A composite selection
    // has several; a single-key selection has one; a whole-analysis (or implicit
    // ROUTE/DOSE_NORM-only) selection has none — those ride their own Route / unit
    // columns and fall back to the analysis label (design B2).
    const compoundKeyAggCols = compoundKeyCols
      .map((k) => agg.col(k))
      .filter((c): c is DG.Column => c !== null);
    const routeCol = agg.col(ROUTE_COL);
    const doseNormCol = agg.col(DOSE_NORM_COL);
    const nCol = agg.col('N');

    // [SF-1, multi-column-group-by] %F + curve thumbnails are keyed by a SINGLE
    // column name, but a composite selection has no single column — passing
    // groupByCols[0] would pool %F/curve across the un-selected dimension while the
    // Group label shows the full composite (BUG-08-class silent cross-dose pooling).
    // Fix: for ≥2 real keys, synthesize a composite key COLUMN (compositeGroupKey,
    // U+0001) onto a working results clone AND the observations frame, and pass THAT
    // column's name — so %F/curve are confined to the full composite cell. A single
    // real key passes the real column directly (no synthesis, byte-identical to the
    // pre-v8 path); whole-analysis passes null (matched by route within the analysis).
    let fpctFrame = resultsDf;
    let compoundColName: string | null = null;
    if (compoundKeyCols.length === 1)
      compoundColName = compoundKeyCols[0];
    else if (compoundKeyCols.length >= 2) {
      compoundColName = COMPOSITE_COMPOUND_COL;
      const rCols = compoundKeyCols.map((k) => resultsDf.col(k)!);
      fpctFrame = resultsDf.clone();
      fpctFrame.columns.add(DG.Column.fromStrings(COMPOSITE_COMPOUND_COL,
        Array.from({length: resultsDf.rowCount}, (_v, i) => compositeGroupKey(rCols, i))));
    }
    // The curve thumbnails key by the SAME compound identity on the observations frame.
    const curveObs = observations === null ? null :
      keyedObservations(observations, compoundKeyCols, resultsDf, perProfile);
    // v4 (chemist-sar-grid S2): cross-profile %F per (compound × EV route), matched
    // to each EXTRAVASCULAR tier row (IV rows → null). The compound key is the same
    // (composite) key the tier uses, so the %F↔tier join is by identical strings.
    const {fPct: fPctByCompound, noIvRefCompounds} =
      computeFpctByCompound(fpctFrame, recipe, compoundColName);

    // v5 (chemist-sar-grid S4): inline curve thumbnails, keyed by the SAME identity
    // the tier row uses — per (compound × ROUTE) when the analysis carries a route
    // column (so IV and PO arms of a compound get SEPARATE curves, never pooled),
    // else per compound (single-route). Built once here; each tier row looks up its
    // own key below.
    const {routeKeyed: curvesRouteKeyed, byKey: curvesByKey, crossDoseOuters: curvesCrossDose} =
      buildCurvesByCompound(curveObs, compoundColName, label);

    // E4 (study-grid-grouping-columns): per-tier-row candidate-column inventory, keyed by
    // the SAME composite key the aggregate groups on. Each tier row looks up its own
    // subpopulation's inventory below (per-subpopulation scope, R2 finding 3).
    const inventoryByTier = buildKeyInventoryByTier(resultsDf, keys, recipe);

    // v10 (study-grid-passthrough-columns): per-tier-row inventory of the user-selected
    // passthrough (annotation) columns, keyed by the SAME composite tier key. Cap-free +
    // always-shown (RG-D) — SEPARATE candidate list from E4's, shared accumulator. Empty map
    // when the recipe selected no passthrough columns (the common case).
    const passthroughByTier = buildPassthroughByTier(resultsDf, keys, recipe);

    for (let i = 0; i < agg.rowCount; i++) {
      // compoundId = the composite tuple of real compound-key values (compositeGroupKey,
      // U+0001), or the analysis label when whole-analysis / no compound value on this
      // row. isWholeAnalysis = the exact negation — the durable label-INDEPENDENT
      // sentinel the study grid + drill-in read (v7; generalized to a set at v8).
      const anyCompoundValue = compoundKeyAggCols.some((c) => !c.isNone(i));
      const isWholeAnalysis = compoundKeyAggCols.length === 0 || !anyCompoundValue;
      const compoundId = isWholeAnalysis ? label : compositeGroupKey(compoundKeyAggCols, i);
      // Per-key display projection (Variant-2, v8): one {key: friendlyName, value} per
      // real compound key, in pick order; empty for a whole-analysis row. The study
      // grid renders one column per key; the drill-in zips these against groupByCols.
      // Dedup by the SAME typographic canonical id the study grid merges key columns on
      // (`normalizeGroupKey`, first-seen wins): two compound-key columns can carry the SAME
      // `friendlyName` (a casefold/whitespace collision), which the grid collapses into ONE key
      // column — so the persisted projection must too, else the identity label + curve `group`
      // print the value twice (`… · 1 · F · 1 · …`). Display-only: `compoundId`/`curveKey` below
      // are computed independently from the raw aggregate columns, NOT from `groupValues`.
      const groupValues: GroupKeyValue[] = isWholeAnalysis ? [] :
        dedupGroupValuesByCanonical(compoundKeyCols.map((k) => {
          const c = agg.col(k)!;
          const keyLabel = resultsDf.col(k)?.getTag('friendlyName') ?? k;
          return {key: keyLabel, value: c.isNone(i) ? '' : c.getString(i)};
        }));
      const route = (routeCol !== null && !routeCol.isNone(i)) ? routeCol.getString(i) : '';
      const iv = isIntravascular(route);
      // %F: null on IV rows (F ≡ 1 by definition); for an EV row, the compound's
      // EV-strata %F (NE sentinel where it has no IV reference). The join key uses
      // the compound VALUE only when the analysis was grouped by a compound column
      // (else '' — %F was matched by route within the whole analysis).
      const compoundKeyForF = compoundColName !== null ? compoundId : '';
      const fPct = iv ? null : (fPctByCompound.get(`${compoundKeyForF}${KEY_SEP}${route}`) ?? NA_MEANCV);

      // Per-group warnings, PLUS the honest NE explanation on an EV row whose
      // compound had a no-IV-reference stratum (blank %F → explained on the Flags
      // tooltip, not silently dropped — rule 14). Only for EV rows (IV %F is null
      // by definition, not "not estimable").
      const aggKey = keyStr((k) => agg.col(k), i);
      const nEx = excludedByGroup.get(aggKey) ?? 0;
      const wmap = perGroupWarnings.get(aggKey) ?? new Map<string, string>();
      const warnings: CompoundWarning[] = Array.from(wmap, ([code, message]) => ({code, message}));
      if (!iv && noIvRefCompounds.has(compoundKeyForF) && !warnings.some((w) => w.code === NO_IV_REFERENCE_CODE))
        warnings.push({code: NO_IV_REFERENCE_CODE, message: NO_IV_REFERENCE_MSG});

      // E3: this row's pooled thumbnail (one series), keyed per (compound × ROUTE) when
      // route-split so it is its OWN route's curve. Relabel the single series' `group` to
      // the row's group label (its `groupValues` joined by ` · `, or `''` for a
      // whole-analysis row) — the SAME label the v8→v9 migration writes on a preserved
      // curve, so fresh and migrated agree and no stale dose value leaks. Disclose a
      // ≥2-dose-cluster pool via CROSS_DOSE_POOL (BUG-08 honesty by disclosure).
      const curveKey = curvesRouteKeyed ? `${compoundId}${KEY_SEP}${route}` : compoundId;
      const rawCurve = curvesByKey.get(curveKey) ?? null;
      const groupLabel = groupValues.map((g) => g.value).join(' · ');
      const curve: CompoundCurve | null = rawCurve === null ? null :
        {series: rawCurve.series.map((s) => ({...s, group: groupLabel}))};
      if (curvesCrossDose.has(curveKey) && !warnings.some((w) => w.code === CROSS_DOSE_POOL_CODE))
        warnings.push({code: CROSS_DOSE_POOL_CODE, message: CROSS_DOSE_POOL_MSG});

      compoundTier.push({
        compoundId,
        isWholeAnalysis,
        groupValues,
        // E4: this row's per-subpopulation candidate-column inventory (present-but-ungrouped
        // keys → the grid's counted-join disclosure). Keyed by the SAME composite key.
        keyInventory: inventoryByTier.get(aggKey) ?? [],
        // v10: this row's per-subpopulation passthrough (annotation) inventory (cap-free,
        // always-shown; structures render bare-when-unique). Keyed by the SAME composite key.
        passthrough: passthroughByTier.get(aggKey) ?? [],
        n: (nCol !== null && !nCol.isNone(i)) ? (nCol.get(i) as number) : 0,
        route,
        doseNorm: (doseNormCol !== null && !doseNormCol.isNone(i)) ?
          doseNormCol.getString(i) : 'absolute',
        cmax: meanCvAt(agg, 'CMAX', i),
        // v5 (chemist-sar-grid S4): Tmax carry-through (aggregate already present)
        // for the row-click drill-in scalar card — not a study-grid column.
        tmax: meanCvAt(agg, 'TMAX', i),
        aucInf: meanCvAt(agg, 'AUC_INF', i, nEx),
        thalf: meanCvAt(agg, 'HALF_LIFE', i, nEx),
        cl: meanCvAt(agg, 'CL', i, nEx),
        vz: meanCvAt(agg, 'VZ', i, nEx),
        // v3 (chemist-sar-grid): MRT + Vss are already aggregated (DEFAULT_AGG_PARAMS
        // includes 'MRT','VSS'), so this is a pure carry-through — MEAN_MRT/CV_PCT_MRT/
        // N_MRT and MEAN_VSS/… already exist on `agg`. Vss is NaN for non-IV (sci-comp
        // gate) → meanCvAt yields NaN mean → the study grid drops it for EV cases.
        mrt: meanCvAt(agg, 'MRT', i, nEx),
        vss: meanCvAt(agg, 'VSS', i, nEx),
        // v4 (chemist-sar-grid S2): %F lead-trio leg (cross-profile, persisted —
        // the study grid sees only the cache). null on IV rows / no pair.
        fPct,
        pctExtrap: meanCvAt(agg, 'PCT_EXTRAP', i, nEx),
        // v5 (chemist-sar-grid S4) + E3 (study-grid-grouping-columns): this row's semi-log
        // geomean thumbnail — now ONE pooled series (the inner dose-split is dropped),
        // relabelled to the row's group. Keyed per (compound × ROUTE) when route-split, so
        // this row's curve is its OWN route's — never a route-mixed pool.
        curve,
        warnings,
      });
    }
  }

  // nCompounds = distinct compound identities (NOT row count — a mixed-route
  // compound contributes ≥2 rows but is one compound).
  const nCompounds = new Set(compoundTier.map((r) => r.compoundId)).size;

  return {
    analysisTier: {
      nCompounds,
      nSubjects,
      route,
      dose,
      status,
      flagSummary,
    },
    compoundTier,
  };
}
