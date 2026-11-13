/**
 * Dose proportionality (Smith power model) — cross-dose regression primitive
 * (FR-410..413 / UC-03, P0).
 *
 * **Why a NEW module, not `aggregateNcaResults` (mirror of `bioavailability.ts`
 * / `toxicokinetics.ts`).** The native aggregator does *within-group* arithmetic
 * (mean / SD / %CV / N) via `groupBy().aggregate()`; it cannot fit a slope.
 * Dose proportionality is intrinsically a *cross-dose* regression — `ln(AUC) =
 * α + β·ln(Dose)` across every subject of a Route×Sex stratum — so it is a
 * separate primitive. `aggregateNcaResults` is untouched.
 *
 * **Correct layering (not vendored math — Rule 5 / AC-9).** The OLS slope + CI
 * and the ANOVA/ANCOVA come from sci-comp (`stats.linearFit`, `stats.runAncova`,
 * `stats.oneWayAnova`). This module owns only the *domain orchestration*: design
 * assembly, the Smith critical region, verdict-mode classification, the
 * minimum-N advisory, completeness/NA handling, and the warning stream. No
 * normal equations or t-quantile are re-derived here.
 *
 * **Scientific contract (the honesty core — thesis #4):**
 *  - **Individual-subject fit is the default** (df = N − 2); the group-mean path
 *    (df = K − 2) is a *labelled* fallback only for composite/no-replicate
 *    designs (AD-1). The 03 oracle proves a group-mean fit at 3 doses is
 *    structurally inconclusive (df = 1) — it is never the default.
 *  - **Three verdict modes, never conflated** (AD-6): *conclusive*
 *    (proportional / sub / supra), *df-driven inconclusive* (too few levels /
 *    composite — no per-subject AUC), *IIV-driven inconclusive* (adequate df but
 *    between-animal variability widens the CI past the acceptance region).
 *  - **< 3 dose levels → no verdict** (`SMITH_INAPPLICABLE`, AD-3): a 2-dose
 *    point estimate with no CI is not a proportionality call.
 *  - **The critical region is recomputed from the realized dose ratio and θ**
 *    (AD-5), never frozen: `CR = [1 + ln(θL)/ln(r), 1 + ln(θH)/ln(r)]`.
 *  - **DOSE_NORM design guard** (AD-10): a stratum mixing absolute / per-kg dose
 *    refuses the fit (`MIXED_DOSE_NORM_FIT`) rather than fit a meaningless slope.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import {f64Column, numAt} from './adapter-utils';
import {occasionNumber} from './occasion';
import {stats} from '@datagrok-libraries/sci-comp';

// --- Defaults (cited, conservative) -------------------------------------

/** Smith equivalence bounds (the FDA/EMA average-bioequivalence interval). */
export const DEFAULT_THETA_LOW = 0.8;
export const DEFAULT_THETA_HIGH = 1.25;
/** Smith CI level = 1 − 2α, α = 0.05 → 90% two-sided (the equivalence default). */
export const DEFAULT_CI_LEVEL = 0.90;
/** Practical upper bound for the minimum-N root-find (animals per group). A
 *  balanced design always has a *mathematically* feasible n (Sxx(n) = n·C → ∞),
 *  so the advisory is a *practical* search; no n ≤ cap ⇒ `NO_FEASIBLE_N` (AD-6,
 *  SY-F1 fix). */
export const DEFAULT_MIN_N_SEARCH_CAP = 50;
/** Within-dose-group residual-SD max/min ratio that flags heteroscedasticity
 *  (AD-8; a descriptive threshold — a formal Bartlett/Levene test has near-zero
 *  power at n = 4/group). Calibrated against the reference data at build. */
export const HETEROSCEDASTICITY_SD_RATIO = 4;
/** Realized-vs-balanced Sxx floor below which an NA drop is flagged as biased
 *  toward the dose extremes (AD-6, R2-F2). */
export const EXTREME_DOSE_DROP_SXX_RATIO = 0.75;
/** Minimum completeness for an individual-fit Smith verdict: ≥ this many dose
 *  levels, each retaining ≥ {@link MIN_ANIMALS_PER_LEVEL} animals (AD-6). */
export const MIN_DOSE_LEVELS = 3;
export const MIN_ANIMALS_PER_LEVEL = 2;
/** Sex×dose slope-homogeneity test level (AD-7). The test itself — an ANCOVA
 *  homogeneity-of-regression-slopes F-test fit as `ln(AUC) ~ ln(Dose) + Sex +
 *  ln(Dose):Sex`, where the interaction term's significance IS the test —
 *  follows Huitema, *The Analysis of Covariance and Alternatives*, 2nd ed.
 *  (Wiley, 2011), §6.6 pp. 144–148 and the general-linear-model form §7.6
 *  p. 178; heterogeneous slopes → per-stratum fits (Huitema ch. 11), the AD-7
 *  fallback. The 0.05→0.10 relaxation is a model-selection convention, not a
 *  confirmatory level: Bendel & Afifi, "Comparison of Stopping Rules in Forward
 *  Stepwise Regression", JASA 72(357):46–53 (1977), §6 p. 52 + Table 3 (0.05 is
 *  too stringent for term retention; recommended 0.10–0.25, with α=0.10 the
 *  backward-elimination optimum). α = 0.10 RATIFIED by USER-DECISION 2026-06-14;
 *  citation pinned 2026-06-26 (DATA-GAP-19). NB Huitema's worked example tests
 *  this assumption at 0.05 and neither source prescribes 0.10 for *this* test —
 *  0.10 is the relaxed model-selection convention applied to the interaction-
 *  retention decision, not a value lifted verbatim from one paper. (The earlier
 *  Gough-1995/Hummel-2009 plan was a mis-citation: those papers ground the power
 *  model and the WLS/heteroscedasticity remedy, not an interaction-test α.) */
export const SEX_INTERACTION_ALPHA = 0.10;

// --- Public types --------------------------------------------------------

/** AUC integration basis for the headline fit (AD-4). */
export type AucBasis = 'aucinf' | 'auclast' | 'auctau';

/** The parameter identity shown in the verdict table's `PARAMETER` column in per-occasion
 *  mode (guardrail "one basis per column": the basis IS the parameter identity, never a
 *  row label). */
export const DP_PARAMETER_ID: Readonly<Record<AucBasis | 'cmax', string>> = {
  aucinf: 'AUC∞', auctau: 'AUC0–τ', auclast: 'AUClast', cmax: 'Cmax',
};

/** The app's λz acceptance limit on %AUCextrap, reused as the Day-1 AUC∞ validity limit (F4):
 *  "valid AUC∞" means the same on the results grid and in dose proportionality. */
export const AUCINF_VALID_EXTRAP_PCT = 20;

/** Note on a verdict row whose class differs on the other printed basis (thesis #4). */
export const VERDICT_DIFFERS_NOTE = 'verdict differs on the other basis';

/** The honest-uncertainty mode of a verdict (AD-6). */
export type VerdictMode = 'conclusive' | 'df-driven' | 'iiv-driven' | 'inapplicable';

/** The dose-proportionality call for one stratum × parameter. */
export type DoseProportionalityVerdict =
  | 'proportional'
  | 'sub-proportional'
  | 'supra-proportional'
  | 'inconclusive-too-few-levels' // df-driven: remedy is more dose levels
  | 'inconclusive-composite' // df-driven: remedy is serial sampling
  | 'inconclusive-underpowered' // iiv-driven: remedy is more animals / lower IIV
  | 'insufficient-completeness' // dn-only; no fit
  | 'inapplicable'; // < 3 dose levels / refused

export type DoseProportionalityWarningCode =
  | 'SMITH_INAPPLICABLE'
  | 'INSUFFICIENT_COMPLETENESS'
  | 'EXTREME_DOSE_DROP'
  | 'AUCINF_NA_STARVED'
  | 'MIXED_DOSE_NORM_FIT'
  | 'HETEROSCEDASTICITY_SUSPECTED'
  | 'SEX_DOSE_CONFOUNDED'
  | 'SEX_STRATUM_INSUFFICIENT'
  | 'SEX_BASIS_DIFFERS'
  | 'NO_FEASIBLE_N'
  | 'MIN_N_ADVISORY_UNAVAILABLE';

export interface DoseProportionalityWarning {
  readonly code: DoseProportionalityWarningCode;
  readonly message: string;
  /** The stratum key the warning applies to (e.g. `'ROUTE=IV|SEX=F'`), or '' for
   *  a table-level warning. */
  readonly stratum: string;
}

/**
 * Short, human-readable label for a dose-proportionality advisory code — for the
 * summary toast (the transient balloon can't hold the per-stratum sentences,
 * which repeat the same code across strata). The full text stays in
 * `result.warnings`; the verdict grid carries the per-stratum VERDICT/MODE.
 */
export function shortDoseProportionalityNote(code: string): string {
  switch (code) {
  case 'SMITH_INAPPLICABLE': return 'too few dose levels';
  case 'INSUFFICIENT_COMPLETENESS': return 'insufficient completeness';
  case 'EXTREME_DOSE_DROP': return 'NA drop biased to dose extremes';
  case 'AUCINF_NA_STARVED': return 'AUCinf NA (try AUClast basis)';
  case 'MIXED_DOSE_NORM_FIT': return 'mixed dose normalization';
  case 'HETEROSCEDASTICITY_SUSPECTED': return 'heteroscedasticity';
  case 'SEX_DOSE_CONFOUNDED': return 'sex × dose confounded';
  case 'SEX_STRATUM_INSUFFICIENT': return 'sex stratum insufficient';
  case 'SEX_BASIS_DIFFERS': return 'sexes on different AUC bases';
  case 'NO_FEASIBLE_N': return 'underpowered (no feasible N)';
  case 'MIN_N_ADVISORY_UNAVAILABLE': return 'min-N advisory N/A';
  default: return code;
  }
}

/** Minimum-N study-design advisory for an underpowered (IIV-driven) stratum. */
export interface MinNAdvisory {
  /** Whether a conclusive verdict is reachable at n ≤ search cap. */
  readonly feasible: boolean;
  /** Smallest animals-per-group that makes the slope CI fit inside the CR; null
   *  when not feasible (≤ cap) or not applicable (unbalanced / non-IIV). */
  readonly nStar: number | null;
  /** Residual σ on log-AUC that drives the advisory. */
  readonly sigma: number;
  /** True when the advisory could not run (unbalanced design — AD-6 scope). */
  readonly unavailable: boolean;
}

/** Per-stratum × parameter structured fit (feeds the table and the P4 charts). */
export interface DoseProportionalityFit {
  readonly route: string;
  readonly sex: string;
  /** The occasion (study day / period) of this fit; '' when the analysis is not per occasion. */
  readonly occasion: string;
  readonly parameter: 'AUC' | 'Cmax';
  readonly basis: AucBasis | null; // null for Cmax
  /** Per-occasion mode: true when this row carries the stratum's verdict; false when it is
   *  shown beside it (β + CI only). Always true outside per-occasion mode. */
  readonly carriesVerdict: boolean;
  /** Per-dose-level N in ascending dose order, e.g. `3/3/2` (attrition disclosure, R2 F7). */
  readonly nByDose: string;
  /** Free-text notes for the row (basis fallback, verdict differs, attrition). */
  readonly notes: ReadonlyArray<string>;
  readonly beta: number;
  readonly intercept: number;
  readonly slopeCI: readonly [number, number];
  readonly rSquared: number;
  readonly doseRatio: number; // realized r = max/min
  readonly criticalRegion: readonly [number, number];
  readonly verdict: DoseProportionalityVerdict;
  readonly mode: VerdictMode;
  readonly n: number; // subjects in the fit
  readonly df: number;
  readonly nDoseLevels: number;
  readonly isGroupMeanFallback: boolean;
  readonly minNAdvisory: MinNAdvisory | null;
  /** Raw per-subject points for the log-log scatter (AD-11). */
  readonly points: ReadonlyArray<{readonly dose: number; readonly value: number}>;
}

/**
 * Secondary dose-normalized-AUC one-way ANOVA across dose groups (AD-9 / FR-412).
 * **Offered, not led** — a non-significant result is *absence of evidence*, not
 * proof of proportionality; the equivalence-framed power-model CI is primary.
 */
export interface DnAucAnova {
  readonly route: string;
  readonly sex: string;
  readonly occasion: string;
  readonly fStatistic: number | null;
  readonly pValue: number | null;
  readonly nGroups: number;
  readonly n: number;
}

/** Per-route sex slope-homogeneity result (AD-7). */
export interface SexHomogeneity {
  readonly route: string;
  readonly occasion: string;
  readonly fStatistic: number | null;
  readonly pValue: number | null;
  readonly homogeneous: boolean;
  /** Pooled common dose slope (only when homogeneous), df = N − 3. */
  readonly pooledSlope: number | null;
  readonly pooledSlopeCI: readonly [number, number] | null;
  readonly pooledDf: number | null;
}

export interface DoseProportionalityOptions {
  /** AUCinf column. Default `'AUC_INF'`. */
  readonly aucInfCol?: string;
  /** AUClast column. Default `'AUC_LAST'`. */
  readonly aucLastCol?: string;
  /** Cmax column. Default `'CMAX'`. */
  readonly cmaxCol?: string;
  /** Per-subject dose column (carried metadata, e.g. `'Dose_mg_kg'`). Required. */
  readonly doseCol: string;
  /** Route column. Default `'ROUTE'`. When absent the analysis collapses to a
   *  single route stratum. */
  readonly routeCol?: string;
  /** Sex column (mandatory stratifier for the headline fit). Default `'SEX'`. */
  readonly sexCol?: string;
  /** AUC integration basis. Default `'aucinf'` (AD-4). */
  readonly basis?: AucBasis;
  /** Smith equivalence bounds. Default 0.8 / 1.25 (AD-5). */
  readonly thetaLow?: number;
  readonly thetaHigh?: number;
  /** Two-sided CI level. Default 0.90 (AD-5). */
  readonly ciLevel?: number;
  /** Practical cap for the minimum-N root-find. Default 50 (AD-6). */
  readonly minNSearchCap?: number;
  /**
   * Per-occasion mode (tk-study-day-stratification F2): the occasion column of `results`
   * (one value per profile). When set with ≥ 2 levels, every stratum is route × sex ×
   * occasion and the AUC basis follows the occasion rule (F4) instead of `basis`. Absent →
   * the single-basis behaviour above, unchanged.
   */
  readonly occasionCol?: string;
  /** Occasion levels in display order (from `resolveOccasion`); the first is "Day 1". */
  readonly occasionLevels?: ReadonlyArray<string>;
  /** AUC0–τ column (see `appendAucTauColumn`); absent → no AUC0–τ basis, AUClast fallback. */
  readonly aucTauCol?: string;
  /** %AUCextrap column for the Day-1 AUC∞ validity guard. Default `'PCT_EXTRAP'`. */
  readonly pctExtrapCol?: string;
  /** Animal per results row (from `resolveOccasion`); null → no animal key. Drives the
   *  paired Δβ and the DN table's Animal column. */
  readonly animals?: ReadonlyArray<string> | null;
  /** Animals excluded from per-animal outputs (failing verification). */
  readonly excludedAnimals?: ReadonlySet<string>;
  /** Subject column of `results`. Default `'SUBJECT'`. */
  readonly subjectCol?: string;
  /** TLAST column, used to compare sampling spans for Δβ suppression and for the Day-1 AUC∞
   *  guard's λz ≤ τ test. Default `'TLAST'`. */
  readonly tlastCol?: string;
  /** The declared dosing interval τ (h); NaN / absent = unknown. */
  readonly tauH?: number;
  /** Every occasion is a single dose (a confirmed `crossover`): no AUC0–τ, the first-occasion
   *  AUC rule on every occasion, and Δβ is labelled a period effect. Default false. */
  readonly singleDoseOccasions?: boolean;
}

/** One time-dependence row: Δβ = β(later) − β(first) on one common basis (F3a). */
export interface DeltaBetaRow {
  readonly route: string;
  readonly sex: string;
  readonly first: string;
  readonly later: string;
  /** The common basis both slopes use, or null when suppressed. */
  readonly basis: AucBasis | null;
  readonly method: 'paired' | 'two-fits' | 'suppressed';
  readonly deltaBeta: number;
  readonly ci: readonly [number, number];
  readonly df: number;
  readonly nFirst: number;
  readonly nLater: number;
  readonly nPaired: number;
  /** N per dose group on each occasion, ascending dose (`3/3/3`) — attrition disclosure (R2 F7). */
  readonly nByDoseFirst: string;
  readonly nByDoseLater: string;
  readonly notes: ReadonlyArray<string>;
}

export interface DoseProportionalityResult {
  /** One row per {route, sex, parameter, basis}: β, CI, r, CR, verdict, mode,
   *  n, df, levels, R². */
  readonly table: DG.DataFrame;
  /** Structured per-stratum fits (feed the P4 charts). */
  readonly fits: DoseProportionalityFit[];
  /** Per-subject dose-normalized exposure (AD-9 bar). */
  readonly dnTable: DG.DataFrame;
  /** Per-route sex slope-homogeneity (when both sexes present & orthogonal). */
  readonly sexHomogeneity: SexHomogeneity[];
  /** Secondary dn-AUC one-way ANOVA per Route×Sex (AD-9 — offered, not led). */
  readonly dnAnova: DnAucAnova[];
  /** Active AUC basis (surfaced in the grid header). */
  readonly basis: AucBasis;
  readonly warnings: DoseProportionalityWarning[];
  /** Per-occasion mode only: Δβ rows (empty otherwise) and their table (null otherwise). */
  readonly deltaBeta: DeltaBetaRow[];
  readonly deltaTable: DG.DataFrame | null;
}

// --- Pure helpers (exported for unit tests) -----------------------------

/**
 * Smith critical region for the power-model slope at dose ratio `r`:
 * `[1 + ln(θL)/ln(r), 1 + ln(θH)/ln(r)]` (AD-5). NaN pair when `r ≤ 1`.
 */
export function computeCriticalRegion(
  r: number, thetaLow: number, thetaHigh: number,
): [number, number] {
  if (!(r > 1)) return [NaN, NaN];
  const lnr = Math.log(r);
  return [1 + Math.log(thetaLow) / lnr, 1 + Math.log(thetaHigh) / lnr];
}

/**
 * Classify the slope CI against the critical region. `conclusive` is true only
 * when the CI lies entirely inside (proportional) or entirely on one side
 * (sub / supra) of the CR.
 */
export function classifyAgainstCR(
  ci: readonly [number, number], cr: readonly [number, number],
): {conclusive: boolean; direction: 'within' | 'below' | 'above' | 'straddle'} {
  const [ciLo, ciHi] = ci;
  const [crLo, crHi] = cr;
  if (!Number.isFinite(ciLo) || !Number.isFinite(ciHi) || !Number.isFinite(crLo) || !Number.isFinite(crHi))
    return {conclusive: false, direction: 'straddle'};
  if (ciLo >= crLo && ciHi <= crHi) return {conclusive: true, direction: 'within'};
  if (ciHi < crLo) return {conclusive: true, direction: 'below'};
  if (ciLo > crHi) return {conclusive: true, direction: 'above'};
  return {conclusive: false, direction: 'straddle'};
}

/**
 * Minimum-N root-find (AD-6, balanced designs only). Smallest integer
 * `n ∈ [2, cap]` for which `t(0.95, nK−2)·σ/√(n·sxxPerAnimal) < crHalf`, where
 * `sxxPerAnimal = Σ_k (ln dose_k − mean ln dose)²` over the K distinct levels.
 * Returns the boundary (`feasible: false`) when no n ≤ cap satisfies it.
 */
export function minNAdvisory(
  sigma: number, sxxPerAnimal: number, nLevels: number, crHalf: number, cap: number,
): {feasible: boolean; nStar: number | null} {
  if (!(sigma > 0) || !(sxxPerAnimal > 0) || !(crHalf > 0) || nLevels < 2)
    return {feasible: false, nStar: null};
  for (let n = 2; n <= cap; n++) {
    const df = n * nLevels - 2;
    if (df <= 0) continue;
    const tCrit = stats.distributions.studentTInv(0.95, df);
    const ciHalf = tCrit * sigma / Math.sqrt(n * sxxPerAnimal);
    if (ciHalf < crHalf) return {feasible: true, nStar: n};
  }
  return {feasible: false, nStar: null};
}

// --- Internal helpers ----------------------------------------------------

function strAt(col: DG.Column | null, row: number): string {
  return col === null || col.isNone(row) ? '' : col.getString(row);
}

function sampleSd(values: ReadonlyArray<number>): number {
  const n = values.length;
  if (n < 2) return NaN;
  let sum = 0;
  for (const v of values) sum += v;
  const m = sum / n;
  let ss = 0;
  for (const v of values) {const d = v - m; ss += d * d;}
  return Math.sqrt(ss / (n - 1));
}

/** A subject's (dose, parameter value) point plus whether the value is usable. */
interface SubjectPoint {dose: number; value: number; usable: boolean; animal: string; tlast: number;}

/** N per dose level in ascending dose order (`3/3/2`), dose-0 rows excluded. */
function nByDoseOf(points: ReadonlyArray<{dose: number}>): string {
  const per = new Map<number, number>();
  for (const p of points) if (p.dose > 0) per.set(p.dose, (per.get(p.dose) ?? 0) + 1);
  return [...per.keys()].sort((a, b) => a - b).map((d) => per.get(d)).join('/');
}

/** One Route×Sex(×Occasion) × parameter design before fitting. */
interface StratumDesign {
  route: string;
  sex: string;
  occasion: string;
  parameter: 'AUC' | 'Cmax';
  basis: AucBasis | null;
  carriesVerdict: boolean;
  notes: string[];
  /** All subjects with a finite dose (usable flags whether the value is finite). */
  subjects: SubjectPoint[];
  /** Distinct DOSE_NORM bases seen over the stratum's rows. */
  doseNorms: Set<string>;
}

interface ParamSpec {parameter: 'AUC' | 'Cmax'; col: DG.Column | null; basis: AucBasis | null;}

// =====================================================================
// Fit pipeline for one Route×Sex × parameter design
// =====================================================================

function fitStratum(
  design: StratumDesign,
  options: Required<Pick<DoseProportionalityOptions,
    'thetaLow' | 'thetaHigh' | 'ciLevel' | 'minNSearchCap'>>,
  warningSink: DoseProportionalityWarning[],
): DoseProportionalityFit {
  const stratumKey = `ROUTE=${design.route}|SEX=${design.sex}` +
    (design.occasion !== '' ? `|OCCASION=${design.occasion}` : '');
  const paramKey = `${stratumKey}|${design.parameter}` +
    (design.occasion !== '' && design.basis !== null ? `(${DP_PARAMETER_ID[design.basis]})` : '');
  const thetaLow = options.thetaLow;
  const thetaHigh = options.thetaHigh;
  // A fit shown BESIDE the verdict row must not repeat the verdict row's warnings.
  const warnings: DoseProportionalityWarning[] = design.carriesVerdict ? warningSink : [];
  const common = {
    route: design.route, sex: design.sex, occasion: design.occasion, parameter: design.parameter,
    basis: design.basis, carriesVerdict: design.carriesVerdict, notes: design.notes,
  };

  const degenerate = (
    verdict: DoseProportionalityVerdict, mode: VerdictMode,
    nLevels: number, n: number,
  ): DoseProportionalityFit => ({
    ...common, nByDose: nByDoseOf(design.subjects),
    beta: NaN, intercept: NaN, slopeCI: [NaN, NaN], rSquared: NaN,
    doseRatio: NaN, criticalRegion: [NaN, NaN], verdict, mode,
    n, df: NaN as unknown as number, nDoseLevels: nLevels,
    isGroupMeanFallback: false, minNAdvisory: null, points: [],
  });

  // --- AD-10: DOSE_NORM design guard. Refuse a mixed-basis fit.
  if (design.doseNorms.size > 1) {
    warnings.push({code: 'MIXED_DOSE_NORM_FIT', stratum: paramKey,
      message: `${paramKey}: mixes absolute / per-kg dose normalization; ` +
        `the power-model fit is refused (a mixed-basis ln(Dose) slope is meaningless).`});
    return degenerate('inapplicable', 'inapplicable', 0, design.subjects.length);
  }

  const usable = design.subjects.filter((s) => s.usable && s.dose > 0 && Number.isFinite(s.value) && s.value > 0);
  const removed = design.subjects.filter((s) => s.dose > 0 && !(s.usable && Number.isFinite(s.value) && s.value > 0));

  // Distinct dose levels and per-level animal counts (after NA removal).
  const perLevel = new Map<number, number>();
  for (const s of usable) perLevel.set(s.dose, (perLevel.get(s.dose) ?? 0) + 1);
  const levels = [...perLevel.keys()].sort((a, b) => a - b);
  const K = levels.length;
  const N = usable.length;

  // --- AD-4 / AC-10: if the AUCinf basis is NA-starved (λz/AUCinf NA across most
  // subjects, e.g. steady-state disease data), recommend the AUClast basis. This
  // is checked BEFORE the degenerate-branch returns below, because a fully
  // NA-starved stratum collapses to < 3 usable levels and would otherwise exit
  // via SMITH_INAPPLICABLE with no AUClast hint.
  if (design.parameter === 'AUC' && design.basis === 'aucinf' &&
      removed.length > 0 && removed.length >= usable.length) {
    warnings.push({code: 'AUCINF_NA_STARVED', stratum: paramKey,
      message: `${paramKey}: AUCinf is NA across most subjects (${removed.length} dropped vs ` +
        `${usable.length} usable). Re-run on the AUClast basis if λz/AUCinf is widely unfit.`});
  }

  // --- AD-3: < 3 dose levels → withhold the Smith verdict.
  if (K < MIN_DOSE_LEVELS) {
    warnings.push({code: 'SMITH_INAPPLICABLE', stratum: paramKey,
      message: `≥ ${MIN_DOSE_LEVELS} dose levels are required for the Smith dose-proportionality ` +
        `criterion; ${paramKey} spans ${K} level(s). Showing dose-normalized exposure only.`});
    return degenerate('inapplicable', 'inapplicable', K, N);
  }

  const maxAnimals = Math.max(...perLevel.values());
  const isGroupMeanFallback = maxAnimals === 1; // one value per dose level → composite

  // --- AD-6: completeness rule (individual path only; composite is exempt).
  if (!isGroupMeanFallback) {
    const levelsWithEnough = levels.filter((d) => (perLevel.get(d) ?? 0) >= MIN_ANIMALS_PER_LEVEL).length;
    if (levelsWithEnough < MIN_DOSE_LEVELS) {
      // (AUCINF_NA_STARVED — the "recommend AUClast" hint — is raised by the
      // NA-starvation pre-check above, which fires only when AUCinf NA is the
      // actual cause; a genuinely sparse design with valid AUCinf gets no
      // misleading AUClast suggestion.)
      const range = K > 0 ? `${levels[0]}–${levels[K - 1]}` : 'n/a';
      warnings.push({code: 'INSUFFICIENT_COMPLETENESS', stratum: paramKey,
        message: `${paramKey}: after NA removal, fewer than ${MIN_DOSE_LEVELS} dose levels ` +
          `retain ≥ ${MIN_ANIMALS_PER_LEVEL} animals (realized: ${K} levels, N=${N}, ` +
          `dose range ${range}). No Smith verdict; dose-normalized exposure only.`});
      return degenerate('insufficient-completeness', 'df-driven', K, N);
    }
  }

  // --- The fit: ln(value) ~ ln(dose). Individual subjects (or composite means).
  const lnDose: number[] = [];
  const lnVal: number[] = [];
  for (const s of usable) {lnDose.push(Math.log(s.dose)); lnVal.push(Math.log(s.value));}
  const fit = stats.linearFit(lnDose, lnVal, {ciLevel: options.ciLevel});

  const r = levels[K - 1] / levels[0];
  const cr = computeCriticalRegion(r, thetaLow, thetaHigh);
  const ci: [number, number] = [fit.slopeCI[0], fit.slopeCI[1]];
  const cls = classifyAgainstCR(ci, cr);

  // --- Verdict + mode.
  let verdict: DoseProportionalityVerdict;
  let mode: VerdictMode;
  if (cls.conclusive) {
    mode = 'conclusive';
    verdict = cls.direction === 'within' ? 'proportional' :
      cls.direction === 'below' ? 'sub-proportional' : 'supra-proportional';
  } else if (isGroupMeanFallback) {
    mode = 'df-driven';
    verdict = 'inconclusive-composite';
  } else {
    mode = 'iiv-driven';
    verdict = 'inconclusive-underpowered';
  }

  // --- AD-6: EXTREME_DOSE_DROP — flag a removal biased toward the dose extremes.
  if (removed.length > 0) flagExtremeDoseDrop(usable, removed, levels, paramKey, warnings);

  // --- AD-8: heteroscedasticity detection (suppressed if any level < 3 residuals).
  if (!isGroupMeanFallback)
    detectHeteroscedasticity(usable, fit.intercept, fit.slope, paramKey, warnings);

  // --- AD-6: minimum-N advisory (IIV-driven only; balanced designs only).
  let advisory: MinNAdvisory | null = null;
  if (mode === 'iiv-driven') {
    advisory = buildMinNAdvisory(
      usable, levels, perLevel, fit, cr, options.minNSearchCap, paramKey, warnings);
  }

  return {
    ...common, nByDose: nByDoseOf(usable),
    beta: fit.slope, intercept: fit.intercept, slopeCI: ci, rSquared: fit.rSquared,
    doseRatio: r, criticalRegion: cr, verdict, mode,
    n: N, df: fit.df, nDoseLevels: K, isGroupMeanFallback, minNAdvisory: advisory,
    points: usable.map((s) => ({dose: s.dose, value: s.value})),
  };
}

function flagExtremeDoseDrop(
  usable: SubjectPoint[], removed: SubjectPoint[], retainedLevels: number[],
  paramKey: string, warnings: DoseProportionalityWarning[],
): void {
  // (a) Realized vs balanced-reallocation Sxx on retained subjects.
  const lnDoses = usable.map((s) => Math.log(s.dose));
  const meanLn = lnDoses.reduce((a, b) => a + b, 0) / lnDoses.length;
  let sxxRealized = 0;
  for (const l of lnDoses) {const d = l - meanLn; sxxRealized += d * d;}
  const lnLevels = retainedLevels.map((d) => Math.log(d));
  const meanLnLevel = lnLevels.reduce((a, b) => a + b, 0) / lnLevels.length;
  let levelSpread = 0;
  for (const l of lnLevels) {const d = l - meanLnLevel; levelSpread += d * d;}
  const sxxBalancedRef = (usable.length / retainedLevels.length) * levelSpread;
  const ratio = sxxBalancedRef > 0 ? sxxRealized / sxxBalancedRef : 1;

  // (b) ≥ 50% of removed subjects at the (full-design) min or max dose.
  const allDoses = [...usable, ...removed].map((s) => s.dose);
  const minDose = Math.min(...allDoses);
  const maxDose = Math.max(...allDoses);
  const removedAtExtreme = removed.filter((s) => s.dose === minDose || s.dose === maxDose).length;
  const extremeFraction = removed.length > 0 ? removedAtExtreme / removed.length : 0;

  if (ratio < EXTREME_DOSE_DROP_SXX_RATIO || extremeFraction >= 0.5) {
    warnings.push({code: 'EXTREME_DOSE_DROP', stratum: paramKey,
      message: `${paramKey}: NA removal is biased toward the dose extremes ` +
        `(retained Sxx ${(ratio * 100).toFixed(0)}% of the balanced reference; ` +
        `${(extremeFraction * 100).toFixed(0)}% of dropped subjects at min/max dose). ` +
        `The realized log-dose leverage is degraded — interpret the slope CI with care.`});
  }
}

function detectHeteroscedasticity(
  usable: SubjectPoint[], intercept: number, slope: number,
  paramKey: string, warnings: DoseProportionalityWarning[],
): void {
  const byLevel = new Map<number, number[]>();
  for (const s of usable) {
    const resid = Math.log(s.value) - (intercept + slope * Math.log(s.dose));
    let arr = byLevel.get(s.dose);
    if (arr === undefined) {arr = []; byLevel.set(s.dose, arr);}
    arr.push(resid);
  }
  // Degenerate-SD guard (R2-F3): every group must keep ≥ 3 residuals.
  const sds: number[] = [];
  for (const arr of byLevel.values()) {
    if (arr.length < 3) return; // suppress — SD too unstable to act on
    sds.push(sampleSd(arr));
  }
  const finite = sds.filter((v) => Number.isFinite(v) && v > 0);
  if (finite.length < 2) return;
  const ratio = Math.max(...finite) / Math.min(...finite);
  if (ratio >= HETEROSCEDASTICITY_SD_RATIO) {
    warnings.push({code: 'HETEROSCEDASTICITY_SUSPECTED', stratum: paramKey,
      message: `${paramKey}: within-dose-group residual-SD ratio ${ratio.toFixed(1)} ` +
        `≥ ${HETEROSCEDASTICITY_SD_RATIO} — variance is not constant across dose. ` +
        `The OLS slope CI may be optimistic; a weighted (WLS) fit is the remedy.`});
  }
}

function buildMinNAdvisory(
  usable: SubjectPoint[], levels: number[], perLevel: Map<number, number>,
  fit: stats.LinearFitResult, cr: readonly [number, number], cap: number,
  paramKey: string, warnings: DoseProportionalityWarning[],
): MinNAdvisory {
  // Balanced check: every retained dose level carries the same animal count.
  const counts = levels.map((d) => perLevel.get(d) ?? 0);
  const balanced = counts.every((c) => c === counts[0]);
  if (!balanced) {
    warnings.push({code: 'MIN_N_ADVISORY_UNAVAILABLE', stratum: paramKey,
      message: `${paramKey}: the study-design (minimum-N) advisory is available for ` +
        `balanced designs only; this stratum has unequal animals per dose level.`});
    return {feasible: false, nStar: null, sigma: NaN, unavailable: true};
  }
  // σ_IIV = residual σ on log-AUC, derived from slopeSe and Sxx (σ = se·√Sxx).
  const lnDoses = usable.map((s) => Math.log(s.dose));
  const meanLn = lnDoses.reduce((a, b) => a + b, 0) / lnDoses.length;
  let sxxRealized = 0;
  for (const l of lnDoses) {const d = l - meanLn; sxxRealized += d * d;}
  const sigma = fit.slopeSe * Math.sqrt(sxxRealized);

  // Sxx per added animal = Σ_k (ln dose_k − mean ln dose_k)² over distinct levels.
  const lnLevels = levels.map((d) => Math.log(d));
  const meanLnLevel = lnLevels.reduce((a, b) => a + b, 0) / lnLevels.length;
  let sxxPerAnimal = 0;
  for (const l of lnLevels) {const d = l - meanLnLevel; sxxPerAnimal += d * d;}

  const crHalf = (cr[1] - cr[0]) / 2;
  const {feasible, nStar} = minNAdvisory(sigma, sxxPerAnimal, levels.length, crHalf, cap);
  if (!feasible) {
    warnings.push({code: 'NO_FEASIBLE_N', stratum: paramKey,
      message: `${paramKey}: underpowered by between-animal variability (σ = ${sigma.toFixed(3)} ` +
        `on log-AUC); a conclusive verdict needs impractically many animals (> ${cap}/group) — ` +
        `lower IIV, not a wider dose range.`});
  }
  return {feasible, nStar, sigma, unavailable: false};
}

// =====================================================================
// AD-7 — Sex slope-homogeneity via runAncova (per route)
// =====================================================================

function computeSexHomogeneity(
  byRoute: Map<string, StratumDesign[]>,
  ciLevel: number,
  warnings: DoseProportionalityWarning[],
): SexHomogeneity[] {
  const out: SexHomogeneity[] = [];
  for (const designs of byRoute.values()) {
    const route = designs[0].route;
    const occasion = designs[0].occasion;
    // Per occasion the same sentence repeats per day, so the stratum names the day.
    const stratum = occasion !== '' ? `ROUTE=${route}, ${occasion}` : `ROUTE=${route}`;
    // Pool the verdict-carrying AUC designs of this route (× occasion) across sexes.
    const aucDesigns = designs.filter((d) => d.parameter === 'AUC' && d.carriesVerdict);
    const lnDose: number[] = [];
    const lnAuc: number[] = [];
    const sexCode: number[] = [];
    const poolDose: number[] = [];
    const sexes = new Set<string>();
    for (const d of aucDesigns) {
      const code = d.sex.toUpperCase().startsWith('F') ? 0 : d.sex.toUpperCase().startsWith('M') ? 1 : -1;
      for (const s of d.subjects) {
        if (!s.usable || !(s.dose > 0) || !Number.isFinite(s.value) || !(s.value > 0)) continue;
        sexes.add(d.sex);
        if (code < 0) continue;
        lnDose.push(Math.log(s.dose));
        lnAuc.push(Math.log(s.value));
        sexCode.push(code);
        poolDose.push(s.dose);
      }
    }
    // AD-7 (SY-F2): exactly two sex categories, deterministically encoded F→0/M→1.
    const codeSet = new Set(sexCode);
    if (sexes.size !== 2 || codeSet.size !== 2) {
      if (sexes.size > 1) {
        warnings.push({code: 'SEX_STRATUM_INSUFFICIENT', stratum,
          message: `${stratum}: sex slope-homogeneity needs exactly two encodable sex ` +
            `categories (F/M); found ${sexes.size}. Per-stratum fits stand; no sex-interaction claim.`});
      }
      continue;
    }
    // Per occasion the AUC basis carrying the verdict is decided per sex (F4), so the two sexes
    // can carry different AUCs (e.g. F on AUC∞, M on AUC0–τ). A sex×dose interaction across two
    // AUC definitions would read a basis difference as a sex effect: withhold the claim.
    const bases = new Set(aucDesigns.map((d) => d.basis));
    if (bases.size > 1) {
      warnings.push({code: 'SEX_BASIS_DIFFERS', stratum,
        message: `${stratum}: the sexes carry the verdict on different AUC bases ` +
          `(${[...bases].map((b) => DP_PARAMETER_ID[b ?? 'cmax']).join(' vs ')}); ` +
          'per-stratum fits stand; no sex-interaction claim.'});
      continue;
    }
    // AD-7 F8 guard: the interaction model is valid only when Sex ⊥ Dose. A
    // sex-stratified dose-level design (a dose level carrying only one sex)
    // confounds sex-slope heterogeneity with dose-level imbalance — the ANCOVA
    // interaction term aliases the dose effect. Withhold the sex claim.
    const sexesAtDose = new Map<number, Set<number>>();
    for (let i = 0; i < poolDose.length; i++) {
      let set = sexesAtDose.get(poolDose[i]);
      if (set === undefined) {set = new Set(); sexesAtDose.set(poolDose[i], set);}
      set.add(sexCode[i]);
    }
    const confounded = [...sexesAtDose.values()].some((s) => s.size < 2);
    if (confounded) {
      warnings.push({code: 'SEX_DOSE_CONFOUNDED', stratum,
        message: `${stratum}: Sex is confounded with Dose (a dose level carries only one ` +
          `sex), so the sex×dose interaction aliases the dose effect. Per-stratum fits stand; ` +
          `no sex-interaction claim.`});
      continue;
    }
    const ancova = stats.runAncova(lnAuc, lnDose, sexCode, {alpha: SEX_INTERACTION_ALPHA});
    // AD-7 (SY-F12): runAncova null (pooled N < k+2) → per-stratum fallback.
    if (ancova === null) {
      warnings.push({code: 'SEX_STRATUM_INSUFFICIENT', stratum,
        message: `${stratum}: too few pooled subjects for the sex-interaction ANCOVA; ` +
          `per-stratum fits stand; no sex-interaction claim.`});
      continue;
    }
    // Singular interaction model → the homogeneity F-test could not run
    // (`runAncova` returns pValue null and defaults `homogeneous` to true). Do
    // NOT report a pooled slope under unverified homogeneity — that is exactly
    // the false-confidence the 3-mode surface forbids. Fall back per-stratum.
    if (ancova.slopeHomogeneity.pValue === null) {
      warnings.push({code: 'SEX_STRATUM_INSUFFICIENT', stratum,
        message: `${stratum}: the sex×dose interaction test could not be evaluated ` +
          `(singular design); homogeneity is unverified. Per-stratum fits stand; no pooled slope.`});
      continue;
    }
    const homogeneous = ancova.slopeHomogeneity.homogeneous;
    // When homogeneous, the pooled common dose slope (runAncova's covariate slope)
    // carries df = N − 3 (intercept + sex offset + slope). Build its CI here.
    let pooledSlope: number | null = null;
    let pooledSlopeCI: [number, number] | null = null;
    let pooledDf: number | null = null;
    if (homogeneous) {
      const dfPooled = lnAuc.length - 3;
      if (dfPooled > 0) {
        pooledSlope = ancova.slope.estimate;
        pooledDf = dfPooled;
        const t = stats.distributions.studentTInv(1 - (1 - ciLevel) / 2, dfPooled);
        const half = t * ancova.slope.se;
        pooledSlopeCI = [pooledSlope - half, pooledSlope + half];
      }
    }
    out.push({
      route, occasion,
      fStatistic: ancova.slopeHomogeneity.fStatistic,
      pValue: ancova.slopeHomogeneity.pValue,
      homogeneous,
      pooledSlope, pooledSlopeCI, pooledDf,
    });
  }
  return out;
}

// =====================================================================
// Public entry point
// =====================================================================

/**
 * Compute the dose-proportionality (Smith power model) tables from a per-subject
 * NCA results frame carrying AUC/Cmax + the dose + route + sex metadata
 * (typically the post-`appendGroupingColumns` frame).
 *
 * @throws if a required column (dose / sex / the active AUC basis / Cmax) is
 *   absent from `results`.
 */
export function computeDoseProportionality(
  results: DG.DataFrame,
  options: DoseProportionalityOptions,
): DoseProportionalityResult {
  const basis = options.basis ?? 'aucinf';
  const aucCol = basis === 'aucinf' ? (options.aucInfCol ?? 'AUC_INF') : (options.aucLastCol ?? 'AUC_LAST');
  const cmaxName = options.cmaxCol ?? 'CMAX';
  const doseName = options.doseCol;
  const routeName = options.routeCol ?? 'ROUTE';
  const sexName = options.sexCol ?? 'SEX';
  const opts = {
    thetaLow: options.thetaLow ?? DEFAULT_THETA_LOW,
    thetaHigh: options.thetaHigh ?? DEFAULT_THETA_HIGH,
    ciLevel: options.ciLevel ?? DEFAULT_CI_LEVEL,
    minNSearchCap: options.minNSearchCap ?? DEFAULT_MIN_N_SEARCH_CAP,
  };

  for (const name of [aucCol, cmaxName, doseName, sexName]) {
    if (results.col(name) === null)
      throw new Error(`computeDoseProportionality: required column '${name}' not found in results.`);
  }

  const dose = results.col(doseName)!;
  const sex = results.col(sexName)!;
  const route = results.col(routeName); // optional
  const doseNorm = results.col('DOSE_NORM');
  const aucColumn = results.col(aucCol)!;
  const cmaxColumn = results.col(cmaxName)!;
  const n = results.rowCount;

  // Per-occasion mode (F2): ≥ 2 occasion levels on a frame that carries the occasion column.
  const occCol = options.occasionCol !== undefined ? results.col(options.occasionCol) : null;
  const levels = occCol !== null ? (options.occasionLevels ?? distinctStrings(occCol)) : [];
  const perOccasion = occCol !== null && levels.length >= 2;
  const animals = options.animals ?? null;
  const excluded = options.excludedAnimals ?? new Set<string>();
  const tlastCol = results.col(options.tlastCol ?? 'TLAST');
  const occasionOf = (i: number): string => perOccasion ? strAt(occCol, i) : '';

  // The AUC bases fitted per stratum: one (the run's basis) outside per-occasion mode; the
  // occasion rule's set inside it (F4).
  const aucInf = results.col(options.aucInfCol ?? 'AUC_INF');
  const aucLast = results.col(options.aucLastCol ?? 'AUC_LAST');
  const aucTau = options.aucTauCol !== undefined ? results.col(options.aucTauCol) : null;
  const colOfBasis = (b: AucBasis): DG.Column | null => b === 'aucinf' ? aucInf : b === 'auclast' ? aucLast : aucTau;
  const singleDoseOccasions = options.singleDoseOccasions ?? false;
  const tauAvailable = !singleDoseOccasions && aucTau !== null && anyFinite(aucTau);
  const basisPlan = perOccasion ?
    occasionBasisPlan(results, levels, occasionOf, dose, sex, route, {
      aucInf, pctExtrap: results.col(options.pctExtrapCol ?? 'PCT_EXTRAP'), tlast: tlastCol,
      tauH: options.tauH ?? NaN, aucTau, singleDoseOccasions,
    }) : null;

  // Build the per (route, sex[, occasion], parameter[, basis]) designs.
  const designs = new Map<string, StratumDesign>();
  const designKey = (rt: string, sx: string, p: string): string => `${rt}${sx}${p}`;
  for (let i = 0; i < n; i++) {
    const d = numAt(dose, i);
    if (!Number.isFinite(d)) continue;
    const rt = route !== null ? strAt(route, i) : '';
    const sx = strAt(sex, i);
    const oc = occasionOf(i);
    const dn = strAt(doseNorm, i);
    const cmaxSpec = {spec: {parameter: 'Cmax' as const, col: cmaxColumn, basis: null}, carries: true, notes: []};
    const plan: {spec: ParamSpec; carries: boolean; notes: string[]}[] = basisPlan === null ?
      [{spec: {parameter: 'AUC', col: aucColumn, basis}, carries: true, notes: []}, cmaxSpec] :
      [...(basisPlan.get(`${rt}|${sx}|${oc}`) ?? []).map((p) => ({
        spec: {parameter: 'AUC' as const, col: colOfBasis(p.basis), basis: p.basis}, carries: p.carries,
        notes: p.notes,
      })), cmaxSpec];
    const an = animals !== null ? animals[i] : '';
    for (const {spec: ps, carries, notes} of plan) {
      const key = designKey(rt, sx, `${oc}|${ps.parameter}|${ps.basis ?? ''}`);
      let des = designs.get(key);
      if (des === undefined) {
        des = {route: rt, sex: sx, occasion: oc, parameter: ps.parameter, basis: ps.basis,
          carriesVerdict: carries, notes: [...notes], subjects: [], doseNorms: new Set()};
        designs.set(key, des);
      }
      const v = ps.col !== null ? numAt(ps.col, i) : NaN;
      des.subjects.push({dose: d, value: v, usable: Number.isFinite(v), animal: an,
        tlast: tlastCol !== null ? numAt(tlastCol, i) : NaN});
      if (dn !== '') des.doseNorms.add(dn);
    }
  }

  const warnings: DoseProportionalityWarning[] = [];
  const fits: DoseProportionalityFit[] = [];
  const occRank = (o: string): number => levels.indexOf(o);
  const basisRank = (b: AucBasis | null): number => b === null ? 9 : ['aucinf', 'auctau', 'auclast'].indexOf(b);
  const designList = [...designs.values()].sort((a, b) =>
    a.route.localeCompare(b.route) || a.sex.localeCompare(b.sex) ||
    occRank(a.occasion) - occRank(b.occasion) ||
    a.parameter.localeCompare(b.parameter) || basisRank(a.basis) - basisRank(b.basis));
  for (const des of designList) fits.push(fitStratum(des, opts, warnings));
  if (perOccasion) markVerdictDiffers(fits);

  // Sex slope-homogeneity, per route (× occasion) (AD-7).
  const byRoute = new Map<string, StratumDesign[]>();
  for (const des of designList) {
    const k = `${des.route}|${des.occasion}`;
    let arr = byRoute.get(k);
    if (arr === undefined) {arr = []; byRoute.set(k, arr);}
    arr.push(des);
  }
  const sexHomogeneity = computeSexHomogeneity(byRoute, opts.ciLevel, warnings);
  const dnAnova = computeDnAucAnova(designs);

  const table = buildVerdictTable(fits, basis, perOccasion);
  const dnAucName = perOccasion && aucInf !== null ? aucInf.name : aucCol;
  const dnTable = buildDnTable(results, {doseName, routeName, sexName, cmaxName, aucCol: dnAucName,
    route: route !== null, subjectName: options.subjectCol ?? 'SUBJECT',
    occasionName: perOccasion ? occCol!.name : null, animals});

  const deltaBeta = perOccasion ?
    computeDeltaBeta(designList, levels, animals !== null, excluded, opts.ciLevel, tauAvailable,
      singleDoseOccasions) : [];
  const deltaTable = perOccasion ? buildDeltaTable(deltaBeta) : null;

  return {table, fits, dnTable, sexHomogeneity, dnAnova, basis, warnings, deltaBeta, deltaTable};
}

function distinctStrings(col: DG.Column): string[] {
  const out: string[] = [];
  for (let i = 0; i < col.length; i++) {
    const v = strAt(col, i);
    if (v !== '' && !out.includes(v)) out.push(v);
  }
  return out;
}

function anyFinite(col: DG.Column): boolean {
  for (let i = 0; i < col.length; i++) if (Number.isFinite(numAt(col, i))) return true;
  return false;
}

/**
 * The occasion rule (F4), per (route, sex, occasion) stratum: which AUC fits are shown and
 * which one carries the verdict.
 * - First occasion: AUC∞ carries the verdict only when it is valid for the WHOLE design (every
 *   dosed profile has an AUC∞ with %AUCextrap ≤ 20 % and, τ known, its λz fitted on points
 *   ≤ τ — the auto window ends at TLAST, so TLAST ≤ τ; after a second dose the terminal phase
 *   is no longer single-dose) — never per-animal removal; otherwise AUC0–τ does. AUC0–τ is
 *   always shown beside.
 * - A later occasion: AUC0–τ carries the verdict; AUC∞ is shown beside.
 * - AUC0–τ unavailable for the stratum (τ unknown, or τ not spanned by EVERY dosed profile of
 *   the stratum): AUClast carries it, labelled as the fallback — decided for the whole stratum,
 *   like the AUC∞ guard, so the fit never runs on the subset of animals sampled long enough.
 * - Single-dose occasions (a confirmed `crossover`: every period is one dose): every occasion
 *   takes the first-occasion rule without AUC0–τ (a dosing-interval quantity) — AUC∞ when valid
 *   for the whole stratum, else AUClast.
 */
function occasionBasisPlan(
  results: DG.DataFrame, levels: ReadonlyArray<string>, occasionOf: (i: number) => string,
  dose: DG.Column, sex: DG.Column, route: DG.Column | null,
  c: {aucInf: DG.Column | null; pctExtrap: DG.Column | null; tlast: DG.Column | null; tauH: number;
    aucTau: DG.Column | null; singleDoseOccasions: boolean},
): Map<string, {basis: AucBasis; carries: boolean; notes: string[]}[]> {
  const singleDose = c.singleDoseOccasions;
  const tauH = singleDose ? NaN : c.tauH;
  const valid = new Map<string, boolean>();
  const tauDosed = new Map<string, {n: number; spanned: number}>();
  for (let i = 0; i < results.rowCount; i++) {
    const d = numAt(dose, i);
    if (!Number.isFinite(d)) continue;
    const key = `${route !== null ? strAt(route, i) : ''}|${strAt(sex, i)}|${occasionOf(i)}`;
    if (!valid.has(key)) {valid.set(key, true); tauDosed.set(key, {n: 0, spanned: 0});}
    if (d <= 0) continue;
    const inf = c.aucInf !== null ? numAt(c.aucInf, i) : NaN;
    const ext = c.pctExtrap !== null ? numAt(c.pctExtrap, i) : NaN;
    const tl = c.tlast !== null ? numAt(c.tlast, i) : NaN;
    const lzWithinTau = !Number.isFinite(tauH) || (Number.isFinite(tl) && tl <= tauH * (1 + 1e-6));
    const ok = Number.isFinite(inf) && inf > 0 && Number.isFinite(ext) && ext <= AUCINF_VALID_EXTRAP_PCT &&
      lzWithinTau;
    valid.set(key, valid.get(key)! && ok);
    const td = tauDosed.get(key)!;
    td.n++;
    const tv = c.aucTau !== null ? numAt(c.aucTau, i) : NaN;
    if (Number.isFinite(tv) && tv > 0) td.spanned++;
  }
  const out = new Map<string, {basis: AucBasis; carries: boolean; notes: string[]}[]>();
  for (const [key, allValid] of valid) {
    const occ = key.substring(key.lastIndexOf('|') + 1);
    const firstRule = singleDose || occ === levels[0];
    const infCarries = firstRule && allValid;
    const infNote = firstRule && !allValid ?
      [`AUC∞ not valid for every animal (%AUCextrap > ${AUCINF_VALID_EXTRAP_PCT} % or not estimated)`] : [];
    const td = tauDosed.get(key)!;
    const tauSpanned = !singleDose && c.aucTau !== null && td.n > 0 && td.spanned === td.n;
    let fallback: string[] = [];
    if (!infCarries) {
      fallback = singleDose ? ['AUClast: single-dose periods, AUC∞ not valid'] :
        c.aucTau === null ? ['AUClast: τ unknown — AUC0–τ unavailable'] :
          [`AUClast: τ not spanned by ${td.n - td.spanned} of ${td.n} profiles — AUC0–τ unavailable`];
    }
    out.set(key, tauSpanned ?
      [{basis: 'aucinf', carries: infCarries, notes: infNote}, {basis: 'auctau', carries: !infCarries, notes: []}] :
      [{basis: 'aucinf', carries: infCarries, notes: infNote},
        {basis: 'auclast', carries: !infCarries, notes: fallback}]);
  }
  return out;
}

/** The verdict class, for "verdict differs on the other basis". */
function verdictClass(v: DoseProportionalityVerdict): string {
  return v === 'proportional' || v === 'sub-proportional' || v === 'supra-proportional' ? v : 'inconclusive';
}

/** Mark the verdict-carrying AUC row of a stratum whose shown-beside AUC row has another class. */
function markVerdictDiffers(fits: DoseProportionalityFit[]): void {
  const key = (f: DoseProportionalityFit): string => `${f.route}|${f.sex}|${f.occasion}`;
  for (let i = 0; i < fits.length; i++) {
    const f = fits[i];
    if (f.parameter !== 'AUC' || !f.carriesVerdict || f.mode === 'inapplicable') continue;
    const other = fits.find((g) => g !== f && g.parameter === 'AUC' && key(g) === key(f) &&
      g.mode !== 'inapplicable' && Number.isFinite(g.beta));
    if (other !== undefined && verdictClass(other.verdict) !== verdictClass(f.verdict))
      fits[i] = {...f, notes: [...f.notes, VERDICT_DIFFERS_NOTE]};
  }
}

/**
 * Δβ = β(later) − β(first) per route × sex, on ONE common basis (F3a, SF-DELTA-BETA-BASIS):
 * AUC0–τ. Without AUC0–τ, AUClast is the common basis only when the two occasions' sampling
 * spans are equal (median tlast within 5 %); otherwise the row is suppressed — a slope
 * difference across two AUC definitions is not a time effect.
 * - Animal key: the exact paired-difference regression d = ln AUC(later) − ln AUC(first) on
 *   ln dose (`stats.linearFit`), complete-case; df = paired − 2.
 * - No animal key: two independent fits, SE = √(SE₁² + SE₂²), Welch–Satterthwaite df.
 * No colour, flag or p-value (the F3c cue waits for the null-effect fixtures).
 */
function computeDeltaBeta(
  designs: ReadonlyArray<StratumDesign>, levels: ReadonlyArray<string>, haveAnimals: boolean,
  excluded: ReadonlySet<string>, ciLevel: number, tauAvailable: boolean, singleDoseOccasions: boolean,
): DeltaBetaRow[] {
  const out: DeltaBetaRow[] = [];
  const strata = new Map<string, StratumDesign[]>();
  for (const d of designs) {
    if (d.parameter !== 'AUC') continue;
    const k = `${d.route}|${d.sex}`;
    strata.set(k, [...(strata.get(k) ?? []), d]);
  }
  const usable = (s: SubjectPoint): boolean => s.usable && s.dose > 0 && Number.isFinite(s.value) && s.value > 0;
  const tlastMedian = (d: StratumDesign): number => {
    const t = d.subjects.map((s) => s.tlast).filter(Number.isFinite).sort((a, b) => a - b);
    return t.length > 0 ? t[Math.floor(t.length / 2)] : NaN;
  };
  const numericOrder = levels.every((l) => Number.isFinite(occasionNumber(l)));
  for (const ds of strata.values()) {
    const route = ds[0].route;
    const sex = ds[0].sex;
    const first = levels[0];
    const pick = (occ: string, b: AucBasis): StratumDesign | undefined =>
      ds.find((d) => d.occasion === occ && d.basis === b);
    for (const later of levels.slice(1)) {
      const notes: string[] = [];
      // ≥ 3 occasions or a first-appearance order: the row label carries the pair (R1 F3).
      if (levels.length > 2 || !numericOrder) notes.push(`${later} vs ${first}`);
      // A crossover's periods are single doses: the slope difference is a period effect.
      if (singleDoseOccasions) notes.push('period effect (single-dose periods), not accumulation');
      let basis: AucBasis | null = null;
      if (tauAvailable && pick(first, 'auctau') !== undefined && pick(later, 'auctau') !== undefined)
        basis = 'auctau';
      else {
        const a = pick(first, 'auclast');
        const b = pick(later, 'auclast');
        const t1 = a !== undefined ? tlastMedian(a) : NaN;
        const t2 = b !== undefined ? tlastMedian(b) : NaN;
        if (Number.isFinite(t1) && Number.isFinite(t2) && Math.abs(t1 - t2) <= 0.05 * Math.max(t1, t2)) {
          basis = 'auclast';
          notes.push('AUClast on both occasions (equal sampling spans)');
        }
      }
      const base = {route, sex, first, later};
      if (basis === null) {
        out.push({...base, basis: null, method: 'suppressed', deltaBeta: NaN, ci: [NaN, NaN], df: NaN,
          nFirst: 0, nLater: 0, nPaired: 0, nByDoseFirst: '', nByDoseLater: '',
          notes: ['not comparable: different sampling spans', ...notes]});
        continue;
      }
      const u1 = pick(first, basis)!.subjects.filter(usable);
      const u2 = pick(later, basis)!.subjects.filter(usable);
      if (haveAnimals) {
        const byAnimal = new Map<string, SubjectPoint>();
        for (const s of u1) if (!excluded.has(s.animal)) byAnimal.set(s.animal, s);
        const x: number[] = [];
        const y: number[] = [];
        const pairedDoses: number[] = [];
        for (const s of u2) {
          const p = excluded.has(s.animal) ? undefined : byAnimal.get(s.animal);
          if (p === undefined || p.dose !== s.dose) continue;
          x.push(Math.log(s.dose));
          y.push(Math.log(s.value) - Math.log(p.value));
          pairedDoses.push(s.dose);
        }
        notes.push(...attritionNote(u1, u2, pairedDoses));
        const n = {nFirst: u1.length, nLater: u2.length, nPaired: x.length,
          nByDoseFirst: nByDoseOf(u1), nByDoseLater: nByDoseOf(u2)};
        if (x.length < 3 || new Set(x).size < 2) {
          out.push({...base, basis, method: 'paired', deltaBeta: NaN, ci: [NaN, NaN], df: NaN, ...n,
            notes: ['too few paired animals', ...notes]});
          continue;
        }
        const fit = stats.linearFit(x, y, {ciLevel});
        out.push({...base, basis, method: 'paired', deltaBeta: fit.slope,
          ci: [fit.slopeCI[0], fit.slopeCI[1]], df: fit.df, ...n, notes});
        continue;
      }
      const f1 = fitLog(u1, ciLevel);
      const f2 = fitLog(u2, ciLevel);
      const n = {nFirst: u1.length, nLater: u2.length, nPaired: 0,
        nByDoseFirst: nByDoseOf(u1), nByDoseLater: nByDoseOf(u2)};
      if (f1 === null || f2 === null) {
        out.push({...base, basis, method: 'two-fits', deltaBeta: NaN, ci: [NaN, NaN], df: NaN, ...n,
          notes: ['too few profiles', ...notes]});
        continue;
      }
      const v1 = f1.se * f1.se;
      const v2 = f2.se * f2.se;
      const se = Math.sqrt(v1 + v2);
      const df = (v1 + v2) ** 2 / (v1 * v1 / f1.df + v2 * v2 / f2.df);
      const t = stats.distributions.studentTInv(1 - (1 - ciLevel) / 2, df);
      const db = f2.slope - f1.slope;
      out.push({...base, basis, method: 'two-fits', deltaBeta: db, ci: [db - t * se, db + t * se], df, ...n,
        notes: ['conservative (two independent fits)', ...notes]});
    }
  }
  return out;
}

/** OLS ln(value) ~ ln(dose), with sci-comp's own slope SE. */
function fitLog(points: ReadonlyArray<SubjectPoint>, ciLevel: number): {slope: number; se: number; df: number} | null {
  if (points.length < 3 || new Set(points.map((p) => p.dose)).size < 2) return null;
  const fit = stats.linearFit(points.map((p) => Math.log(p.dose)), points.map((p) => Math.log(p.value)), {ciLevel});
  return {slope: fit.slope, se: fit.slopeSe, df: fit.df};
}

/** Attrition disclosure (R2 F7): `paired N < occasion N`, plus `dropped animals concentrated
 *  at <dose>` when every dropped animal sits at one dose level. */
function attritionNote(u1: ReadonlyArray<SubjectPoint>, u2: ReadonlyArray<SubjectPoint>,
  pairedDoses: ReadonlyArray<number>): string[] {
  const nPaired = pairedDoses.length;
  if (nPaired >= u1.length && nPaired >= u2.length) return [];
  const count = (doses: ReadonlyArray<number>): Map<number, number> => {
    const m = new Map<number, number>();
    for (const d of doses) m.set(d, (m.get(d) ?? 0) + 1);
    return m;
  };
  const cp = count(pairedDoses);
  const droppedAt = new Set<number>();
  for (const u of [u1, u2])
    for (const [d, k] of count(u.map((s) => s.dose))) if (k > (cp.get(d) ?? 0)) droppedAt.add(d);
  let note = `paired N ${nPaired} < occasion N ${u1.length} / ${u2.length}`;
  if (droppedAt.size === 1) note += `; dropped animals concentrated at ${[...droppedAt][0]}`;
  return [note];
}

/**
 * Secondary dn-AUC one-way ANOVA across dose groups, per Route×Sex (AD-9). Uses
 * the AUC design's per-subject dn-AUC; `null` F/p when the ANOVA is undefined
 * (< 2 dose groups, no within-group variation — `oneWayAnova` returns null).
 */
function computeDnAucAnova(designs: Map<string, StratumDesign>): DnAucAnova[] {
  const out: DnAucAnova[] = [];
  for (const des of designs.values()) {
    if (des.parameter !== 'AUC' || !des.carriesVerdict) continue;
    const dn: number[] = [];
    const grp: number[] = [];
    const levelCode = new Map<number, number>();
    for (const s of des.subjects) {
      if (!s.usable || !(s.dose > 0) || !Number.isFinite(s.value) || !(s.value > 0)) continue;
      let code = levelCode.get(s.dose);
      if (code === undefined) {code = levelCode.size; levelCode.set(s.dose, code);}
      dn.push(s.value / s.dose);
      grp.push(code);
    }
    const anova = dn.length > 0 ? stats.oneWayAnova(dn, grp) : null;
    out.push({
      route: des.route, sex: des.sex, occasion: des.occasion,
      fStatistic: anova?.fStatistic ?? null,
      pValue: anova?.pValue ?? null,
      nGroups: levelCode.size,
      n: dn.length,
    });
  }
  return out.sort((a, b) => a.route.localeCompare(b.route) || a.sex.localeCompare(b.sex));
}

// --- Output table builders ----------------------------------------------

const VERDICT_LABEL: Record<DoseProportionalityVerdict, string> = {
  'proportional': 'Proportional',
  'sub-proportional': 'Sub-proportional',
  'supra-proportional': 'Supra-proportional',
  'inconclusive-too-few-levels': 'Inconclusive (too few dose levels)',
  'inconclusive-composite': 'Inconclusive (composite — no per-subject AUC)',
  'inconclusive-underpowered': 'Inconclusive (underpowered by IIV)',
  'insufficient-completeness': 'Insufficient completeness',
  'inapplicable': 'Inapplicable (< 3 dose levels)',
};

function buildVerdictTable(
  fits: ReadonlyArray<DoseProportionalityFit>, basis: AucBasis, perOccasion = false,
): DG.DataFrame {
  if (perOccasion) return buildOccasionVerdictTable(fits);
  const m = fits.length;
  const routeOut = new Array<string>(m);
  const sexOut = new Array<string>(m);
  const paramOut = new Array<string>(m);
  const basisOut = new Array<string>(m);
  const betaOut = new Float64Array(m);
  const ciLoOut = new Float64Array(m);
  const ciHiOut = new Float64Array(m);
  const rOut = new Float64Array(m);
  const crLoOut = new Float64Array(m);
  const crHiOut = new Float64Array(m);
  const verdictOut = new Array<string>(m);
  const modeOut = new Array<string>(m);
  const nOut = new Float64Array(m);
  const dfOut = new Float64Array(m);
  const levelsOut = new Float64Array(m);
  const r2Out = new Float64Array(m);
  const minNOut = new Float64Array(m);

  for (let i = 0; i < m; i++) {
    const f = fits[i];
    routeOut[i] = f.route;
    sexOut[i] = f.sex;
    paramOut[i] = f.parameter;
    basisOut[i] = f.parameter === 'AUC' ? (f.basis ?? basis) : '';
    betaOut[i] = f.beta;
    ciLoOut[i] = f.slopeCI[0];
    ciHiOut[i] = f.slopeCI[1];
    rOut[i] = f.doseRatio;
    crLoOut[i] = f.criticalRegion[0];
    crHiOut[i] = f.criticalRegion[1];
    verdictOut[i] = VERDICT_LABEL[f.verdict] + (f.isGroupMeanFallback ? ' [group-mean df-limited]' : '');
    modeOut[i] = f.mode;
    nOut[i] = f.n;
    dfOut[i] = Number.isFinite(f.df) ? f.df : NaN;
    levelsOut[i] = f.nDoseLevels;
    r2Out[i] = f.rSquared;
    minNOut[i] = f.minNAdvisory?.nStar ?? NaN;
  }

  const cols: DG.Column[] = [
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'ROUTE', routeOut),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'SEX', sexOut),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'PARAMETER', paramOut),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'BASIS', basisOut),
    f64Column('BETA', betaOut).setTag('label', 'β (power slope)'),
    f64Column('CI_LOW', ciLoOut),
    f64Column('CI_HIGH', ciHiOut),
    f64Column('DOSE_RATIO', rOut).setTag('label', 'r (max/min)'),
    f64Column('CR_LOW', crLoOut),
    f64Column('CR_HIGH', crHiOut),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'VERDICT', verdictOut),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'MODE', modeOut),
    f64Column('N', nOut),
    f64Column('DF', dfOut),
    f64Column('N_LEVELS', levelsOut),
    f64Column('R_SQUARED', r2Out),
    f64Column('MIN_N_PER_GROUP', minNOut).setTag('label', 'min n/group for a verdict'),
  ];
  const table = DGNS.DataFrame.fromColumns(cols);
  table.name = 'Dose Proportionality';
  return table;
}

/**
 * The per-occasion verdict table: one row per route × sex × occasion × parameter identity
 * (`AUC∞`, `AUC0–τ`, `AUClast`, `Cmax` — the basis IS the parameter, architect S7). `VERDICT`
 * is filled only on the row that carries the stratum's verdict; a row shown beside it carries
 * β + CI only. `N_BY_DOSE` and `NOTE` disclose attrition, fallbacks and a basis disagreement.
 */
function buildOccasionVerdictTable(fits: ReadonlyArray<DoseProportionalityFit>): DG.DataFrame {
  const s = (f: (x: DoseProportionalityFit) => string, name: string): DG.Column =>
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, name, fits.map(f));
  const v = (f: (x: DoseProportionalityFit) => number, name: string): DG.Column =>
    f64Column(name, Float64Array.from(fits, f));
  const table = DGNS.DataFrame.fromColumns([
    s((f) => f.route, 'ROUTE'),
    s((f) => f.sex, 'SEX'),
    s((f) => f.occasion, 'OCCASION'),
    s((f) => f.parameter === 'Cmax' ? DP_PARAMETER_ID.cmax : DP_PARAMETER_ID[f.basis ?? 'aucinf'], 'PARAMETER'),
    v((f) => f.beta, 'BETA').setTag('label', 'β (power slope)'),
    v((f) => f.slopeCI[0], 'CI_LOW'),
    v((f) => f.slopeCI[1], 'CI_HIGH'),
    v((f) => f.doseRatio, 'DOSE_RATIO').setTag('label', 'r (max/min)'),
    v((f) => f.criticalRegion[0], 'CR_LOW'),
    v((f) => f.criticalRegion[1], 'CR_HIGH'),
    s((f) => f.carriesVerdict ? VERDICT_LABEL[f.verdict] + (f.isGroupMeanFallback ? ' [group-mean df-limited]' : '') :
      '', 'VERDICT'),
    s((f) => f.carriesVerdict ? f.mode : '', 'MODE'),
    v((f) => f.n, 'N'),
    s((f) => f.nByDose, 'N_BY_DOSE').setTag('label', 'N per dose'),
    v((f) => Number.isFinite(f.df) ? f.df : NaN, 'DF'),
    v((f) => f.nDoseLevels, 'N_LEVELS'),
    v((f) => f.rSquared, 'R_SQUARED'),
    v((f) => f.minNAdvisory?.nStar ?? NaN, 'MIN_N_PER_GROUP').setTag('label', 'min n/group for a verdict'),
    s((f) => f.notes.join('; '), 'NOTE'),
  ]);
  table.name = 'Dose Proportionality';
  return table;
}

/** The Δβ table (F3a): one row per route × sex × later occasion. */
function buildDeltaTable(rows: ReadonlyArray<DeltaBetaRow>): DG.DataFrame {
  const s = (f: (x: DeltaBetaRow) => string, name: string): DG.Column =>
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, name, rows.map(f));
  const v = (f: (x: DeltaBetaRow) => number, name: string): DG.Column => f64Column(name, Float64Array.from(rows, f));
  const table = DGNS.DataFrame.fromColumns([
    s((r) => r.route, 'ROUTE'),
    s((r) => r.sex, 'SEX'),
    s((r) => r.first, 'OCCASION_FIRST'),
    s((r) => r.later, 'OCCASION_LATER'),
    s((r) => r.basis === null ? '' : DP_PARAMETER_ID[r.basis], 'PARAMETER'),
    s((r) => r.method, 'METHOD'),
    v((r) => r.deltaBeta, 'DELTA_BETA').setTag('label', 'Δβ (later − first)'),
    v((r) => r.ci[0], 'CI_LOW'),
    v((r) => r.ci[1], 'CI_HIGH'),
    v((r) => r.df, 'DF'),
    v((r) => r.nFirst, 'N_FIRST'),
    v((r) => r.nLater, 'N_LATER'),
    v((r) => r.nPaired, 'N_PAIRED'),
    s((r) => r.nByDoseFirst, 'N_BY_DOSE_FIRST').setTag('label', 'N per dose (first)'),
    s((r) => r.nByDoseLater, 'N_BY_DOSE_LATER').setTag('label', 'N per dose (later)'),
    s((r) => r.notes.join('; '), 'NOTE'),
  ]);
  table.name = 'Dose Proportionality over time';
  return table;
}

function buildDnTable(
  results: DG.DataFrame,
  cols: {doseName: string; routeName: string; sexName: string; cmaxName: string; aucCol: string; route: boolean;
    subjectName: string; occasionName: string | null; animals: ReadonlyArray<string> | null},
): DG.DataFrame {
  const dose = results.col(cols.doseName)!;
  const sex = results.col(cols.sexName)!;
  const route = cols.route ? results.col(cols.routeName) : null;
  const auc = results.col(cols.aucCol)!;
  const cmax = results.col(cols.cmaxName)!;
  const subject = results.col(cols.subjectName);
  const occasion = cols.occasionName !== null ? results.col(cols.occasionName) : null;
  const n = results.rowCount;

  const routeOut: string[] = [];
  const sexOut: string[] = [];
  const subjectOut: string[] = [];
  const animalOut: string[] = [];
  const occasionOut: string[] = [];
  const doseOut: number[] = [];
  const dnAucOut: number[] = [];
  const dnCmaxOut: number[] = [];
  for (let i = 0; i < n; i++) {
    const d = numAt(dose, i);
    if (!Number.isFinite(d) || d <= 0) continue;
    routeOut.push(strAt(route, i));
    sexOut.push(strAt(sex, i));
    subjectOut.push(strAt(subject, i));
    if (cols.animals !== null) animalOut.push(cols.animals[i]);
    if (occasion !== null) occasionOut.push(strAt(occasion, i));
    doseOut.push(d);
    const a = numAt(auc, i);
    const c = numAt(cmax, i);
    dnAucOut.push(Number.isFinite(a) ? a / d : NaN);
    dnCmaxOut.push(Number.isFinite(c) ? c / d : NaN);
  }

  // G9: the subject (and, when known, the animal and the occasion) make each row identifiable.
  const out: DG.Column[] = [
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'ROUTE', routeOut),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'SEX', sexOut),
    ...(occasion !== null ? [DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'OCCASION', occasionOut)] : []),
    ...(subject !== null ? [DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'SUBJECT', subjectOut)] : []),
    ...(cols.animals !== null && occasion !== null ?
      [DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'ANIMAL', animalOut)] : []),
    f64Column('DOSE', doseOut),
    f64Column('DN_AUC', dnAucOut).setTag('label', 'AUC / Dose'),
    f64Column('DN_CMAX', dnCmaxOut).setTag('label', 'Cmax / Dose'),
  ];
  const table = DGNS.DataFrame.fromColumns(out);
  table.name = 'Dose-Normalized Exposure';
  return table;
}
