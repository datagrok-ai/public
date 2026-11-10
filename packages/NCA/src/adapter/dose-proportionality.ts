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
export type AucBasis = 'aucinf' | 'auclast';

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
  readonly parameter: 'AUC' | 'Cmax';
  readonly basis: AucBasis | null; // null for Cmax
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
  readonly fStatistic: number | null;
  readonly pValue: number | null;
  readonly nGroups: number;
  readonly n: number;
}

/** Per-route sex slope-homogeneity result (AD-7). */
export interface SexHomogeneity {
  readonly route: string;
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
interface SubjectPoint {dose: number; value: number; usable: boolean;}

/** One Route×Sex × parameter design before fitting. */
interface StratumDesign {
  route: string;
  sex: string;
  parameter: 'AUC' | 'Cmax';
  basis: AucBasis | null;
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
  warnings: DoseProportionalityWarning[],
): DoseProportionalityFit {
  const stratumKey = `ROUTE=${design.route}|SEX=${design.sex}`;
  const paramKey = `${stratumKey}|${design.parameter}`;
  const thetaLow = options.thetaLow;
  const thetaHigh = options.thetaHigh;

  const degenerate = (
    verdict: DoseProportionalityVerdict, mode: VerdictMode,
    nLevels: number, n: number,
  ): DoseProportionalityFit => ({
    route: design.route, sex: design.sex, parameter: design.parameter, basis: design.basis,
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
    route: design.route, sex: design.sex, parameter: design.parameter, basis: design.basis,
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
  for (const [route, designs] of byRoute) {
    // Pool the AUC designs of this route across sexes.
    const aucDesigns = designs.filter((d) => d.parameter === 'AUC');
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
        warnings.push({code: 'SEX_STRATUM_INSUFFICIENT', stratum: `ROUTE=${route}`,
          message: `ROUTE=${route}: sex slope-homogeneity needs exactly two encodable sex ` +
            `categories (F/M); found ${sexes.size}. Per-stratum fits stand; no sex-interaction claim.`});
      }
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
      warnings.push({code: 'SEX_DOSE_CONFOUNDED', stratum: `ROUTE=${route}`,
        message: `ROUTE=${route}: Sex is confounded with Dose (a dose level carries only one ` +
          `sex), so the sex×dose interaction aliases the dose effect. Per-stratum fits stand; ` +
          `no sex-interaction claim.`});
      continue;
    }
    const ancova = stats.runAncova(lnAuc, lnDose, sexCode, {alpha: SEX_INTERACTION_ALPHA});
    // AD-7 (SY-F12): runAncova null (pooled N < k+2) → per-stratum fallback.
    if (ancova === null) {
      warnings.push({code: 'SEX_STRATUM_INSUFFICIENT', stratum: `ROUTE=${route}`,
        message: `ROUTE=${route}: too few pooled subjects for the sex-interaction ANCOVA; ` +
          `per-stratum fits stand; no sex-interaction claim.`});
      continue;
    }
    // Singular interaction model → the homogeneity F-test could not run
    // (`runAncova` returns pValue null and defaults `homogeneous` to true). Do
    // NOT report a pooled slope under unverified homogeneity — that is exactly
    // the false-confidence the 3-mode surface forbids. Fall back per-stratum.
    if (ancova.slopeHomogeneity.pValue === null) {
      warnings.push({code: 'SEX_STRATUM_INSUFFICIENT', stratum: `ROUTE=${route}`,
        message: `ROUTE=${route}: the sex×dose interaction test could not be evaluated ` +
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
      route,
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

  const params: ParamSpec[] = [
    {parameter: 'AUC', col: aucColumn, basis},
    {parameter: 'Cmax', col: cmaxColumn, basis: null},
  ];

  // Build the per (route, sex, parameter) designs.
  const designs = new Map<string, StratumDesign>();
  const designKey = (rt: string, sx: string, p: 'AUC' | 'Cmax'): string => `${rt}${sx}${p}`;
  for (let i = 0; i < n; i++) {
    const d = numAt(dose, i);
    if (!Number.isFinite(d)) continue;
    const rt = route !== null ? strAt(route, i) : '';
    const sx = strAt(sex, i);
    const dn = strAt(doseNorm, i);
    for (const ps of params) {
      const key = designKey(rt, sx, ps.parameter);
      let des = designs.get(key);
      if (des === undefined) {
        des = {route: rt, sex: sx, parameter: ps.parameter, basis: ps.basis, subjects: [], doseNorms: new Set()};
        designs.set(key, des);
      }
      const v = ps.col !== null ? numAt(ps.col, i) : NaN;
      des.subjects.push({dose: d, value: v, usable: Number.isFinite(v)});
      if (dn !== '') des.doseNorms.add(dn);
    }
  }

  const warnings: DoseProportionalityWarning[] = [];
  const fits: DoseProportionalityFit[] = [];
  const designList = [...designs.values()].sort((a, b) =>
    a.route.localeCompare(b.route) || a.sex.localeCompare(b.sex) ||
    a.parameter.localeCompare(b.parameter));
  for (const des of designList) fits.push(fitStratum(des, opts, warnings));

  // Sex slope-homogeneity, per route (AD-7).
  const byRoute = new Map<string, StratumDesign[]>();
  for (const des of designList) {
    let arr = byRoute.get(des.route);
    if (arr === undefined) {arr = []; byRoute.set(des.route, arr);}
    arr.push(des);
  }
  const sexHomogeneity = computeSexHomogeneity(byRoute, opts.ciLevel, warnings);
  const dnAnova = computeDnAucAnova(designs);

  const table = buildVerdictTable(fits, basis);
  const dnTable = buildDnTable(results, {doseName, routeName, sexName, cmaxName, aucCol, route: route !== null});

  return {table, fits, dnTable, sexHomogeneity, dnAnova, basis, warnings};
}

/**
 * Secondary dn-AUC one-way ANOVA across dose groups, per Route×Sex (AD-9). Uses
 * the AUC design's per-subject dn-AUC; `null` F/p when the ANOVA is undefined
 * (< 2 dose groups, no within-group variation — `oneWayAnova` returns null).
 */
function computeDnAucAnova(designs: Map<string, StratumDesign>): DnAucAnova[] {
  const out: DnAucAnova[] = [];
  for (const des of designs.values()) {
    if (des.parameter !== 'AUC') continue;
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
      route: des.route, sex: des.sex,
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

function buildVerdictTable(fits: ReadonlyArray<DoseProportionalityFit>, basis: AucBasis): DG.DataFrame {
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

function buildDnTable(
  results: DG.DataFrame,
  cols: {doseName: string; routeName: string; sexName: string; cmaxName: string; aucCol: string; route: boolean},
): DG.DataFrame {
  const dose = results.col(cols.doseName)!;
  const sex = results.col(cols.sexName)!;
  const route = cols.route ? results.col(cols.routeName) : null;
  const auc = results.col(cols.aucCol)!;
  const cmax = results.col(cols.cmaxName)!;
  const n = results.rowCount;

  const routeOut: string[] = [];
  const sexOut: string[] = [];
  const doseOut: number[] = [];
  const dnAucOut: number[] = [];
  const dnCmaxOut: number[] = [];
  for (let i = 0; i < n; i++) {
    const d = numAt(dose, i);
    if (!Number.isFinite(d) || d <= 0) continue;
    routeOut.push(strAt(route, i));
    sexOut.push(strAt(sex, i));
    doseOut.push(d);
    const a = numAt(auc, i);
    const c = numAt(cmax, i);
    dnAucOut.push(Number.isFinite(a) ? a / d : NaN);
    dnCmaxOut.push(Number.isFinite(c) ? c / d : NaN);
  }

  const out: DG.Column[] = [
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'ROUTE', routeOut),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'SEX', sexOut),
    f64Column('DOSE', doseOut),
    f64Column('DN_AUC', dnAucOut).setTag('label', 'AUC / Dose'),
    f64Column('DN_CMAX', dnCmaxOut).setTag('label', 'Cmax / Dose'),
  ];
  const table = DGNS.DataFrame.fromColumns(out);
  table.name = 'Dose-Normalized Exposure';
  return table;
}
