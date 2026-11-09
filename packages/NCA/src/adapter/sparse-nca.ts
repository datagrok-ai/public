/**
 * Sparse / destructive-sampling NCA orchestration (UC-04 / FR-301..306, P0).
 *
 * **Why a NEW module, not `aggregateNcaResults` (mirror of `dose-proportionality.ts`
 * / `bioavailability.ts` / `toxicokinetics.ts`).** The native aggregator does
 * *within-group* arithmetic on a *per-subject results* frame. Sparse NCA has no
 * per-subject result — each animal contributes one (destructive) or a few (batch)
 * raw samples, so the analysis is a *cross-animal composite* built from the raw
 * observations: pool by nominal time → composite mean profile → closed-form AUC
 * variance. `aggregateNcaResults` cannot build a composite or call the Holder
 * variance, so this is a separate view-invoked primitive. The dense pipeline
 * (controller, aggregator) is untouched (AD-13).
 *
 * **Correct layering (not vendored math — Rule 5).** All sparse statistics come
 * from sci-comp v0.9.0: `nca.sparseAuc` (Holder A1+A3 variance + Nedelman-Jia df),
 * `nca.buildCompositeProfile`, `nca.summarizeBootstrap` (BCa, min-n gated). The
 * derived parameters (Cmax, t½, CL, Vz, MRT, Vss, AUCinf) come from running the
 * **dense** `nca.computeNca` on the composite mean profile — exactly the
 * "build composite, then NCA on the composite" workflow the PKNCA reference
 * script uses. No NCA math is re-derived here; this module owns only the domain
 * orchestration (grouping, composite assembly, composite-BW resolution, the
 * honest-uncertainty surface, the warning stream).
 *
 * **Scientific contract (the honesty core — thesis #4):**
 *  - **Closed-form AUC is the headline** (AD-1): `AUClast` + the Holder/Bailer SE
 *    + Satterthwaite-t CI is the regulatory interval (NFR-04/05). The bootstrap
 *    AUC CI is *exploratory only* (AD-1, SF-2) and bounded to the nonlinear
 *    parameters as a primary use.
 *  - **Topology from the data** (AD-2): destructive / batch / serial is derived
 *    from the r_ij matrix and cross-checked against the declared label.
 *  - **Per-kg composite** (nca-cl-bw-normalised-units D-d, amending uc-04 SF-1): a
 *    `/kg` dose gives CL / Vz / Vss = dose per kg ÷ AUC of the composite (no weight
 *    needed; `SPARSE_PER_KG_WEIGHT_BASED` states the assumption); the absolute
 *    companions use the group-mean body weight (`SPARSE_COMPOSITE_BW`).
 *  - **Honest failure** (AD-6): the bootstrap is suppressed at small n; the
 *    closed-form CI is never gated (Holder holds to n=2).
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import {f64Column, numAt, ensureBoolBlq} from './adapter-utils';
import {nca} from '@datagrok-libraries/sci-comp';

import {DEFAULT_RULES, DEFAULT_ACCEPTANCE_RSQUARED} from './defaults';
import {evaluateAcceptance, reasonsCellValue, tagReasonPillsColumn, ACCEPT_REASONS_COL} from './acceptance';
import type {AcceptanceReading, AcceptanceThresholds} from './acceptance';
import {routeStringToCode, ROUTE_STRING_TO_CODE} from './profile-inputs';
import {ROUTE_CODE_TO_STRING} from './results-writer';
import {detectMixedCompoundGroups, mixedCompoundMessage, findColByPattern} from './column-autodetect';
import {
  CONC_TO_MG_L, BW_TO_KG, TIME_TO_HOUR, DOSE_TO_MG, CANONICAL_UNIT_TAGS, isPerKgDose,
  normaliseClVolumes, compositeDoseWeight,
} from './unit-conversion';
import type {BodyWeightContext, BodyWeightUnit} from './unit-conversion';
import {computeDecisionFingerprint} from './fingerprint';
import type {CompositeUnitInput} from './review-units';

// --- Public types --------------------------------------------------------

/** Which closed-form variance produced the AUC CI. */
export type AucCiMethod = 'Holder' | 'Bailer' | 'none';

export interface SparseNcaOptions {
  /** Concentration column. Default `'Conc'`. */
  readonly concCol?: string;
  /** Nominal-time column (the composite binning key). Required. */
  readonly nominalTimeCol: string;
  /** Animal-ID column. Optional — absent ⇒ destructive (Bailer) + a warning. */
  readonly animalIdCol?: string;
  /** BLQ-flag column (0/1 or bool). Optional — absent ⇒ all measurable. */
  readonly blqCol?: string;
  /** Dose column in mg/kg (per animal). Required for CL/Vz/Vss (per kg; no weight needed). */
  readonly doseCol?: string;
  /**
   * Absolute dose column in mg (per animal). **Opt-in BW bypass (AD-8):** when
   * set, CL/Vss use the group-mean absolute mg dose directly and the group-mean
   * body-weight normalization is skipped — for studies that carry an exact mg
   * dose and want to avoid the composite-BW convention. Takes precedence over
   * `doseCol` × group-mean BW.
   */
  readonly doseMgCol?: string;
  /** Body-weight column (kg). Optional — enables the absolute companions (`CL_ABS` …) via the
   *  group-mean BW, and the per-kg values of an absolute-mg dose. */
  readonly bodyWeightCol?: string;
  /** Sampling-type label column (declared topology; cross-checked). Optional. */
  readonly samplingTypeCol?: string;
  /**
   * Study-level declared topology, set when the user picks a non-`(auto)` Sampling
   * type in the mapping form. Overrides the per-group column-derived label and is
   * routed to `sparseAuc`'s existing AD-2 cross-check, so a wrong declaration
   * surfaces `SPARSE_TOPOLOGY_MISMATCH` (UC-04 AC-A2). `(auto)` ⇒ undefined ⇒ the
   * `samplingTypeCol` / `dominantLabel` path is used unchanged (zero regression).
   */
  readonly declaredTopology?: nca.SamplingTopology;
  /** Dose-group column the analysis strata are keyed on. Default = `doseCol`. */
  readonly groupCol?: string;
  /** Route label column or scalar. Default `'PO'` (extravascular). */
  readonly route?: string;
  /** Concentration unit (informational; CL/Vss inherit it). Default `'mg/L'`. */
  readonly concUnit?: nca.ConcentrationUnit;
  /** Scalar LLOQ in concentration units. Default 0. */
  readonly lloq?: number;
  /** λz acceptance thresholds for the composite ACCEPTANCE verdict (U2). Default: the shipped
   *  `R_acc` default, the default λz min points and %AUCextrap warn threshold. */
  readonly acceptance?: AcceptanceThresholds;
  /** Two-sided CI level. Default 0.95. */
  readonly ciLevel?: number;
  /** Bootstrap iterations. Default 1000. */
  readonly bootstrapIterations?: number;
  /** Bootstrap master seed. Default 12345. */
  readonly bootstrapSeed?: number;
  /** Calibrated bootstrap min-n per timepoint. Default 3 (sci-comp default). */
  readonly bootstrapMinN?: number;
  /** Compute on these rows of `observations` only (default: every row). The groups' `rows`
   *  index `observations` either way — {@link computeSparseStrata} passes one route's rows. */
  readonly rows?: ReadonlyArray<number>;
  /** The route qualifier of the group keys (`ROUTE=<route>|DOSE=…`), set only when the analysis
   *  holds more than one route (BUG-19 AD-3); absent / null ⇒ `DOSE=…`. */
  readonly routeKey?: string | null;
  /** Declared unit of the dose-group column, for the `DOSE_LEVEL` header. Default `'mg/kg'`. */
  readonly doseLevelUnit?: string;
}

export type SparseNcaWarningCode =
  | nca.SparseWarningCode
  | 'SPARSE_COMPOSITE_BW'
  | 'SPARSE_NO_DOSE_NO_CL'
  | 'SPARSE_BOOTSTRAP_N_TOO_SMALL'
  | 'SPARSE_BOOTSTRAP_AUC_EXPLORATORY'
  // BUG-19 data notes (AD-3 / AD-7 / AD-8 / S-4).
  | 'SPARSE_LLOQ_VARIES'
  | 'SPARSE_CONC_MISSING_DROPPED'
  | 'SPARSE_DOSE_SPELLING_MERGED'
  | 'SPARSE_BLQ_FLAG_WITH_VALUE'
  // nca-cl-bw-normalised-units AD-5 / AD-5b: the per-kg composite.
  | 'SPARSE_PER_KG_WEIGHT_BASED'
  | 'PER_KG_NOT_FORMED'
  | 'COMPOSITE_BW_PARTIAL'
  | 'COMPOSITE_NOT_WEIGHT_BASED';

export interface SparseNcaWarning {
  readonly code: SparseNcaWarningCode;
  readonly message: string;
  /** Group key the warning applies to, or '' for a table-level note. */
  readonly group: string;
}

/** A bootstrap interval reduced to the fields the view needs. */
export interface SparseBootstrapCi {
  readonly estimate: number;
  readonly ciLow: number;
  readonly ciHigh: number;
  readonly method: 'BCa' | 'percentile' | 'none';
  readonly suppressed: boolean;
  readonly suppressReason: string | null;
}

export interface SparseGroupResult {
  readonly groupKey: string;
  /** Route of the composite (`ROUTE_CODE_TO_STRING` spelling). */
  readonly route: string;
  /** The dose level as the source spells it (BUG-19 AD-3, S-6). `key` is the raw dose-group
   *  string the group key carries (`''` with no dose column); `value` its parsed number (NaN for a
   *  text label such as "Low") — the dose label, the sort and the decision fingerprint read it. */
  readonly doseLevel: {readonly key: string; readonly value: number};
  /** Member rows of the composite, indexing the frame passed to `computeSparseNca` (every member
   *  row, including a row AD-7 dropped from the pool). */
  readonly rows: ReadonlyArray<number>;
  /** The dose the CL / Vz / Vss compute reads, in canonical mg/kg; NaN for an absolute dose. */
  readonly doseMgKg: number;
  /** Group-mean absolute dose (mg): the absolute dose column's mean, or mg/kg × group-mean BW;
   *  NaN when neither is known. */
  readonly doseMg: number;
  readonly meanBodyWeightKg: number | null;
  /** How the ENGINE dose of this composite was formed (nca-cl-bw-normalised-units AD-5): `per_kg`
   *  for a `/kg` dose (the engine runs on mg/kg), `absolute` for an absolute-mg dose column. */
  readonly engineDoseNorm: 'absolute' | 'per_kg';
  /** The engine's own CL (L/h/kg on `per_kg`, L/h on `absolute`) — the implied-CL plausibility
   *  check reads THIS with {@link engineDoseNorm}, never the reported mL/min/kg `cl` (architect N1). */
  readonly engineCl: number;
  /** The weight the engine CL is divided by for the per-kg value on the `absolute` path
   *  (`compositeDoseWeight.perKgDivisor`), null when there is none — the plausibility check's BW. */
  readonly enginePerKgBodyWeightKg: number | null;
  readonly topology: nca.SamplingTopology;
  // Closed-form AUC headline.
  readonly aucLast: number;
  readonly aucSe: number;
  readonly aucDf: number;
  readonly aucCi: readonly [number, number];
  readonly aucCiMethod: AucCiMethod;
  // Derived (dense computeNca on the composite).
  readonly cmax: number;
  readonly tmax: number;
  readonly aucInf: number;
  readonly lambdaZ: number;
  readonly halfLife: number;
  /** Per kg: CL in mL/min/kg, Vz / Vss in L/kg (`normaliseClVolumes`); NaN when not formed. */
  readonly cl: number;
  readonly vz: number;
  readonly mrt: number;
  readonly vss: number;
  /** Absolute companions (L/h, L, L): the engine value on an absolute dose; per kg × group-mean BW
   *  on a `/kg` dose (NaN there without a weight). */
  readonly clAbs: number;
  readonly vzAbs: number;
  readonly vssAbs: number;
  readonly pctExtrap: number;
  // Fragility inputs for the honest-uncertainty surface (AD-9).
  /** Terminal-phase R² of the composite λz fit (`NaN` when λz not estimable). */
  readonly lambdaZR2: number;
  /** Adjusted R², span ratio and point count of the SAME composite λz fit, carried from
   *  `dense.provenance.lambdaZ` (already computed — no new math). They feed the composite's
   *  ACCEPTANCE verdict through the shared `evaluateAcceptance` predicate, so a pooled cohort is
   *  judged by the same criteria as a dense profile (nca-workflow-stages F-S4.8 / U2). `NaN`
   *  when λz was not estimable. */
  readonly lambdaZAdjR2: number;
  readonly spanRatio: number;
  readonly lambdaZNPoints: number;
  /** The composite's λz acceptance reading (verdict, flags, reasons) against the analysis's
   *  thresholds — the SAME predicate as the per-subject ACCEPTANCE column. */
  readonly acceptance: AcceptanceReading;
  /** Animals pooled into the composite, and the fewest / most sampled at any nominal time —
   *  the pooled-cohort cue (AC-S4.8.0: "composite of {n} animals ({min}–{max} per timepoint)"). */
  readonly nAnimals: number;
  readonly maxTimepointN: number;
  /** Fewest animals sampled at any nominal time (drives the n<4 fragility flag). */
  readonly minTimepointN: number;
  /** A singleton timepoint's variance was modeled (Nedelman var-vs-mean). */
  readonly varianceModeled: boolean;
  /** The bootstrap was suppressed (min-n gate / degenerate resamples). */
  readonly bootstrapSuppressed: boolean;
  // Bootstrap (nonlinear primary + exploratory AUC).
  readonly bootstrapAucExploratory: SparseBootstrapCi;
  readonly bootstrapCmax: SparseBootstrapCi;
  readonly bootstrapHalfLife: SparseBootstrapCi;
  readonly bootstrapCl: SparseBootstrapCi;
  readonly composite: nca.CompositeProfileResult;
  /** Raw per-animal points for the VR-07 scatter overlay. */
  readonly rawPoints: ReadonlyArray<{nominalTime: number; conc: number; blq: boolean}>;
  readonly warnings: SparseNcaWarning[];
}

export interface SparseNcaResult {
  readonly groups: SparseGroupResult[];
  /** One row per dose group: AUC headline + CI + method + derived + flags. */
  readonly parameterTable: DG.DataFrame;
  /** One row per (group × nominal time): mean/SD/%CV/n/%BLQ. */
  readonly compositeTable: DG.DataFrame;
  readonly warnings: SparseNcaWarning[];
}

/** What `computeSparseNca` returns: the result plus the dose groups that got no composite — every
 *  member sample had neither a concentration nor a BLQ flag (AD-7 dropped them all). Never silent:
 *  `computeSparseStrata` reports them as blocked strata (adding the route). */
export type SparseNcaRouteResult = SparseNcaResult & {
  readonly skipped: ReadonlyArray<Omit<SparseBlockedStratum, 'route'>>;
};

// --- Internal helpers ----------------------------------------------------

function boolAt(col: DG.Column | null, row: number): boolean {
  if (col === null || col.isNone(row)) return false;
  const v = col.get(row);
  return v === true || v === 1;
}

/** Display form of a nominal time or a concentration in a note. */
function fmtNum(v: number): string {
  return Number.isFinite(v) ? String(Number(v.toPrecision(4))) : 'NA';
}

/**
 * The ONE builder of a sparse composite's group key (BUG-19 AD-3; architect reuse 2): `DOSE=<raw>`,
 * `ROUTE=<route>|DOSE=<raw>` for an analysis with more than one route, `all` without a dose column.
 * Read by `computeSparseNca` and by the pipeline save's flag keys, so the two cannot drift.
 */
export function sparseCompositeKey(route: string | null, doseKey: string | null): string {
  const dose = doseKey === null ? 'all' : `DOSE=${doseKey}`;
  return route === null ? dose : `ROUTE=${route}|${dose}`;
}

/** One composite dose group: its raw key, parsed value, member rows and the spellings it pooled. */
export interface SparseDoseGroup {
  readonly key: string;
  readonly value: number;
  readonly rows: number[];
  readonly spellings: string[];
}

/**
 * The composite dose groups of `rows` (BUG-19 AD-3, plan review R1-F4). When every non-empty dose
 * value parses as a number, rows group by that NUMBER — the pipeline's rule — so "10" and "10.0"
 * are one composite, keyed by the lexicographically smallest spelling (independent of row order,
 * R2-N3). Otherwise (a text label such as "Low") rows group by the raw string. An empty dose cell
 * is its own group with key `''`. Groups come in first-seen order.
 */
export function sparseDoseGroups(col: DG.Column, rows: ReadonlyArray<number>): SparseDoseGroup[] {
  const rawAt = (i: number): string => (col.isNone(i) ? '' : col.getString(i));
  const numeric = rows.every((i) => rawAt(i).trim() === '' || Number.isFinite(numAt(col, i)));
  const byKey = new Map<string, {key: string; value: number; rows: number[]; spellings: string[]}>();
  for (const i of rows) {
    const raw = rawAt(i);
    const empty = raw.trim() === '';
    const value = empty ? NaN : numAt(col, i);
    const id = empty ? '\u0001empty' : (numeric ? String(value) : raw);
    let g = byKey.get(id);
    if (g === undefined) byKey.set(id, g = {key: raw, value, rows: [], spellings: []});
    g.rows.push(i);
    if (!g.spellings.includes(raw)) g.spellings.push(raw);
    if (raw < g.key) g.key = raw;
  }
  return [...byKey.values()];
}

/** Build the flat sci-comp SparseInput for one group's rows. An unflagged missing concentration is
 *  dropped from the pool (AD-7, user decision 2026-10-16 — PKNCA with the row removed); a
 *  BLQ-flagged row is always kept (`set-zero` reads it as 0). */
function buildSparseInput(
  rows: ReadonlyArray<number>,
  cols: {nom: DG.Column; conc: DG.Column; blq: DG.Column | null; animal: DG.Column | null},
  lloq: number,
): {
  input: nca.SparseInput; raw: Array<{nominalTime: number; conc: number; blq: boolean}>; kept: number[];
  dropped: number[]; flaggedWithValue: number;
} {
  const kept: number[] = [];
  const dropped: number[] = [];
  let flaggedWithValue = 0;
  for (const r of rows) {
    const c = numAt(cols.conc, r);
    const isBlq = boolAt(cols.blq, r);
    if (!Number.isFinite(c) && !isBlq) {dropped.push(r); continue;}
    if (isBlq && c > 0) flaggedWithValue++;
    kept.push(r);
  }
  const m = kept.length;
  const nominalTime = new Float64Array(m);
  const conc = new Float64Array(m);
  const blqMask = new Uint8Array(m);
  const animalId = cols.animal !== null ? new Int32Array(m) : null;
  const animalCode = new Map<string, number>();
  const raw: Array<{nominalTime: number; conc: number; blq: boolean}> = [];
  for (let i = 0; i < m; i++) {
    const r = kept[i];
    nominalTime[i] = numAt(cols.nom, r);
    conc[i] = numAt(cols.conc, r);
    const isBlq = boolAt(cols.blq, r);
    blqMask[i] = isBlq ? 1 : 0;
    if (animalId !== null) {
      const a = cols.animal!.isNone(r) ? `__r${r}` : cols.animal!.getString(r);
      let code = animalCode.get(a);
      if (code === undefined) {code = animalCode.size + 1; animalCode.set(a, code);}
      animalId[i] = code;
    }
    raw.push({nominalTime: nominalTime[i], conc: conc[i], blq: isBlq});
  }
  return {input: {nominalTime, conc, blqMask, lloq, animalId}, raw, kept, dropped, flaggedWithValue};
}

/** The `SPARSE_CONC_MISSING_DROPPED` sentence: the dropped count per nominal time, and every time
 *  left with no sample (AD-7 / R1-F1 / AC-7.3 / AC-7.6 — a late-time cluster stays visible). */
function droppedMessage(
  cols: {nom: DG.Column; conc: DG.Column; blq: DG.Column | null},
  rows: ReadonlyArray<number>, dropped: ReadonlyArray<number>,
): string {
  const count = (list: ReadonlyArray<number>): Map<number, number> => {
    const m = new Map<number, number>();
    for (const r of list) m.set(numAt(cols.nom, r), (m.get(numAt(cols.nom, r)) ?? 0) + 1);
    return m;
  };
  const byTime = count(dropped);
  const all = count(rows);
  // Measurable = kept, not BLQ-flagged, concentration > 0 (what can carry Tlast / the λz window).
  const measurable = count(rows.filter((r) => !boolAt(cols.blq, r) && numAt(cols.conc, r) > 0));
  const times = [...byTime.keys()].sort((a, b) => a - b);
  const gone = times.filter((t) => byTime.get(t) === all.get(t));
  const onlyBlq = times.filter((t) => byTime.get(t) !== all.get(t) && (measurable.get(t) ?? 0) === 0);
  const list = (ts: number[]): string => ts.map((t) => `${fmtNum(t)} h`).join(', ');
  return `${dropped.length} sample(s) with no concentration and no BLQ flag dropped from the composite: ` +
    `${times.map((t) => `${byTime.get(t)} at ${fmtNum(t)} h`).join(', ')}.` +
    (gone.length > 0 ? ` No sample left at ${list(gone)}: that time is not in the composite.` : '') +
    // AC-7.6: a time left with BLQ samples only can no longer be Tlast or in the λz window.
    (onlyBlq.length > 0 ? ` No measurable sample left at ${list(onlyBlq)}.` : '');
}

/**
 * λz rules for the composite NCA. The default `minRSquared: 0.85` is a hard gate
 * that *withholds* λz/AUCinf/CL when the terminal fit is poor — correct for a
 * dense per-subject profile, but wrong for a sparse composite: PKNCA reports the
 * best-fit λz on the composite even at low R² (e.g. the mouse 10 mg/kg composite
 * has a noisy terminal blip → R²≈0.40, yet PKNCA ships cl.obs≈0.076). Per thesis
 * #4 the honest move is to *ship the fragile value and flag it* (the FRAGILITY
 * `R²<0.90` flag), not hide it — so the composite path accepts the best-fit λz
 * (`minRSquared: 0`). Dense UC-01/02/03 rules are untouched.
 */
const SPARSE_COMPOSITE_RULES: nca.NcaRules = {
  ...DEFAULT_RULES,
  lambdaZ: {...DEFAULT_RULES.lambdaZ, minRSquared: 0},
};

/**
 * Run the dense `computeNca` on the composite mean profile of a sparse input —
 * the "build composite, then NCA on the composite" workflow (PKNCA parity). The
 * composite means feed a single pseudo-profile; `doseMg` is the dose per kg on a `/kg`
 * composite (CL/Vz/Vss then come out per kg) or the group-mean mg on an absolute one.
 */
function compositeNca(
  input: nca.SparseInput, doseMg: number, bwKg: number | null,
  routeCode: nca.RouteCode, concUnit: nca.ConcentrationUnit,
  prebuilt?: nca.CompositeProfileResult,
): nca.ComputeResult {
  // Reuse the composite already built by `sparseAuc` for the point estimate
  // (architect SIMPLIFY-1); the bootstrap path rebuilds per resample.
  const comp = prebuilt ?? nca.buildCompositeProfile(input, 'set-zero');
  const tps = comp.timepoints;
  const time = Float64Array.from(tps.map((t) => t.time));
  // Composite mean per timepoint; a fully-BLQ-majority timepoint contributes 0.
  const conc = Float64Array.from(tps.map((t) => (t.n > 0 && t.nBlq / t.n > 0.5 ? 0 : t.mean)));
  const blqMask = Uint8Array.from(tps.map((t) => (t.mean > 0 ? 0 : 1)));
  const profile: nca.ProfileInputs = {
    time, conc, blqMask, lloq: input.lloq,
    dose: doseMg, doseUnits: 'mg', concentrationUnits: concUnit, timeUnits: 'h',
    route: routeCode, infusionDuration: null, bodyWeight: bwKg,
  };
  return nca.computeNca(profile, SPARSE_COMPOSITE_RULES);
}

function toBootstrapCi(s: nca.BootstrapSummary): SparseBootstrapCi {
  return {
    estimate: s.estimate, ciLow: s.ci[0], ciHigh: s.ci[1],
    method: s.ciMethod, suppressed: s.suppressed, suppressReason: s.suppressReason,
  };
}

// --- Public entry point --------------------------------------------------

/**
 * Compute sparse / destructive-sampling NCA from a raw observations DataFrame.
 *
 * @throws if a required column (nominal time) is absent, or if a dose group holds animals of
 *   ≥ 2 compounds (a column named like a compound id varies within the group — BUG-07 pattern;
 *   the message is `mixedCompoundMessage`). The view shows the message as a banner.
 */
export function computeSparseNca(
  observations: DG.DataFrame,
  options: SparseNcaOptions,
): SparseNcaRouteResult {
  const concName = options.concCol ?? 'Conc';
  const nomName = options.nominalTimeCol;
  const groupName = options.groupCol ?? options.doseCol;
  const lloq = options.lloq ?? 0;
  const ciLevel = options.ciLevel ?? 0.95;
  const concUnit = options.concUnit ?? 'mg/L';
  const routeCode = routeStringToCode(options.route ?? 'PO');
  const iterations = options.bootstrapIterations ?? 1000;
  const seed = options.bootstrapSeed ?? 12345;
  const minN = options.bootstrapMinN ?? nca.DEFAULT_MIN_N_PER_TIMEPOINT;
  const thresholds: AcceptanceThresholds = options.acceptance ?? {
    lambdaAcceptanceRSquared: DEFAULT_ACCEPTANCE_RSQUARED,
    lambdaMinPoints: DEFAULT_RULES.lambdaZ.minPoints,
    extrapWarnPct: DEFAULT_RULES.extrapWarnPct,
  };

  const nomCol = observations.col(nomName);
  const concCol = observations.col(concName);
  if (nomCol === null) throw new Error(`computeSparseNca: nominal-time column '${nomName}' not found.`);
  if (concCol === null) throw new Error(`computeSparseNca: concentration column '${concName}' not found.`);
  const animalCol = options.animalIdCol !== undefined ? observations.col(options.animalIdCol) : null;
  const blqCol = options.blqCol !== undefined ? observations.col(options.blqCol) : null;
  const doseCol = options.doseCol !== undefined ? observations.col(options.doseCol) : null;
  const doseMgCol = options.doseMgCol !== undefined ? observations.col(options.doseMgCol) : null;
  const bwCol = options.bodyWeightCol !== undefined ? observations.col(options.bodyWeightCol) : null;
  const typeCol = options.samplingTypeCol !== undefined ? observations.col(options.samplingTypeCol) : null;
  const groupCol = groupName !== undefined ? observations.col(groupName) : null;

  const allRows = options.rows ?? Array.from({length: observations.rowCount}, (_, i) => i);
  const routeKey = options.routeKey ?? null;

  // --- Group rows by the dose level (single stratum if no group column) — AD-3.
  const doseGroups: SparseDoseGroup[] = groupCol !== null ? sparseDoseGroups(groupCol, allRows) :
    [{key: '', value: NaN, rows: [...allRows], spellings: []}];
  const keyOfRow = new Map<number, string>();
  for (const g of doseGroups) for (const r of g.rows) keyOfRow.set(r, groupCol !== null ? g.key : 'all');
  // BUG-07 pattern: the composite key has no compound dimension, so a dose group holding animals
  // of two compounds would average them into one composite — silently wrong. Refuse instead.
  const mixed = detectMixedCompoundGroups(observations, (i) => keyOfRow.get(i) ?? null,
    groupName !== undefined ? [groupName] : []);
  if (mixed.length > 0) {
    throw new Error(mixedCompoundMessage(mixed[0],
      groupCol !== null ? `Dose group '${mixed[0].group}'` : 'The study (no dose-group column)'));
  }

  const cols = {nom: nomCol, conc: concCol, blq: blqCol, animal: animalCol};
  const groups: SparseGroupResult[] = [];
  const skipped: Array<Omit<SparseBlockedStratum, 'route'>> = [];

  // Numeric on the parsed dose level, a text label (NaN) last, then the raw key (S-6) — the
  // lexical fallback no longer orders "10" before "2".
  const sortedGroups = [...doseGroups].sort((a, b) => {
    const an = Number.isNaN(a.value); const bn = Number.isNaN(b.value);
    if (an !== bn) return an ? 1 : -1;
    return (an ? 0 : a.value - b.value) || a.key.localeCompare(b.key);
  });
  for (const g of sortedGroups) {
    const groupKey = sparseCompositeKey(routeKey, groupCol !== null ? g.key : null);
    const {input, raw, kept, dropped, flaggedWithValue} = buildSparseInput(g.rows, cols, lloq);
    if (kept.length === 0) {
      skipped.push({doseKey: g.key, rows: g.rows,
        reason: 'no sample of this dose group has a concentration or a BLQ flag'});
      continue;
    }
    // Data notes that change the pool, so they stay in the decision fingerprint (AD-9).
    const dataNotes: SparseNcaWarning[] = [];
    if (g.spellings.length > 1) {
      dataNotes.push({code: 'SPARSE_DOSE_SPELLING_MERGED', group: groupKey,
        message: `Dose values ${g.spellings.map((s) => `"${s}"`).join(', ')} are one dose level, ` +
          `pooled into one composite keyed "${g.key}".`});
    }
    if (dropped.length > 0) {
      dataNotes.push({code: 'SPARSE_CONC_MISSING_DROPPED', group: groupKey,
        message: droppedMessage(cols, g.rows, dropped)});
    }
    if (flaggedWithValue > 0) {
      dataNotes.push({code: 'SPARSE_BLQ_FLAG_WITH_VALUE', group: groupKey,
        message: `${flaggedWithValue} BLQ-flagged sample(s) carry a reported concentration; the BLQ rule ` +
          'reads them as below the LLOQ (0).'});
    }

    // Declared topology label for the AD-2 cross-check. A user-supplied study-level
    // `declaredTopology` (non-`(auto)` Sampling type) takes precedence over the
    // per-group column-derived most-common label.
    const declared = options.declaredTopology ??
      (typeCol !== null ? dominantLabel(typeCol, g.rows) : undefined);

    // --- Closed-form AUC (the headline).
    const sa = nca.sparseAuc(input, {ciLevel, blqRule: 'set-zero', declaredTopology: declared});
    const aucCiMethod: AucCiMethod = Number.isFinite(sa.ci[0]) ?
      (sa.topology === 'destructive' ? 'Bailer' : 'Holder') : 'none';

    // --- Dose for CL/Vz/Vss (nca-cl-bw-normalised-units AD-5, amending uc-04 SF-1).
    // Over the POOLED rows only (`kept`): an animal whose sample AD-7 dropped contributes no
    // concentration, so it contributes no body weight or dose either — PKNCA with the row removed
    // (user decision at the BUG-19 review gate, 2026-10-16).
    //  - `/kg` dose: the engine runs on the dose per kg with no weight (the dense no-BW
    //    pass-through) → per-kg CL = D_kg / AUC, no weight needed; CL_ABS = per kg × group-mean BW.
    //  - absolute-mg dose column: the engine runs on the group-mean mg dose; per kg = mean dose /
    //    mean BW / AUC over the SAME members (`compositeDoseWeight`, shared with the NONMEM pooler).
    const meanBwKg = bwCol !== null ? groupMean(bwCol, kept) : null;
    const firstDose = doseCol !== null ? numAt(doseCol, kept[0]) : NaN;
    const doseMgKg = Number.isFinite(firstDose) ? firstDose : NaN;
    const absolute = doseMgCol !== null ? compositeDoseWeight(kept.map((r) => ({
      dose: numAt(doseMgCol, r), bw: bwCol !== null ? numAt(bwCol, r) : NaN}))) : null;
    const usingAbsolute = absolute !== null && Number.isFinite(absolute.meanDose);
    const engineDoseNorm: 'absolute' | 'per_kg' = usingAbsolute ? 'absolute' : 'per_kg';
    const engineDose = usingAbsolute ? absolute!.meanDose : doseMgKg;
    const doseMg = usingAbsolute ? absolute!.meanDose :
      ((meanBwKg !== null && Number.isFinite(doseMgKg)) ? doseMgKg * meanBwKg : NaN);
    // The weight the per-kg rule reads: the group-mean BW on a `/kg` dose (it only forms CL_ABS),
    // the same-member divisor on an absolute dose.
    const perKgBw: BodyWeightContext = {
      bodyWeightKg: usingAbsolute ? absolute!.perKgDivisor : meanBwKg,
      bodyWeightRangeKg: null, perKgDose: !usingAbsolute,
    };

    const warnings: SparseNcaWarning[] = [...sa.warnings.map((w) => ({
      code: w.code, message: w.message, group: groupKey,
    })), ...dataNotes];
    if (usingAbsolute) {
      for (const w of absolute!.warnings)
        warnings.push({code: w.code as SparseNcaWarningCode, message: w.message, group: groupKey});
      if (absolute!.perKgDivisor === null) {
        warnings.push({code: 'PER_KG_NOT_FORMED', group: groupKey,
          message: 'No body weight in this group: CL, Vz and Vss per kg are not reported for the ' +
            'absolute-mg dose. The absolute values are.'});
      }
    } else if (Number.isFinite(doseMgKg)) {
      warnings.push({code: 'SPARSE_PER_KG_WEIGHT_BASED', group: groupKey,
        message: 'Per-kg composite CL, Vz and Vss assume every animal received the same dose per kg ' +
          '(dose per kg ÷ AUC of the mean concentration curve).'});
      if (meanBwKg !== null) {
        warnings.push({
          code: 'SPARSE_COMPOSITE_BW',
          message: `CL abs., Vz abs. and Vss abs. use the group-mean body weight (${meanBwKg.toFixed(3)} kg) — ` +
            `the composite pools animals of differing BW (Phoenix/PKNCA convention). ` +
            `Supply an absolute mg dose (doseMgCol) to bypass body-weight normalization.`,
          group: groupKey,
        });
      }
    } else {
      warnings.push({
        code: 'SPARSE_NO_DOSE_NO_CL',
        message: 'No dose for this group; CL, Vz and Vss are not reported ' +
          '(AUC/Cmax/t½ are unaffected).',
        group: groupKey,
      });
    }

    // --- Derived params from the composite (dense computeNca on the composite mean).
    // Reuse `sa.composite` (already built by sparseAuc) — no second composite build. A per-kg run
    // carries no weight (the engine's per-kg pass-through); an absolute run keeps the group mean.
    const engineBw = usingAbsolute ? meanBwKg : null;
    const dense = compositeNca(input, engineDose, engineBw, routeCode, concUnit, sa.composite);
    const dv = dense.values;
    const norm = normaliseClVolumes(dv, engineDoseNorm, perKgBw);
    const lz = dense.provenance.lambdaZ;
    const lambdaZR2 = lz !== null ? lz.rSquared : NaN;
    const lambdaZAdjR2 = lz !== null ? lz.adjRSquared : NaN;
    const spanRatio = lz !== null ? lz.spanRatio : NaN;
    const lambdaZNPoints = lz !== null ? lz.pointsUsed.length : NaN;
    const acceptance = evaluateAcceptance(
      {spanRatio, adjR2: lambdaZAdjR2, nPoints: lambdaZNPoints, pctExtrap: dv.pctExtrap}, thresholds);
    const maxTimepointN = sa.composite.timepoints.reduce((acc, t) => Math.max(acc, t.n), 0);
    const nAnimals = animalCol !== null ?
      new Set(kept.map((r) => (animalCol.isNone(r) ? '' : animalCol.getString(r)))).size :
      kept.length;
    // Fragility inputs (AD-9): fewest animals at any nominal time; modeled variance.
    const minTimepointN = sa.composite.timepoints.reduce(
      (acc, t) => Math.min(acc, t.n), Number.POSITIVE_INFINITY);
    const varianceModeled = sa.warnings.some((w) => w.code === 'SPARSE_VARIANCE_MODELED');

    // --- Bootstrap: nonlinear primary (Cmax, t½, CL) + exploratory AUC (SF-2).
    const bootOpts = {iterations, masterSeed: seed, ciLevel, minNPerTimepoint: minN};
    const denseVals = (r: nca.SparseInput): nca.ParameterValues =>
      compositeNca(r, engineDose, engineBw, routeCode, concUnit).values;
    const bAuc = nca.summarizeBootstrap(input, (r) => nca.sparseAuc(r).auc, bootOpts);
    const bCmax = nca.summarizeBootstrap(input, (r) => denseVals(r).cmax, bootOpts);
    const bHl = nca.summarizeBootstrap(input, (r) => denseVals(r).halfLife, bootOpts);
    // The CL interval shares the point estimate's basis and unit (per kg, mL/min/kg — AD-5).
    const bCl = nca.summarizeBootstrap(input,
      (r) => normaliseClVolumes(denseVals(r), engineDoseNorm, perKgBw).cl, bootOpts);

    if (bAuc.suppressed) {
      warnings.push({
        code: 'SPARSE_BOOTSTRAP_N_TOO_SMALL',
        message: `Bootstrap suppressed: ${bAuc.suppressReason}. The closed-form ` +
          `(Holder/Bailer) CI stands as the AUC interval.`,
        group: groupKey,
      });
    } else {
      warnings.push({
        code: 'SPARSE_BOOTSTRAP_AUC_EXPLORATORY',
        message: 'The bootstrap AUC CI is exploratory — not for regulatory use. The ' +
          'closed-form Holder/Bailer CI is the AUC headline (NFR-04/05).',
        group: groupKey,
      });
    }

    groups.push({
      groupKey, route: ROUTE_CODE_TO_STRING[routeCode], doseLevel: {key: g.key, value: g.value}, rows: g.rows,
      doseMgKg, doseMg, meanBodyWeightKg: meanBwKg, topology: sa.topology,
      engineDoseNorm, engineCl: dv.cl, enginePerKgBodyWeightKg: usingAbsolute ? absolute!.perKgDivisor : null,
      aucLast: sa.auc, aucSe: sa.se, aucDf: sa.df, aucCi: sa.ci, aucCiMethod,
      cmax: dv.cmax, tmax: dv.tmax, aucInf: dv.aucInf, lambdaZ: dv.lambdaZ,
      halfLife: dv.halfLife, cl: norm.cl, vz: norm.vz, mrt: dv.mrt, vss: norm.vss,
      clAbs: norm.clAbs, vzAbs: norm.vzAbs, vssAbs: norm.vssAbs,
      pctExtrap: dv.pctExtrap, lambdaZR2, lambdaZAdjR2, spanRatio, lambdaZNPoints, acceptance,
      nAnimals, maxTimepointN,
      minTimepointN: Number.isFinite(minTimepointN) ? minTimepointN : 0,
      varianceModeled, bootstrapSuppressed: bAuc.suppressed,
      bootstrapAucExploratory: toBootstrapCi(bAuc),
      bootstrapCmax: toBootstrapCi(bCmax),
      bootstrapHalfLife: toBootstrapCi(bHl),
      bootstrapCl: toBootstrapCi(bCl),
      composite: sa.composite, rawPoints: raw, warnings,
    });
  }

  return {...assembleSparseResult(groups, options.doseLevelUnit ?? 'mg/kg'), skipped};
}

/**
 * The result tables over composite groups — one route's (`computeSparseNca`) or several routes'
 * merged (`computeSparseForRecipe`). The same builders either way, so a merged result keeps the
 * one schema the presentation maps and the BUG-12 guard check.
 */
export function assembleSparseResult(groups: SparseGroupResult[], doseLevelUnit: string): SparseNcaResult {
  return {
    groups,
    parameterTable: buildParameterTable(groups, doseLevelUnit),
    compositeTable: buildCompositeTable(groups, doseLevelUnit),
    warnings: groups.flatMap((g) => g.warnings),
  };
}

function dominantLabel(col: DG.Column, rows: number[]): nca.SamplingTopology | undefined {
  const counts = new Map<string, number>();
  for (const r of rows) {
    if (col.isNone(r)) continue;
    const v = col.getString(r).trim().toLowerCase();
    if (v === '') continue;
    counts.set(v, (counts.get(v) ?? 0) + 1);
  }
  let best = ''; let bestN = 0;
  for (const [k, c] of counts) if (c > bestN) {best = k; bestN = c;}
  if (best === 'destructive' || best === 'batch' || best === 'serial') return best;
  return undefined;
}

function groupMean(col: DG.Column, rows: number[]): number | null {
  let s = 0; let c = 0;
  for (const r of rows) {const v = numAt(col, r); if (Number.isFinite(v)) {s += v; c++;}}
  return c > 0 ? s / c : null;
}

// --- Output table builders ----------------------------------------------

/**
 * Short, human-readable label for a sparse advisory code — shared by the toast
 * summary (view) and the persistent NOTES column, so both stay concise and
 * consistent. The full sentence lives in `result.warnings` (audit / provenance).
 */
export function shortSparseNote(code: string): string {
  switch (code) {
  case 'SPARSE_COMPOSITE_BW': return 'CL abs. on group-mean BW';
  case 'SPARSE_NO_DOSE_NO_CL': return 'CL/Vss N/A (no dose)';
  case 'SPARSE_PER_KG_WEIGHT_BASED': return 'per kg assumes weight-based dosing';
  case 'PER_KG_NOT_FORMED': return 'per-kg CL N/A (no BW)';
  case 'COMPOSITE_BW_PARTIAL': return 'per kg on animals with BW';
  case 'COMPOSITE_NOT_WEIGHT_BASED': return 'dose not ∝ BW (per kg approx.)';
  case 'SPARSE_BOOTSTRAP_AUC_EXPLORATORY': return 'bootstrap AUC exploratory';
  case 'SPARSE_BOOTSTRAP_N_TOO_SMALL': return 'bootstrap suppressed (low N)';
  case 'SPARSE_TOPOLOGY_MISMATCH': return 'topology ≠ declared label';
  case 'SPARSE_ANIMAL_ID_ABSENT_ASSUMED_DESTRUCTIVE': return 'no animal-ID → destructive';
  case 'SPARSE_VARIANCE_MODELED': return 'variance modeled (n=1)';
  case 'SPARSE_TERMINAL_OVEREST': return 'terminal AUC may overestimate';
  case 'SPARSE_DF_UNAVAILABLE': return 'df unavailable';
  case 'SPARSE_LLOQ_VARIES': return 'LLOQ varies (largest used)';
  case 'SPARSE_CONC_MISSING_DROPPED': return 'missing conc dropped';
  case 'SPARSE_DOSE_SPELLING_MERGED': return 'dose spellings merged';
  case 'SPARSE_BLQ_FLAG_WITH_VALUE': return 'BLQ flag with a value → 0';
  default: return code;
  }
}

/** Distinct short notes for one group, `;`-joined (persistent NOTES column). */
function groupNotes(g: SparseGroupResult): string {
  const seen = new Set<string>();
  const out: string[] = [];
  for (const w of g.warnings) {
    const s = shortSparseNote(w.code);
    if (!seen.has(s)) {seen.add(s); out.push(s);}
  }
  return out.join('; ');
}

/** Per-parameter fragility flags (AD-9): the honest-uncertainty surface. Returns
 *  a `;`-joined flag string, or `ok` when none fire.
 *
 *  **No λz-fit token** (scientist-ratified): the
 *  former hard-coded plain-`R²<0.90` token is replaced by the composite's ACCEPTANCE verdict
 *  (adj-R² vs `R_acc`, span, N points, %AUCextrap — the same predicate as a dense profile), so
 *  there is one λz acceptance criterion, not two that can disagree. The remaining flags are
 *  sampling-design uncertainty, informational only. */
function fragilityFlags(g: SparseGroupResult): string {
  const flags: string[] = [];
  if (g.minTimepointN > 0 && g.minTimepointN < 4) flags.push('n<4');
  if (!Number.isFinite(g.aucDf) || g.aucDf < 3) flags.push('df<3');
  if (g.varianceModeled) flags.push('variance-modeled');
  if (g.bootstrapSuppressed) flags.push('bootstrap-suppressed');
  return flags.length > 0 ? flags.join('; ') : 'ok';
}

/** Bootstrap AUC CI rendered as ONE subordinate string (AC-C2 / R2 NF-R2-2):
 *  the suppressed state is explicit text, never a blank or numeric CI. */
function bootstrapAucCiText(g: SparseGroupResult): string {
  const b = g.bootstrapAucExploratory;
  if (b.suppressed) return 'Bootstrap CI unavailable (N too small)';
  if (!Number.isFinite(b.ciLow) || !Number.isFinite(b.ciHigh)) return 'n/a';
  return `${b.ciLow.toFixed(3)}–${b.ciHigh.toFixed(3)} (${b.method})`;
}

/** The `DOSE_LEVEL` column (AD-4 / S-5): the dose level in the DECLARED unit of the dose-group
 *  column, headed `Dose (<unit>)` by its `label` tag (the presentation map reads the tag). */
function doseLevelColumn(values: ArrayLike<number>, unit: string): DG.Column {
  return f64Column('DOSE_LEVEL', values).setTag('label', `Dose (${unit})`);
}

/** Sparse parameter-table columns tagged with a unit (the per-kg family + companions + BW). */
const SPARSE_UNIT_TAGGED: ReadonlySet<string> =
  new Set(['CL', 'VZ', 'VSS', 'CL_ABS', 'VZ_ABS', 'VSS_ABS', 'BW_KG']);

function buildParameterTable(groups: ReadonlyArray<SparseGroupResult>, doseLevelUnit: string): DG.DataFrame {
  const m = groups.length;
  const mk = (): Float64Array => new Float64Array(m);
  const doseOut = mk(); const aucOut = mk(); const seOut = mk(); const dfOut = mk();
  const ciLoOut = mk(); const ciHiOut = mk(); const cmaxOut = mk(); const tmaxOut = mk();
  const aucInfOut = mk(); const hlOut = mk(); const clOut = mk(); const vzOut = mk();
  const mrtOut = mk(); const vssOut = mk(); const extrapOut = mk(); const r2Out = mk();
  const clAbsOut = mk(); const vzAbsOut = mk(); const vssAbsOut = mk(); const bwOut = mk();
  const adjR2Out = mk(); const acceptOut = new Array<string>(m); const reasonsOut = new Array<string>(m);
  const methodOut = new Array<string>(m); const topoOut = new Array<string>(m);
  const fragilityOut = new Array<string>(m); const bootCiOut = new Array<string>(m);
  const notesOut = new Array<string>(m);

  groups.forEach((g, i) => {
    doseOut[i] = g.doseLevel.value; aucOut[i] = g.aucLast; seOut[i] = g.aucSe; dfOut[i] = g.aucDf;
    ciLoOut[i] = g.aucCi[0]; ciHiOut[i] = g.aucCi[1]; cmaxOut[i] = g.cmax; tmaxOut[i] = g.tmax;
    aucInfOut[i] = g.aucInf; hlOut[i] = g.halfLife; clOut[i] = g.cl; vzOut[i] = g.vz;
    mrtOut[i] = g.mrt; vssOut[i] = g.vss; extrapOut[i] = g.pctExtrap; r2Out[i] = g.lambdaZR2;
    clAbsOut[i] = g.clAbs; vzAbsOut[i] = g.vzAbs; vssAbsOut[i] = g.vssAbs;
    bwOut[i] = g.meanBodyWeightKg ?? NaN;
    adjR2Out[i] = g.lambdaZAdjR2; acceptOut[i] = g.acceptance.verdict;
    reasonsOut[i] = reasonsCellValue(g.acceptance);
    methodOut[i] = g.aucCiMethod; topoOut[i] = g.topology;
    fragilityOut[i] = fragilityFlags(g); bootCiOut[i] = bootstrapAucCiText(g);
    notesOut[i] = groupNotes(g);
  });

  // The closed-form AUClast + Holder/Bailer CI is the headline (leftmost, numeric).
  // The bootstrap AUC CI is the LAST column, a single annotated string — visually
  // and positionally SUBORDINATE so it can't be read as a co-primary regulatory
  // interval (AC-C2 / AD-1 / R1-plan NB-1). FRAGILITY surfaces the per-parameter
  // honest-uncertainty flags (AD-9 / thesis #4).
  const cols: DG.Column[] = [
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'ROUTE', groups.map((g) => g.route)),
    doseLevelColumn(doseOut, doseLevelUnit),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'TOPOLOGY', topoOut),
    f64Column('AUC_LAST', aucOut).setTag('label', 'AUClast (closed-form, headline)'),
    f64Column('AUC_SE', seOut),
    f64Column('AUC_DF', dfOut),
    f64Column('AUC_CI_LOW', ciLoOut),
    f64Column('AUC_CI_HIGH', ciHiOut),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'AUC_CI_METHOD', methodOut),
    f64Column('CMAX', cmaxOut),
    f64Column('TMAX', tmaxOut),
    f64Column('AUC_INF', aucInfOut),
    f64Column('HALF_LIFE', hlOut),
    f64Column('CL', clOut),
    f64Column('VZ', vzOut),
    f64Column('MRT', mrtOut),
    f64Column('VSS', vssOut),
    f64Column('PCT_AUCEXTRAP', extrapOut),
    f64Column('LAMBDAZ_R2', r2Out).setTag('label', 'λz R²'),
    f64Column('LAMBDAZ_ADJ_R2', adjR2Out).setTag('label', 'adj. R² (λz)'),
    // The composite's λz ACCEPTANCE verdict (U2) — the same predicate and vocabulary as the
    // per-subject results grid's ACCEPTANCE column.
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'ACCEPTANCE', acceptOut),
    // One pill per violated criterion — the same column as the per-subject grid's `Violated`.
    tagReasonPillsColumn(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, ACCEPT_REASONS_COL, reasonsOut)),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'FRAGILITY', fragilityOut),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'AUC_BOOT_CI', bootCiOut)
      .setTag('label', 'AUC bootstrap CI — exploratory, not for regulatory use'),
    // Persistent advisory stream (thesis #4): the full sentences are in
    // result.warnings; this column keeps the short labels visible after the
    // transient toast is gone.
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'NOTES', notesOut),
    // Absolute companions + the group-mean weight (AD-5, architect N4) — the dense frame's names,
    // so the λz-dependent set and Exclude λz reach them. BW_KG is not λz-dependent.
    f64Column('CL_ABS', clAbsOut).setTag('label', 'CL abs.'),
    f64Column('VZ_ABS', vzAbsOut).setTag('label', 'Vz abs.'),
    f64Column('VSS_ABS', vssAbsOut).setTag('label', 'Vss abs.'),
    f64Column('BW_KG', bwOut).setTag('label', 'Group-mean BW'),
  ];
  // The per-kg family and its companions carry the dense frame's unit tags, from the one map
  // (AD-5). The other sparse columns stay as they were.
  for (const col of cols)
    if (SPARSE_UNIT_TAGGED.has(col.name)) col.setTag('units', CANONICAL_UNIT_TAGS[col.name]);
  const table = DGNS.DataFrame.fromColumns(cols);
  table.name = 'Sparse NCA';
  return table;
}

function buildCompositeTable(groups: ReadonlyArray<SparseGroupResult>, doseLevelUnit: string): DG.DataFrame {
  const routeOut: string[] = [];
  const doseOut: number[] = []; const timeOut: number[] = []; const nOut: number[] = [];
  const nBlqOut: number[] = []; const meanOut: number[] = []; const sdOut: number[] = [];
  const cvOut: number[] = []; const pctBlqOut: number[] = [];
  for (const g of groups) {
    for (const tp of g.composite.timepoints) {
      routeOut.push(g.route);
      doseOut.push(g.doseLevel.value); timeOut.push(tp.time); nOut.push(tp.n); nBlqOut.push(tp.nBlq);
      meanOut.push(tp.mean); sdOut.push(tp.sd);
      // %CV is undefined at mean=0 (0/0). Render it blank (FLOAT_NULL), not the
      // raw text "NaN" — BLQ%=100 / N=BLQ already convey the all-BLQ state.
      cvOut.push(Number.isFinite(tp.cvPct) ? tp.cvPct : DGNS.FLOAT_NULL);
      pctBlqOut.push(tp.n > 0 ? 100 * tp.nBlq / tp.n : NaN);
    }
  }
  const cols: DG.Column[] = [
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'ROUTE', routeOut),
    doseLevelColumn(doseOut, doseLevelUnit),
    f64Column('NOMINAL_TIME', timeOut),
    f64Column('N', nOut),
    f64Column('N_BLQ', nBlqOut),
    f64Column('MEAN_CONC', meanOut),
    f64Column('SD', sdOut),
    f64Column('CV_PCT', cvOut).setTag('label', '%CV'),
    f64Column('PCT_BLQ', pctBlqOut).setTag('label', '%BLQ'),
  ];
  const table = DGNS.DataFrame.fromColumns(cols);
  table.name = 'Sparse Composite Profile';
  return table;
}

// --- The app's sparse panel compute, view-free ----------------------------
// Moved from `views/nca-mapping-pane.ts` (nca-pipeline-save-rework F1): the headless pipeline
// computes the SAME composites the app computes on open, to derive the saved stage — so the
// compute lives here and both callers share it (review 5A: "move the needed part into the
// adapter rather than copy it"). Since BUG-19 both the app and `ncaCompute` run the ONE
// canonicalization + strata function below (`computeSparseStrata`).

/** The declared units of a sparse compute's columns — REQUIRED on every path (BUG-19, guardrail
 *  "every path that calls the engine canonicalizes every unit"): a missing unit is a compile error,
 *  never a silent hours default. */
export interface SparseUnits {
  readonly time: nca.TimeUnit;
  readonly concentration: nca.ConcentrationUnit;
  /** Unit of the dose column: per-kg (`mg/kg`, `µg/kg`) or absolute (`mg`, `µg`, `ng`). */
  readonly dose: nca.DoseUnit;
  readonly bodyWeight: BodyWeightUnit;
}

/** What {@link computeSparseStrata} reads: role column names, units, LLOQ, route source. */
export interface SparseStrataSpec {
  readonly concCol: string;
  readonly nominalTimeCol: string;
  readonly animalIdCol: string | null;
  /** BLQ flag column: bool, 0/1 or text (`Y`, `BLQ`, `<LLOQ`, … — read by `ensureBoolBlq`). */
  readonly blqCol: string | null;
  /** The RAW dose column: it keys the composites; a canonical copy feeds CL / Vz / Vss. */
  readonly doseCol: string | null;
  readonly bodyWeightCol: string | null;
  readonly lloq: {readonly colName: string | null; readonly scalar: number};
  /** Per-row route codes (`DERIVED_ROUTE_CODE_COL`, int32, null cell = the scalar) or none. */
  readonly route: {readonly codesCol: string | null; readonly scalar: string};
  readonly units: SparseUnits;
  readonly samplingTypeCol?: string | null;
  readonly declaredTopology?: nca.SamplingTopology;
  readonly acceptance: AcceptanceThresholds;
  readonly bootstrapIterations?: number;
}

/** One computed composite: a route × dose-level cell. `rows` index the input frame. */
export interface SparseStratumCell {
  readonly route: string;
  readonly doseKey: string;
  readonly rows: ReadonlyArray<number>;
  readonly group: SparseGroupResult;
}

/** A route × dose-level cell that was not computed, with the reason. */
export interface SparseBlockedStratum {
  readonly route: string;
  readonly doseKey: string;
  readonly rows: ReadonlyArray<number>;
  readonly reason: string;
}

export type SparseStrata =
  | {readonly ok: false; readonly reason: string}
  | {readonly ok: true; readonly cells: SparseStratumCell[]; readonly blocked: SparseBlockedStratum[];
     /** The declared dose unit (the `DOSE_LEVEL` header). */
     readonly doseLevelUnit: string};

/** Why a sparse IV-infusion stratum is not computed (one wording for the app and the pipeline). */
export const SPARSE_INFUSION_BLOCKED_REASON =
  'sparse infusion is not supported (no infusion duration in the composite engine)';

/** A derived Float64 column `name` = `src` × `factor` (a platform null stays null) — the dense
 *  `convertToCanonical` arithmetic: multiply by the factor, so sparse and dense times are the same
 *  floats (plan review R1-F7). */
function addScaled(df: DG.DataFrame, src: string, name: string, factor: number): string {
  const c = df.col(src)!;
  df.columns.add(DGNS.Column.fromFloat64Array(name, Float64Array.from({length: df.rowCount},
    (_, i) => c.isNone(i) ? DGNS.FLOAT_NULL : numAt(c, i) * factor)));
  return name;
}

/** Up to `max` distinct non-empty cells of `col` that do not parse as a finite number. */
function nonNumericTokens(col: DG.Column, max = 3): string[] {
  const out: string[] = [];
  for (let i = 0; i < col.length && out.length < max; i++) {
    if (col.isNone(i)) continue;
    const s = col.getString(i).trim();
    if (s !== '' && !Number.isFinite(numAt(col, i)) && !out.includes(s)) out.push(s);
  }
  return out;
}

/**
 * THE sparse composite compute, shared by the app (`computeSparseForRecipe`) and `ncaCompute`
 * (BUG-19 AD-2). On a clone of `df` it converts every unit to the canonical basis (conc → mg/L,
 * nominal time → h, BW → kg, dose → mg/kg or absolute mg; LLOQ = the largest value in the LLOQ
 * column, else the scalar, × the conc factor), splits the rows by route code, refuses an
 * IV-infusion stratum, and runs one `computeSparseNca` per remaining route with the RAW dose column
 * as the group key. Keys carry the route only when the rows hold more than one route code,
 * counting blocked strata (R2-N5). Returns one cell per route × dose level; `rows` index `df`
 * (the clone is row-aligned — AD-7 drops rows inside the engine-input builder, never here).
 *
 * @throws what `computeSparseNca` throws (a dose group pooling two compounds — BUG-07).
 */
export function computeSparseStrata(df: DG.DataFrame, spec: SparseStrataSpec): SparseStrata {
  const u = spec.units;
  const fConc = CONC_TO_MG_L[u.concentration];
  if (!Number.isFinite(fConc))
    return {ok: false, reason: `sparse: molar concentration unit '${u.concentration}' is not supported`};
  const concSrc = df.col(spec.concCol);
  if (concSrc === null) return {ok: false, reason: `concentration column '${spec.concCol}' not found`};
  if (df.col(spec.nominalTimeCol) === null)
    return {ok: false, reason: `nominal-time column '${spec.nominalTimeCol}' not found`};
  // AD-7 / R1-F1: a text token ("BLQ", "ND", "<0.1") never reaches the pool as a number.
  const bad = nonNumericTokens(concSrc);
  if (bad.length > 0) {
    return {ok: false, reason: `concentration column '${spec.concCol}' is not numeric: ` +
      bad.map((t) => `"${t}"`).join(', ')};
  }
  const scalarCode = ROUTE_STRING_TO_CODE[spec.route.scalar];

  const work = df.clone();
  const conc = addScaled(work, spec.concCol, '~SPARSE_CONC_MG_L', fConc);
  const nominal = addScaled(work, spec.nominalTimeCol, '~SPARSE_NOMINAL_H', TIME_TO_HOUR[u.time]);
  const has = (name: string | null): name is string => name !== null && work.col(name) !== null;
  const bw = has(spec.bodyWeightCol) ?
    addScaled(work, spec.bodyWeightCol, '~SPARSE_BW_KG', BW_TO_KG[u.bodyWeight]) : undefined;
  const perKg = isPerKgDose(u.dose);
  const dose = has(spec.doseCol) ?
    addScaled(work, spec.doseCol, perKg ? '~SPARSE_DOSE_MG_KG' : '~SPARSE_DOSE_MG', DOSE_TO_MG[u.dose]) : undefined;
  const blq = has(spec.blqCol) ? ensureBoolBlq(work, spec.blqCol) : null;

  // LLOQ (S-4, R1-F6): the LARGEST value in the column — deterministic and conservative. Inert
  // under the shipped `set-zero` composite rules, so it changes no number; a varying column is
  // disclosed by a note that names the value used.
  let lloq = spec.lloq.scalar;
  let lloqNote: string | null = null;
  if (has(spec.lloq.colName)) {
    const c = work.col(spec.lloq.colName)!;
    const vals = new Set<number>();
    for (let i = 0; i < c.length; i++) {const v = numAt(c, i); if (Number.isFinite(v)) vals.add(v);}
    if (vals.size > 0) lloq = [...vals].reduce((a, b) => Math.max(a, b), -Infinity);
    if (vals.size > 1) {
      lloqNote = `The LLOQ column holds ${vals.size} distinct values; the composite uses the largest, ` +
        `${fmtNum(lloq * fConc)} mg/L.`;
    }
  }
  lloq = Number.isFinite(lloq) ? lloq * fConc : 0;

  // Route strata, in first-seen order. A blank route cell takes the scalar (the pipeline's
  // `defaultRoute` rule).
  const codes = spec.route.codesCol !== null ? work.col(spec.route.codesCol) : null;
  const byRoute = new Map<number, number[]>();
  for (let i = 0; i < work.rowCount; i++) {
    const code = codes !== null && !codes.isNone(i) ? Number(codes.get(i)) : scalarCode;
    if (code === undefined) return {ok: false, reason: `route '${spec.route.scalar}' is not a known route`};
    let list = byRoute.get(code);
    if (list === undefined) byRoute.set(code, list = []);
    list.push(i);
  }
  const multiRoute = byRoute.size > 1;
  const rawDose = has(spec.doseCol) ? work.col(spec.doseCol) : null;

  const cells: SparseStratumCell[] = [];
  const blocked: SparseBlockedStratum[] = [];
  for (const [code, rows] of byRoute) {
    const route = ROUTE_CODE_TO_STRING[code];
    if (route === 'IV-infusion') {
      const groups = rawDose !== null ? sparseDoseGroups(rawDose, rows) : [{key: '', rows}];
      for (const g of groups)
        blocked.push({route, doseKey: g.key, rows: g.rows, reason: SPARSE_INFUSION_BLOCKED_REASON});
      continue;
    }
    const res = computeSparseNca(work, {
      concCol: conc,
      nominalTimeCol: nominal,
      ...(has(spec.animalIdCol) ? {animalIdCol: spec.animalIdCol} : {}),
      ...(blq !== null ? {blqCol: blq} : {}),
      ...(dose !== undefined ? (perKg ? {doseCol: dose} : {doseMgCol: dose}) : {}),
      ...(rawDose !== null ? {groupCol: spec.doseCol!} : {}),
      ...(bw !== undefined ? {bodyWeightCol: bw} : {}),
      ...(has(spec.samplingTypeCol ?? null) ? {samplingTypeCol: spec.samplingTypeCol!} : {}),
      ...(spec.declaredTopology !== undefined ? {declaredTopology: spec.declaredTopology} : {}),
      route,
      concUnit: 'mg/L',
      lloq,
      acceptance: spec.acceptance,
      ...(spec.bootstrapIterations !== undefined ? {bootstrapIterations: spec.bootstrapIterations} : {}),
      rows,
      routeKey: multiRoute ? route : null,
      doseLevelUnit: u.dose,
    });
    for (const g of res.groups) {
      if (lloqNote !== null) g.warnings.push({code: 'SPARSE_LLOQ_VARIES', message: lloqNote, group: g.groupKey});
      cells.push({route, doseKey: g.doseLevel.key, rows: g.rows, group: g});
    }
    for (const s of res.skipped) blocked.push({route, doseKey: s.doseKey, rows: s.rows, reason: s.reason});
  }
  return {ok: true, cells, blocked, doseLevelUnit: u.dose};
}

/** The sparse panel's mapping for one compute — what the app derives from its Setup inputs. */
export interface SparseComputeMapping {
  subjectName: string | null; timeName: string | null; concName: string | null;
  /** The RAW dose column (the composite key); its unit is `units.dose`. */
  doseColName: string | null;
  bwColName: string | null;
  /** REQUIRED: the declared units of time / conc / dose / BW (BUG-19 — no hours default). */
  units: SparseUnits;
  /** Per-row route codes (`DERIVED_ROUTE_CODE_COL`) or null, plus the route scalar (S-1). */
  route: {codesCol: string | null; scalar: string};
  /** LLOQ column or scalar (in `units.concentration`). */
  lloq: {colName: string | null; scalar: number};
  /** The MAPPED BLQ-flag role; null ⇒ the BLQ name heuristic (AD-8). */
  blqColName: string | null;
  /** Explicit nominal-time role from the Sparse-design section; null ⇒ heuristic. */
  nominalTimeName?: string | null;
  /** Explicit animal-ID role; null ⇒ Subject-then-heuristic fallback. */
  animalIdName?: string | null;
  /** Declared topology from the Sampling-type choice; null ⇒ derive-from-data. */
  declaredTopology?: 'destructive' | 'batch' | 'serial' | null;
  /** The live λz acceptance thresholds — the composite ACCEPTANCE verdict reads the SAME
   *  `R_acc` as the per-subject grid (nca-workflow-stages F-S4.8 / U2). */
  acceptance: AcceptanceThresholds;
  /** Bootstrap iterations override (the Confirm dry run uses a short run — no bootstrap
   *  field enters a decision fingerprint). Default: the sparse module's own default. */
  bootstrapIterations?: number;
}

/** True when the frame is a sparse / destructive composite design: a distinct
 *  nominal-time column exists AND ≥ 2 animals share at least one nominal time
 *  (pooling present). A dense serial design (one animal per profile, full time
 *  course) does not pool and is handled by the normal pipeline. */
export function isSparseDesign(df: DG.DataFrame, nominalName: string, animalName: string): boolean {
  const nom = df.col(nominalName);
  const animal = df.col(animalName);
  if (nom === null || animal === null) return false;
  const perTime = new Map<number, Set<string>>();
  const n = df.rowCount;
  for (let i = 0; i < n; i++) {
    if (nom.isNone(i)) continue;
    const t = Number(nom.get(i));
    const a = animal.isNone(i) ? `__r${i}` : animal.getString(i);
    let s = perTime.get(t);
    if (s === undefined) {s = new Set<string>(); perTime.set(t, s);}
    s.add(a);
  }
  for (const s of perTime.values()) if (s.size >= 2) return true;
  return false;
}

/**
 * The COMPUTE half of the sparse panel — no view, no dock, no notice (nca-workflow-stages
 * F-S3.5): shared by the app's sparse panel, its Confirm dry run and the pipeline save, so all
 * compute the SAME composites. `result` is null when the frame is not a sparse design or the
 * compute failed (`error` then carries the message).
 */
export function computeSparseForRecipe(
  df: DG.DataFrame,
  mapping: SparseComputeMapping,
): {result: SparseNcaResult | null; blocked: SparseBlockedStratum[]; error: string | null; elapsedMs: number} {
  const {subjectName, concName} = mapping;
  // The nominal-time role is the sparse trigger; it must be DISTINCT from the
  // observed time used by the dense pipeline. An explicit mapping (Sparse-design
  // section) overrides the name heuristic; null falls back to the heuristic so a
  // conventionally-named dataset (toggle off) keeps the prior behaviour.
  const nominalName = mapping.nominalTimeName ?? findColByPattern(df, /nominal.*time/i);
  const animalName = mapping.animalIdName ?? subjectName ??
    findColByPattern(df, /animal|mouse|subject/i);
  if (nominalName === null || animalName === null || concName === null ||
      !isSparseDesign(df, nominalName, animalName))
    return {result: null, blocked: [], error: null, elapsedMs: 0};
  // AD-8: the MAPPED BLQ role, else the name heuristic (no dataset loses its flag).
  const blqName = mapping.blqColName ?? findColByPattern(df, /^blq$|blq.?flag|^blq_/i);
  const samplingTypeName = findColByPattern(df, /sampling.?type|^stype$/i);

  const t0 = (typeof performance !== 'undefined' ? performance.now() : 0);
  let strata: SparseStrata;
  try {
    strata = computeSparseStrata(df, {
      concCol: concName, nominalTimeCol: nominalName, animalIdCol: animalName, blqCol: blqName,
      doseCol: mapping.doseColName, bodyWeightCol: mapping.bwColName, lloq: mapping.lloq,
      route: mapping.route, units: mapping.units, samplingTypeCol: samplingTypeName,
      // User-declared topology (Sampling-type choice) overrides the column-derived
      // label and is cross-checked against the data (SPARSE_TOPOLOGY_MISMATCH).
      ...(mapping.declaredTopology != null ? {declaredTopology: mapping.declaredTopology} : {}),
      acceptance: mapping.acceptance,
      ...(mapping.bootstrapIterations !== undefined ? {bootstrapIterations: mapping.bootstrapIterations} : {}),
    });
  } catch (e) {
    return {result: null, blocked: [], error: String((e as Error)?.message ?? e), elapsedMs: 0};
  }
  if (!strata.ok) return {result: null, blocked: [], error: strata.reason, elapsedMs: 0};
  const result = assembleSparseResult(strata.cells.map((c) => c.group), strata.doseLevelUnit);
  const elapsedMs = (typeof performance !== 'undefined' ? performance.now() : 0) - t0;
  return {result, blocked: strata.blocked, error: null, elapsedMs};
}

/** Bootstrap-derived warning codes are STOCHASTIC-adjacent (resampling) — kept out of every
 *  decision fingerprint (AD-7 determinism premise: the Confirm dry run must predict exactly). */
const SPARSE_BOOTSTRAP_CODE = /^SPARSE_BOOTSTRAP_/;
/** Informational codes no composite number depends on (BUG-19 AD-9): also kept out of the
 *  fingerprint. The rule: a code enters the fingerprint only if the composite's numbers depend on
 *  it — so `SPARSE_CONC_MISSING_DROPPED`, `SPARSE_DOSE_SPELLING_MERGED` and (user decision (a),
 *  the rows are zeroed) `SPARSE_BLQ_FLAG_WITH_VALUE` stay in. */
const SPARSE_NUMBER_FREE_CODES: ReadonlySet<string> = new Set(['SPARSE_LLOQ_VARIES']);
// (The per-kg notes — nca-cl-bw-normalised-units AD-5 — are left out by `computeDecisionFingerprint`
// itself, `BASIS_NOTE_CODES`, for dense profiles and composites alike.)

/** UC-04 group composites as review units, each with its decision fingerprint (members + fit
 *  summary + verdict + the warning codes the numbers depend on). Empty when the sparse panel did
 *  not run. */
export function sparseCompositeUnits(
  result: SparseNcaResult | null,
): Array<CompositeUnitInput & {fingerprint: string}> {
  if (result === null) return [];
  return result.groups.map((g) => {
    const fingerprint = computeDecisionFingerprint({
      points: g.rawPoints.map((p) => [String(p.nominalTime), p.blq ? 'BLQ' : String(p.conc)] as const),
      // The raw dose level (S-6): `String(value)` is the pre-BUG-19 `String(doseMgKg)` byte for
      // byte on a per-kg study, and stays raw when the compute dose is canonicalized.
      dose: String(g.doseLevel.value), route: '', pointsUsedTimes: [],
      lambdaZ: g.lambdaZ, adjR2: g.lambdaZAdjR2, aucInf: g.aucInf, pctExtrap: g.pctExtrap,
      verdict: g.acceptance.verdict, flagCodes: g.acceptance.flagCodes,
      warnings: g.warnings
        .filter((w) => !SPARSE_BOOTSTRAP_CODE.test(w.code) && !SPARSE_NUMBER_FREE_CODES.has(w.code))
        .map((w) => ({code: w.code, severity: 'info'})),
    });
    return {
      groupKey: g.groupKey, acceptance: g.acceptance,
      measurableTimepoints: g.composite.timepoints.filter((t) => t.mean > 0).length,
      fingerprint,
    };
  });
}
