/**
 * Toxicokinetics (TK) aggregates — repeat-dose accumulation ratio, steady-state
 * theoretical ceiling, NOAEL exposure/safety margin, and sex fold-difference
 * (FR-420..424).
 *
 * **Why a NEW module, not `aggregateNcaResults` (mirror of `bioavailability.ts`).**
 * `aggregateNcaResults` does *within-group* arithmetic (mean / SD / %CV / N) via
 * the native `groupBy().aggregate()`. The accumulation ratio is intrinsically a
 * *cross-StudyDay* ratio (Day N ÷ Day 1 of the *same* dose×sex cell) — the native
 * aggregator cannot divide cell A by cell B. Exactly the reason `bioavailability.ts`
 * (a cross-route ratio) is also separate. `aggregateNcaResults` is untouched (it is
 * the #11 contract this cycle *consumes*; per-row-route AD-6 precedent).
 *
 * **Correct layering (not vendored math).** The per-profile primitives this
 * consumes (AUCτ via AUClast, Cmax, λz, %extrap) are sci-comp's. The TK aggregates
 * are cross-profile ratios over DataFrame strata — DataFrame orchestration, which
 * CLAUDE.md places in this package. No core math is re-derived.
 *
 * **Ratified scientific contract (USER-DECISION D1..D4, blueprint gate):**
 *  - **Rac(AUC) Formula A (default, D1)** = AUCτ,ss / AUCτ,Day1 (CDISC SEND RACAUC).
 *    AUCτ = AUClast under the `tlast == τ` guard (D4). Formula B
 *    (AUCτ,ss / AUCinf,Day1, WinNonlin) is **opt-in only** where the Day-1 λz
 *    quality gates pass (λz present **and** adj.R² ≥ 0.90 **and** %AUCextrap < 20%);
 *    otherwise the cell **falls back to A and is labelled** — never N/A, never a
 *    number resting on an unreliable tail. The classical "Accumulation Index" is
 *    the *theoretical* `1/(1−e^(−λz·τ))` ({@link accumulationIndexTheoretical},
 *    `RAC_THEO`), NOT Formula B — Formula B is labelled `"B (AUCτ,ss/AUCinf,Day1)"`.
 *  - **τ is the dosing interval**, NOT a StudyDay span (R1/F-3). The caller resolves
 *    it (user study-design input, defaulting to `max(NominalTime_hr)`); this module
 *    receives the number.
 *  - **RAC_THEO (D1/P2)** computed **only** where the Day-1 λz passes the same gates;
 *    NaN-with-tooltip otherwise (honest uncertainty, thesis #4). Exact for IV-bolus
 *    1-compartment; for oral/EV it *approximates* AUCinf(Day1)/AUCτ(Day1).
 *  - **Exposure margin (D2)** = exposure-based default (rat NOAEL AUCτ ÷ *empirical*
 *    human AUC, ICH M3(R2)) — NaN when no human AUC is provided (never invented) —
 *    plus HED dose-based (BSA/Km, FDA 2005) and an allometrically-scaled CL
 *    projection, all **per-sex with the worst-case flagged**. NOAEL is a user input,
 *    never auto-derived from PK.
 *  - **Sex fold (D-/G8)** = `max(M,F)/min(M,F)` (≥ 1) + a direction label; carries
 *    N and %CV; **no significance claim** at composite N.
 *  - **Mixed-frame guards inherited** (per-row-route / cl-vz-convention): a cell
 *    mixing `CL_BASIS` or `DOSE_NORM` is **refused** (NaN + warning), not divided.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import {f64Column, numAt} from './adapter-utils';
import {BW_UNIT, CL_PER_KG_UNIT, L_PER_H_TO_ML_PER_MIN} from './unit-conversion';

// --- Gate / default constants (cited, conservative) ----------------------

/** λz-acceptance bar for Formula B and RAC_THEO (FR-223; matches the WinNonlin
 *  default; formal regulatory citation tracked as RG-5 — direction conservative). */
export const ADJ_R2_GATE = 0.90;
/** %AUCextrap(Day1) bar for Formula B and RAC_THEO (FR-226). */
export const PCT_EXTRAP_GATE = 20;
/** Default Rac(Cmax)/Rac(AUC) divergence bounds — an UNCITED screening heuristic
 *  (R1/F-6), exposed as an option, not a hard constant. */
export const DEFAULT_DIVERGENCE_BOUNDS: readonly [number, number] = [0.7, 1.43];
/** Default observed-vs-theoretical RAC_SS_FLAG bounds on `RAC_AUC / RAC_THEO` —
 *  an UNCITED screening heuristic (exposed as an option, not a hard constant;
 *  same treatment as {@link DEFAULT_DIVERGENCE_BOUNDS}). */
export const DEFAULT_SS_FLAG_BOUNDS: readonly [number, number] = [0.8, 1.25];
/** FDA 2005 Table 1 Km factors (BSA/Km HED scaling). */
export const DEFAULT_RAT_KM = 6;
export const DEFAULT_HUMAN_KM = 37;
/** Reference dose (mg/kg) for the projected-margin rat CL (the lowest non-vehicle
 *  dose; CL is ≈ dose-independent here). */
export const DEFAULT_CL_REF_DOSE_MG_KG = 10;
/** Standard reference rat weight (allometric anchor) and exponent: CL/kg scales
 *  with BW^(−0.25) (CL scales with BW^0.75). */
export const DEFAULT_RAT_BW_KG = 0.25;
export const ALLOMETRIC_EXPONENT = -0.25;

// --- Warning channel (mirrors BioavailabilityWarning) --------------------

export type TkWarningCode =
  | 'RAC_TLAST_NE_TAU'
  | 'RAC_TLAST_UNVERIFIED'
  | 'RAC_DIVERGENCE'
  | 'MIXED_CL_BASIS_AGGREGATE'
  | 'MIXED_DOSE_NORM_AGGREGATE'
  | 'NO_NOAEL_INPUT'
  | 'NO_HUMAN_AUC_INPUT'
  | 'NO_CL_REFERENCE_FOR_SEX'
  | 'MARGIN_NON_ALLOMETRIC'
  // nca-cl-bw-normalised-units D-m: the margin's CL input contract.
  | 'MARGIN_CL_NOT_PER_KG'
  | 'MARGIN_CL_ROWS_DROPPED';

/** One cell-level warning surfaced alongside a TK aggregate table. */
export interface TkWarning {
  readonly code: TkWarningCode;
  readonly message: string;
  /** The cell key the warning applies to (e.g. `'SEX=F|DOSE=100'`), or '' when
   *  the warning is table-level (e.g. NO_NOAEL_INPUT). */
  readonly cell: string;
}

// --- Shared helpers ------------------------------------------------------

/** Arithmetic mean over finite values; NaN if none (matches the
 *  `aggregateNcaResults` / VAL-01 arithmetic-mean convention). */
function mean(values: ReadonlyArray<number>): number {
  let sum = 0;
  let k = 0;
  for (const v of values) if (Number.isFinite(v)) {sum += v; k++;}
  return k === 0 ? NaN : sum / k;
}

/** Distinct non-null string values of `col` over the given rows (mixed-basis
 *  detection). Returns an empty set when the column is absent. */
function distinctOverRows(col: DG.Column | null, rows: ReadonlyArray<number>): Set<string> {
  const seen = new Set<string>();
  if (col === null) return seen;
  for (const r of rows) if (!col.isNone(r)) seen.add(col.getString(r));
  return seen;
}

// =====================================================================
// AD-1 / AD-2 — Accumulation ratio + theoretical steady-state ceiling
// =====================================================================

export interface AccumulationOptions {
  /** Dosing interval τ (hours). REQUIRED — τ is the dosing interval, NOT a
   *  StudyDay span (R1/F-3). The caller resolves it (user input, defaulting to
   *  `max(NominalTime_hr)`). */
  readonly tau: number;
  /** Per-subject dose column (the within-cell pivot key). Required. */
  readonly doseCol: string;
  /** StudyDay column (the cross-day pivot dimension). Default `'StudyDay'`. */
  readonly studyDayCol?: string;
  /** Stratification keys — sex is mandatory; append covariates as needed.
   *  Default `['SEX']`. */
  readonly stratifyCols?: ReadonlyArray<string>;
  /** Cmax column. Default `'CMAX'`. */
  readonly cmaxCol?: string;
  /** AUCτ column (AUClast under the tlast==τ guard). Default `'AUC_LAST'`. */
  readonly aucTauCol?: string;
  /** AUCinf column (Formula B denominator). Default `'AUC_INF'`. */
  readonly aucInfCol?: string;
  /** Day-1 λz column (Formula B gate + RAC_THEO). Default `'LAMBDA_Z'`. */
  readonly lambdaZCol?: string;
  /** Day-1 adj.R² column (gate). Default `'LAMBDAZ_ADJ_R2'`. */
  readonly adjR2Col?: string;
  /** Day-1 %AUCextrap column (gate). Default `'PCT_EXTRAP'`. */
  readonly pctExtrapCol?: string;
  /** Per-profile last-sampling-time column for the `tlast == τ` guard (D4). When
   *  absent from the frame the guard cannot run and emits RAC_TLAST_UNVERIFIED.
   *  Default `'TLAST'`. */
  readonly tlastCol?: string;
  /** Absolute tolerance (hours) for `tlast == τ`. Default `0.01`. */
  readonly tlastTolHr?: number;
  /** Baseline day. Default = min observed StudyDay. */
  readonly day1?: number;
  /** Steady-state day ("Day N"). Default = max observed StudyDay (data-driven,
   *  never literal 28). */
  readonly dayN?: number;
  /** Rac formula: `'A'` (default, AUCτ/AUCτ) or `'B'` (opt-in, AUCτ,ss/AUCinf,Day1
   *  where the Day-1 λz gate passes — else falls back to A). */
  readonly formula?: 'A' | 'B';
  /** Divergence bounds for ratio-of-ratios (Cmax-Rac ÷ AUC-Rac). Default
   *  {@link DEFAULT_DIVERGENCE_BOUNDS}. */
  readonly divergenceBounds?: readonly [number, number];
  /** Bounds for the observed-vs-theoretical RAC_SS_FLAG: `[lo, hi]` on
   *  `RAC_AUC / RAC_THEO`. Below `lo` ⇒ "not-yet-SS / autoinduction"; above `hi`
   *  ⇒ "time-dependent kinetics"; in-band ⇒ "consistent with theoretical". An
   *  UNCITED screening heuristic (exposed, not a hard constant — same treatment
   *  as {@link DEFAULT_DIVERGENCE_BOUNDS}). Default {@link DEFAULT_SS_FLAG_BOUNDS}. */
  readonly ssFlagBounds?: readonly [number, number];
}

export interface AccumulationResult {
  /** Per-cell table: stratify cols, `DOSE`, `RAC_AUC`, `RAC_CMAX`, `RAC_FORMULA`
   *  (per-cell label), `RAC_THEO`, `RAC_SS_FLAG` (text), `RAC_DIVERGENCE` (bool). */
  readonly table: DG.DataFrame;
  /** Resolved Day-1 / Day-N used (data-driven). */
  readonly day1: number;
  readonly dayN: number;
  readonly warnings: TkWarning[];
}

/**
 * Closed-form steady-state accumulation index (the WinNonlin/Phoenix
 * "Accumulation Index"): `1 / (1 − e^(−λz·τ))`. Exact only for IV-bolus
 * 1-compartment kinetics; for oral/EV it approximates AUCinf(Day1)/AUCτ(Day1).
 * Returns NaN for a non-positive / non-finite λz or τ.
 */
export function accumulationIndexTheoretical(lambdaZ: number, tau: number): number {
  if (!Number.isFinite(lambdaZ) || !Number.isFinite(tau) || lambdaZ <= 0 || tau <= 0)
    return NaN;
  const denom = 1 - Math.exp(-lambdaZ * tau);
  return denom > 0 ? 1 / denom : NaN;
}

/** True when a Day-1 cell's λz is reliable enough for Formula B / RAC_THEO. */
function dayOneGatePasses(lambdaZ: number, adjR2: number, pctExtrap: number): boolean {
  return Number.isFinite(lambdaZ) && lambdaZ > 0 &&
    Number.isFinite(adjR2) && adjR2 >= ADJ_R2_GATE &&
    Number.isFinite(pctExtrap) && pctExtrap < PCT_EXTRAP_GATE;
}

/**
 * Compute per-cell accumulation ratios (Day N ÷ Day 1) from a per-subject NCA
 * results frame stratified by `{…stratify, dose}` across StudyDay.
 *
 * @throws if a required column (dose / studyDay / cmax / aucτ / any stratify
 *   key) is absent from `results`.
 */
export function computeAccumulationRatios(
  results: DG.DataFrame,
  options: AccumulationOptions,
): AccumulationResult {
  const tau = options.tau;
  const doseName = options.doseCol;
  const dayName = options.studyDayCol ?? 'StudyDay';
  const stratifyCols = options.stratifyCols ?? ['SEX'];
  const cmaxName = options.cmaxCol ?? 'CMAX';
  const aucTauName = options.aucTauCol ?? 'AUC_LAST';
  const aucInfName = options.aucInfCol ?? 'AUC_INF';
  const lzName = options.lambdaZCol ?? 'LAMBDA_Z';
  const adjR2Name = options.adjR2Col ?? 'LAMBDAZ_ADJ_R2';
  const pctExtrapName = options.pctExtrapCol ?? 'PCT_EXTRAP';
  const tlastName = options.tlastCol ?? 'TLAST';
  const tlastTol = options.tlastTolHr ?? 0.01;
  const formula = options.formula ?? 'A';
  const [divLo, divHi] = options.divergenceBounds ?? DEFAULT_DIVERGENCE_BOUNDS;
  const [ssLo, ssHi] = options.ssFlagBounds ?? DEFAULT_SS_FLAG_BOUNDS;

  const required = [doseName, dayName, cmaxName, aucTauName, ...stratifyCols];
  for (const name of required) {
    if (results.col(name) === null)
      throw new Error(`computeAccumulationRatios: required column '${name}' not found in results.`);
  }

  const doseCol = results.col(doseName)!;
  const dayCol = results.col(dayName)!;
  const cmaxCol = results.col(cmaxName)!;
  const aucTauCol = results.col(aucTauName)!;
  const aucInfCol = results.col(aucInfName);
  const lzCol = results.col(lzName);
  const adjR2Col = results.col(adjR2Name);
  const pctExtrapCol = results.col(pctExtrapName);
  const tlastCol = results.col(tlastName);
  const clBasisCol = results.col('CL_BASIS');
  const doseNormCol = results.col('DOSE_NORM');
  const stratCols = stratifyCols.map((s) => results.col(s)!);
  const n = results.rowCount;

  // Resolve day1 / dayN data-driven.
  let minDay = Infinity;
  let maxDay = -Infinity;
  for (let i = 0; i < n; i++) {
    const d = numAt(dayCol, i);
    if (!Number.isFinite(d)) continue;
    if (d < minDay) minDay = d;
    if (d > maxDay) maxDay = d;
  }
  const day1 = options.day1 ?? minDay;
  const dayN = options.dayN ?? maxDay;

  // Group rows into cells keyed by {…stratify, dose}; record per-day row lists.
  interface Cell {
    stratumVals: string[];
    dose: number;
    day1Rows: number[];
    dayNRows: number[];
    allRows: number[];
  }
  const cells = new Map<string, Cell>();
  for (let i = 0; i < n; i++) {
    const dose = numAt(doseCol, i);
    if (!Number.isFinite(dose)) continue;
    const day = numAt(dayCol, i);
    const stratumVals = stratifyCols.map((s, j) => stratCols[j].isNone(i) ? '' : stratCols[j].getString(i));
    const key = `${stratumVals.join('|')}|D=${dose}`;
    let cell = cells.get(key);
    if (cell === undefined) {
      cell = {stratumVals, dose, day1Rows: [], dayNRows: [], allRows: []};
      cells.set(key, cell);
    }
    cell.allRows.push(i);
    if (day === day1) cell.day1Rows.push(i);
    else if (day === dayN) cell.dayNRows.push(i);
  }

  const cellList = [...cells.values()].sort((a, b) =>
    a.stratumVals.join('|').localeCompare(b.stratumVals.join('|')) || a.dose - b.dose);

  const warnings: TkWarning[] = [];
  const m = cellList.length;
  const stratOut = stratifyCols.map(() => new Array<string>(m));
  const doseOut = new Float64Array(m);
  const racAucOut = new Float64Array(m);
  const racCmaxOut = new Float64Array(m);
  const formulaOut = new Array<string>(m);
  const racTheoOut = new Float64Array(m);
  const ssFlagOut = new Array<string>(m);
  const divergenceOut = new Array<boolean>(m);

  for (let r = 0; r < m; r++) {
    const cell = cellList[r];
    const cellKey = stratifyCols.map((s, j) => `${s}=${cell.stratumVals[j]}`).join('|') + `|DOSE=${cell.dose}`;
    for (let j = 0; j < stratifyCols.length; j++) stratOut[j][r] = cell.stratumVals[j];
    doseOut[r] = cell.dose;

    // --- Mixed-basis guard (AC-5): refuse a cell that mixes CL_BASIS/DOSE_NORM.
    const mixedClBasis = distinctOverRows(clBasisCol, cell.allRows).size > 1;
    const mixedDoseNorm = distinctOverRows(doseNormCol, cell.allRows).size > 1;
    if (mixedClBasis || mixedDoseNorm) {
      racAucOut[r] = NaN; racCmaxOut[r] = NaN; racTheoOut[r] = NaN;
      formulaOut[r] = 'refused (mixed basis)';
      ssFlagOut[r] = '';
      divergenceOut[r] = false;
      if (mixedClBasis) {
        warnings.push({code: 'MIXED_CL_BASIS_AGGREGATE', cell: cellKey,
          message: `Cell ${cellKey} mixes CL/CL/F basis; Rac refused (NaN) — not pooled.`});
      }
      if (mixedDoseNorm) {
        warnings.push({code: 'MIXED_DOSE_NORM_AGGREGATE', cell: cellKey,
          message: `Cell ${cellKey} mixes absolute/per-kg dose normalization; Rac refused (NaN) — not pooled.`});
      }
      continue;
    }

    // --- Cmax ratio (always computable; CMAX is never NA here).
    const cmaxD1 = mean(cell.day1Rows.map((i) => numAt(cmaxCol, i)));
    const cmaxDN = mean(cell.dayNRows.map((i) => numAt(cmaxCol, i)));
    racCmaxOut[r] = (Number.isFinite(cmaxD1) && cmaxD1 > 0) ? cmaxDN / cmaxD1 : NaN;

    // --- tlast == τ guard (D4): both Day-1 and Day-N profiles must sample to τ.
    let tlastOk = true;
    if (tlastCol === null) {
      tlastOk = false; // cannot verify
      warnings.push({code: 'RAC_TLAST_UNVERIFIED', cell: cellKey,
        message: `No '${tlastName}' column; cannot verify tlast == τ (${tau} h) for ${cellKey}. ` +
          `Rac(AUC) assumes AUClast ≡ AUCτ — verify the sampling reaches τ.`});
    } else {
      const tlastVals = [...cell.day1Rows, ...cell.dayNRows].map((i) => numAt(tlastCol, i));
      for (const tl of tlastVals)
        if (!Number.isFinite(tl) || Math.abs(tl - tau) > tlastTol) {tlastOk = false; break;}

      if (!tlastOk && tlastVals.length > 0) {
        warnings.push({code: 'RAC_TLAST_NE_TAU', cell: cellKey,
          message: `Cell ${cellKey}: tlast ≠ τ (${tau} h); sampling stops before the dosing interval. ` +
            `Rac(AUC) = NaN (AUClast would under-report AUCτ).`});
      }
    }

    // --- AUC ratio (Formula A default; B opt-in where Day-1 gate passes).
    const aucTauD1 = mean(cell.day1Rows.map((i) => numAt(aucTauCol, i)));
    const aucTauDN = mean(cell.dayNRows.map((i) => numAt(aucTauCol, i)));

    // Day-1 λz quality (single composite cell → first Day-1 row drives the gate).
    const d1Row = cell.day1Rows.length > 0 ? cell.day1Rows[0] : -1;
    const lzD1 = d1Row >= 0 && lzCol !== null ? numAt(lzCol, d1Row) : NaN;
    const adjR2D1 = d1Row >= 0 && adjR2Col !== null ? numAt(adjR2Col, d1Row) : NaN;
    const extrapD1 = d1Row >= 0 && pctExtrapCol !== null ? numAt(pctExtrapCol, d1Row) : NaN;
    const gatePasses = dayOneGatePasses(lzD1, adjR2D1, extrapD1);

    let racAuc = NaN;
    let label: string;
    if (formula === 'B') {
      if (gatePasses) {
        const aucInfD1 = aucInfCol !== null ? mean(cell.day1Rows.map((i) => numAt(aucInfCol, i))) : NaN;
        if (Number.isFinite(aucTauDN) && Number.isFinite(aucInfD1) && aucInfD1 > 0)
          racAuc = aucTauDN / aucInfD1;
        label = 'B (AUCτ,ss/AUCinf,Day1)';
      } else {
        // Fall back to A — never N/A, never a number on an unreliable tail.
        if (Number.isFinite(aucTauDN) && Number.isFinite(aucTauD1) && aucTauD1 > 0)
          racAuc = aucTauDN / aucTauD1;
        label = 'A (fallback; B gate failed)';
      }
    } else {
      if (Number.isFinite(aucTauDN) && Number.isFinite(aucTauD1) && aucTauD1 > 0)
        racAuc = aucTauDN / aucTauD1;
      label = 'A (AUCτ/AUCτ)';
    }
    if (!tlastOk && tlastCol !== null) racAuc = NaN; // hard guard only when verifiable
    racAucOut[r] = racAuc;
    formulaOut[r] = label;

    // --- RAC_THEO (AD-2): gated closed form from Day-1 λz; else NaN.
    const racTheo = gatePasses ? accumulationIndexTheoretical(lzD1, tau) : NaN;
    racTheoOut[r] = racTheo;

    // --- Steady-state interpretation flag (text, not a number).
    if (Number.isFinite(racTheo) && Number.isFinite(racAuc) && racTheo > 0) {
      const ratio = racAuc / racTheo;
      if (ratio > ssHi) ssFlagOut[r] = 'observed ≫ theoretical — possible time-dependent kinetics';
      else if (ratio < ssLo) ssFlagOut[r] = 'observed ≪ theoretical — not-yet-SS or autoinduction';
      else ssFlagOut[r] = 'consistent with theoretical';
    } else
      ssFlagOut[r] = ''; // λz unfit / R² < 0.90 / extrap ≥ 20% — Rac_theo not reliable


    // --- Divergence flag (Cmax-Rac vs AUC-Rac, configurable bounds).
    let diverges = false;
    if (Number.isFinite(racCmaxOut[r]) && Number.isFinite(racAuc) && racAuc > 0) {
      const rr = racCmaxOut[r] / racAuc;
      diverges = rr < divLo || rr > divHi;
      if (diverges) {
        warnings.push({code: 'RAC_DIVERGENCE', cell: cellKey,
          message: `Cell ${cellKey}: Rac(Cmax)/Rac(AUC) = ${rr.toFixed(2)} outside [${divLo}, ${divHi}] ` +
            `— possible absorption/Tmax shift at steady state (screening heuristic).`});
      }
    }
    divergenceOut[r] = diverges;
  }

  const cols: DG.Column[] = [];
  for (let j = 0; j < stratifyCols.length; j++)
    cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, stratifyCols[j], stratOut[j]));
  cols.push(f64Column('DOSE', doseOut));
  const racAucCol = f64Column('RAC_AUC', racAucOut);
  racAucCol.setTag('label', `Rac(AUC) Day${dayN}/Day${day1}`);
  cols.push(racAucCol);
  const racCmaxColOut = f64Column('RAC_CMAX', racCmaxOut);
  racCmaxColOut.setTag('label', `Rac(Cmax) Day${dayN}/Day${day1}`);
  cols.push(racCmaxColOut);
  cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'RAC_FORMULA', formulaOut));
  const racTheoColOut = f64Column('RAC_THEO', racTheoOut);
  racTheoColOut.setTag('label', 'Accumulation Index (theoretical)');
  cols.push(racTheoColOut);
  cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'RAC_SS_FLAG', ssFlagOut));
  cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.BOOL, 'RAC_DIVERGENCE', divergenceOut));

  const table = DGNS.DataFrame.fromColumns(cols);
  table.name = 'Day N / Day 1';
  return {table, day1, dayN, warnings};
}

// =====================================================================
// AD-3 — NOAEL exposure / safety margin (per-sex, worst-case flagged)
// =====================================================================

export interface ExposureMarginOptions {
  /** NOAEL dose (mg/kg). REQUIRED — never auto-derived from PK. Absent ⇒
   *  NO_NOAEL_INPUT and no margin computed. */
  readonly noaelMgKg?: number;
  /** Clinical (human) dose, mg. Required for HED and projected margins. */
  readonly humanDoseMg: number;
  /** Human body weight, kg (clinical dose-per-kg denominator). */
  readonly humanBodyWeightKg: number;
  /** Empirical human AUC (same AUCτ units as the rat). When provided, the
   *  exposure-based default is computed; absent ⇒ MARGIN_EXPOSURE = NaN +
   *  NO_HUMAN_AUC_INPUT. */
  readonly humanAuc?: number;
  /** Steady-state day for the rat NOAEL AUCτ. Default = max observed StudyDay. */
  readonly dayN?: number;
  /** Reference dose (mg/kg) for the projected-margin rat CL. Default
   *  {@link DEFAULT_CL_REF_DOSE_MG_KG}. */
  readonly clRefDoseMgKg?: number;
  /** Day for the rat CL reference. Default = min observed StudyDay (Day 1). */
  readonly clRefDay?: number;
  /** Reference rat weight (kg) for allometric scaling when the CL-reference rows carry no body
   *  weight. Default {@link DEFAULT_RAT_BW_KG}. When they do, the per-sex cohort mean of
   *  {@link bodyWeightCol} over those rows is the anchor (nca-cl-bw-normalised-units D-m). */
  readonly ratBodyWeightKg?: number;
  /** Body-weight column (kg) of the results frame. Default `'BW_KG'`. */
  readonly bodyWeightCol?: string;
  /** Disable allometric CL scaling (labels the projection "non-allometric" and
   *  raises MARGIN_NON_ALLOMETRIC). Default false (allometric is the default). */
  readonly disableAllometry?: boolean;
  /** Km factors for HED. Default 6 (rat) / 37 (human), FDA 2005. */
  readonly ratKm?: number;
  readonly humanKm?: number;
  readonly doseCol: string;
  readonly studyDayCol?: string;
  readonly stratifyCols?: ReadonlyArray<string>;
  readonly aucTauCol?: string;
  readonly clCol?: string;
}

export interface ExposureMarginResult {
  /** Per-sex table: stratify cols, `RAT_NOAEL_AUC`, `MARGIN_EXPOSURE`,
   *  `MARGIN_HED`, `MARGIN_PROJECTED`, `PROJECTION_ALLOMETRIC` (bool),
   *  `WORST_CASE` (bool — the smaller finite MARGIN_PROJECTED's sex), and the
   *  allometric anchor `RAT_REF_BW_KG` + `RAT_REF_BW_SOURCE` (`cohort mean of N profiles` /
   *  `default`; both blank when no projection is formed). */
  readonly table: DG.DataFrame | null;
  readonly warnings: TkWarning[];
}

/**
 * Compute the per-sex NOAEL exposure/safety margin trichotomy (exposure-based
 * default + HED + allometric CL projection) from a per-subject NCA results frame.
 *
 * @throws if a required column (dose / studyDay / aucτ / any stratify key) is
 *   absent from `results`.
 */
export function computeExposureMargin(
  results: DG.DataFrame,
  options: ExposureMarginOptions,
): ExposureMarginResult {
  const warnings: TkWarning[] = [];
  if (options.noaelMgKg === undefined || !Number.isFinite(options.noaelMgKg)) {
    warnings.push({code: 'NO_NOAEL_INPUT', cell: '',
      message: 'No NOAEL input (mg/kg); exposure/safety margin not computed. ' +
        'NOAEL is a user input, never auto-derived from PK.'});
    return {table: null, warnings};
  }
  const noael = options.noaelMgKg;
  const doseName = options.doseCol;
  const dayName = options.studyDayCol ?? 'StudyDay';
  const stratifyCols = options.stratifyCols ?? ['SEX'];
  const aucTauName = options.aucTauCol ?? 'AUC_LAST';
  const clName = options.clCol ?? 'CL';
  const ratKm = options.ratKm ?? DEFAULT_RAT_KM;
  const humanKm = options.humanKm ?? DEFAULT_HUMAN_KM;
  const clRefDose = options.clRefDoseMgKg ?? DEFAULT_CL_REF_DOSE_MG_KG;
  const defaultRatBw = options.ratBodyWeightKg ?? DEFAULT_RAT_BW_KG;
  const allometric = !options.disableAllometry;

  const required = [doseName, dayName, aucTauName, ...stratifyCols];
  for (const name of required) {
    if (results.col(name) === null)
      throw new Error(`computeExposureMargin: required column '${name}' not found in results.`);
  }
  const doseCol = results.col(doseName)!;
  const dayCol = results.col(dayName)!;
  const aucTauCol = results.col(aucTauName)!;
  // The margin's CL input contract (nca-cl-bw-normalised-units D-m): the ratified formula is in
  // L/h/kg, so the column must SAY it is per kg. `CL` is per kg in mL/min/kg on every surface; an
  // untagged or absolute column is refused rather than read as per kg (E6: an absolute CL read
  // as per kg understated the margin ~4.5× on the NONMEM 06 path).
  const clColRaw = results.col(clName);
  const clUnit = clColRaw?.getTag('units');
  const clToLhKg = clUnit === CL_PER_KG_UNIT ? 1 / L_PER_H_TO_ML_PER_MIN : (clUnit === 'L/h/kg' ? 1 : null);
  if (clColRaw !== null && clToLhKg === null) {
    warnings.push({code: 'MARGIN_CL_NOT_PER_KG', cell: '',
      message: `Column '${clName}' is not a per-kg clearance (units: ${clUnit == null || clUnit === '' ?
        'none' : clUnit}); MARGIN_PROJECTED needs CL in mL/min/kg or L/h/kg and is not computed.`});
  }
  const clCol = clToLhKg === null ? null : clColRaw;
  const bwCol = results.col(options.bodyWeightCol ?? 'BW_KG');
  const stratCols = stratifyCols.map((s) => results.col(s)!);
  const n = results.rowCount;

  let minDay = Infinity;
  let maxDay = -Infinity;
  for (let i = 0; i < n; i++) {
    const d = numAt(dayCol, i);
    if (!Number.isFinite(d)) continue;
    if (d < minDay) minDay = d;
    if (d > maxDay) maxDay = d;
  }
  const dayN = options.dayN ?? maxDay;
  const clRefDay = options.clRefDay ?? minDay;

  const humanDosePerKg = options.humanDoseMg / options.humanBodyWeightKg;
  const hedMgKg = noael * (ratKm / humanKm);
  const marginHed = humanDosePerKg > 0 ? hedMgKg / humanDosePerKg : NaN; // sex-identical
  if (!allometric) {
    warnings.push({code: 'MARGIN_NON_ALLOMETRIC', cell: '',
      message: 'Allometric CL scaling disabled — MARGIN_PROJECTED uses rat CL/kg directly ' +
        `(≈${Math.pow(70 / 0.25, 0.25).toFixed(1)}× margin inflation vs the allometric default).`});
  }

  // Build the ordered list of stratum values present (sex worst-case ranking).
  const strata = new Map<string, string[]>();
  for (let i = 0; i < n; i++) {
    const vals = stratifyCols.map((s, j) => stratCols[j].isNone(i) ? '' : stratCols[j].getString(i));
    strata.set(vals.join('|'), vals);
  }
  const stratumList = [...strata.values()].sort((a, b) => a.join('|').localeCompare(b.join('|')));

  const m = stratumList.length;
  const stratOut = stratifyCols.map(() => new Array<string>(m));
  const ratNoaelAucOut = new Float64Array(m);
  const marginExposureOut = new Float64Array(m);
  const marginHedOut = new Float64Array(m);
  const marginProjectedOut = new Float64Array(m);
  const allometricOut = new Array<boolean>(m);
  const ratRefBwOut = new Float64Array(m);
  const ratRefBwSourceOut = new Array<string>(m);

  let humanAucMissingWarned = false;
  for (let r = 0; r < m; r++) {
    const vals = stratumList[r];
    for (let j = 0; j < stratifyCols.length; j++) stratOut[j][r] = vals[j];
    const cellKey = stratifyCols.map((s, j) => `${s}=${vals[j]}`).join('|');
    const matches = (i: number, dose: number, day: number): boolean => {
      if (numAt(doseCol, i) !== dose || numAt(dayCol, i) !== day) return false;
      for (let j = 0; j < stratifyCols.length; j++) {
        const v = stratCols[j].isNone(i) ? '' : stratCols[j].getString(i);
        if (v !== vals[j]) return false;
      }
      return true;
    };

    // rat NOAEL AUCτ at Day N for this sex.
    const noaelAucVals: number[] = [];
    for (let i = 0; i < n; i++) if (matches(i, noael, dayN)) noaelAucVals.push(numAt(aucTauCol, i));
    const ratNoaelAuc = mean(noaelAucVals);
    ratNoaelAucOut[r] = ratNoaelAuc;

    // (1) Exposure-based default (preferred, ICH M3(R2)).
    if (options.humanAuc !== undefined && Number.isFinite(options.humanAuc) && options.humanAuc > 0)
      marginExposureOut[r] = ratNoaelAuc / options.humanAuc;
    else {
      marginExposureOut[r] = NaN;
      if (!humanAucMissingWarned) {
        humanAucMissingWarned = true;
        warnings.push({code: 'NO_HUMAN_AUC_INPUT', cell: '',
          message: 'No empirical human AUC provided — the exposure-based margin (the ' +
            'ratified default) is NaN. HED + projection are shown instead (DG-1).'});
      }
    }

    // (2) HED dose-based (sex-identical).
    marginHedOut[r] = marginHed;

    // (3) CL-projected (allometric default), per sex. The rat CL is per kg (L/h/kg); a reference
    // row with no per-kg CL is left out and the count is stated (`n used of N`). The allometric
    // anchor is the cohort's own mean weight over the SAME rows when they carry one, else the
    // labelled default. The estimator stays the ratified one (mean rat CL × one factor); against a
    // mean of per-animal projections it differs only through the BW spread inside (BW/W_h)^0.25,
    // under 1 % for a ±20 % spread (0.8^0.25 = 0.946, 1.2^0.25 = 1.047, mean 0.997).
    const clRefVals: number[] = [];
    const bwRefVals: number[] = [];
    let nRef = 0;
    if (clCol !== null) {
      for (let i = 0; i < n; i++) {
        if (!matches(i, clRefDose, clRefDay)) continue;
        nRef++;
        const v = numAt(clCol, i) * clToLhKg!;
        if (!Number.isFinite(v)) continue;
        clRefVals.push(v);
        if (bwCol !== null) bwRefVals.push(numAt(bwCol, i));
      }
    }
    if (nRef > clRefVals.length && clRefVals.length > 0) {
      warnings.push({code: 'MARGIN_CL_ROWS_DROPPED', cell: cellKey,
        message: `${cellKey}: MARGIN_PROJECTED uses the rat CL of ${clRefVals.length} of ${nRef} ` +
          `Day ${clRefDay} ${clRefDose} mg/kg profiles (the others have no per-kg CL).`});
    }
    const cohortBw = mean(bwRefVals);
    const nBw = bwRefVals.filter((v) => Number.isFinite(v)).length;
    const ratBw = Number.isFinite(cohortBw) && cohortBw > 0 ? cohortBw : defaultRatBw;
    const allomFactor = allometric ? Math.pow(options.humanBodyWeightKg / ratBw, ALLOMETRIC_EXPONENT) : 1;
    const ratCl = mean(clRefVals);
    const projected = Number.isFinite(ratCl) && ratCl > 0;
    // The anchor is stated only for a projection that was formed: with no usable CL there is no
    // projection, and a source would claim a number that does not exist (thesis #4). `n` counts the
    // CL-reference profiles (a pooled NONMEM cell is one profile, not one animal).
    ratRefBwOut[r] = projected ? ratBw : NaN;
    ratRefBwSourceOut[r] = !projected ? '' : Number.isFinite(cohortBw) && cohortBw > 0 ?
      `cohort mean of ${nBw} ${nBw === 1 ? 'profile' : 'profiles'}` : 'default';
    if (projected) {
      const clHuman = ratCl * allomFactor;
      const projectedHumanAuc = humanDosePerKg / clHuman;
      marginProjectedOut[r] = Number.isFinite(projectedHumanAuc) && projectedHumanAuc > 0 ?
        ratNoaelAuc / projectedHumanAuc : NaN;
    } else {
      marginProjectedOut[r] = NaN;
      warnings.push({code: 'NO_CL_REFERENCE_FOR_SEX', cell: cellKey,
        message: `No valid Day ${clRefDay} CL at ${clRefDose} mg/kg for ${cellKey}` +
          (nRef > 0 ? ` (0 of ${nRef} profiles have a per-kg CL)` : '') +
          `; MARGIN_PROJECTED is NaN (present ≠ valid — unreliable λz-derived CL excluded).`});
    }
    allometricOut[r] = allometric;
  }

  // Worst-case = the sex with the smallest finite MARGIN_PROJECTED.
  const worstOut = new Array<boolean>(m).fill(false);
  let worstIdx = -1;
  let worstVal = Infinity;
  for (let r = 0; r < m; r++) {
    const v = marginProjectedOut[r];
    if (Number.isFinite(v) && v < worstVal) {worstVal = v; worstIdx = r;}
  }
  if (worstIdx >= 0) worstOut[worstIdx] = true;

  const cols: DG.Column[] = [];
  for (let j = 0; j < stratifyCols.length; j++)
    cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, stratifyCols[j], stratOut[j]));
  cols.push(f64Column('RAT_NOAEL_AUC', ratNoaelAucOut));
  const expCol = f64Column('MARGIN_EXPOSURE', marginExposureOut);
  expCol.setTag('label', 'Exposure margin (AUCτ-based)');
  cols.push(expCol);
  const hedColOut = f64Column('MARGIN_HED', marginHedOut);
  hedColOut.setTag('label', 'HED margin (BSA/Km)');
  cols.push(hedColOut);
  const projCol = f64Column('MARGIN_PROJECTED', marginProjectedOut);
  projCol.setTag('label', allometric ? 'Projected margin (allometric)' : 'Projected margin (non-allometric)');
  cols.push(projCol);
  cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.BOOL, 'PROJECTION_ALLOMETRIC', allometricOut));
  cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.BOOL, 'WORST_CASE', worstOut));
  cols.push(f64Column('RAT_REF_BW_KG', ratRefBwOut).setTag('units', BW_UNIT));
  cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'RAT_REF_BW_SOURCE', ratRefBwSourceOut));

  const table = DGNS.DataFrame.fromColumns(cols);
  table.name = 'ExposureMargin';
  return {table, warnings};
}

// =====================================================================
// AD-4 — Sex fold-difference (over the aggregateNcaResults output)
// =====================================================================

export interface SexFoldOptions {
  /** Sex column in the aggregated frame. Default `'SEX'`. */
  readonly sexCol?: string;
  /** Parameter base names to fold (the `MEAN_<P>` columns). Default: every
   *  `MEAN_<P>` present in the aggregated frame. */
  readonly params?: ReadonlyArray<string>;
}

export interface SexFoldResult {
  /** Per-cell (non-sex group combination) table: the non-sex group columns, then
   *  per param `F_MEAN_<P>` / `M_MEAN_<P>` (+ N and CV% when present),
   *  `FOLD_<P>` (= max/min ≥ 1) and `FOLD_DIR_<P>` ('F>M' | 'M>F' | 'M=F').
   *  No significance column (composite N — thesis #4 / G9). */
  readonly table: DG.DataFrame;
}

/**
 * Derive sex fold-differences from an {@link aggregateNcaResults} output that was
 * grouped by (… , Sex). `FOLD = max(M,F)/min(M,F)` (≥ 1) plus a direction label;
 * carries the N and %CV the fold rests on. No p-value / significance claim.
 *
 * @throws if the sex column is absent from `aggregated`.
 */
export function computeSexFold(
  aggregated: DG.DataFrame,
  options: SexFoldOptions = {},
): SexFoldResult {
  const sexName = options.sexCol ?? 'SEX';
  const sexCol = aggregated.col(sexName);
  if (sexCol === null)
    throw new Error(`computeSexFold: sex column '${sexName}' not found in aggregated frame.`);

  // Discover params from MEAN_<P> columns unless given.
  const meanPrefix = 'MEAN_';
  const allParams: string[] = [];
  for (const c of aggregated.columns as Iterable<DG.Column>)
    if (c.name.startsWith(meanPrefix)) allParams.push(c.name.slice(meanPrefix.length));
  const params = (options.params ?? allParams).filter((p) => aggregated.col(`${meanPrefix}${p}`) !== null);

  // Non-sex group columns = string/int columns that are neither sex nor a
  // derived stat (MEAN_/SD_/CV_PCT_/N_/N). These define the fold cells.
  const isStat = (name: string): boolean =>
    name.startsWith('MEAN_') || name.startsWith('SD_') || name.startsWith('CV_PCT_') ||
    name.startsWith('N_') || name === 'N';
  const groupCols: DG.Column[] = [];
  for (const c of aggregated.columns as Iterable<DG.Column>)
    if (c.name !== sexName && !isStat(c.name)) groupCols.push(c);

  const nRows = aggregated.rowCount;
  // Group rows by the non-sex key; collect per-sex row index.
  interface Cell {groupVals: string[]; bySex: Map<string, number>;}
  const cells = new Map<string, Cell>();
  for (let i = 0; i < nRows; i++) {
    const groupVals = groupCols.map((c) => c.isNone(i) ? '' : c.getString(i));
    const key = groupVals.join('|');
    let cell = cells.get(key);
    if (cell === undefined) {cell = {groupVals, bySex: new Map()}; cells.set(key, cell);}
    cell.bySex.set(sexCol.isNone(i) ? '' : sexCol.getString(i), i);
  }
  const cellList = [...cells.values()].sort((a, b) => a.groupVals.join('|').localeCompare(b.groupVals.join('|')));
  const m = cellList.length;

  const groupOut = groupCols.map(() => new Array<string>(m));
  // Per param: F_MEAN, M_MEAN, (F_N, M_N), (F_CV, M_CV), FOLD, FOLD_DIR.
  const fMean = new Map<string, Float64Array>();
  const mMean = new Map<string, Float64Array>();
  const fN = new Map<string, Float64Array>();
  const mN = new Map<string, Float64Array>();
  const fCv = new Map<string, Float64Array>();
  const mCv = new Map<string, Float64Array>();
  const fold = new Map<string, Float64Array>();
  const foldDir = new Map<string, string[]>();
  const hasN = new Map<string, boolean>();
  const hasCv = new Map<string, boolean>();
  for (const p of params) {
    fMean.set(p, new Float64Array(m)); mMean.set(p, new Float64Array(m));
    fN.set(p, new Float64Array(m)); mN.set(p, new Float64Array(m));
    fCv.set(p, new Float64Array(m)); mCv.set(p, new Float64Array(m));
    fold.set(p, new Float64Array(m)); foldDir.set(p, new Array<string>(m));
    hasN.set(p, aggregated.col(`N_${p}`) !== null);
    hasCv.set(p, aggregated.col(`CV_PCT_${p}`) !== null);
  }

  for (let r = 0; r < m; r++) {
    const cell = cellList[r];
    for (let j = 0; j < groupCols.length; j++) groupOut[j][r] = cell.groupVals[j];
    const fRow = cell.bySex.get('F');
    const mRow = cell.bySex.get('M');
    for (const p of params) {
      const meanCol = aggregated.col(`${meanPrefix}${p}`)!;
      const fVal = fRow !== undefined ? numAt(meanCol, fRow) : NaN;
      const mVal = mRow !== undefined ? numAt(meanCol, mRow) : NaN;
      fMean.get(p)![r] = fVal;
      mMean.get(p)![r] = mVal;
      if (hasN.get(p)) {
        const nCol = aggregated.col(`N_${p}`)!;
        fN.get(p)![r] = fRow !== undefined ? numAt(nCol, fRow) : NaN;
        mN.get(p)![r] = mRow !== undefined ? numAt(nCol, mRow) : NaN;
      }
      if (hasCv.get(p)) {
        const cvCol = aggregated.col(`CV_PCT_${p}`)!;
        fCv.get(p)![r] = fRow !== undefined ? numAt(cvCol, fRow) : NaN;
        mCv.get(p)![r] = mRow !== undefined ? numAt(cvCol, mRow) : NaN;
      }
      if (Number.isFinite(fVal) && Number.isFinite(mVal) && fVal > 0 && mVal > 0) {
        const hi = Math.max(fVal, mVal);
        const lo = Math.min(fVal, mVal);
        fold.get(p)![r] = hi / lo;
        foldDir.get(p)![r] = fVal === mVal ? 'M=F' : (fVal > mVal ? 'F>M' : 'M>F');
      } else {
        fold.get(p)![r] = NaN;
        foldDir.get(p)![r] = '';
      }
    }
  }

  const cols: DG.Column[] = [];
  for (let j = 0; j < groupCols.length; j++)
    cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, groupCols[j].name, groupOut[j]));
  for (const p of params) {
    cols.push(f64Column(`F_MEAN_${p}`, fMean.get(p)!));
    cols.push(f64Column(`M_MEAN_${p}`, mMean.get(p)!));
    if (hasN.get(p)) {
      cols.push(f64Column(`F_N_${p}`, fN.get(p)!));
      cols.push(f64Column(`M_N_${p}`, mN.get(p)!));
    }
    if (hasCv.get(p)) {
      cols.push(f64Column(`F_CV_PCT_${p}`, fCv.get(p)!));
      cols.push(f64Column(`M_CV_PCT_${p}`, mCv.get(p)!));
    }
    const foldColOut = f64Column(`FOLD_${p}`, fold.get(p)!);
    foldColOut.setTag('label', `Fold ${p} (max/min)`);
    cols.push(foldColOut);
    cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, `FOLD_DIR_${p}`, foldDir.get(p)!));
  }

  const table = DGNS.DataFrame.fromColumns(cols);
  table.name = 'SexFold';
  return {table};
}
