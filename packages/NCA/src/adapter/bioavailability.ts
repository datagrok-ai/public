/**
 * Bioavailability (F%) — cross-stratum reference-ratio primitive.
 *
 * **Why this is a NEW module, not `aggregateNcaResults`.** `aggregateNcaResults`
 * does *within-group* arithmetic (mean / SD / %CV / N) via the native
 * `groupBy().aggregate()`. F% is intrinsically a *cross-group* ratio: each
 * extravascular stratum's dose-normalized AUCinf divided by the **IV arm's**
 * dose-normalized AUCinf of the **same sex** (+ optional covariates). The native
 * aggregator cannot divide stratum A by stratum B, so this is a separate
 * primitive; `aggregateNcaResults` is untouched (no regression — AD-6 / AC-D4).
 *
 * **Correct layering (not vendored math).** The per-profile primitives this
 * consumes (AUCinf, dose-normalization) are sci-comp's. F% itself is a
 * cross-profile aggregate over DataFrame strata — DataFrame orchestration, which
 * CLAUDE.md places in this package. No core math is re-derived.
 *
 * **Scientific contract (AD-6):**
 *  - Denominator is **sex-stratified** (+ optional covariate keys). Pooling
 *    across sex reintroduces a real ~36% error (F iv_dn_AUC = 1.874 vs M = 2.549
 *    in the rat reference). Mandatory.
 *  - **Uncapped** — F% > 100 is real (SC/F/25 = 117.5%); never clamped.
 *  - A stratum with **no IV arm** → `F_pct = NaN` + `NO_IV_REFERENCE_FOR_STRATUM`
 *    warning; never a pooled / cross-stratum fallback.
 *  - Reference statistic is `denominatorStat` (`'mean'` arithmetic by default —
 *    matches the validated fixture; `'gmean'` geometric is the log-normal/PKNCA
 *    alternative, R1-SF-5). Surfaced so a ratification flip needs no re-impl.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import {f64Column, numAt} from './adapter-utils';

/** Central-tendency statistic for the IV reference denominator (R1-SF-5). */
export type DenominatorStat = 'mean' | 'gmean';

export interface BioavailabilityOptions {
  /** AUCinf column in the results frame. Default `'AUC_INF'`. */
  readonly aucInfCol?: string;
  /** Per-subject dose column (carried metadata, e.g. `'Dose_mg_kg'`). Required. */
  readonly doseCol: string;
  /** Route column (writer output strings). Default `'ROUTE'`. Intravascular =
   *  values starting with `'IV'` (IV-bolus / IV-infusion). */
  readonly routeCol?: string;
  /** Stratification keys for matching the IV reference — sex is mandatory;
   *  append study covariates as needed. Default `['SEX']`. */
  readonly stratifyCols?: ReadonlyArray<string>;
  /** Central-tendency statistic for the IV denominator. Default `'mean'`. */
  readonly denominatorStat?: DenominatorStat;
  /**
   * The occasion column (tk-study-day-stratification F2), when it is one of `stratifyCols`
   * (same-occasion pairing of a repeat-dose design). An extravascular cell whose IV reference
   * exists only on another occasion takes that occasion's IV with a
   * `CROSS_OCCASION_IV_REFERENCE` warning naming both occasions (P6), and the table gains an
   * `IV_OCCASION` column.
   */
  readonly occasionCol?: string;
}

/** One stratum-level warning surfaced alongside the F% table. */
export interface BioavailabilityWarning {
  readonly code: 'NO_IV_REFERENCE_FOR_STRATUM' | 'CROSS_OCCASION_IV_REFERENCE';
  readonly message: string;
  /** The stratum key that lacked an IV reference (e.g. `'SEX=M'`). */
  readonly stratum: string;
}

export interface BioavailabilityResult {
  /**
   * Stratum-level F% table: `ROUTE`, the stratify columns, `DOSE`, `DN_AUC`
   * (stratum central-tendency dose-normalized AUCinf), `IV_DN_AUC` (the matched
   * IV reference), `F_PCT`. One row per non-IV `{route, …stratify, dose}` cell.
   */
  readonly table: DG.DataFrame;
  /** Active denominator statistic (surfaced in the grid header). */
  readonly denominatorStat: DenominatorStat;
  /** Strata that had no IV reference (their `F_PCT` is NaN). */
  readonly warnings: BioavailabilityWarning[];
}

/** True when a route string denotes an intravascular case (IV-bolus / IV-infusion)
 *  — the IV/EV split %F is defined on. This is the SINGLE owner of the IV/EV
 *  classification for bioavailability: a row is "EV" (an F% is defined for it)
 *  exactly when this returns false. Exported so the compound-tier %F join
 *  (nca-analysis-open.ts) classifies rows by the identical rule (no parallel
 *  copy that could drift if the route convention changes). */
export function isIntravascular(routeStr: string): boolean {
  return routeStr.toUpperCase().startsWith('IV');
}

function centralTendency(values: ReadonlyArray<number>, stat: DenominatorStat): number {
  if (values.length === 0) return NaN;
  if (stat === 'gmean') {
    let sumLn = 0;
    let k = 0;
    for (const v of values)
      if (v > 0) {sumLn += Math.log(v); k++;}

    return k === 0 ? NaN : Math.exp(sumLn / k);
  }
  let sum = 0;
  for (const v of values) sum += v;
  return sum / values.length;
}

/**
 * Compute the stratum-level bioavailability (F%) table from a per-subject NCA
 * results frame that carries route + dose + the stratification metadata
 * (typically the post-`appendGroupingColumns` frame).
 *
 * @throws if a required column (AUCinf / dose / route / any stratify key) is
 *   absent from `results`.
 */
export function computeBioavailability(
  results: DG.DataFrame,
  options: BioavailabilityOptions,
): BioavailabilityResult {
  const aucName = options.aucInfCol ?? 'AUC_INF';
  const routeName = options.routeCol ?? 'ROUTE';
  const doseName = options.doseCol;
  const stratifyCols = options.stratifyCols ?? ['SEX'];
  const stat = options.denominatorStat ?? 'mean';

  const required = [aucName, routeName, doseName, ...stratifyCols];
  for (const name of required) {
    if (results.col(name) === null)
      throw new Error(`computeBioavailability: required column '${name}' not found in results.`);
  }

  const aucCol = results.col(aucName)!;
  const routeCol = results.col(routeName)!;
  const doseCol = results.col(doseName)!;
  const stratCols = stratifyCols.map((s) => results.col(s)!);
  const n = results.rowCount;

  const stratumKey = (row: number): string =>
    stratifyCols.map((s, j) => `${s}=${stratCols[j].isNone(row) ? '' : stratCols[j].getString(row)}`).join('|');

  // Pass 1: IV reference — collect dn_AUCinf per stratum over intravascular rows.
  const ivByStratum = new Map<string, number[]>();
  // Pass 1b: non-IV cells keyed by {route|stratum|dose}.
  interface Cell {
    route: string; stratumVals: string[]; dose: number; dn: number[];
  }
  const cells = new Map<string, Cell>();

  for (let i = 0; i < n; i++) {
    // `numAt` is the isNone-first accessor (hoisted, adapter-utils.ts) — the
    // fourth copy of this idiom used to be spelled inline right here. Same
    // semantics: a platform-null cell arrives as NaN and is dropped by the
    // `isFinite` guard on the next line.
    const auc = numAt(aucCol, i);
    const dose = numAt(doseCol, i);
    if (!Number.isFinite(auc) || !Number.isFinite(dose) || dose <= 0) continue;
    const dn = auc / dose;
    if (!Number.isFinite(dn)) continue;

    const routeStr = routeCol.isNone(i) ? '' : routeCol.getString(i);
    const sKey = stratumKey(i);
    if (isIntravascular(routeStr)) {
      let arr = ivByStratum.get(sKey);
      if (arr === undefined) {arr = []; ivByStratum.set(sKey, arr);}
      arr.push(dn);
      continue;
    }
    // Non-IV cell.
    const cellKey = `${routeStr}|${sKey}|${dose}`;
    let cell = cells.get(cellKey);
    if (cell === undefined) {
      cell = {
        route: routeStr,
        stratumVals: stratifyCols.map((s, j) => stratCols[j].isNone(i) ? '' : stratCols[j].getString(i)),
        dose,
        dn: [],
      };
      cells.set(cellKey, cell);
    }
    cell.dn.push(dn);
  }

  const ivRefByStratum = new Map<string, number>();
  for (const [sKey, arr] of ivByStratum)
    ivRefByStratum.set(sKey, centralTendency(arr, stat));

  // P6: the occasion-free stratum key → the IV strata on other occasions, in first-seen order.
  const occIdx = options.occasionCol !== undefined ? stratifyCols.indexOf(options.occasionCol) : -1;
  const keyOf = (vals: ReadonlyArray<string>, dropOcc: boolean): string =>
    stratifyCols.map((s, j) => dropOcc && j === occIdx ? `${s}=*` : `${s}=${vals[j]}`).join('|');
  const ivByOccFree = new Map<string, {occ: string; key: string}[]>();
  if (occIdx >= 0) {
    for (const sKey of ivByStratum.keys()) {
      const vals = sKey.split('|').map((p) => p.substring(p.indexOf('=') + 1));
      const free = keyOf(vals, true);
      ivByOccFree.set(free, [...(ivByOccFree.get(free) ?? []), {occ: vals[occIdx], key: sKey}]);
    }
  }
  const ivOccOut: string[] = [];

  // Build output rows (deterministic order: route, stratum, dose).
  const cellList = [...cells.values()].sort((a, b) =>
    a.route.localeCompare(b.route) ||
    a.stratumVals.join('|').localeCompare(b.stratumVals.join('|')) ||
    a.dose - b.dose);

  const warnings: BioavailabilityWarning[] = [];
  const warnedStrata = new Set<string>();

  const routeOut = new Array<string>(cellList.length);
  const stratOut = stratifyCols.map(() => new Array<string>(cellList.length));
  const doseOut = new Float64Array(cellList.length);
  const dnOut = new Float64Array(cellList.length);
  const ivDnOut = new Float64Array(cellList.length);
  // F-S5.2 (D11 "N beside each parameter"): how many subjects feed each side of the ratio — an
  // Exclude λz decision drops a subject's AUC∞, and the N makes that visible instead of silently
  // moving %F.
  const nOut = new Int32Array(cellList.length);
  const ivNOut = new Int32Array(cellList.length);
  const fpctOut = new Float64Array(cellList.length);

  for (let r = 0; r < cellList.length; r++) {
    const cell = cellList[r];
    let sKey = stratifyCols.map((s, j) => `${s}=${cell.stratumVals[j]}`).join('|');
    const dn = centralTendency(cell.dn, stat);
    let ivRef = ivRefByStratum.get(sKey);
    if (occIdx >= 0) {
      let ivOcc = cell.stratumVals[occIdx];
      if (ivRef === undefined) {
        const other = ivByOccFree.get(keyOf(cell.stratumVals, true))?.[0];
        if (other !== undefined) {
          const ownKey = sKey;
          sKey = other.key;
          ivRef = ivRefByStratum.get(sKey);
          ivOcc = other.occ;
          if (!warnedStrata.has(`${ownKey}>${sKey}`)) {
            warnedStrata.add(`${ownKey}>${sKey}`);
            warnings.push({
              code: 'CROSS_OCCASION_IV_REFERENCE',
              message: `No IV arm on ${options.occasionCol} ${cell.stratumVals[occIdx]} for ${ownKey}; ` +
                `%F uses the IV reference of ${options.occasionCol} ${other.occ} — exposure may differ ` +
                'between occasions (accumulation, time-dependent clearance).',
              stratum: ownKey,
            });
          }
        }
      }
      ivOccOut.push(ivRef !== undefined ? ivOcc : '');
    }

    routeOut[r] = cell.route;
    for (let j = 0; j < stratifyCols.length; j++) stratOut[j][r] = cell.stratumVals[j];
    doseOut[r] = cell.dose;
    dnOut[r] = dn;
    nOut[r] = cell.dn.length;
    ivNOut[r] = ivByStratum.get(sKey)?.length ?? 0;

    if (ivRef === undefined || !Number.isFinite(ivRef)) {
      ivDnOut[r] = NaN;
      fpctOut[r] = NaN; // uncapped; NaN = no IV reference
      if (!warnedStrata.has(sKey)) {
        warnedStrata.add(sKey);
        warnings.push({
          code: 'NO_IV_REFERENCE_FOR_STRATUM',
          message:
            `No IV arm for stratum ${sKey}; %F is undefined (NaN) for its ` +
            `extravascular subjects. Provide an IV reference for this stratum.`,
          stratum: sKey,
        });
      }
    } else {
      ivDnOut[r] = ivRef;
      fpctOut[r] = (dn / ivRef) * 100; // uncapped (F% > 100 is real)
    }
  }

  const cols: DG.Column[] = [DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'ROUTE', routeOut)];
  for (let j = 0; j < stratifyCols.length; j++)
    cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, stratifyCols[j], stratOut[j]));
  cols.push(f64Column('DOSE', doseOut));
  cols.push(f64Column('DN_AUC', dnOut));
  cols.push(DGNS.Column.fromInt32Array('N', nOut));
  cols.push(f64Column('IV_DN_AUC', ivDnOut));
  cols.push(DGNS.Column.fromInt32Array('IV_N', ivNOut));
  if (occIdx >= 0) cols.push(DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'IV_OCCASION', ivOccOut));
  const fpctCol = f64Column('F_PCT', fpctOut);
  fpctCol.setTag('label', stat === 'gmean' ? '%F (geom.)' : '%F (arith.)');
  cols.push(fpctCol);

  const table = DGNS.DataFrame.fromColumns(cols);
  table.name = 'Bioavailability';
  return {table, denominatorStat: stat, warnings};
}
