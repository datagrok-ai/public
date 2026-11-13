/**
 * Unit conversion at the DG.DataFrame ↔ nca.ProfileInputs boundary.
 *
 * **Canonical internal units (AD-2 — fixture-constrained by NFR-04).** Time
 * `h`, concentration `mg/L`, dose `mg`, body weight `kg`. The PKNCA reference
 * fixtures (theoph in h/mg/L; indometh in h/mcg/mL ≡ mg/L) pin these via
 * NFR-04 tolerance assertions; any alternative would force reference-suite
 * recalibration.
 *
 * **Locus (AD-1).** Rescaling happens once, at the adapter↔core boundary,
 * inside `scaleProfileInputs`. Sci-comp stays unit-agnostic. The pipeline-
 * bypass invariant (AD-12 / G9): all `nca.ProfileInputs` MUST be constructed
 * via `extractProfileInputs` → `scaleProfileInputs`. Direct construction
 * bypasses canonical rescaling and is forbidden. See `profile-inputs.ts`
 * JSDoc + architecture spec invariant #7.
 *
 * **Type safety (E5 / F1).** `scaleProfileInputs` parameter types come from
 * the sci-comp unit type literals (`nca.TimeUnit`, `nca.ConcentrationUnit`,
 * `nca.DoseUnit`); passing a non-canonical string is a compile-time error.
 *
 * **`BodyWeightUnit` placement (architect S-2).** Defined here (leaf module);
 * re-exported via `adapter/index.ts`; `profile-inputs.ts` imports the type
 * from `../adapter` (not from this module directly) to avoid a circular
 * import (this module is a sibling-leaf of `profile-inputs.ts`; the import
 * direction is extractor → scaler).
 *
 * **No side effects, no DG imports.** Testable in isolation; ~22 tests in
 * `tests/adapter-unit-conversion-tests.ts`.
 */

import {nca} from '@datagrok-libraries/sci-comp';

/** Body weight unit (AD-6). v1 supports `kg` and `g`; `lb` is forward-deferred
 *  per AD-6 known limitation (a 15-lb beagle reported as 15 with bwUnit='kg'
 *  is not caught by the species-absent heuristic at the kg-15 boundary). */
export type BodyWeightUnit = 'kg' | 'g';

/** Time → hours scale factor (multiply input by this). */
export const TIME_TO_HOUR: Readonly<Record<nca.TimeUnit, number>> = Object.freeze({
  'h': 1,
  'min': 1 / 60,
  's': 1 / 3600,
  'd': 24,
});

/** Concentration → mg/L scale factor (mass-based only — molar `µmol/L` /
 *  `nmol/L` are rejected by `UNITS_MISMATCH` validator pre-conversion). */
export const CONC_TO_MG_L: Readonly<Record<nca.ConcentrationUnit, number>> = Object.freeze({
  'mg/L': 1,
  'µg/mL': 1, // numerically identical
  'µg/L': 1 / 1000,
  'ng/mL': 1 / 1000,
  // Molar entries listed for type completeness; multiplication by these
  // numbers is never reached — `validateScalarsForMolar` hard-blocks earlier.
  // NaN keeps any erroneous bypass numerically obvious instead of silently
  // returning identity.
  'µmol/L': NaN,
  'nmol/L': NaN,
});

/** Dose → mg base scale factor. Per-kg variants (`mg/kg`, `µg/kg`) are
 *  combined with `bodyWeight_kg` inside `scaleProfileInputs`. */
export const DOSE_TO_MG: Readonly<Record<nca.DoseUnit, number>> = Object.freeze({
  'mg': 1,
  'µg': 1 / 1000,
  'ng': 1 / 1e6,
  'mg/kg': 1, // × bodyWeight_kg at call site
  'µg/kg': 1 / 1000, // × bodyWeight_kg at call site
});

/** A per-kg dose unit (`mg/kg`, `µg/kg`) — the only units the scaler multiplies by body weight. */
export const isPerKgDose = (u: nca.DoseUnit): boolean => u === 'mg/kg' || u === 'µg/kg';

/** Body weight → kg scale factor. */
export const BW_TO_KG: Readonly<Record<BodyWeightUnit, number>> = Object.freeze({
  'kg': 1,
  'g': 1 / 1000,
});

/** CL per kg, the reporting unit of the `CL` column on every surface (D-b). */
export const CL_PER_KG_UNIT = 'mL/min/kg';
/** Vz / Vss per kg, the reporting unit of the `VZ` / `VSS` columns (D-b). */
export const VOLUME_PER_KG_UNIT = 'L/kg';
/** Absolute CL, the unit of the `CL_ABS` companion (canonical engine unit). */
export const CL_ABS_UNIT = 'L/h';
/** Absolute volume, the unit of the `VZ_ABS` / `VSS_ABS` companions. */
export const VOLUME_ABS_UNIT = 'L';
/** Body weight, the unit of the `BW_KG` column (canonical). */
export const BW_UNIT = 'kg';
/** L/h/kg → mL/min/kg (× 1000 mL/L ÷ 60 min/h). The one scale factor of the per-kg CL. */
export const L_PER_H_TO_ML_PER_MIN = 1000 / 60;

/**
 * The package's reporting-basis contract for CL / VZ / VSS (nca-cl-bw-normalised-units
 * AD-8): `per-kg-1` = per-kg in mL/min/kg and L/kg, absolute companions `CL_ABS` / `VZ_ABS`
 * / `VSS_ABS`. Persisted with every freshly computed summary cache (EMS `results_basis`) and
 * carried in-band on the pipeline output (`RESULTS_BASIS`, `nca.resultsBasis`). A cache stamped
 * with another value (or none) is stale. Bump it when the meaning of a reported CL / volume
 * column changes again.
 */
export const RESULTS_BASIS = 'per-kg-1';

/**
 * Canonical `column.tags.units` strings, keyed by UPPERCASE column name.
 * Numeric columns get a tag; string / boolean / count / JSON columns are
 * intentionally absent from this map (their `column.tags.units` stays
 * `undefined`). `CL` / `VZ` / `VSS` mean per-kg on every row (the column identity,
 * nca-cl-bw-normalised-units AD-1), so their unit is uniform and lives here like
 * every other column's; the absolute values ride the `_ABS` companions.
 *
 * Format: Unicode-dot for v1 (`mg/L·h`, `mg/L·h²`, etc.). CDISC SEND PP
 * ASCII-asterisk translation (`mg/L·h` → `mg*h/L`) is deferred to Phase 3
 * export — a 1-line `.replace(...)` at SEND-export time; no storage change here.
 */
export const CANONICAL_UNIT_TAGS: Readonly<Record<string, string>> = Object.freeze({
  CMAX: 'mg/L',
  C0_EXTRAP: 'mg/L',
  TMAX: 'h',
  AUC_LAST: 'mg/L·h',
  AUC_INF: 'mg/L·h',
  PCT_EXTRAP: '%',
  LAMBDA_Z: '1/h',
  HALF_LIFE: 'h',
  CL: CL_PER_KG_UNIT,
  VZ: VOLUME_PER_KG_UNIT,
  VSS: VOLUME_PER_KG_UNIT,
  CL_ABS: CL_ABS_UNIT,
  VZ_ABS: VOLUME_ABS_UNIT,
  VSS_ABS: VOLUME_ABS_UNIT,
  BW_KG: BW_UNIT,
  // FR-200 moment / lag params (sci-comp 0.7.0). AUMC = first-moment area
  // (conc·time²).
  AUMC_LAST: 'mg/L·h²',
  AUMC_INF: 'mg/L·h²',
  MRT: 'h',
  TLAG: 'h',
  PCT_AUMC_EXTRAP: '%',
  LAMBDAZ_TFIRST: 'h',
  LAMBDAZ_TLAST: 'h',
  LAMBDAZ_INTERCEPT: 'mg/L',
});

/**
 * What the per-kg rule needs to know about one profile's body weight
 * (nca-cl-bw-normalised-units AD-1 / AD-2). Built by the extractor; `NO_BODY_WEIGHT`
 * when no weight is mapped.
 */
export interface BodyWeightContext {
  /** The engine's body weight (kg) — the profile's first row, as `scaleProfileInputs` used it. */
  readonly bodyWeightKg: number | null;
  /** Lowest and highest non-missing body weight (kg) over the profile's rows; `null` when none. */
  readonly bodyWeightRangeKg: {readonly min: number; readonly max: number} | null;
  /** The dose was given per kg (`mg/kg`, `µg/kg`): the per-kg value is exact whatever weight
   *  the engine multiplied by, so a varying weight does not refuse it. */
  readonly perKgDose: boolean;
}

/** No body weight mapped (a per-kg dose then reaches the engine per kg; an absolute dose has no per-kg value). */
export const NO_BODY_WEIGHT: BodyWeightContext = Object.freeze({
  bodyWeightKg: null, bodyWeightRangeKg: null, perKgDose: false,
});

/** The engine's three dose-normalised outputs (canonical: L/h, L, L — or L/h/kg, L/kg, L/kg on the per-kg path). */
export interface EngineClVolumes {
  readonly cl: number;
  readonly vz: number;
  readonly vss: number;
}

/** The seven reported values of {@link normaliseClVolumes}. NaN = cannot be formed. */
export interface NormalisedClVolumes {
  /** CL per kg, mL/min/kg. */
  readonly cl: number;
  /** Vz per kg, L/kg. */
  readonly vz: number;
  /** Vss per kg, L/kg. */
  readonly vss: number;
  /** Absolute CL, L/h (NaN on the per-kg path: no weight to multiply back). */
  readonly clAbs: number;
  /** Absolute Vz, L. */
  readonly vzAbs: number;
  /** Absolute Vss, L. */
  readonly vssAbs: number;
  /** The weight the per-kg value used (kg), NaN when there is none. */
  readonly bwKg: number;
  /** `PER_KG_NOT_FORMED` (per-kg refused) or `BW_VARIES_WITHIN_PROFILE` (absolute companions rest on the first row). */
  readonly warnings: ReadonlyArray<nca.ParameterWarning>;
}

const fmtKg = (v: number): string => Number(v.toPrecision(4)).toString();

/**
 * Form the reported CL / Vz / Vss of one profile from the engine values — the ONE place the
 * per-kg rule lives (nca-cl-bw-normalised-units AD-1, guardrail "one canonicalization function,
 * every caller"). The engine and its inputs are untouched: this is a division by the mapped
 * weight and a constant scale, like {@link BW_TO_KG}.
 *
 * | engine dose | weight | `CL`/`VZ`/`VSS` | `_ABS` | `BW_KG` |
 * |---|---|---|---|---|
 * | `per_kg` | none (every dense per-kg run) | engine value | NaN | NaN |
 * | `per_kg` | known (a sparse `/kg` composite's group mean) | engine value | per-kg × BW | BW |
 * | `absolute` | constant within the profile (or a `/kg` dose) | abs ÷ BW | engine value | BW |
 * | `absolute`, absolute-mg dose | varies within the profile | NaN + `PER_KG_NOT_FORMED` | engine | first-row BW |
 * | `absolute` | none | NaN + `PER_KG_NOT_FORMED` | engine value | NaN |
 *
 * A `/kg` dose with a varying weight keeps the exact per-kg value (the weight cancels) and warns
 * `BW_VARIES_WITHIN_PROFILE`: the absolute companions rest on the first row's weight (AD-2).
 * A weight ≤ 0 or non-finite is treated as absent. NaN engine values stay NaN.
 */
export function normaliseClVolumes(
  engine: EngineClVolumes, doseNorm: 'absolute' | 'per_kg', bw: BodyWeightContext,
): NormalisedClVolumes {
  const k = L_PER_H_TO_ML_PER_MIN;
  const w = bw.bodyWeightKg;
  const hasW = w !== null && Number.isFinite(w) && w > 0;
  if (doseNorm === 'per_kg') {
    // A dense per-kg run has no weight by construction (the scaler multiplies a `/kg` dose by any
    // weight it has). A sparse `/kg` composite runs per kg on purpose and may know its group-mean
    // weight (AD-5): then the absolute values are the per-kg ones × that weight.
    return {cl: engine.cl * k, vz: engine.vz, vss: engine.vss,
      clAbs: hasW ? engine.cl * w! : NaN, vzAbs: hasW ? engine.vz * w! : NaN,
      vssAbs: hasW ? engine.vss * w! : NaN, bwKg: hasW ? w! : NaN, warnings: []};
  }
  const abs = {clAbs: engine.cl, vzAbs: engine.vz, vssAbs: engine.vss};
  if (!hasW) {
    return {cl: NaN, vz: NaN, vss: NaN, ...abs, bwKg: NaN, warnings: [{
      code: 'PER_KG_NOT_FORMED', severity: 'info',
      message: 'No body weight: CL, Vz and Vss per kg are not available. The absolute values ' +
        '(CL abs., Vz abs., Vss abs.) are reported.',
    }]};
  }
  const r = bw.bodyWeightRangeKg;
  const varies = r !== null && r.min !== r.max;
  const range = varies ? `${fmtKg(r!.min)}–${fmtKg(r!.max)} kg` : '';
  if (varies && !bw.perKgDose) {
    return {cl: NaN, vz: NaN, vss: NaN, ...abs, bwKg: w!, warnings: [{
      code: 'PER_KG_NOT_FORMED', severity: 'warning',
      message: `Body weight varies within the profile (${range}); it is not known which weight to ` +
        'normalise by, so CL, Vz and Vss per kg are not reported. The absolute values are.',
    }]};
  }
  const warnings: nca.ParameterWarning[] = varies ? [{
    code: 'BW_VARIES_WITHIN_PROFILE', severity: 'info',
    message: `Body weight varies within the profile (${range}). CL, Vz and Vss per kg are exact ` +
      `(the dose is per kg); the absolute values use the first row's weight (${fmtKg(w!)} kg).`,
  }] : [];
  return {cl: engine.cl / w! * k, vz: engine.vz / w!, vss: engine.vss / w!, ...abs, bwKg: w!, warnings};
}

/** One pooled member's absolute dose and body weight (NaN = missing), in any consistent units. */
export interface CompositeMember {
  readonly dose: number;
  readonly bw: number;
}

/** The dose and the per-kg divisor of one composite (pooled) profile — {@link compositeDoseWeight}. */
export interface CompositeDoseWeight {
  /** Mean dose over the members with a dose — the engine's dose (unchanged by this rule). */
  readonly meanDose: number;
  /** Mean body weight over the members with BOTH a dose and a weight; null when none has both. */
  readonly meanBw: number | null;
  /** The weight the engine's absolute CL / volumes are divided by so that per kg =
   *  mean(dose) / mean(BW) / AUC over the SAME members (equals `meanBw` when every member has a
   *  weight); null when no member has both. */
  readonly perKgDivisor: number | null;
  readonly n: number;
  readonly nBw: number;
  /** `COMPOSITE_BW_PARTIAL` and / or `COMPOSITE_NOT_WEIGHT_BASED`. */
  readonly warnings: ReadonlyArray<nca.ParameterWarning>;
}

/** Decimals of a recorded value, read from its shortest representation at single precision (a
 *  float32-stored 0.2441 still reads 4, not 17). */
function recordedDecimals(v: number): number {
  const s = String(Number(Math.abs(v).toPrecision(7)));
  if (s.includes('e')) return 0;
  const i = s.indexOf('.');
  return i < 0 ? 0 : s.length - i - 1;
}

/**
 * The dose and per-kg divisor of a composite pooling animals of different weights
 * (nca-cl-bw-normalised-units AD-5 / AD-5b) — ONE helper for the sparse pooler and the NONMEM
 * cell pooler, so the two cannot drift. Per kg = mean(dose) / mean(BW) / AUC over the members
 * that have both (`COMPOSITE_BW_PARTIAL` when some lack a weight). That identity equals
 * dose-per-kg / AUC exactly when the dosing was weight-based, which is tested from the data, with
 * no threshold: every member's dose / BW must agree within the propagated rounding of the recorded
 * values, `½·10^−d(dose) / dose + ½·10^−d(BW) / BW` (`COMPOSITE_NOT_WEIGHT_BASED` otherwise).
 */
export function compositeDoseWeight(members: ReadonlyArray<CompositeMember>): CompositeDoseWeight {
  let doseSum = 0; let doseN = 0;
  const both: Array<{dose: number; bw: number}> = [];
  for (const m of members) {
    if (!Number.isFinite(m.dose)) continue;
    doseSum += m.dose; doseN++;
    if (Number.isFinite(m.bw) && m.bw > 0) both.push({dose: m.dose, bw: m.bw});
  }
  const meanDose = doseN > 0 ? doseSum / doseN : NaN;
  if (both.length === 0)
    return {meanDose, meanBw: null, perKgDivisor: null, n: doseN, nBw: 0, warnings: []};
  const meanBw = both.reduce((a, b) => a + b.bw, 0) / both.length;
  const meanDoseBoth = both.reduce((a, b) => a + b.dose, 0) / both.length;
  const perKgDivisor = meanDoseBoth > 0 ? meanBw * meanDose / meanDoseBoth : meanBw;
  const warnings: nca.ParameterWarning[] = [];
  if (both.length < doseN) {
    warnings.push({code: 'COMPOSITE_BW_PARTIAL', severity: 'warning',
      message: `Per-kg CL, Vz and Vss use the ${both.length} of ${doseN} animals with a body weight ` +
        `(mean ${fmtKg(meanBw)} kg).`});
  }
  // Weight-based check: dose / BW the same for every member within the recorded precision.
  // A zero-dose (control) composite has nothing to check.
  const q = both.filter((b) => b.dose > 0).map((b) => ({
    r: b.dose / b.bw,
    tol: 0.5 * 10 ** -recordedDecimals(b.dose) / b.dose + 0.5 * 10 ** -recordedDecimals(b.bw) / b.bw,
  }));
  let weightBased = true;
  for (let i = 0; i < q.length && weightBased; i++) {
    for (let j = i + 1; j < q.length; j++) {
      if (Math.abs(q[i].r - q[j].r) > (q[i].tol + q[j].tol) * Math.max(q[i].r, q[j].r)) {
        weightBased = false; break;
      }
    }
  }
  if (!weightBased) {
    warnings.push({code: 'COMPOSITE_NOT_WEIGHT_BASED', severity: 'warning',
      message: 'Doses are not proportional to body weight in this composite; per-kg CL = ' +
        'mean dose / mean body weight / AUC is an approximation.'});
  }
  return {meanDose, meanBw, perKgDivisor, n: doseN, nBw: both.length, warnings};
}

/**
 * Concentration-unit conversion factor. Single factor source — also consumed
 * by the LLOQ rescaling inside `scaleProfileInputs` (AD-7 same-factor
 * invariant) and by the LLOQ sanity validator's understanding of the
 * pre-rescaled scale.
 *
 * Returns `NaN` for molar units (`µmol/L`/`nmol/L`); `validateScalarsForMolar`
 * runs upstream and prevents these reaching the conversion site.
 */
export function convertConcentrationFactor(unit: nca.ConcentrationUnit): number {
  return CONC_TO_MG_L[unit];
}

/**
 * Output of {@link scaleProfileInputs}. The rescaled arrays/scalars are
 * fresh allocations (typed-array copies), never views into the caller's
 * buffers. `canonicalScalars` reports the canonical units the rescaled
 * fields now carry — always `h` / `mg/L` / `mg`.
 *
 * `warnings` carries non-blocking heuristic advisories surfaced during
 * rescaling: `BW_UNIT_HEURISTIC_WARNING` (species-absent fallback, AD-6),
 * `INFUSION_DURATION_HEURISTIC_WARNING` (AD-10 2× margin), `PER_KG_NO_BW`
 * (nca-cl-vz-convention AD-2: per-kg `/kg` dose without body weight). They
 * propagate into `ComputeResult.provenance.warnings` so they surface in the
 * `WARNINGS_JSON` results column.
 *
 * **`doseNormalization` discriminant (nca-cl-vz-convention AD-2).** Reports
 * whether the rescaled `dose` is absolute canonical `mg` (`'absolute'`) or a
 * per-kg dose left in `mg/(animal·kg)` because no body weight was available
 * (`'per_kg'`). It is per-profile interpretation metadata: the writer carries
 * it into the `DOSE_NORM` result column and {@link normaliseClVolumes} reads it to
 * form the per-kg `CL`/`VZ`/`VSS` and the absolute companions (the unit tags are
 * fixed, nca-cl-bw-normalised-units AD-1). It cannot be re-derived downstream, so it
 * must travel on this struct alongside `route`.
 */
export interface ScaledProfileInputs {
  readonly time: Float64Array;
  readonly conc: Float64Array;
  readonly dose: number;
  readonly lloq: Float64Array | number;
  readonly bodyWeight: number | null;
  readonly infusionDuration: number | null;
  /**
   * `'absolute'` — `dose` is canonical `mg` (BW present, or a non-`/kg` unit).
   * `'per_kg'` — `dose` is per-kg (a `/kg` unit with no body weight; the `×BW`
   * step was skipped per AD-2). Drives the per-kg CL/Vz unit tags downstream.
   */
  readonly doseNormalization: 'absolute' | 'per_kg';
  readonly canonicalScalars: {
    readonly timeUnits: nca.TimeUnit; // 'h'
    readonly concentrationUnits: nca.ConcentrationUnit; // 'mg/L'
    readonly doseUnits: nca.DoseUnit; // 'mg'
  };
  readonly warnings: ReadonlyArray<nca.ParameterWarning>;
}

/** Threshold for the species-absent BW heuristic (AD-6 / R1-finding-2).
 *  BW > 15 with bwUnit='kg' fires the warning — catches mouse (20-40g
 *  reported as 30 → silent 1000× error class), rat (250 → still fires
 *  for safety), and tolerates a 15-kg beagle as a non-blocking advisory.
 *  Threshold = 15 was selected to catch the dominant Persona 1 (Discovery
 *  DMPK) rodent silent-error pattern; 50 (the original synthesis-time
 *  default) missed mice. Lb edge case is a known v1 limitation. */
export const BW_HEURISTIC_THRESHOLD_KG = 15;

/** Margin factor for the infusion-duration heuristic (AD-10 / R2-4).
 *  `infusionDuration > max(time) × INFUSION_DURATION_HEURISTIC_MARGIN`
 *  fires the warning. 2× margin avoids false positives on legitimate
 *  continuous / extended-infusion designs that intentionally span the
 *  observation window. */
export const INFUSION_DURATION_HEURISTIC_MARGIN = 2;

/**
 * Raw-concentration plausibility ceiling for the P1 conc-mis-scale guard
 * (value ratified by a DMPK scientist). When the canonical-default (`mg/L`) concentration's
 * per-profile positive median exceeds this, the data is almost certainly a
 * ng/mL→mg/L 1000× mis-scale (`CONC_MAGNITUDE_HEURISTIC_WARNING`).
 *
 * **Value = 500 mg/L, calibrated empirically vs `files/*.csv` (rule 16).**
 * It sits in the gap between the highest *legitimate* in-repo profile median
 * (06_tk_serial D300 Day28 = 298.3 mg/L) and the lowest mis-scale we must catch
 * (03_rat_iv 5 mg/kg mapped-as-mg/L = 895.2 mg/L), below the ~1000 mg/L
 * aqueous-solubility ceiling. **Documented structural miss (AC-3, R2-F2):** the
 * mid/low-dose + PO/SC mis-scales whose median falls under 500 (03_rat_iv
 * 1 mg/kg ≤ 388.6; 03_rat_po 5 mg/kg ≤ 420.5) do NOT fire here — they are the
 * irreducible overlap zone P2 (implied-CL cross-check) reaches instead. This is
 * a catch-egregious ceiling, not a complete discriminant. */
export const CONC_PLAUSIBILITY_CEILING_MG_L = 500;

/** Minimum positive-observation count below which the P1 conc guard does NOT
 *  fire (nca-conc-dose-misscale-guard AD-1, R1 S-F3). A *structural* floor, not
 *  a fixture-tuned threshold: a "median" of 1–2 positive points is dominated by
 *  a single value and is statistically meaningless (exactly the single-value
 *  mis-computation that flipped the research headline). 3 is the minimum count
 *  at which a median is not just one value. {@link positiveMedian} returns
 *  `null` below it ⇒ no fire. */
export const CONC_GUARD_MIN_POSITIVE_COUNT = 3;

/**
 * Median of the **finite, strictly-positive** values of `arr`, or `null` when
 * fewer than {@link CONC_GUARD_MIN_POSITIVE_COUNT} such values exist.
 *
 * **Positive-only by design (nca-conc-dose-misscale-guard AD-1, METHOD-10).**
 * Excludes NaN, ±Inf, `0`, and negatives — which drops BLQ-substituted zeros
 * and any t=0 zero. The result is a peak-biased central tendency that *aids*
 * mis-scale detection (a profile dominated by sub-LLOQ zeros should not dilute
 * the magnitude signal of its measurable peak) and is documented as such; it is
 * a guard-specific heuristic statistic, NOT the validated NCA median (which has
 * no consumable sci-comp export — see RG-conc-median-upstream).
 *
 * Median = middle element (odd count) or arithmetic mean of the two middle
 * (even count) over the ascending-sorted positive subset. Not barrel-exported
 * (package-internal — test access via direct `../adapter/unit-conversion`
 * import only; NOT re-exported via `adapter/index.ts`, architect gate — single
 * in-package caller).
 */
export function positiveMedian(arr: Float64Array): number | null {
  const pos: number[] = [];
  for (let i = 0; i < arr.length; i++) {
    const v = arr[i];
    if (Number.isFinite(v) && v > 0) pos.push(v);
  }
  if (pos.length < CONC_GUARD_MIN_POSITIVE_COUNT) return null;
  pos.sort((a, b) => a - b);
  const mid = pos.length >> 1;
  return pos.length % 2 === 1 ? pos[mid] : (pos[mid - 1] + pos[mid]) / 2;
}

/**
 * Type narrowing the `scaleProfileInputs` first argument to just the unit
 * fields it actually consumes. `profile-inputs.ts` passes a full
 * `ProfileScalars`; this narrows it so `unit-conversion.ts` does NOT depend
 * on `profile-inputs.ts` (would create a circular import — see
 * architect-gate S-2).
 */
export interface UnitScalars {
  readonly timeUnits: nca.TimeUnit;
  readonly concentrationUnits: nca.ConcentrationUnit;
  readonly doseUnits: nca.DoseUnit;
}

/**
 * Rescale one profile's raw inputs to canonical units (`h` / `mg/L` / `mg` /
 * `kg`). Pure function: inputs are never mutated; outputs are fresh
 * allocations.
 *
 * Inputs:
 * - `scalars` — at minimum `timeUnits` / `concentrationUnits` / `doseUnits`
 *   (type-safe via the sci-comp unit literals; wrong-unit is a compile-time
 *   error).
 * - `time`, `conc` — raw `Float64Array` slices (from `extractProfileInputs`).
 * - `dose` — raw scalar dose; per-kg variants are × `bodyWeight_kg` here.
 * - `lloq` — per-row array OR scalar; rescaled by the **same factor as
 *   concentration** (AD-7 — RATIFIED).
 * - `bodyWeight` — raw scalar in `bodyWeightUnit`; rescaled to kg.
 * - `infusionDuration` — raw scalar; rescaled by the same factor as time
 *   (AD-10).
 * - `bodyWeightUnit` — `'kg'` or `'g'` (default `'kg'` at the caller).
 *
 * Output `ScaledProfileInputs` carries canonical arrays/scalars +
 * `canonicalScalars: {timeUnits:'h', concentrationUnits:'mg/L', doseUnits:'mg'}`
 * + any heuristic warnings.
 *
 * NaN / ±Infinity propagate through Float64 arithmetic — never silently
 * substituted.
 */
export function scaleProfileInputs(
  scalars: UnitScalars,
  time: Float64Array,
  conc: Float64Array,
  dose: number,
  lloq: Float64Array | number,
  bodyWeight: number | null,
  infusionDuration: number | null,
  bodyWeightUnit: BodyWeightUnit,
): ScaledProfileInputs {
  const warnings: nca.ParameterWarning[] = [];

  // Body weight first — dose rescaling for per-kg (mg/kg, µg/kg) units consumes the
  // canonical (kg) value.
  const bwFactor = BW_TO_KG[bodyWeightUnit];
  const bodyWeightKg: number | null = bodyWeight === null ? null : bodyWeight * bwFactor;

  // Species-absent heuristic (AD-6 / R1-finding-2). Strictly-greater
  // comparison at the threshold so an exact 15-kg beagle does NOT fire
  // — the boundary case is acceptable for human dispatch; values above 15
  // (mouse-as-kg, rat-as-kg, dog-edge) trigger the non-blocking warning.
  if (bodyWeightKg !== null && bodyWeightUnit === 'kg' && bodyWeight !== null &&
      bodyWeight > BW_HEURISTIC_THRESHOLD_KG) {
    warnings.push({
      code: 'BW_UNIT_HEURISTIC_WARNING',
      severity: 'warning',
      message:
        `Body weight ${bodyWeight} kg exceeds species-absent threshold ` +
        `${BW_HEURISTIC_THRESHOLD_KG} kg. If the source data is grams ` +
        `reported as kg (common rodent CRO format), CL / Vz will be off ` +
        `by 1000×. Verify the BW unit picker.`,
    });
  }

  // Time → hours.
  const timeFactor = TIME_TO_HOUR[scalars.timeUnits];
  const timeOut = new Float64Array(time.length);
  if (timeFactor === 1)
    timeOut.set(time);
  else
    for (let i = 0; i < time.length; i++) timeOut[i] = time[i] * timeFactor;


  // Concentration → mg/L.
  const concFactor = convertConcentrationFactor(scalars.concentrationUnits);
  const concOut = new Float64Array(conc.length);
  if (concFactor === 1)
    concOut.set(conc);
  else
    for (let i = 0; i < conc.length; i++) concOut[i] = conc[i] * concFactor;

  // P1 conc-mis-scale guard (nca-conc-dose-misscale-guard AD-1). Fires ONLY
  // when the user left the canonical default `mg/L` — i.e. they did NOT pick a
  // non-canonical unit (ng/mL etc.). When a non-canonical unit was picked the
  // data was already rescaled ÷1000 here and the high magnitude is gone (AC-4):
  // a correctly-labelled ng/mL profile must never fire. The per-profile
  // positive median (peak-biased; null below the count floor) separates a
  // legitimate high-dose TK study (≤298 mg/L in-repo) from a ng/mL→mg/L 1000×
  // mis-scale; above the ceiling the inflated concentration would silently
  // corrupt CL/Vz/AUC while λz/t½ (scale-invariant) look fine. Non-blocking —
  // legit 298 mg/L in-repo forbids a hard block (BW precedent, not LLOQ).
  if (scalars.concentrationUnits === 'mg/L') {
    const concMedian = positiveMedian(concOut);
    if (concMedian !== null && concMedian > CONC_PLAUSIBILITY_CEILING_MG_L) {
      warnings.push({
        code: 'CONC_MAGNITUDE_HEURISTIC_WARNING',
        severity: 'warning',
        message:
          `Concentration median ${concMedian.toPrecision(4)} mg/L exceeds the ` +
          `plausibility ceiling ${CONC_PLAUSIBILITY_CEILING_MG_L} mg/L. If the ` +
          `source data is ng/mL reported as mg/L (common bioanalytical CRO ` +
          `format), CL / Vz / AUC will be off by 1000× while λz / t½ look ` +
          `correct. Verify the concentration unit picker.`,
      });
    }
  }

  // LLOQ — same factor as concentration (AD-7 ratified). Array OR scalar.
  let lloqOut: Float64Array | number;
  if (typeof lloq === 'number')
    lloqOut = concFactor === 1 ? lloq : lloq * concFactor;
  else {
    if (concFactor === 1) {
      lloqOut = new Float64Array(lloq.length);
      lloqOut.set(lloq);
    } else {
      const out = new Float64Array(lloq.length);
      for (let i = 0; i < lloq.length; i++) out[i] = lloq[i] * concFactor;
      lloqOut = out;
    }
  }

  // Dose → mg. Per-kg variants (mg/kg, µg/kg) consume canonical `bodyWeightKg`
  // to form an absolute mg dose when body weight is available. When a `/kg`
  // dose unit arrives WITHOUT a body weight (`bodyWeightKg === null`), this is
  // NOT a defect — it is the ratified discovery-PK per-kg path
  // (nca-cl-vz-convention AD-2): the `×BW` step is skipped, the dose stays
  // per-kg (`mg/(animal·kg)` after the base factor), `doseNormalization` is
  // marked `'per_kg'`, and a `PER_KG_NO_BW` advisory is pushed so the
  // downstream CL/Vz units (L/h/kg, L/kg) are auditable. The
  // `MISSING_BODY_WEIGHT` validator no longer blocks this case.
  const baseDoseFactor = DOSE_TO_MG[scalars.doseUnits];
  let doseOut: number;
  let doseNormalization: 'absolute' | 'per_kg' = 'absolute';
  // INVARIANT (nca-conc-dose-misscale-guard R2-F1): `mg/kg` / `µg/kg` are the
  // ONLY per-kg entries in DOSE_TO_MG, so this branch is the sole producer of
  // `doseNormalization === 'per_kg'`, and per_kg ⟺ BW-absent by construction
  // (the `bodyWeightKg === null` sub-branch below). The P2 implied-CL guard
  // (cl-plausibility.ts) and `normaliseClVolumes` RELY on this: they divide CL by
  // body weight only on the `'absolute'` path, where the engine's CL is in L/h, and
  // read a `'per_kg'` CL as is — a single normalization, never double. A future
  // `/kg` dose unit added to DOSE_TO_MG
  // MUST extend this conditional in lockstep, or the per_kg⟺BW-absent invariant
  // silently breaks and P2 could double-normalize.
  if (isPerKgDose(scalars.doseUnits)) {
    if (bodyWeightKg === null) {
      doseOut = dose * baseDoseFactor; // per-kg dose, no ×BW
      doseNormalization = 'per_kg';
      // Per kg is the reported basis everywhere (nca-cl-bw-normalised-units D-a/c), so the
      // advisory names what is missing — the absolute companions — not the basis.
      warnings.push({
        code: 'PER_KG_NO_BW',
        severity: 'warning',
        message:
          `No body weight: absolute CL, Vz and Vss (CL abs., Vz abs., Vss abs.) are not ` +
          `available; the per-kg values are reported. If the dose column ('${scalars.doseUnits}') ` +
          `is actually absolute, re-map the dose unit.`,
      });
    } else
      doseOut = dose * baseDoseFactor * bodyWeightKg;
  } else
    doseOut = dose * baseDoseFactor;


  // Infusion duration — same factor as time (AD-10 / E10).
  let infusionOut: number | null;
  if (infusionDuration === null)
    infusionOut = null;
  else {
    infusionOut = timeFactor === 1 ? infusionDuration : infusionDuration * timeFactor;

    // 2× heuristic on the ORIGINAL (pre-rescale) max(time) — since
    // infusionDuration is rescaled by the same factor, the ratio is
    // invariant under rescaling. Computing on canonical values keeps
    // the check unit-independent.
    let maxTimeCanonical = -Infinity;
    for (let i = 0; i < timeOut.length; i++) {
      const v = timeOut[i];
      if (!Number.isNaN(v) && v > maxTimeCanonical) maxTimeCanonical = v;
    }
    if (Number.isFinite(maxTimeCanonical) && Number.isFinite(infusionOut) &&
        infusionOut > maxTimeCanonical * INFUSION_DURATION_HEURISTIC_MARGIN) {
      warnings.push({
        code: 'INFUSION_DURATION_HEURISTIC_WARNING',
        severity: 'warning',
        message:
          `Infusion duration ${infusionOut} h exceeds ${INFUSION_DURATION_HEURISTIC_MARGIN}× ` +
          `the observation window (max time ${maxTimeCanonical} h) — verify units. ` +
          `The infusion duration is assumed to share the time picker's unit.`,
      });
    }
  }

  return {
    time: timeOut,
    conc: concOut,
    dose: doseOut,
    lloq: lloqOut,
    bodyWeight: bodyWeightKg,
    infusionDuration: infusionOut,
    doseNormalization,
    canonicalScalars: {
      timeUnits: 'h',
      concentrationUnits: 'mg/L',
      doseUnits: 'mg',
    },
    warnings,
  };
}
