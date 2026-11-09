/**
 * Profile-inputs extractor — builds `nca.ProfileInputs` for one profile by
 * copying the right slice of an observations DataFrame into fresh
 * `Float64Array` / `Uint8Array` instances.
 *
 * **Copy, never view.** Sci-comp keeps references to incoming `Float64Array` /
 * `Uint8Array` buffers (see `nca/core/compute-nca.ts`). Handing it a
 * DataFrame's `getRawData()` view would couple core lifetime to Datagrok GC
 * — the spec §10 boundary discipline rule rules this out. Copies cost
 * microseconds at NCA sizes (< 100 points per profile).
 *
 * **Float32 widening (spec §6.1).** When a numeric input column is
 * Float32-backed (e.g. created via `Column.float()` or `addNewFloat()`), the
 * extractor widens by copying into a fresh `Float64Array` and pushes a
 * `FLOAT32_WIDENED` warning into the returned `warnings` array. PKNCA
 * bit-identity (NFR-04) is not guaranteed for widened profiles; the warning
 * makes the loss auditable. This is a soft warning, NOT a `ValidationError`
 * — the validator set stays at the 5 codes per spec §3.1 row 5.
 *
 * **Route mapping (cycle-state phase_3_plan resolved_decisions #6).** The
 * adapter is the single source of truth for the string→`RouteCode` mapping;
 * the core never sees route strings. Unknown route strings throw at the
 * adapter layer — they are not part of the 5-code surface. A structured
 * `BAD_ROUTE` code can be added later when column-mapping UI (Phase 5)
 * surfaces this case.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';
import {nca} from '@datagrok-libraries/sci-comp';

import type {ColumnContract} from './column-contract';
import type {ProfileIndex} from './profile-index';
import type {BodyWeightContext, BodyWeightUnit} from './unit-conversion';
import {BW_TO_KG, isPerKgDose, scaleProfileInputs} from './unit-conversion';

/**
 * Scalar fallbacks used when the corresponding optional column is absent
 * from the contract or the DataFrame. The (Phase 4) controller is
 * responsible for collecting these from the column-mapping dialog and
 * passing them in.
 *
 * `lloq` may be a per-row column (when present in the contract) OR a scalar
 * carried here; the same field name keeps the call-site shape symmetric.
 */
export interface ProfileScalars {
  readonly dose: number;
  readonly route: nca.RouteCode;
  readonly doseUnits: nca.DoseUnit;
  readonly concentrationUnits: nca.ConcentrationUnit;
  readonly timeUnits: nca.TimeUnit;
  readonly bodyWeight: number | null;
  /**
   * Body weight unit (`'kg'` | `'g'`). Default `'kg'`; `'g'` for rodent
   * studies reporting grams. Drives `scaleProfileInputs` rescaling to
   * canonical kg (AD-6). Triggers the species-absent heuristic warning
   * (`BW_UNIT_HEURISTIC_WARNING`) when `bodyWeight > 15` with
   * `bodyWeightUnit === 'kg'` — catches the mouse-as-kg silent 1000× error.
   */
  readonly bodyWeightUnit: BodyWeightUnit;
  readonly infusionDuration: number | null;
  /**
   * Lower limit of quantification scalar (in `concentrationUnits` pre-rescale).
   *
   * **"0 = unset" convention.** The mapping pane defaults `lloq` to `0` when
   * the user has not entered a value (the same form-level convention used by
   * `bodyWeight` and `infusionDuration`, except those collapse to `null` via
   * `scalarOrNull`; `lloq` keeps the `number` type for legacy reasons). The
   * scalar value is consumed by `extractProfileInputs` only when no LLOQ
   * column is mapped (`lloqCol === null` path). It does NOT drive BLQ
   * inference — BLQ-mask construction is gated on `blqFlagCol !== null` OR
   * `lloqCol !== null`, never on the scalar value (see `profile-inputs.ts`
   * lines 145-159). Therefore `lloq=0` cannot silently affect downstream
   * BLQ classification.
   *
   * Sci-comp `set-half-lloq` would compute `lloq/2 = 0`, but only for rows
   * that the precomputed `blqMask` marked BLQ. With no `blqFlagCol` and no
   * `lloqCol`, the mask is all-zeros and `lloq=0` is inert.
   *
   * The LLOQ sanity validator (`validateLloqSanity`) treats `lloq=0` with
   * no `lloqCol` mapped as "not specified" and skips per-profile sanity
   * conditions (a/c/d/e). This matches the form convention.
   */
  readonly lloq: number;
}

/**
 * Adapter→core route mapping (string → `RouteCode`). Locked-in by spec §7.1
 * row Route and cycle-state phase_3_plan resolved_decisions #6.
 */
export const ROUTE_STRING_TO_CODE: Readonly<Record<string, nca.RouteCode>> = {
  'IV-bolus': nca.ROUTE_IV_BOLUS,
  'IV-infusion': nca.ROUTE_IV_INFUSION,
  'PO': nca.ROUTE_PO,
  'SC': nca.ROUTE_SC,
  'IM': nca.ROUTE_IM,
  'other': nca.ROUTE_OTHER,
};

/**
 * Translate a route string to its integer `RouteCode`. Throws on unknown
 * input — UC-01 has no `BAD_ROUTE` `ValidationCode`, so the validator
 * cannot surface this; it is a defect upstream (column-mapping dialog
 * should constrain the string to the closed set, Phase 5).
 */
export function routeStringToCode(route: string): nca.RouteCode {
  const code = ROUTE_STRING_TO_CODE[route];
  if (code === undefined) {
    const known = Object.keys(ROUTE_STRING_TO_CODE).join(', ');
    throw new Error(`Unknown route string '${route}'; expected one of ${known}.`);
  }
  return code;
}

/**
 * Build a single profile's `nca.ProfileInputs` plus any extractor warnings.
 *
 * Inputs:
 * - `df` — observations DataFrame (post-validation; the caller should have
 *   already run `validateObservations` and `validateUnits` and aborted on
 *   errors).
 * - `index` — `ProfileIndex` for the one profile to extract.
 * - `contract` — the column contract (tells the extractor which optional
 *   columns are present and what they're called in the DataFrame).
 * - `scalars` — fallbacks for dose/route/units/body weight/lloq/infusion
 *   duration when the corresponding column is absent, plus `bodyWeightUnit`
 *   (`'kg'` | `'g'`).
 *
 * Returns `{inputs, warnings}`. The controller (Phase 4) merges `warnings`
 * into the downstream `ComputeResult.provenance.warnings` per
 * `phase_3_plan.acceptance` row 5.
 *
 * **Pipeline-bypass invariant (AD-12 / G9 / architecture spec invariant #7).**
 * All `nca.ProfileInputs` MUST be constructed via `extractProfileInputs` →
 * `scaleProfileInputs`. Direct construction of `nca.ProfileInputs` from raw
 * (non-canonical) values bypasses the canonical rescaling step and is
 * forbidden. `nca.ProfileInputs` is a plain TypeScript interface (sci-comp
 * public), so this is enforced by code review + the architecture spec
 * invariant rather than at compile time. Future hardening via a named-builder
 * refactor is a separate cycle if/when bypass risk materialises empirically.
 *
 * **Canonical units invariant (AD-1 / AD-2).** The returned `inputs.timeUnits`
 * / `inputs.concentrationUnits` / `inputs.doseUnits` are ALWAYS `'h'` /
 * `'mg/L'` / `'mg'` (canonical), regardless of the user's mapping-pane
 * selections. `scaleProfileInputs` is the single rescaling point.
 */
export function extractProfileInputs(
  df: DG.DataFrame,
  index: ProfileIndex,
  contract: ColumnContract,
  scalars: ProfileScalars,
): {
  inputs: nca.ProfileInputs;
  warnings: nca.ParameterWarning[];
  /** Per-profile dose-normalization basis from `scaleProfileInputs`
   *  (nca-cl-vz-convention AD-2): `'absolute'` (canonical mg) or `'per_kg'`
   *  (a `/kg` dose with no body weight → CL/Vz reported per kg). Threaded to
   *  the results writer like `route`. */
  doseNormalization: 'absolute' | 'per_kg';
  /** The profile's body weight as the per-kg rule needs it (nca-cl-bw-normalised-units AD-2):
   *  the engine's first-row weight, the range over the profile's non-missing weights, and
   *  whether the dose was per kg. Threaded to the results writer beside `doseNormalization`. */
  bodyWeight: BodyWeightContext;
} {
  const warnings: nca.ParameterWarning[] = [];
  const rows = index.observationRows;
  const obs = contract.observations;

  const timeCol = mustGetCol(df, obs.time.name);
  const concCol = mustGetCol(df, obs.concentration.name);
  const rawTime = extractFloat64(timeCol, rows, warnings);
  const rawConc = extractFloat64(concCol, rows, warnings);

  const blqFlagCol = obs.blqFlag !== undefined ? df.col(obs.blqFlag.name) : null;
  const lloqCol = obs.lloq !== undefined ? df.col(obs.lloq.name) : null;

  let blqMask: Uint8Array;
  if (blqFlagCol !== null)
    blqMask = extractUint8(blqFlagCol, rows);
  else if (lloqCol !== null) {
    // Spec §7.1: BLQ default is "inferred: conc === 0 treated as BLQ when LLOQ
    // column is present, else all 0s". The LLOQ-column gate is intentional —
    // a global LLOQ scalar does not justify per-row BLQ inference.
    // Inference uses the RAW (pre-rescale) concentration: zero is zero in
    // any unit, so the BLQ mask is rescale-invariant.
    blqMask = new Uint8Array(rows.length);
    for (let i = 0; i < rows.length; i++)
      blqMask[i] = (rawConc[i] === 0) ? 1 : 0;
  } else
    blqMask = new Uint8Array(rows.length);


  const rawLloq: Float64Array | number = lloqCol !== null ?
    extractFloat64(lloqCol, rows, warnings) :
    scalars.lloq;

  const firstRow = rows[0];

  const rawDose = obs.dose !== undefined && df.col(obs.dose.name) !== null ?
    readNumericCell(df.col(obs.dose.name)!, firstRow) :
    scalars.dose;

  // Route per-row read: route column is rawType 'int32' per ColumnContract
  // (Phase 5 column-mapping dialog encodes the user's string column to
  // RouteCode before the contract is built). When the cell is INT_NULL,
  // `readNumericCell` returns NaN; `NaN | 0 = 0 = ROUTE_IV_BOLUS` would
  // silently mislabel CL/F output. Fall back to the dialog's scalar route
  // when the cell is null.
  let route: nca.RouteCode;
  if (obs.route !== undefined && df.col(obs.route.name) !== null) {
    const raw = readNumericCell(df.col(obs.route.name)!, firstRow);
    route = Number.isNaN(raw) ? scalars.route : ((raw | 0) as nca.RouteCode);
  } else
    route = scalars.route;


  // Per-subject body weight (mirrors the per-row route / infusion reads below).
  // A blank/null cell on a mapped BW column reads as NaN via `readNumericCell`
  // (FLOAT_NULL→NaN), so it MUST be coerced to `null` — not passed through. The
  // `null` is the only value `scaleProfileInputs` treats as BW-absent (its
  // `bodyWeight === null` check, unit-conversion.ts), which routes a mg/kg dose
  // down the ratified per-kg path (DOSE_NORM='per_kg', PER_KG_NO_BW). Passing the
  // raw NaN instead leaks `bodyWeightKg = NaN` (≠ null): the per-kg branch is
  // skipped, DOSE_NORM is mislabeled 'absolute', and the dose becomes
  // dose×factor×NaN = NaN — silently corrupting every dose-dependent parameter
  // for that subject. (blank-BW-dosenorm-fix; honest per_kg fallback — we do NOT
  // fabricate an absolute dose from the global scalar weight for an unknown BW.)
  let rawBodyWeight: number | null = scalars.bodyWeight;
  // Range over the profile's NON-MISSING weights (AD-2, R2 N-2): a blank cell is ignored, never
  // read as a different weight. The first-row read above stays the engine's weight.
  let bwRangeRaw: {min: number; max: number} | null = null;
  if (obs.bodyWeight !== undefined && df.col(obs.bodyWeight.name) !== null) {
    const bwCol = df.col(obs.bodyWeight.name)!;
    const raw = readNumericCell(bwCol, firstRow);
    rawBodyWeight = Number.isFinite(raw) ? raw : null;
    for (const r of rows) {
      const v = readNumericCell(bwCol, r);
      if (!Number.isFinite(v)) continue;
      bwRangeRaw = bwRangeRaw === null ? {min: v, max: v} :
        {min: Math.min(bwRangeRaw.min, v), max: Math.max(bwRangeRaw.max, v)};
    }
  } else if (rawBodyWeight !== null && Number.isFinite(rawBodyWeight))
    bwRangeRaw = {min: rawBodyWeight, max: rawBodyWeight};

  // Per-subject infusion duration (Slice A2 / BD-2): when an infusionDuration
  // column is on the contract (the event-driven `~InfusionDuration` column),
  // read this profile's value at its first row. A null/NaN cell (RATE=0 subject)
  // falls back to the dialog scalar — exactly mirroring the per-row `route` read
  // above, so a single global scalar no longer overwrites every subject's
  // duration. Absent column ⇒ pure scalar path (dense studies byte-unchanged).
  let rawInfusionDuration: number | null = scalars.infusionDuration;
  if (obs.infusionDuration !== undefined && df.col(obs.infusionDuration.name) !== null) {
    const raw = readNumericCell(df.col(obs.infusionDuration.name)!, firstRow);
    if (Number.isFinite(raw)) rawInfusionDuration = raw;
  }

  // Single rescaling point (AD-1 / AD-12 pipeline-bypass invariant). Outputs
  // canonical arrays/scalars + heuristic warnings (BW_UNIT_HEURISTIC_WARNING,
  // INFUSION_DURATION_HEURISTIC_WARNING).
  const scaled = scaleProfileInputs(
    scalars,
    rawTime,
    rawConc,
    rawDose,
    rawLloq,
    rawBodyWeight,
    rawInfusionDuration,
    scalars.bodyWeightUnit,
  );
  for (const w of scaled.warnings) warnings.push(w);

  const inputs: nca.ProfileInputs = {
    time: scaled.time,
    conc: scaled.conc,
    blqMask,
    lloq: scaled.lloq,
    dose: scaled.dose,
    doseUnits: scaled.canonicalScalars.doseUnits,
    concentrationUnits: scaled.canonicalScalars.concentrationUnits,
    timeUnits: scaled.canonicalScalars.timeUnits,
    route,
    infusionDuration: scaled.infusionDuration,
    bodyWeight: scaled.bodyWeight,
  };
  const bwFactor = BW_TO_KG[scalars.bodyWeightUnit];
  const bodyWeight: BodyWeightContext = {
    bodyWeightKg: scaled.bodyWeight,
    bodyWeightRangeKg: bwRangeRaw === null ? null :
      {min: bwRangeRaw.min * bwFactor, max: bwRangeRaw.max * bwFactor},
    perKgDose: isPerKgDose(scalars.doseUnits),
  };
  return {inputs, warnings, doseNormalization: scaled.doseNormalization, bodyWeight};
}

/**
 * Extract a per-row Float64 slice for one profile from a Datagrok column.
 *
 * Discriminator switch on `getRawData()` constructor:
 * - `Float64Array` — fast path, copy with `FLOAT_NULL → NaN` translation.
 * - `Float32Array` — widen-copy + push `FLOAT32_WIDENED` warning.
 * - `Int32Array`   — int → float copy with `INT_NULL → NaN`.
 * - other          — throw (unsupported buffer type for a numeric column).
 *
 * Output is a freshly-allocated `Float64Array` of `rows.length`; no view into
 * the DataFrame's storage is retained.
 */
function extractFloat64(
  col: DG.Column,
  rows: Int32Array,
  warnings: nca.ParameterWarning[],
): Float64Array {
  const raw = col.getRawData();
  const out = new Float64Array(rows.length);
  if (raw instanceof Float64Array) {
    for (let i = 0; i < rows.length; i++) {
      const v = raw[rows[i]];
      out[i] = (v === DGNS.FLOAT_NULL) ? NaN : v;
    }
  } else if (raw instanceof Float32Array) {
    for (let i = 0; i < rows.length; i++) {
      const v = raw[rows[i]];
      out[i] = (v === DGNS.FLOAT_NULL) ? NaN : v;
    }
    warnings.push({
      code: 'FLOAT32_WIDENED',
      severity: 'warning',
      message: `Column ${col.name} widened from Float32 to Float64; PKNCA bit-identity not guaranteed (NFR-04)`,
    });
  } else if (raw instanceof Int32Array) {
    for (let i = 0; i < rows.length; i++) {
      const v = raw[rows[i]];
      out[i] = (v === DGNS.INT_NULL) ? NaN : v;
    }
  } else
    throw new Error(`extractFloat64: column '${col.name}' has unsupported raw buffer type '${raw.constructor.name}'.`);

  return out;
}

/**
 * Extract a per-row Uint8 mask for one profile from a BOOL column. BOOL
 * storage in Datagrok is a `Uint32Array` bit array (not directly indexable as
 * bytes), so we go through the per-cell accessor — at NCA profile sizes (< 100
 * points) the loop is microseconds.
 */
function extractUint8(col: DG.Column, rows: Int32Array): Uint8Array {
  const out = new Uint8Array(rows.length);
  for (let i = 0; i < rows.length; i++)
    out[i] = col.get(rows[i]) ? 1 : 0;
  return out;
}

function readNumericCell(col: DG.Column, row: number): number {
  const raw = col.getRawData();
  const v = raw[row];
  if (raw instanceof Float64Array || raw instanceof Float32Array)
    return (v === DGNS.FLOAT_NULL) ? NaN : v;
  if (raw instanceof Int32Array)
    return (v === DGNS.INT_NULL) ? NaN : v;
  return v;
}

function mustGetCol(df: DG.DataFrame, name: string): DG.Column {
  const col = df.col(name);
  if (col === null) {
    throw new Error(
      `extractProfileInputs: required column '${name}' not present in DataFrame ` +
      `(should have been caught by validators).`,
    );
  }
  return col;
}
