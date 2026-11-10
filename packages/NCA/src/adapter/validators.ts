/**
 * Hand-written validators for the UC-01 observations DataFrame.
 *
 * Scope (UC-01): the 5 codes locked in spec §3.1 row 5 — `MISSING_COLUMN`,
 * `WRONG_TYPE`, `NEGATIVE_TIME`, `NON_NUMERIC_CONC`, `DUPLICATE_SUBJECT_TIME`.
 * Interface §2 lists three additional codes (`INVALID_ROUTE_CODE`,
 * `UNITS_MISMATCH`, `INFUSION_DURATION_REQUIRED`); they are out of UC-01 scope
 * (route-string lookup throws at the adapter layer per cycle-state Phase 3
 * resolved_decisions #6, infusion duration is not required for the UC-01
 * theoph/indometh fixtures, and units are entered in the mapping form rather
 * than parsed from the DataFrame).
 *
 * `NON_NUMERIC_CONC` interpretation: within a Datagrok numeric column the only
 * non-finite, non-NaN value is `±Infinity`; spec §7.1 row Concentration requires
 * "finite or NaN", so `±Infinity` is the sole violation that is not already
 * covered by `WRONG_TYPE` (which fires when the user maps a non-numeric column
 * to concentration).
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import type {
  BooleanRawArray,
  ColumnContract,
  ColumnRawType,
  ColumnSpec,
  IntegerRawArray,
  NumericRawArray,
  StringRawArray,
  ValidationError,
} from './column-contract';
import type {ProfileScalars} from './profile-inputs';
import type {ProfileIndex} from './profile-index';
import {detectAnalyteColumns, detectMixedCompoundGroups} from './column-autodetect';
import {stringifyKey} from './adapter-utils';

type AnyColumnSpec = ColumnSpec<NumericRawArray | IntegerRawArray | BooleanRawArray | StringRawArray>;

// `rawType='float64'` accepts BOTH FLOAT and INT — extractor (profile-inputs.ts
// extractFloat64) translates INT_NULL→NaN and widens INT columns to Float64 at
// extraction, so an INT-backed numeric column (e.g. an integer-coded time axis
// from a legacy CSV) is a valid input. Closes NOTE-04 (Phase 3 design decision
// #3: extractor permissive on numeric source type; validator was stricter than
// the extractor, blocking realistic legacy-CSV shapes).
const COLUMN_TYPE_BY_RAW: Record<ColumnRawType, ReadonlyArray<string>> = {
  float64: [DGNS.COLUMN_TYPE.FLOAT, DGNS.COLUMN_TYPE.INT],
  int32: [DGNS.COLUMN_TYPE.INT],
  uint8: [DGNS.COLUMN_TYPE.BOOL],
  string: [DGNS.COLUMN_TYPE.STRING],
};

/**
 * Validate one observations DataFrame against a UC-01 column contract. Returns
 * an empty array when valid; otherwise a flat list of every failure (the
 * controller (Phase 4) renders these in the column-mapping error panel — spec
 * §6.1 row 2 — and blocks compute).
 *
 * The validator is allocation-light: each column is scanned at most twice
 * (once for type/missing checks, once for row-level scans), and the
 * `DUPLICATE_SUBJECT_TIME` pass uses `Map<string, Set<number>>` rather than
 * sorting. Costs are O(n_rows · n_required_columns).
 *
 * Validators never throw; structured `ValidationError[]` is the only output
 * channel. Conditions that would have thrown (e.g. unknown route string) live
 * in `extractProfileInputs` instead — they are not part of the 5-code surface.
 */
export function validateObservations(
  df: DG.DataFrame,
  contract: ColumnContract,
): ValidationError[] {
  const errs: ValidationError[] = [];
  const obs = contract.observations;

  const subjectErrs = checkColumnSpec(df, obs.subjectId);
  errs.push(...subjectErrs);
  const timeErrs = checkColumnSpec(df, obs.time);
  errs.push(...timeErrs);
  const concErrs = checkColumnSpec(df, obs.concentration);
  errs.push(...concErrs);

  for (const optional of [obs.dose, obs.route, obs.bodyWeight, obs.blqFlag, obs.lloq,
    obs.infusionDuration]) {
    if (optional !== undefined && df.col(optional.name) !== null)
      errs.push(...checkColumnSpec(df, optional));
  }

  if (timeErrs.length === 0)
    errs.push(...scanNegativeTime(df, obs.time.name));

  if (concErrs.length === 0)
    errs.push(...scanNonNumericConc(df, obs.concentration.name));

  if (subjectErrs.length === 0 && timeErrs.length === 0)
    errs.push(...scanDuplicateSubjectTime(df, obs.subjectId.name, obs.time.name));

  return errs;
}

function checkColumnSpec(
  df: DG.DataFrame,
  spec: AnyColumnSpec,
): ValidationError[] {
  const col = df.col(spec.name);
  if (col === null) {
    if (spec.required) {
      return [{
        code: 'MISSING_COLUMN',
        message: `Required column '${spec.name}' is missing from the input DataFrame.`,
        column: spec.name,
      }];
    }
    return [];
  }
  const expected = COLUMN_TYPE_BY_RAW[spec.rawType];
  if (!expected.includes(col.type)) {
    const expectedStr = expected.length === 1 ? expected[0] : expected.join(' | ');
    return [{
      code: 'WRONG_TYPE',
      message: `Column '${spec.name}' has type '${col.type}'; expected '${expectedStr}' (raw '${spec.rawType}').`,
      column: spec.name,
    }];
  }
  return [];
}

function scanNegativeTime(df: DG.DataFrame, name: string): ValidationError[] {
  const errs: ValidationError[] = [];
  const col = df.col(name);
  if (col === null) return errs;
  const raw = col.getRawData();
  const n = df.rowCount;
  for (let i = 0; i < n; i++) {
    const v = raw[i];
    if (v === DGNS.FLOAT_NULL) continue;
    if (Number.isNaN(v)) continue;
    if (v < 0) {
      errs.push({
        code: 'NEGATIVE_TIME',
        message: `Negative time '${v}' at row ${i} in column '${name}'.`,
        column: name,
        row: i,
      });
    }
  }
  return errs;
}

function scanNonNumericConc(df: DG.DataFrame, name: string): ValidationError[] {
  const errs: ValidationError[] = [];
  const col = df.col(name);
  if (col === null) return errs;
  const raw = col.getRawData();
  const n = df.rowCount;
  for (let i = 0; i < n; i++) {
    const v = raw[i];
    if (v === DGNS.FLOAT_NULL) continue;
    if (Number.isFinite(v) || Number.isNaN(v)) continue;
    errs.push({
      code: 'NON_NUMERIC_CONC',
      message: `Non-finite, non-NaN concentration '${v}' at row ${i} in column '${name}'.`,
      column: name,
      row: i,
    });
  }
  return errs;
}

/**
 * Validate unit-load-bearing pre-conditions on `(df, contract, scalars)`
 * before extract / compute. Returns an empty array when the inputs are
 * consistent with the unit semantics ratified by `nca-unit-semantics`;
 * otherwise a flat list of `ValidationError`s the controller surfaces in
 * the mapping-pane error panel.
 *
 * Four validators (E6 / AD-4 / AD-8 / AD-9):
 *  1. {@link validateScalarsForDose} — `/kg` dose-unit hook. No longer blocks
 *     (nca-cl-vz-convention P-A / AD-1): `/kg` dose without body weight is the
 *     valid per-kg path, advised non-blockingly via the `PER_KG_NO_BW` warning
 *     on `provenance.warnings`, not `MISSING_BODY_WEIGHT`. Returns `[]`.
 *  2. {@link validateScalarsForMolar} — molar concentration unit ⇒
 *     `UNITS_MISMATCH` with Path C `details.kind = 'molar-not-supported'`.
 *  3. {@link validateMultiAnalyte} — multi-analyte CSV ⇒
 *     `MULTI_ANALYTE_NOT_SUPPORTED`.
 *  4. {@link validateLloqSanity} — per-profile AD-8 sanity conditions
 *     (a/c/d/e) ⇒ `LLOQ_OUT_OF_RANGE`.
 *
 * Controller invokes this AFTER `validateObservations` and AFTER
 * `buildProfileIndices` (the LLOQ validator needs the indices) and BEFORE
 * `extractProfileInputs`. Non-empty result ⇒ `{status:'validation-failed'}`.
 */
export function validateUnits(
  df: DG.DataFrame,
  contract: ColumnContract,
  scalars: ProfileScalars,
  indices: ReadonlyArray<ProfileIndex>,
): ValidationError[] {
  const errs: ValidationError[] = [];
  errs.push(...validateScalarsForDose(contract, scalars));
  errs.push(...validateScalarsForMolar(scalars));
  errs.push(...validateMultiAnalyte(df, contract.observations.subjectId.name));
  errs.push(...validateLloqSanity(df, contract, scalars, indices));
  return errs;
}

/**
 * `/kg` dose-unit body-weight check.
 *
 * **No longer a hard block (nca-cl-vz-convention AD-1, ratified P-A).** A `/kg`
 * dose unit (`mg/kg`, `µg/kg`) WITHOUT a body-weight source used to emit
 * `MISSING_BODY_WEIGHT` and block compute. It no longer does. The ratified
 * discovery-PK convention reports per-kg `CL` (mL/min/kg) and `Vz` (L/kg) natively
 * in that case (PKNCA `pk.calc.clearance` pass-through precedent). The advisory
 * is not the validator's job: `scaleProfileInputs` skips the `×BW` step, marks
 * `doseNormalization='per_kg'`, and pushes a non-blocking `PER_KG_NO_BW` warning
 * onto `provenance.warnings` (→ `WARNINGS_JSON`). (The pre-Run per-kg confirmation
 * and its `requiresPerKgConfirmation` predicate were retired when per kg became the
 * reported basis everywhere.) This
 * validator therefore returns `[]` for every dose unit — it is kept as the named
 * hook for any future dose-scalar rule.
 *
 * "Body-weight source present" (column path takes precedence at extract time
 * per `profile-inputs.ts`; a positive `scalars.bodyWeight` is the scalar path)
 * no longer changes the outcome — both the present and absent cases are now
 * valid — but the distinction still governs which path `scaleProfileInputs`
 * takes (absolute vs per-kg).
 */
export function validateScalarsForDose(
  _contract: ColumnContract,
  _scalars: ProfileScalars,
): ValidationError[] {
  // AD-1: per-kg `/kg` dose without body weight is a valid (per-kg) path, not a
  // block. No dose-scalar condition currently rejects; the relaxed advisory
  // rides `provenance.warnings` (PER_KG_NO_BW) from `scaleProfileInputs`.
  return [];
}

/**
 * Molar concentration unit (`µmol/L` / `nmol/L`) ⇒ `UNITS_MISMATCH` with
 * Path C `details.kind = 'molar-not-supported'`. Molar→mass conversion
 * requires drug molecular weight, which is not part of the v1 column
 * contract; the architectural pin for MW lives at `ProfileScalars`
 * (forward-deferred to UC-05).
 *
 * Reusing the canonical `UNITS_MISMATCH` from sci-comp interface §2 over
 * a bespoke `UNITS_MOLAR_UNSUPPORTED` code (AD-4 Path C) satisfies CLAUDE.md
 * rule 5 — molar IS unit-incompatibility — and keeps the discriminator
 * extensible for future unit-mismatch sub-kinds (e.g. mixed-unit cohort
 * detection, forward-deferred per E1).
 */
export function validateScalarsForMolar(scalars: ProfileScalars): ValidationError[] {
  const u = scalars.concentrationUnits;
  if (u !== 'µmol/L' && u !== 'nmol/L') return [];
  return [{
    code: 'UNITS_MISMATCH',
    message:
      `Concentration unit '${u}' is molar; mass-based canonical (mg/L) ` +
      `conversion requires the drug molecular weight, which is not part of ` +
      `the v1 column contract. Convert your CSV to a mass-based unit (e.g. ` +
      `ng/mL, µg/L) and re-upload.`,
    column: 'concentration',
    details: {kind: 'molar-not-supported', column: 'concentration', unit: u},
  }];
}

/**
 * Multi-analyte CSV ⇒ `MULTI_ANALYTE_NOT_SUPPORTED`. Hard-block before
 * partial-fill autodetect produces confusing per-profile results.
 *
 * Detection is delegated to {@link detectAnalyteColumns} in
 * `column-autodetect.ts` (single source of truth; tested independently).
 * Two signals trigger: (a) ≥ 2 columns matching analyte-name regexes
 * (`Parent_*`, `Metabolite_*`, `M\d+_Conc`), (b) per-row STRING column
 * named `analyte`/`compound` with cardinality ≥ 2.
 *
 * AD-9 / E6 / scientist memo §3.
 *
 * **Cassette (BUG-07 interim, 2026-10-05).** A third signal: a column named like a compound
 * identifier (`CompoundID`, `Analyte Name`, `Drug`, `Test Article`, …) whose value VARIES within
 * a subject ({@link detectMixedCompoundGroups} keyed on `subjectName`). Profiles are keyed on the
 * subject alone, so that subject's compounds would pool into ONE profile and every parameter
 * would be wrong with no signal. Same code, `details.kind = 'cassette-shared-subject'`; a column
 * signal (a)/(b) already reported is not reported twice. The root fix — a composite
 * (subject × compound × occasion) profile key — is OQ-2, a separate cycle.
 */
export function validateMultiAnalyte(df: DG.DataFrame, subjectName: string): ValidationError[] {
  const errs: ValidationError[] = [];
  const detection = detectAnalyteColumns(df);
  if (detection.multiAnalyte) {
    errs.push({
      code: 'MULTI_ANALYTE_NOT_SUPPORTED',
      message:
        `v1 supports single-analyte NCA only. Detected: ${detection.evidence.join('; ')}. ` +
        `Filter your CSV to one analyte (e.g., parent or M1) and re-upload.`,
      details: {evidence: detection.evidence},
    });
  }
  const subj = df.col(subjectName);
  if (subj === null) return errs;
  const subjectOf = (i: number): string | null => subj.isNone(i) ? null : subj.getString(i);
  for (const c of detectMixedCompoundGroups(df, subjectOf, [subjectName])) {
    if (detection.columns.includes(c.column)) continue;
    errs.push({
      code: 'MULTI_ANALYTE_NOT_SUPPORTED',
      message:
        `Subject '${c.group}' carries ${c.values.length} values of '${c.column}' ` +
        `(${c.values.join(', ')}) — a cassette (one subject dosed with several compounds). ` +
        'Profiles are built per subject, so these compounds would be pooled into one profile ' +
        `and every parameter would be wrong. Not supported yet: run one '${c.column}' value ` +
        'per analysis (filter the table), or give each compound its own subject ids.',
      column: c.column,
      details: {kind: 'cassette-shared-subject', column: c.column, subject: c.group, values: [...c.values]},
    });
  }
  return errs;
}

/**
 * Per-profile LLOQ sanity validator (AD-8). Runs PRE-extract — therefore
 * BLQ inference (`profile-inputs.ts` lines 122-133 'conc===0 inferred BLQ
 * when LLOQ column present') has NOT yet applied. "Non-BLQ" is defined by
 * the unified invariant: `NOT (blqFlag[i] is true/1) AND conc[i] > 0`. The
 * `conc > 0` conjunct applies in BOTH branches (with and without a blqFlag
 * column) and closes two semantically-distinct zero-encoding cases under a
 * single rule:
 *   1. **BLQ-coded-as-zero** (Watson / Thermo Fisher LIMS encode BLQ as
 *      literal zero with no flag column — DATA-GAP-13 closure).
 *   2. **Predose t=0=0 with explicit blqFlag=false** (composite / serial
 *      datasets — "no drug administered yet" is semantically distinct from
 *      "below detection"; the scientist would rightly object to flagging
 *      it BLQ. BUG-04 closure, Ed validation observation #8).
 *
 * Emits `LLOQ_OUT_OF_RANGE` per failing profile on any of:
 *  - **(a)** `LLOQ > min(non-BLQ conc) × 1.5` — LLOQ above the lowest
 *    measured value.
 *  - **(c)** `LLOQ > max(conc) × 10` — LLOQ wildly above all measurements.
 *    **Skipped for single-observation profiles** (`rows.length === 1`,
 *    sparse-destructive design) — a mouse sampled at the terminal tail
 *    with sub-LLOQ conc is normal toxicology TK, not a data error
 *    (BUG-03 / Ed validation #7).
 *  - **(d)** `LLOQ ≤ 0` — data error.
 *  - **(e)** all conc < LLOQ — entire profile BLQ (flag, don't silently
 *    zero-out). **Skipped for single-observation profiles** for the same
 *    sparse-destructive reason as (c). For 2+ row profiles (serial /
 *    composite), conditions (c) and (e) still surface real LLOQ / unit
 *    errors.
 *
 * Condition (b) `LLOQ < min(non-BLQ conc) / 1000` was dropped from v1 per
 * R2-finding-4 (high-sensitivity LC-MS/MS routinely exceeds the 1000 ratio
 * on sensitive compounds — uncited hard-block risk). Forward-deferred to
 * REGISTRY stream `unit-semantics-lloq-ratio-threshold-research`.
 *
 * The LLOQ scalar/column lives in pre-conversion units (same as
 * concentration); per AD-7 the same-factor rescaling invariant ensures the
 * comparison is unit-consistent post-rescale, but for the validator we work
 * in pre-conversion space — only the relative magnitudes matter, and they
 * are invariant under the same-factor rescaling.
 */
export function validateLloqSanity(
  df: DG.DataFrame,
  contract: ColumnContract,
  scalars: ProfileScalars,
  indices: ReadonlyArray<ProfileIndex>,
): ValidationError[] {
  const errs: ValidationError[] = [];
  const obs = contract.observations;
  const concCol = df.col(obs.concentration.name);
  if (concCol === null) return errs;
  const concRaw = concCol.getRawData();

  const lloqCol = obs.lloq !== undefined ? df.col(obs.lloq.name) : null;
  const lloqRaw = lloqCol !== null ? lloqCol.getRawData() : null;
  const lloqScalar = scalars.lloq;

  const blqFlagCol = obs.blqFlag !== undefined ? df.col(obs.blqFlag.name) : null;

  for (const index of indices) {
    const rows = index.observationRows;
    if (rows.length === 0) continue;

    // Gather conc values for this profile, FLOAT_NULL → skip.
    const concVals: number[] = [];
    const nonBlqConc: number[] = [];
    let i = 0;
    for (; i < rows.length; i++) {
      const r = rows[i];
      const c = concRaw[r];
      if (c === DGNS.FLOAT_NULL) continue;
      if (Number.isNaN(c)) continue;
      concVals.push(c);
      // Non-BLQ classification — unified invariant `NOT(flagged BLQ) AND
      // conc > 0`. See method JSDoc for the two zero-encoding cases this
      // closes (BLQ-as-zero / predose t=0=0). The `conc > 0` guard applies
      // in both branches.
      let isBlq: boolean;
      if (blqFlagCol !== null) {
        // BOOL storage is per-cell; use `get` to read true/false → 0/1.
        const v = blqFlagCol.get(r);
        const flaggedBlq = v === true || v === 1;
        isBlq = flaggedBlq || !(c > 0); // BUG-04: predose t=0=0 not non-BLQ
      } else
        isBlq = !(c > 0); // DATA-GAP-13: BLQ-coded-as-zero (no flag column)

      if (!isBlq) nonBlqConc.push(c);
    }

    if (concVals.length === 0) continue;

    const lloqForProfile: number = lloqRaw !== null ?
      // Per-row LLOQ — use the first row's value as the profile-level
      // sanity reference. (UC-01 lloq is conventionally constant per
      // profile / per assay run; using row 0 is the cheapest defensible
      // sample. Per-row variation, if it ever occurs, is a different
      // validator class — out of v1 scope.)
      (lloqRaw[rows[0]] === DGNS.FLOAT_NULL ? lloqScalar : (lloqRaw[rows[0]] as number)) :
      lloqScalar;

    // No LLOQ specified — skip sanity. The form's "0 = unset" convention
    // (see ProfileScalars.lloq JSDoc) means a 0 scalar with no LLOQ column
    // mapped indicates the user did not provide an LLOQ. Skipping is the
    // correct disambiguation: condition (d) "LLOQ ≤ 0" must not fire on
    // unset, only on user-provided ≤0 — but since `lloq=0` cannot reach
    // sci-comp's BLQ comparisons (BLQ inference is gated on lloqCol, not
    // scalar; see profile-inputs.ts:145-159), the science risk is bounded:
    // no path exists where a scalar `lloq=0` silently affects AUC/CL.
    if (lloqForProfile === 0 && lloqCol === null) continue;

    // (d) LLOQ ≤ 0
    if (lloqForProfile <= 0) {
      errs.push({
        code: 'LLOQ_OUT_OF_RANGE',
        message:
          `Subject '${index.subjectKey}': LLOQ = ${lloqForProfile} is ≤ 0 ` +
          `(condition d). Enter a positive LLOQ in the mapping pane.`,
        details: {
          kind: 'non-positive',
          subject: index.subjectKey,
          lloq: lloqForProfile,
        },
      });
      continue;
    }

    // Conditions (c) and (e) are sampling-design-aware: skipped for
    // single-observation (sparse-destructive) profiles where an all-BLQ
    // / wildly-off-LLOQ shape is a legitimate workflow outcome — a mouse
    // sampled at the terminal tail with conc < LLOQ is normal toxicology
    // TK, not a data error (Ed validation #7 / BUG-03). The compute layer
    // already handles such profiles per-subject (`status: 'failed'`); the
    // hard-block here would kill an entire batch on one legitimate row.
    // Condition (d) "LLOQ ≤ 0" still applies (data error regardless of
    // sampling design); condition (a) is self-protected (sparse all-BLQ
    // has `nonBlqConc.length === 0` → its check is skipped naturally).
    // For serial / composite profiles (`rows.length > 1`), (c) and (e)
    // still surface real LLOQ / unit errors.
    if (rows.length > 1) {
      let maxC = -Infinity;
      for (let j = 0; j < concVals.length; j++)
        if (concVals[j] > maxC) maxC = concVals[j];

      // (c) LLOQ > max(conc) × 10
      if (lloqForProfile > maxC * 10) {
        errs.push({
          code: 'LLOQ_OUT_OF_RANGE',
          message:
            `Subject '${index.subjectKey}': LLOQ = ${lloqForProfile} exceeds ` +
            `max(conc) × 10 = ${maxC * 10} (condition c). LLOQ is wildly above ` +
            `all measurements — verify the LLOQ value and concentration unit.`,
          details: {
            kind: 'above-max-times-ten',
            subject: index.subjectKey,
            lloq: lloqForProfile,
            maxConc: maxC,
          },
        });
        continue;
      }

      // (e) all conc < LLOQ
      let allBelow = true;
      for (let j = 0; j < concVals.length; j++)
        if (concVals[j] >= lloqForProfile) {allBelow = false; break;}

      if (allBelow) {
        errs.push({
          code: 'LLOQ_OUT_OF_RANGE',
          message:
            `Subject '${index.subjectKey}': all ${concVals.length} measured ` +
            `concentrations are < LLOQ ${lloqForProfile} (condition e). The ` +
            `entire profile is below limit-of-quantification — verify the LLOQ ` +
            `value or check whether the wrong subject was uploaded.`,
          details: {
            kind: 'all-below-lloq',
            subject: index.subjectKey,
            lloq: lloqForProfile,
          },
        });
        continue;
      }
    }

    // (a) LLOQ > min(non-BLQ conc) × 1.5
    if (nonBlqConc.length > 0) {
      let minNonBlq = Infinity;
      for (let j = 0; j < nonBlqConc.length; j++)
        if (nonBlqConc[j] < minNonBlq) minNonBlq = nonBlqConc[j];

      if (lloqForProfile > minNonBlq * 1.5) {
        errs.push({
          code: 'LLOQ_OUT_OF_RANGE',
          message:
            `Subject '${index.subjectKey}': LLOQ = ${lloqForProfile} exceeds ` +
            `min(non-BLQ conc) × 1.5 = ${minNonBlq * 1.5} (condition a). LLOQ ` +
            `must be at or below the lowest measured non-BLQ value — verify ` +
            `the LLOQ value or check whether BLQ rows are correctly flagged.`,
          details: {
            kind: 'above-min-non-blq',
            subject: index.subjectKey,
            lloq: lloqForProfile,
            minNonBlqConc: minNonBlq,
          },
        });
      }
    }
  }

  return errs;
}

/**
 * Sparse-design topology cross-check (UC-04 / AD-2, AD-11). Mapping-time early
 * feedback complementing the authoritative compute-time check in sci-comp
 * `sparseAuc` (which derives the topology from the full r_ij matrix).
 *
 * Emits **warning-severity** `SPARSE_TOPOLOGY_MISMATCH` (non-blocking — the
 * caller adds it to the diagnostics stream, NOT the blocking `ValidationError[]`
 * gate, per the route-code precedent) when the declared sampling type is
 * `destructive` but the animal × nominal-time table shows an animal sampled at
 * ≥ 2 distinct nominal times (i.e. the design is actually batch/serial). A wrong
 * "destructive" label would otherwise yield a too-narrow Bailer SE on paired
 * data (FM-1).
 *
 * No-op when the sparse roles (`animalId`, `nominalTime`, `samplingType`) are not
 * all mapped — a dense analysis or an under-specified sparse mapping is silent
 * here (the compute layer still guards via `SPARSE_ANIMAL_ID_ABSENT_…`).
 *
 * Cheap and focused: a single pass building `animal → set<nominal time>`; it does
 * NOT recompute the Holder covariance (that is `sparseAuc`'s job — rule 5).
 */
export function validateSparseDesign(
  df: DG.DataFrame,
  contract: ColumnContract,
): ValidationError[] {
  const obs = contract.observations;
  if (obs.animalId === undefined || obs.nominalTime === undefined || obs.samplingType === undefined)
    return [];
  const animalCol = df.col(obs.animalId.name);
  const nomCol = df.col(obs.nominalTime.name);
  const typeCol = df.col(obs.samplingType.name);
  if (animalCol === null || nomCol === null || typeCol === null) return [];

  const nomRaw = nomCol.getRawData();
  const n = df.rowCount;

  // Declared destructive iff the dominant label is "destructive".
  let destructiveRows = 0;
  let labeledRows = 0;
  for (let i = 0; i < n; i++) {
    const label = typeCol.isNone(i) ? '' : typeCol.getString(i).trim().toLowerCase();
    if (label === '') continue;
    labeledRows++;
    if (label === 'destructive') destructiveRows++;
  }
  if (labeledRows === 0 || destructiveRows !== labeledRows) return [];

  // Per-animal distinct nominal times.
  const byAnimal = new Map<string, Set<number>>();
  for (let i = 0; i < n; i++) {
    const a = stringifyKey(animalCol.get(i));
    const t = nomRaw[i];
    if (t === DGNS.FLOAT_NULL || Number.isNaN(t)) continue;
    let set = byAnimal.get(a);
    if (set === undefined) {set = new Set<number>(); byAnimal.set(a, set);}
    set.add(t);
  }
  let spanning = 0;
  for (const set of byAnimal.values()) if (set.size >= 2) spanning++;
  if (spanning === 0) return [];

  return [{
    code: 'SPARSE_TOPOLOGY_MISMATCH',
    message:
      `Sampling type is declared 'destructive', but ${spanning} animal(s) are ` +
      `sampled at ≥ 2 nominal times — the design is actually batch/serial. The ` +
      `data-derived topology will be used (the Holder covariance SE, not a ` +
      `too-narrow Bailer SE). Verify the SamplingType / animal-ID mapping.`,
    column: obs.samplingType.name,
    details: {kind: 'declared-destructive-but-paired', spanningAnimals: spanning},
  }];
}

function scanDuplicateSubjectTime(
  df: DG.DataFrame,
  subjectCol: string,
  timeCol: string,
): ValidationError[] {
  const errs: ValidationError[] = [];
  const subj = df.col(subjectCol);
  const tcol = df.col(timeCol);
  if (subj === null || tcol === null) return errs;
  const timeRaw = tcol.getRawData();
  const n = df.rowCount;
  const seen = new Map<string, Map<number, number>>();
  for (let i = 0; i < n; i++) {
    const sk = stringifyKey(subj.get(i));
    const t = timeRaw[i];
    let bucket = seen.get(sk);
    if (bucket === undefined) {
      bucket = new Map<number, number>();
      seen.set(sk, bucket);
    }
    const prev = bucket.get(t);
    if (prev !== undefined) {
      errs.push({
        code: 'DUPLICATE_SUBJECT_TIME',
        message: `Duplicate (subject='${sk}', time=${t}) at rows ${prev} and ${i}.`,
        column: timeCol,
        row: i,
      });
    } else
      bucket.set(t, i);
  }
  return errs;
}

