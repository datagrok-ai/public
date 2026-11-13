/**
 * ColumnContract — type-only description of the columns the adapter expects in
 * a UC-01 observations DataFrame, plus the structured validation error shape
 * the (Phase 3) validators will emit.
 *
 * Scope (UC-01, single-compound, single-occasion, single-analyte): the
 * observations table is the only required DataFrame. Dose and route may be
 * per-row columns OR scalars entered in the column-mapping dialog (Phase 5);
 * the contract describes only the per-row column case. Body weight is
 * per-subject (constant within a subject) when present; the validator will
 * enforce constancy in Phase 3.
 *
 * Deferred (per OQ-2 in TODO.md and PRD §10): separate dosing DataFrame,
 * multi-analyte profile keys, multi-occasion (Day 1 / Day N) routing. The
 * interface §2 sketch shows the long-term multi-table shape.
 *
 * No runtime validators here; this module is types only. Hand-written
 * validators land in `validators.ts` (Phase 3).
 */

/**
 * Raw runtime type of a column. Used by the (Phase 3) validators to check the
 * DataFrame column matches the contract, and by the (Phase 3) extractors to
 * pick the right typed-array constructor when copying into ProfileInputs.
 *
 * - `string`  — subjectId (subjects may be string identifiers)
 * - `float64` — time, concentration, dose, body weight, LLOQ
 * - `int32`   — route code (RouteCode encoded as 0..5; see sci-comp `nca.RouteCode`)
 * - `uint8`   — BLQ flag (0 measurable, 1 below LOQ)
 *
 * Interface v2.1 §2 also lists `'float32' | 'uint32'`; UC-01 doesn't need
 * those, so they're omitted here. The set extends naturally if a later phase
 * needs them.
 */
export type ColumnRawType = 'string' | 'float64' | 'int32' | 'uint8';

/**
 * Raw-array category aliases — mirror interface §2 nomenclature so the
 * `ColumnSpec<TRaw>` generic constraint reads identically to the interface
 * sketch. UC-01 narrows each alias to the single typed-array variant the
 * spec needs; later phases may widen (e.g. add `Float32Array` to
 * `NumericRawArray`) without breaking the constraint shape.
 */
export type NumericRawArray = Float64Array;
export type IntegerRawArray = Int32Array;
export type BooleanRawArray = Uint8Array;
/** UC-01 addition: subject IDs may be strings (spec §7.1 row 1). */
export type StringRawArray = string[];

/**
 * One column's contract: name in the DataFrame, whether it's required, its
 * runtime raw type, and an optional conceptual unit string for semantic notes.
 *
 * The `TRaw` generic mirrors the interface §2 sketch — it type-tags a
 * `ColumnSpec` with its runtime raw-array category so a (Phase 3) typed
 * extraction API can return the correct `Float64Array` / `Int32Array` /
 * `Uint8Array` / `string[]` without a runtime cast on the caller side. The
 * generic is a phantom (no field of the interface references `TRaw`); the
 * `rawType` discriminator carries the runtime information.
 */
// `TRaw` is a phantom generic (no field references it) — kept to align with
// interface §2's typed-extraction sketch (Phase 2 design decision #1). The
// runtime extractor uses `rawType` as the discriminator; the type parameter
// tags ColumnSpec instances at the type level so a future typed-extractor
// API can return the matching typed-array variant without a runtime cast.
// eslint-disable-next-line @typescript-eslint/no-unused-vars
export interface ColumnSpec<TRaw extends NumericRawArray | IntegerRawArray | BooleanRawArray | StringRawArray> {
  /** Column name as it appears in the source DataFrame (post-mapping). */
  readonly name: string;
  /** True if the column must be present in the DataFrame for compute to proceed. */
  readonly required: boolean;
  /** Raw runtime type of the DataFrame column. */
  readonly rawType: ColumnRawType;
  /**
   * Optional conceptual unit string (e.g. `'h'`, `'ng/mL'`, `'mg/kg'`). The
   * adapter records the unit chosen by the user in the column-mapping dialog
   * (Phase 5) and propagates it into ProfileInputs.{timeUnits, ...}.
   */
  readonly units?: string;
}

/**
 * Schema for the observations DataFrame the UC-01 adapter consumes.
 *
 * Required:
 * - `subjectId` — string or integer identifier; rows with the same subjectId
 *   form one profile (Phase 3 profile-index builder).
 * - `time` — non-negative numeric, sorted ascending within each subject; no
 *   duplicate `(subject, time)` pairs.
 * - `concentration` — non-negative numeric or NaN; non-numeric values yield a
 *   `NON_NUMERIC_CONC` error.
 *
 * Optional (each is the per-row column variant; the column-mapping dialog
 * supports a scalar fallback for `dose`, `route`, `lloq`):
 * - `dose` — per-row administered dose; falls back to the dialog's scalar
 *   when absent.
 * - `route` — per-row route code (`RouteCode` 0..5); falls back to the
 *   dialog's scalar.
 * - `bodyWeight` — per-subject body weight in kg, used for weight-normalised
 *   CL/Vz. Validator (Phase 3) enforces constancy within a subject.
 * - `blqFlag` — 0/1 BLQ indicator; absent → all 0s (Phase 3 extractor).
 * - `lloq` — per-row lower limit of quantification; falls back to the
 *   dialog's scalar (default 0).
 */
export interface ColumnContract {
  readonly observations: {
    readonly subjectId: ColumnSpec<StringRawArray | IntegerRawArray>;
    readonly time: ColumnSpec<NumericRawArray>;
    readonly concentration: ColumnSpec<NumericRawArray>;
    readonly dose?: ColumnSpec<NumericRawArray>;
    readonly route?: ColumnSpec<IntegerRawArray>;
    readonly bodyWeight?: ColumnSpec<NumericRawArray>;
    readonly blqFlag?: ColumnSpec<BooleanRawArray>;
    readonly lloq?: ColumnSpec<NumericRawArray>;
    /**
     * Per-row infusion duration (h), used by the IV-infusion −T_inf/2 MRT/Vss
     * correction. Falls back to the dialog's `infusionDuration` scalar when
     * absent or null at a profile's first row. The event-driven (NONMEM)
     * transform emits this as the `~InfusionDuration` column so each subject's
     * AMT/RATE-derived duration reaches the pipeline independently — a single
     * global scalar cannot represent a multi-subject infusion dataset where
     * subjects have different RATE (Slice A2 / BD-2). Optional: a dense
     * single-route study omits it and is byte-unchanged.
     */
    readonly infusionDuration?: ColumnSpec<NumericRawArray>;

    // ── Sparse / destructive-sampling roles (UC-04, AD-11) ────────────────
    // The single-occasion contract above carries only an observed `time` and a
    // per-subject identifier. Sparse NCA needs three further roles, all
    // OPTIONAL (a dense UC-01/02/03 study omits them and is byte-unchanged):
    /**
     * Animal / subject-pool identifier for the r_ij covariance (AD-2). When
     * present, the overlap matrix is derived from the animal × nominal-time
     * table (deduplicated per dose group — IDs are commonly reused across
     * groups in batch designs). When **absent**, the design is treated as
     * destructive (r_ij = 0, Bailer) and `sparseAuc` emits
     * `SPARSE_ANIMAL_ID_ABSENT_ASSUMED_DESTRUCTIVE` — never a silent batch SE.
     * Distinct from `subjectId`: in a destructive design every observation is
     * its own subject, but a *batch* animal recurs across nominal times.
     */
    readonly animalId?: ColumnSpec<StringRawArray | IntegerRawArray>;
    /**
     * Nominal (scheduled) time, the binning key the composite means are built
     * on — distinct from the observed `time` (which drifts per animal). Sparse
     * pooling is by nominal time; dense AUC integration stays on observed time.
     */
    readonly nominalTime?: ColumnSpec<NumericRawArray>;
    /**
     * Declared sampling type (`destructive` / `batch` / `serial`). The *label*
     * only — `sparseAuc` derives the real topology from r_ij and cross-checks
     * it against this, warning `SPARSE_TOPOLOGY_MISMATCH` on disagreement (AD-2).
     */
    readonly samplingType?: ColumnSpec<StringRawArray>;
  };
}

/**
 * Structured validation error code emitted by the hand-written validators.
 *
 * Five UC-01 codes (PRD §3.1 row 5):
 *  `MISSING_COLUMN`, `WRONG_TYPE`, `NEGATIVE_TIME`, `NON_NUMERIC_CONC`,
 *  `DUPLICATE_SUBJECT_TIME`.
 *
 * Four codes added by `nca-unit-semantics` (E2 / AD-4):
 *  - `UNITS_MISMATCH` — canonical from sci-comp interface §2; v1 use: molar
 *    concentration (Path C hybrid with `details.kind = 'molar-not-supported'`).
 *    Architecture invariant #3 amended — the "out of UC-01 scope" phrasing
 *    is retired.
 *  - `MISSING_BODY_WEIGHT` — **deprecated, no longer emitted** (nca-cl-vz-convention
 *    P-A / AD-1). `/kg` dose without BW is now the valid per-kg path (CL in
 *    mL/min/kg, Vz in L/kg — nca-cl-bw-normalised-units) carrying a `PER_KG_NO_BW`
 *    warning instead of a block.
 *    The enum member is kept for schema/back-compat; do not remove.
 *  - `MULTI_ANALYTE_NOT_SUPPORTED` — multi-analyte CSV detected (e.g.
 *    `Parent_*` and `Metabolite_*` columns together, or a per-row
 *    `analyte`/`compound` categorical column with cardinality ≥ 2).
 *  - `LLOQ_OUT_OF_RANGE` — LLOQ sanity validator (AD-8 conditions a/c/d/e).
 *
 * Six codes added by `per-row-route` (AD-5). These register the full route
 * diagnostic vocabulary; their **severity is mixed** (see per-code notes).
 * Error-severity codes are emitted as blocking `ValidationError`s (the
 * controller returns `validation-failed`); warning-severity codes are NOT
 * blocking and ride the per-profile `provenance.warnings` stream
 * (`nca.ParameterWarning.code` is an open union, so the same code string
 * surfaces in `WARNINGS_JSON`) or the F% stratum table:
 *  - `UNMAPPABLE_ROUTE_COLUMN` — **error**. Column-quality pre-scan: < 50% of
 *    the route column's unique non-null values map to a known synonym (catches
 *    SEX mapped as route). Blocks.
 *  - `UNMAPPABLE_ROUTE` — **error** for an IV-adjacent unmappable string
 *    (O-6 mis-augmentation guard, blocks); **warning** for a clearly-non-IV
 *    unknown string (→ `other`, computes).
 *  - `ROUTE_SCALAR_COLUMN_CONFLICT` — **warning**. A mapped route column and a
 *    non-null scalar route disagree (defensive/programmatic path; the UI
 *    disables the scalar when a column is mapped, so it is not the steady
 *    state). Column wins.
 *  - `ROUTE_INCONSISTENT_WITHIN_SUBJECT` — **warning**. A subject's rows carry
 *    more than one route (G5); the profile resolves to the first row + warning.
 *  - `NO_IV_REFERENCE_FOR_STRATUM` — **warning**. F% requested for a stratum
 *    with no IV arm → `F_pct = NaN`.
 *  - `ROUTE_INFUSION_DURATION_MISSING` — **warning**. A profile resolves to
 *    `IV-infusion` but no positive infusion duration is present. Computes
 *    (sci-comp v0.6 applies extravascular augmentation and ignores duration);
 *    metadata-completeness flag.
 *
 * Two codes added by `uc-04-sparse-nca` (AD-2 / AD-11). Both are
 * **warning-severity** — they ride the non-blocking diagnostics stream (the
 * sparse pipeline surfaces them via `WARNINGS_JSON`; `validateSparseDesign`
 * returns them for early mapping-time feedback), never the blocking
 * `ValidationError[]` gate:
 *  - `SPARSE_TOPOLOGY_MISMATCH` — the declared sampling type disagrees with the
 *    animal × nominal-time table (e.g. label "destructive" but an animal spans
 *    ≥ 2 nominal times). The data-derived topology wins.
 *  - `SPARSE_ANIMAL_ID_ABSENT_ASSUMED_DESTRUCTIVE` — no animal-ID role mapped on
 *    a sparse analysis; the design is assumed destructive (Bailer). Emitted by
 *    `sparseAuc` (compute-time) and mirrored here for the shared code vocabulary.
 *
 * **Count grows 5 → 9 → 15 → 17** (5 preserved + 4 unit-semantics + 6 route + 2 sparse).
 */
export type ValidationCode =
  | 'MISSING_COLUMN'
  | 'WRONG_TYPE'
  | 'NEGATIVE_TIME'
  | 'NON_NUMERIC_CONC'
  | 'DUPLICATE_SUBJECT_TIME'
  | 'UNITS_MISMATCH'
  // deprecated — no longer emitted (nca-cl-vz-convention P-A / AD-1); kept for
  // schema/enum back-compat. See JSDoc above.
  | 'MISSING_BODY_WEIGHT'
  | 'MULTI_ANALYTE_NOT_SUPPORTED'
  | 'LLOQ_OUT_OF_RANGE'
  | 'UNMAPPABLE_ROUTE_COLUMN'
  | 'UNMAPPABLE_ROUTE'
  | 'ROUTE_SCALAR_COLUMN_CONFLICT'
  | 'ROUTE_INCONSISTENT_WITHIN_SUBJECT'
  | 'NO_IV_REFERENCE_FOR_STRATUM'
  | 'ROUTE_INFUSION_DURATION_MISSING'
  | 'SPARSE_TOPOLOGY_MISMATCH'
  | 'SPARSE_ANIMAL_ID_ABSENT_ASSUMED_DESTRUCTIVE'
  // Two codes added by `nca-nonmem-event-driven` (AD-7). Both **error-severity**
  // (blocking) — the event-driven transform rejects an SS profile it cannot honor
  // rather than emit a wrong AUCτ. Unvalidated by oracle (DATA-GAP-27):
  //  - `NONMEM_SS_UNSUPPORTED` — SS=1 without an inter-dose interval (II), or an
  //    unrecognised SS flag. Reject (no silent wrong steady-state AUCτ).
  //  - `NONMEM_SS_SUPERPOSITION_UNSUPPORTED` — SS=2 (steady state by
  //    superposition over a non-uniform prior schedule). Always rejected — the
  //    SS=1 honor path would compute a wrong AUCτ on a non-uniform-τ history.
  | 'NONMEM_SS_UNSUPPORTED'
  | 'NONMEM_SS_SUPERPOSITION_UNSUPPORTED';

/**
 * Structured validation error. Returned as `ValidationError[]` by the
 * validators; if non-empty, the controller blocks compute and the UI
 * renders the list in an error panel below the column-mapping form
 * (spec §6.1).
 *
 * Optional `column` / `row` locate the offending cell when applicable; e.g.
 * a `MISSING_COLUMN` error names the column but omits the row, while a
 * `NEGATIVE_TIME` error names both.
 *
 * Optional `details` is a Path C discriminator (AD-4 / E2): `UNITS_MISMATCH`
 * with `details.kind = 'molar-not-supported'` distinguishes the v1 molar
 * case from future per-row / mixed-cohort unit mismatches. `unknown` is
 * deliberately loose to avoid locking the discriminator shape ahead of F5
 * (mixed-unit cohort detection, forward-deferred per E1).
 */
export interface ValidationError {
  readonly code: ValidationCode;
  readonly message: string;
  readonly column?: string;
  readonly row?: number;
  readonly details?: Readonly<Record<string, unknown>>;
}
