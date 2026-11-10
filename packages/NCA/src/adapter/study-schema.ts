/**
 * NCA Project entity schema (v1) — the serialized form of a named analysis
 * workspace, plus the forward-migration registry and the source-signature
 * primitive used for drift detection.
 *
 * **Scope.** This module is pure data + pure functions — no DG view state, no
 * persistence I/O (that is `study-coordinator.ts`), no compute (that is
 * `runNcaPipeline`). It owns three things:
 *
 * 1. The TYPE of a persisted project (`NcaStudy`) and its parts, including the
 *    EXHAUSTIVE {@link AnalysisRecipe} — every mapping-pane input that flows into
 *    a computed parameter. A
 *    slim `{mapping, rules}` recipe would silently revert the NONMEM / sparse
 *    mode toggles on replay and change the compute PATH (different Cmax / AUCinf /
 *    CL / Vz with no error) — so the recipe carries the mode toggles, the
 *    confirmed NONMEM route, the unit pickers, and the full `NcaRules`.
 *
 * 2. The {@link SourceSignature} type and the synchronous FNV-1a string hash. Persistence is
 *    `nca-repository.ts` (the EMS schema `nca`); the file format and its migration ladder are gone.
 */

import type * as DG from 'datagrok-api/dg';
import type {NcaRules, RouteCode} from './index';
import {INSTALLED_ENGINE_VERSION} from './engine-version';
import {RESULTS_BASIS} from './unit-conversion';

/* Field history of the retired `study.json` format (v1..v16). The file format, its
 *  `schemaVersion` and its migration ladder are gone (nca-ems-storage-migration: rows in the EMS
 *  schema `nca`, versioned by the platform; json columns carry `recipe_version`). The notes below
 *  stay because they explain why each field of the types exists.
 *
 *  v2 (study-unified-grid): the per-compound summary tier gained honest
 *  per-(compound × route × dose-norm) shape — `route`/`doseNorm`/`pctExtrap` on
 *  {@link CompoundTierRow}, per-param `n` on {@link MeanCv} — and the recipe
 *  gained a persisted {@link AnalysisRecipe.groupByCol}.
 *
 *  v3 (chemist-sar-grid): {@link CompoundTierRow} gained `mrt` + `vss` (MeanCv) —
 *  the chemist SAR lead trio (AUC∞ / MRT / %F) needs MRT in the cross-compound
 *  tier, and Vss is the IV-relevant steady-state volume. Pure projection
 *  additions (the aggregator already computes both); the v2→v3 migration
 *  back-fills the two new fields with the "no estimate" sentinel and PRESERVES
 *  every existing row (no data loss on read). *
 *  v4 (chemist-sar-grid S2): {@link CompoundTierRow} gained `fPct` (MeanCv | null)
 *  — the third leg of the chemist SAR lead trio (AUC∞ / MRT / %F). %F is a
 *  CROSS-profile artifact (EV dose-normalized AUC∞ ÷ the sex-matched IV
 *  reference), so — unlike the within-group aggregates — it cannot be recovered
 *  from the cached tier alone and MUST be persisted (the study grid sees only the
 *  cache, not the results frame). `null` on IV rows / where no IV reference
 *  exists; the v3→v4 migration back-fills `fPct: null` on every existing row and
 *  PRESERVES all other fields (same safe pattern as v2→v3 — no clear, no data
 *  loss). *
 *  v5 (chemist-sar-grid S4): {@link CompoundTierRow} gained `curve`
 *  ({@link CompoundCurve} | null) — the inline per-compound semi-log geometric-mean
 *  thumbnail (one polyline series per dose group, BUG-08 no cross-dose pooling) —
 *  and `tmax` (MeanCv) for the row-click drill-in scalar card. The study overview
 *  is cache-only (zero recompute at render), so the thumbnail's curve points must
 *  be PERSISTED (they cannot be re-derived from the cached tier — no observations
 *  at render time). Both are built at compute time in `buildSummaryCache` (curve
 *  reuses `mean-sd-profile.ts`; `tmax` is a pure aggregate carry-through — `TMAX`
 *  is already in `DEFAULT_AGG_PARAMS`). Like v2→v3 / v3→v4, the v4→v5 migration
 *  BACK-FILLS the two new fields (`curve: null`, `tmax` = the "no estimate"
 *  sentinel) and PRESERVES every existing row — no data loss. *
 *  v6 (chemist-sar-grid S4 polish): {@link CompoundTierRow.warnings} changed from
 *  `string[]` (bare codes) to {@link CompoundWarning}`[]` (`{code, message}`) so the
 *  study grid's Flags tooltip can show sci-comp's own human-readable `message`
 *  instead of the raw code — the message was previously dropped at cache-build time.
 *  The v5→v6 migration maps each old code `c` to `{code: c, message: c}` (the code is
 *  the only text available for already-cached rows; a recompute repopulates the real
 *  message). Every other field is preserved — no data loss. *
 *  v7 (group-col-stale-on-rename): {@link CompoundTierRow} gained `isWholeAnalysis`
 *  (boolean) — the durable, label-INDEPENDENT sentinel for the study grid's "whole
 *  analysis → —" case. It replaces a `compoundId === AnalysisEntry.label` string
 *  compare that broke on rename (the label mutates, the compute-time-frozen
 *  `compoundId` does not, so the stale pre-rename label leaked into the Group column).
 *  Set at compute time in `buildSummaryCache` (the exact negation of the compoundId
 *  real-value branch). The v6→v7 migration BACK-FILLS it from the immutable persisted
 *  `groupByCol` (null / `'ROUTE'` / `'DOSE_NORM'` → whole-analysis), NOT the label,
 *  and PRESERVES every existing row — no data loss. *
 *  v8 (multi-column-group-by): the persisted grouping key widened from a scalar
 *  {@link AnalysisRecipe.groupByCol} (`string | null`) to a set
 *  {@link AnalysisRecipe.groupByCols} (`readonly string[]`) so a scientist can group
 *  by a composite `(Sex, Dose, …)` tuple, and {@link CompoundTierRow} gained
 *  `groupValues` ({@link GroupKeyValue}`[]`) — the per-key `{key, value}` projection
 *  the study grid renders as one column per selected key (Variant-2, user decision
 *  2026-07-09; per-dimension filter/sort). The v7→v8 migration is a LOSSLESS
 *  singleton widen: `groupByCol` `col` → `groupByCols` `[col]`, `null` → `[]`, and
 *  each tier row's `groupValues` back-fills from the immutable scalar `compoundId`
 *  (`[{key: <groupByCol>, value: <compoundId>}]` for a real group, `[]` for a
 *  whole-analysis row). A singleton array yields an identical grouping, so the
 *  cached `compoundTier` stays valid and needs no re-projection. *
 *  v9 (study-grid-grouping-columns): {@link CompoundTierRow} gained `keyInventory`
 *  ({@link KeyInventoryEntry}`[]`) — the per-tier-row, subpopulation-scoped inventory of
 *  a candidate grouping column's distinct values WITH per-value subject counts. It lets
 *  the study grid render a key an analysis has in its data but did NOT group by as its
 *  counted composition (`F (15); M (15)`) instead of an ambiguous `—` (D1). Scoped to the
 *  ROW's subpopulation (not the analysis) so a route/dose-norm-split whole-analysis case
 *  discloses each row's OWN composition. **The migration ALSO nulls a stale multi-series
 *  `curve`** (`series.length > 1`): the study grid's inline thumbnail no longer auto-splits
 *  by dose (E3 pools per real group + discloses a cross-dose pool), so a pre-fix dose-split
 *  curve is semantically stale and cannot be re-derived from the scalar cache → nulled to
 *  blank until recompute (the v4→v5 baked-curve-null precedent). A `series.length ≤ 1` curve
 *  is ALREADY the pooled single series E3 produces, so it is PRESERVED — but its stale
 *  `series[0].group` (a dose VALUE) is REWRITTEN to the row's group label (its `groupValues`
 *  joined, or `''` for a whole-analysis row) so a preserved curve matches what a fresh
 *  recompute would label. The v8→v9 migration back-fills `keyInventory: []` on every row and
 *  PRESERVES all other fields. *
 *  v10 (study-grid-passthrough-columns): user-selected source columns (SMILES → structures,
 *  chemical name, series, project, …) are carried into the cache-only study grid as passthrough
 *  (annotation) columns. Two additive shape changes: {@link AnalysisRecipe} gained
 *  `passthroughCols` ({@link PassthroughColumnSpec}`[]` — `{name, isStructure}`, the user's
 *  dialog pick; `isStructure` captured ONCE from the source frame's live `semType` and PRESERVED
 *  across recompute, never re-derived), and {@link CompoundTierRow} gained `passthrough`
 *  ({@link PassthroughEntry}`[]` — per-subpopulation distinct values + subject counts + the
 *  structure flag, the cache-only render substrate). Unlike E4's `keyInventory`, the passthrough
 *  candidate list is CAP-FREE (a 25-compound cassette's SMILES must survive) and always shown.
 *  The v9→v10 migration back-fills `passthrough: []` on every tier row AND `passthroughCols: []`
 *  on every analysis, PRESERVES all other fields — purely additive (the v7→v8 pattern), so an
 *  empty set reproduces pre-v10 behavior exactly (zero passthrough columns). *
 *  v11 (context-panel-ct-profile): {@link CurveSeries} gained the OPTIONAL parallel band arrays
 *  `geoLower?` / `geoUpper?` — the ×/÷ gSD geometric band already computed by
 *  `buildMeanSdProfile` (`mean-sd-profile.ts:145-148`) but previously DISCARDED by
 *  `buildCurvesByCompound`. Persisting them lets the study-overview context panel draw the band
 *  CACHE-ONLY (no re-derivation) in both single-row and multi-row overlay (BQ-1 band-always).
 *  The fields are OPTIONAL, not back-filled: their ABSENCE on a legacy curve is the honest
 *  "band not yet computed" signal (the panel discloses "recompute to show"), distinct from a
 *  PRESENT band whose element is `FLOAT_NULL` at a point (genuine n=1 / no estimable spread —
 *  honest no-whisker, no disclosure). A required-field + universal-`FLOAT_NULL` back-fill would
 *  make those two cases visually identical (both plain-line) → thesis-#4 false confidence, so
 *  `MIGRATIONS[10]` does NOT back-fill — it only bumps `schemaVersion`. A freshly-built curve
 *  ALWAYS carries both arrays (band computed). Purely additive; every existing field preserved
 *  (the v9→v10 lossless-additive pattern). *
 *  v12 (group-by-duplicate-restore, BUG-11): REPAIR migration — dedupes the two persisted
 *  user-multi-selected SETS, {@link AnalysisRecipe.groupByCols} and
 *  {@link AnalysisRecipe.passthroughCols} (by spec `name`), to the distinct entries they always
 *  meant. Both come from the same stock columns picker, which can commit a name twice. The
 *  `passthroughCols` half is the more dangerous one: its counts are keyed by column name, so a
 *  repeated spec DOUBLES every subject count in the study-grid disclosure (a silent wrong
 *  number), where the `groupByCols` half only doubles labels. Studies
 *  saved before the BUG-11 fix carry a `groupByCols` the mapping pane grew by re-committing
 *  already-selected names (`['ROUTE','Dose_mg_kg','ROUTE','Sex','Dose_mg_kg']` for a 3-key
 *  pick; unbounded on repeat — see `dedupeGroupBySelection` in `nca-mapping-pane.ts`). A
 *  repeated key is a NO-OP partition (grouping by [A,B,A] ≡ [A,B]), so no cached NUMBER is
 *  wrong and this repairs LABELS + the replayed selection only. The cached `summaryCache` is
 *  left UNTOUCHED (unlike the v7→v8 re-projection): its `compoundId`/`curveKey` were built
 *  from the duplicated key list and are self-CONSISTENT (the same doubled string on both
 *  sides of every join), so rewriting `groupByCols` alone cannot desync them, while
 *  re-projecting them here would have to re-derive values this migration does not have. The
 *  next compute regenerates the tier cleanly from the repaired recipe; until then
 *  `dedupGroupValuesByCanonical` (BUG-10) keeps the doubled label off the screen. Preserves
 *  every existing field — no data loss. *
 *  v13 (lambdaz-adjr2-threshold-spec-v2 §7): {@link AnalysisRecipe} gained the OPTIONAL
 *  `acceptanceRSquared` (`R_acc`, the λz adj-R² acceptance threshold), persisted so save/load
 *  locks it with the analysis. `MIGRATIONS[12]` BACK-FILLS `0.9` (the introduction-time default;
 *  a migration is a historical snapshot, so it hardcodes the value rather than tracking the live
 *  `DEFAULT_ACCEPTANCE_RSQUARED`) and preserves every existing field. The field is OPTIONAL so a
 *  reader tolerates a file that lacks it (falls back to the default); `getCurrentRecipe()` always
 *  emits it, so every v13 SAVE carries it. Additive — no data loss (the v9→v10 pattern).
 *
 *  v14 (sci-comp-nca-blq-c0-gate F-N.6): {@link AnalysisEntry} gained `engineVersion` — the
 *  `@datagrok-libraries/sci-comp` version that computed the cached `summaryCache`, stamped on every
 *  compute writeback. The cache was invalidated ONLY on SOURCE drift; a numbers-changing engine
 *  release (0.11.0) would otherwise leave existing studies silently showing pre-fix numbers until
 *  each analysis happened to be reopened (thesis #4). A `null` ("computed by an unknown engine") is
 *  treated as stale by `isStale` (nca-repository.ts) — the "stale — refreshes on open" cue fires and
 *  the on-open recompute stamps the installed version. (The schema migrations and the coordinator's
 *  drift probe were removed by the EMS storage migration.)
 *
 *  v15 (nca-workflow-stages F-S1.7): the analysis entry becomes the durable record of the
 *  analyst's WORK, not only of the recipe and its cached numbers. {@link AnalysisEntry} gained
 *  `confirmed` (the recipe snapshot + fingerprint of the last `Confirm`, so a reopen can land on
 *  the right stage and `Revert changes` works after a reopen — AD-5), `lzOverrides` (the keyed λz
 *  point overrides, identity `(profileKey, timeKey, ordinal)` independent of any frame's row index,
 *  so derived-mode overrides survive a Compute and a reopen — D22 / AD-2), `reviewDecisions` (one
 *  ledger per review unit: the current decision + every superseded one, never obscured — D13,
 *  21 CFR 11.10(e)), `stageCache` (the derived stage for the hub; a CACHE, never truth — D1) and
 *  `origin` (set when a standalone analysis is moved into a study, for link redirect). `summaryCache`
 *  becomes NULLABLE: a never-confirmed analysis has none (a study entry never does — AC-S1.7.2).
 *  {@link MeanCv} gained the optional `nExcluded`. `MIGRATIONS[14]` back-fills empties / nulls;
 *  D23 means no legacy review semantics are designed, but the ladder stays versioned.
 *
 *  v16 (nca-pipeline-save-rework F2, review Q1 + 6d): `origin` became a discriminated union —
 *  `{kind: 'standalone'}` (unchanged) or `{kind: 'pipeline'}` (the run that wrote the entry:
 *  run id / time / user / engine, the mapping and rules it used, its partition identity and the
 *  hash of the recipe it wrote). {@link AnalysisEntry} gained the optional `pipelineFlags` (the
 *  pipeline's per-profile verdict and reasons, kept out of the source CSV) and {@link NcaStudy}
 *  the optional `pipelineKey` (owner + study key of a pipeline-created study, checked against its
 *  derived id). `MIGRATIONS[15]` only bumps the version (the `MIGRATIONS[10]` precedent): every
 *  new field is optional with an "absent ≡ empty" reading, and every existing `origin` is already
 *  `null` or a standalone origin. The bump itself is load-bearing: a v15 reader would read a
 *  pipeline origin as a standalone one, and `SchemaTooNewError` made it refuse instead. */

// ─── Recipe (AD-5, exhaustive) ──────────────────────────────────────────────

/** How an optional scalar-or-column input is sourced. Mirrors the mapping-pane
 *  `<role>SourceInput` choice (`'value'|'column'`) plus `'none'` for the roles
 *  that can be omitted entirely (infusion / body weight / LLOQ). */
export type SourceMode = 'none' | 'value' | 'column';

/**
 * Every column / scalar mapping-pane input that affects a computed parameter.
 * Field names mirror the pane's input variables so {@link AnalysisRecipe} replay
 * is a one-to-one assignment (`nca-mapping-pane.ts` AD-5 replay block).
 *
 * Columns are stored by NAME (string), resolved against the source DataFrame at
 * replay time — a name that no longer exists is dropped to null by the picker,
 * exactly as a stale autodetect would be.
 */
export interface RecipeColumns {
  readonly subject: string | null;
  readonly time: string | null;
  readonly concentration: string | null;

  readonly doseSource: SourceMode;
  readonly doseCol: string | null;
  readonly doseScalar: number;

  readonly routeSource: SourceMode;
  readonly routeScalar: string;
  readonly routeCol: string | null;

  readonly infusionSource: SourceMode;
  readonly infusionDurationCol: string | null;
  readonly infusionDurationScalar: number;

  readonly bodyWeightSource: SourceMode;
  readonly bodyWeightCol: string | null;
  readonly bodyWeightScalar: number;

  readonly blqFlagCol: string | null;

  readonly lloqSource: SourceMode;
  readonly lloqCol: string | null;
  readonly lloqScalar: number;

  readonly sexCol: string | null;
  /** Raw-value → M/F/U designation map for the sex column. */
  readonly sexValueMap: Readonly<Record<string, string>>;
}

/** Declared sampling-topology sentinel — `'(auto)'` is serialized VERBATIM,
 *  never as null (R1 SY-F3: `null` "no input" and `'(auto)'` "declared
 *  auto-detect" have different semantics; `null` is not in the input's items). */
export type SamplingType = '(auto)' | 'destructive' | 'batch' | 'serial';

/** Mode toggles that change the COMPUTE PATH (not just labeling). Reverting any
 *  of these on replay silently changes the parameters. */
export interface RecipeModes {
  /** Event-driven NONMEM ingestion (EVID/AMT/CMT/RATE/MDV …). */
  readonly nonmem: boolean;
  /**
   * Resolved NONMEM route when RATE/CMT are ambiguous (route-confirm dialog
   * result, `nca-mapping-pane.ts:439`). `null` when not confirmed / not NONMEM.
   * R1 SY-F1: determines CL-vs-CL/F labeling + IV t-corrections; must round-trip
   * so replay does not re-trigger the dialog (R2-F1 deserialization).
   */
  readonly confirmedRoute: RouteCode | null;
  readonly sparse: {
    readonly on: boolean;
    readonly nominalTimeCol: string | null;
    readonly animalIdCol: string | null;
    readonly samplingType: SamplingType;
  };
}

/** One user-selected passthrough (annotation) column carried into the study grid
 *  (v10, study-grid-passthrough-columns). `name` is the source-column name (also the
 *  first-seen friendly header, before the grid's typographic merge); `isStructure`
 *  is the persisted structure designation — captured ONCE at the add-to-study dialog
 *  from the source DataFrame's live `semType` (`=== DG.SEMTYPE.MOLECULE`) and PRESERVED
 *  across recompute (LBC-2). It MUST NOT be re-derived from the results frame:
 *  `buildPerSubjectColumn` (`grouping-columns.ts`) rebuilds each carried column via
 *  `Column.fromFloat64Array`/`fromList`, which copies NEITHER tags NOR `semType`, so a
 *  SMILES column reaching the results frame has lost its MOLECULE semType — re-detecting
 *  there would silently always yield `false` and defeat structure rendering. Keeping
 *  `{name, isStructure}` as one unit (not two parallel `string[]` lists) avoids drift. */
export interface PassthroughColumnSpec {
  readonly name: string;
  readonly isStructure: boolean;
}

/** Source units for the four unit pickers — rescaled to canonical h / mg/L /
 *  mg / kg before compute. */
export interface RecipeUnits {
  readonly time: string;
  readonly concentration: string;
  readonly dose: string;
  readonly bodyWeight: string;
}

/**
 * The full replayable recipe — the lossless seam between the mapping pane
 * (`getCurrentRecipe()` emits) and lazy hydration (`mountAnalysisIntoView`
 * replays). `ncaRules` is typed `NcaRules` (R2-F5 — NOT `unknown`/plain object,
 * so the replay boundary is compiler-checked) and is the output of the pane's
 * `buildRules()` closure, not a parallel serializer.
 */
export interface AnalysisRecipe {
  readonly columns: RecipeColumns;
  readonly modes: RecipeModes;
  readonly ncaRules: NcaRules;
  readonly units: RecipeUnits;
  /**
   * The results-frame column(s) the analysis-view "Group by" ribbon selected
   * (v2 scalar `groupByCol`; widened to a set at v8, multi-column-group-by).
   * Persisted so the cache-only study grid can regroup without recompute. `[]` =
   * no grouping (whole-analysis tier); one entry = single-column grouping; ≥2 =
   * composite `(Sex, Dose, …)` grouping keyed by the tuple. Order is the user's
   * pick order and is preserved (it drives the per-key column order in the study
   * grid + the composite key order). There is no auto-resolvable "compound column"
   * (`RecipeColumns` has no analyte role); the group-by choice is the single honest
   * source, so it rides the recipe.
   */
  readonly groupByCols: readonly string[];
  /**
   * User-selected source columns carried into the cache-only study grid as
   * passthrough (annotation) columns (v10, study-grid-passthrough-columns) — the
   * SAR-grid plumbing (SMILES → structures, chemical name, series, …). Chosen at
   * add-to-study time (a DIALOG-time pick, NOT a mapping-pane input), so
   * `getCurrentRecipe()` emits `[]` by default and `openAddToStudyDialog` overrides
   * it with the user's selection when constructing the {@link AnalysisEntry}. Each
   * entry carries the persisted structure designation ({@link PassthroughColumnSpec}),
   * captured once at the dialog and preserved across recompute. `[]` = no passthrough
   * columns (pre-v10 behavior exactly).
   */
  readonly passthroughCols: readonly PassthroughColumnSpec[];
  /**
   * The λz adj-R² ACCEPTANCE threshold (`R_acc`, lambdaz-adjr2-threshold-spec-v2 §3.2/§7),
   * persisted so save/load locks it with the analysis. It is NOT a member of `ncaRules`
   * (SAP/presentation, thesis #1 — sci-comp never reads it), so it rides the recipe as its
   * own scalar field rather than inside `NcaRules`.
   *
   * **Optional for reader-tolerance:** a pre-persistence recipe (schema < 13) has no such
   * field; the reader falls back to the nca-studio default (`applyRecipe` / the input-init
   * both `?? DEFAULT_ACCEPTANCE_RSQUARED`). `MIGRATIONS[12]` back-fills it so a migrated file
   * is complete; `getCurrentRecipe()` always emits it, so every v13+ save carries it.
   */
  readonly acceptanceRSquared?: number;
}

// ─── Source pointer (decision e — pointer, never embedded data) ──────────────

/** `stored` = a d42 snapshot in the analysis row's `source_file` (every persisted analysis —
 *  nca-ems-storage-migration SD-4); `tableId` = the in-session launch binding, never persisted;
 *  `platformPath` / `dbQuery` = reserved kinds with no writer today. */
export type SourceKind = 'stored' | 'platformPath' | 'tableId' | 'dbQuery';

export interface SourcePointer {
  readonly kind: SourceKind;
  /** Snapshot file path (`System:DomainFiles/nca/analysis/<uuid>/source.d42`) / platform file path /
   *  DG table id / DB query text. */
  readonly ref: string;
}

// ─── Source signature (decision g — drift detection) ─────────────────────────

export type SignatureKind = 'rowcount-mtime' | 'content-hash' | 'computed-at-only';

export interface SourceSignature {
  readonly kind: SignatureKind;
  readonly rowCount: number;
  /** Platform-file mtime (ISO) when kind === 'rowcount-mtime'. */
  readonly mtime: string | null;
  /** FNV-1a content hash (hex) when kind === 'content-hash'. */
  readonly hash: string | null;
  /** Compute timestamp (ISO) — always present; the only signal for
   *  'computed-at-only' (DB-query sources with no cheap probe). */
  readonly computedAt: string;
}

// ─── Summary cache (decision h — two-tier overview, summary-only) ────────────

/** One semi-log polyline for the per-compound curve (v5, chemist-sar-grid S4;
 *  regrouped E3, study-grid-grouping-columns 2026-07-13). Post-E3 a curve is a
 *  SINGLE POOLED series per (compound × route) with a `CROSS_DOSE_POOL` disclosure
 *  when it spans ≥2 dose levels — it does NOT split per dose group (the pre-E3
 *  "one series per dose group / never pool" model is retired; pooling is disclosed,
 *  not avoided by splitting). `group` is the row's group label (its `groupValues`
 *  joined by ` · `, or `''` for a whole-analysis row). `time`/`geoMean` are the
 *  paired points of the geometric-mean concentration-time curve (the SAME geomean
 *  the analysis's own mean±SD chart plots — `mean-sd-profile.ts` — just small),
 *  sorted ascending by time, all observed/nominal timepoints (NO lossy
 *  downsample: the honest-thumbnail contract). Only measurable points ride here
 *  (a geometric mean is undefined for ≤ 0, so all-BLQ timepoints are absent).
 *
 *  `geoLower?`/`geoUpper?` (v11, context-panel-ct-profile) are the OPTIONAL ×/÷ gSD
 *  band — same length/index as `time`/`geoMean` when present, each element `FLOAT_NULL`
 *  where the band is genuinely inestimable (n=1 / degenerate spread, matching
 *  `mean-sd-profile.ts` null handling). A freshly-built curve ALWAYS carries both
 *  (band computed); ABSENT arrays mark a legacy pre-v11 curve (band not-yet-computed
 *  — the panel discloses "recompute to show"). The thumbnail ignores these columns;
 *  the context panel draws them as whiskers, cache-only. */
export interface CurveSeries {
  readonly group: string;
  readonly time: ReadonlyArray<number>;
  readonly geoMean: ReadonlyArray<number>;
  readonly geoLower?: ReadonlyArray<number>;
  readonly geoUpper?: ReadonlyArray<number>;
}

/** The persisted inline curve for one compound-tier row (v5, chemist-sar-grid S4;
 *  regrouped E3). A small, cache-only representation of the compound's semi-log
 *  geometric-mean profile — post-E3 a SINGLE POOLED {@link CurveSeries} per
 *  (compound × route), cross-dose pooling DISCLOSED via `CROSS_DOSE_POOL` (not
 *  avoided by per-dose splitting). It is PERSISTED because the study overview never
 *  recomputes at render (cache-only, AC-C3), so no observations exist at render time
 *  to derive it from. `null` / empty `series` = no measurable curve (no observations
 *  plumbed, or every point BLQ) — the thumbnail cell renders blank, never faked (rule 14). */
export interface CompoundCurve {
  readonly series: ReadonlyArray<CurveSeries>;
}

/** One per-case warning carried into the tier (v6, chemist-sar-grid S4 polish). The
 *  `code` is the stable identifier used for logic (e.g. NO_IV_REFERENCE detection);
 *  the `message` is sci-comp's own human-readable text, shown in the study grid's
 *  Flags tooltip. For a pre-v6 cached row the migration sets `message = code` (the
 *  only text then available); a recompute repopulates the real message. */
export interface CompoundWarning {
  readonly code: string;
  readonly message: string;
}

/** One selected grouping key's value for a compound-tier row (v8,
 *  multi-column-group-by). `key` is the friendly display label the study grid uses
 *  as the per-key column header (the results column's `friendlyName`, or the column
 *  name when it has none); `value` is this row's value for that key. Emitted in the
 *  user's pick order (matching {@link AnalysisRecipe.groupByCols}), EXCLUDING the
 *  implicit ROUTE / DOSE_NORM sub-keys (they ride their own Route / unit columns).
 *  A whole-analysis row carries an EMPTY `groupValues` (every per-key cell renders
 *  `—`). The drill-in reconstructs the row's composite tuple by zipping the
 *  analysis's real-compound `groupByCols` with these values in order. */
export interface GroupKeyValue {
  readonly key: string;
  readonly value: string;
}

/** One distinct value of a candidate grouping column WITH its subject count (v9,
 *  study-grid-grouping-columns). `count` is the number of SUBJECTS in the tier row's
 *  subpopulation carrying `value` — the honest-uncertainty signal (thesis #4): an
 *  `F (27); M (3)` imbalance must not read like `F (15); M (15)`. */
export interface KeyValueCount {
  readonly value: string;
  readonly count: number;
}

/** The per-tier-row inventory of ONE candidate grouping column an analysis has in its
 *  data but did NOT group by (v9, study-grid-grouping-columns). `key` is the friendly
 *  display label (the results column's `friendlyName`, or its name); the study grid
 *  normalizes it via `normalizeGroupKey` at merge time — the canonical id is NOT persisted
 *  (a single source of truth: `normalizeGroupKey` alone, no drift risk). `values` are this
 *  ROW's subpopulation's distinct non-null values + subject counts. The grid renders a
 *  present-but-ungrouped key's cell as `values` joined-with-counts; `—` only when the key
 *  is genuinely absent from the row's subpopulation (no inventory entry). Scoped
 *  per-tier-row (NOT per-analysis) so each subpopulation discloses its OWN composition. */
export interface KeyInventoryEntry {
  readonly key: string;
  readonly values: readonly KeyValueCount[];
}

/** The per-tier-row inventory of ONE user-selected passthrough (annotation) column (v10,
 *  study-grid-passthrough-columns). Modeled on {@link KeyInventoryEntry} but carries the
 *  persisted structure flag AND is CAP-FREE / always-shown (a passthrough column an analysis
 *  carries always renders, however many distinct values it holds — the E4 cardinality cap does
 *  NOT apply). `key` is the source column's friendly label (the grid normalizes it via
 *  `normalizeGroupKey` at merge time); `values` are this ROW's subpopulation's distinct non-null
 *  values + subject counts (reuses {@link KeyValueCount}). `isStructure` is copied verbatim from
 *  the recipe's {@link PassthroughColumnSpec} (single source of truth, never re-detected). The
 *  grid renders it BARE when `values.length === 1` (a lone SMILES → a 2D molecule), a counted
 *  join when many-text, and `—` + count when many-structure or absent (thesis #4 honesty). */
export interface PassthroughEntry {
  readonly key: string;
  readonly values: readonly KeyValueCount[];
  readonly isStructure: boolean;
}

export interface MeanCv {
  readonly mean: number;
  readonly cvPct: number;
  /** Per-parameter contributor count = the aggregator's `N_P` (number of
   *  subjects in the group that produced a fittable value for THIS parameter,
   *  na.rm). Distinct from {@link CompoundTierRow.n} (the group's subject
   *  count): a partial group has `n < CompoundTierRow.n`. Honest-uncertainty
   *  (thesis #4) — a mean of 2 is not a mean of 5 (v2, study-unified-grid). */
  readonly n: number;
  /** Profiles in the group whose λz-dependent value was withheld by an Exclude λz review
   *  decision (v15, nca-workflow-stages F-S5.2). `n` already excludes them; this lets the study
   *  grid say `n = a of b` and read `not reportable` when every contributor was excluded — a
   *  state distinct from "not estimated". Absent ≡ 0 (λz-independent params never carry it). */
  readonly nExcluded?: number;
}

/** The analysis-level line of a summary. Its compute time is not part of it: the stored analysis
 *  keeps `results_computed_at`, set only by a write that changes the results
 *  ({@link AnalysisEntry.resultsComputedAt}, ems-modification decision 2). */
export interface AnalysisTierSummary {
  readonly nCompounds: number;
  readonly nSubjects: number;
  readonly route: string;
  readonly dose: string;
  readonly status: 'ok' | 'warn' | 'error';
  readonly flagSummary: string;
}

export interface CompoundTierRow {
  readonly compoundId: string;
  /**
   * True when this case IS the whole analysis (no grouping, or the group-by is an
   * implicit route/dose-norm sub-key) — i.e. `compoundId` fell back to the analysis
   * label at compute time (v7, group-col-stale-on-rename). The study grid shows `—`
   * for these and the drill-in resolves them by this flag. It is the durable,
   * label-INDEPENDENT sentinel: an earlier design detected the whole-analysis case
   * by string-comparing `compoundId` against the LIVE `AnalysisEntry.label`, which
   * broke on rename (the label mutates, the frozen `compoundId` does not, so the
   * stale label leaked into the Group column). Never string-compare against the
   * label again — read this flag.
   */
  readonly isWholeAnalysis: boolean;
  /**
   * Per-selected-key values for this row (v8, multi-column-group-by) — the
   * display/persistence projection the study grid renders as one column per key
   * (Variant-2). One {@link GroupKeyValue} per user-selected grouping column, in
   * pick order, EXCLUDING the implicit ROUTE / DOSE_NORM sub-keys. EMPTY for a
   * whole-analysis row (`isWholeAnalysis === true`). The internal composite dedup /
   * curve key stays on `compoundId` (a `compositeGroupKey` U+0001 join); these are
   * the split-out values the grid + drill-in read per dimension.
   */
  readonly groupValues: readonly GroupKeyValue[];
  /**
   * Per-tier-row inventory of candidate grouping columns present in this row's
   * subpopulation but NOT grouped by (v9, study-grid-grouping-columns). One
   * {@link KeyInventoryEntry} per candidate column (Sex / Route / Dose / an unmapped
   * low-cardinality covariate), each holding that column's distinct values + SUBJECT
   * counts WITHIN this row's `groupByCols × route × dose-norm` partition. The study grid
   * uses it to render a present-but-ungrouped key as its counted composition
   * (`F (15); M (15)`) instead of an ambiguous `—` (D1, thesis #4). EMPTY when the row's
   * subpopulation carries no candidate column (or on a pre-v9 migrated row until recompute).
   * Scoped per-tier-row (NOT `SummaryCache`) so each subpopulation discloses its own
   * composition — a route-split whole-analysis case shows each route row's OWN counts.
   */
  readonly keyInventory: readonly KeyInventoryEntry[];
  /**
   * Per-tier-row inventory of the user-selected passthrough (annotation) columns this
   * analysis carries (v10, study-grid-passthrough-columns). One {@link PassthroughEntry}
   * per selected column present in this row's subpopulation, each holding that column's
   * distinct values + SUBJECT counts WITHIN this row's `groupByCols × route × dose-norm`
   * partition, plus the persisted structure flag. Unlike {@link keyInventory}, the
   * candidate list is CAP-FREE (the explicit user pick — a 25-distinct SMILES column is
   * kept) and ALWAYS shown (not only when ungrouped). EMPTY when the analysis selected no
   * passthrough columns (the common case) or on a pre-v10 migrated row until recompute.
   * Scoped per-tier-row (a coarse pooled row discloses ALL its distinct values with counts,
   * never one silently picked — thesis #4).
   */
  readonly passthrough: readonly PassthroughEntry[];
  /** Group subject count (the aggregator's `N` = group row count). */
  readonly n: number;
  /**
   * Route string of this route-homogeneous case (v2 — the results `ROUTE`
   * value, e.g. `'IV-bolus'`/`'PO'`). Every case is route-homogeneous by
   * construction (route is a group sub-key), so it 1:1-determines the CL/Vz
   * basis (`CL` vs `CL/F`) the study grid routes on at display time.
   */
  readonly route: string;
  /**
   * Dose-normalization basis of this unit-homogeneous case (v2 — the results
   * `DOSE_NORM` value, `'absolute'`/`'per_kg'`) — how the ENGINE dose was formed. It no
   * longer selects the CL/Vz unit: `cl`/`vz`/`vss` are per kg (mL/min/kg, L/kg) on every row
   * (nca-cl-bw-normalised-units AD-1/AD-8).
   */
  readonly doseNorm: string;
  readonly cmax: MeanCv;
  /** Time-to-Cmax mean/CV/n (v5, chemist-sar-grid S4) — a pure aggregate
   *  carry-through (`TMAX` is already in `DEFAULT_AGG_PARAMS`), persisted so the
   *  row-click drill-in scalar card can show Tmax from cache (no recompute). Not a
   *  study-grid column; card-only. */
  readonly tmax: MeanCv;
  readonly aucInf: MeanCv;
  readonly thalf: MeanCv;
  readonly cl: MeanCv;
  readonly vz: MeanCv;
  /** Mean residence time mean/CV/n (v3, chemist-sar-grid) — route-uniform
   *  (sci-comp emits one MRT column for all routes; `AUMC∞/AUC∞ − T_inf/2`). One
   *  of the chemist SAR lead trio (AUC∞ / MRT / %F). */
  readonly mrt: MeanCv;
  /** Volume of distribution at steady state mean/CV/n (v3, chemist-sar-grid) —
   *  IV-only (sci-comp gates `vss` to NaN for non-IV routes: an extravascular Vss
   *  would be Vss/F confounded by absorption). The study grid routes it like Vz
   *  (shown for intravascular cases; genuinely N/A — dropped — for EV). */
  readonly vss: MeanCv;
  /**
   * Absolute bioavailability %F mean/CV/n (v4, chemist-sar-grid S2) — the third
   * leg of the chemist SAR lead trio (AUC∞ / MRT / %F). CROSS-profile: the
   * mean/%CV/n over this compound's EXTRAVASCULAR strata F%, each F% = 100 ×
   * dose-normalized AUC∞(EV) ÷ the sex-matched IV reference
   * (`computeBioavailability`, sex-stratified denominator MANDATORY). `null` on
   * IV rows (F ≡ 1 by definition, never shown) AND where no paired IV+EV arm
   * exists at all; a `{mean: NaN, …}` sentinel where the compound is EV but has
   * no IV reference (NE — computed, not estimable; never faked, rule 14). The
   * study grid's all-empty guard omits the whole %F column when no case carries
   * a finite %F (Decision 2: no paired arm anywhere → no column). */
  readonly fPct: MeanCv | null;
  /** %AUCextrap mean/CV/n (v2) — the fragility flag (>20% = unreliable extrap)
   *  the study grid conditionally formats (design B4). */
  readonly pctExtrap: MeanCv;
  /** Inline per-compound curve (v5, chemist-sar-grid S4; regrouped E3) — the compound's
   *  semi-log geometric-mean profile, a SINGLE POOLED {@link CurveSeries} per (compound ×
   *  route) with a `CROSS_DOSE_POOL` disclosure when it spans ≥2 dose levels (post-E3: pooling
   *  is disclosed, NOT avoided by per-dose splitting). Carries the optional v11 ×/÷ gSD band
   *  (`geoLower`/`geoUpper`) so the context panel can draw whiskers cache-only. `null` when no
   *  observations were plumbed into the compute (e.g. a pre-v5 recompute, or a compute path
   *  with no dense frame) or the compound has no measurable points — the thumbnail cell then
   *  renders blank (never faked, rule 14). */
  readonly curve: CompoundCurve | null;
  /** Per-case warnings (v6: `{code, message}` — was bare `string[]` codes). The
   *  study grid shows each `message` on the Flags tooltip; `code` drives logic. */
  readonly warnings: ReadonlyArray<CompoundWarning>;
}

export interface SummaryCache {
  readonly analysisTier: AnalysisTierSummary;
  readonly compoundTier: ReadonlyArray<CompoundTierRow>;
}

// ─── Analysis + project ──────────────────────────────────────────────────────

/**
 * One analysis inside a project. Extends {@link AnalysisRecipe} so the recipe
 * fields (`columns`/`modes`/`ncaRules`/`units`) sit FLAT on the entry (matching
 * the AD-6 JSON shape) while `getCurrentRecipe()` can still return the
 * `AnalysisRecipe` subset as a compiler-checked unit.
 */
export interface AnalysisEntry extends AnalysisRecipe {
  readonly id: string;
  readonly label: string;
  readonly source: SourcePointer;
  readonly sourceSignature: SourceSignature;
  /** The cached two-tier summary, or `null` for an analysis that was never confirmed (v15 —
   *  a standalone analysis saved at Setup has no numbers). **A study entry always has a
   *  non-null cache** (AC-S1.7.2): Move to study is offered only at `stage ≠ setup`, and an
   *  un-confirmed edit keeps the last confirmed cache. Readers still narrow a null to the
   *  recompute placeholder rather than dereference. */
  readonly summaryCache: SummaryCache | null;
  /**
   * The `@datagrok-libraries/sci-comp` version that computed `summaryCache`
   * (v14), stamped on every compute writeback from `INSTALLED_ENGINE_VERSION`.
   * `null` = computed by an unknown engine (a pre-v14 file, back-filled by the
   * migration) — treated as STALE by the coordinator exactly like a stamp that
   * differs from the installed version: the cached numbers may predate a
   * numbers-changing release. OPTIONAL on the type so a test double that
   * constructs an entry keeps compiling (absent ≡ `null` ≡ unknown engine);
   * every path that WRITES an entry (`updateAnalysisSummary`,
   * `updateAnalysisEntry`, add-to-study) stamps it.
   */
  readonly engineVersion?: string | null;
  /**
   * The reporting-basis contract `summaryCache` was computed under (`RESULTS_BASIS`,
   * nca-cl-bw-normalised-units AD-8) — stamped ONLY where a cache is freshly computed and written;
   * every other write round-trips it. `null` / absent = written before the per-kg basis → STALE
   * (`isStale`), and the study grid blanks its CL / Vz / Vss cells.
   */
  readonly resultsBasis?: string | null;
  /** When the stored results last CHANGED (ems-modification decision 2) — read-only, set by the
   *  repository on read; a reopen that recomputes the same numbers does not move it. Absent / null
   *  = no results, or not read. */
  readonly resultsComputedAt?: string | null;
  /** Last `Confirm` (v15): the recipe snapshot, its recipe fingerprint, and when.
   *  `null` = never confirmed → the analysis is at Setup. Absent ≡ `null`. */
  readonly confirmed?: ConfirmedRecipe | null;
  /** Keyed λz point overrides (v15, AD-2). Absent ≡ `[]`. */
  readonly lzOverrides?: ReadonlyArray<LzOverride>;
  /** Review-decision ledgers, one per review unit that was ever decided (v15). Absent ≡ `[]`. */
  readonly reviewDecisions?: ReadonlyArray<ReviewLedger>;
  /** Derived stage cached for the hub (v15). Written only together with the state it derives
   *  from (AD-8). Read only through `formatStageBadge` + `isStale`. Absent ≡ `null`. */
  readonly stageCache?: StageCache | null;
  /** Where the entry came from (v15/v16): a standalone analysis moved into a study (F-S2.5 —
   *  the old `analyses/<label>` link redirects here), or a pipeline run that wrote it (v16).
   *  Absent / `null` ≡ made in the app. Readers narrow on `kind`. */
  readonly origin?: AnalysisOrigin | null;
  /** The pipeline's per-profile verdict and reasons from the run in `origin` (v16, review Q1).
   *  Kept here, not in the source CSV, so the CSV stays the delivered data. Absent ≡ `[]`. */
  readonly pipelineFlags?: ReadonlyArray<PipelineFlag>;
  /** The stored count of profiles with a pipeline-only reason (FIELD-11) — read-only, set by the
   *  repository on read, so a reader without the flags themselves (the study overview) can still
   *  show it. Absent → count {@link pipelineFlags}. */
  readonly pipelineFlagCount?: number;
}

/** One profile's pipeline verdict (v16). `profileKey` is the results-row `SUBJECT` for a dense
 *  profile and the app's composite unit key (`sparse:DOSE=<raw dose>`) for a sparse composite,
 *  so it joins the app's own review units. `reasons` routed the profile to `Review`; `info` is
 *  shown but routes nothing. */
export interface PipelineFlag {
  readonly profileKey: string;
  readonly verdict: 'Accepted' | 'Review';
  readonly reasons: ReadonlyArray<string>;
  readonly info: ReadonlyArray<string>;
  /** The subset of `reasons` the app's own λz acceptance check does NOT raise — a CRO incident
   *  comment, a screening heuristic, a unit-error heuristic. The app shows and counts only these
   *  (approved design 2026-10-15): the acceptance reasons are already in its ACCEPTANCE verdict
   *  and review queue, so repeating them would duplicate `N not reviewed`. */
  readonly pipelineOnly: ReadonlyArray<string>;
}

/** The {@link AnalysisRecipe} subset of an entry (the flat recipe fields). Moved from
 *  `views/nca-analysis-open.ts` (nca-pipeline-save-rework F1): the pipeline's "untouched" test
 *  hashes this projection, so it must be reachable without the view layer. */
export function recipeOf(entry: AnalysisEntry): AnalysisRecipe {
  return {
    columns: entry.columns,
    modes: entry.modes,
    ncaRules: entry.ncaRules,
    units: entry.units,
    groupByCols: entry.groupByCols,
    // v10: the persisted passthrough selection travels ON the recipe so a replay recompute
    // rebuilds the annotation columns from the authoritative dialog-time pick (LBC-2 —
    // never re-derived). getCurrentRecipe() emits [] (the pane has no passthrough input), so
    // the recompute call site (mountAnalysisIntoView) re-injects entry.passthroughCols.
    passthroughCols: entry.passthroughCols,
    // v13: the adj-R² acceptance threshold rides the recipe back into the pane on reopen, so
    // the input restores the SAVED value rather than the default. Like every field here, it
    // is a SEPARATE enumeration from getCurrentRecipe — omitting it reset R_acc on reopen
    // (user QA 2026-08-25). `undefined` on a pre-v13 entry → the pane's `?? DEFAULT` fallback.
    acceptanceRSquared: entry.acceptanceRSquared,
  };
}

/** What a NEW {@link AnalysisEntry} is built from ({@link makeAnalysisEntry}). */
export interface MakeAnalysisEntryArgs {
  readonly id: string;
  readonly label: string;
  readonly source: SourcePointer;
  readonly sourceSignature: SourceSignature;
  /** The full recipe — its `passthroughCols` IS the entry's annotation selection (one source). */
  readonly recipe: AnalysisRecipe;
  /** `null` = never confirmed (a Setup-stage standalone save, D1). */
  readonly summaryCache: SummaryCache | null;
  /** The v15 workflow state the entry starts with. */
  readonly workflow: {
    readonly confirmed: ConfirmedRecipe | null;
    readonly lzOverrides: ReadonlyArray<LzOverride>;
    readonly reviewDecisions: ReadonlyArray<ReviewLedger>;
    readonly stageCache: StageCache | null;
  };
  readonly origin: AnalysisOrigin | null;
  /** v16: the pipeline's per-profile flags (pipeline entries only). */
  readonly pipelineFlags?: ReadonlyArray<PipelineFlag>;
}

/**
 * The ONE builder of a new {@link AnalysisEntry} (nca-pipeline-save-rework F1, review 7a): Add to
 * study, the first standalone save, the pipeline save and the source-replaced save-beside all call
 * it, so a field added to the schema is set in one place. It SPREADS the recipe rather than listing
 * its fields — a per-field list is how `acceptanceRSquared` was once dropped on save (user QA
 * 2026-08-25). The engine stamp follows the cache: an entry with numbers is stamped with the
 * installed sci-comp version, one without (`summaryCache === null`) carries `null`.
 * `moveStandaloneIntoStudy` does not use it: a move COPIES an existing entry (`{...entry}`), so
 * every field travels by construction.
 */
export function makeAnalysisEntry(a: MakeAnalysisEntryArgs): AnalysisEntry {
  return {
    ...a.recipe,
    id: a.id,
    label: a.label,
    source: a.source,
    sourceSignature: a.sourceSignature,
    summaryCache: a.summaryCache,
    engineVersion: a.summaryCache === null ? null : INSTALLED_ENGINE_VERSION,
    resultsBasis: a.summaryCache === null ? null : RESULTS_BASIS,
    confirmed: a.workflow.confirmed,
    lzOverrides: a.workflow.lzOverrides,
    reviewDecisions: a.workflow.reviewDecisions,
    stageCache: a.workflow.stageCache,
    origin: a.origin,
    ...(a.pipelineFlags !== undefined ? {pipelineFlags: a.pipelineFlags} : {}),
  };
}

// ─── Workflow stage + review (v15, nca-workflow-stages) ──────────────────────

/** The three workflow stages (D7). `complete` = results exist and nothing is left to review
 *  (D18 — derived, never clicked). */
export type WorkflowStage = 'setup' | 'review' | 'complete';

/** The hub's cached copy of `deriveStage` output (v15). A CACHE: the live view always
 *  re-derives; the hub reads this through `formatStageBadge` with the staleness overlay. */
export interface StageCache {
  readonly stage: WorkflowStage;
  /** Review units still `To review`. */
  readonly pending: number;
  /** Review units the design makes `Not applicable` (e.g. per-animal rows of a destructive
   *  design — SF-3). Disclosed, never counted as pending. */
  readonly notApplicable: number;
}

/** The recipe the analyst last confirmed (v15, AD-5). The SNAPSHOT is persisted, not only its
 *  fingerprint, so `Revert changes` can re-apply it after a reopen. */
export interface ConfirmedRecipe {
  readonly recipe: AnalysisRecipe;
  readonly fingerprint: string;
  /** ISO timestamp of the Confirm. */
  readonly at: string;
}

/** One manual λz point override (v15, AD-2). Identity is FRAME-INDEPENDENT:
 *  `profileKey` = the results-row `SUBJECT` string; `timeKey` = the observation time in the
 *  canonical `getString` form (CSV round-trip safe); `ordinal` = its index among the same
 *  profile's rows at that same time (replicate samples). `concKey` is the concentration in the
 *  same canonical form, checked on reopen: a changed value reverts the override to auto. */
export interface LzOverride {
  readonly profileKey: string;
  readonly timeKey: string;
  readonly ordinal: number;
  readonly concKey: string;
  readonly action: 'force-in' | 'force-out';
  /** The analyst's display name as recorded (what the UI shows). */
  readonly user: string;
  /** The analyst's platform user id (IR-12 / Part 11 attribution: a name can change or be shared,
   *  an id cannot). Absent on a record made without a platform user (headless, tests). */
  readonly userId?: string;
  /** ISO timestamp. */
  readonly at: string;
}

/** The three review actions that need a reason (D11/D12). A point edit is not one of them —
 *  the edit is its own justification. */
export type ReviewAction = 'accept' | 'acknowledge' | 'exclude';

/** One chosen reason. `code` is the stable registry id (`review-reasons.ts`) or `OTHER`;
 *  `text` is the text shown at decision time (free text for `OTHER`). */
export interface ReviewReason {
  readonly code: string;
  readonly text: string;
}

/** Why a decision stopped being current. `reset` = the analyst pressed Reset; `auto-reset` =
 *  the fit changed since the decision (fingerprint mismatch on compute / reopen — D13);
 *  `replaced` = a λz point edit after the decision, or a new decision on the same unit. */
export type SupersededBy = 'reset' | 'auto-reset' | 'replaced';

/** One review decision (v15). `fingerprint` is the decision fingerprint of what the analyst
 *  saw; a later mismatch resets the decision (D13). */
export interface ReviewDecisionRecord {
  readonly action: ReviewAction;
  readonly reasons: ReadonlyArray<ReviewReason>;
  /** The analyst's display name as recorded (what the UI shows). */
  readonly user: string;
  /** The analyst's platform user id, as on {@link LzOverride.userId}. */
  readonly userId?: string;
  /** ISO timestamp. */
  readonly at: string;
  readonly fingerprint: string;
  /** Set on a history record only. */
  readonly supersededBy?: SupersededBy;
  /** ISO timestamp, set on a history record only. */
  readonly supersededAt?: string;
}

/** Every decision ever recorded for one review unit (v15). `current === null` = the unit is
 *  back to `To review`; `history` keeps each superseded record, oldest first — a superseded
 *  decision is never obscured (21 CFR 11.10(e), D13). */
export interface ReviewLedger {
  readonly unitKey: string;
  readonly current: ReviewDecisionRecord | null;
  readonly history: ReadonlyArray<ReviewDecisionRecord>;
}

/** Where an entry came from (v16 union): moved from the standalone store, or written by a
 *  pipeline run. */
export type AnalysisOrigin = StandaloneOrigin | PipelineOrigin;

/** Where a study entry came from when it was moved out of the standalone Analyses store. */
export interface StandaloneOrigin {
  readonly kind: 'standalone';
  /** The standalone container's label at move time (the `analyses/<label>` link segment). */
  readonly label: string;
  readonly containerId: string;
}

/** The pipeline run that wrote an entry (v16, review Q1 + 6d) — run provenance, and what the
 *  re-run's "untouched" test compares against (`pipelineUntouched`, `pipeline-save.ts`). */
export interface PipelineOrigin {
  readonly kind: 'pipeline';
  readonly runId: string;
  /** ISO timestamp of the run. */
  readonly runAt: string;
  /** The runner's login. */
  readonly user: string;
  /** The sci-comp version the run computed with. */
  readonly engineVersion: string;
  /** The mapping the run used (the parsed `ncaCompute` mapping, verbatim). */
  readonly mapping: Readonly<Record<string, unknown>>;
  /** The rules JSON the run used, verbatim. */
  readonly rules: string;
  /** The partition identity this entry belongs to — its label before `uniqueLabel`,
   *  route-qualified for a sparse analysis. A re-run updates the entry whose partition matches. */
  readonly partition: string;
  /** `computeRecipeHash(recipeOf(entry))` of the recipe the run wrote. */
  readonly recipeHash: string;
  /** `RECIPE_HASH_VERSION` the hash was computed with; another version cannot be compared. */
  readonly recipeHashVersion: number;
  /** The id of the pipeline entry this one superseded when the sampling design changed between
   *  runs (dense ↔ sparse). The old entry is left in place; both rows say so in the study grid. */
  readonly supersedes?: string;
}

/** What the caller may do with a study — the server's `~can_*` (nca-ems-storage-migration; the
 *  in-JSON permission lists are gone, the platform enforces row security). */
export interface StudyAccessFlags {
  readonly canEdit: boolean;
  readonly canDelete: boolean;
  readonly canShare: boolean;
}

/** A study (or a standalone container — `kind: 'standalone'`, one analysis) with its analyses. */
export interface NcaStudy {
  readonly id: string;
  readonly label: string;
  readonly kind: 'study' | 'standalone';
  readonly createdAt: string;
  readonly updatedAt: string;
  readonly access: StudyAccessFlags;
  readonly analyses: ReadonlyArray<AnalysisEntry>;
  /** Set on a pipeline-created study: the owner and delivered study key (its EMS business key). */
  readonly pipelineKey?: {readonly owner: string; readonly studyKey: string};
}

// ─── FNV-1a content hash (AC-B6 — synchronous, deterministic, order-stable) ──

/** FNV-1a 32-bit over a UTF-16-code-unit string. Synchronous and deterministic
 *  (R2-F4: NOT Web Crypto, which is async and would cascade). The exact function
 *  matters less than that it is fixed, deterministic, synchronous, order-stable. */
export function fnv1a(s: string): string {
  let h = 0x811c9dc5;
  for (let i = 0; i < s.length; i++) {
    h ^= s.charCodeAt(i);
    // FNV prime 16777619, kept in 32-bit via Math.imul.
    h = Math.imul(h, 0x01000193);
  }
  // Unsigned 32-bit hex.
  return (h >>> 0).toString(16).padStart(8, '0');
}

// The source signature of a stored snapshot is `snapshotSignature` (nca-repository.ts): a full
// content hash written with the snapshot (SD-19). The sampled `computeContentHash`, the mtime
// signature and `hasDrifted` were removed with the AppData store (SD-4): a stored snapshot is
// immutable, so it cannot drift — only the engine can (`isStale`).
