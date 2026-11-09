/**
 * Column auto-detection for the docked NCA mapping pane.
 *
 * Pure-function detector that maps `DG.DataFrame` columns to the 8 mapping-pane
 * fields (subject / time / concentration / dose / bodyWeight / blqFlag / lloq /
 * route) by exact (case-insensitive) name match against per-field synonym
 * tables, gated by a type filter. Used by `mountMappingPane` to pre-populate
 * inputs at view construction (silently, and only at initial mount).
 *
 * **Boundary invariants (do not regress):**
 *
 * 1. **No mutation.** Detector reads `df` but never writes column tags, never
 *    sets `column.semType`, never adds columns. semType emission is a separate
 *    forward-deferred research stream.
 * 2. **No side effects.** Two calls on the same DF return reference-equal
 *    results (same `DG.Column` objects). No I/O. No state.
 * 3. **Determinism.** Synonym tables are `const` and frozen at module-load.
 * 4. **Bucket structure is the disambiguation tie-breaker.** The `actualTime`
 *    and `nominalTime` buckets are NOT positionally-resolved within one flat
 *    synonym list; they are structurally separate. Bare `time` IS in the
 *    `actualTime` bucket — it does not fall through. PRD FR-106 mandates actual
 *    time for AUC calculation; nominal time for grouping. Auto-filling the
 *    Time *calculation* input from a `nominal_time` column would silently shift
 *    AUC in sparse-sampling studies (UC-04 / UC-06).
 * 5. **BLQ INT-0/1 acceptance.** NONMEM `MDV`, SAS exports, and most CRO
 *    deliveries encode BLQ as INT `0`/`1`. The BOOL-only filter would reject
 *    the dominant encoding. `isBoolLikeNumeric` widens the filter to BOOL or
 *    INT-with-cardinality-≤-2-and-values-⊆-{0,1}; symmetric with the
 *    mapping-pane BLQ picker filter.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

/** Heuristic resolve of a role column by name regex. Returns the source column name
 *  (verbatim) or null. Moved from `views/nca-mapping-pane.ts` (nca-pipeline-save-rework F1):
 *  the sparse compute that moved to `sparse-nca.ts` resolves its BLQ / sampling-type roles
 *  with it. */
export function findColByPattern(df: DG.DataFrame, pattern: RegExp): string | null {
  for (const c of df.columns.names()) if (pattern.test(c)) return c;
  return null;
}

/**
 * Result of `detectMapping`. Each field holds the `DG.Column` reference if
 * exactly one synonym-and-type-matching candidate was found in the input
 * DataFrame; otherwise `undefined`.
 *
 * The `time` field is filled from the `actualTime` bucket only (bare 'time' IS
 * in this bucket); `nominalTime` never promotes to the calculation slot.
 */
export interface DetectedMapping {
  readonly subject?: DG.Column;
  readonly time?: DG.Column;
  readonly concentration?: DG.Column;
  readonly dose?: DG.Column;
  readonly bodyWeight?: DG.Column;
  readonly blqFlag?: DG.Column;
  readonly lloq?: DG.Column;
  /** Per-row route column (string or int). Filled when exactly one synonym-and-
   *  type-matching candidate exists; the mapping pane encodes it to RouteCode
   *  via the tolerant normalizer (per-row-route AD-1). */
  readonly route?: DG.Column;

  /**
   * Unit string parsed from a unit-suffixed *matched* column name
   * (e.g. `Conc_ng_mL` → `'ng/mL'`, `Dose_mg_kg` → `'mg/kg'`, `Time_hr` →
   * `'h'`, `BodyWeight_g` → `'g'`). Present ONLY when the matched column
   * carries an unambiguous unit token that maps to a value the mapping-pane
   * dropdown actually offers; absent otherwise (the pane keeps its canonical
   * default). Closes GAP-11: without this, mapping `Conc_ng_mL` and leaving
   * the Concentration dropdown at canonical mg/L silently mis-scales
   * CL/Vz/AUC by 1000× (λz/t½ are scale-invariant, so the error is invisible
   * in the "shape" parameters). Each unit field is independent of its column
   * field being filled here — the token is parsed from the matched column's
   * own name, so it appears only alongside its column match.
   */
  readonly concentrationUnit?: string;
  readonly doseUnit?: string;
  readonly timeUnit?: string;
  readonly bodyWeightUnit?: string;
}

/** Per-field synonym tables (lowercased — match is `name.toLowerCase()`).
 *  `Object.freeze` enforces synthesis §1.4 rule 7 (determinism) at runtime —
 *  the `ReadonlyArray` type prevents mutation through typed references but
 *  not through `any`-cast access; freezing makes the invariant load-bearing.
 */
const SUBJECT_SYNONYMS: ReadonlyArray<string> = Object.freeze([
  'subject', 'subj', 'subjectid', 'subjid', 'usubjid', 'animal', 'participant', 'pid',
  // Preclinical species-id column names (GAP-14): the canonical mouse/rat fixtures
  // use `MouseID` / `AnimalNo`, which the clinical-leaning list above missed —
  // leaving Subject unmapped on first Compute (the sparse UC-04 collapse). A bare
  // `id` token is intentionally excluded: `uniqueMatch` is exact-equality, so a
  // generic column literally named `ID` (a row/record/sample id, not the subject)
  // would silently fill Subject. Suffixed names like `StudyID`/`GroupID` are
  // already safe under exact-match — listing only specific `*id` forms keeps the
  // detector from over-matching.
  'mouseid', 'ratid', 'dogid', 'monkeyid', 'animalid', 'animalno',
  'animal_id', 'animal_no',
]);
// `actualTime` bucket members (PRD FR-106 canonical for AUC calculation). The
// `nominalTime` bucket (`nominal_time`, `nomtime`, `tpt`) is intentionally NOT
// represented here — nominal columns flow through the picker dropdown for
// manual selection, never through this detector, so the calculation slot
// cannot be silently filled from a nominal-time column.
const ACTUAL_TIME_SYNONYMS: ReadonlyArray<string> = Object.freeze([
  'actual_time', 'acttime', 'sampling_time', 'time', 'time_h', 'time_min',
  // Unit-suffixed preclinical / CRO variants (e.g. rat dataset `Time_hr`). The
  // bucket stays ACTUAL-only — no nominal-time names — so nominal time never
  // silently fills the AUC-calculation slot.
  'time_hr', 'time_hrs', 'time_hour', 'time_hours', 'time_minutes',
  'actual_time_hr', 'actualtime_hr', 'sampling_time_hr',
]);
const CONCENTRATION_SYNONYMS: ReadonlyArray<string> = Object.freeze([
  'conc', 'concentration', 'dv', 'aval', 'cobs', 'dv_obs',
  // Unit-suffixed preclinical / CRO variants, mass-based units (e.g. rat
  // dataset `Conc_ng_mL`). Molar units (nmol/L, …) are intentionally absent —
  // the engine rejects molar downstream, so auto-mapping them would only lead
  // to a deferred error; the user picks those manually.
  'conc_ng_ml', 'conc_ug_ml', 'conc_µg_ml', 'conc_mg_ml', 'conc_mg_l', 'conc_ug_l',
  'concentration_ng_ml', 'concentration_ug_ml', 'concentration_mg_l',
]);
const DOSE_SYNONYMS: ReadonlyArray<string> = Object.freeze([
  'dose', 'amt', 'dose_mg', 'dose_amount', 'total_dose',
  // Unit-suffixed preclinical / CRO variants, incl. per-kg (e.g. rat dataset
  // `Dose_mg_kg`). The unit-conversion layer canonicalizes per-kg → mg via
  // body weight; the validator hard-blocks per-kg dose without a BW source.
  'dose_mg_kg', 'dose_mgkg', 'dose_ug_kg', 'dose_µg_kg', 'dose_mg_l',
]);
const BODY_WEIGHT_SYNONYMS: ReadonlyArray<string> = Object.freeze([
  'weight', 'wt', 'bw', 'bodyweight', 'body_weight', 'bw_kg', 'weight_kg',
  // Unit-suffixed preclinical / CRO variants (e.g. rat dataset `BodyWeight_g`).
  'bodyweight_g', 'body_weight_g', 'bw_g', 'weight_g',
  'bodyweight_kg', 'body_weight_kg',
]);
const BLQ_SYNONYMS: ReadonlyArray<string> = Object.freeze(
  ['blq', 'bloq', 'belowlloq', 'lloq_flag']);
// Route column synonyms (per-row-route AD-1). `exroute` is the CDISC SEND
// EX-domain route variable; `route`/`routecode` are the common CRO headers.
const ROUTE_SYNONYMS: ReadonlyArray<string> = Object.freeze(
  ['route', 'routecode', 'route_code', 'exroute', 'dosing_route', 'admin_route',
    'administration_route', 'route_of_administration']);
const LLOQ_SYNONYMS: ReadonlyArray<string> = Object.freeze(
  ['lloq', 'lod', 'qlim', 'dl', 'detection_limit']);

/**
 * Per-field unit-suffix dictionaries: `[lowercased name suffix, unit string]`.
 * The unit strings are EXACTLY the mapping-pane dropdown values (TIME / CONC /
 * DOSE unit choices + body-weight `kg`/`g`) so the pane can assign them
 * directly. `parseUnitToken` returns the FIRST entry whose suffix the
 * lowercased column name ends with, so entries are ordered longest-first
 * within each overlap class (e.g. `_mg_kg` before `_mg`; `_kg` before `_g`)
 * to prevent a shorter suffix from shadowing a longer one.
 *
 * Conservative by construction (GAP-11 "only set on an unambiguous token"):
 * tokens with no supported dropdown value are deliberately ABSENT, so the
 * pane keeps its canonical default rather than guessing. Notably omitted:
 *  - `_mg_ml` (`mg/mL` concentration): synonym `conc_mg_ml` matches the
 *    column, but mg/mL is not a Concentration dropdown choice → leave
 *    canonical. (Distinct from `_mg_l` → `mg/L`, which IS a choice and IS
 *    present below — the suffix-vs-unit distinction is load-bearing here.)
 *  - molar concentration units (engine rejects molar downstream — never
 *    auto-selected here);
 *  - `mg/L` dose (synonym `dose_mg_l` exists but mg/L is not a Dose choice).
 *
 * Frozen at module-load (determinism, mirrors the synonym tables).
 */
const CONC_UNIT_SUFFIXES: ReadonlyArray<readonly [string, string]> = Object.freeze([
  ['_ng_ml', 'ng/mL'],
  ['_µg_ml', 'µg/mL'], ['_ug_ml', 'µg/mL'],
  ['_µg_l', 'µg/L'], ['_ug_l', 'µg/L'],
  ['_mg_l', 'mg/L'],
] as ReadonlyArray<readonly [string, string]>);
const DOSE_UNIT_SUFFIXES: ReadonlyArray<readonly [string, string]> = Object.freeze([
  ['_mg_kg', 'mg/kg'], ['_mgkg', 'mg/kg'],
  ['_µg_kg', 'µg/kg'], ['_ug_kg', 'µg/kg'],
  ['_mg', 'mg'],
] as ReadonlyArray<readonly [string, string]>);
const TIME_UNIT_SUFFIXES: ReadonlyArray<readonly [string, string]> = Object.freeze([
  ['_minutes', 'min'], ['_min', 'min'],
  ['_hours', 'h'], ['_hour', 'h'], ['_hrs', 'h'], ['_hr', 'h'], ['_h', 'h'],
] as ReadonlyArray<readonly [string, string]>);
const BODY_WEIGHT_UNIT_SUFFIXES: ReadonlyArray<readonly [string, string]> = Object.freeze([
  ['_kg', 'kg'], ['_g', 'g'],
] as ReadonlyArray<readonly [string, string]>);

/**
 * Parse a unit string from a matched column name by suffix lookup. Returns the
 * first dictionary entry whose lowercased suffix the (lowercased) name ends
 * with, or `undefined` when no supported unit token is present. Pure.
 */
function parseUnitToken(
  name: string,
  dict: ReadonlyArray<readonly [string, string]>,
): string | undefined {
  const lower = name.toLowerCase();
  for (const [suffix, unit] of dict)
    if (lower.endsWith(suffix)) return unit;
  return undefined;
}

const STRING_OR_INT = (c: DG.Column): boolean =>
  c.type === DGNS.COLUMN_TYPE.STRING || c.type === DGNS.COLUMN_TYPE.INT;
const NUMERIC = (c: DG.Column): boolean =>
  c.type === DGNS.COLUMN_TYPE.FLOAT || c.type === DGNS.COLUMN_TYPE.INT;

/**
 * Predicate used by BOTH the autodetect BLQ field (AD-4) AND the mapping-pane
 * BLQ picker dropdown filter (AD-7) — consistent affordance: the detector
 * never accepts a class the manual picker rejects, and vice versa.
 *
 * Accepts BOOL unconditionally, INT only when every non-null value is `0` or
 * `1`. Loop-only (no dependency on `DG.Column.stats.uniqueCount`, which may or
 * may not include nulls per Datagrok platform conventions — FM-7 from R2
 * synthesis review). Single pass over the INT raw buffer; O(N) where N is the
 * column row count.
 *
 * Empty / all-null INT column returns `false` (no observable 0/1 pattern).
 */
export function isBoolLikeNumeric(col: DG.Column): boolean {
  if (col.type === DGNS.COLUMN_TYPE.BOOL) return true;
  if (col.type !== DGNS.COLUMN_TYPE.INT) return false;
  const raw = col.getRawData() as Int32Array;
  let saw0 = false;
  let saw1 = false;
  for (let i = 0; i < raw.length; i++) {
    const v = raw[i];
    if (v === DGNS.INT_NULL) continue;
    if (v === 0) {saw0 = true; continue;}
    if (v === 1) {saw1 = true; continue;}
    return false;
  }
  return saw0 || saw1;
}

/**
 * Map columns to mapping-pane fields by name + type filter. Pure: no
 * mutation of `df`. Idempotent: two calls return reference-equal results.
 *
 * For each field, scans columns left-to-right and collects every column whose
 * lowercased name matches any synonym AND whose type passes the field filter.
 * If exactly 1 candidate → use it; if 0 or ≥ 2 → leave the field empty.
 *
 * The `time` field (NCA calculation input, PRD FR-106 canonical) is resolved
 * from the `actualTime` bucket only:
 *
 * - 1 actualTime candidate (any nominalTime state) → `time` fills with the
 *   actualTime candidate. `nominalTime` columns remain available in the
 *   picker dropdown for manual selection (this detector does not touch them).
 * - 0 actualTime candidates → `time` stays empty. Nominal-time columns do NOT
 *   fall through to the calculation slot (it would silently shift AUC).
 * - ≥ 2 actualTime candidates → `time` stays empty (cardinality rule applies
 *   inside the bucket; ambiguous actualTime does NOT fall through to
 *   nominalTime).
 */
/**
 * Result of {@link detectAnalyteColumns}. `multiAnalyte` is `true` when ANY
 * of the regex / synonym / cardinality signals fire; `evidence` is a list
 * of human-readable strings naming the detected signal (consumed by the
 * `MULTI_ANALYTE_NOT_SUPPORTED` error message in `validators.ts`).
 */
export interface AnalyteDetection {
  readonly multiAnalyte: boolean;
  readonly evidence: ReadonlyArray<string>;
  /** The columns behind `evidence` — lets a column-level check skip a column already reported. */
  readonly columns: ReadonlyArray<string>;
}

/** Regex patterns identifying analyte-named numeric columns. Case-insensitive
 *  match against `column.name`. `Parent_*` / `Metabolite_*` are the dominant
 *  CRO LC-MS/MS multi-analyte column-naming conventions; `M[0-9]_Conc` covers
 *  metabolite-by-number variants. Frozen to keep determinism (the synonym
 *  vocabulary is closed for v1; widening is a deliberate future change). */
const ANALYTE_REGEXES: ReadonlyArray<RegExp> = Object.freeze([
  /^parent_.+$/i,
  /^metabolite_.+$/i,
  /^m\d+_conc$/i,
]);

/** Exact name match (case-insensitive) for analyte categorical columns. */
const ANALYTE_CATEGORICAL_NAMES: ReadonlyArray<string> = Object.freeze(['analyte', 'compound']);

/**
 * Detect multi-analyte structure in `df`. Two orthogonal signals:
 *
 *  (a) **Wide form** — ≥ 2 columns whose names match one of {@link
 *      ANALYTE_REGEXES} (`Parent_*`, `Metabolite_*`, `M\d+_Conc`), OR a
 *      `parent_*` AND `metabolite_*` pair (count ≥ 2 across the wide signal
 *      space).
 *  (b) **Long form** — a STRING column named `analyte` / `compound`
 *      (case-insensitive) whose distinct non-null values are ≥ 2 (per-row
 *      analyte categorical).
 *
 *  Either signal sets `multiAnalyte = true` with evidence appended.
 *
 *  Pure: no mutation; idempotent on the same `df`.
 *
 *  AD-9 / E6 / scientist memo §3. Existing scan-once pattern (rule-5 reuse
 *  of the column-autodetect file).
 */
export function detectAnalyteColumns(df: DG.DataFrame): AnalyteDetection {
  const cols = df.columns.toList();
  const evidence: string[] = [];
  const flagged: string[] = [];

  const wideMatches: string[] = [];
  for (const col of cols) {
    for (const re of ANALYTE_REGEXES) {
      if (re.test(col.name)) {
        wideMatches.push(col.name);
        break;
      }
    }
  }
  if (wideMatches.length >= 2) {
    evidence.push(`multiple analyte-named columns: ${wideMatches.join(', ')}`);
    flagged.push(...wideMatches);
  }

  for (const col of cols) {
    if (!ANALYTE_CATEGORICAL_NAMES.includes(col.name.toLowerCase())) continue;
    if (col.type !== DGNS.COLUMN_TYPE.STRING) continue;
    // Distinct-value count via Set over the row values; cheap at NCA sizes
    // (< 10k rows) and avoids platform-specific `stats.uniqueCount` semantics
    // (FM-7 in nca-mapping-autofill R2).
    const distinct = new Set<string>();
    const n = df.rowCount;
    for (let i = 0; i < n; i++) {
      const v = col.get(i);
      if (v !== null && v !== undefined && v !== '')
        distinct.add(String(v));
      if (distinct.size >= 2) break;
    }
    if (distinct.size >= 2) {
      evidence.push(`per-row analyte categorical column '${col.name}' with ≥ 2 distinct values`);
      flagged.push(col.name);
    }
  }

  return {multiAnalyte: evidence.length > 0, evidence, columns: flagged};
}

/**
 * Industry names of a per-row compound / analyte IDENTIFIER column (BUG-07 interim, DATA-GAP-14):
 * `CompoundID`, `Compound_Name`, `CmpdID`, `Analyte Name`, `AnalyteID`, `Drug`, `DrugName`,
 * `Test Article`, `Article`, `TestItem`, `Molecule`, `MoleculeID`, `Substance`, `Chemical`,
 * `Parent` / `Metabolite` (categorical), `DGK_ID` (the Datagrok demo compound id), … — matched on
 * the name with spaces / underscores / hyphens removed, case-insensitive. The name alone never
 * blocks: a group must also hold ≥ 2 of its values ({@link detectMixedCompoundGroups}).
 */
const COMPOUND_ID_STEMS = 'compound|cmpd|analyte|drug|testarticle|testitem|article|molecule|substance|chemical|' +
  'parent|metabolite|dgk';
const COMPOUND_ID_NAME = new RegExp(`^(${COMPOUND_ID_STEMS})(id|name|code|no)?$`);

/** True when `name` reads as a compound / analyte identifier column (see {@link COMPOUND_ID_NAME}). */
export function isCompoundIdColumnName(name: string): boolean {
  return COMPOUND_ID_NAME.test(name.toLowerCase().replace(/[\s_-]+/g, ''));
}

/** One pooling group (a subject, a dose group, a composite cell) that holds ≥ 2 values of a
 *  recognised compound-identifier column. */
export interface MixedCompoundGroup {
  readonly column: string;
  /** The group's key as the pooling code forms it (a subject id, `10`, `D10_F_Day28`, …). */
  readonly group: string;
  /** The distinct values seen in that group (≥ 2), in first-seen order. */
  readonly values: ReadonlyArray<string>;
}

/**
 * Pooling-across-compounds guard (BUG-07 pattern, 2026-10-05): every place that POOLS rows by a
 * key — a profile per subject, a sparse composite per dose group, a NONMEM composite per
 * dose×sex×day cell — silently mixes compounds when that key has no compound dimension. For each
 * column named like a compound identifier ({@link isCompoundIdColumnName}), report the first
 * group (`keyOf(row)`) that holds ≥ 2 of its values. Rows with a missing compound value, or a
 * `null` key, are ignored. `skipColumns` = columns that ARE the pooling key (never reported
 * against themselves).
 *
 * Name AND structure on purpose — the structural test alone was tried and reverted (2026-07-14):
 * a per-observation flag (`BLQ`) or a mislabelled covariate (`Sex`) varies within a subject just
 * like a cassette compound id, so a name-blind scan false-positives on clean data. A compound
 * column that is CONSTANT within each group (a batch screen with disjoint subjects —
 * `05_multicompound_batch_screen.csv`) is fine and is not reported.
 */
export function detectMixedCompoundGroups(
  df: DG.DataFrame, keyOf: (row: number) => string | null,
  skipColumns: ReadonlyArray<string> = [],
): MixedCompoundGroup[] {
  const out: MixedCompoundGroup[] = [];
  for (const col of df.columns.toList()) {
    if (skipColumns.includes(col.name) || !isCompoundIdColumnName(col.name)) continue;
    const byGroup = new Map<string, string[]>();
    for (let i = 0; i < df.rowCount; i++) {
      if (col.isNone(i)) continue;
      const v = col.getString(i);
      if (v === '') continue;
      const g = keyOf(i);
      if (g === null) continue;
      const seen = byGroup.get(g);
      if (seen === undefined) byGroup.set(g, [v]);
      else if (!seen.includes(v)) {
        seen.push(v);
        out.push({column: col.name, group: g, values: seen});
        break;
      }
    }
  }
  return out;
}

/**
 * The analyst-facing explanation of a pooling group that mixes compounds — ONE wording for every
 * pooling surface. `groupLabel` names the group that pools (`Dose group '10'`,
 * `Composite cell 'D10_F_Day28'`).
 */
export function mixedCompoundMessage(m: MixedCompoundGroup, groupLabel: string): string {
  return `${groupLabel} pools ${m.values.length} values of '${m.column}' ` +
    `(${m.values.join(', ')}) — different compounds would be averaged into one composite and ` +
    `every parameter would be wrong. Not supported yet: run one '${m.column}' value per ` +
    'analysis (filter the table).';
}

/**
 * Result of {@link detectEventDrivenColumns}. All column refs are optional — a
 * NONMEM file legitimately omits some roles (`03_rat_multiroute_nonmem.csv` has
 * no `ADDL`/`II`/`SS`). `isEventDriven` is the trigger: the mapping pane
 * *proposes* the Event-driven (NONMEM) mode when it is true (never silent —
 * AD-2). This all-optional shape is deliberately distinct from the
 * required-present `EventContract` the transform consumes (architect SIMPLIFY-2).
 */
export interface EventDrivenDetection {
  readonly isEventDriven: boolean;
  readonly evid?: DG.Column;
  readonly amt?: DG.Column;
  readonly cmt?: DG.Column;
  readonly rate?: DG.Column;
  readonly mdv?: DG.Column;
  readonly addl?: DG.Column;
  readonly ii?: DG.Column;
  readonly ss?: DG.Column;
}

/** NONMEM data-item synonyms (exact case-insensitive match, frozen — mirrors the
 *  detector convention). Single-token each: the NONMEM column names are
 *  standardized data items, not free CRO headers. */
const EVENT_SYNONYMS: Readonly<Record<keyof Omit<EventDrivenDetection, 'isEventDriven'>, ReadonlyArray<string>>> =
  Object.freeze({
    evid: ['evid'], amt: ['amt'], cmt: ['cmt'], rate: ['rate'],
    mdv: ['mdv'], addl: ['addl'], ii: ['ii'], ss: ['ss'],
  });

/** Valid NONMEM EVID values (`0` obs, `1` dose, `2` other, `3` reset, `4`
 *  reset+dose). The trigger requires ≥ 1 row carrying one of these (R2 N-5). */
const VALID_EVID_VALUES: ReadonlySet<number> = new Set([0, 1, 2, 3, 4]);

/**
 * Detect NONMEM event-driven layout in `df` (AD-2, P-1). Returns the resolved
 * column refs for the NONMEM data items plus `isEventDriven`.
 *
 * **Trigger (three-part conjunction, R2 N-5):** an `EVID` column **and** an
 * `AMT` column are present **and** ≥ 1 row carries a numeric `EVID ∈ {0..4}`.
 * The conjunction avoids a false trigger on a coincidental "Event ID" / Phoenix
 * `AMT`-as-nominal column (a single-column trigger is too weak).
 *
 * Pure: no mutation; idempotent on the same `df`. The id / time / concentration
 * roles are detected by {@link detectMapping} (DV ∈ `CONCENTRATION_SYNONYMS`,
 * TIME ∈ `ACTUAL_TIME_SYNONYMS`, ID ∈ `SUBJECT_SYNONYMS`); the pane bridges both
 * into an `EventContract` at the confirm step.
 */
export function detectEventDrivenColumns(df: DG.DataFrame): EventDrivenDetection {
  const cols = df.columns.toList();
  const find = (role: keyof typeof EVENT_SYNONYMS): DG.Column | undefined =>
    uniqueMatch(cols, EVENT_SYNONYMS[role], NUMERIC);

  const evid = find('evid');
  const amt = find('amt');

  let hasValidEvid = false;
  if (evid !== undefined) {
    const raw = evid.getRawData();
    const n = df.rowCount;
    for (let i = 0; i < n; i++) {
      const v = raw[i];
      if (v === DGNS.INT_NULL || v === DGNS.FLOAT_NULL) continue;
      if (VALID_EVID_VALUES.has(Math.round(v))) {hasValidEvid = true; break;}
    }
  }
  const isEventDriven = evid !== undefined && amt !== undefined && hasValidEvid;

  const out: {-readonly [K in keyof EventDrivenDetection]: EventDrivenDetection[K]} = {isEventDriven};
  if (evid !== undefined) out.evid = evid;
  if (amt !== undefined) out.amt = amt;
  for (const role of ['cmt', 'rate', 'mdv', 'addl', 'ii', 'ss'] as const) {
    const c = find(role);
    if (c !== undefined) out[role] = c;
  }
  return out;
}

export function detectMapping(df: DG.DataFrame): DetectedMapping {
  const cols = df.columns.toList();

  const subject = uniqueMatch(cols, SUBJECT_SYNONYMS, STRING_OR_INT);
  const time = uniqueMatch(cols, ACTUAL_TIME_SYNONYMS, NUMERIC);
  const concentration = uniqueMatch(cols, CONCENTRATION_SYNONYMS, NUMERIC);
  const dose = uniqueMatch(cols, DOSE_SYNONYMS, NUMERIC);
  const bodyWeight = uniqueMatch(cols, BODY_WEIGHT_SYNONYMS, NUMERIC);
  const blqFlag = uniqueMatch(cols, BLQ_SYNONYMS, isBoolLikeNumeric);
  const lloq = uniqueMatch(cols, LLOQ_SYNONYMS, NUMERIC);
  const route = uniqueMatch(cols, ROUTE_SYNONYMS, STRING_OR_INT);

  const detected: {-readonly [K in keyof DetectedMapping]: DetectedMapping[K]} = {};
  if (subject !== undefined) detected.subject = subject;
  if (time !== undefined) detected.time = time;
  if (concentration !== undefined) detected.concentration = concentration;
  if (dose !== undefined) detected.dose = dose;
  if (bodyWeight !== undefined) detected.bodyWeight = bodyWeight;
  if (blqFlag !== undefined) detected.blqFlag = blqFlag;
  if (lloq !== undefined) detected.lloq = lloq;
  if (route !== undefined) detected.route = route;

  // Parse the unit token from each matched column's own name (GAP-11). Each
  // is independent of the others; emit only when an unambiguous, dropdown-
  // supported token is present so the pane keeps its canonical default
  // otherwise. Unit fields are set, never overwritten by a column being empty.
  if (concentration !== undefined) {
    const u = parseUnitToken(concentration.name, CONC_UNIT_SUFFIXES);
    if (u !== undefined) detected.concentrationUnit = u;
  }
  if (dose !== undefined) {
    const u = parseUnitToken(dose.name, DOSE_UNIT_SUFFIXES);
    if (u !== undefined) detected.doseUnit = u;
  }
  if (time !== undefined) {
    const u = parseUnitToken(time.name, TIME_UNIT_SUFFIXES);
    if (u !== undefined) detected.timeUnit = u;
  }
  if (bodyWeight !== undefined) {
    const u = parseUnitToken(bodyWeight.name, BODY_WEIGHT_UNIT_SUFFIXES);
    if (u !== undefined) detected.bodyWeightUnit = u;
  }
  return detected;
}

function uniqueMatch(
  cols: ReadonlyArray<DG.Column>,
  synonyms: ReadonlyArray<string>,
  typeFilter: (c: DG.Column) => boolean,
): DG.Column | undefined {
  let match: DG.Column | undefined;
  for (const col of cols) {
    if (!synonyms.includes(col.name.toLowerCase())) continue;
    if (!typeFilter(col)) continue;
    if (match !== undefined) return undefined;
    match = col;
  }
  return match;
}

