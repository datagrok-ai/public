/**
 * Study unified grid builder (study-unified-grid AD-3) — the pure
 * `NcaStudy → DG.DataFrame` projection the L2 study overview renders as a native
 * `DG.Grid`. Concatenates every analysis's cached `compoundTier` rows into one
 * frame, tagged with the source **analysis label** (provenance) and a
 * back-pointer `(analysisId, compoundId, route)` for the row-click drill-in.
 *
 * **Computes nothing (thesis #1).** Every value is read straight off the cached
 * {@link CompoundTierRow} (which `buildSummaryCache` derived via the validated
 * `aggregateNcaResults`). This module only *projects* — union of parameter
 * columns, drop-all-empty, basis-routing (CL vs CL/F by route), the
 * internal→display label map, and warnings/NE carry-through. No mean/CV/n is
 * recomputed here.
 *
 * **Basis-routing (AC-3b).** Each case is route-homogeneous, so its `route`
 * 1:1-determines the CL/Vz basis: an intravascular case's `cl`/`vz` land in the
 * `CL`/`Vz` columns, an extravascular case's in `CL/F`/`Vz/F`. The unused column
 * is null (genuinely not-applicable — never a "missing" blank).
 *
 * **Vocabulary firewall (CT-9).** Internal ids (`MEAN_*`/`CV_PCT_*`/`N_*`/
 * `CL_BASIS`/`PCT_EXTRAP`/`DOSE_NORM`) never reach an output column name — the
 * display-label map + the `nca-param`/`nca-stat` grouping tags are the only
 * things the grid header sees.
 */

import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import {f64Column} from './adapter-utils';
import {nca} from '@datagrok-libraries/sci-comp';

import {ROUTE_STRING_TO_CODE} from './profile-inputs';
// One-way edge: workflow-stage no longer imports from this module (review: cycle removed).
import {
  formatNotReviewed, STALE_ENGINE_FLAG, formatPipelineFlags, pipelineFlagCount, PER_ANIMAL_SPARSE_NOTE,
  formatSupersedes, formatSupersededBy,
} from './workflow-stage';
import type {
  NcaStudy, AnalysisEntry, CompoundTierRow, CompoundCurve, MeanCv, GroupKeyValue,
} from './study-schema';
import {isBasisStale, isCurvePending} from './nca-repository';
import type {LoadedStudy} from './nca-repository';
import {CL_PER_KG_UNIT, VOLUME_PER_KG_UNIT} from './unit-conversion';

/** Provenance / back-pointer column display names. */
export const ANALYSIS_COL = 'Analysis';
/** Inline per-compound curve thumbnail column (v5, chemist-sar-grid S4). A
 *  `DG.COLUMN_TYPE.DATA_FRAME` column — each cell is a small frame of the
 *  compound's semi-log geomean curve (Time / Conc / Series), rendered in-cell by
 *  the overview via the native `plot.line` (no custom renderer). `null` cell where
 *  a case has no cached curve (placeholder row / pre-v5 recompute). */
export const CURVE_COL = 'C–t';
/** Display token for a per-key group cell that does not apply to a row: the row is
 *  whole-analysis (empty `groupValues`) AND its subpopulation carries no inventory value
 *  for that key, or the key is genuinely absent from the row's analysis. `—` (not blank)
 *  keeps an inapplicable cell visibly distinct from a genuinely-empty value. **The generic
 *  `Group` column was RETIRED**: the grid ALWAYS emits one named column per grouping key
 *  (per-dimension filter / sort, the reference-tool convention). */
export const NO_GROUP = '—';
export const ROUTE_COL = 'Route';
/** Warnings column display name (renamed from "Flags", user-directed 2026-07-07) —
 *  the case's human-readable warning text (v6: NE + warning messages, NEWLINE-joined
 *  per message; '' when clean). The text is the cell VALUE so it survives table
 *  export (replaced the render-only glyph + the separate hidden column). The const
 *  stays `FLAGS_COL` internally (identifier only; the display header is "Warnings"). */
export const FLAGS_COL = 'Warnings';
/** Hidden back-pointer to the source analysis (AD-4 hides it; used by row-click
 *  to re-open + hydrate the case's profile). Internal `__`-prefixed name (never shown)
 *  so a user grouping key named `analysisId` cannot CASE-collide with it at
 *  `DataFrame.fromColumns` (DG column names are case-insensitive). */
export const ANALYSIS_ID_COL = '__analysisId';
/** Hidden back-pointer to the case's composite group identity (v8,
 *  multi-column-group-by) — the tier row's `compoundId` (a `compositeGroupKey`
 *  U+0001 join for a composite selection). The row-click drill-in reads THIS (not a
 *  display column), so it works uniformly regardless of the visible per-key columns.
 *  Internal `__`-prefixed name (never shown) so a user grouping key named `CompoundID`
 *  (the common compound-screen case) gets a CLEAN `CompoundID` column instead of
 *  CASE-colliding with this back-pointer at `DataFrame.fromColumns`
 *  (study-grid-grouping-columns regression fix, 2026-07-13). */
export const COMPOUND_ID_COL = '__compoundId';

/** Column tag: the parameter a Mean/%CV/n column belongs to (drives the
 *  friendly `${stat} ${param}` header + fragility detection), e.g. `'Cmax'`,
 *  `'CL/F'`. */
export const PARAM_TAG = 'nca-param';
/** Column tag: the statistic role (`'Mean'` | `'%CV'` | `'n'`). */
export const STAT_TAG = 'nca-stat';
/** Column tag marking a passthrough (annotation) column (v10, study-grid-passthrough-columns).
 *  Value is `'structure'` (a MOLECULE-designated column — `configureStudyGrid` sets
 *  `semType = MOLECULE` on it) or `'text'`. LOAD-BEARING (probe PP-3): `studyGridVisibleOrder`
 *  classes ANY non-known / non-param column as a grouping key, so WITHOUT this discriminator a
 *  passthrough column would leak into the grouping-key zone. The tag places it in its OWN zone
 *  (structure passthrough first, after the grouping keys, before Route). Mirrors {@link PARAM_TAG}. */
export const PASSTHROUGH_TAG = 'nca-passthrough';
/** {@link PASSTHROUGH_TAG} value for a structure-designated passthrough column. */
export const PASSTHROUGH_STRUCT = 'structure';
/** {@link PASSTHROUGH_TAG} value for a plain-text passthrough column. */
export const PASSTHROUGH_TEXT = 'text';

/** Placeholder flag text for a v1-migrated analysis not yet recomputed in the
 *  v2 shape (AC-2b / R1-F7): visible, never silently absent, never a wrong
 *  collapsed row. */
export const RECOMPUTE_FLAG = 'recompute to include per-compound rows';

/** Re-exported from its home in `workflow-stage.ts` (the ONE wording module) — kept importable
 *  from here for the existing study-grid readers. */
export {STALE_ENGINE_FLAG};

/** True when a route string denotes an intravascular case (→ true CL/Vz, not
 *  CL/F·Vz/F). Route 1:1-determines the CL basis (results-writer). */
function isIntravascular(route: string): boolean {
  const code = ROUTE_STRING_TO_CODE[route];
  return code === nca.ROUTE_IV_BOLUS || code === nca.ROUTE_IV_INFUSION;
}

/** Canonical display units for the headline params (sci-comp canonicalizes to these before
 *  compute — so they are stable constants). CL / Vz / Vss are per kg on every cached row since
 *  nca-cl-bw-normalised-units (AD-8), so their units are fixed too. */
const FIXED_UNITS: Readonly<Record<string, string>> = {
  'Cmax': 'mg/L',
  'AUC∞': 'mg·h/L',
  't½': 'h',
  'MRT': 'h',
  'Tmax': 'h',
  '%AUCextrap': '%',
  'CL': CL_PER_KG_UNIT,
  'CL/F': CL_PER_KG_UNIT,
  'Vz': VOLUME_PER_KG_UNIT,
  'Vz/F': VOLUME_PER_KG_UNIT,
  'Vss': VOLUME_PER_KG_UNIT,
};

/** The study-grid parameter groups whose cached value depends on the reporting basis. On a row of
 *  a basis-stale analysis ({@link isBasisStale}) they are blank: an old absolute mean must not sit
 *  under a mL/min/kg header (architect B3). Every other cached number is basis-independent. */
export const BASIS_DEPENDENT_PARAMS: ReadonlySet<string> = new Set(['CL', 'CL/F', 'Vz', 'Vz/F', 'Vss']);

/** `df.temp` key on the study-grid frame: per grid row, 1 when the row's analysis is basis-stale
 *  and its {@link BASIS_DEPENDENT_PARAMS} cells were blanked. Read via {@link studyCellBasisStale}. */
export const STUDY_BASIS_STALE_TEMP = 'ncaStudyBasisStale';

/** Was this study-grid cell blanked because its analysis was cached under an older reporting basis
 *  (until someone with edit access reopens it)? False for every other cell. */
export function studyCellBasisStale(df: DG.DataFrame, colName: string, row: number): boolean {
  const label = df.col(colName)?.getTag(PARAM_TAG) ?? '';
  if (!BASIS_DEPENDENT_PARAMS.has(label) || row < 0) return false;
  const flags = df.temp[STUDY_BASIS_STALE_TEMP] as Uint8Array | null | undefined;
  return (flags?.[row] ?? 0) === 1;
}

/** `df.temp` key on the study-grid frame: per grid row, 1 when the row's curve is still being
 *  loaded (ems-modification step 5 — the overview reads the result rows without curves and loads
 *  them after the first render). Its C–t cell is then empty for now, not "no curve". Read via
 *  {@link studyCurvePending}. */
export const STUDY_CURVE_PENDING_TEMP = 'ncaStudyCurvePending';

/** Is the C–t thumbnail of grid row `row` still being loaded? */
export function studyCurvePending(df: DG.DataFrame, row: number): boolean {
  const flags = df.temp[STUDY_CURVE_PENDING_TEMP] as Uint8Array | null | undefined;
  return row >= 0 && (flags?.[row] ?? 0) === 1;
}

/** One projected grid row (flattened case or a placeholder). */
interface GridRow {
  readonly analysisLabel: string;
  readonly analysisId: string;
  readonly compoundId: string;
  readonly route: string;
  readonly iv: boolean;
  /** The row's analysis cached its CL / Vz / Vss under an older reporting basis (blanked). */
  readonly basisStale: boolean;
  /** The row's curve is still being loaded (see {@link STUDY_CURVE_PENDING_TEMP}). */
  readonly curvePending: boolean;
  /** null on a placeholder row (no per-compound data yet). */
  readonly row: CompoundTierRow | null;
  /** Flags cell = the case's human-readable flag text (v6): `NE` (unfit) and/or the
   *  warning messages, `; `-joined; `RECOMPUTE_FLAG` on a placeholder; '' when clean.
   *  The text lives IN the cell (not a render-only glyph) so it survives table export
   *  (user-directed 2026-07-07). */
  readonly flags: string;
}

/** `df.temp` key on the study-grid frame: per parameter label, the per-row count of profiles whose
 *  λz-dependent value was withheld by an Exclude λz decision (`MeanCv.nExcluded`, F-S5.2). Only
 *  labels with at least one excluded profile are present. Read via {@link studyCellExclusion}. */
export const STUDY_EXCLUDED_N_TEMP = 'ncaStudyExcludedN';

/**
 * The Exclude λz disclosure of one study-grid cell (F-S5.2): for a λz-dependent parameter group's
 * Mean / %CV / n cell whose row had profiles excluded by review — `n` (the profiles that
 * contribute), `total` (`n` + excluded — the `n = a of b` of the tooltip) and `notReportable`
 * (EVERY contributor excluded: the empty Mean / %CV reads `not reportable`, never the not-estimated
 * blank). null when nothing on that row was excluded.
 */
export function studyCellExclusion(
  df: DG.DataFrame, colName: string, row: number,
): {n: number; total: number; notReportable: boolean} | null {
  const col = df.col(colName);
  const label = col?.getTag(PARAM_TAG) ?? '';
  const stat = col?.getTag(STAT_TAG) ?? '';
  if (col === null || label === '' || stat === '' || row < 0) return null;
  const map = df.temp[STUDY_EXCLUDED_N_TEMP] as Record<string, Int32Array> | null | undefined;
  const ex = map?.[label]?.[row] ?? 0;
  if (!(ex > 0)) return null;
  const nCol = df.col(`${label} n`);
  const n = nCol !== null && !nCol.isNone(row) ? Number(nCol.get(row)) : 0;
  return {n, total: n + ex, notReportable: n === 0 && stat !== 'n'};
}

/** The Warnings cell text for a case (v6): the compact `NE` token (when the headline
 *  exposure AUC∞ is unfit) followed by each warning's human-readable message. Joined
 *  by NEWLINE — NOT `; ` — because the messages THEMSELVES contain `; ` (e.g.
 *  "…re-map the dose unit; otherwise verify…"), so a `; ` delimiter would tear one
 *  message across bullets. Newline is a true per-message boundary (messages carry no
 *  newlines). '' when clean. The text is the cell VALUE (exportable), not a
 *  render-only glyph — so "save table" carries the warning info. */
function caseFlagsText(row: CompoundTierRow): string {
  const parts: string[] = [];
  // v15 (F-S5.2 / D11): Exclude λz decisions withheld λz-dependent values of some profiles — say
  // how many, and name the state when EVERY profile was withheld ("not reportable", distinct from
  // the not-estimated `NE`). The per-parameter `n` columns already exclude them.
  const nEx = row.thalf.nExcluded ?? row.aucInf.nExcluded ?? 0;
  if (nEx > 0) {
    parts.push(nEx >= row.n ?
      'λz-dependent parameters not reportable (every profile excluded by review)' :
      `λz excluded for ${nEx} of ${row.n} profiles — t½, AUC∞, CL, Vz, MRT not reportable for them`);
  }
  if (Number.isNaN(row.aucInf.mean) && nEx < row.n) parts.push('NE');
  for (const w of row.warnings) parts.push(w.message);
  return parts.join('\n');
}

/**
 * Typographic canonical id for a grouping-key label (study-grid-grouping-columns E2):
 * **casefold + strip WHITESPACE only** (punctuation PRESERVED). Merges cross-analysis
 * key columns that differ only in case/whitespace into ONE column:
 * `route`/`ROUTE`/`'Route '` → `'route'`; `'Dose (mg/kg)'`/`'Dose(mg/kg)'` → `'dose(mg/kg)'`.
 * Punctuation is KEPT so signed/qualified labels stay DISTINCT — `'Day -1'` → `'day-1'` vs
 * `'Day 1'` → `'day1'` (a CDISC baseline-vs-dosing-day distinction; stripping `-` would
 * wrongly merge two real groups, R1 finding 1) and `'Dose'` → `'dose'` vs `'Dose (mg/kg)'`
 * → `'dose(mg/kg)'` (genuinely different). Semantic synonyms (`Sex`↔`Gender`) also stay
 * distinct — synonym reconciliation is OUT of scope (typographic merge only). Exported +
 * pure so the union-key and the join-match read from ONE canonical-id source (no regex
 * drift → no silent fragmentation) and the typographic boundary is unit-testable in isolation.
 */
export function normalizeGroupKey(label: string): string {
  return label.toLowerCase().replace(/\s+/g, '');
}

/**
 * A row's {@link GroupKeyValue} list with canonical-duplicate keys removed (first-seen wins) —
 * the SAME merge {@link buildStudyGrid}'s key columns apply via {@link normalizeGroupKey}. Two
 * entries collapse when their keys share a canonical id (the `friendlyName`-collision case where
 * two result columns carry the same friendly label). Applied at `groupValues` CONSTRUCTION so the
 * persisted projection matches how every consumer already reads it (the grid merges to one key
 * column; the identity label / curve `series.group` must print the value ONCE). Display-only — the
 * numeric identity (`compoundId`, computed independently via `compositeGroupKey`) is untouched.
 */
export function dedupGroupValuesByCanonical(
  groupValues: readonly GroupKeyValue[],
): GroupKeyValue[] {
  const seen = new Set<string>();
  const out: GroupKeyValue[] = [];
  for (const kv of groupValues) {
    const canonical = normalizeGroupKey(kv.key);
    if (seen.has(canonical)) continue;
    seen.add(canonical);
    out.push(kv);
  }
  return out;
}

/**
 * The per-key cell value for a row and a merged key column, by canonical id (E4
 * 3-step precedence):
 *  1. the row is GROUPED by this key (a `groupValues` entry whose canonical id matches) →
 *     that entry's value;
 *  2. else the key is present-but-ungrouped in this row's subpopulation (a `keyInventory`
 *     entry whose canonical id matches) → its distinct values joined WITH counts
 *     (`F (15); M (15)`) — the honest-composition disclosure (thesis #4). A whole-analysis /
 *     pooled row lands here (empty `groupValues` → step 1 never matches), which is exactly
 *     where covariate-imbalance disclosure matters most;
 *  3. else `—` — the key is genuinely absent from this row's subpopulation.
 * NO `isWholeAnalysis` short-circuit (R1 finding 2b) — a pooled row discloses its counted
 * composition, not an unconditional `—`. Per-ROW inventory (R2 finding 3) — a
 * route/dose-norm-split row shows ITS OWN subpopulation's counts, never an analysis-wide one.
 */
function keyCellValue(row: CompoundTierRow | null, canonical: string): string {
  if (row === null) return NO_GROUP; // placeholder row (no per-compound data yet)
  const grouped = row.groupValues.find((kv) => normalizeGroupKey(kv.key) === canonical);
  if (grouped !== undefined) return grouped.value;
  const inv = row.keyInventory.find((e) => normalizeGroupKey(e.key) === canonical);
  if (inv !== undefined && inv.values.length > 0)
    return inv.values.map((vc) => `${vc.value} (${vc.count})`).join('; ');
  return NO_GROUP;
}

/**
 * The cell value for a PASSTHROUGH (annotation) column, by canonical id (v10,
 * study-grid-passthrough-columns). DISTINCT from {@link keyCellValue} — a passthrough cell shows
 * the BARE value when unique (so a lone SMILES renders as a 2D molecule; a name renders as plain
 * text), NOT a counted join. Multiplicity rules (thesis #4 honesty, user-approved 2026-07-14):
 *  - exactly 1 distinct value → the bare value;
 *  - >1 distinct, TEXT column → the counted-join disclosure `A (3); B (2)` (reuses the
 *    keyInventory join form — visible, honest);
 *  - >1 distinct, STRUCTURE column → `— (N)` (the conservative default: a MOLECULE-semType cell
 *    cannot render "N structures" text or one depiction standing for a pool of different
 *    compounds — thesis #4 — so it discloses the distinct COUNT, which survives table export and
 *    surfaces on hover, rather than a misleading single structure);
 *  - absent from this row's subpopulation → `—`.
 * Per-ROW (the row's own `passthrough` inventory) and per-entry `isStructure` (the persisted
 * designation), so a coarse pooled row honestly discloses its composition.
 */
function passthroughCellValue(row: CompoundTierRow | null, canonical: string): string {
  if (row === null) return NO_GROUP; // placeholder row (no per-compound data yet)
  const e = row.passthrough.find((p) => normalizeGroupKey(p.key) === canonical);
  if (e === undefined || e.values.length === 0) return NO_GROUP;
  if (e.values.length === 1) return e.values[0].value; // bare (lone SMILES → molecule; name → text)
  // many structures → bare '—' cell (a MOLECULE-semType cell can't host a text count); the distinct
  // COUNT is disclosed on the cell TOOLTIP by configureStudyGrid (user-approved 2026-07-14, RG-B).
  if (e.isStructure) return NO_GROUP;
  return e.values.map((vc) => `${vc.value} (${vc.count})`).join('; '); // many text → counted join
}

/**
 * The flattened (analysis, tier-row) sequence the study grid renders, in grid-row order: one entry
 * per `compoundTier` row, plus ONE placeholder entry (`row === null`) per analysis whose
 * `compoundTier` is empty. The SINGLE source of grid-row order — {@link buildStudyGrid} maps it to
 * grid rows, and the overview's `configureStudyGrid` tooltip maps a grid row index back to its tier
 * row (so a many-structure passthrough cell, shown as a bare `—`, can disclose its distinct count on
 * hover). Keeping both consumers on one order function prevents the row-index drift a re-derived
 * flattening would risk.
 */
export function studyGridRowTiers(
  study: NcaStudy,
): ReadonlyArray<{a: AnalysisEntry; r: CompoundTierRow | null; i: number}> {
  const out: {a: AnalysisEntry; r: CompoundTierRow | null; i: number}[] = [];
  for (const a of study.analyses) {
    // A null cache (v15 — never confirmed) narrows to the recompute placeholder row, exactly
    // like an empty tier: a study row never dereferences a missing cache (AC-S1.7.2).
    const tier = a.summaryCache?.compoundTier ?? [];
    if (tier.length === 0) out.push({a, r: null, i: -1});
    else tier.forEach((r, i) => out.push({a, r, i}));
  }
  return out;
}

/**
 * Project a study's cached summaries into the unified grid DataFrame. Row count
 * = Σ compoundTier lengths, plus one placeholder row per analysis whose
 * `compoundTier` is empty (v1-migrated / not-yet-computed). Empty parameter
 * groups (no case carries a value) are dropped. Returns an empty-but-shaped
 * frame (the provenance columns only) when the study has no analyses.
 */
export function buildStudyGrid(
  study: LoadedStudy,
  isStale: (a: AnalysisEntry) => boolean = () => false,
): DG.DataFrame {
  // A row of an engine-stale analysis carries STALE_ENGINE_FLAG as the FIRST line of its Flags
  // cell (v14 F-N.6). The predicate is the coordinator's `isStale` — injected, not re-derived,
  // so every reader of cached rows shares one definition; the default (never stale) keeps the
  // frame-shape callers/tests unchanged.
  // v15 (nca-workflow-stages F-S6.4 / D14 / D19): an analysis added with profiles still to review
  // carries `N not reviewed` as the NEXT line — same mechanism and position rule as the stale mark
  // (stale first, then not-reviewed, then the case's own flags). Wording from ONE helper.
  // v16 (nca-pipeline-save-rework, approved design 2026-10-15): then `N pipeline flags` (the
  // pipeline's own reasons — AC-9.2), the sparse per-animal disclosure (AC-5.8b), and the
  // supersede pair on both rows of a design flip (AC-4.2b). Every wording from ONE helper.
  const labelById = new Map(study.analyses.map((x) => [x.id, x.label]));
  const supersededBy = new Map<string, string>();
  for (const x of study.analyses) {
    const by = x.origin?.kind === 'pipeline' ? x.origin.supersedes : undefined;
    if (by !== undefined && labelById.has(by)) supersededBy.set(by, x.label);
  }
  const withStale = (a: AnalysisEntry, flags: string): string => {
    const supersedes = a.origin?.kind === 'pipeline' && a.origin.supersedes !== undefined ?
      labelById.get(a.origin.supersedes) : undefined;
    const by = supersededBy.get(a.id);
    const lines = [
      isStale(a) ? STALE_ENGINE_FLAG : '', formatNotReviewed(a.stageCache),
      // The overview reads the stored count, not the flags (ems-modification step 5).
      formatPipelineFlags(a.pipelineFlagCount ?? pipelineFlagCount(a.pipelineFlags)),
      // Guarded: hand-built entries (fixtures, the pre-v16 grid tests) may carry no `modes`.
      a.modes?.sparse?.on === true ? PER_ANIMAL_SPARSE_NOTE : '',
      supersedes !== undefined ? formatSupersedes(supersedes) : '',
      by !== undefined ? formatSupersededBy(by) : '',
      flags,
    ].filter((l) => l !== '');
    return lines.join('\n');
  };
  // Grid rows in the shared studyGridRowTiers order (single source of row order — the overview's
  // passthrough tooltip maps a grid row index back to its tier row via the same function).
  const rows: GridRow[] = studyGridRowTiers(study).map(({a, r, i}) => r === null ?
    // Placeholder: visible, flagged to recompute, never a wrong collapsed row.
    {
      analysisLabel: a.label, analysisId: a.id,
      compoundId: a.label, route: a.summaryCache?.analysisTier.route ?? '',
      iv: isIntravascular(a.summaryCache?.analysisTier.route ?? ''),
      basisStale: false, curvePending: false, row: null, flags: withStale(a, RECOMPUTE_FLAG),
    } : {
      analysisLabel: a.label, analysisId: a.id,
      compoundId: r.compoundId, route: r.route, iv: isIntravascular(r.route),
      basisStale: isBasisStale({resultsBasis: a.resultsBasis, hasCache: a.summaryCache !== null}),
      curvePending: isCurvePending(study, a.id, i), row: r,
      // v6: the human-readable flag text (NE + sci-comp's messages) IN the cell.
      flags: withStale(a, caseFlagsText(r)),
    });

  const n = rows.length;
  // Grouping-key columns:
  // ALWAYS one named column per grouping key — the generic `Group` column and the
  // `perKeyMode ≥ 2` gate are RETIRED. Key columns are the UNION of the grouping keys any
  // row carries, MERGED across analyses on a TYPOGRAPHIC canonical id (`normalizeGroupKey`)
  // so `route`/`ROUTE`/`Route ` fold into ONE column; the header is the first-seen friendly
  // label. A study where NO analysis groups → the union is empty → zero key columns. Flat
  // columns only — an embedded `DG.Viewer.grid` ignores `.columnGroups` spanning headers
  // (guardrails footgun); adjacency + the per-column header carry the grouping, and each key
  // stays independently filterable/sortable in the native grid.
  interface KeyCol {readonly canonical: string; readonly label: string; colName: string;}
  const keyCols: KeyCol[] = [];
  const canonicalIndex = new Map<string, KeyCol>();
  for (const r of rows) {
    for (const kv of r.row?.groupValues ?? []) {
      const canonical = normalizeGroupKey(kv.key);
      if (!canonicalIndex.has(canonical)) {
        const kc: KeyCol = {canonical, label: kv.key, colName: kv.key};
        canonicalIndex.set(canonical, kc);
        keyCols.push(kc);
      }
    }
  }

  // Passthrough (annotation) columns (v10, study-grid-passthrough-columns): the UNION of the
  // passthrough keys any row carries, MERGED across analyses on the SAME typographic canonical id
  // as the grouping keys (E2 — `SMILES`/`smiles` fold into one column; `SMILES` vs
  // `Canonical_SMILES` stay separate, the accepted boundary), first-seen friendly label + isStructure.
  // Ordered STRUCTURE-FIRST then source order (user-approved 2026-07-14) — a stable sort keeps
  // first-seen order within each group. Tagged PASSTHROUGH_TAG so the overview's ordering /
  // semType logic places + designates them (PP-3). Empty when no analysis selected passthrough columns.
  interface PassCol {
    readonly canonical: string; readonly label: string; readonly isStructure: boolean; colName: string;
  }
  const passCols: PassCol[] = [];
  const passIndex = new Map<string, PassCol>();
  for (const r of rows) {
    for (const e of r.row?.passthrough ?? []) {
      const canonical = normalizeGroupKey(e.key);
      if (!passIndex.has(canonical)) {
        const pc: PassCol = {canonical, label: e.key, isStructure: e.isStructure, colName: e.key};
        passIndex.set(canonical, pc);
        passCols.push(pc);
      }
    }
  }
  passCols.sort((a, b) => (a.isStructure ? 0 : 1) - (b.isStructure ? 0 : 1));

  const cols: DG.Column[] = [
    DGNS.Column.fromStrings(ANALYSIS_COL, rows.map((r) => r.analysisLabel)),
  ];
  // Collision-free column names: a user key label could clash with a provenance /
  // parameter column (`DataFrame.fromColumns` throws on a duplicate). Reserve the fixed
  // names first, then finalize each key column's display name (first-seen friendly label).
  const used = new Set<string>([
    ANALYSIS_COL, ROUTE_COL, CURVE_COL, FLAGS_COL, ANALYSIS_ID_COL, COMPOUND_ID_COL]);
  for (const kc of keyCols) {kc.colName = uniqueColName(kc.label, used); used.add(kc.colName);}
  for (const kc of keyCols)
    cols.push(DGNS.Column.fromStrings(kc.colName, rows.map((r) => keyCellValue(r.row, kc.canonical))));
  // Passthrough columns after the grouping keys (structure-first). Collision-safe names (case-
  // insensitive, guards a SMILES/name label clashing with a key or provenance column). Each
  // tagged PASSTHROUGH_TAG=structure|text so the overview zones + semType-designates it (PP-3);
  // the tag is the load-bearing discriminator that keeps it OUT of the grouping-key zone.
  for (const pc of passCols) {pc.colName = uniqueColName(pc.label, used); used.add(pc.colName);}
  for (const pc of passCols) {
    const col = DGNS.Column.fromStrings(
      pc.colName, rows.map((r) => passthroughCellValue(r.row, pc.canonical)));
    col.setTag(PASSTHROUGH_TAG, pc.isStructure ? PASSTHROUGH_STRUCT : PASSTHROUGH_TEXT);
    cols.push(col);
  }
  cols.push(DGNS.Column.fromStrings(ROUTE_COL, rows.map((r) => r.route)));

  // ── Parameter groups (Mean | %CV | n), dropping all-empty. ────────────────
  const excludedN: Record<string, Int32Array> = {};
  // valueOf returns the MeanCv to place in THIS group for a row, or null (blank
  // = not-applicable, e.g. the wrong CL basis, or a placeholder row). A basis-dependent group
  // is blank on a basis-stale row (B3) — the row keeps its stale Warnings line.
  const addGroup = (label: string, valueOf: (r: GridRow) => MeanCv | null): void => {
    const means = new Float64Array(n);
    const cvs = new Float64Array(n);
    const ns = new Int32Array(n);
    const exN = new Int32Array(n); // F-S5.2: profiles withheld by Exclude λz, per row
    let anyExcluded = false;
    let anyMean = false;
    let anyBlanked = false;
    const basisDependent = BASIS_DEPENDENT_PARAMS.has(label);
    for (let i = 0; i < n; i++) {
      const own = valueOf(rows[i]);
      const blanked = basisDependent && rows[i].basisStale;
      if (blanked && own !== null && Number.isFinite(own.mean)) anyBlanked = true;
      const mv = blanked ? null : own;
      const mean = (mv !== null && Number.isFinite(mv.mean)) ? mv.mean : NaN;
      means[i] = Number.isNaN(mean) ? DGNS.FLOAT_NULL : mean;
      cvs[i] = (mv !== null && Number.isFinite(mv.cvPct)) ? mv.cvPct : DGNS.FLOAT_NULL;
      ns[i] = (mv !== null && mv.n > 0) ? mv.n : DGNS.INT_NULL;
      exN[i] = mv?.nExcluded ?? 0;
      if (exN[i] > 0) anyExcluded = true;
      if (!Number.isNaN(mean)) anyMean = true;
    }
    // Drop the whole group if no case carries a value (hide-all-empty, AC-3a). A group whose
    // only values were blanked as basis-stale stays: its empty cells carry the stale tooltip,
    // and a vanished CL column would not say why it is gone (B3, all-stale study).
    if (!anyMean && !anyBlanked) return;

    const meanCol = f64Column(`${label} Mean`, means);
    const cvCol = f64Column(`${label} %CV`, cvs);
    const nCol = DGNS.Column.fromInt32Array(`${label} n`, ns);
    // One unit per parameter since every cached row reports one basis (AD-8).
    const unit = FIXED_UNITS[label];
    if (unit !== undefined) meanCol.setTag('units', unit);
    cvCol.setTag('units', '%');
    meanCol.setTag(PARAM_TAG, label); meanCol.setTag(STAT_TAG, 'Mean');
    cvCol.setTag(PARAM_TAG, label); cvCol.setTag(STAT_TAG, '%CV');
    nCol.setTag(PARAM_TAG, label); nCol.setTag(STAT_TAG, 'n');
    cols.push(meanCol, cvCol, nCol);
    if (anyExcluded) excludedN[label] = exN;
  };

  // A standalone single-value parameter column (design B3: `%AUCextrap` is "its
  // own single column (no %CV/n triple)"). No STAT_TAG → excluded from the
  // Mean|%CV|n column groups; the header shows the param label directly.
  const addSingle = (
    label: string, unit: string, valueOf: (r: GridRow) => number | null,
  ): void => {
    const vals = new Float64Array(n);
    let anyVal = false;
    for (let i = 0; i < n; i++) {
      const v = valueOf(rows[i]);
      const finite = v !== null && Number.isFinite(v);
      vals[i] = finite ? (v as number) : DGNS.FLOAT_NULL;
      if (finite) anyVal = true;
    }
    if (!anyVal) return;
    const col = f64Column(label, vals);
    col.setTag('units', unit);
    col.setTag(PARAM_TAG, label);
    cols.push(col);
  };

  const active = (r: GridRow): boolean => r.row !== null;
  // Visible-param order per design B3: Cmax, AUC∞, %AUCextrap, t½, then the
  // basis-routed CL/CL/F, Vz/Vz/F pairs. %AUCextrap is a single column.
  addGroup('Cmax', (r) => active(r) ? r.row!.cmax : null);
  // Tmax (context-panel-ct-profile B, user-directed): carried on the grid frame so the row-click
  // "NCA parameters" Forms viewer (bound to this frame) can display it. It is HIDDEN in the grid —
  // `studyGridVisibleOrder` never lists it, so the SAR scan is unaffected — but present in the
  // frame for the form. `TMAX` is already a per-row aggregate on the tier (`tmax`).
  addGroup('Tmax', (r) => active(r) ? r.row!.tmax : null);
  addGroup('AUC∞', (r) => active(r) ? r.row!.aucInf : null);
  addSingle('%AUCextrap', '%', (r) => active(r) ? r.row!.pctExtrap.mean : null);
  addGroup('t½', (r) => active(r) ? r.row!.thalf : null);
  // Basis-routed CL / Vz: intravascular → CL/Vz, extravascular → CL/F·Vz/F.
  addGroup('CL', (r) => (active(r) && r.iv) ? r.row!.cl : null);
  addGroup('CL/F', (r) => (active(r) && !r.iv) ? r.row!.cl : null);
  addGroup('Vz', (r) => (active(r) && r.iv) ? r.row!.vz : null);
  addGroup('Vz/F', (r) => (active(r) && !r.iv) ? r.row!.vz : null);
  // v3 (chemist-sar-grid): MRT (route-uniform) + Vss (IV-only, routed like Vz —
  // sci-comp gives NaN for EV, so the group self-drops on an all-EV screen). No
  // `/F` counterpart for Vss (an apparent steady-state volume is not defined).
  // Visible-order/chemist-priority (lead trio AUC∞/MRT/%F) is applied at the grid
  // layer (S3), not here — this projection only makes the columns available.
  addGroup('MRT', (r) => active(r) ? r.row!.mrt : null);
  // v4 (chemist-sar-grid S2): %F — the third leg of the SAR lead trio. Cross-
  // profile mean/%CV/n (persisted on the tier; the study grid sees only the
  // cache). `fPct` is null on IV rows AND where no paired IV+EV arm exists; the
  // addGroup all-empty guard omits the WHOLE column when no case carries a finite
  // %F (Decision 2: all-IV / unpaired screen → no %F column — e.g. fixture 05).
  // Unit '%' (dimensionless percentage); no basis-routing (%F is route-uniform
  // over its EV strata).
  addGroup('%F', (r) => active(r) ? (r.row!.fPct ?? null) : null);
  addGroup('Vss', (r) => (active(r) && r.iv) ? r.row!.vss : null);

  // Visible Flags text (NE + warning messages, in-cell so it exports) + hidden
  // back-pointer.
  cols.push(DGNS.Column.fromStrings(FLAGS_COL, rows.map((r) => r.flags)));
  cols.push(DGNS.Column.fromStrings(ANALYSIS_ID_COL, rows.map((r) => r.analysisId)));
  // Hidden composite-identity back-pointer (v8): the row-click drill-in reads this,
  // not a display column, so it is display-form-independent (single `Group` or per-key).
  cols.push(DGNS.Column.fromStrings(COMPOUND_ID_COL, rows.map((r) => r.compoundId)));

  const df = DGNS.DataFrame.fromColumns(cols);
  df.name = study.label;
  df.temp[STUDY_EXCLUDED_N_TEMP] = excludedN; // F-S5.2 — read by studyCellExclusion
  df.temp[STUDY_BASIS_STALE_TEMP] = Uint8Array.from(rows, (r) => r.basisStale ? 1 : 0); // AD-8 B3
  df.temp[STUDY_CURVE_PENDING_TEMP] = Uint8Array.from(rows, (r) => r.curvePending ? 1 : 0);

  // v5 (chemist-sar-grid S4): the inline curve thumbnail as a DATA_FRAME column —
  // each cell holds the compound's small semi-log geomean frame, rendered in-cell
  // by the overview's native `plot.line` (frontend-ui-gate Rule 0: copies the
  // platform in-cell chart pattern, `compute-utils/fitting-view.ts`). Built here so
  // the projection stays the single place that shapes the grid frame. `null` cell
  // where a case has no cached curve (placeholder / pre-v5) or its curve is still loading
  // (STUDY_CURVE_PENDING_TEMP). Only added when at least one case carries — or will carry — a
  // curve (hide-all-empty parity with the param groups).
  if (rows.some((r) => r.curvePending || (r.row?.curve?.series.length ?? 0) > 0)) {
    df.columns.addNew(CURVE_COL, DGNS.COLUMN_TYPE.DATA_FRAME).init((i: number) => {
      const curve = rows[i].row?.curve ?? null;
      return (curve !== null && curve.series.length > 0) ? buildCurveFrame(curve) : null;
    });
  }
  return df;
}

/** A column name unique against `used` — appends ` (2)`, ` (3)`, … on collision.
 *  Guards the per-key grouping columns (a user key's friendly label could collide with a
 *  provenance / parameter column name; `DataFrame.fromColumns` throws on a duplicate).
 *  **Case-INSENSITIVE** (study-grid-grouping-columns regression fix): DG column names are
 *  case-insensitive, so `'CompoundID'` and `'compoundId'` collide at `fromColumns` — a
 *  case-sensitive guard would let a differently-cased key crash the grid. */
function uniqueColName(name: string, used: ReadonlySet<string>): string {
  const usedLc = new Set<string>();
  for (const u of used) usedLc.add(u.toLowerCase());
  if (!usedLc.has(name.toLowerCase())) return name;
  for (let k = 2; ; k++) {
    const candidate = `${name} (${k})`;
    if (!usedLc.has(candidate.toLowerCase())) return candidate;
  }
}

/** Placeholder `Series` category for a thumbnail line whose group label is empty
 *  (a whole-analysis pooled curve — its `series.group` is `''` since v9/E3). The
 *  overview's `plot.line` colors lines by the `Series` CATEGORY, and an EMPTY category
 *  renders BLACK; substituting a non-empty placeholder restores the standard first-palette
 *  line color. The label is never shown (the sparkline's legend/hover are off), so any
 *  stable non-empty token works (study-grid-grouping-columns render fix, 2026-07-13). */
const THUMB_SERIES_DEFAULT = 'curve';

/** Column accumulator for {@link buildCurveFrame} / {@link buildOverlayCurveFrame} — the ONE
 *  frame shape both the grid thumbnail and the context panel consume, so correspondence between
 *  them holds by construction (F1, context-panel-ct-profile). */
interface CurveFrameCols {
  readonly time: number[]; readonly conc: number[];
  readonly geoLower: number[]; readonly geoUpper: number[]; readonly series: string[];
}

/** Append one curve's points to the accumulator. `labelOverride` (F3 overlay) forces the
 *  `Series` category to a row-identity label; when null the series' own `group` is used (mapped
 *  to {@link THUMB_SERIES_DEFAULT} when empty — a whole-analysis pooled curve's `''` group would
 *  otherwise render black). The optional v11 band (`geoLower`/`geoUpper`) is emitted when the
 *  series carries it, else `FLOAT_NULL` at every point (a legacy band-absent curve → no whiskers,
 *  the panel discloses "recompute to show"; guarded with `!== undefined` — a plain-TS optional,
 *  NOT the `view.temp`/`getTag` null-vs-undefined dart-boundary footgun). */
function pushCurve(acc: CurveFrameCols, curve: CompoundCurve, labelOverride: string | null): void {
  for (const s of curve.series) {
    const label = labelOverride ?? (s.group === '' ? THUMB_SERIES_DEFAULT : s.group);
    const hasBand = s.geoLower !== undefined && s.geoUpper !== undefined;
    for (let i = 0; i < s.time.length; i++) {
      const t = s.time[i];
      const c = s.geoMean[i];
      // A tuple with no time or no concentration is NOT a plottable point — drop
      // it whole (time / conc / series / band together, so the arrays stay
      // parallel) rather than emitting a placeholder. This makes the LEGACY
      // persisted path behave exactly like the fresh one, which already skips
      // such points at `nca-analysis-open.ts:330-332`.
      //
      // Why the type annotation is not enough: the study cache is JSON
      // (`study-coordinator.ts`, no replacer) and JSON cannot represent NaN, so
      // a pre-guard persisted curve round-trips an absent value as `null` while
      // its declared type stays `number`. `Float64Array.from([null])` would have
      // yielded a real `0.0` — a datum at zero on a semi-log axis, not a gap —
      // and `f64Column` rejects it outright (AC-4b). Skipping is the honest
      // answer to both.
      if (typeof t !== 'number' || !Number.isFinite(t)) continue;
      if (typeof c !== 'number' || !Number.isFinite(c)) continue;
      acc.time.push(t); acc.conc.push(c); acc.series.push(label);
      acc.geoLower.push(bandEdge(hasBand ? s.geoLower![i] : undefined));
      acc.geoUpper.push(bandEdge(hasBand ? s.geoUpper![i] : undefined));
    }
  }
}

/** One persisted band edge, normalised onto the platform null. An absent band is
 *  already representable (`FLOAT_NULL` → no whisker, "recompute to show"), so a
 *  legacy `null` / `undefined` / `NaN` edge maps onto that existing meaning
 *  instead of tripping `f64Column`'s non-numeric guard. A real edge (including a
 *  pre-existing `FLOAT_NULL`) passes through untouched. */
function bandEdge(v: number | null | undefined): number {
  return (typeof v === 'number' && !Number.isNaN(v)) ? v : DGNS.FLOAT_NULL;
}

/** Materialize the accumulator into a DataFrame. Columns `Time` / `Conc` (= geoMean) /
 *  `geoLower` / `geoUpper` / `Series`. The thumbnail (`THUMB_OPTS`: y=`Conc`, split by `Series`)
 *  ignores the two band columns; the context panel binds whiskers to them. */
function framifyCurve(acc: CurveFrameCols): DG.DataFrame {
  return DGNS.DataFrame.fromColumns([
    f64Column('Time', acc.time),
    f64Column('Conc', acc.conc),
    f64Column('geoLower', acc.geoLower),
    f64Column('geoUpper', acc.geoUpper),
    DGNS.Column.fromStrings('Series', acc.series),
  ]);
}

/**
 * The shared curve→frame builder (F1, context-panel-ct-profile) for ONE compound curve: one row
 * per (series × point), columns `Time` / `Conc` / `geoLower` / `geoUpper` / `Series`. Consumed by
 * BOTH the grid thumbnail (`buildStudyGrid`) and the study-overview context panel — a SINGLE frame
 * source, so the panel corresponds to the thumbnail by construction (the bug this cycle fixes).
 * `seriesLabel` overrides the `Series` category (unused by the thumbnail; the overlay uses
 * {@link buildOverlayCurveFrame}). Post-E3 a curve is a single pooled series; its `group` is `''`
 * for a whole-analysis row → mapped to {@link THUMB_SERIES_DEFAULT} (an empty category renders
 * black on `plot.line`). Exported so both consumers share it (rule 5; correspondence by construction).
 */
export function buildCurveFrame(curve: CompoundCurve, seriesLabel?: string): DG.DataFrame {
  const acc: CurveFrameCols = {time: [], conc: [], geoLower: [], geoUpper: [], series: []};
  pushCurve(acc, curve, seriesLabel ?? null);
  return framifyCurve(acc);
}

/**
 * The F3 overlay assembler (context-panel-ct-profile): concatenate several selected rows' cached
 * curves into ONE frame, each series relabelled by its row identity (`Analysis · group · route`),
 * so a single semi-log axis carries every overlaid curve coloured/keyed by row. Reuses the SAME
 * {@link pushCurve}/{@link framifyCurve} emit as the single-curve {@link buildCurveFrame} — the
 * overlay and the single-row panel cannot drift in column shape. Post-E3 each row's curve is
 * exactly one series (LBC-2), so one entry contributes one line.
 */
export function buildOverlayCurveFrame(
  entries: ReadonlyArray<{curve: CompoundCurve; label: string}>,
): DG.DataFrame {
  const acc: CurveFrameCols = {time: [], conc: [], geoLower: [], geoUpper: [], series: []};
  for (const e of entries) {
    // AC-F3.2 (LBC-2): post-E3 each row's curve is EXACTLY ONE pooled series; a `series.length > 1`
    // curve is stale migrated data (the v8→v9 migration nulls those, so this cannot occur under
    // current paths — this is a defense-in-depth integrity assert). LOG the anomaly rather than
    // treat it as a normal multi-series case. Rendering is safe regardless: `pushCurve` forces the
    // row-identity `label` onto every point, so even a multi-series curve collapses to ONE labelled
    // line here — never silently rendered as multiple series.
    if (e.curve.series.length > 1) {
      console.warn(
        `buildOverlayCurveFrame: overlay row '${e.label}' carries ${e.curve.series.length} curve ` +
        'series (post-E3 expects exactly 1 — stale migrated data?); overlaying as ONE row-labelled ' +
        'series, not multiple. Data-integrity anomaly.');
    }
    pushCurve(acc, e.curve, e.label);
  }
  return framifyCurve(acc);
}
