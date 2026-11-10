/**
 * Results view — configures the per-subject results DataFrame for presentation
 * and applies the grid-side profile.
 *
 * **Column counts are never restated here.** The schema size is
 * `FIXED_RESULT_COLUMNS.length` (drift-guarded in the writer) and the coloured
 * set is `Object.keys(buildColourRules(...))`. This header used to say
 * "24-column" while the writer emitted 31, and three separate readers have
 * miscounted the banded set — hence the rule.
 *
 * **Background colour coding with default text.** Anomalies are highlighted
 * by tinting the CELL BACKGROUND; the text stays in Datagrok's default grid
 * colour (good readability on every tier). Three-tier pastel palette:
 * - ok (mint green)        — STATUS = ok, AND every parameter value in the
 *                            normal range (PCT_EXTRAP ≤ 20, adj R² > 0.90,
 *                            N pts > 3). Unified happy-path background.
 * - warning (peach)        — STATUS = partial, adj R² 0.85–0.90 (borderline),
 *                            N pts = 3.
 * - critical (soft red)    — STATUS = failed, PCT_EXTRAP > 20, adj R² < 0.85
 *                            (poor fit — λz unreliable).
 *
 * **Column header tooltips explain the colours.** Each coloured column has a
 * rich tooltip wired via `grid.onCellTooltip` that names the parameter and
 * lists its highlighting rules with colour swatches. No legend lives outside
 * the table — the rule is discoverable at the column header.
 *
 * **JSON audit columns hidden by default.** WARNINGS_JSON + PROVENANCE_JSON
 * are audit-trail surfaces; the user can re-enable them via Datagrok's
 * column-chooser (right-click on a column header). All 22 columns remain in
 * the DataFrame for downstream `.d42` round-trip and export.
 *
 * **Per-row CL/Vz basis (per-row-route AD-4).** The CL/VZ numeric columns now
 * carry **neutral** `column.tags.label` (`'CL'`/`'Vz'`) — the F-corrected
 * interpretation is the per-row `CL_BASIS` STRING column (`'CL'`/`'CL/F'`). The
 * column tooltips direct the reader to `CL_BASIS`; the old all-or-nothing
 * `'CL/F'`/`'Vz/F'` frame-level labels are removed.
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';
// Value namespace — needed to CREATE the acceptance columns
// ({@link writeAcceptanceColumns}). The type-only `DG` alias above stays for the
// type positions; `DGNS` mirrors the writer / grouping-columns idiom
// (`Column.fromList` + `COLUMN_TYPE`).
import * as DGNS from 'datagrok-api/dg';

import {DEFAULT_RULES, DEFAULT_ACCEPTANCE_RSQUARED, numAt, isLzServiceColumn,
  formatOverrideTimestamp, evaluateAcceptance, fmtThreshold, SPAN_RATIO_POOR,
  ACCEPT_ACCEPTED, ACCEPT_FLAGGED, ACCEPT_NA, isLzDependentParam, formatReasons,
  REVIEW_ACTION_LABEL, isSparseLzDependentCol, REVIEW_COL, REVIEW_EXCLUDED,
  LZ_REVIEW_INFO_COL, REVIEW_TO_REVIEW, upsertStringColumn,
  ACCEPT_REASONS_COL, ACCEPTANCE_FLAG_CODES, criterionShortLabel, reasonsCellValue, ACCEPTANCE_COL, FLAG_COLS,
  tagReasonPillsColumn, BW_UNIT} from '../adapter';
import type {
  AcceptanceReading, AcceptanceFlagCode, FlagState, ReviewDecisionRecord, PipelineFlag, AnalysisOrigin, ReviewInfo,
} from '../adapter';
import {richTipEl} from './rich-tooltip';
import type {RichTipSpec} from './rich-tooltip';
import {STATUS_ACCEPTED, STATUS_FLAGGED, STATUS_NOT_ESTIMATED} from './status-palette';

/** ARGB integers (Datagrok colour format — alpha-first 32-bit). Pastel-tier
 *  background colours; default text on every tier gives high readability.
 *  Most `setConditional` rulesets cover the entire numeric domain — values
 *  outside an anomaly band get the OK background explicitly.
 *
 *  **`LAMBDAZ_NPOINTS` deliberately does NOT** (2026-07-30): its sub-minimum band
 *  was removed because the only value that reaches it is `DG.INT_NULL`, and
 *  painting a cell that holds no value assigns it a severity it cannot have. Full
 *  domain coverage is a convention here, never a DG requirement — treat it as the
 *  default, not an invariant, and drop it wherever a band's only occupant is a
 *  null sentinel. See the band list in `buildColourRules`. */
const COLOR_OK = 0xFFDCF3E7; // mint green
const COLOR_WARNING = 0xFFFFECDB; // peach / soft yellow
const COLOR_CRITICAL = 0xFFEB6767; // soft red

// The per-subject results grid colours the cell TEXT instead of the background FILL (cleaner to
// read). The text uses the SAME verdict
// hues the facet tiles use for their lines/markers ({@link buildFacetVerdictColors}), so the
// grid reads in the same colours as the plots — declared ONCE here and referenced by BOTH, so
// grid and facet cannot drift. The band LOGIC is unchanged:
// configureResultsGrid's onCellPrepare finds a cell's tier through the SAME buildColourRules
// bands and maps the matched fill → its verdict hue via `fillToTextColour`, while
// configureResultsDataFrame no longer applies the background bands. Legend + tooltips (from
// buildColourRules directly) are untouched.
//
// The Radix "teal-success" verdict triad — the SINGLE source is `status-palette.ts`
// (user-directed 2026-09-09, superseding the DG categorical-palette sourcing of 2026-09-02).
// Grid text, the facet lines/markers, the header-tooltip legend swatches, and the λz-editor
// header ALL read these three constants (grid text via `fillToTextColour`), so nothing can
// drift from the palette. See status-palette.ts for the rationale (equal-lightness, CVD-safe
// teal instead of green) and the source report.
const HUE_ACCEPTED = STATUS_ACCEPTED; // Radix teal-11 #008573
const HUE_FLAGGED = STATUS_FLAGGED; // Radix amber-11 #AB6400
const HUE_NOT_ESTIMATED = STATUS_NOT_ESTIMATED; // Radix red-11 #CE2C31
/** Muted text of an excluded (not reportable) λz-dependent value — the neutral grey-4 used for
 *  muted text elsewhere in the app (no severity colour: excluded is a decision, not a verdict). */
export const EXCLUDED_TEXT = 0xFF9497A0;
/** Link blue (`--blue-1`) of a clickable "To review" cell. */
const REVIEW_LINK_TEXT = 0xFF2083D5;
/** Map a banded cell's FILL tier (COLOR_OK/WARNING/CRITICAL) → the verdict hue the facet uses,
 *  so the results-grid cell TEXT matches the plots (user-directed 2026-09-02). */
function fillToTextColour(fill: number): number {
  return fill === COLOR_WARNING ? HUE_FLAGGED :
    fill === COLOR_CRITICAL ? HUE_NOT_ESTIMATED : HUE_ACCEPTED;
}

/** Prepend `bold` to a grid cell's font, PRESERVING its size/family. ONE idiom shared by the
 *  results grid (STATUS/ACCEPTANCE) and the aggregated grid's fragile-N cells (nca-mapping-pane.ts
 *  imports this — no re-typed twin, per rule 5). Idempotent (never double-bolds); falls back to a
 *  sane default when the base font is empty. */
export function boldGridCell(gc: DG.GridCell): void {
  const base = gc.style.font;
  gc.style.font = (base !== null && base !== undefined && base.length > 0) ?
    (/bold/i.test(base) ? base : `bold ${base}`) : 'bold 13px sans-serif';
}

/** Numeric thresholds that drive both the cell-background colour rules and
 *  the column-header tooltip legends. Sourced from the four method-settings
 *  inputs that are currently surface-visible in the mapping pane's `Rules`
 *  section (`extrapWarnPct`, `extrapErrorPct`, `λz min adj-R²`,
 *  `λz min points`) so the grid presentation tracks the user-edited
 *  rules — not a stale set of hardcoded constants. */
export interface ResultsColorThresholds {
  readonly extrapWarnPct: number;
  readonly extrapErrorPct: number;
  /** The sci-comp point-selection GATE (`nca.LambdaZStrategy.minRSquared`) — the
   *  minimum fit quality the core admits during λz window selection. Projected from
   *  `ncaRules.lambdaZ.minRSquared`. It is NO LONGER the acceptance threshold (that is
   *  `lambdaAcceptanceRSquared`); it survives here only as the reference gate the
   *  MRT/Vss fragility tint reads and as a validation floor for `R_acc`. */
  readonly lambdaMinRSquared: number;
  /** The λz adj-R² ACCEPTANCE threshold (`R_acc`, lambdaz-adjr2-threshold-spec-v2 §3.2):
   *  a finished fit with `adjR² < lambdaAcceptanceRSquared` is Flagged. The value is
   *  applied VERBATIM — no derived recomputation (the removed `goodR`).
   *
   *  **The ONE field of this interface that is NOT a projection from `nca.NcaRules`.**
   *  Every sibling is sourced from `ncaRules.*` (the AD-4 invariant below), but `R_acc`
   *  is an nca-studio SAP/presentation threshold with NO sci-comp counterpart (thesis #1),
   *  so it is sourced from its own mapping-pane input, not from the rules object. There is
   *  therefore no `ncaRules.lambdaZ.acceptanceRSquared` a maintainer could accidentally
   *  wire — the exact AD-4 collision vector is structurally absent for this field. */
  readonly lambdaAcceptanceRSquared: number;
  readonly lambdaMinPoints: number;
  /** %AUMCextrap soft-warning threshold (FR-200). Drives the PCT_AUMC_EXTRAP
   *  band and the MRT/VSS fragility tint. House heuristic — no error tier. */
  readonly extrapWarnPctAumc: number;
}

/** Default thresholds derived from {@link DEFAULT_RULES}. Used by
 *  {@link openResultsView} and as the default argument of
 *  {@link configureResultsDataFrame} / {@link configureResultsGrid} so
 *  call sites that have no Rules context still get a sensible legend. */
export const DEFAULT_RESULTS_THRESHOLDS: ResultsColorThresholds = {
  extrapWarnPct: DEFAULT_RULES.extrapWarnPct,
  extrapErrorPct: DEFAULT_RULES.extrapErrorPct,
  lambdaMinRSquared: DEFAULT_RULES.lambdaZ.minRSquared,
  // NOT a DEFAULT_RULES projection — the acceptance threshold has no sci-comp
  // counterpart (see the interface field's contract). Its own nca-studio default.
  lambdaAcceptanceRSquared: DEFAULT_ACCEPTANCE_RSQUARED,
  lambdaMinPoints: DEFAULT_RULES.lambdaZ.minPoints,
  extrapWarnPctAumc: DEFAULT_RULES.extrapWarnPctAumc,
};

/**
 * λz span-ratio display thresholds, in half-lives (RG-C, user-ratified
 * 2026-08-11): `< 2` poor · `2–4` borderline · `≥ 4` adequate.
 *
 * **AD-4 — these are deliberately NOT fields of {@link ResultsColorThresholds},
 * and that omission is a structural guard rather than an oversight.** Every field
 * of that interface is a *projection from* `nca.NcaRules` — an
 * invariant the code states at its construction site (`nca-mapping-pane.ts`,
 * "Project from the just-built NcaRules … per the Rules→presentation projection
 * invariant") and that `DEFAULT_RESULTS_THRESHOLDS` follows below — with ONE
 * documented exception, `lambdaAcceptanceRSquared` (`R_acc`): it is an nca-studio
 * SAP threshold with no `nca.NcaRules` counterpart, so it is sourced from its own
 * mapping-pane input, not projected. That exception is SAFE from this very failure
 * mode precisely because there is no `ncaRules.lambdaZ.acceptanceRSquared` to
 * accidentally wire — see the field's own contract. The span thresholds below are
 * different: a `minSpanRatio` DOES exist on `nca.NcaRules`, so keeping them as
 * module constants is what prevents the wire-up. That
 * invariant is exactly the collision vector: adding a span field there would
 * oblige the next maintainer, *correctly following the established rule*, to
 * source it from `ncaRules.lambdaZ.minSpanRatio` — thereby setting sci-comp's
 * `minSpanRatio` as a side effect of tidying.
 *
 * That would not change a single number (sci-comp's contract: the threshold is
 * read at one site, only to push a `LAMBDAZ_LOW_SPAN` warning, and is absent
 * from λz window selection — "emitting the warning is all it does; it never
 * affects the fit"). What it WOULD do is inject unrequested warnings into
 * `provenance.warnings` → `WARNINGS_JSON` → the compute-notices banner, and
 * corrupt the audit log's `ruleHash` so the audit trail asserts a computation
 * rule the scientist never set. That is a correctness-of-provenance defect, not
 * a parity risk — and keeping these constants module-local means there is no
 * field to project, so the invariant that would have caused it now prevents it.
 *
 * Neither constant is named `minSpanRatio`, for the same reason (AC-3.5).
 */
// `SPAN_RATIO_POOR` (= 2) now lives in `adapter/acceptance.ts` beside the SPAN criterion it
// also drives (nca-workflow-stages AD-4) — same AD-4 guard, imported above.
const SPAN_RATIO_BORDERLINE = 4;
/** DG has no open-ended lower band, so `≥ 4` needs an upper sentinel. 100
 *  half-lives is beyond anything physically observable (the widest span in our
 *  entire reference corpus is 10.94), and — like the `AUC_DF` band's `upTo` —
 *  it is a rendering artifact that `legendText` deliberately never shows, where
 *  it would read as a real ceiling. */
const SPAN_RATIO_LEGEND_CEILING = 100;

// ──────────────────────────────────────────────────────────────────────────
// λz acceptance verdict + flags (nca-grid-status-spec).
//
// A single ROLL-UP verdict answering "does this profile need an analyst
// decision?", plus four per-criterion flag columns for filter / count. The
// STATUS column answers a DIFFERENT question — "did the algorithm run?" — and
// on a real frame reads `ok` for every subject while span-ratio 0.69, adj-R²
// 0.85–0.90 and %AUMCextrap tints mark profiles that clearly need review. That
// gap (thesis #4: surface honest uncertainty) is what these columns close.
//
// **Thresholds TRACK LIVE RULES** — the flags read the SAME `ResultsColorThresholds`
// object and the SAME module-local `SPAN_RATIO_POOR` that drive the colour bands,
// so a flag and its column's band can never disagree (user decision).
//
//   • `RSQ_FLAG`: fires at the λz adj-R² ACCEPTANCE threshold `R_acc`
//     (`lambdaAcceptanceRSquared`), the value the user enters, applied VERBATIM —
//     `adjR² < R_acc` → Flagged (lambdaz-adjr2-threshold-spec-v2 §4). `R_acc` is a
//     SEPARATE field from the core point-selection gate `minRSquared`, and is NOT
//     derived from it: the earlier `goodR = minR + …` formula that conflated the two
//     roles was removed (spec §2), so the applied criterion equals what was set. The
//     `R_acc ≥ minR` consistency constraint (§3.2) is enforced at input.
//   • `NEXT_FLAG`: spec `N < 3`, but sci-comp REJECTS a window below `minPoints`, so a
//     literal `< 3` flag is structurally dead (an accepted fit can never trip it). Uses
//     `≤ minPoints`, matching {@link readFragility} and the peach `N = minP` band (a fit
//     AT the floor is fragile).
// ──────────────────────────────────────────────────────────────────────────

/** The results-frame column names for the verdict + the four flag columns (defined in the
 *  adapter, `acceptance.ts`, so the passthrough gate can reserve them; re-exported for the views). */
export {ACCEPTANCE_COL, FLAG_COLS};
/** The `Violated` pills column beside ACCEPTANCE (acceptance-reason-pills) — written in the SAME
 *  pass as the verdict, so a pill and the verdict it explains cannot disagree. */
export {ACCEPT_REASONS_COL};

/** Hidden helper column carrying the problem-first sort rank (0 = most severe,
 *  leads the grid). NEVER shown in either grid profile ({@link resultsGridVisibleOrder}
 *  skips it, {@link applyResultsColumnProfile} then paints `visible = false`), and
 *  written PRESENTATION-side in the SAME pass as ACCEPTANCE ({@link writeAcceptanceColumns}),
 *  from the same thresholds — so the rank and the verdict it encodes cannot drift.
 *  `_`-prefixed to read as internal; excluded from any future frame-level export the
 *  way the JSON audit columns are. */
export const PROBLEM_RANK_COL = '_PROBLEM_RANK';

// ──────────────────────────────────────────────────────────────────────────
// λz review state (nca-workflow-stages F-S4.1 / F-S4.5; D15; approved design 2026-09-26).
// ──────────────────────────────────────────────────────────────────────────

/** The `Review` column: one neutral word per profile — whether the analyst still has to look at
 *  it and, if they decided, what (D15). Written PRESENTATION-side beside ACCEPTANCE (not a writer
 *  column, so `FIXED_RESULT_COLUMNS` is untouched). NOT colour-banded: the verdict colour already
 *  lives in ACCEPTANCE, and review state is a separate axis. Defined in the adapter (the study
 *  summary cache reads it too) and re-exported here for the views. */
export {REVIEW_COL, REVIEW_EXCLUDED};

/** Frame tag: `'true'` when the Review column is visible in the COMPACT profile (Review λz
 *  stage). Set by the view on every stage change; read by {@link resultsGridVisibleOrder} so a
 *  Compact/Full toggle can never undo the stage scoping. Full always shows the column. */
export const REVIEW_COMPACT_TAG = '.nca-review-compact';

/** Frame tag: `'false'` when the scientist switched the `Violated` pills column off (Settings).
 *  Read by {@link resultsGridVisibleOrder}, so the column stays hidden in BOTH profiles and a
 *  Compact/Full toggle cannot bring it back. Any other value leaves the usual rule in force. */
export const REASON_PILLS_TAG = '.nca-reason-pills';

/** Show / hide the `Violated` column of the Sparse NCA grid: the Settings switch AND at least one
 *  Flagged composite — the same "say it only where it can" rule as the per-subject compact grid.
 *  Grid-side only. */
export function applySparseReasonPills(grid: DG.Grid, show: boolean): void {
  const gc = grid.col(ACCEPT_REASONS_COL);
  const df = grid.dataFrame;
  if (gc !== null) gc.visible = show && df != null && anyReason(df);
}

// ── Pipeline reasons (nca-pipeline-save-rework F10, approved design 2026-10-15) ──────────

/** The `Pipeline` column: the pipeline run's OWN reasons for a profile (the ones the λz acceptance
 *  check does not raise — `PipelineFlag.pipelineOnly`). Neutral text, no hue (the verdict colour
 *  stays on ACCEPTANCE). Present only when the analysis carries at least one such reason. */
export const PIPELINE_COL = 'PIPELINE';
/** Hidden `_lz*` detail of the `Pipeline` cell (JSON {@link PipelineCellInfo}) for its tooltip. */
export const LZ_PIPELINE_INFO_COL = '_lz_pipeline';

/** The pipeline run behind an analysis, as the view shows it: who/when + the reasons per unit key
 *  (`SUBJECT`, or the composite unit key `sparse:<group>`). */
export interface PipelineReasonsView {
  readonly runAt: string;
  readonly user: string;
  readonly byKey: ReadonlyMap<string, ReadonlyArray<string>>;
}

/** One `Pipeline` cell's tooltip detail. */
interface PipelineCellInfo {
  readonly runAt: string;
  readonly user: string;
  readonly reasons: ReadonlyArray<string>;
}

/** The view model of an entry's pipeline flags, or null when it carries no pipeline-only reason
 *  (then no column, no pane line — nothing the app does not already show). */
export function pipelineReasonsViewOf(
  flags: ReadonlyArray<PipelineFlag> | undefined, origin: AnalysisOrigin | null | undefined,
): PipelineReasonsView | null {
  if (origin?.kind !== 'pipeline') return null;
  const byKey = new Map<string, ReadonlyArray<string>>();
  for (const f of flags ?? []) if ((f.pipelineOnly ?? []).length > 0) byKey.set(f.profileKey, f.pipelineOnly);
  return byKey.size === 0 ? null : {runAt: origin.runAt, user: origin.user, byKey};
}

/** The `Pipeline` cell text: the first reason, `(+N)` for the rest. */
export function pipelineCellText(reasons: ReadonlyArray<string>): string {
  return reasons.length === 0 ? '' : reasons[0] + (reasons.length > 1 ? ` (+${reasons.length - 1})` : '');
}

/** The Review-pane line for one unit: `Pipeline (run <local time>): <reasons · …>`. */
export function pipelineReasonsLine(pv: PipelineReasonsView, reasons: ReadonlyArray<string>): string {
  return `Pipeline (run ${formatOverrideTimestamp(pv.runAt)}): ${reasons.join(' · ')}`;
}

/** The `Pipeline` cell tooltip: the run (local time · user) and one bullet per reason. */
function pipelineTipSpec(info: PipelineCellInfo): RichTipSpec {
  const who = info.user === '' ? 'unknown user' : info.user;
  return {
    title: `Pipeline run ${formatOverrideTimestamp(info.runAt)} · ${who}`,
    summary: 'Reasons the pipeline sent this profile to review that the λz acceptance check does not raise.',
    extra: info.reasons.map((r) => ui.divText(`• ${r}`)),
  };
}

/** The `Pipeline` cell tooltip of frame row `i`, or null (blank cell / no detail). */
export function pipelineTipAt(df: DG.DataFrame | null, i: number): RichTipSpec | null {
  const c = df?.col(LZ_PIPELINE_INFO_COL) ?? null;
  if (c === null || c.isNone(i) || c.getString(i) === '') return null;
  try {
    return pipelineTipSpec(JSON.parse(c.getString(i)) as PipelineCellInfo);
  } catch {
    return null;
  }
}

/**
 * Write (or remove) the `Pipeline` column + its hidden detail on a results frame, row `i` keyed by
 * `keyAt(i)`. With no pipeline-only reason on this frame the columns are REMOVED, so the column
 * appears only when it carries something (the compact rule). Returns true when the frame's columns
 * changed (written — `upsertStringColumn` re-creates them — or removed); the caller then re-applies
 * the column profile.
 */
export function writePipelineColumns(
  df: DG.DataFrame, pv: PipelineReasonsView | null, keyAt: (i: number) => string,
): boolean {
  const n = df.rowCount;
  const texts: string[] = new Array(n).fill('');
  const infos: string[] = new Array(n).fill('');
  let any = false;
  if (pv !== null) {
    for (let i = 0; i < n; i++) {
      const reasons = pv.byKey.get(keyAt(i));
      if (reasons === undefined || reasons.length === 0) continue;
      texts[i] = pipelineCellText(reasons);
      infos[i] = JSON.stringify({runAt: pv.runAt, user: pv.user, reasons} satisfies PipelineCellInfo);
      any = true;
    }
  }
  const had = df.col(PIPELINE_COL) !== null;
  if (!any) {
    if (had) df.columns.remove(PIPELINE_COL);
    if (df.col(LZ_PIPELINE_INFO_COL) !== null) df.columns.remove(LZ_PIPELINE_INFO_COL);
    return had;
  }
  upsertStringColumn(df, PIPELINE_COL, texts);
  upsertStringColumn(df, LZ_PIPELINE_INFO_COL, infos);
  df.col(PIPELINE_COL)!.setTag('friendlyName', DISPLAY_HEADERS[PIPELINE_COL]);
  return true;
}

function reviewInfoAt(df: DG.DataFrame | null, i: number): ReviewInfo | null {
  const c = df?.col(LZ_REVIEW_INFO_COL) ?? null;
  if (c === null || c.isNone(i)) return null;
  const t = c.getString(i);
  if (t === '') return null;
  try {return JSON.parse(t) as ReviewInfo;} catch {return null;}
}

/** Who/when line of a decision, local time (the ✎ convention). */
function decidedLine(d: ReviewDecisionRecord): string {
  const who = d.user === '' ? 'unknown user' : d.user;
  return `${REVIEW_ACTION_LABEL[d.action]} — ${who} · ${formatOverrideTimestamp(d.at)}`;
}

/** The Exclude λz decision of results row `i`, or null — the ONE read of "is this profile
 *  excluded" for every per-subject surface that marks a value (the grid cell, the context-panel
 *  NCA parameters readout — AC-S5.1.4), so the two cannot disagree. */
export function excludedDecisionAt(df: DG.DataFrame | null, i: number): ReviewDecisionRecord | null {
  const info = reviewInfoAt(df, i);
  return info?.state === REVIEW_EXCLUDED ? info.decision : null;
}

/** The rich tooltip of an EXCLUDED λz-dependent value (F-S4.5) — ONE wording for the dense
 *  per-subject grid and the UC-04 composite table, so the two surfaces cannot drift. */
export function excludedCellTipSpec(decision: ReviewDecisionRecord, composite = false): RichTipSpec {
  return {
    title: 'excluded (λz)',
    summary: composite ?
      'Not reportable — this λz-dependent value of the pooled composite is excluded from the report.' :
      'Not reportable — this λz-dependent value is left out of every statistic ' +
      '(grouped summary, dose proportionality, %F, TK, study grid).',
    detail: `Reasons: ${formatReasons(decision.reasons)}`,
    note: decidedLine(decision),
  };
}

/** Paint one excluded λz-dependent cell (F-S4.5): the value kept, muted, marked ` (excl.)`,
 *  one line truncated at the END so the number stays whole. Display-only — the VALUE, sort,
 *  filter and band are untouched. Shared by the dense grid and the composite table. */
function paintExcludedCell(gc: DG.GridCell, text: string): void {
  gc.customText = `${text} (excl.)`;
  gc.style.textColor = EXCLUDED_TEXT;
  gc.style.textWrap = 'ellipsis';
}

/**
 * Widen the given columns so an excluded value reads WHOLE as `{value} (excl.)` (manual QA
 * 2026-09-30: in a narrow column DG's ellipsis cut the number itself — `).908…l.)`). Only widens,
 * never narrows; a no-op when no row is excluded. Width from the longest excluded cell's text at
 * the grid's ~7 px/char + padding. Shared by the per-subject grid and the composite table.
 */
export function fitExcludedColumns(
  grid: DG.Grid, cols: ReadonlyArray<string>, isExcluded: (row: number) => boolean,
): void {
  const df = grid.dataFrame;
  if (df == null) return;
  const rows: number[] = [];
  for (let i = 0; i < df.rowCount; i++) if (isExcluded(i)) rows.push(i);
  if (rows.length === 0) return;
  for (const name of cols) {
    const col = df.col(name);
    const gc = grid.columns.byName(name);
    if (col === null || gc === null) continue;
    let longest = 0;
    for (const i of rows) if (!col.isNone(i)) longest = Math.max(longest, `${col.getString(i)} (excl.)`.length);
    if (longest === 0) continue;
    const need = longest * 7 + 16;
    if (gc.width < need) gc.width = need;
  }
}

/** Rich tooltip of one Review cell (D15): the state, then why (category / reasons), then who and
 *  when. The Not-applicable BLQ-only wording is the AC-S1.3.1a disclosure, verbatim. */
export function buildReviewTooltipSpec(df: DG.DataFrame, i: number): RichTipSpec | null {
  const info = reviewInfoAt(df, i);
  if (info === null) return null;
  const lines: string[] = [];
  if (info.state === 'Not applicable') {
    lines.push('Not applicable — reviewed at the group composite.');
    if (info.blqOnly) {
      const times = [...new Set(info.blqTimes)].join(', ');
      lines.push(`This animal's sample(s) are BLQ at ${times === '' ? 'every sampled time' : times}.`);
    }
    return {title: 'Review', summary: lines.join(' ')};
  }
  if (info.state === 'Not needed')
    return {title: 'Review', summary: 'Not needed — the λz fit meets every acceptance criterion.'};
  const why = info.category === 'no-exposure' ?
    'No exposure — no measurable concentration (all BLQ / zero / missing).' :
    info.category === 'not-estimated' ? 'Not estimated — λz could not be fitted.' :
      'Flagged — at least one λz acceptance criterion is out of bounds.';
  if (info.decision === null)
    return {title: 'To review', summary: why, detail: 'Click to review this profile.'};
  return {
    title: info.state, summary: why,
    detail: `Reasons: ${formatReasons(info.decision.reasons)}`,
    note: decidedLine(info.decision),
  };
}

/** Results-frame service column carrying the per-subject λz manual-override
 *  indicator: the analyst label (from `_lz_user`) when the subject is under
 *  override, `''` otherwise. Hidden from BOTH grid profiles (an `_lz*` service
 *  column — {@link resultsGridVisibleOrder} skips it, {@link isLzServiceColumn}),
 *  read only by the ACCEPTANCE cell's ✏️ glyph + tooltip (F7 / AC-7.1). Lives on
 *  the RESULTS frame keyed by SUBJECT — NOT the source df where the `_lz_action`
 *  truth lives — because the render hooks bind to the results grid; the controller
 *  projects one onto the other via {@link writeSubjectOverrideFlags}. A cell is
 *  non-empty iff the subject is overridden; an unknown analyst renders as
 *  {@link OVERRIDE_UNKNOWN_USER} so the ✏️ still shows. */
export const LZ_OVERRIDDEN_COL = '_lz_overridden';

/** Results-frame service column carrying the per-subject λz-override TIMESTAMP (raw
 *  `_lz_time` ISO of the same representative point {@link LZ_OVERRIDDEN_COL} takes the
 *  analyst from), projected by {@link writeSubjectOverrideFlags} / {@link setSubjectOverrideFlag}.
 *  Read only by the ACCEPTANCE cell's ✎ tooltip (the gray who/when line), formatted at
 *  read time via `formatOverrideTimestamp`. An `_lz*` service column — hidden from both
 *  grid profiles by the same {@link isLzServiceColumn} path as {@link LZ_OVERRIDDEN_COL}. */
export const LZ_OVERRIDDEN_TIME_COL = '_lz_overridden_time';

/** Fallback analyst label for an override whose `_lz_user` was blank (no shell user
 *  resolvable when it was made) — keeps the flag cell non-empty so the ✏️ shows. */
export const OVERRIDE_UNKNOWN_USER = 'unknown user';

/** ACCEPTANCE verdict vocabulary (closed set, spec §2). `Not estimated` is the
 *  third, spec-implied state: when λz was not fitted every criterion is `n/a`, so
 *  the profile can be neither accepted nor flagged. It is the MOST severe outcome —
 *  none of the terminal-phase parameters (t½, CL, Vz, AUC₀–∞, MRT, Vss) exist — so
 *  it is banded CRITICAL (red), distinct from `Flagged` (params exist but questionable). */
// Verdict strings ACCEPT_ACCEPTED / ACCEPT_FLAGGED / ACCEPT_NA live in adapter/acceptance.ts
// (imported above) so the non-UI review classifier reads the same vocabulary.

/** Per-criterion flag glyphs (spec §3 tri-state). DISTINCT strings, deliberately:
 *  the spec makes "distinguishing `false` from `n/a`" mandatory, and a nullable
 *  DG bool column blurs exactly that (a null bool and a `false` bool render alike).
 *  A string column with three glyphs renders them distinctly and still filters /
 *  counts as three categories. `✗` = criterion violated (the eye-catching one, a
 *  problem), `✓` = met, `n/a` = not evaluated (λz not estimated). */
const FLAG_VIOLATED = '✗';
const FLAG_MET = '✓';
const FLAG_NA = 'n/a';

// ──────────────────────────────────────────────────────────────────────────
// Colour bands — ONE source for the DG rule and the tooltip legend.
//
// A tinted cell whose header tooltip does not explain the tint is a defect
// (MRT / Vss shipped that way; user QA 2026-07-30). It happened because the
// colour rule and the legend that describes it were written in two places,
// each restating the threshold as its own literal. Everything below exists so
// a threshold is written ONCE, as a number, and BOTH the DG range key
// (`'<=20'`) and the human legend (`'≤ 20% — …'`) are derived from it.
// ──────────────────────────────────────────────────────────────────────────

/** A band's numeric predicate. `none` = the band is painted by a cross-column
 *  `onCellPrepare` hook, so it has no per-column DG rule — only a legend. */
type BandOp =
  | {readonly op: 'le' | 'lt' | 'gt' | 'eq'; readonly v: number}
  /** `≥ v` to the scientist, `v-upTo` to DG — which has no open-ended lower
   *  band, so an upper sentinel is required. `upTo` is a rendering artifact and
   *  is deliberately NOT shown in the legend: "≥ 3 — adequate" is the rule, "3
   *  to 100" is an implementation detail that would read as a real ceiling. */
  | {readonly op: 'ge'; readonly v: number; readonly upTo: number}
  | {readonly op: 'range'; readonly lo: number; readonly hi: number}
  | {readonly op: 'cat'; readonly v: string}
  | {readonly op: 'none'};

/** The thresholds a band asserts to the SCIENTIST — the numbers that must appear
 *  in both the DG rule and the legend. Excludes `ge.upTo` (a DG-only sentinel).
 *  Drives the correspondence test in results-view-presentation-tests. */
function semanticThresholds(b: BandOp): number[] {
  switch (b.op) {
  case 'le': case 'lt': case 'gt': case 'eq': case 'ge': return [b.v];
  case 'range': return [b.lo, b.hi];
  case 'cat': case 'none': return [];
  }
}

/** One colour band: the predicate, the fill, and the prose after the em-dash. */
interface ColourBand {
  readonly band: BandOp;
  readonly color: number;
  /** Meaning, e.g. `'acceptable extrapolation'`. The threshold itself is NOT
   *  repeated here — it is rendered from {@link BandOp} by {@link legendText}. */
  readonly gloss: string;
  /** Unit suffix on the threshold in the legend (`'%'`); omitted for bare numbers. */
  readonly unit?: string;
  /** Rule-only band: registered with DG but NOT listed in the legend. For
   *  DEFENSIVE bands that paint a value the engine cannot emit (e.g. λz points
   *  BELOW the configured minimum — sci-comp rejects such fits) and so carry no
   *  meaning a scientist needs to read. Never use this to hide a band a real
   *  value can land in: that is the defect this whole file guards against. */
  readonly hidden?: boolean;
}

/** A column's colour rule: how it is painted, and its bands. */
interface ColourRule {
  /** `conditional` → numeric `setConditional`; `categorical` → `setCategorical`;
   *  `crossColumn` → painted by {@link configureResultsGrid}'s `onCellPrepare`
   *  from SIBLING columns, so `configureResultsDataFrame` registers nothing. */
  readonly kind: 'conditional' | 'categorical' | 'crossColumn';
  readonly bands: ReadonlyArray<ColourBand>;
}

/** True when `v` falls in `op`'s predicate — the READ-side mirror of {@link dgRangeKey}
 *  (which WRITES the DG rule). Used by the results grid's text-colour `onCellPrepare` to find
 *  a cell's tier without re-declaring any threshold (they stay single-sourced in
 *  {@link buildColourRules}). `none` never matches — it is only a cross-column marker. */
function bandContains(op: BandOp, v: unknown): boolean {
  if (op.op === 'cat') return typeof v === 'string' && v === op.v;
  if (op.op === 'none') return false;
  if (typeof v !== 'number' || !Number.isFinite(v)) return false;
  switch (op.op) {
  case 'le': return v <= op.v;
  case 'lt': return v < op.v;
  case 'gt': return v > op.v;
  case 'eq': return v === op.v;
  case 'ge': return v >= op.v && v <= op.upTo;
  case 'range': return v >= op.lo && v <= op.hi;
  }
}

/** The FILL colour of the FIRST band `v` lands in, or null (unbanded value / no match). One
 *  evaluator over {@link buildColourRules}' bands, so the text tier a cell gets and the legend
 *  line it reads come from the same source and cannot disagree. */
function matchBandColour(rule: ColourRule, v: unknown): number | null {
  for (const cb of rule.bands)
    if (bandContains(cb.band, v)) return cb.color;
  return null;
}

/** The TEXT colour the results grid paints `value` in for `column` — the SAME band → text
 *  mapping the grid's `onCellPrepare` uses — or null (column not banded, no band matches, or a
 *  cross-column tint). Exported so the λz context-panel readout colours its criterion values
 *  EXACTLY as the grid does (user-directed 2026-09-02). */
export function gridBandTextColour(
  column: string, value: unknown,
  thresholds: ResultsColorThresholds = DEFAULT_RESULTS_THRESHOLDS,
): number | null {
  const rule = buildColourRules(thresholds)[column];
  if (rule === undefined || rule.kind === 'crossColumn') return null;
  const fill = matchBandColour(rule, value);
  return fill === null ? null : fillToTextColour(fill);
}

// `fmtThreshold` (BOTH the DG range key and the legend text; strips IEEE-754 noise such as
// `0.9500000000000001`) now lives in adapter/acceptance.ts, shared with the acceptance
// reasons — imported above.

/** The DG conditional/categorical rule key for a band, or `null` when the band
 *  carries no per-column rule (cross-column tint). */
function dgRangeKey(b: BandOp): string | null {
  switch (b.op) {
  case 'le': return `<=${fmtThreshold(b.v)}`;
  case 'lt': return `<${fmtThreshold(b.v)}`;
  case 'gt': return `>${fmtThreshold(b.v)}`;
  case 'eq': return fmtThreshold(b.v);
  case 'ge': return `${fmtThreshold(b.v)}-${fmtThreshold(b.upTo)}`;
  case 'range': return `${fmtThreshold(b.lo)}-${fmtThreshold(b.hi)}`;
  case 'cat': return b.v;
  case 'none': return null;
  }
}

/** The legend line for a band: `"≤ 20% — acceptable extrapolation"`. A
 *  cross-column band has no threshold prefix of its own — its gloss states the
 *  sibling conditions in full. */
function legendText(cb: ColourBand): string {
  const u = cb.unit ?? '';
  const b = cb.band;
  const prefix =
    b.op === 'le' ? `≤ ${fmtThreshold(b.v)}${u}` :
      b.op === 'lt' ? `< ${fmtThreshold(b.v)}${u}` :
        b.op === 'gt' ? `> ${fmtThreshold(b.v)}${u}` :
          b.op === 'ge' ? `≥ ${fmtThreshold(b.v)}${u}` :
            b.op === 'eq' ? `= ${fmtThreshold(b.v)}${u}` :
              b.op === 'range' ? `${fmtThreshold(b.lo)}–${fmtThreshold(b.hi)}${u}` :
                b.op === 'cat' ? b.v : null;
  return prefix === null ? cb.gloss : `${prefix} — ${cb.gloss}`;
}

/** Build the DG colour-rule object (`range key → argb`) for a rule's bands.
 *  Bands with no key (cross-column) are skipped. */
function dgColorRules(rule: ColourRule): Record<string, number> {
  const out: Record<string, number> = {};
  for (const cb of rule.bands) {
    const key = dgRangeKey(cb.band);
    if (key !== null) out[key] = cb.color;
  }
  return out;
}

/** One colour band, flattened into its two projections for inspection. */
export interface ColourBandProjection {
  /** Which grid the band belongs to — `results` (live Rules thresholds) or one of
   *  the {@link FIXED_GRID_SPECS} tabs (own fixed thresholds; see TD-06). */
  readonly grid: 'results' | ConfigurableGrid;
  readonly column: string;
  /** What DG is told to paint. `null` for a cross-column band, which has no
   *  per-column rule (painted by `onCellPrepare` from sibling columns). */
  readonly ruleKey: string | null;
  /** What the scientist is told it means. `null` for a hidden defensive band. */
  readonly legend: string | null;
  /** The numbers this band asserts to the scientist — must appear in BOTH
   *  projections. Excludes DG-only sentinels (see {@link semanticThresholds}). */
  readonly thresholds: ReadonlyArray<number>;
}

/**
 * Every colour band in the package, across all three grids that have one.
 *
 * Exported so a pure test can assert the two projections agree — the invariant
 * the single source exists to hold, and the one nothing checked when MRT / Vss
 * shipped tinted with a silent tooltip. Needs no DG Grid.
 */
export function colourBandProjections(
  thresholds: ResultsColorThresholds = DEFAULT_RESULTS_THRESHOLDS,
): ReadonlyArray<ColourBandProjection> {
  const out: ColourBandProjection[] = [];
  const push = (
    grid: ColourBandProjection['grid'], column: string, rule: ColourRule,
  ): void => {
    for (const cb of rule.bands) {
      out.push({
        grid, column,
        ruleKey: rule.kind === 'crossColumn' ? null : dgRangeKey(cb.band),
        legend: cb.hidden === true ? null : legendText(cb),
        thresholds: semanticThresholds(cb.band),
      });
    }
  };
  for (const [column, rule] of Object.entries(buildColourRules(thresholds)))
    push('results', column, rule);
  // EVERY spec-driven grid, from the one registry — not a hand-picked subset.
  // Only the two sparse tabs declare bands today; iterating all five means a
  // `colour` added to the F% or dose-proportionality spec tomorrow is covered by
  // this guard on the day it lands, with no second list to remember to update.
  for (const [gridId, spec] of Object.entries(FIXED_GRID_SPECS)) {
    for (const [column, d] of Object.entries(spec))
      if (d.colour !== undefined) push(gridId as ConfigurableGrid, column, d.colour);
  }
  return out;
}

/** The legend rows for a rule, ready to drop into `RichTipSpec.extra`. The swatch shows the
 *  SAME verdict hue the cell TEXT now uses ({@link fillToTextColour}, not the pale fill), so
 *  the header-tooltip legend matches what the scientist sees in the cells (text-colour trial,
 *  user-directed 2026-09-02). */
function legendRows(rule: ColourRule): HTMLElement[] {
  const shown = rule.bands.filter((cb) => cb.hidden !== true);
  return [legendSpacer(), ...shown.map((cb) => legendRow(fillToTextColour(cb.color), legendText(cb)))];
}

/**
 * Every coloured column of the per-subject results grid, as ONE declaration
 * consumed by three places: {@link configureResultsDataFrame} (registers the DG
 * rules), {@link buildHeaderTooltips} (renders the legend), and
 * {@link missingPresentationEntries} (asserts no column is coloured without a
 * legend). Rebuilt per call because every threshold is live — the user edits
 * them in the mapping pane's Rules section and re-runs Compute.
 *
 * `MRT` / `VSS` are `crossColumn`: their tint is decided by SIBLING columns
 * (%AUMCextrap, adj. R², N points), so no per-column DG rule can express it and
 * `configureResultsGrid` paints them in `onCellPrepare`. They still declare a
 * band here — that is precisely how their legend stops going missing.
 */
function buildColourRules(
  t: ResultsColorThresholds,
): Readonly<Record<string, ColourRule>> {
  const wp = t.extrapWarnPct;
  const ep = t.extrapErrorPct;
  const minR = t.lambdaMinRSquared;
  // The λz adj-R² ACCEPTANCE threshold (R_acc) — the single line the band and the
  // RSQ acceptance flag both read (lambdaz-adjr2-threshold-spec-v2 §3.2/§4). Applied
  // verbatim; there is no derived `goodR` anymore (§2).
  const rAcc = t.lambdaAcceptanceRSquared;
  const minP = t.lambdaMinPoints;
  const wpa = t.extrapWarnPctAumc;

  // The one sentence MRT and Vss share. Reads the SAME thresholds object the
  // `onCellPrepare` predicate reads, so the text cannot describe a rule the
  // code does not apply.
  const fragileGloss =
    `fragile — %AUMCextrap > ${fmtThreshold(wpa)}%, ` +
    `adj. R² < ${fmtThreshold(minR)}, or N points ≤ ${fmtThreshold(minP)}`;
  // Single-band rule: `onCellPrepare` paints ONLY the fragile branch, leaving
  // every other cell at the grid's default background. A COLOR_OK row here
  // would claim a mint fill that is never drawn.
  const crossColumnFragile: ColourRule = {
    kind: 'crossColumn',
    bands: [{band: {op: 'none'}, color: COLOR_WARNING, gloss: fragileGloss}],
  };

  return {
    STATUS: {kind: 'categorical', bands: [
      {band: {op: 'cat', v: 'ok'}, color: COLOR_OK,
        gloss: 'all parameters computed successfully'},
      {band: {op: 'cat', v: 'partial'}, color: COLOR_WARNING,
        gloss: 'results produced but some parameters could not be derived'},
      {band: {op: 'cat', v: 'failed'}, color: COLOR_CRITICAL,
        gloss: 'computation did not produce results'},
    ]},

    // The λz ACCEPTANCE roll-up (nca-grid-status-spec §2). Categorical, like
    // STATUS — and banding it is what puts it in the compact profile: the
    // mandatory-visible floor is DERIVED from `Object.keys(buildColourRules(...))`
    // (AD-2), so a verdict the scientist must act on is compact-visible with no
    // edit to the profile logic. ALL THREE states are banded: `Not estimated`
    // (`ACCEPT_NA`) is CRITICAL/red on purpose — a λz-less profile has NONE of the
    // terminal-phase parameters, the most severe outcome, so the scientist must SEE
    // it, not read a blank cell (an earlier draft left it unbanded/blank; changed on
    // user request). Glosses hold NO digits, so the legend-vs-rule number-
    // correspondence guard passes vacuously for a categorical band.
    ACCEPTANCE: {kind: 'categorical', bands: [
      {band: {op: 'cat', v: ACCEPT_ACCEPTED}, color: COLOR_OK,
        gloss: 'all λz acceptance criteria met'},
      {band: {op: 'cat', v: ACCEPT_FLAGGED}, color: COLOR_WARNING,
        gloss: 'one or more λz acceptance criteria violated — analyst review needed'},
      {band: {op: 'cat', v: ACCEPT_NA}, color: COLOR_CRITICAL,
        gloss: 'λz not fitted — no terminal-phase parameters, criteria not evaluated'},
    ]},

    // Three-tier band driven by both Rules → Extrapolation inputs.
    //
    // Edge case: `wp == ep` produces a single-point middle range (e.g. "20-20").
    // DG renders this as an effectively-empty band — values at exactly `wp` fall
    // to CRITICAL via `>${ep}`. Acceptable since the inputs (both `min:0,
    // max:100`) have no cross-field ordering guard and an interactive user is
    // unlikely to set them equal; documented here so the failure mode is
    // discoverable.
    PCT_EXTRAP: {kind: 'conditional', bands: [
      {band: {op: 'le', v: wp}, color: COLOR_OK, unit: '%',
        gloss: 'acceptable extrapolation'},
      {band: {op: 'range', lo: wp, hi: ep}, color: COLOR_WARNING, unit: '%',
        gloss: 'borderline; review'},
      {band: {op: 'gt', v: ep}, color: COLOR_CRITICAL, unit: '%',
        gloss: 'AUC₀–∞ considered unreliable'},
    ]},

    // Two-tier band at the λz adj-R² ACCEPTANCE threshold `R_acc`
    // (lambdaz-adjr2-threshold-spec-v2 §4). A finished fit passes acceptance iff
    // `adjR² ≥ R_acc` — the boundary is inclusive on the "good" side (§4) — so the OK
    // band is `[R_acc, 1]` and the fail band is `< R_acc`. The threshold is applied
    // VERBATIM: no derived `goodR` (§2, removed), so the number painted equals the
    // number the user entered.
    //
    // The fail tint is WARNING (peach), matching the Flagged verdict a violated RSQ
    // criterion rolls up to — NOT CRITICAL/red, which is reserved for the
    // Not-estimated (λz-less) verdict. There is NO separate "below the core gate"
    // tier: sci-comp rejects any window with `adjR² < minRSquared` during selection,
    // so a reported fit never carries an adjR² below the gate (a failed fit writes
    // FLOAT_NULL, which no band paints). `R_acc ≥ minR` is enforced at input, so the
    // OK band can never start below the gate.
    LAMBDAZ_ADJ_R2: {kind: 'conditional', bands: [
      {band: {op: 'ge', v: rAcc, upTo: 1}, color: COLOR_OK,
        gloss: 'meets acceptance threshold'},
      {band: {op: 'lt', v: rAcc}, color: COLOR_WARNING,
        gloss: 'below acceptance threshold — analyst review'},
    ]},

    // NO band below `minP` — the omission is DELIBERATE and load-bearing.
    //
    // A `{op:'lt', v:minP}` band used to sit here, painted COLOR_OK and marked
    // `hidden: true`, justified as "cannot occur — sci-comp rejects sub-minimum
    // fits… no real value lands here". That justification was false. A real value
    // lands there on EVERY failed λz fit: `results-writer.ts` writes
    // `DG.INT_NULL` (= -2147483648) into this column when `lz === null`, and
    // -2147483648 < minP. So a row where the terminal phase could not be fitted at
    // all was tinted mint-green "robust fit", and the band that did it was
    // deliberately withheld from the legend.
    //
    // That is precisely what {@link ColourBand.hidden}'s own contract forbids
    // ("Never use this to hide a band a real value can land in: that is the defect
    // this whole file guards against"), and it regresses the invariant closed by
    // `48998b9` — explain every tinted cell by binding the colour to its legend.
    // Found by user QA, 2026-07-30.
    //
    // A cell that holds no value takes the most conservative option: NO
    // COLOUR. An unpainted cell reports "no fit" honestly; any tint would assign a
    // severity to a value that does not exist. Do NOT re-add a band here to "cover
    // the whole domain" — domain completeness is not a DG requirement we have ever
    // verified, and it is what produced the defect.
    LAMBDAZ_NPOINTS: {kind: 'conditional', bands: [
      {band: {op: 'gt', v: minP}, color: COLOR_OK, gloss: 'robust fit'},
      {band: {op: 'eq', v: minP}, color: COLOR_WARNING,
        gloss: 'minimum allowed; fit may be fragile'},
    ]},

    // Two-tier band on the Rules → %AUMCextrap threshold (FR-200). Deliberately
    // NOT three-tier like PCT_EXTRAP: there is no error-level %AUMCextrap cutoff
    // (no regulatory basis, and the "2-3×" critical ratio is illustrative not
    // validated — research AD-8), so the core emits only a 'warning'-severity
    // AUMC_EXTRAP_HIGH. The band mirrors that single threshold honestly.
    PCT_AUMC_EXTRAP: {kind: 'conditional', bands: [
      {band: {op: 'le', v: wpa}, color: COLOR_OK, unit: '%', gloss: 'acceptable'},
      {band: {op: 'gt', v: wpa}, color: COLOR_WARNING, unit: '%',
        gloss: 'MRT / Vss fragile (house heuristic; no regulatory cutoff)'},
    ]},

    // Three-tier band on the λz window's span in half-lives (RG-C, user-ratified
    // 2026-08-11). Mirrors the LAMBDAZ_ADJ_R2 good / borderline / poor shape, so
    // it needs no new machinery — one more `conditional` entry.
    //
    // WHY THIS COLUMN EARNS A BAND. Span ratio catches a failure mode adj-R²
    // structurally cannot: a 3-point terminal phase over a fraction of one
    // half-life fits a near-perfect straight line, so adj-R² reads "good" while
    // the extrapolated tail — and therefore AUC∞, t½, CL, Vz, MRT and Vss — rests
    // on almost no elimination. Our own reference corpus contains exactly that
    // case: indometh subject 1 pairs a span of 0.685 half-lives with an adj-R² of
    // 0.994 (measured, not illustrative). Banding it is what stops that row
    // reading green.
    //
    // The two thresholds are DELIBERATELY module-local constants and are NOT
    // fields of ResultsColorThresholds — see AD-4 at their declaration.
    SPAN_RATIO: {kind: 'conditional', bands: [
      {band: {op: 'ge', v: SPAN_RATIO_BORDERLINE, upTo: SPAN_RATIO_LEGEND_CEILING},
        color: COLOR_OK, gloss: 'adequate span'},
      {band: {op: 'range', lo: SPAN_RATIO_POOR, hi: SPAN_RATIO_BORDERLINE},
        color: COLOR_WARNING,
        // The conditional form of the finding is load-bearing and must survive
        // into the tooltip: the ≥ 4-half-life bar is reported UNDER HIGH
        // MEASUREMENT VARIABILITY, not unconditionally. An earlier draft restated
        // it as a flat threshold — in the one place a scientist reads it as the
        // tool's own assertion. The citation itself lives in the header tooltip
        // summary (a year in this string would also collide with the
        // legend-vs-rule number guard, which is a second reason to site it there).
        gloss: 'borderline — under high measurement variability, reliable window ' +
          'selection needs a span of about 4 half-lives'},
      {band: {op: 'lt', v: SPAN_RATIO_POOR}, color: COLOR_CRITICAL,
        // NOT "PKNCA's default": PKNCA's `half.life.R` codes no span threshold —
        // its `min.span.ratio` is the same opt-in, diagnostic-only shape as
        // sci-comp's. The 2 is a real and widely-used convention, so the VALUE
        // stays; what changes is that the legend no longer credits a library with
        // a default it does not ship.
        gloss: 'poor — the terminal phase spans under 2 half-lives ' +
          '(below the conventional minimum)'},
    ]},

    MRT: crossColumnFragile,
    VSS: crossColumnFragile,
  };
}

/** Columns hidden in the default grid presentation (spec §7.2 rows 21-22). */
const JSON_COLUMNS = ['WARNINGS_JSON', 'PROVENANCE_JSON'] as const;

/** Display headers. Set via the DG
 *  `friendlyName` column tag — the grid renders this in the header instead
 *  of the underlying column name. Most headers carry no unit (the header
 *  tooltip states it); the per-kg family and its absolute companions DO —
 *  `CL (mL/min/kg)`, `CL abs. (L/h)` — because their basis is what a DMPK
 *  scientist must not misread (nca-cl-bw-normalised-units, user-directed
 *  2026-10-23).
 *
 *  `CL` / `VZ` / `VSS` and the `_ABS` companions are not in this map — their
 *  headers are the writer-set `tags.label` plus the column's own `units` tag
 *  ({@link headerWithUnit}), applied in {@link configureResultsDataFrame}. */
/** Datagrok number formats. DG `'0.000'` is used as a practical
 *  proxy for "3 significant digits" — exact 3-sig-figs would
 *  require a custom cell renderer; in the UC-01 numeric range (1–1000) the
 *  approximation matches.
 *
 *  **An omission here does NOT mean "renders raw".** A Float64 column with no
 *  format falls through to DG's own default, which is TEN decimals: `t_lag = 0`
 *  shipped as `0.0000000000` (user QA 2026-08-14). The comment this one replaces
 *  claimed the opposite — "renders at the actual protocol units … no trailing
 *  zeros" — describing an intent the platform never implemented. `TMAX` still
 *  carries that exemption and still renders at ten decimals: the same defect,
 *  left standing only because this change was scoped to `TLAG` by the user.
 *
 *  Time-valued columns are always in HOURS, never the mapped source unit — the
 *  adapter scales through `TIME_TO_HOUR` in `convertToCanonical` before
 *  `computeNca` sees a point. So a time column's decimal budget is set by the
 *  finest sampling interval expressed in h, not by an unknown unit. */
const NUMBER_FORMATS: Readonly<Record<string, string>> = {
  CMAX: '0.000',
  C0_EXTRAP: '0.000', // mirrors CMAX — same quantity, same units (mg/L)
  AUC_LAST: '0.000',
  AUC_INF: '0.000',
  PCT_EXTRAP: '0.0', // 1 decimal per spec
  HALF_LIFE: '0.000',
  LAMBDA_Z: '0.000',
  CL: '0.000',
  VZ: '0.000',
  CL_ABS: '0.000',
  VZ_ABS: '0.000',
  VSS_ABS: '0.000',
  BW_KG: '0.000',
  AUMC_LAST: '0.000',
  AUMC_INF: '0.000',
  MRT: '0.000',
  VSS: '0.000',
  PCT_AUMC_EXTRAP: '0.0', // 1 decimal, mirrors PCT_EXTRAP
  // Sampling-grid value, not a computed estimate: `t_lag` is the last time with
  // C = 0 before the first measurable concentration, so it can only ever BE a
  // point of the protocol's own time grid. Three decimals is that grid's ceiling
  // (5 min = 0.083 h; everything coarser needs fewer), and `#` drops the rest —
  // a 0.25 / 1 / 24 h study reads "0.25" / "1" / "24", and no lag reads "0", not
  // "0.000". Fixed decimals would be false precision on a quoted time, and the
  // ten-decimal DG default was the reported defect. User decision 2026-08-14.
  TLAG: '0.###',
  LAMBDAZ_R2: '0.000', // 3 decimals per spec
  LAMBDAZ_ADJ_R2: '0.000', // 3 decimals per spec
  // 2 decimals: the band edges are whole half-lives, so a third decimal would be
  // false precision on a diagnostic read against "2" and "4".
  SPAN_RATIO: '0.00',
  LAMBDAZ_TFIRST: '0.00',
  LAMBDAZ_TLAST: '0.00',
  LAMBDAZ_INTERCEPT: '0.000',
  // TMAX: no format (raw protocol value)
  // LAMBDAZ_NPOINTS: integer, no format needed
};

const DISPLAY_HEADERS: Readonly<Record<string, string>> = {
  SUBJECT: 'Subject',
  ROUTE: 'Route',
  CL_BASIS: 'CL basis',
  DOSE_NORM: 'Dose norm.',
  // 'Computation', not 'Status' (nca-grid-status-spec §1): the column reflects the
  // technical fact of computation (ok/partial/failed), renamed to separate it from
  // the ACCEPTANCE verdict below. Display-only — the underlying column NAME stays
  // STATUS (writer schema + the forceRepaint nudge target both key off it).
  STATUS: 'Computation',
  // λz acceptance verdict + flags (nca-grid-status-spec §2/§3). The shared 'Flag:'
  // prefix stands in for a native column group: an embedded DG.Viewer.grid IGNORES
  // `.columnGroups`, so the group heading rides in the friendly names instead.
  ACCEPTANCE: 'Acceptance',
  [ACCEPT_REASONS_COL]: 'Violated',
  REVIEW: 'Review',
  PIPELINE: 'Pipeline',
  SPAN_FLAG: 'Flag: span',
  RSQ_FLAG: 'Flag: R²',
  NEXT_FLAG: 'Flag: N',
  EXTRAP_FLAG: 'Flag: %extrap',
  CMAX: 'C_max',
  // `C_0`, not the raw schema name and not WinNonlin's bare `C0` — follows the
  // underscore convention every sibling here uses (C_max / t_max / t_lag).
  // The "estimate, not observation / IV-bolus-only" caveat lives in the header
  // tooltip, which has room to state it honestly (user decision 2026-07-29).
  C0_EXTRAP: 'C_0',
  TMAX: 't_max',
  AUC_LAST: 'AUC_0–t',
  AUC_INF: 'AUC_0–∞',
  PCT_EXTRAP: '%AUCextrap',
  HALF_LIFE: 't_½',
  LAMBDA_Z: 'λz',
  AUMC_LAST: 'AUMC_0–t',
  AUMC_INF: 'AUMC_0–∞',
  MRT: 'MRT',
  TLAG: 't_lag',
  PCT_AUMC_EXTRAP: '%AUMCextrap',
  BW_KG: `BW (${BW_UNIT})`,
  LAMBDAZ_R2: 'R² (λz)',
  LAMBDAZ_ADJ_R2: 'adj. R² (λz)',
  LAMBDAZ_NPOINTS: 'N points (λz)',
  SPAN_RATIO: 'Span ratio (λz)',
  LAMBDAZ_TFIRST: 't_first (λz)',
  LAMBDAZ_TLAST: 't_last (λz)',
  LAMBDAZ_INTERCEPT: 'ln(C₀) intercept',
  AUC_METHOD: 'AUC method',
  COMPENSATED: 'Compensated',
};

/** Presentation-only column order.
 *  Differs from the underlying writer schema only on the HALF_LIFE / LAMBDA_Z
 *  pair (writer ships LAMBDA_Z before HALF_LIFE; presentation shows HALF_LIFE
 *  first as the human-meaningful disposition parameter). Underlying DataFrame
 *  schema stays untouched — applies to grid layout only. */
const PRESENTATION_COLUMN_ORDER: ReadonlyArray<string> = [
  // ACCEPTANCE rides right after STATUS: the two verdicts (did it run / can it be
  // trusted) read together at the front. It is compact-visible via the derived
  // banded floor; STATUS stays for the same reason.
  // `Violated` (the reason pills) rides right after the verdict it explains.
  'SUBJECT', 'ROUTE', 'CL_BASIS', 'DOSE_NORM', 'STATUS', 'ACCEPTANCE', ACCEPT_REASONS_COL, REVIEW_COL,
  PIPELINE_COL,
  'CMAX', 'C0_EXTRAP', 'TMAX', 'AUC_LAST', 'AUC_INF', 'PCT_EXTRAP',
  'HALF_LIFE', 'LAMBDA_Z', 'CL', 'VZ',
  // FR-200 moment / lag params after VZ (AD-13 display order).
  'AUMC_LAST', 'AUMC_INF', 'MRT', 'VSS', 'TLAG', 'PCT_AUMC_EXTRAP',
  // Absolute companions + the weight (nca-cl-bw-normalised-units D-a/c): neither banded nor in
  // COMPACT_CORE, so Full only (non-parameters live in Full).
  'CL_ABS', 'VZ_ABS', 'VSS_ABS', 'BW_KG',
  'LAMBDAZ_R2', 'LAMBDAZ_ADJ_R2', 'LAMBDAZ_NPOINTS', 'SPAN_RATIO',
  // The four λz acceptance flags sit CONTIGUOUS, right after the λz diagnostics
  // they summarize, so the shared 'Flag:' header prefix reads as one group. They
  // are Full-only (non-banded → they fall through the compact partition), for
  // filtering / counting rather than at-a-glance reading (ACCEPTANCE carries that).
  'SPAN_FLAG', 'RSQ_FLAG', 'NEXT_FLAG', 'EXTRAP_FLAG',
  'LAMBDAZ_TFIRST', 'LAMBDAZ_TLAST', 'LAMBDAZ_INTERCEPT',
  'AUC_METHOD', 'COMPENSATED',
  'WARNINGS_JSON', 'PROVENANCE_JSON',
];

/** Columns whose grid header comes from the writer-set `label` tag rather than
 *  {@link DISPLAY_HEADERS}, because it is route-dependent (`CL` vs `CL/F`, `Vz`
 *  vs `Vz/F`) and cannot be a static string; the header appends the column's `units`
 *  tag (`CL (mL/min/kg)`, `CL abs. (L/h)` — user-directed 2026-10-23). ONE list,
 *  consumed by both the header-hoisting loop in {@link configureResultsDataFrame} and
 *  {@link HEADER_EXEMPT} — two hand-kept copies of these names would be the very
 *  drift this module's guard exists to prevent. */
const ROUTE_LABELED_COLUMNS: ReadonlyArray<string> = ['CL', 'VZ', 'VSS', 'CL_ABS', 'VZ_ABS', 'VSS_ABS'];

/** A grid header with the column's unit: `CL (mL/min/kg)`. The unit is the column's own `units` tag
 *  (one source — `CANONICAL_UNIT_TAGS`), never a restated literal; no tag → the bare label. */
export function headerWithUnit(label: string, col: DG.Column): string {
  const unit = col.getTag('units');
  return unit !== null && unit !== undefined && unit !== '' ? `${label} (${unit})` : label;
}

/** Columns intentionally absent from {@link DISPLAY_HEADERS}: the route-labelled
 *  trio above, plus the two JSON audit columns (hidden by default). Every OTHER
 *  presentation column must carry a static friendly name. */
const HEADER_EXEMPT: ReadonlyArray<string> = [...ROUTE_LABELED_COLUMNS, ...JSON_COLUMNS];

/** Columns intentionally absent from {@link NUMBER_FORMATS}. Non-numeric columns
 *  (string / bool / the JSON tail) plus two numerics: `LAMBDAZ_NPOINTS` is an
 *  integer count, and `TMAX` is a survivor of the retired "renders raw" belief
 *  documented on {@link NUMBER_FORMATS} — it is exempt, and it does render at
 *  DG's ten-decimal default. `TLAG` left this list on 2026-08-14 and took a
 *  `'0.###'` format; `TMAX` is the identical case awaiting the same call.
 *  Every OTHER presentation column must carry an explicit format, or it renders
 *  at DG's default precision beside a 3-decimal sibling. */
const FORMAT_EXEMPT: ReadonlyArray<string> = [
  'SUBJECT', 'ROUTE', 'CL_BASIS', 'DOSE_NORM', 'STATUS',
  // λz acceptance verdict + flags — all STRING, no number format.
  ACCEPTANCE_COL, ACCEPT_REASONS_COL, ...FLAG_COLS, REVIEW_COL, PIPELINE_COL,
  'TMAX', 'LAMBDAZ_NPOINTS', 'AUC_METHOD', 'COMPENSATED',
  ...JSON_COLUMNS,
];

/** Per-subject grid column widths, in px. An OVERRIDE map, not a layout: a
 *  column absent here keeps DG's content-derived default, which is right for the
 *  3-decimal parameter block. Applied once per `configureResultsGrid` (i.e. per
 *  Compute) rather than per profile toggle, so a manual resize survives until the
 *  next recompute — the same semantics {@link applyFixedGridSpec} gives the
 *  sparse / dose-proportionality tabs.
 *
 *  `TLAG` is here because DG sized it against the ten-decimal `0.0000000000`
 *  render it no longer produces; its real content is `0` / `0.083` / `0.25` under
 *  a 5-character header. 56 px mirrors the sparse tab's `Tmax` — same quantity
 *  class, same header length. */
const GRID_WIDTHS: Readonly<Record<string, number>> = {
  TLAG: 56,
  // Three pills (`span < 2` `N ≤ 3` `%extrap > 20%`) — the most measured on the fixtures; narrower,
  // the Tags renderer cuts a pill to `N…` and hides the reason (QA 2026-10-19).
  [ACCEPT_REASONS_COL]: 170,
};

/** Presentation columns that carry no friendly header / no header tooltip / no
 *  number format / no colour legend. All four arrays are empty in a healthy build. */
interface MissingPresentationEntries {
  readonly noHeader: ReadonlyArray<string>;
  readonly noTooltip: ReadonlyArray<string>;
  readonly noFormat: ReadonlyArray<string>;
  /** Columns this module COLOURS whose header tooltip carries no legend rows —
   *  the scientist sees a tint and hovering does not say what it means. */
  readonly noLegend: ReadonlyArray<string>;
}

/**
 * Drift guard for the per-subject results grid's presentation maps.
 *
 * `C0_EXTRAP` reached scientists as a raw `C0_EXTRAP` header carrying DG's
 * generic min/avg/max/std stats tooltip: it was added to the writer schema and
 * to {@link PRESENTATION_COLUMN_ORDER}, but to none of the presentation maps.
 * The writer already drift-guards its own side (`FIXED_RESULT_COLUMNS` vs the
 * built `cols` array); this is the missing presentation-side counterpart — a
 * pure function needing no DG Grid, asserted empty by
 * `results-view-presentation-tests`. A column added to the order without a
 * header, a tooltip and a number format now fails the suite instead of
 * shipping unexplained.
 *
 * `noLegend` is the FOURTH list, added 2026-07-30 after MRT and Vss shipped
 * tinted-but-unexplained. It closes a different hole from the first three: those
 * catch a column with NO tooltip, this catches a column whose tooltip exists,
 * reads well, and simply never mentions the colour on the cell. The three
 * original lists were all green while the defect was live — a tooltip's
 * PRESENCE was checked, its CONTENT was not.
 *
 * Each list has its own explicit exemption set ({@link HEADER_EXEMPT},
 * {@link FORMAT_EXEMPT}) — the point of the guard is to separate a DELIBERATE
 * omission from an ACCIDENTAL one, which a silent `undefined` fallback cannot
 * do. `noLegend` needs no exemption list: it is driven by
 * {@link buildColourRules}, so a column is checked exactly when it is coloured.
 *
 * @returns Presentation columns missing a header / tooltip / format / legend.
 */
export function missingPresentationEntries(): MissingPresentationEntries {
  const tooltips = buildHeaderTooltips(DEFAULT_RESULTS_THRESHOLDS);
  const coloured = buildColourRules(DEFAULT_RESULTS_THRESHOLDS);
  const headerExempt = new Set<string>(HEADER_EXEMPT);
  const formatExempt = new Set<string>(FORMAT_EXEMPT);
  return {
    noHeader: PRESENTATION_COLUMN_ORDER.filter(
      (n) => !headerExempt.has(n) && !(n in DISPLAY_HEADERS)),
    noTooltip: PRESENTATION_COLUMN_ORDER.filter((n) => !(n in tooltips)),
    noFormat: PRESENTATION_COLUMN_ORDER.filter(
      (n) => !formatExempt.has(n) && !(n in NUMBER_FORMATS)),
    noLegend: Object.keys(coloured).filter(
      (n) => (tooltips[n]?.extra?.length ?? 0) === 0),
  };
}

// ──────────────────────────────────────────────────────────────────────────
// Compact / full column profile (F1) + inapplicable-column suppression (F5)
// ──────────────────────────────────────────────────────────────────────────

/** `view.temp` key holding the results grid's active profile (`true` = Full).
 *  Deliberately on the VIEW, not the frame: `swapOrDockResults` replaces the
 *  DataFrame on every Compute, so the study grid's `df.temp` mechanism (where
 *  the frame persists) would silently reset the scientist's choice each run. */
export const RESULTS_EXTENDED_TAG = 'ncaResultsExtended';

/** The VR-20 scan set plus identity — what a Discovery DMPK lead ranks on.
 *  Deliberately NOT the whole story: the mandatory-visible floor below adds
 *  every colour-banded column, and the union is the compact profile. */
const COMPACT_CORE: ReadonlyArray<string> = [
  'SUBJECT', 'CMAX', 'TMAX', 'AUC_LAST', 'AUC_INF', 'HALF_LIFE', 'CL', 'VZ',
];

/**
 * The columns the compact profile may NEVER hide, **derived** from the one
 * declaration that paints them.
 *
 * This is the load-bearing decision of the feature (AD-2). A colour-banded
 * column that the default view hides is a tinted cue the scientist cannot see —
 * the honest-uncertainty argument (thesis #4) rendered invisible by a column
 * profile. The obvious implementation is to write the banded set down here; the
 * record says not to. **Three independent readers have restated that set and
 * produced three different wrong answers**: the GAP-23 proposal said 6, the
 * architecture doc said 4, and the truth at the time was 7 (it is 8 now that
 * `SPAN_RATIO` is banded). Deriving it makes the whole defect class
 * structurally unreachable, and a ninth band added later stays visible with no
 * edit here.
 *
 * The band-key SET is threshold-independent (thresholds move the edges, not
 * which columns are banded), so passing the defaults purely to enumerate keys is
 * sound — the same trick {@link missingPresentationEntries} already uses so that
 * "a column is checked exactly when it is coloured."
 */
function bandedColumns(): ReadonlySet<string> {
  return new Set(Object.keys(buildColourRules(DEFAULT_RESULTS_THRESHOLDS)));
}

/**
 * Schema columns the compact profile never shows because they are not
 * PARAMETERS — they describe the study, not its result (user-directed
 * 2026-08-13).
 *
 * `ROUTE` is study design: which arm a subject was on. Useful, but it is a PROXY
 * for the thing that changes how a number must be read, and that thing has its
 * own column — see {@link COMPACT_WHEN_VARYING}.
 *
 * `STATUS` is deliberately NOT in this list even though it is not a parameter
 * either: it is colour-banded, and the AD-2 floor keeps a banded column visible
 * because a tint the scientist cannot see is thesis-#4 uncertainty made invisible.
 * It is, however, hidden by the SEPARATE value-based {@link statusSuppressed} rule
 * when every row is `ok` — legal under the same "band can't fire" exception as VSS,
 * NOT via this list. `DOSE_NORM`, `AUC_METHOD` and `COMPENSATED` need no entry —
 * they are neither in {@link COMPACT_CORE} nor banded, so they were already Full-only.
 */
const COMPACT_EXCLUDED_NON_PARAMETERS: ReadonlyArray<string> = ['ROUTE'];

/**
 * Non-parameter columns the compact profile shows **only when they carry more
 * than one distinct value**.
 *
 * **Why `CL_BASIS` earns the exception and `ROUTE` does not.** `CL` and `Vz`
 * carry a NEUTRAL header whenever a frame mixes routes — a column tag cannot
 * express a per-row property, so the writer deliberately refuses to label the
 * whole column `CL/F`. `CL_BASIS` is that per-row label. With it hidden on a
 * mixed-route frame there was nothing on screen separating true CL (intravascular,
 * F = 1) from apparent CL/F (extravascular, includes bioavailability) — values
 * that differ by ~1/F. A DMPK scientist ranking compounds by CL could read one as
 * the other with no cue at all: thesis #4's exact failure mode, on the parameter
 * the primary persona ranks on.
 *
 * The `> 1 distinct` gate is what keeps this from undoing the parameters-only
 * rule: on a single-route study `CL_BASIS` repeats one value and stays hidden, so
 * the column returns **only** on the frames where the ambiguity is real. That is
 * also precisely the case RG-J existed for — subsuming it wholesale removed the
 * cue exactly where it was load-bearing.
 */
const COMPACT_WHEN_VARYING: ReadonlyArray<string> = ['CL_BASIS'];

/** Distinct NON-NULL values in a column. Nulls are excluded deliberately (RG-J):
 *  counting them would make a single-route frame with one failed row look like
 *  two distinct bases and force the column on for the wrong reason. */
function distinctNonNull(col: DG.Column): number {
  const seen = new Set<string>();
  for (let i = 0; i < col.length; i++) {
    if (col.isNone(i)) continue;
    seen.add(col.getString(i));
  }
  return seen.size;
}

/** What the frame's `ROUTE` column says about which route-gated parameters can
 *  hold a value at all. `known` is false when there is no route information, in
 *  which case NOTHING is suppressed — an unknown route is not evidence of
 *  inapplicability. */
function routeCoverage(df: DG.DataFrame): {
  known: boolean; anyIvBolus: boolean; anyIv: boolean; anyNonIv: boolean;
} {
  const col = df.col('ROUTE');
  let anyIvBolus = false; let anyIv = false; let anyNonIv = false; let seen = 0;
  if (col !== null) {
    for (let i = 0; i < col.length; i++) {
      if (col.isNone(i)) continue;
      seen++;
      const v = col.getString(i);
      if (v === 'IV-bolus') {anyIvBolus = true; anyIv = true;} else if (v === 'IV-infusion') anyIv = true;
      else anyNonIv = true;
    }
  }
  return {known: seen > 0, anyIvBolus, anyIv, anyNonIv};
}

/**
 * Columns the compact profile suppresses because this frame's ROUTES make them
 * structurally incapable of holding a value — never because they *happen* to be
 * empty (F5 / AC-5.2).
 *
 * **The discriminator is route, and that is the whole point.** An all-null test
 * cannot tell "inapplicable here" from "the computation failed for every
 * subject", and the second case is a finding the scientist must see. So a study
 * whose λz fits all failed keeps its entire λz block on screen holding explicit
 * nulls — that IS the result.
 *
 * The three parameters have **three different** route contracts and must not
 * share one predicate (they nearly did):
 *
 * | Column | Null when | Therefore hide iff |
 * |---|---|---|
 * | `VSS` | route is not IV — bolus **or** infusion both count | no row is IV |
 * | `C0_EXTRAP` | route is not IV-**bolus** — an IV-infusion row is null too | no row is IV-bolus |
 * | `TLAG` | route **is** IV | every row is IV |
 *
 * Collapsing `VSS` and `C0_EXTRAP` into one "no IV rows" rule leaves an
 * **infusion-only** study showing a permanently-null `C0_EXTRAP` — reproducing
 * the exact dead-column defect this function exists to remove.
 *
 * `C0_EXTRAP` additionally has null semantics unlike every other parameter: it
 * is NOT null on a *failed* IV-bolus row, because sci-comp's default
 * back-extrapolation chain ends in `set0`, which returns 0 unconditionally. So
 * it can never be all-null on a frame containing any IV-bolus row, and this
 * route rule is the only thing that can hide it. Do not later "simplify" it into
 * an emptiness test — an emptiness test would never fire.
 */
function routeSuppressed(df: DG.DataFrame): ReadonlySet<string> {
  const out = new Set<string>();
  const r = routeCoverage(df);
  if (!r.known) return out;
  // VSS is BANDED, so hiding it needs the AC-5.5 justification, and it has one
  // structurally rather than by luck: on a frame with no IV rows every VSS cell
  // is null, and the cross-column tint's own first guard is `own === null`
  // (`GridCell.value` is the one read API that honours null). The tint therefore
  // cannot paint on such a frame, so hiding the column removes no cue. The
  // invariant is "never hide a banded column whose band CAN FIRE" — not "never
  // hide a banded column".
  if (!r.anyIv) out.add('VSS');
  // Its absolute companion has the same IV-only contract (architect S5); it is not banded.
  if (!r.anyIv) out.add('VSS_ABS');
  if (!r.anyIvBolus) out.add('C0_EXTRAP');
  if (!r.anyNonIv) out.add('TLAG');
  return out;
}

/**
 * `{STATUS}` when EVERY subject computed cleanly (status === 'ok'), else empty.
 *
 * STATUS (Computation) is banded, so hiding it needs the AC-5.5 justification —
 * and it has it: the column's only NON-happy bands (`partial` → peach, `failed` →
 * red) CANNOT fire when every row is `ok`, and the `ok` band is the happy path,
 * not a cue the scientist must act on. So an all-ok Computation column carries
 * nothing, and hiding it in compact declutters the common case. The moment ANY row
 * is partial/failed the band can fire, so the column returns — the same "never hide
 * a banded column whose band CAN fire" rule that governs VSS. **Compact-only** (Full
 * means full).
 *
 * NOT the forbidden emptiness test: STATUS is never null, so "all ok" is a definite
 * state, not the ambiguous "all-null" {@link routeSuppressed} must never key off.
 * An all-FAILED study keeps STATUS visible (every row `failed` → the band fires) —
 * exactly the case the scientist most needs to see.
 */
function statusSuppressed(df: DG.DataFrame): ReadonlySet<string> {
  const col = df.col('STATUS');
  if (col === null) return new Set();
  for (let i = 0; i < col.length; i++) {
    if (col.isNone(i)) continue;
    if (col.getString(i) !== 'ok') return new Set(); // a non-ok row → the band can fire
  }
  return new Set(['STATUS']);
}

/** Pill fill of the `Violated` column, one per λz criterion, so a reason reads at a glance and down the
 *  column (user-directed 2026-10-21, `acceptance-reason-pill-colors`, superseding the single warm grey
 *  `#F3F1EE`). Light Okabe-Ito-derived tints: colour-blind-safe, and none of them is a verdict hue
 *  (`status-palette.ts`). Fill only — PowerGrid's `Tags` renderer draws no border and writes the
 *  contrast colour (dark) on the fill, so the mockup's text / border colours are not reachable natively.
 *  Keyed by CODE, not label: the label carries the live threshold. */
const REASON_PILL_FILLS: Readonly<Record<AcceptanceFlagCode, number>> = {
  SPAN_LT: 0xFFBFE4F7, // blue
  ADJR2_LT: 0xFFF5DB98, // yellow
  NPTS_LE: 0xFFEDC6DD, // pink
  EXTRAP_GT: 0xFFF3C6A6, // orange
};

/** Pill text size — one step below the grid's 13 px, so the reasons read as subordinate to the
 *  verdict beside them. The renderer draws in the cell font, set per cell in `onCellPrepare`. */
const REASON_PILL_FONT_PX = 12;

/** Paint each `Violated` pill in its criterion's {@link REASON_PILL_FILLS}. `Tags` keys a colour by the
 *  pill TEXT, so each label the thresholds produce gets its code's fill (the threshold moves the text,
 *  hence per call). */
export function applyReasonPillColours(grid: DG.Grid, t: ResultsColorThresholds): void {
  const gc = grid.col(ACCEPT_REASONS_COL);
  if (gc === null) return;
  const colours: Record<string, number> = {};
  for (const code of ACCEPTANCE_FLAG_CODES) colours[criterionShortLabel(code, t)] = REASON_PILL_FILLS[code];
  gc.categoryColors = colours;
}

/** Set a `Violated` cell's font to {@link REASON_PILL_FONT_PX}, PRESERVING its weight/family (the
 *  {@link boldGridCell} idiom). Called from each grid's single `onCellPrepare`. */
function reasonPillFont(gc: DG.GridCell): void {
  const base = gc.style.font;
  gc.style.font = (base !== null && base !== undefined && /\d+px/.test(base)) ?
    base.replace(/\d+px/, `${REASON_PILL_FONT_PX}px`) : `${REASON_PILL_FONT_PX}px Roboto`;
}

/** True when some row of the `Violated` column holds a reason (≥ 1 Flagged profile). */
function anyReason(df: DG.DataFrame): boolean {
  const col = df.col(ACCEPT_REASONS_COL);
  if (col === null) return false;
  for (let i = 0; i < col.length; i++)
    if (!col.isNone(i) && col.getString(i) !== '') return true;
  return false;
}

/** True when a compact grid shows / hides `Violated` against its current content — after a λz
 *  override flipped the first (or last) profile to / from Flagged. The caller then re-applies the
 *  profile; checked first so an ordinary override keeps the scientist's column order. */
export function reasonsVisibilityStale(grid: DG.Grid, extended: boolean): boolean {
  const df = grid.dataFrame;
  const gc = grid.col(ACCEPT_REASONS_COL);
  if (extended || df == null || gc === null || df.getTag(REASON_PILLS_TAG) === 'false') return false;
  return gc.visible !== anyReason(df);
}

/**
 * The ordered list of VISIBLE column names for the results grid.
 *
 * Pure and grid-free so the partition is unit-testable without a live render —
 * the `applyStudyGridHeaders` lesson, where header logic that needed a real grid
 * to exercise shipped "undefined undefined" past every automated test.
 *
 * Order: `SUBJECT`, then (Full only) the carried metadata columns, then the
 * profile's parameters in presentation order.
 *
 * **Compact is SUBJECT + PARAMETERS** (user-directed 2026-08-13). Nothing here
 * removes data: this function decides grid-side VISIBILITY only, the DataFrame
 * always keeps every column, and the scientist can re-enable any of them
 * individually through Datagrok's own column chooser (right-click a header) —
 * which is what makes "hidden by default" an honest default rather than a
 * suppression.
 *
 * @param extended `false` = compact (identity + parameters: the scan set plus
 *   every fragility cue), `true` = every presentation column. Full means full:
 *   the route-based suppression above is a DEFAULT-view courtesy, not a filter
 *   the scientist cannot switch off.
 */
export function resultsGridVisibleOrder(df: DG.DataFrame, extended: boolean): string[] {
  const present = new Set(df.columns.names());
  const order: string[] = [];
  const push = (name: string): void => {if (present.has(name)) order.push(name);};

  push('SUBJECT');

  // Carried metadata (study id, species, strain, sex, dose, matrix, anticoagulant,
  // sampling type, SMILES, …) — appended by `appendGroupingColumns`, so they are
  // exactly the present columns the fixed schema does not know about.
  //
  // FULL PROFILE ONLY (user-directed 2026-08-13, replacing "always visible in
  // both"). On a real rat cassette this block was NINE columns before the first
  // parameter — and seven of them (`StudyID`, `Species`, `Strain`, `Matrix`,
  // `Anticoagulant`, `SamplingType`, `Smiles`) carried ONE value for the whole
  // study, so they pushed C_max off the first screen to say nothing. The prior
  // rule read them as "covariates the scientist grouped by"; in practice
  // `appendGroupingColumns` carries every constant-per-subject column, grouped-on
  // or not, so most of them were never chosen by anyone.
  //
  // Unknown names are still KEPT rather than dropped — they simply ride the Full
  // profile now, and the column chooser reaches any one of them individually.
  const schema = new Set<string>(PRESENTATION_COLUMN_ORDER);
  if (extended) {
    for (const name of df.columns.names()) {
      // `_PROBLEM_RANK` is a hidden sort helper, not a result column — it must not
      // ride the Full profile the way other unknown (carried-metadata) columns do.
      // Same for any `_lz*` service column (e.g. `_lz_overridden`, F7): it feeds the
      // ✏️ render hook, it is not a scientist-facing column.
      if (!schema.has(name) && name !== 'SUBJECT' && name !== PROBLEM_RANK_COL &&
          !isLzServiceColumn(name)) order.push(name);
    }
  }

  const json = new Set<string>(JSON_COLUMNS);
  const banded = bandedColumns();
  // Compact-only suppression: route-structural (VSS/C0/TLAG) + STATUS when every
  // row computed cleanly. Full shows everything.
  const suppressed = extended ? new Set<string>() :
    new Set<string>([...routeSuppressed(df), ...statusSuppressed(df)]);
  const nonParameter = new Set<string>(COMPACT_EXCLUDED_NON_PARAMETERS);
  const whenVarying = new Set<string>(COMPACT_WHEN_VARYING);

  const inCompact = (name: string): boolean => {
    // Route-structural suppression wins over everything below: a column that
    // cannot hold a value on this frame is noise in a scan view.
    if (suppressed.has(name)) return false;
    // Study description, not result — see COMPACT_EXCLUDED_NON_PARAMETERS.
    if (nonParameter.has(name)) return false;
    // A non-parameter that changes how a PARAMETER must be read, shown only where
    // it actually changes — see COMPACT_WHEN_VARYING.
    if (whenVarying.has(name)) {
      const col = df.col(name);
      return col !== null && distinctNonNull(col) > 1;
    }
    // `C0_EXTRAP` / `TLAG` are conditional-MAIN: route-gated, but genuinely part
    // of the scan set wherever the route can produce them. Reaching here means
    // `routeSuppressed` did NOT hide them, so this frame's routes can populate
    // them (or carry no route information at all, in which case nothing is
    // suppressed — an unknown route is not evidence of inapplicability).
    if (name === 'C0_EXTRAP' || name === 'TLAG') return true;
    // Stage-scoped (D15): compact-visible on Review λz only (frame tag set by the view).
    if (name === REVIEW_COL) return df.getTag(REVIEW_COMPACT_TAG) === 'true';
    // The pipeline run's own reasons (nca-pipeline-save-rework F10, approved design 2026-10-15): the
    // column exists only when it carries a reason (`writePipelineColumns`), so present ⇒ compact.
    if (name === PIPELINE_COL) return true;
    // The λz reason pills (acceptance-reason-pills): shown only when some profile is Flagged — on
    // a frame with no violated criterion every cell is empty and the column says nothing.
    if (name === ACCEPT_REASONS_COL) return anyReason(df);
    if (COMPACT_CORE.includes(name)) return true;
    return banded.has(name); // AD-2: derived floor, never a restated list
  };

  // The Settings switch for the reason pills overrides both profiles.
  const pillsOff = df.getTag(REASON_PILLS_TAG) === 'false';
  for (const name of PRESENTATION_COLUMN_ORDER) {
    if (name === 'SUBJECT' || json.has(name)) continue;
    if (name === ACCEPT_REASONS_COL && pillsOff) continue;
    if (!extended && !inCompact(name)) continue;
    if (extended && suppressed.has(name)) continue;
    push(name);
  }
  return order;
}

/**
 * Apply the compact/full profile to the grid: order the visible columns per
 * {@link resultsGridVisibleOrder}, then hide everything else.
 *
 * The FRAME is untouched — only grid-side order and visibility — so sorting a
 * hidden-then-revealed column still works, and a frame-level export still
 * carries every column (measured: a grid-aware `toCsv` honours visibility, a
 * frame-level one does not).
 *
 * Idempotent, and safe to re-run after `grid.dataFrame = newDf` or after
 * `appendGroupingColumns` adds columns to the frame in place.
 */
export function applyResultsColumnProfile(grid: DG.Grid, extended: boolean): void {
  const df = grid.dataFrame;
  if (df === null || df === undefined) return;
  const visible = resultsGridVisibleOrder(df, extended);
  const visibleSet = new Set(visible);
  // setOrder gets EVERY column name (visible first, then the rest) so DG keeps
  // the frame complete while the profile drives the left-to-right order.
  const hidden = df.columns.names().filter((n) => !visibleSet.has(n));
  grid.columns.setOrder([...visible, ...hidden]);
  for (const name of df.columns.names()) {
    const gc = grid.col(name);
    if (gc !== null) gc.visible = visibleSet.has(name);
  }
  fitResultsRowHeight(grid, visible);
  forceRepaintAfterVisibilityChange(grid, visible);
}

/** The platform grid's default row height (px), measured on the results grid with no structure column. */
const DEFAULT_RESULTS_ROW_HEIGHT = 28;
/** `grid.temp` key: the row height the platform chose for a structure column, kept to restore in Full. */
const STRUCT_ROW_HEIGHT_KEY = 'ncaStructRowHeight';

/**
 * Size the results rows for what is VISIBLE. A carried `Smiles` column is re-detected as `Molecule`
 * by the platform (Chem's detector; `buildPerSubjectColumn` strips the semType, the detector sets
 * it again) and the grid then raises EVERY row to the molecule height (measured 109 px) even while
 * the column is hidden in Compact — tall, mostly empty rows with the pills stuck at the top
 * (user-reported 2026-11-03). Compact (no visible structure) → the default height; Full with a
 * visible structure → the platform's molecule height, remembered when it raised it.
 * Any height above the default is remembered as that height, so rows an analyst enlarged by hand
 * return to the default on the next profile apply (Confirm, Compact/Full) — accepted: the grid
 * keeps no user row height anywhere else either.
 */
function fitResultsRowHeight(grid: DG.Grid, visible: readonly string[]): void {
  const df = grid.dataFrame;
  const current = grid.props.rowHeight;
  if (current > DEFAULT_RESULTS_ROW_HEIGHT) grid.temp[STRUCT_ROW_HEIGHT_KEY] = current;
  const structureShown = visible.some((n) => df.col(n)?.semType === DGNS.SEMTYPE.MOLECULE);
  if (!structureShown) {
    if (current !== DEFAULT_RESULTS_ROW_HEIGHT) grid.props.rowHeight = DEFAULT_RESULTS_ROW_HEIGHT;
    return;
  }
  const remembered = grid.temp[STRUCT_ROW_HEIGHT_KEY] as number | null | undefined;
  if (typeof remembered === 'number' && current < remembered) grid.props.rowHeight = remembered;
}

/** `grid.temp` key: the subscription to the CURRENT frame's semantic-type detection. */
const SEMTYPE_SUB_KEY = 'ncaSemTypeSub';

/**
 * Re-fit the row height when the platform detects a semantic type on the results frame after the
 * grid was configured (the detection that sets `Smiles` to `Molecule` is not ours and may come
 * later). One subscription per frame: a re-Compute swaps the frame, so the previous one is dropped.
 */
function watchSemanticTypes(grid: DG.Grid): void {
  (grid.temp[SEMTYPE_SUB_KEY] as {unsubscribe(): void} | null | undefined)?.unsubscribe();
  const df = grid.dataFrame;
  grid.temp[SEMTYPE_SUB_KEY] = df.onSemanticTypeDetected.subscribe(() => {
    // Deferred a tick: the grid raises its row height in response to the same detection. The
    // visible set is read from the grid, so a Compact/Full toggle since subscribing is honoured.
    setTimeout(() => {
      if (grid.dataFrame !== df) return;
      fitResultsRowHeight(grid, df.columns.names().filter((n) => grid.col(n)?.visible === true));
    }, 0);
  });
}

/**
 * Make the grid actually PAINT the visibility change just written above.
 *
 * **This exists because the state was right and the canvas was not.** A freshly
 * docked results grid rendered only `Subject` and `Status`; touching a column
 * border — any resize — made the remaining fifteen appear at once, unchanged.
 * That rules out this package's partition logic, and the measurements agree: the
 * production sequence reproduced end-to-end on the reported all-IV-bolus indometh
 * frame reports 17 columns, `visible === true` on every one, widths 69/63/64/…/111
 * summing to 1264 px in a 1145 px host, and no unresolved `GridColumn`
 * (`results-grid-platform-tests`, "the compact profile RENDERS its columns"). The
 * defect is a stale canvas, so the fix has to be a repaint, not a partition change.
 *
 * Two steps, and the second is not redundant:
 *   1. `Grid.invalidate()` — the platform's own documented "causes the grid to
 *      repaint". Tried first because it is the intended call and has no side
 *      effects.
 *   2. A real one-pixel width change, applied and then reverted on the next
 *      macrotask. The reported symptom IS that a visibility-only change failed to
 *      reach the canvas while a geometry change succeeded, so this reproduces the
 *      one input known to work rather than trusting that (1) covers it.
 *
 * The revert is deferred rather than inline because setting a width to `w + 1`
 * and back to `w` within one tick can coalesce into no change at all — which
 * would make this whole function a no-op.
 *
 * Deferring, though, is exactly what makes the naive version DRIFT, and the guard
 * in `results-grid-platform-tests` caught it at 64 px → 66 px: `runCompute`
 * applies the profile several times in one tick (configure, again after
 * `appendGroupingColumns`, again per toggle), so call 2 reads the width call 1 had
 * already bumped, and each deferred revert restores its own inflated capture.
 * Hence {@link pendingWidthNudge}: the FIRST nudge on a grid records the true
 * base, later nudges in the same tick cancel the pending revert and reuse it, and
 * the column lands back exactly where it started however many times this runs.
 *
 * If a later platform build makes (1) sufficient on its own, DELETE (2) — do not
 * leave both standing on the assumption that more nudging is safer.
 */
function forceRepaintAfterVisibilityChange(
  grid: DG.Grid, visible: ReadonlyArray<string>,
): void {
  try {
    grid.invalidate();
  } catch {
    // Older platform build without the call — step 2 still applies.
  }
  // `STATUS` by preference when it is visible, else whatever leads the visible
  // order (STATUS is now hidden in compact on an all-`ok` frame — statusSuppressed —
  // so the fallback is a real path, not a degenerate one). Never a hidden column:
  // nudging one changes nothing on screen.
  const target = visible.includes('STATUS') ? 'STATUS' : (visible[0] ?? null);
  if (target === null) return;
  const gc = grid.col(target);
  if (gc === null) return;

  const pending = pendingWidthNudge.get(grid) ?? null;
  let base = gc.width;
  if (pending !== null && pending.column === target) {
    // Same column, still mid-nudge: its pending revert holds the TRUE base and
    // `gc.width` is the inflated value. Take the base and cancel the revert.
    window.clearTimeout(pending.timer);
    base = pending.base;
  }
  gc.width = base + 1;
  const timer = window.setTimeout(() => {
    if (pendingWidthNudge.get(grid)?.timer === timer) pendingWidthNudge.delete(grid);
    try {
      gc.width = base;
    } catch {
      // Grid torn down between the two writes — nothing left to restore.
    }
  }, 0);
  pendingWidthNudge.set(grid, {timer, base, column: target});
}

/** In-flight width nudges, per grid — see {@link forceRepaintAfterVisibilityChange}.
 *  A `WeakMap` so a closed grid takes its entry with it. */
const pendingWidthNudge = new WeakMap<
  DG.Grid, {readonly timer: number; readonly base: number; readonly column: string}
>();

/**
 * Register colour-coding rules on the results DataFrame columns. Mutates
 * column metadata; safe to call once per DataFrame, idempotent.
 *
 * Rules colour only anomalies — "ok" / normal values stay default-colour.
 * Returns the same DataFrame for fluent chaining.
 */
export function configureResultsDataFrame(
  df: DG.DataFrame,
  thresholds: ResultsColorThresholds = DEFAULT_RESULTS_THRESHOLDS,
): DG.DataFrame {
  // Derive the λz ACCEPTANCE verdict + flag columns FIRST, so the band loop below
  // finds ACCEPTANCE and attaches its categorical band (the loop skips a column
  // that is not yet present). Uses the SAME `thresholds` the bands read, so a flag
  // and its column's band cannot disagree.
  writeAcceptanceColumns(df, thresholds);

  // TRIAL (user-directed 2026-09-02): the results grid colours cell TEXT, NOT the background
  // fill — so NO `col.meta.colors` background band is registered on this frame. The per-cell
  // text colour is applied in `configureResultsGrid`'s `onCellPrepare` (gc.style.textColor via
  // matchBandColour), reading the SAME `buildColourRules` bands, so the tier a cell shows and
  // the header legend still come from one source. The legend + header tooltips (derived from
  // buildColourRules directly) are unchanged. To restore fills: re-register the bands here via
  // `col.meta.colors.setCategorical/setConditional(dgColorRules(rule))`.

  // Display headers — static map for every column except CL/VZ (which carry
  // a route-derived label set by the writer; we hoist that to friendlyName).
  for (const [name, friendly] of Object.entries(DISPLAY_HEADERS)) {
    const col = df.col(name);
    if (col !== null) col.setTag('friendlyName', friendly);
  }
  // CL / VZ / VSS carry a writer-set `tags.label` ('CL'/'CL/F', 'Vz', 'Vss')
  // hoisted to the grid header; they are intentionally absent from DISPLAY_HEADERS
  // (the same list backs HEADER_EXEMPT — one source, no hand-kept second copy).
  for (const name of ROUTE_LABELED_COLUMNS) {
    const col = df.col(name);
    if (col !== null) {
      const label = col.tags['label'];
      if (label) col.setTag('friendlyName', headerWithUnit(label, col));
    }
  }

  // Number formats per spec.
  for (const [name, fmt] of Object.entries(NUMBER_FORMATS)) {
    const col = df.col(name);
    if (col !== null) col.meta.format = fmt;
  }

  return df;
}

/**
 * Make a docked result grid read-only: these grids present COMPUTED output
 * (NCA parameters, TK aggregates, F%, dose-proportionality verdicts, sparse
 * results). Hand-editing a cell would desync the display from the validated
 * computation with no recompute path — an audit-integrity hazard (the math is the source of truth,
 * the grid only presents it).
 * `allowEdit = false` keeps selection, sorting, and tooltips while blocking cell
 * edits. Re-applied by each `configure*Grid` (idempotent) so it survives a
 * `grid.dataFrame = newDf` swap on re-Compute.
 */
function setGridReadOnly(grid: DG.Grid): void {
  grid.setOptions({allowEdit: false});
}

/** One row's MRT / Vss fragility verdict plus the live values behind it. */
export interface FragilityReading {
  /** Whether the cross-column tint paints on this row. */
  readonly fragile: boolean;
  /** One phrase per condition that ACTUALLY fired, carrying this row's value
   *  and the threshold it crossed. Empty when nothing fired. */
  readonly reasons: ReadonlyArray<string>;
}

/**
 * Evaluate the MRT / Vss cross-column fragility conditions at one row.
 *
 * **The single source for both the paint and its explanation.** The tint used to
 * be decided here-ish and described in `buildColourRules`; that split is exactly
 * how MRT and Vss shipped tinted with a tooltip that never mentioned the colour.
 * The legend still states the RULE (from the same thresholds object); this adds
 * the row's own numbers, which is what a scientist standing on a tinted cell
 * actually needs.
 *
 * All three siblings are read through `numAt`, which is `isNone`-first by
 * construction. A bare `Column.get()` would return the raw `DG.FLOAT_NULL`
 * sentinel — `2.68e-34`, a FINITE POSITIVE number — so `adjR2 < minR` would be
 * TRUE on a row where λz never fit, reporting "fragile fit, adj. R² = 0.00"
 * where there was no fit at all. That exact break shipped once already, a few
 * lines from here.
 *
 * **Exported for inspection, not for reuse.** `grid.onCellTooltip` hands its
 * callback to the platform and gives nothing back, so the tooltip's *content* is
 * unreachable from a test through the grid. Returning a plain reading — data,
 * not a closure — is what makes AC-4.5's both-directions assertion possible, and
 * it is the same move that made the header tooltips checkable (`RichTipSpec`
 * replacing `() => HTMLElement`, whose stated lesson was precisely "content is
 * only checkable if the artifact is data"). Product code has exactly two callers,
 * both in this file: the paint and the tooltip.
 */
export function readFragility(
  pctAumcCol: DG.Column | null,
  adjR2Col: DG.Column | null,
  nPtsCol: DG.Column | null,
  row: number,
  t: ResultsColorThresholds,
): FragilityReading {
  const pct = pctAumcCol !== null ? numAt(pctAumcCol, row) : NaN;
  const adjR2 = adjR2Col !== null ? numAt(adjR2Col, row) : NaN;
  const nPts = nPtsCol !== null ? numAt(nPtsCol, row) : NaN;
  const highExtrap = Number.isFinite(pct) && pct > t.extrapWarnPctAumc;
  const lowR2 = Number.isFinite(adjR2) && adjR2 < t.lambdaMinRSquared;
  // `<= minP`: a fit AT the floor is fragile; one somehow below it is more so
  // (defensive — sci-comp rejects sub-minimum fits, but do not rely on it).
  const fewPts = Number.isFinite(nPts) && nPts <= t.lambdaMinPoints;

  const reasons: string[] = [];
  if (highExtrap) {
    reasons.push(`%AUMCextrap = ${pct.toFixed(1)}% ` +
      `(> ${fmtThreshold(t.extrapWarnPctAumc)}%)`);
  }
  if (lowR2) {
    reasons.push(`adj. R² = ${adjR2.toFixed(3)} ` +
      `(< ${fmtThreshold(t.lambdaMinRSquared)})`);
  }
  if (fewPts) {
    reasons.push(`N points (λz) = ${nPts.toFixed(0)} ` +
      `(≤ ${fmtThreshold(t.lambdaMinPoints)})`);
  }
  return {fragile: highExtrap || lowR2 || fewPts, reasons};
}

// `FlagState` / `AcceptanceReading` are re-exported from adapter/acceptance.ts.
export type {AcceptanceReading};

/**
 * Evaluate the four λz acceptance criteria at one row (nca-grid-status-spec §2/§3).
 *
 * **Sibling to {@link readFragility}, and the same discipline.** All four source
 * columns are read through `numAt`, which is `isNone`-first — a bare `Column.get()`
 * would return the finite-positive `DG.FLOAT_NULL` sentinel and make `adjR2 < minR`
 * TRUE on a row where λz never fit, reporting a violated criterion where there was
 * no fit at all. Each criterion is `n/a` exactly when ITS source value is null, so a
 * failed λz fit yields four `n/a` flags and the `—` verdict.
 *
 * **Thresholds track live Rules** (user decision): `t` is the same
 * {@link ResultsColorThresholds} the colour bands read, and `SPAN_RATIO_POOR` is
 * the same module constant the SPAN_RATIO band reads (AD-4 — span is NOT a Rules
 * field). So a flag and its column's band cannot disagree. `NEXT_FLAG` uses
 * `≤ minPoints` (not the spec's literal `< minPoints`, which never fires since
 * sci-comp rejects sub-minimum fits) — matching `readFragility` and the peach
 * `N = minP` band. `RSQ_FLAG` fires at `lambdaAcceptanceRSquared` (`R_acc`) — the
 * acceptance threshold the user enters, applied verbatim (`adjR² < R_acc`), the SAME
 * `R_acc` the LAMBDAZ_ADJ_R2 band reads. It is a separate field from the core gate
 * `minRSquared` and NOT derived from it (the removed `goodR`;
 * lambdaz-adjr2-threshold-spec-v2 §2/§4).
 *
 * **Exported for inspection AND reuse.** Three product callers, all in this file:
 * {@link writeAcceptanceColumns} (the column values), the `ACCEPTANCE` cell-tooltip in
 * {@link configureResultsGrid} (the reasons), and {@link buildFacetVerdictDescriptions}
 * (the facet caption). One source feeding all three is what stops any of them describing a
 * verdict the column does not hold.
 */
export function readAcceptance(
  spanCol: DG.Column | null,
  adjR2Col: DG.Column | null,
  nPtsCol: DG.Column | null,
  pctExtrapCol: DG.Column | null,
  row: number,
  t: ResultsColorThresholds,
): AcceptanceReading {
  // Column-reading wrapper (nca-workflow-stages AD-4): the pure predicate lives in
  // adapter/acceptance.ts; this reads the four sources through `numAt` (isNone-first,
  // so a FLOAT_NULL cell is NaN = `n/a`, never a finite 2.68e-34) and delegates.
  return evaluateAcceptance({
    spanRatio: spanCol !== null ? numAt(spanCol, row) : NaN,
    adjR2: adjR2Col !== null ? numAt(adjR2Col, row) : NaN,
    nPoints: nPtsCol !== null ? numAt(nPtsCol, row) : NaN,
    pctExtrap: pctExtrapCol !== null ? numAt(pctExtrapCol, row) : NaN,
  }, t);
}

/** Flag-state → grid glyph. `violated` draws the eye; `met`/`n/a` are quiet. */
const FLAG_GLYPH: Readonly<Record<FlagState, string>> = {
  'violated': FLAG_VIOLATED, 'met': FLAG_MET, 'n/a': FLAG_NA,
};

/**
 * Problem-first sort rank for one results row (item 1 of "surface problem
 * profiles first"). Lower = more severe = higher in the grid. User-ratified
 * order (2026-08-25): failed → Not estimated → Flagged → Accepted.
 *
 * **STATUS is PRIMARY, ACCEPTANCE secondary — and that is not redundant with the
 * verdict alone.** sci-comp couples the two axes: `failed` (no Cmax at all) and
 * `partial` (Cmax/AUClast obtained but λz not fit) both, in practice, carry the
 * `Not estimated` verdict (λz not fit ⇒ every criterion n/a), so the verdict axis
 * cannot separate a genuine engine failure (rank 0) from a computed-but-λz-less
 * profile (rank 1) — only STATUS can. Testing STATUS FIRST also makes the rank
 * robust in the rare/degenerate case a `partial` row's verdict were to differ: it
 * still lands in bucket 1. On an `ok` row λz WAS fit, so its verdict is only ever
 * Accepted / Flagged; the `Not estimated` branch below is therefore defensive.
 */
export function problemSortRank(status: string | null, verdict: string): number {
  if (status === 'failed') return 0;
  if (status === 'partial') return 1;
  // `Not estimated` is defensive here — unreachable on a real `ok` row, whose λz was fit.
  if (verdict === ACCEPT_NA) return 1;
  if (verdict === ACCEPT_FLAGGED) return 2;
  if (verdict === ACCEPT_ACCEPTED) return 3;
  return 4; // unknown verdict → sinks to the bottom
}

/**
 * Compute the ACCEPTANCE verdict + four λz acceptance-flag columns and add them
 * to the per-subject results frame (nca-grid-status-spec §2/§3).
 *
 * **Presentation-side on purpose.** The values depend on thresholds that live at
 * presentation (the live Rules projection — {@link ResultsColorThresholds}), so
 * this runs from {@link configureResultsDataFrame} rather than the writer, and the
 * writer's `FIXED_RESULT_COLUMNS` drift-guard is untouched. The columns are still
 * real DataFrame columns, so they filter / count / export like any other (spec §3
 * requires the reason to be stored, not tooltip-only).
 *
 * Reads the four source columns the criteria are computed from; a degenerate frame
 * missing them yields all-`n/a` flags and the `—` verdict rather than throwing.
 * Returns `df` for fluent chaining.
 */
export function writeAcceptanceColumns(
  df: DG.DataFrame,
  t: ResultsColorThresholds = DEFAULT_RESULTS_THRESHOLDS,
): DG.DataFrame {
  const spanCol = df.col('SPAN_RATIO');
  const adjR2Col = df.col('LAMBDAZ_ADJ_R2');
  const nPtsCol = df.col('LAMBDAZ_NPOINTS');
  const pctExtrapCol = df.col('PCT_EXTRAP');
  // STATUS drives the problem-first rank's top two buckets (failed / partial); read
  // it here so the rank is written in THIS same pass — one source, no drift.
  const statusCol = df.col('STATUS');

  const n = df.rowCount;
  const acceptance = new Array<string>(n);
  const reasons = new Array<string>(n);
  const spanF = new Array<string>(n);
  const rsqF = new Array<string>(n);
  const nextF = new Array<string>(n);
  const extrapF = new Array<string>(n);
  // Problem-first sort rank (item 1). Single-digit strings so the existing
  // `upsertStringColumn` helper carries it and a plain ascending grid sort orders
  // it correctly ('0' < '1' < … < '4'). NOTE the single-digit assumption: lexicographic
  // order equals numeric order only for 0–9, so if `problemSortRank` ever grows past
  // 9 buckets this must zero-pad (else '10' would sort before '2').
  const rank = new Array<string>(n);
  for (let i = 0; i < n; i++) {
    const a = readAcceptance(spanCol, adjR2Col, nPtsCol, pctExtrapCol, i, t);
    // All three verdicts are VISIBLE, banded values — including `Not estimated`
    // (red), which the scientist must see rather than read as an empty cell.
    acceptance[i] = a.verdict;
    reasons[i] = reasonsCellValue(a);
    spanF[i] = FLAG_GLYPH[a.flags.span];
    rsqF[i] = FLAG_GLYPH[a.flags.rsq];
    nextF[i] = FLAG_GLYPH[a.flags.next];
    extrapF[i] = FLAG_GLYPH[a.flags.extrap];
    const st = statusCol !== null && !statusCol.isNone(i) ? statusCol.getString(i) : null;
    rank[i] = String(problemSortRank(st, a.verdict));
  }
  upsertStringColumn(df, ACCEPTANCE_COL, acceptance);
  upsertStringColumn(df, ACCEPT_REASONS_COL, reasons);
  tagReasonPillsColumn(df.col(ACCEPT_REASONS_COL)!);
  upsertStringColumn(df, FLAG_COLS[0], spanF);
  upsertStringColumn(df, FLAG_COLS[1], rsqF);
  upsertStringColumn(df, FLAG_COLS[2], nextF);
  upsertStringColumn(df, FLAG_COLS[3], extrapF);
  upsertStringColumn(df, PROBLEM_RANK_COL, rank);
  return df;
}

/**
 * Re-derive the ACCEPTANCE verdict / flag / problem-rank cells IN PLACE from the
 * (already-patched) numeric λz columns — for the interactive λz override (AC-7.3).
 *
 * Unlike {@link writeAcceptanceColumns} (which REMOVES + re-adds the columns via
 * `upsertStringColumn`, reordering them to the frame tail and dropping the grid's
 * banding/position config), this writes each cell with `col.set`, so the bound grid
 * re-renders the changed verdict without any column churn. Requires the columns to
 * already exist (a full Compute ran) — a no-op otherwise. Reuses {@link readAcceptance}
 * so the patched verdict cannot disagree with the facet caption (rule 5). Does NOT
 * re-sort: PROBLEM_RANK values update but the grid order is re-asserted only on the
 * next full Compute (D10).
 */
export function refreshAcceptanceColumnsInPlace(
  df: DG.DataFrame, t: ResultsColorThresholds = DEFAULT_RESULTS_THRESHOLDS,
): void {
  const acc = df.col(ACCEPTANCE_COL);
  if (acc === null) return; // no full Compute yet — nothing to refresh
  const spanCol = df.col('SPAN_RATIO');
  const adjR2Col = df.col('LAMBDAZ_ADJ_R2');
  const nPtsCol = df.col('LAMBDAZ_NPOINTS');
  const pctExtrapCol = df.col('PCT_EXTRAP');
  const statusCol = df.col('STATUS');
  const reasonsCol = df.col(ACCEPT_REASONS_COL);
  const f0 = df.col(FLAG_COLS[0]);
  const f1 = df.col(FLAG_COLS[1]);
  const f2 = df.col(FLAG_COLS[2]);
  const f3 = df.col(FLAG_COLS[3]);
  const rankCol = df.col(PROBLEM_RANK_COL);
  const n = df.rowCount;
  for (let i = 0; i < n; i++) {
    const a = readAcceptance(spanCol, adjR2Col, nPtsCol, pctExtrapCol, i, t);
    acc.set(i, a.verdict);
    reasonsCol?.set(i, reasonsCellValue(a));
    f0?.set(i, FLAG_GLYPH[a.flags.span]);
    f1?.set(i, FLAG_GLYPH[a.flags.rsq]);
    f2?.set(i, FLAG_GLYPH[a.flags.next]);
    f3?.set(i, FLAG_GLYPH[a.flags.extrap]);
    const st = statusCol !== null && !statusCol.isNone(i) ? statusCol.getString(i) : null;
    rankCol?.set(i, String(problemSortRank(st, a.verdict)));
  }
}

/** Ensure a blank string service column named `name` exists on the results frame
 *  (all-blank when freshly added). Idempotent; returns the column. Shared by the ✎
 *  override flag ({@link LZ_OVERRIDDEN_COL}) and its timestamp sibling
 *  ({@link LZ_OVERRIDDEN_TIME_COL}). */
function ensureBlankStringColumn(df: DG.DataFrame, name: string): DG.Column {
  const existing = df.col(name);
  if (existing !== null) return existing;
  return df.columns.add(
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, name, new Array<string>(df.rowCount).fill('')));
}

/** Display string stored in {@link LZ_OVERRIDDEN_COL} for an overridden subject: the
 *  analyst label, or {@link OVERRIDE_UNKNOWN_USER} when blank so the cell stays
 *  non-empty (the overridden test). */
function overrideFlagText(user: string): string {
  return user !== '' ? user : OVERRIDE_UNKNOWN_USER;
}

/**
 * Project the controller's per-subject λz-override → analyst map onto the results
 * frame's {@link LZ_OVERRIDDEN_COL} (F7 / AC-7.1). Full pass, called once per Compute
 * after the results frame is built. Map MEMBERSHIP is the overridden test (the value
 * may legitimately be `''`); a subject absent from `bySubject` is cleared to `''` —
 * an override that was reset away no longer shows the ✏️.
 */
export function writeSubjectOverrideFlags(
  df: DG.DataFrame, bySubject: ReadonlyMap<string, string>,
  timeBySubject?: ReadonlyMap<string, string>,
): void {
  const col = ensureBlankStringColumn(df, LZ_OVERRIDDEN_COL);
  const timeCol = timeBySubject !== undefined ?
    ensureBlankStringColumn(df, LZ_OVERRIDDEN_TIME_COL) : null;
  const subjCol = df.col('SUBJECT');
  if (subjCol === null) return;
  const n = df.rowCount;
  for (let i = 0; i < n; i++) {
    const key = subjCol.isNone(i) ? '' : subjCol.getString(i);
    const ov = bySubject.has(key);
    col.set(i, ov ? overrideFlagText(bySubject.get(key)!) : '');
    timeCol?.set(i, ov && timeBySubject!.has(key) ? timeBySubject!.get(key)! : '');
  }
}

/**
 * Update ONE subject's {@link LZ_OVERRIDDEN_COL} cell after a single-subject override
 * toggle (F7) — the fast counterpart to {@link writeSubjectOverrideFlags}. `user` is
 * the analyst label (`''` allowed), or `null` to clear the ✏️ when the subject
 * reverted to auto. Caller must repaint the grid (`grid.invalidate()`) so the
 * ACCEPTANCE cell's `onCellPrepare` re-reads this hidden sibling.
 */
export function setSubjectOverrideFlag(
  df: DG.DataFrame, subjectKey: string, user: string | null, time?: string | null,
): void {
  const col = ensureBlankStringColumn(df, LZ_OVERRIDDEN_COL);
  const timeCol = time !== undefined ?
    ensureBlankStringColumn(df, LZ_OVERRIDDEN_TIME_COL) : null;
  const subjCol = df.col('SUBJECT');
  if (subjCol === null) return;
  const n = df.rowCount;
  for (let i = 0; i < n; i++) {
    if (!subjCol.isNone(i) && subjCol.getString(i) === subjectKey) {
      col.set(i, user === null ? '' : overrideFlagText(user));
      timeCol?.set(i, (user === null || time == null) ? '' : time);
      return;
    }
  }
}

/**
 * Facet ordering key per subject, so the small-multiples concentration-time plot
 * follows the SAME problem-first order as the results grid (item 2). Reads the
 * SUBJECT + hidden {@link PROBLEM_RANK_COL} columns of the results frame and returns
 * `subjectKey → rank * 1e6 + rowIndex`: rank as the primary key, the frame's own
 * (subject-lex, determinism-contract) row order as the tiebreak — matching the grid's
 * `[rank, SUBJECT]` sort. The key is `getString`, exactly how the facet builder keys
 * its subjects (`String(raw)`), so the lookup lines up for every non-null subject.
 * A degenerate frame missing either column yields an empty map (facet keeps its
 * natural order).
 */
export function buildFacetSubjectOrder(resultsDf: DG.DataFrame): ReadonlyMap<string, number> {
  const out = new Map<string, number>();
  const subjCol = resultsDf.col('SUBJECT');
  const rankCol = resultsDf.col(PROBLEM_RANK_COL);
  if (subjCol === null || rankCol === null) return out;
  const n = resultsDf.rowCount;
  for (let i = 0; i < n; i++) {
    const key = subjCol.isNone(i) ? '' : subjCol.getString(i);
    const r = rankCol.isNone(i) ? 4 : Number(rankCol.getString(i));
    if (!out.has(key)) out.set(key, r * 1e6 + i);
  }
  return out;
}

/**
 * Per-subject λz acceptance verdict text for the facet small-multiples' viewer
 * DESCRIPTION caption (facet-plot-title-verdict-separators). Sibling to
 * {@link buildFacetSubjectOrder}: reads the SAME `outcome.results` frame and keys by
 * `SUBJECT.getString(i)`, exactly how the facet builder keys its tiles (`String(raw)`),
 * so the lookup lines up per subject.
 *
 * **REUSES {@link readAcceptance}** — the SAME function whose values populate the
 * ACCEPTANCE column and its cell tooltip — so a tile's caption cannot disagree with
 * the grid's verdict (no re-derivation). Returns, per subject:
 *   - Accepted      → `'Accepted'`
 *   - Flagged       → the violated criteria joined by `' · '` (this row's OWN values,
 *                     the same phrasing the ACCEPTANCE tooltip lists — the verdict
 *                     word is dropped here because the caption's job is to name WHY,
 *                     user decision this spike)
 *   - Not estimated → `'Not estimated'` (λz not fit / computation failed)
 *
 * A degenerate frame missing SUBJECT yields an empty map, so a tile with no entry
 * shows no caption rather than a wrong one.
 */
export function buildFacetVerdictDescriptions(
  resultsDf: DG.DataFrame,
  thresholds: ResultsColorThresholds = DEFAULT_RESULTS_THRESHOLDS,
  nonContiguous?: ReadonlySet<string>,
): ReadonlyMap<string, string> {
  const out = new Map<string, string>();
  const subjCol = resultsDf.col('SUBJECT');
  if (subjCol === null) return out;
  // Same source columns readAcceptance reads in configureResultsGrid — one set of
  // names, so the caption and the grid cell evaluate the identical criteria.
  const spanCol = resultsDf.col('SPAN_RATIO');
  const adjR2Col = resultsDf.col('LAMBDAZ_ADJ_R2');
  const nPtsCol = resultsDf.col('LAMBDAZ_NPOINTS');
  const pctExtrapCol = resultsDf.col('PCT_EXTRAP');
  const n = resultsDf.rowCount;
  for (let i = 0; i < n; i++) {
    const key = subjCol.isNone(i) ? '' : subjCol.getString(i);
    if (out.has(key)) continue;
    const a = readAcceptance(spanCol, adjR2Col, nPtsCol, pctExtrapCol, i, thresholds);
    // Flagged → the reasons alone (user decision); Accepted / Not estimated → the
    // verdict word itself (a.verdict is already 'Accepted' / 'Not estimated').
    const base = a.verdict === ACCEPT_FLAGGED ? a.reasons.join(' · ') : a.verdict;
    // Decision 2: a non-contiguous manual λz selection appends a caption note, so the
    // analyst editing on the tile sees it immediately (the same fact is in the
    // WARNINGS_JSON stream). Additive to whatever verdict text stands.
    out.set(key, nonContiguous?.has(key) ? `${base} · non-contiguous λz` : base);
  }
  return out;
}

/**
 * Per-subject facet tile COLOUR keyed to the ACCEPTANCE verdict (experiment
 * `facet-lambdaz-line-points`): the tile reads at a glance as the same three-way
 * disposition as the grid's ACCEPTANCE cell — Accepted → teal, Flagged → amber,
 * Not estimated → red, the Radix "teal-success" triad from {@link ./status-palette}
 * (via `HUE_ACCEPTED/FLAGGED/NOT_ESTIMATED`), NOT the pale `COLOR_OK/WARNING/CRITICAL`
 * cell-band tints which read too faint as a line/marker hue. Sibling to
 * {@link buildFacetVerdictDescriptions}: same
 * `outcome.results` frame, same `readAcceptance`, same `SUBJECT` keying, so a tile's
 * colour cannot disagree with its verdict caption or its grid cell. A degenerate
 * frame missing SUBJECT yields an empty map (caller falls back to `STATUS_NEUTRAL`
 * — a tile with no verdict must not wear a colour that reads as one).
 */
export function buildFacetVerdictColors(
  resultsDf: DG.DataFrame,
  thresholds: ResultsColorThresholds = DEFAULT_RESULTS_THRESHOLDS,
): ReadonlyMap<string, number> {
  const out = new Map<string, number>();
  const subjCol = resultsDf.col('SUBJECT');
  if (subjCol === null) return out;
  const spanCol = resultsDf.col('SPAN_RATIO');
  const adjR2Col = resultsDf.col('LAMBDAZ_ADJ_R2');
  const nPtsCol = resultsDf.col('LAMBDAZ_NPOINTS');
  const pctExtrapCol = resultsDf.col('PCT_EXTRAP');
  const n = resultsDf.rowCount;
  for (let i = 0; i < n; i++) {
    const key = subjCol.isNone(i) ? '' : subjCol.getString(i);
    if (out.has(key)) continue;
    const a = readAcceptance(spanCol, adjR2Col, nPtsCol, pctExtrapCol, i, thresholds);
    // Accepted → teal, Flagged → amber, Not estimated → red. The SAME constants the
    // results-grid cell TEXT uses (`fillToTextColour`), so grid and facet cannot drift.
    const color = a.verdict === ACCEPT_ACCEPTED ? HUE_ACCEPTED :
      a.verdict === ACCEPT_FLAGGED ? HUE_FLAGGED : HUE_NOT_ESTIMATED;
    out.set(key, color);
  }
  return out;
}

/** Per-subject λz terminal-line parameters for the facet tiles (spike
 *  `facet-lambdaz-line-points`, VR-10). The line is `ln(C) = intercept − λz·t`
 *  over the terminal window `[tFirst, tLast]` — every field is already a fixed
 *  results-frame column (`results-writer.ts` FIXED_RESULT_COLUMNS). */
export interface FacetLambdaZParams {
  readonly lambdaZ: number;
  readonly intercept: number;
  readonly tFirst: number;
  readonly tLast: number;
}

/**
 * Per-subject λz line parameters for the facet small-multiples (VR-10 static
 * overlay). Sibling to {@link buildFacetVerdictDescriptions} / {@link buildFacetSubjectOrder}:
 * reads the SAME `outcome.results` frame and keys by `SUBJECT.getString(i)` (first-seen
 * wins), exactly how the facet builder keys its tiles, so the lookup lines up per subject.
 * A subject whose fit failed carries `NaN` params (no line drawn); a degenerate frame
 * missing any source column yields an empty map.
 */
export function buildFacetLambdaZParams(
  resultsDf: DG.DataFrame,
): ReadonlyMap<string, FacetLambdaZParams> {
  const out = new Map<string, FacetLambdaZParams>();
  const subjCol = resultsDf.col('SUBJECT');
  const lzCol = resultsDf.col('LAMBDA_Z');
  const interceptCol = resultsDf.col('LAMBDAZ_INTERCEPT');
  const tFirstCol = resultsDf.col('LAMBDAZ_TFIRST');
  const tLastCol = resultsDf.col('LAMBDAZ_TLAST');
  if (subjCol === null || lzCol === null || interceptCol === null ||
      tFirstCol === null || tLastCol === null) return out;
  const n = resultsDf.rowCount;
  for (let i = 0; i < n; i++) {
    const key = subjCol.isNone(i) ? '' : subjCol.getString(i);
    if (out.has(key)) continue;
    out.set(key, {
      lambdaZ: numAt(lzCol, i),
      intercept: numAt(interceptCol, i),
      tFirst: numAt(tFirstCol, i),
      tLast: numAt(tLastCol, i),
    });
  }
  return out;
}

/**
 * Apply the default problem-first row sort to a results grid (item 1): the hidden
 * {@link PROBLEM_RANK_COL} (0 = most severe, leads) then SUBJECT — the
 * determinism-contract tiebreak (user decision 2026-08-25). Re-asserted on every
 * Compute ({@link configureResultsGrid} calls it, and the mapping-pane re-applies it
 * after the grouping-column pass), so a scientist's manual header sort survives a
 * Compact/Full toggle but a fresh Compute returns to problem-first. Guarded: a
 * degenerate frame without the rank column is left unsorted rather than throwing.
 * Sorting BY a hidden column is a supported, tested property of this grid
 * (`results-grid-platform-tests` AC-1.2).
 */
export function applyResultsDefaultSort(grid: DG.Grid): void {
  const df = grid.dataFrame;
  if (df === null || df === undefined || df.col(PROBLEM_RANK_COL) === null) return;
  grid.sort([PROBLEM_RANK_COL, 'SUBJECT']);
}

/** The MRT / Vss cross-column fragility tooltip spec — ONE source shared by the results grid cell
 *  tooltip and the λz context-panel value tooltip ({@link paramValueTooltipSpec}), so the "fragile
 *  estimate" wording cannot drift between the two (architect review). `column` is 'MRT' or 'VSS'. */
function fragilityTooltipSpec(column: string, f: FragilityReading): RichTipSpec {
  return {
    title: `${column === 'MRT' ? 'MRT' : 'Vss'} — fragile estimate`,
    summary: f.reasons.join(' · '),
    detail: 'Moment parameters carry double λz sensitivity (the 1/λz² tail), so a fragile ' +
      'terminal fit moves them further than it moves CL or Vz.',
  };
}

/**
 * The ACCEPTANCE cell tooltip SPEC for row `i` — the SINGLE source shared by the results grid's
 * `onCellTooltip` and the λz context-panel status line, so the two cannot disagree (user-directed
 * 2026-09-02). Flagged lists the violated criteria (this row's values), Accepted confirms, Not
 * estimated says λz was not fitted; a manual override appends the analyst.
 */
export function buildAcceptanceTooltipSpec(
  resultsDf: DG.DataFrame, i: number,
  thresholds: ResultsColorThresholds = DEFAULT_RESULTS_THRESHOLDS,
): RichTipSpec {
  const a = readAcceptance(
    resultsDf.col('SPAN_RATIO'), resultsDf.col('LAMBDAZ_ADJ_R2'),
    resultsDf.col('LAMBDAZ_NPOINTS'), resultsDf.col('PCT_EXTRAP'), i, thresholds);
  const oc = resultsDf.col(LZ_OVERRIDDEN_COL);
  const ov = oc !== null && !oc.isNone(i) ? oc.getString(i) : '';
  const tc = resultsDf.col(LZ_OVERRIDDEN_TIME_COL);
  const when = tc !== null && !tc.isNone(i) ? formatOverrideTimestamp(tc.getString(i)) : '';
  const overrideExtra: HTMLElement[] = ov !== '' ? [
    ui.divText(`✏️ Manual λz override — ${ov}`, 'nca-richtip-summary'),
    ...(when !== '' ? [ui.divText(when, 'nca-richtip-summary')] : []),
  ] : [];
  if (a.verdict === ACCEPT_FLAGGED) {
    return {
      title: 'Acceptance — flagged',
      summary: 'Violated:',
      extra: [...a.reasons.map((r) => ui.divText(`• ${r}`, 'nca-richtip-summary')), ...overrideExtra],
    };
  }
  const summary = a.verdict === ACCEPT_ACCEPTED ?
    'All λz acceptance criteria met.' : 'λz not estimated — criteria not evaluated.';
  return {title: 'Acceptance', summary, extra: overrideExtra.length > 0 ? overrideExtra : undefined};
}

/**
 * The `C0_EXTRAP` CELL tooltip spec for row `i`, or `null` when the cell has nothing beyond the
 * header to say (DG's default then applies). ONE extra line — `back-extrapolated: 20.6 % of
 * AUC₀–∞` — read from `PROVENANCE_JSON.c0.pctAucBackExtrap` (sci-comp ≥ 0.11.0, the Phoenix
 * `AUC_%Back_Ext` analogue), shown ONLY when the value was estimated (`method !== 'observed'`:
 * a measured dose-time sample has nothing extrapolated — the line is omitted rather than
 * printing `0 %`) and the share exists (omitted on a `partial` row, where there is no AUC₀–∞).
 * No column, no band, no threshold: a diagnostic the scientist reads on the cell that carries
 * the number (GROK-20960). Exported for the presentation test.
 */
export function buildC0TooltipSpec(
  resultsDf: DG.DataFrame, i: number, header: RichTipSpec,
): RichTipSpec | null {
  const pc = resultsDf.col('PROVENANCE_JSON');
  if (pc === null || pc.isNone(i)) return null;
  interface C0Json {method?: unknown; pctAucBackExtrap?: unknown}
  const readC0 = (): C0Json | null => {
    try {
      const parsed = JSON.parse(pc.getString(i)) as {c0?: C0Json | null};
      return parsed.c0 ?? null;
    } catch {
      return null;
    }
  };
  const c0 = readC0();
  if (c0 === null || c0.method === 'observed') return null;
  const pct = c0.pctAucBackExtrap;
  if (typeof pct !== 'number' || !Number.isFinite(pct)) return null;
  return {
    ...header,
    extra: [ui.divText(`back-extrapolated: ${pct.toFixed(1)} % of AUC₀–∞`, 'nca-richtip-summary')],
  };
}

/** The columns that rest on the profile's FIRST-row body weight when the weight varies (AD-2). */
const BW_VARIES_COMPANIONS: ReadonlySet<string> = new Set(['CL_ABS', 'VZ_ABS', 'VSS_ABS', 'BW_KG']);

/**
 * The tooltip SPEC of an absolute companion / `BW_KG` cell whose profile's body weight varies
 * (nca-cl-bw-normalised-units AD-2): the header explanation plus the profile's own
 * `BW_VARIES_WITHIN_PROFILE` note from `WARNINGS_JSON` (the range and the weight used), so the
 * scientist reads it on the number it qualifies. `null` (DG's default) when the note is absent.
 * Exported for the presentation test.
 */
export function buildBwVariesTooltipSpec(
  resultsDf: DG.DataFrame, column: string, i: number, header: RichTipSpec | undefined,
): RichTipSpec | null {
  if (!BW_VARIES_COMPANIONS.has(column)) return null;
  const wc = resultsDf.col('WARNINGS_JSON');
  if (wc === null || wc.isNone(i)) return null;
  let note: string | null = null;
  try {
    const parsed = JSON.parse(wc.getString(i)) as Array<{code?: unknown; message?: unknown}>;
    const w = Array.isArray(parsed) ? parsed.find((x) => x?.code === 'BW_VARIES_WITHIN_PROFILE') : undefined;
    note = typeof w?.message === 'string' ? w.message : null;
  } catch {
    return null;
  }
  if (note === null) return null;
  return {...(header ?? {title: column, summary: ''}), extra: [ui.divText(note, 'nca-richtip-summary')]};
}

/**
 * The per-VALUE tooltip SPEC for `column` at row `i` — what to say about THIS value, distinct from
 * the header tooltip that explains the parameter. For MRT / VSS: the cross-column fragility (this
 * row's reasons), matching the grid cell. For a banded column: the LEGEND of the tier the value
 * falls in (why it is coloured so). Otherwise the header explanation. User-directed 2026-09-02
 * (λz context-panel value tooltips).
 */
export function paramValueTooltipSpec(
  resultsDf: DG.DataFrame, column: string, i: number,
  thresholds: ResultsColorThresholds = DEFAULT_RESULTS_THRESHOLDS,
): RichTipSpec {
  const header = buildHeaderTooltips(thresholds)[column];
  const title = header?.title ?? column;
  if (column === 'MRT' || column === 'VSS') {
    const c = resultsDf.col(column);
    const own = c !== null && !c.isNone(i) ? numAt(c, i) : NaN;
    if (Number.isFinite(own)) {
      const f = readFragility(
        resultsDf.col('PCT_AUMC_EXTRAP'), resultsDf.col('LAMBDAZ_ADJ_R2'),
        resultsDf.col('LAMBDAZ_NPOINTS'), i, thresholds);
      if (f.fragile) return fragilityTooltipSpec(column, f);
    }
    return header ?? {title, summary: ''};
  }
  const rule = buildColourRules(thresholds)[column];
  if (rule !== undefined && rule.kind !== 'crossColumn') {
    const c = resultsDf.col(column);
    const raw = c !== null && !c.isNone(i) ? numAt(c, i) : NaN;
    const band = rule.bands.find((cb) => bandContains(cb.band, raw));
    if (band !== undefined) return {title, summary: legendText(band)};
  }
  return header ?? {title, summary: ''};
}

/**
 * Apply grid-side presentation to a results Grid:
 *   1. Hide the JSON audit columns by default.
 *   2. Install rich header tooltips that explain each coloured column's rules.
 *
 * Background colour coding works through Datagrok's default fill rendering —
 * no per-column toggle is required; `setConditional` / `setCategorical` calls
 * on the DataFrame side already paint the cell backgrounds.
 *
 * Caller must have already invoked {@link configureResultsDataFrame} on the
 * Grid's DataFrame so the colour rules exist before this function runs.
 *
 * Safe to re-invoke (e.g. after `existingGrid.dataFrame = newDf` on re-Compute):
 *   - `grid.col(name).visible = false` is idempotent.
 *   - `grid.onCellTooltip(handler)` REPLACES the previous handler, not appends.
 */
export function configureResultsGrid(
  grid: DG.Grid,
  thresholds: ResultsColorThresholds = DEFAULT_RESULTS_THRESHOLDS,
  extended = false,
): void {
  setGridReadOnly(grid);
  // Order + visibility in ONE write-point (the two drifted apart on the study
  // grid when they were set in separate places). The JSON audit columns are not
  // in either profile's visible order, so they stay hidden in both — unchanged
  // behaviour, now a consequence of the profile rather than a separate loop.
  applyResultsColumnProfile(grid, extended);
  watchSemanticTypes(grid);

  // Narrow overrides for columns DG's content-based auto-size gets wrong. Same
  // idiom as `applyFixedGridSpec`; guarded, since a degenerate frame may lack the
  // column. Placed after the profile so it also covers a column the profile just
  // revealed, and idempotent across the re-Compute `dataFrame =` swap.
  for (const [name, w] of Object.entries(GRID_WIDTHS)) {
    const gc = grid.col(name);
    if (gc !== null) gc.width = w;
  }

  // F2 — pin the subject so identity survives horizontal scroll. Requirement (1),
  // and load-bearing rather than a nicety: the compact profile is a ~40%
  // reduction, NOT a no-scroll guarantee (a real cassette with Compound/Dose/Sex
  // still lands around 20 columns), so a scientist reading CL at the right edge
  // needs to know whose CL it is.
  //
  // Verified live on THIS host before shipping (the embedded `Viewer.grid`, which
  // is the only one the app reaches): pinning relocates a middle column from
  // index 12 to 1, survives the re-Compute `dataFrame =` swap, and leaves every
  // banded column's per-row colours byte-identical across pin/unpin. That last
  // check is why this is safe to place beside the colour bands at all — DG has a
  // documented pin() formatting-loss defect, fixed in 1.23 / PowerGrid 1.5.0
  // (installed here: 1.27.7 / 1.7.0).
  //
  // Guarded, never a bare assertion: a degenerate result set may lack the column.
  // Idempotent — this function re-runs on every Compute.
  grid.col('SUBJECT')?.pin();

  // `Violated` pill fill — re-keyed per Compute, since the labels carry the live thresholds.
  applyReasonPillColours(grid, thresholds);

  // Default problem-first row order (item 1): failed → Not estimated → Flagged →
  // Accepted, lex-by-subject within a block. Orthogonal to the pin above (row order
  // vs frozen column) and re-asserted on every Compute (Q2). Placed after pin so the
  // frozen SUBJECT rides at index 0 regardless of the sort.
  applyResultsDefaultSort(grid);

  const tooltips = buildHeaderTooltips(thresholds);
  const df = grid.dataFrame;
  const pctAumcCol = df?.col('PCT_AUMC_EXTRAP') ?? null;
  const adjR2Col = df?.col('LAMBDAZ_ADJ_R2') ?? null;
  const nPtsCol = df?.col('LAMBDAZ_NPOINTS') ?? null;
  /** The row's fragility reading — ONE function, read by both the paint below
   *  and the tooltip above, from the same `thresholds` object. That is what
   *  makes it impossible for the tooltip to describe a rule the code does not
   *  apply, which is the invariant this module exists to hold. */
  const fragilityAt = (i: number): FragilityReading =>
    readFragility(pctAumcCol, adjR2Col, nPtsCol, i, thresholds);

  // ONE handler. `grid.onCellTooltip` REPLACES the previous handler rather than
  // appending (documented above), so registering a second one for the cell
  // branch would silently destroy every header tooltip in this grid.
  grid.onCellTooltip((cell, x, y) => {
    // F4 — the MRT / Vss fragility tint explains itself with THIS ROW's actual
    // numbers. The header legend can only state the rule ("adj. R² < 0.85");
    // standing on a tinted cell, the question is which condition fired here and
    // by how much. Placed before the header early-return, which returned false
    // for every table cell.
    if (cell.isTableCell) {
      const nm = cell.tableColumn?.name;
      const i = cell.cell.rowIndex;
      if (i < 0) return false;
      // Review — the state, why, reasons, who/when (D15), inside this ONE handler.
      if (nm === REVIEW_COL) {
        if (df === null) return false;
        const spec = buildReviewTooltipSpec(df, i);
        if (spec === null) return false;
        ui.tooltip.show(richTipEl(spec), x, y);
        return true;
      }
      // Pipeline — the run's own reasons + who/when (F10), inside this ONE handler.
      if (nm === PIPELINE_COL) {
        const spec = pipelineTipAt(df, i);
        if (spec === null) return false;
        ui.tooltip.show(richTipEl(spec), x, y);
        return true;
      }
      // An EXCLUDED profile's λz-dependent value (F-S4.5): kept visible, marked, explained.
      if (nm != null && isLzDependentParam(nm) && df !== null) {
        const info = reviewInfoAt(df, i);
        if (info?.state === REVIEW_EXCLUDED && info.decision !== null) {
          ui.tooltip.show(richTipEl(excludedCellTipSpec(info.decision)), x, y);
          return true;
        }
      }
      // ACCEPTANCE — a Flagged cell explains WHICH criteria fired, with this row's
      // own values (spec §2). A `—` cell says why it could not be judged; an
      // Accepted cell falls through to DG's default (nothing to explain).
      // EVERY ACCEPTANCE cell gets a tooltip (not just Flagged): Flagged lists the
      // violated criteria with this row's values, Accepted confirms it passed, a
      // blank (null) cell says λz was not evaluated.
      // `Violated` pills explain themselves with the SAME tip as the verdict (the full criteria);
      // an empty cell (nothing violated) keeps DG's default.
      if (nm === ACCEPTANCE_COL ||
          (nm === ACCEPT_REASONS_COL && cell.cell.value != null && cell.cell.value !== '')) {
        // ONE source for this spec — {@link buildAcceptanceTooltipSpec} — reused verbatim by the
        // λz context-panel status line, so the grid cell and the panel tooltips cannot disagree
        // (user-directed 2026-09-02). It reads the SAME readAcceptance + `_lz_overridden`.
        if (df === null) return false;
        ui.tooltip.show(richTipEl(buildAcceptanceTooltipSpec(df, i, thresholds)), x, y);
        return true;
      }
      // C0_EXTRAP — the header explanation plus ONE per-row line (the back-extrapolated share
      // of AUC₀–∞) when the value was estimated; a measured / partial row keeps DG's default.
      if (nm === 'C0_EXTRAP') {
        if (df === null) return false;
        const spec = buildC0TooltipSpec(df, i, tooltips.C0_EXTRAP);
        if (spec === null) return false;
        ui.tooltip.show(richTipEl(spec), x, y);
        return true;
      }
      // Absolute companions + BW_KG — the profile's varying-weight note (AD-2), on the values it
      // qualifies; a profile without the note keeps DG's default.
      if (nm != null && BW_VARIES_COMPANIONS.has(nm)) {
        if (df === null) return false;
        const spec = buildBwVariesTooltipSpec(df, nm, i, tooltips[nm]);
        if (spec === null) return false;
        ui.tooltip.show(richTipEl(spec), x, y);
        return true;
      }
      if (nm !== 'MRT' && nm !== 'VSS') return false;
      // Same emptiness guard as the paint: an untinted cell keeps DG's default
      // behaviour rather than gaining a tooltip that explains a tint it has not
      // got. `GridCell.value` is the one read API that honours null.
      const own = cell.cell.value;
      if (own === null || (typeof own === 'number' && !Number.isFinite(own))) return false;
      const f = fragilityAt(i);
      if (!f.fragile) return false;
      // ONE source with the context-panel value tooltip — {@link fragilityTooltipSpec} (architect
      // review dedup): the "fragile estimate" wording cannot drift between grid cell and panel.
      ui.tooltip.show(richTipEl(fragilityTooltipSpec(nm, f)), x, y);
      return true;
    }
    if (!cell.isColHeader)
      return false;
    const name = cell.tableColumn?.name;
    if (name === undefined || name === null)
      return false;
    const spec = tooltips[name];
    if (spec === undefined)
      return false;
    ui.tooltip.show(richTipEl(spec), x, y);
    return true;
  });

  // Cell colouring lives in ONE `onCellPrepare` (DG REPLACES, not appends — a second
  // registration would silently wipe this one, taking the ✏️ + header/fragility handling with
  // it). TRIAL (user-directed 2026-09-02): banded columns colour the cell TEXT, not the
  // background — so `configureResultsDataFrame` registers no `col.meta.colors` fill and the
  // tiers are applied here via `gc.style.textColor`. The text tiers, the MRT/VSS legend text
  // in `buildColourRules`, and the per-row tooltip all read the SAME `thresholds`; the
  // per-column tint and the fragility tint now share `buildColourRules` / `readFragility`
  // respectively, so a cell's colour and its explanation cannot disagree.
  //
  // AC-B3 (P4/P5): MRT/VSS carry DOUBLE λz sensitivity (the 1/λz² extrapolation tail), so a
  // fragile terminal fit hits them harder than CL/Vz — hence the cross-column Warning tint
  // (now on the TEXT) when %AUMCextrap is high OR the λz fit is fragile.
  const colourRules = buildColourRules(thresholds);
  grid.onCellPrepare((gc) => {
    if (!gc.isTableCell)
      return;
    const nm = gc.gridColumn?.name;
    const i = gc.cell.rowIndex;
    if (i < 0 || nm == null)
      return;
    // "To review" is a LINK (user-directed 2026-09-30): clicking it starts the review of that profile
    // (the view wires the click). Link blue, not a severity hue — the verdict colour stays on
    // ACCEPTANCE. Every other Review state keeps the default neutral text.
    if (nm === REVIEW_COL) {
      if (gc.cell.value === REVIEW_TO_REVIEW) gc.style.textColor = REVIEW_LINK_TEXT;
      return;
    }
    if (nm === ACCEPT_REASONS_COL) {
      reasonPillFont(gc);
      return;
    }
    // F-S4.5 (approved design): an EXCLUDED profile's λz-dependent value stays visible but reads
    // as not reportable — muted text + ` (excl.)` via `customText` (display-only: the VALUE, sort,
    // filter and band are untouched — the ✎ precedent). Returns before the band colour so the tier
    // hue never paints a value that is not reported.
    if (isLzDependentParam(nm) && gc.cell.value !== null && df !== null) {
      const rc = df.col(REVIEW_COL);
      if (rc !== null && !rc.isNone(i) && rc.getString(i) === REVIEW_EXCLUDED) {
        const v = gc.cell.value;
        if (typeof v === 'number' && Number.isFinite(v)) {
          const col = df.col(nm);
          // One line, truncated at the END (user-directed 2026-09-30): the default wrapped the marker
          // mid-token in narrow columns (`1.961 (` / `excl.)`).
          paintExcludedCell(gc, col !== null ? col.getString(i) : String(v));
          return;
        }
      }
    }
    // Manual-override indicator (F7 / AC-7.1). A subject whose λz fit was hand-edited shows a
    // pencil AFTER its verdict, WITHOUT changing the ACCEPTANCE cell VALUE (`customText` is
    // display-only; the value-keyed logic still reads the bare verdict). The truth is the hidden
    // `_lz_overridden` sibling. The glyph is the MONOCHROME pencil `✎` (U+270E), NOT the emoji
    // `✏️` (U+270F + FE0F): an emoji renders in its own multicolour glyph and IGNORES the cell
    // text colour, whereas `✎` inherits `gc.style.textColor`, so the marker takes the verdict
    // colour of the cell font (user-directed 2026-09-02). Falls THROUGH to the per-column text
    // colour below — ACCEPTANCE is itself a (categorical) banded column.
    if (nm === ACCEPTANCE_COL) {
      const oc = df?.col(LZ_OVERRIDDEN_COL) ?? null;
      const ov = oc !== null && !oc.isNone(i) ? oc.getString(i) : '';
      if (ov !== '') gc.customText = `${gc.cell.value ?? ''} ✎`;
    }
    // Cross-column MRT/VSS fragility → TEXT colour (was a background tint). An empty cell
    // (VSS null on non-IV rows) carries nothing to flag — the guard that also makes hiding an
    // all-null VSS legitimate under the compact profile (the tint can never paint).
    if (nm === 'MRT' || nm === 'VSS') {
      const own = gc.cell.value;
      if (own === null || (typeof own === 'number' && !Number.isFinite(own)))
        return;
      if (fragilityAt(i).fragile)
        gc.style.textColor = HUE_FLAGGED;
      return;
    }
    // Per-column banded columns (STATUS / ACCEPTANCE / %extrap / adj-R² / N / span /
    // %AUMCextrap): colour the cell TEXT by tier, from the SAME bands as the legend
    // (matchBandColour → fillToTextColour). A non-banded column matches nothing → default text.
    const rule = colourRules[nm];
    if (rule === undefined || rule.kind === 'crossColumn')
      return;
    const fill = matchBandColour(rule, gc.cell.value);
    if (fill !== null)
      gc.style.textColor = fillToTextColour(fill);
    // Bold the two VERDICT columns — Computation (STATUS) + Acceptance — so the status word
    // stands out (user-directed 2026-09-02). Same bold idiom as the aggregated grid's fragile-N
    // cells. The numeric banded columns stay regular weight.
    if (nm === 'STATUS' || nm === ACCEPTANCE_COL)
      boldGridCell(gc);
  });
}

// ---------------------------------------------------------------------------
// Aggregated ("By group") grid presentation.
//
// The aggregator (results-aggregator.ts) emits MEAN_<P> / SD_<P> / N_<P> /
// CV_PCT_<P> columns plus a group-size `N`. DG's groupBy carries the group-key
// columns' tags through (so ROUTE / CL_BASIS / DOSE_NORM already show their
// per-subject friendlyName), but the freshly-built stat columns have no tags —
// hence raw "MEAN_CMAX" headers and only the default stats tooltip. This block
// gives the stat columns friendly headers, mirrored number formats, and header
// tooltips explaining each one.
// ---------------------------------------------------------------------------

/** Per-parameter display name + one-line gloss for the aggregated grid.
 *  Kept separate from the per-subject `buildHeaderTooltips` text so the validated
 *  per-subject tooltips are not touched (user decision, spike/groupby-grid-presentation).
 *  The `friendly` strings MIRROR {@link DISPLAY_HEADERS} so the aggregated
 *  "Mean C_max" matches the per-subject "C_max"; keep the two in sync if either
 *  changes. Covers DEFAULT_AGG_PARAMS; any param absent here falls back to its raw
 *  name with no gloss. EXPORTED so the study-overview unified grid
 *  (nca-study-overview.ts) — another Mean/%CV/n aggregated grid — draws its
 *  parameter-gloss tooltips from this ONE source rather than a parallel table
 *  (Rule 5; study-grid-header-tooltips spike). The two grids share the gloss
 *  DEFINITION; each renders the friendly token in its own label convention. */
export const AGG_PARAM: Readonly<Record<string, {friendly: string; gloss: string}>> = {
  CMAX: {friendly: 'C_max', gloss: 'maximum observed plasma concentration'},
  TMAX: {friendly: 't_max', gloss: 'time of the maximum observed concentration'},
  AUC_LAST: {friendly: 'AUC_0–t', gloss: 'AUC to the last measurable concentration'},
  AUC_INF: {friendly: 'AUC_0–∞', gloss: 'AUC extrapolated to infinity'},
  PCT_EXTRAP: {friendly: '%AUCextrap', gloss: 'share of AUC₀–∞ obtained by extrapolation'},
  LAMBDA_Z: {friendly: 'λz', gloss: 'terminal elimination rate constant'},
  HALF_LIFE: {friendly: 't_½', gloss: 'terminal half-life, ln(2)/λz'},
  CL: {friendly: 'CL', gloss:
    'clearance per kg of body weight (Dose/AUC₀–∞ ÷ BW, mL/min/kg); CL vs CL/F varies per ' +
    'subject — a group that pools both withholds the mean (shown as NE)'},
  VZ: {friendly: 'Vz', gloss: 'terminal-phase volume of distribution per kg (L/kg)'},
  AUMC_INF: {friendly: 'AUMC_0–∞', gloss: 'first-moment area extrapolated to infinity'},
  MRT: {friendly: 'MRT', gloss: 'mean residence time, AUMC₀–∞/AUC₀–∞ (− T_inf/2 for infusion)'},
  VSS: {friendly: 'Vss', gloss: 'steady-state volume of distribution per kg (IV only, L/kg)'},
  CL_ABS: {friendly: 'CL abs.', gloss: 'absolute clearance, not divided by body weight (L/h)'},
  VZ_ABS: {friendly: 'Vz abs.', gloss: 'absolute terminal-phase volume of distribution (L)'},
  VSS_ABS: {friendly: 'Vss abs.', gloss: 'absolute steady-state volume of distribution (IV only, L)'},
};

/** Aggregate stat prefixes — header word + tooltip sentence + value format.
 *  Prefixes are mutually non-prefixing, so first match wins when parsing. */
const AGG_STATS: ReadonlyArray<{
  prefix: string;
  word: string;
  sentence: (friendly: string) => string;
  format: (param: string) => string | null;
}> = [
  {prefix: 'MEAN_', word: 'Mean',
    sentence: (f) => `Arithmetic mean of ${f} across the subjects in this group.`,
    format: (p) => NUMBER_FORMATS[p] ?? null},
  {prefix: 'SD_', word: 'SD',
    sentence: (f) => `Sample standard deviation of ${f} across the subjects in this group.`,
    format: (p) => NUMBER_FORMATS[p] ?? null},
  {prefix: 'CV_PCT_', word: '%CV',
    sentence: (f) =>
      `Coefficient of variation of ${f}: 100 × SD / mean (unitless; high values ` +
      `flag a heterogeneous group).`,
    format: () => '0.0'},
  {prefix: 'N_', word: 'N',
    sentence: (f) =>
      `Number of subjects contributing a non-null ${f} to this group ` +
      `(≤ the group's N when some subjects lack this parameter).`,
    format: () => null},
];

/** The aggregate-stat tooltip sentence for a stat `word` (`'Mean'` | `'SD'` |
 *  `'%CV'` | `'N'`), phrased for `friendly`, or `null` for an unknown word.
 *  Exported so the study-overview grid (which tags stats as `Mean` / `%CV` / `n`)
 *  can reuse the SAME sentences as the aggregated grid — mapping its lowercase
 *  `n` to `N` at the call site (study-grid-header-tooltips spike, Rule 5). */
export function aggStatSentence(word: string, friendly: string): string | null {
  const stat = AGG_STATS.find((s) => s.word === word);
  return stat === undefined ? null : stat.sentence(friendly);
}

/** Plain header tooltips for the aggregated grid's group-key / size columns.
 *  Carried metadata group columns (compound, dose, sex, …) are self-describing
 *  user column names and intentionally get no tooltip. */
const AGG_GROUP_TOOLTIPS: Readonly<Record<string, {title: string; desc: string}>> = {
  N: {title: 'N', desc:
    'Group size — number of subjects in this group (cohort N). A per-parameter ' +
    'N (e.g. "N C_max") may be smaller when some subjects lack that parameter.'},
  ROUTE: {title: 'Route', desc: 'Route of administration (grouping key).'},
  CL_BASIS: {title: 'CL basis', desc:
    'Clearance basis of the group: "CL" (true, intravascular) or "CL/F" ' +
    '(apparent, extravascular). Group by CL basis to separate true and apparent CL — ' +
    'a group that pools both withholds CL, Vz, Vss, MRT and AUMC (shown as empty/NE).'},
  DOSE_NORM: {title: 'Dose norm.', desc:
    'How the dose of the group was passed to the computation: "absolute" (mg) or "per_kg" ' +
    '(a mg/kg dose with no body weight). CL, Vz and Vss are per kg either way, so a group ' +
    'pooling both reports one mean.'},
  STATUS: {title: 'Status', desc: 'NCA computation status of the subjects in this group.'},
};

/** Parse an aggregate column name into its stat + source parameter, or null
 *  when the name is not a stat column (group key, plain `N`, metadata). */
function parseAggColumn(
  name: string,
): {stat: (typeof AGG_STATS)[number]; param: string} | null {
  for (const stat of AGG_STATS) {
    if (name.startsWith(stat.prefix))
      return {stat, param: name.slice(stat.prefix.length)};
  }
  return null;
}

/** Resolve a stat column's display names, PREFERRING the aggregator-set `label`
 *  tag over the static {@link AGG_PARAM} name. The aggregator writes
 *  `label = "${statWord} ${sourceLabel}"` on MEAN_/SD_/CV_PCT_ of CL/VZ/VSS
 *  (e.g. "Mean CL/F" on an extravascular group) precisely to keep the apparent-
 *  clearance (CL/F) vs true-clearance (CL) distinction in the grouped header —
 *  CL and CL/F differ by ~1/F, so a static "Mean CL" on a CL/F group is a science
 *  mislabel. We strip the `"${statWord} "` prefix back off to recover the bare
 *  parameter name ("CL/F") for the tooltip sentence/gloss. Unset DG tags read
 *  back falsy. `statWord` MUST match the prefix the aggregator wrote (e.g. the
 *  AGG_STATS `word` "%CV" must equal the aggregator's `%CV ` label prefix) — the
 *  exported signature lets a pure unit test pin that contract (GAP-10). */
export function aggDisplayNames(
  statWord: string,
  param: string,
  label: string | null | undefined,
): {bare: string; header: string} {
  const staticFriendly = AGG_PARAM[param]?.friendly ?? param;
  const prefix = `${statWord} `;
  const bare = label && label.startsWith(prefix) ? label.slice(prefix.length) : staticFriendly;
  return {bare, header: `${statWord} ${bare}`};
}

/**
 * Apply presentation to the "By group" aggregated grid: friendly headers and
 * mirrored number formats on the stat columns, plus header tooltips that explain
 * each column (the statistic + the underlying parameter + the bold rule).
 *
 * Idempotent and safe to re-invoke after `grid.dataFrame = newAgg` on re-Compute:
 * the tag/format writes are on the (fresh) DataFrame columns and must be re-run
 * each swap; `grid.onCellTooltip` REPLACES the previous handler, not appends.
 *
 * @param fragileN The caller's fragile-cohort threshold — the SAME constant that
 *   drives this grid's bold-font `onCellPrepare` hook (`FRAGILE_N`, owned by
 *   nca-mapping-pane, which attaches the hook). Passed in rather than imported
 *   because that module imports THIS one; re-declaring the number here would put
 *   the bold rule and the sentence describing it in two places, which is exactly
 *   how MRT / Vss ended up tinted with a tooltip that never said so.
 */
export function configureAggregatedGrid(grid: DG.Grid, fragileN: number): void {
  setGridReadOnly(grid);
  const df = grid.dataFrame;
  // The N / N_* columns render bold below `fragileN` (see the mapping pane's
  // onCellPrepare). Stated on every N header so the emphasis is readable.
  const fragileNote = `Bold when n < ${fmtThreshold(fragileN)} (fragile estimate).`;
  for (const col of df.columns) {
    const parsed = parseAggColumn(col.name);
    if (parsed === null)
      continue;
    const {stat, param} = parsed;
    const {header} = aggDisplayNames(stat.word, param, col.tags['label']);
    col.setTag('friendlyName', header);
    const fmt = stat.format(param);
    if (fmt !== null)
      col.meta.format = fmt;
  }

  grid.onCellTooltip((cell, x, y) => {
    if (!cell.isColHeader)
      return false;
    const name = cell.tableColumn?.name;
    if (name === undefined || name === null)
      return false;

    const parsed = parseAggColumn(name);
    if (parsed !== null) {
      const {stat, param} = parsed;
      const {bare, header} = aggDisplayNames(stat.word, param, df.col(name)?.tags['label']);
      const gloss = AGG_PARAM[param]?.gloss;
      ui.tooltip.show(richTipEl({
        title: header,
        summary: stat.sentence(bare),
        detail: gloss !== undefined ? `${bare} — ${gloss}.` : undefined,
        // Per-parameter contributor counts carry the same bold rule as cohort N.
        note: stat.word === 'N' ? fragileNote : undefined,
      }), x, y);
      return true;
    }

    const group = AGG_GROUP_TOOLTIPS[name];
    if (group !== undefined) {
      ui.tooltip.show(richTipEl({
        title: group.title,
        summary: group.desc,
        note: name === 'N' ? fragileNote : undefined,
      }), x, y);
      return true;
    }
    return false;
  });
}

// ---------------------------------------------------------------------------
// Toxicokinetics aggregate grids (tk-aggregates AD-6) — friendly headers,
// number formats, and header tooltips for the Accumulation / Exposure Margin /
// Sex Comparison panels (freshly-built DataFrames carry only raw column names —
// "F_MEAN_CMAX" etc. — and the default stats tooltip; mirrors
// {@link configureAggregatedGrid}).
// ---------------------------------------------------------------------------

/** Fixed (non-parametric) TK columns → header + tooltip + optional number format
 *  + optional grid width (for headers that truncate at the default width). */
const TK_FIXED: Readonly<Record<string, {header: string; desc: string; format?: string; width?: number}>> = {
  DOSE: {header: 'Dose', desc: 'Dose level of the (Dose × Sex) cell.'},
  RAC_AUC: {header: 'Rac (AUC)', format: '0.000', desc:
    'Accumulation ratio on AUCτ = AUCτ(Day N) / AUCτ(Day 1) — Formula A (CDISC SEND RACAUC).'},
  RAC_CMAX: {header: 'Rac (Cmax)', format: '0.000', desc:
    'Accumulation ratio on Cmax = Cmax(Day N) / Cmax(Day 1).'},
  RAC_FORMULA: {header: 'Rac formula', desc:
    'Formula used: A = AUCτ/AUCτ (default); B = AUCτ,ss/AUCinf,Day1 (opt-in, only where the ' +
    'Day-1 λz gate passes — else labelled fallback to A).'},
  RAC_THEO: {header: 'Rac (theoretical)', format: '0.000', width: 122, desc:
    'Steady-state accumulation index 1/(1−e^(−λz·τ)) from the Day-1 λz. Shown only when λz is ' +
    'present AND adj.R² ≥ 0.90 AND %AUCextrap < 20%; otherwise blank (honest uncertainty).'},
  RAC_SS_FLAG: {header: 'Steady-state', desc:
    'Observed-vs-theoretical reading: ≫ theoretical → time-dependent kinetics; ≪ → not-yet-SS / autoinduction.'},
  RAC_DIVERGENCE: {header: 'Divergence', desc:
    'True when Rac(Cmax) and Rac(AUC) disagree beyond the screening bounds — possible absorption / Tmax shift at SS.'},
  RAT_NOAEL_AUC: {header: 'Rat NOAEL AUCτ', format: '0.0', desc:
    'Rat AUCτ at the NOAEL dose on Day N, per sex (the margin numerator).'},
  MARGIN_EXPOSURE: {header: 'Exposure margin', format: '0.00', desc:
    'Ratified default (ICH M3(R2)): rat NOAEL AUCτ ÷ empirical human AUC. ' +
    'Blank when no measured human AUC is provided.'},
  MARGIN_HED: {header: 'HED margin', format: '0.00', desc:
    'Dose-based (FDA 2005, BSA/Km): HED ÷ clinical dose. No sex term → identical for both sexes.'},
  MARGIN_PROJECTED: {header: 'Projected margin', format: '0.00', desc:
    'CL-projected (allometric by default): rat NOAEL AUCτ ÷ projected human AUC, from the rat CL ' +
    'per kg. Blank when no valid Day-1 λz-derived per-kg CL.'},
  PROJECTION_ALLOMETRIC: {header: 'Allometric', desc:
    'True = projected margin uses allometric CL scaling (default); ' +
    'false = non-allometric (inflated, labelled fallback).'},
  WORST_CASE: {header: 'Worst case', desc:
    'Flags the sex with the smallest (most conservative) projected margin.'},
  RAT_REF_BW_KG: {header: 'Rat ref. BW (kg)', format: '0.000', desc:
    'Rat body weight the allometric projection scales from: the mean body weight of this sex\'s ' +
    'CL-reference profiles when they have one, else the standard 0.25 kg.'},
  RAT_REF_BW_SOURCE: {header: 'Rat ref. BW from', desc:
    'Where the rat reference weight came from: "cohort mean of N profiles" (the CL-reference profiles) ' +
    'or "default" (0.25 kg — no body weight in the data). Blank when no projected margin is formed.'},
};

/** Sex-fold stat prefixes (parametric `<prefix>_<PARAM>` columns). Order matters:
 *  FOLD_DIR_ before FOLD_, and the CV/N prefixes before the bare sex prefixes. */
const TK_FOLD_STATS: ReadonlyArray<{
  prefix: string; label: (p: string) => string; desc: (p: string) => string; format?: string;
}> = [
  {prefix: 'FOLD_DIR_', label: (p) => `${p} direction`,
    desc: (p) => `Which sex is higher for ${p} ("F>M" / "M>F").`},
  {prefix: 'FOLD_', label: (p) => `${p} fold`, format: '0.00',
    desc: (p) => `Sex fold-difference for ${p} = max(M,F) / min(M,F) (≥ 1). No significance claim at composite N.`},
  {prefix: 'F_MEAN_', label: (p) => `${p} mean (F)`,
    desc: (p) => `Female mean ${p} (the fold rests on this).`},
  {prefix: 'M_MEAN_', label: (p) => `${p} mean (M)`,
    desc: (p) => `Male mean ${p} (the fold rests on this).`},
  {prefix: 'F_CV_PCT_', label: (p) => `${p} %CV (F)`, format: '0.0',
    desc: (p) => `Female %CV of ${p}.`},
  {prefix: 'M_CV_PCT_', label: (p) => `${p} %CV (M)`, format: '0.0',
    desc: (p) => `Male %CV of ${p}.`},
  {prefix: 'F_N_', label: (p) => `${p} N (F)`, desc: (p) => `Female N contributing ${p}.`},
  {prefix: 'M_N_', label: (p) => `${p} N (M)`, desc: (p) => `Male N contributing ${p}.`},
];

/** Resolve a TK column's header + tooltip + format (fixed map, then sex-fold
 *  parametric prefixes); null when the column is self-describing (Sex / Dose key). */
function tkColumnDisplay(
  name: string,
): {header: string; desc: string; format?: string; width?: number} | null {
  const fixed = TK_FIXED[name];
  if (fixed !== undefined) return fixed;
  for (const s of TK_FOLD_STATS) {
    if (name.startsWith(s.prefix)) {
      const rawParam = name.slice(s.prefix.length);
      const param = AGG_PARAM[rawParam]?.friendly ?? rawParam;
      return {header: s.label(param), desc: s.desc(param), format: s.format};
    }
  }
  return null;
}

/**
 * Apply presentation to a TK aggregate grid (Accumulation / Exposure Margin /
 * Sex Comparison): friendly headers + number formats + header tooltips. Mirrors
 * {@link configureAggregatedGrid}; idempotent and safe to re-invoke after
 * `grid.dataFrame = newTable` on re-Compute (tags/format are on the fresh frame).
 */
export function configureTkGrid(grid: DG.Grid): void {
  setGridReadOnly(grid);
  const df = grid.dataFrame;
  for (const col of df.columns) {
    const d = tkColumnDisplay(col.name);
    if (d === null) continue;
    col.setTag('friendlyName', d.header);
    if (d.format !== undefined) col.meta.format = d.format;
    // Widen columns whose friendly header truncates at the default width (e.g.
    // "Rac (theoretical)"); mirrors the width handling in `applyFixedGridSpec`.
    if (d.width !== undefined) {
      const gc = grid.columns.byName(col.name);
      if (gc !== null) gc.width = d.width;
    }
  }
  grid.onCellTooltip((cell, x, y) => {
    if (!cell.isColHeader) return false;
    const name = cell.tableColumn?.name;
    if (name === undefined || name === null) return false;
    const d = tkColumnDisplay(name);
    if (d === null) return false;
    ui.tooltip.show(richTipEl({title: d.header, summary: d.desc}), x, y);
    return true;
  });
}

/** Fixed display map for the Bioavailability (F%) grid (per-row-route AD-6).
 *  Columns: `ROUTE`, the stratify column(s) (e.g. `Sex`), `DOSE`, `DN_AUC`,
 *  `IV_DN_AUC`, `F_PCT`. Stratify columns are carried verbatim and
 *  self-describing → no entry here, left untouched. AUC naming mirrors
 *  {@link DISPLAY_HEADERS} (`AUC_0–∞`), with a `dn-` prefix for the
 *  dose-normalized variants. */
const BIOAVAIL_FIXED:
  Readonly<Record<string,
    {header: string; desc: string; format?: string; width?: number}>> = {
      ROUTE: {
        header: 'Route', width: 70,
        desc: 'Extravascular route (PO / SC / IM) of this bioavailability cell.',
      },
      DOSE: {
        header: 'Dose', width: 64,
        desc: 'Dose level of the extravascular stratum, in the mapped dose units.',
      },
      DN_AUC: {
        header: 'dn-AUC_0–∞ (EV)', width: 132,
        desc: 'Dose-normalized AUC_0–∞ of the extravascular stratum: the stratum ' +
        'central-tendency AUC_0–∞ divided by dose.',
        format: '0.000',
      },
      N: {
        header: 'N', width: 44, format: '0',
        desc: 'Extravascular subjects contributing to dn-AUC_0–∞ (EV) — those with a finite AUC_0–∞ ' +
        '(a profile whose λz was excluded by review, or not estimated, does not contribute).',
      },
      IV_DN_AUC: {
        header: 'dn-AUC_0–∞ (IV ref.)', width: 156,
        desc: 'Matched IV-arm dose-normalized AUC_0–∞ — the %F reference ' +
        'denominator for the same stratum. NaN when the dataset has no IV arm ' +
        'for the stratum (never a pooled / cross-stratum fallback).',
        format: '0.000',
      },
      IV_N: {
        header: 'N (IV ref.)', width: 80, format: '0',
        desc: 'IV-arm subjects contributing to the reference dn-AUC_0–∞ of the same stratum.',
      },
      F_PCT: {
        header: '%F', width: 64,
        desc: 'Absolute bioavailability F = 100 × dn-AUC(EV) / dn-AUC(IV ref.). ' +
        'Uncapped — values above 100% are reported, not clamped. NaN when no ' +
        'IV reference exists for the stratum.',
        format: '0.0',
      },
    };

/**
 * Apply presentation to the Bioavailability (F%) grid: friendly headers +
 * number formats + header tooltips. Mirrors {@link configureTkGrid};
 * presentation-only (the F% table from `computeBioavailability` is unchanged).
 * Stratify columns (e.g. `Sex`) are self-describing and left untouched.
 * Idempotent and safe to re-invoke after `grid.dataFrame = newTable` on
 * re-Compute (tags / format live on the fresh frame).
 */
export function configureBioavailabilityGrid(grid: DG.Grid): void {
  // Friendly headers / formats / widths (the long `dn-AUC_0–∞ (…)` labels
  // truncate at the default width) + header tooltips — shared with the dose-
  // proportionality grids via {@link applyFixedGridSpec}.
  applyFixedGridSpec(grid, BIOAVAIL_FIXED);
}

// ── Dose proportionality (Smith power model) grids ─────────────────────────

type FixedGridSpec = Readonly<Record<string,
  {/** Absent ⇒ the column's own `label` tag is the header — for a header that carries a unit only
    *  the producer knows (the sparse `DOSE_LEVEL` column's declared dose unit, BUG-19 S-5). */
   header?: string; desc: string; format?: string; width?: number;
   /** Append the column's own `units` tag to the header (`CL (mL/min/kg)`) — the unit is the
    *  producer's (`CANONICAL_UNIT_TAGS`), never restated here (nca-cl-bw-normalised-units AC-4.4). */
   unitFromTag?: boolean;
   /** Optional colour bands. ONE declaration drives BOTH the DG conditional
    *  ruleset (via {@link dgColorRules}, applied by the grid's `configure*`
    *  function) and the legend rows appended below `desc` in the header tooltip
    *  (via {@link legendRows}, applied by {@link applyFixedGridSpec}).
    *
    *  Previously these two were written separately — the rule as `'<=20'` in the
    *  configure function and the legend as the string `'≤ 20% — …'` here — so
    *  each threshold existed twice with nothing holding the copies together.
    *  They happened to agree; nothing made them. */
   colour?: ColourRule}>>;

/** Per-Route×Sex Smith verdict table (the "Dose Proportionality" tab). */
const DOSE_PROP_FIXED: FixedGridSpec = {
  ROUTE: {header: 'Route', width: 72, desc: 'Route of administration.'},
  SEX: {header: 'Sex', width: 44, desc: 'Sex stratum.'},
  PARAMETER: {header: 'Parameter', width: 84, desc: 'Exposure metric: AUC or Cmax.'},
  BASIS: {header: 'AUC basis', width: 70,
    desc: 'AUC integration basis (AUCinf default, AUClast fallback). Blank for Cmax.'},
  BETA: {header: 'β', width: 52, format: '0.000',
    desc: 'Power-model slope: ln(exposure) = α + β·ln(dose). β = 1 ⇒ dose-proportional.'},
  CI_LOW: {header: '90% CI low', width: 82, format: '0.000',
    desc: 'Lower bound of the 90% confidence interval on β.'},
  CI_HIGH: {header: '90% CI high', width: 86, format: '0.000',
    desc: 'Upper bound of the 90% confidence interval on β.'},
  DOSE_RATIO: {header: 'Dose ratio (r)', width: 96, format: '0.0',
    desc: 'Highest ÷ lowest dose level (max/min).'},
  CR_LOW: {header: 'Accept. low', width: 80, format: '0.000',
    desc: 'Smith acceptance region for β. The 90% CI fully inside ⇒ proportional.'},
  CR_HIGH: {header: 'Accept. high', width: 84, format: '0.000',
    desc: 'Upper bound of the Smith acceptance region for β.'},
  VERDICT: {header: 'Verdict', width: 232,
    desc: 'Proportional / sub- / supra-proportional / inconclusive / inapplicable.'},
  MODE: {header: 'Basis of call', width: 96,
    desc: 'Why: conclusive · df-driven (too few levels / composite) · ' +
      'iiv-driven (underpowered by between-animal variability) · inapplicable.'},
  N: {header: 'N', width: 44, format: '0', desc: 'Subjects in the fit.'},
  DF: {header: 'df', width: 44, format: '0', desc: 'Residual degrees of freedom (N − 2).'},
  N_LEVELS: {header: 'Dose levels', width: 80, format: '0',
    desc: 'Distinct dose levels (≥ 3 required for a verdict).'},
  R_SQUARED: {header: 'R²', width: 52, format: '0.000',
    desc: 'Coefficient of determination of the log-log fit.'},
  MIN_N_PER_GROUP: {header: 'Min n/group', width: 92, format: '0',
    desc: 'Smallest animals/group for a conclusive verdict at this IIV. ' +
      'Blank if already conclusive or infeasible (> 50/group).'},
};

/** Per-subject dose-normalized exposure (the "Dose-Normalized" tab). */
const DOSE_PROP_DN_FIXED: FixedGridSpec = {
  ROUTE: {header: 'Route', width: 72, desc: 'Route of administration.'},
  SEX: {header: 'Sex', width: 44, desc: 'Sex stratum.'},
  DOSE: {header: 'Dose', width: 64, format: '0.###', desc: 'Dose level (per subject), in the mapped units.'},
  DN_AUC: {header: 'AUC/Dose', width: 96, format: '0.000',
    desc: 'Dose-normalized AUC. Flat across doses ⇒ dose-proportional.'},
  DN_CMAX: {header: 'Cmax/Dose', width: 96, format: '0.000',
    desc: 'Dose-normalized Cmax. Flat across doses ⇒ dose-proportional.'},
};

/** Apply friendly headers, formats, widths, colour bands, and header tooltips
 *  from a spec. A spec entry carrying `colour` gets BOTH its DG conditional
 *  ruleset and its tooltip legend from that one declaration, so a band cannot be
 *  painted without being explained. */
function applyFixedGridSpec(
  grid: DG.Grid, spec: FixedGridSpec,
  extra: {
    /** Legend rows for a column banded OUTSIDE the spec (e.g. the sparse grid's adj-R² /
     *  ACCEPTANCE, which read the live `buildColourRules`). */
    readonly headerLegend?: (name: string) => HTMLElement[];
    /** A per-cell rich tooltip (null → DG default) — folded into this grid's ONE handler. */
    readonly cellTip?: (name: string, row: number) => RichTipSpec | null;
    /** Paint the spec's colour bands in the caller's `onCellPrepare` (via `matchBandColour`)
     *  instead of DG column colour-coding. Column coding is drawn OVER `gc.style`, so a band
     *  could not be withheld from one row — the Sparse NCA table needs exactly that for an
     *  excluded composite (manual QA 2026-10-01). Same rules, same background look. */
    readonly bandsInPrepare?: boolean;
  } = {},
): void {
  setGridReadOnly(grid); // covers the F% / dose-proportionality / dn grids
  const df = grid.dataFrame;
  for (const col of df.columns) {
    const d = spec[col.name];
    if (d === undefined) continue;
    col.setTag('friendlyName', fixedHeaderOf(d, col));
    if (d.format !== undefined) col.meta.format = d.format;
    if (d.colour !== undefined && extra.bandsInPrepare !== true)
      col.meta.colors.setConditional(dgColorRules(d.colour));
    if (d.width !== undefined) {
      const gc = grid.columns.byName(col.name);
      if (gc !== null) gc.width = d.width;
    }
  }
  // ONE handler (DG replaces, not appends): header legends + an optional per-cell branch.
  grid.onCellTooltip((cell, x, y) => {
    if (cell.isTableCell && extra.cellTip !== undefined) {
      const nm = cell.tableColumn?.name;
      const i = cell.cell.rowIndex;
      if (nm == null || i < 0) return false;
      const tip = extra.cellTip(nm, i);
      if (tip === null) return false;
      ui.tooltip.show(richTipEl(tip), x, y);
      return true;
    }
    if (!cell.isColHeader) return false;
    const name = cell.tableColumn?.name;
    if (name === undefined || name === null) return false;
    const d = spec[name];
    if (d === undefined) return false;
    // Same `colour` declaration that painted the column above — a band cannot
    // reach the cell without its legend reaching the hover.
    const legend = d.colour !== undefined ? legendRows(d.colour) : (extra.headerLegend?.(name) ?? []);
    ui.tooltip.show(richTipEl({title: fixedHeaderOf(d, cell.tableColumn!), summary: d.desc, extra: legend}), x, y);
    return true;
  });
}

/** A fixed-spec column's header: the spec's, else the column's `label` tag, else its name. */
function fixedHeaderOf(d: FixedGridSpec[string], col: DG.Column): string {
  const base = d.header ?? col.getTag('label') ?? col.name;
  return d.unitFromTag === true ? headerWithUnit(base, col) : base;
}

/** Hide a grid's `ROUTE` column unless it varies (BUG-19 AD-4; the `CL_BASIS` compact rule): on a
 *  single-route analysis — most sparse studies — it says nothing. Re-applied on every configure,
 *  so a recompute's frame swap re-evaluates it. */
function hideRouteUnlessVarying(grid: DG.Grid): void {
  const col = grid.dataFrame?.col('ROUTE') ?? null;
  const gc = grid.columns.byName('ROUTE');
  if (col !== null && gc !== null) gc.visible = col.categories.length > 1;
}

/** Friendly headers + tooltips for the "Dose Proportionality" verdict grid. */
export function configureDoseProportionalityGrid(grid: DG.Grid): void {
  applyFixedGridSpec(grid, DOSE_PROP_FIXED);
}

/** Friendly headers + tooltips for the "Dose-Normalized" exposure grid. */
export function configureDoseProportionalityDnGrid(grid: DG.Grid): void {
  applyFixedGridSpec(grid, DOSE_PROP_DN_FIXED);
}

// ── UC-04 sparse panels — friendly headers + tooltips ──────────────────────
// Units are the canonical basis `computeSparseStrata` converts every input to (BUG-19): conc →
// mg/L, dose → mg or mg/kg, time → h ⇒ AUC mg·h/L, Cmax mg/L, t½/MRT/Tmax h; CL / Vz / Vss and
// their companions carry their unit as a column tag (headers read it — `unitFromTag`). The dose
// LEVEL stays in the declared unit (its header comes from the column's `label` tag).

/** Per-dose closed-form + derived parameters (the "Sparse NCA" tab). */
const SPARSE_NCA_FIXED: FixedGridSpec = {
  ROUTE: {width: 72, header: 'Route',
    desc: 'Route of the composite. An analysis with more than one route has one composite per route × dose.'},
  DOSE_LEVEL: {width: 88, format: '0.###', desc: 'Nominal dose level, in the unit of the dose column.'},
  // Written by the view only when the pipeline run raised a reason of its own (F10).
  PIPELINE: {header: 'Pipeline', width: 180,
    desc: 'Reasons the pipeline run sent this composite to review that the λz acceptance check does ' +
      'not raise — hover a cell for every reason and the run.'},
  TOPOLOGY: {header: 'Topology', width: 84,
    desc: 'Sampling design derived from the animal × nominal-time table: ' +
      'destructive (1 sample/animal) / batch (partial overlap) / serial (full).'},
  AUC_LAST: {header: 'AUClast', width: 80, format: '0.000',
    desc: 'Closed-form composite AUClast (mg·h/L), linear trapezoid from (0, 0) — the regulatory headline. ' +
      'The linear rule can overestimate a convex decline across wide sampling gaps (see Notes).'},
  AUC_SE: {header: 'SE', width: 64, format: '0.000',
    desc: 'Holder/Bailer standard error of AUClast (mg·h/L).'},
  // The Satterthwaite df floor is the sparse path's OWN threshold — it is not
  // the λz `minPoints` rule and does not track Rules (they collide at 3 only by
  // coincidence). Upper bound 100 mirrors the previous explicit ruleset.
  AUC_DF: {header: 'df', width: 48, format: '0.0',
    desc: 'Nedelman-Jia Satterthwaite degrees of freedom for the t-CI.',
    colour: {kind: 'conditional', bands: [
      {band: {op: 'ge', v: 3, upTo: 100}, color: COLOR_OK, gloss: 'adequate'},
      {band: {op: 'lt', v: 3}, color: COLOR_WARNING,
        gloss: 'very few effective observations; wide/unstable CI'},
    ]}},
  AUC_CI_LOW: {header: 'CI low', width: 70, format: '0.000',
    desc: 'Lower 95% CI on AUClast (Student-t, mg·h/L).'},
  AUC_CI_HIGH: {header: 'CI high', width: 72, format: '0.000',
    desc: 'Upper 95% CI on AUClast (mg·h/L).'},
  AUC_CI_METHOD: {header: 'CI method', width: 78,
    desc: 'Closed-form variance used: Bailer (destructive, r_ij = 0) / Holder (batch/serial).'},
  CMAX: {header: 'Cmax', width: 70, format: '0.000', desc: 'Composite peak mean concentration (mg/L).'},
  TMAX: {header: 'Tmax', width: 56, format: '0.##', desc: 'Time of the composite Cmax (h).'},
  AUC_INF: {header: 'AUCinf', width: 76, format: '0.000',
    desc: 'Composite AUC extrapolated to infinity (mg·h/L): the linear AUClast of this row + Clast / λz. ' +
      '%AUCextrap, CL, Vz, MRT and Vss rest on it.'},
  HALF_LIFE: {header: 't½', width: 56, format: '0.000',
    desc: 'Terminal half-life from the composite λz (h).'},
  // Per kg on every composite (nca-cl-bw-normalised-units D-d). The estimator is named (plan-review
  // F-4): a composite value is a ratio of means, not the mean of per-animal values a dense table shows.
  CL: {header: 'CL', unitFromTag: true, width: 104, format: '0.000',
    desc: 'Composite estimate — dose per kg ÷ AUCinf of the mean concentration curve (not a mean of ' +
      'animal values). Assumes every animal received the same dose per kg.'},
  VZ: {header: 'Vz', unitFromTag: true, width: 76, format: '0.000',
    desc: 'Composite estimate of the terminal volume of distribution per kg, from the mean ' +
      'concentration curve (not a mean of animal values).'},
  MRT: {header: 'MRT', width: 64, format: '0.000', desc: 'Mean residence time (h).'},
  VSS: {header: 'Vss', unitFromTag: true, width: 76, format: '0.000',
    desc: 'Composite estimate of the steady-state volume of distribution per kg (IV only), from the ' +
      'mean concentration curve (not a mean of animal values).'},
  CL_ABS: {header: 'CL abs.', unitFromTag: true, width: 88, format: '0.000',
    desc: 'Absolute composite clearance: the per-kg value × the group-mean body weight (a mg/kg ' +
      'dose), or dose / AUCinf (an absolute mg dose). Blank without a body weight on a mg/kg dose.'},
  VZ_ABS: {header: 'Vz abs.', unitFromTag: true, width: 76, format: '0.000',
    desc: 'Absolute composite terminal volume of distribution.'},
  VSS_ABS: {header: 'Vss abs.', unitFromTag: true, width: 80, format: '0.000',
    desc: 'Absolute composite steady-state volume of distribution (IV only).'},
  BW_KG: {header: 'Group-mean BW', unitFromTag: true, width: 140, format: '0.000',
    desc: 'Mean body weight of the animals pooled into the composite (CL abs. uses it).'},
  // Mirrors the dense grid's PCT_EXTRAP bands, but with the sparse path's own
  // FIXED thresholds — the composite NCA deliberately does not follow the Rules
  // inputs (sparse-nca.ts: it accepts the best-fit λz at `minRSquared: 0` and
  // flags rather than hides). See TD-06 in TODO.md.
  PCT_AUCEXTRAP: {header: '%AUCextrap', width: 90, format: '0.0',
    desc: 'Percent of AUCinf from λz extrapolation.',
    colour: {kind: 'conditional', bands: [
      {band: {op: 'le', v: 20}, color: COLOR_OK, unit: '%',
        gloss: 'acceptable extrapolation'},
      {band: {op: 'range', lo: 20, hi: 50}, color: COLOR_WARNING, unit: '%',
        gloss: 'borderline; review'},
      {band: {op: 'gt', v: 50}, color: COLOR_CRITICAL, unit: '%',
        gloss: 'AUCinf considered unreliable'},
    ]}},
  LAMBDAZ_R2: {header: 'λz R²', width: 60, format: '0.000',
    desc: 'R² of the composite terminal-slope (λz) fit.',
    colour: {kind: 'conditional', bands: [
      {band: {op: 'gt', v: 0.9}, color: COLOR_OK, gloss: 'good fit'},
      {band: {op: 'range', lo: 0.85, hi: 0.9}, color: COLOR_WARNING,
        gloss: 'borderline'},
      {band: {op: 'lt', v: 0.85}, color: COLOR_CRITICAL,
        gloss: 'poor fit; λz/CL/t½ are fragile'},
    ]}},
  LAMBDAZ_ADJ_R2: {header: 'adj. R² (λz)', width: 80, format: '0.000',
    desc: 'Adjusted R² of the composite terminal-slope (λz) fit — the value the λz acceptance ' +
      'criterion reads.'},
  ACCEPTANCE: {header: 'λz acceptance', width: 100,
    desc: 'λz acceptance verdict of the composite fit — the same criteria as a dense profile ' +
      '(adj. R², span ratio, N points, %AUCextrap). A pooled-cohort fit: the criteria were ' +
      'calibrated on individual serial profiles.'},
  // 170 px: three pills (`span < 2` `N ≤ 3` `%extrap > 20%`) — the most measured on the fixtures.
  [ACCEPT_REASONS_COL]: {header: 'Violated', width: 170,
    desc: 'The λz acceptance criteria the composite fit violates, one pill each (span ratio, adj. R², ' +
      'N points, %AUCextrap). Hover a cell for the full criteria.'},
  FRAGILITY: {header: 'Fragility', width: 170,
    desc: 'Sampling-design uncertainty flags: n<4 · df<3 · variance-modeled · bootstrap-suppressed. ' +
      'The λz fit quality is judged by the λz acceptance column.'},
  AUC_BOOT_CI: {header: 'AUC bootstrap CI', width: 210,
    desc: 'Exploratory BCa bootstrap CI — NOT for regulatory use. The closed-form ' +
      'CI is the headline; this is shown for completeness only.'},
  NOTES: {header: 'Notes', width: 220,
    desc: 'Per-group advisories (composite BW, exploratory bootstrap, topology mismatch, …).'},
};

/** Per (dose × nominal time) composite profile (the "Composite Profile" tab). */
const SPARSE_COMPOSITE_FIXED: FixedGridSpec = {
  ROUTE: {width: 72, header: 'Route', desc: 'Route of the composite.'},
  DOSE_LEVEL: {width: 88, format: '0.###', desc: 'Dose group, in the unit of the dose column.'},
  NOMINAL_TIME: {header: 'Nominal time (h)', width: 110, format: '0.##', desc: 'Scheduled sampling time.'},
  N: {header: 'N', width: 44, format: '0', desc: 'Animals sampled at this nominal time.'},
  N_BLQ: {header: 'BLQ', width: 48, format: '0', desc: 'Below-LLOQ samples at this time.'},
  MEAN_CONC: {header: 'Mean conc (mg/L)', width: 116, format: '0.000',
    desc: 'Arithmetic mean concentration (BLQ imputed before averaging).'},
  SD: {header: 'SD', width: 70, format: '0.000', desc: 'Sample standard deviation (mg/L).'},
  CV_PCT: {header: '%CV', width: 64, format: '0.0', desc: 'Coefficient of variation.'},
  PCT_BLQ: {header: '%BLQ', width: 64, format: '0.0',
    desc: 'Percent of samples below LLOQ at this nominal time.',
    colour: {kind: 'conditional', bands: [
      {band: {op: 'le', v: 50}, color: COLOR_OK, unit: '%',
        gloss: 'mean rests on measured data'},
      {band: {op: 'gt', v: 50}, color: COLOR_WARNING, unit: '%',
        gloss: 'mean zeroed for the AUClast endpoint'},
    ]}},
};

/**
 * Friendly headers + tooltips + reliability colour bands for the "Sparse NCA"
 * grid. Only the *diagnostic* columns are tinted (the same convention as the
 * main results grid, which colours PCT_EXTRAP / R² / npoints — never the
 * exposure values themselves). The bands encode the thresholds stated in the
 * header tooltips (thesis #4): AUC-extrapolation reliability, λz fit quality,
 * and the Satterthwaite df.
 */
/** One sparse composite's pooled-cohort shape (AC-S4.8.0), aligned with the parameter-table rows. */
export interface SparseCohort {
  readonly nAnimals: number;
  readonly minPerTimepoint: number;
  readonly maxPerTimepoint: number;
  /** The composite's violated acceptance criteria (short form), empty unless Flagged. */
  readonly reasons: ReadonlyArray<string>;
}

/** The pooled-cohort cue that leads every composite verdict surface (AC-S4.8.0; R1/LBC-3): the
 *  verdict is visibly a POOLED fit, never presented as if it were one animal's serial profile.
 *  `min === max` collapses the range to "({n} per timepoint)" (R2/C). */
export function pooledCohortCue(c: SparseCohort): string {
  const per = c.minPerTimepoint === c.maxPerTimepoint ?
    `${c.minPerTimepoint} per timepoint` : `${c.minPerTimepoint}–${c.maxPerTimepoint} per timepoint`;
  return `Pooled-cohort fit — composite of ${c.nAnimals} animals (${per}); criteria calibrated ` +
    'on individual serial profiles.';
}

/**
 * Friendly headers + tooltips + bands for the "Sparse NCA" grid. The spec-declared tints
 * (PCT_AUCEXTRAP / LAMBDAZ_R2 / AUC_DF) come from SPARSE_NCA_FIXED as before. The composite
 * `adj. R² (λz)` and `λz acceptance` columns (nca-workflow-stages F-S4.8 / U2) are banded from
 * the SAME `buildColourRules` declarations as the per-subject grid — two tiers at the live `R_acc`
 * (the value the user entered, applied verbatim), the verdict triad for ACCEPTANCE — painted as TEXT colour like
 * the results grid, with no second band definition. Each ACCEPTANCE cell's tooltip leads with the
 * pooled-cohort cue.
 */
export function configureSparseNcaGrid(
  grid: DG.Grid, thresholds: ResultsColorThresholds = DEFAULT_RESULTS_THRESHOLDS,
  cohorts: ReadonlyArray<SparseCohort> = [],
  /** The Exclude λz decision of the composite on row `i`, or null: its λz-dependent cells then read as
   *  not reportable, exactly as a dense excluded profile's. Read at render time, so a decision
   *  repaints on `invalidate()`. */
  excludedAt: (row: number) => ReviewDecisionRecord | null = () => null,
): void {
  const rules = buildColourRules(thresholds);
  const banded: Record<string, ColourRule | undefined> = {
    LAMBDAZ_ADJ_R2: rules.LAMBDAZ_ADJ_R2, [ACCEPTANCE_COL]: rules[ACCEPTANCE_COL],
  };
  const df = grid.dataFrame;
  // The `Pipeline` cell's hidden tooltip detail is a service column, never a visible one.
  const pipelineInfo = grid.columns.byName(LZ_PIPELINE_INFO_COL);
  if (pipelineInfo !== null) pipelineInfo.visible = false;
  applyFixedGridSpec(grid, SPARSE_NCA_FIXED, {
    bandsInPrepare: true,
    headerLegend: (name) => {const r = banded[name]; return r !== undefined ? legendRows(r) : [];},
    cellTip: (name, i) => {
      if (name === PIPELINE_COL) return pipelineTipAt(df, i);
      if (isSparseLzDependentCol(name)) {
        const d = excludedAt(i);
        return d !== null ? excludedCellTipSpec(d, true) : null;
      }
      // The `Violated` pills share the verdict's tip (the full criteria); an empty cell has none.
      const onPills = name === ACCEPT_REASONS_COL && (df?.col(ACCEPT_REASONS_COL)?.getString(i) ?? '') !== '';
      if ((name !== ACCEPTANCE_COL && !onPills) || i >= cohorts.length) return null;
      const c = cohorts[i];
      const verdict = df?.col(ACCEPTANCE_COL)?.getString(i) ?? '';
      return {
        title: verdict === '' ? 'λz acceptance' : verdict,
        summary: pooledCohortCue(c),
        ...(c.reasons.length > 0 ? {detail: c.reasons.join(' · ')} : {}),
      };
    },
  });
  grid.onCellPrepare((gc) => {
    if (!gc.isTableCell) return;
    const nm = gc.gridColumn?.name;
    if (nm === ACCEPT_REASONS_COL) {
      reasonPillFont(gc);
      return;
    }
    // An EXCLUDED composite's λz-dependent value: kept, muted, marked — the dense F-S4.5 rendering.
    if (nm != null && isSparseLzDependentCol(nm) && typeof gc.cell.value === 'number' &&
        Number.isFinite(gc.cell.value) && excludedAt(gc.cell.rowIndex) !== null) {
      const col = df?.col(nm) ?? null;
      // Returns BEFORE the spec's background band below — a band hue never sits under a value
      // that is not reported (§1 "Review column + excluded values").
      paintExcludedCell(gc, col !== null ? col.getString(gc.cell.rowIndex) : String(gc.cell.value));
      return;
    }
    // The spec-declared bands (%AUCextrap / λz R² / df) as BACKGROUND fills, painted here rather
    // than by column colour-coding (`bandsInPrepare`) so the branch above can withhold them.
    const specRule = nm == null ? undefined : SPARSE_NCA_FIXED[nm]?.colour;
    if (specRule !== undefined) {
      const fill = matchBandColour(specRule, gc.cell.value);
      if (fill !== null) gc.style.backColor = fill;
    }
    const r = nm == null ? undefined : banded[nm];
    if (r === undefined || r.kind === 'crossColumn') return;
    const fill = matchBandColour(r, gc.cell.value);
    if (fill !== null) gc.style.textColor = fillToTextColour(fill);
    if (nm === ACCEPTANCE_COL) boldGridCell(gc);
  });
  applyReasonPillColours(grid, thresholds);
  hideRouteUnlessVarying(grid);
}

/**
 * Friendly headers + tooltips + the BLQ% reliability band for the "Composite
 * Profile" grid. A timepoint with >50% BLQ has its mean zeroed for the AUClast
 * endpoint, so it is flagged.
 */
export function configureSparseCompositeGrid(grid: DG.Grid): void {
  // PCT_BLQ's band is declared in SPARSE_COMPOSITE_FIXED (rule + legend, one source).
  applyFixedGridSpec(grid, SPARSE_COMPOSITE_FIXED);
  hideRouteUnlessVarying(grid);
}

/** The secondary grids that carry a hand-maintained presentation map. */
export type ConfigurableGrid =
  | 'tk' | 'bioavailability' | 'doseProportionality'
  | 'doseProportionalityDn' | 'sparseNca' | 'sparseComposite';

/** The five spec-driven grids, in ONE registry. Read by both completeness
 *  guards — {@link unmappedGridColumns} (is every column mapped?) and
 *  {@link colourBandProjections} (does every band explain itself?) — so neither
 *  carries its own hand-picked list of which grids to check. `tk` is absent: it
 *  resolves through {@link tkColumnDisplay}'s prefix logic, not a flat spec. */
const FIXED_GRID_SPECS: Readonly<Record<
  Exclude<ConfigurableGrid, 'tk'>, FixedGridSpec>> = {
    bioavailability: BIOAVAIL_FIXED,
    doseProportionality: DOSE_PROP_FIXED,
    doseProportionalityDn: DOSE_PROP_DN_FIXED,
    sparseNca: SPARSE_NCA_FIXED,
    sparseComposite: SPARSE_COMPOSITE_FIXED,
  };

/**
 * Columns of `df` that this grid's `configure*Grid` would leave with a raw
 * header and DG's generic min/avg/max/std tooltip.
 *
 * The per-subject results grid gets {@link missingPresentationEntries}, which can
 * diff against a static column list. These grids have no such list — their
 * columns depend on the caller's stratify/covariate choices — so completeness is
 * checked against a REAL computed table instead: run the adapter's own
 * `compute*` on a fixture, hand the result here, assert the return is empty.
 * That is strictly stronger evidence than reading the source (rule 16), and it
 * needs no change to the adapter modules.
 *
 * Each case delegates to the SAME lookup the grid itself uses (`tkColumnDisplay`,
 * the `FixedGridSpec` maps), so the check cannot drift from real behaviour.
 *
 * @param df            A table produced by the corresponding adapter `compute*`.
 * @param grid          Which grid's presentation map to check against.
 * @param selfDescribing Caller-supplied columns that intentionally carry no
 *                      entry — the stratify / covariate columns passed through
 *                      verbatim (e.g. `Sex`, `Dose_mg_kg`), which are named by
 *                      the user and self-explanatory. Anything NOT listed here
 *                      is expected to be mapped.
 */
export function unmappedGridColumns(
  df: DG.DataFrame,
  grid: ConfigurableGrid,
  selfDescribing: ReadonlyArray<string> = [],
): ReadonlyArray<string> {
  const resolves = grid === 'tk' ?
    (n: string): boolean => tkColumnDisplay(n) !== null :
    (n: string): boolean => n in FIXED_GRID_SPECS[grid];
  const exempt = new Set<string>(selfDescribing);
  return df.columns.names().filter((n) => !exempt.has(n) && !resolves(n));
}

/**
 * Open the results DataFrame as a new top-level TableView with full
 * presentation applied (colour rules + grid-side tweaks).
 */
export function openResultsView(df: DG.DataFrame): DG.TableView {
  configureResultsDataFrame(df);
  const view = grok.shell.addTableView(df);
  configureResultsGrid(view.grid);
  return view;
}

// ──────────────────────────────────────────────────────────────────────────
// Header tooltips
// ──────────────────────────────────────────────────────────────────────────

/** ARGB integer (`0xAARRGGBB`) → CSS hex string `#RRGGBB`. Alpha discarded;
 *  swatches in the tooltip are opaque. */
function argbToCssHex(argb: number): string {
  return '#' + (argb & 0x00FFFFFF).toString(16).padStart(6, '0').toUpperCase();
}

/** Small filled square that previews a text colour in the legend. */
function colorSwatch(argb: number): HTMLElement {
  const sw = ui.div();
  sw.style.cssText =
    'display: inline-block;' +
    'width: 10px; height: 10px;' +
    'border-radius: 2px;' +
    `background: ${argbToCssHex(argb)};` +
    'margin-right: 6px;' +
    'vertical-align: middle;' +
    'flex-shrink: 0;';
  return sw;
}

/** One row in the legend: [swatch] [description]. The swatch carries the
 *  colour; the description stays in default tooltip text colour. */
function legendRow(argb: number, text: string): HTMLElement {
  const row = ui.divH([colorSwatch(argb), ui.divText(text)]);
  row.style.cssText = 'align-items: center; margin-top: 2px;';
  return row;
}

/** Spacer between description and legend rows — keeps the tooltip readable
 *  without introducing a section title. */
function legendSpacer(): HTMLElement {
  const s = ui.div();
  s.style.cssText = 'height: 6px;';
  return s;
}

/** Build a plain header tooltip spec: title + one-paragraph description, no
 *  colour legend. Use for columns that have no anomaly highlighting. */
function plainTooltip(title: string, description: string): RichTipSpec {
  return {title, summary: description};
}

/** DataFrame column name → header-tooltip SPEC. Every result column has a
 *  header tooltip: whichever columns {@link buildColourRules} returns carry a
 *  colour legend built from that same declaration — the one that registers the
 *  paint — and the rest get a plain title + description.
 *
 *  **The coloured set is deliberately NOT enumerated here.** This comment used
 *  to name "the seven coloured columns" and went stale the moment `SPAN_RATIO`
 *  was banded — the fourth restatement of that set to be wrong, in the very
 *  file whose column profile derives it precisely so it cannot be (AD-2).
 *  Function (not const)
 *  because the legends embed live Rules thresholds: a new map is needed each
 *  time the user edits a Rules input and Compute is re-run.
 *
 *  Returns SPECS, not rendered elements or lazy factories, so the completeness
 *  guard can INSPECT them — `missingPresentationEntries` asserts that every
 *  coloured column's spec actually carries legend rows. An opaque
 *  `() => HTMLElement` could not be checked, which is why MRT / Vss went tinted
 *  but unexplained for a full release cycle.
 *
 *  The map MUST cover every column in {@link PRESENTATION_COLUMN_ORDER}: a
 *  missing key makes `onCellTooltip` return `false` and DG falls through to its
 *  generic min/avg/max/std stats tooltip, which names no parameter and explains
 *  nothing. That is how `C0_EXTRAP` shipped unexplained — it was added to the
 *  order + the writer schema but to none of the presentation maps. Enforced by
 *  {@link missingPresentationEntries} (pinned in results-view-presentation-tests). */
export function buildHeaderTooltips(
  thresholds: ResultsColorThresholds,
): Readonly<Record<string, RichTipSpec>> {
  // TEXT only. Legends are attached below by KEY, not written here — a column
  // cannot be given the wrong column's legend, and cannot be silently missed.
  const text: Record<string, RichTipSpec> = {
    STATUS: {
      title: 'Computation',
      summary:
        'Technical computation status for the subject\'s profile — whether the ' +
        'algorithm produced results (ok / partial / failed). Distinct from ' +
        'Acceptance, which asks whether the λz fit can be trusted.',
    },

    // Banded (categorical) → its legend is attached by the loop below, from the
    // same declaration that paints the cell.
    ACCEPTANCE: {
      title: 'Acceptance',
      summary:
        'λz acceptance verdict (SAP criteria). Flagged if span / adj. R² / ' +
        'N points / %AUCextrap is out of bounds; Not estimated when λz could not be fitted.',
    },

    [ACCEPT_REASONS_COL]: {
      title: 'Violated',
      summary: 'The λz acceptance criteria a Flagged profile violates, one pill each.',
      detail: 'Hover a cell for the full criteria.',
    },

    REVIEW: {
      title: 'Review',
      summary:
        'Whether the analyst still has to look at this profile, and what they decided: To review · ' +
        '✓ Accepted* (accepted as is, with a reason) · Acknowledged · Excluded (λz-dependent values ' +
        'not reportable) · Not needed (Accepted fit) · Not applicable (reviewed at the group composite).',
      detail: 'Hover a cell for the reasons, who and when.',
    },

    PIPELINE: {
      title: 'Pipeline',
      summary:
        'Reasons the pipeline run sent this profile to review that the λz acceptance check does not ' +
        'raise — a CRO incident comment, a screening heuristic (exposure outlier, terminal rise) or a ' +
        'unit-error heuristic. Shown only when the run raised one.',
      detail: 'Hover a cell for every reason and the run that raised it.',
    },

    PCT_EXTRAP: {
      title: '%AUCextrap',
      summary: 'Share of AUC₀–∞ obtained by extrapolation past the last measurable concentration.',
    },

    LAMBDAZ_ADJ_R2: {
      title: 'adj. R² (λz)',
      summary:
        'Adjusted R² of the log-linear regression used to estimate the ' +
        'terminal elimination rate constant λz.',
      detail:
        `Acceptance: a fit with adj. R² < ${fmtThreshold(thresholds.lambdaAcceptanceRSquared)} ` +
        '(the adj-R² acceptance threshold) is flagged. The threshold is applied as ' +
        'entered — it is a separate setting from the λz min adj-R² point-selection gate.',
    },

    LAMBDAZ_NPOINTS: {
      title: 'N points (λz)',
      summary: 'Number of terminal-phase points included in the λz regression.',
    },

    SPAN_RATIO: {
      title: 'Span ratio (λz)',
      // The last clause is the whole reason this column sits BESIDE adj. R² rather
      // than instead of it, and it is the minimum that survived the tooltip trim
      // (user-directed, restored on review 2026-08-13). Our own corpus carries the
      // case it describes: indometh subject 1 pairs a span of 0.685 half-lives
      // with an adj. R² of 0.994 — a near-perfect line fitted over almost no
      // observed elimination. Without it a scientist has two λz diagnostics side
      // by side and no stated reason for the second.
      summary:
        'How many terminal half-lives the λz regression window actually spans: ' +
        '(t_last − t_first) / t_½. Diagnostic only — it never changes the fit or ' +
        'any derived parameter. A high adj. R² does not rule out a short span ' +
        '(Noe 2020).',
      // The `detail` paragraph was REMOVED (user-directed 2026-08-13, tooltip too
      // tall). It explained the mechanism — a short span leaves AUC_0–∞, t_½, CL,
      // Vz, MRT and Vss resting on very little observed elimination, which a high
      // adj. R² cannot rule out — and attributed the ≥ 4 figure to Noe (2020).
      //
      // What survives without it: the three legend rows below still carry both
      // thresholds AND the conditional form of the 4-half-life bar ("under high
      // measurement variability"), so no number is stated more strongly than the
      // source supports. What is GONE from the product: the adj-R²-cannot-see-this
      // argument, i.e. the reason this column sits beside adj. R² at all, and the
      // Noe attribution — recorded in TODO.md rather than silently dropped.
      //
      // Do NOT move the citation into a legend gloss to "keep" it: a year there
      // trips the legend-vs-rule number correspondence guard in
      // results-view-presentation-tests, which is why it lived up here.
    },

    PCT_AUMC_EXTRAP: {
      title: '%AUMCextrap',
      summary:
        'Share of AUMC₀–∞ obtained by extrapolation past the last measurable ' +
        'point. Always ≥ %AUCextrap (the moment tail is time-weighted) — a ' +
        'non-regulatory fragility flag for MRT and Vss, which are highlighted ' +
        'when this is high.',
    },

    // Plain tooltips — no anomaly colouring, just title + description.
    SUBJECT: plainTooltip('Subject', 'Subject identifier.'),
    ROUTE: plainTooltip('Route', 'Route of administration.'),
    CL_BASIS: plainTooltip(
      'CL basis',
      'Per-subject clearance basis: "CL" (true clearance, intravascular, F = 1) ' +
      'or "CL/F" (apparent clearance, extravascular, includes bioavailability F). ' +
      'Follows this subject\'s route; also governs Vz vs Vz/F.',
    ),
    DOSE_NORM: plainTooltip(
      'Dose norm.',
      'How the dose was passed to the computation: "absolute" (mg — an absolute dose, or a ' +
      'mg/kg dose × body weight) or "per_kg" (a mg/kg dose with no body weight). CL, Vz and Vss ' +
      'are reported per kg either way.',
    ),
    // The four per-criterion λz acceptance flags (nca-grid-status-spec §3). Plain
    // tooltips: they are NOT colour-banded (ACCEPTANCE carries the tint), so no
    // legend — the glyph vocabulary and the threshold are stated in prose instead.
    // Thresholds are the SAME live Rules the colour bands read.
    SPAN_FLAG: plainTooltip(
      'Flag: span',
      'λz acceptance flag — span ratio < 2 (the terminal window covers under two ' +
      'half-lives). ✗ violated · ✓ met · n/a when λz was not estimated.',
    ),
    RSQ_FLAG: plainTooltip(
      'Flag: R²',
      `λz acceptance flag — adj. R² < ${fmtThreshold(thresholds.lambdaAcceptanceRSquared)} ` +
      '(the adj-R² acceptance threshold, applied as entered — a separate setting from the ' +
      'λz min adj-R² point-selection gate). ✗ violated · ✓ met · n/a when λz was not estimated.',
    ),
    NEXT_FLAG: plainTooltip(
      'Flag: N',
      'λz acceptance flag — N points ≤ the Rules minimum (default 3; a fit AT the ' +
      'floor is fragile). ✗ violated · ✓ met · n/a when λz was not estimated.',
    ),
    EXTRAP_FLAG: plainTooltip(
      'Flag: %extrap',
      'λz acceptance flag — %AUCextrap above the Rules warn threshold (default ' +
      '20%). ✗ violated · ✓ met · n/a when λz was not estimated.',
    ),
    CMAX: plainTooltip('C_max', 'Maximum observed plasma concentration.'),
    // Sized to its siblings (C_max is one line). The method chain, the library
    // it ports, the AUC/AUMC insertion and the λz peak-exclusion are all TRUE
    // and all deliberately absent — a header hover carries what changes a
    // reading, not the reference text (that lives in the writer/controller
    // comments). What survives is the
    // one fact a scientist cannot recover from the grid: the number may not be
    // a measurement, and the row does not say. The `0.000` clause stays because
    // that fallback is the case that LOOKS like real data. User-trimmed twice
    // (2026-07-29, BUG-12) — read the TODO entry before growing this again.
    C0_EXTRAP: {
      title: 'C_0',
      summary: 'Concentration at the dose time. IV bolus only.',
      detail:
        'Measured when a t = 0 sample exists, otherwise estimated — or 0.000 ' +
        'when there is nothing to estimate from. Which one applies is recorded ' +
        'per row in the provenance column.',
    },
    TMAX: plainTooltip('t_max', 'Time of the maximum observed concentration.'),
    AUC_LAST: plainTooltip(
      'AUC_0–t',
      'AUC from time 0 to the last measurable concentration. ' +
    'Computed via the linear-up / log-down trapezoidal method.',
    ),
    AUC_INF: plainTooltip(
      'AUC_0–∞',
      'AUC extrapolated to infinity: AUC_0–t + C_last / λz.',
    ),
    HALF_LIFE: plainTooltip('t_½', 'Terminal half-life: ln(2) / λz.'),
    LAMBDA_Z: plainTooltip('λz', 'Terminal elimination rate constant.'),
    // The basis sentence is SELF-CONTAINED, not a pointer. It used to read "see
    // the row columns: (1) CL basis …", and when the compact profile briefly moved
    // `CL_BASIS` to Full that instruction pointed at a column the default view was
    // hiding. `CL_BASIS` is now compact-visible wherever it VARIES
    // (COMPACT_WHEN_VARYING), so the pointer would be correct again on exactly the
    // frames that need it — but stating the rule costs nothing and survives any
    // future profile change, which the pointer demonstrably did not.
    CL: plainTooltip(
      'CL (mL/min/kg)',
      'Clearance per kg of body weight: Dose / AUC_0–∞ ÷ body weight. Basis: intravascular ' +
    'rows are true CL (F = 1); extravascular rows are apparent CL/F, which ' +
    'includes bioavailability — the two differ by ~1/F. This header stays neutral ' +
    'when a study mixes routes; the per-row basis is the CL basis column, shown ' +
    'whenever it varies. Blank when an absolute dose has no body weight to divide by ' +
    '(CL abs. is then reported).',
    ),
    VZ: plainTooltip(
      'Vz (L/kg)',
      'Volume of distribution during the terminal phase, per kg of body weight. Basis: ' +
    'intravascular rows are true Vz; extravascular rows are apparent Vz/F, which includes ' +
    'bioavailability. This header stays neutral when a study mixes routes; the per-row ' +
    'basis is the CL basis column, shown whenever it varies.',
    ),
    // Absolute companions + the weight (nca-cl-bw-normalised-units D-a/c) — Full profile only.
    CL_ABS: plainTooltip(
      'CL abs. (L/h)',
      'Absolute clearance (Dose / AUC_0–∞), not divided by body weight — the same CL vs CL/F ' +
    'basis as CL. Blank when a mg/kg dose has no body weight to multiply it by.',
    ),
    VZ_ABS: plainTooltip(
      'Vz abs. (L)',
      'Absolute terminal volume of distribution, not divided by body weight. Blank when a ' +
    'mg/kg dose has no body weight.',
    ),
    VSS_ABS: plainTooltip(
      'Vss abs. (L)',
      'Absolute steady-state volume of distribution (IV only), not divided by body weight. ' +
    'Blank when a mg/kg dose has no body weight.',
    ),
    BW_KG: plainTooltip(
      'BW (kg)',
      'The body weight CL, Vz and Vss are divided by (the profile\'s first row, in kg). ' +
    'Blank when no body weight is mapped.',
    ),
    AUMC_LAST: plainTooltip(
      'AUMC_0–t',
      'Area under the first-moment curve (∫ t·C dt) to the last measurable ' +
    'concentration. Reported even when λz did not fit (AUMC_0–∞ is then blank).',
    ),
    AUMC_INF: plainTooltip(
      'AUMC_0–∞',
      'First-moment area to infinity: AUMC_0–t + the two-term tail ' +
    '(t_last·C_last/λz + C_last/λz²). The basis for MRT and Vss.',
    ),
    // MRT / VSS carry the CROSS-COLUMN fragility tint (painted in
    // `configureResultsGrid`'s `onCellPrepare` from %AUMCextrap / adj. R² /
    // N points). Both shipped tinted with a plain tooltip that never mentioned
    // the colour — user QA 2026-07-30. Their legend is attached by key below,
    // from the same declaration that drives the paint, and the guard fails the
    // build if either loses it.
    MRT: {
      title: 'MRT',
      summary:
        'Mean residence time = AUMC_0–∞ / AUC_0–∞ (single column for every route). ' +
        'For IV infusion the infusion duration is corrected out (− T_inf/2). For ' +
        'extravascular dosing this is the absorption-inclusive MRT (MRT_iv + MAT).',
    },
    VSS: {
      title: 'Vss (L/kg)',
      summary:
        'Volume of distribution at steady state = Dose·AUMC/AUC² (= CL·MRT), per kg of body ' +
        'weight. IV only — blank for extravascular data (an oral "Vss" is Vss/F confounded ' +
        'by absorption). ' +
        'Valid only under linear, time-invariant disposition: under saturable ' +
        '(Michaelis–Menten) kinetics the NCA Vss is dose-dependent and not the ' +
        'true steady-state volume.',
    },
    TLAG: plainTooltip(
      't_lag',
      'Absorption lag time: the last time with C = 0 before the first ' +
    'measurable concentration. Extravascular only — blank for IV routes.',
    ),
    LAMBDAZ_R2: plainTooltip(
      'R² (λz)',
      'Coefficient of determination of the log-linear λz regression.',
    ),
    LAMBDAZ_TFIRST: plainTooltip(
      't_first (λz)',
      'Time of the first point included in the λz regression.',
    ),
    LAMBDAZ_TLAST: plainTooltip(
      't_last (λz)',
      'Time of the last point included in the λz regression.',
    ),
    LAMBDAZ_INTERCEPT: plainTooltip(
      'ln(C₀) intercept',
      'Natural-log intercept of the log-linear λz regression.',
    ),
    AUC_METHOD: plainTooltip(
      'AUC method',
      'Trapezoidal method used for AUC integration ' +
    '(linear-up / log-down — standard for extravascular dosing).',
    ),
    COMPENSATED: plainTooltip(
      'Compensated',
      'Whether compensated summation was applied (BLQ adjustments, ' +
    'missing-point handling, etc.).',
    ),
    WARNINGS_JSON: plainTooltip(
      'Warnings',
      'NCA engine warnings emitted while computing this subject\'s profile ' +
    '(precision widening, λz-fit notes, BLQ fraction, etc.).',
    ),
    PROVENANCE_JSON: plainTooltip(
      'Provenance',
      'Calculation provenance: indices of points used for λz, ' +
      'excluded BLQ values, and other auditable metadata.',
    ),
  };

  // Attach each colour legend to its OWN column, keyed off the rule map rather
  // than named at each tooltip. Driving the loop from `buildColourRules` means a
  // coloured column cannot be skipped (there is no per-tooltip opt-in to forget)
  // and cannot receive a sibling's legend (the key is the same on both sides).
  // The completeness guard still runs — but it now backstops this loop instead
  // of being the only thing standing between a tint and its explanation.
  for (const [name, rule] of Object.entries(buildColourRules(thresholds))) {
    const spec = text[name];
    if (spec !== undefined) text[name] = {...spec, extra: legendRows(rule)};
  }
  return text;
}
