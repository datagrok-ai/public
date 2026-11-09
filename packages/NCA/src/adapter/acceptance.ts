/**
 * λz acceptance predicate — the PURE core of the results-grid `ACCEPTANCE` verdict
 * (nca-grid-status-spec §2/§3), moved out of `views/results-view.ts` so non-UI code
 * (the review-unit classifier, the sparse composite verdict, the decision fingerprint,
 * the hub stage) reads ONE predicate (nca-workflow-stages AD-4 / F-S1.1).
 *
 * **Move, not rewrite.** The logic is byte-identical to the former `readAcceptance`
 * body; `readAcceptance` in results-view is now a column-reading wrapper that feeds
 * this function through `numAt` (AC-S1.1.1: every existing verdict, reason string and
 * `_PROBLEM_RANK` is unchanged).
 *
 * No `grok` / `ui` imports — this module must stay importable headless.
 */

/** Poor span-ratio threshold, in half-lives (RG-C, user-ratified 2026-08-11). Shared by
 *  the SPAN_RATIO colour band and the SPAN criterion, so a band and its flag cannot
 *  disagree. **Deliberately NOT a field of `ResultsColorThresholds`** (AD-4 of the
 *  results-grid spec): `nca.NcaRules` has a `minSpanRatio`, and a threshold field would
 *  invite projecting it from there. Never name it `minSpanRatio` (AC-3.5). */
export const SPAN_RATIO_POOR = 2;

/** Verdict strings — the cell VALUES of the ACCEPTANCE column. */
export const ACCEPT_ACCEPTED = 'Accepted';
export const ACCEPT_FLAGGED = 'Flagged';
export const ACCEPT_NA = 'Not estimated';

export type AcceptanceVerdict = typeof ACCEPT_ACCEPTED | typeof ACCEPT_FLAGGED | typeof ACCEPT_NA;

/** One criterion's state. `n/a` when its source value is null (λz not estimated), which
 *  is distinct from `met` (an evaluated pass). */
export type FlagState = 'violated' | 'met' | 'n/a';

/** Stable identifier of a VIOLATED criterion (AC-S1.1.2). Hashed by the decision
 *  fingerprint and used to intersect bulk-review reasons — never displayed. The display
 *  phrase lives in `reasons` and may be copy-edited without moving any fingerprint. */
export const ACCEPTANCE_FLAG_CODES = ['SPAN_LT', 'ADJR2_LT', 'NPTS_LE', 'EXTRAP_GT'] as const;
export type AcceptanceFlagCode = typeof ACCEPTANCE_FLAG_CODES[number];

/** The four thresholds the predicate reads. A structural subset of
 *  `ResultsColorThresholds`, so the view passes its live Rules projection unchanged. */
export interface AcceptanceThresholds {
  readonly lambdaAcceptanceRSquared: number;
  readonly lambdaMinPoints: number;
  readonly extrapWarnPct: number;
}

/** The four source values of one profile. `NaN` (or any non-finite) = null. */
export interface AcceptanceInputs {
  readonly spanRatio: number;
  readonly adjR2: number;
  readonly nPoints: number;
  readonly pctExtrap: number;
}

export interface AcceptanceReading {
  readonly verdict: AcceptanceVerdict;
  readonly flags: {
    readonly span: FlagState; readonly rsq: FlagState;
    readonly next: FlagState; readonly extrap: FlagState;
  };
  /** One short-form phrase per VIOLATED criterion; empty unless Flagged. */
  readonly reasons: ReadonlyArray<string>;
  /** Stable codes of the violated criteria, in the same order as `reasons`. */
  readonly flagCodes: ReadonlyArray<AcceptanceFlagCode>;
  /** The pill label of each violated criterion ({@link criterionShortLabel}), same order. */
  readonly shortReasons: ReadonlyArray<string>;
}

/** Format a threshold for display AND for a DG range key. A bare `${n}` leaks IEEE-754
 *  noise (measured 2026-07-30: `0.9500000000000001` reached both a tooltip and a
 *  conditional range key). Rounding to 6 decimals moves an edge by ≤ 1e-7 of a unit,
 *  far below any threshold a user can enter. */
export function fmtThreshold(n: number): string {
  return Number.isInteger(n) ? String(n) : String(Number(n.toFixed(6)));
}

/** Short-form reason text per criterion code: criterion + threshold crossed, WITHOUT the
 *  row's measured value (the colour band already carries the magnitude). The single
 *  source for both the verdict reasons and the review picker's criterion reasons. */
export function criterionReasonText(code: AcceptanceFlagCode, t: AcceptanceThresholds): string {
  switch (code) {
  case 'SPAN_LT': return `span ratio < ${fmtThreshold(SPAN_RATIO_POOR)}`;
  case 'ADJR2_LT': return `adj. R² < ${fmtThreshold(t.lambdaAcceptanceRSquared)}`;
  case 'NPTS_LE': return `N points (λz) ≤ ${fmtThreshold(t.lambdaMinPoints)}`;
  case 'EXTRAP_GT': return `%AUCextrap > ${fmtThreshold(t.extrapWarnPct)}%`;
  }
}

/** Pill label per criterion code — the `Violated` column (acceptance-reason-pills D4). Shorter than
 *  {@link criterionReasonText} so four pills fit one grid cell; the words are those of the four
 *  `Flag:` column headers (`span` / `R²` / `N` / `%extrap`). The full phrase stays the tooltip. */
export function criterionShortLabel(code: AcceptanceFlagCode, t: AcceptanceThresholds): string {
  switch (code) {
  case 'SPAN_LT': return `span < ${fmtThreshold(SPAN_RATIO_POOR)}`;
  case 'ADJR2_LT': return `R² < ${fmtThreshold(t.lambdaAcceptanceRSquared)}`;
  case 'NPTS_LE': return `N ≤ ${fmtThreshold(t.lambdaMinPoints)}`;
  case 'EXTRAP_GT': return `%extrap > ${fmtThreshold(t.extrapWarnPct)}%`;
  }
}

/** The `Violated` column: one pill per violated λz criterion, beside ACCEPTANCE. */
export const ACCEPT_REASONS_COL = 'ACCEPT_REASONS';
/** The results-frame column names for the verdict + the four flag columns. Added
 *  PRESENTATION-side by the results view's `writeAcceptanceColumns`, NOT by the writer — the
 *  thresholds that decide their values live at presentation (the live Rules projection), so the
 *  writer and its `FIXED_RESULT_COLUMNS` drift-guard stay out of it. Here (moved from
 *  `views/results-view.ts`, nca-hub-demo-studies) so the passthrough gate's reserved set can fold
 *  them in without the adapter importing a view. */
export const ACCEPTANCE_COL = 'ACCEPTANCE';
export const FLAG_COLS = ['SPAN_FLAG', 'RSQ_FLAG', 'NEXT_FLAG', 'EXTRAP_FLAG'] as const;

/** Separates the pill labels in a `Violated` cell. No label contains it (thresholds format with
 *  `.`, see {@link fmtThreshold}). */
export const REASONS_SEPARATOR = ',';

/** The `Violated` cell value of one reading: its pill labels joined, `''` when nothing is violated. */
export function reasonsCellValue(a: AcceptanceReading): string {
  return a.shortReasons.join(REASONS_SEPARATOR);
}

/** The column tags that make Datagrok draw a `Violated` cell as pills: the platform `Tags` cell
 *  renderer (PowerGrid) splitting on {@link REASONS_SEPARATOR}, which also gives the native
 *  multi-value filter (precedent: BiostructureViewer `prolif.ts`). Tag names are the literal
 *  values of `DG.TAGS.MULTI_VALUE_SEPARATOR` / `DG.TAGS.CELL_RENDERER`, so this module needs no
 *  runtime `datagrok-api` import and stays headless. */
export function tagReasonPillsColumn<C extends {setTag(name: string, value: string): unknown; semType: string}>(
  col: C,
): C {
  col.setTag('.multi-value-separator', REASONS_SEPARATOR);
  col.setTag('cell.renderer', 'Tags');
  col.semType = 'Tags';
  return col;
}

/**
 * Evaluate the four λz acceptance criteria of one profile.
 *
 * - RSQ fires at `R_acc` applied VERBATIM (`adjR² < R_acc`; lambdaz-adjr2-threshold-spec-v2 §4).
 * - NEXT uses `≤ minPoints` (a fit AT the floor is fragile; `< minPoints` never fires since
 *   sci-comp rejects sub-minimum fits).
 * - Each criterion is `n/a` exactly when ITS value is null, so a failed λz fit yields four
 *   `n/a` flags and the `Not estimated` verdict.
 */
export function evaluateAcceptance(v: AcceptanceInputs, t: AcceptanceThresholds): AcceptanceReading {
  const rAcc = t.lambdaAcceptanceRSquared;
  const state = (finite: boolean, violated: boolean): FlagState =>
    !finite ? 'n/a' : (violated ? 'violated' : 'met');
  const spanState = state(Number.isFinite(v.spanRatio), v.spanRatio < SPAN_RATIO_POOR);
  const rsqState = state(Number.isFinite(v.adjR2), v.adjR2 < rAcc);
  const nextState = state(Number.isFinite(v.nPoints), v.nPoints <= t.lambdaMinPoints);
  const extrapState = state(Number.isFinite(v.pctExtrap), v.pctExtrap > t.extrapWarnPct);

  const flagCodes: AcceptanceFlagCode[] = [];
  if (spanState === 'violated') flagCodes.push('SPAN_LT');
  if (rsqState === 'violated') flagCodes.push('ADJR2_LT');
  if (nextState === 'violated') flagCodes.push('NPTS_LE');
  if (extrapState === 'violated') flagCodes.push('EXTRAP_GT');
  const reasons = flagCodes.map((c) => criterionReasonText(c, t));

  const states: FlagState[] = [spanState, rsqState, nextState, extrapState];
  const verdict: AcceptanceVerdict = states.every((s) => s === 'n/a') ? ACCEPT_NA :
    (states.includes('violated') ? ACCEPT_FLAGGED : ACCEPT_ACCEPTED);
  return {
    verdict,
    flags: {span: spanState, rsq: rsqState, next: nextState, extrap: extrapState},
    reasons, flagCodes,
    shortReasons: flagCodes.map((c) => criterionShortLabel(c, t)),
  };
}
