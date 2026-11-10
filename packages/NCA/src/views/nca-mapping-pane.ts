/**
 * NCA mapping pane — docked replacement for the UC-01 modal column-mapping
 * dialog. Lives on the left side of the analysis table view at dock_ratio=0.25
 * and stays open across Compute clicks (re-Compute overwrites the results grid
 * in place).
 *
 * @see archived `column-mapping-dialog.ts` (UC-01 Phase 5) — form / validator
 *      rendering / contract builder are preserved verbatim; only the container
 *      changes (dialog → docked pane).
 *
 * **Layout (spec §4 decision 3, §5.4 step 5; right area re-cut 2026-08-13):**
 * - This pane: `dockManager.dock(paneRoot, 'left', null, 'Column mapping', 0.25)`
 *   — full height; takes 25% of the table view width (mapping-pane-unit-dual-mode —
 *   restored from the UX-02 0.18 so the dual-mode unit pickers show inline by default;
 *   was 0.3 → 0.25 → 0.18 → 0.25).
 * - Main grid (the platform's built-in `view.grid`): sibling on the right; its
 *   `DockNode` is retrieved via `view.dockManager.findNode(view.grid.root)`
 *   and held as `mainGridNode` for the results-dock target.
 * - Results grid (after first Compute): `FILL` onto `mainGridNode` — a TAB of the
 *   source-data grid, not a split, so the tables row keeps its full height. Every
 *   companion TABLE joins that same group.
 * - All CHARTS live in a second group split off the top of it (`up`, 0.5), in tab
 *   order `Mean ± SD` | `Combined` | `Facet` | `Composite`. That group has no
 *   permanent member, so it is addressed through whichever of its tabs is still
 *   attached — see {@link CHARTS_ZONE_TAG} before touching any dock call here.
 *
 * **Re-Compute (spec §4 decision 5):** swap in place via
 * `existingGrid.dataFrame = results`. The DockNode and the underlying viewer
 * stay alive; only the bound DataFrame changes. Dock ratio + position
 * preserved by the platform.
 *
 * **No Cancel button (spec §4 explicit deviation from UC-01 dialog).** Closing
 * happens via the dock manager (right-click pane header → Close) or by
 * closing the entire view.
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import type * as DG from 'datagrok-api/dg';
import * as DGNS from 'datagrok-api/dg';

import {NcaController} from '../controller';
import type {PerProfileResult, OverrideAction} from '../controller';
import type {
  ColumnContract, ColumnRawType, ProfileScalars, ValidationError,
  TimeUnit, ConcentrationUnit, DoseUnit, NcaRules, BodyWeightUnit, RouteCode,
} from '../adapter';
import type {
  AnalysisRecipe, RecipeColumns, RecipeModes, RecipeUnits, SamplingType, SourceMode,
  SourcePointer, PassthroughColumnSpec, LzOverride, ReviewLedger, CompositeUnitInput,
  ConfirmedRecipe, StageCache, StageResult, ReviewUnit, LzAction, ReviewDecisionRecord,
} from '../adapter';
import {openAddToStudyDialog, updateAnalysisInPlace, saveNewStandaloneAnalysis} from './nca-add-to-study';
import type {SavedAnalysis} from './nca-add-to-study';
import {updateNcaStudyPath, updateNcaAnalysisPath} from './nca-studies-router';
// Runtime-only (the hint's study link, at click time) — same lazy-cycle safety as below.
import {openStudyOverview, UNGROUPED_SECTION} from './nca-studies-list';
// Runtime-only import (used inside promoteToStudyBacked, never at module init) — the cyclic edge
// nca-mapping-pane → nca-analysis-open is safe lazily (analysis-update-in-place F4/D7).
import {makeStudyBackedOnComputed, VIEW_ONLY_TIP} from './nca-analysis-open';
import {createWorkflowStrip, reachableStages} from './workflow-strip';
import {keepBrowseInFront, NCA_RESULTS_PRESENTED_EVENT} from './shell-panels';
import type {Subscription} from 'rxjs';
import {renderReviewPane} from './review-pane';
import {pooledCohortCue} from './results-view';
import type {ViewStage} from './workflow-strip';
import {
  DEFAULT_RULES, DEFAULT_ACCEPTANCE_RSQUARED,
  ACCEPTANCE_RSQUARED_MIN, ACCEPTANCE_RSQUARED_MAX,
  ROUTE_STRING_TO_CODE, detectMapping, isBoolLikeNumeric,
  appendGroupingColumns, aggregateNcaResults, buildMeanSdProfile,
  normalizeRouteColumn, DERIVED_ROUTE_CODE_COL,
  resolveAggregateGroupKeys, computeBioavailability, KEY_SEP, ROUTE_CODE_TO_STRING,
  blankMixedBasisFamily, IV_EV_INCOMPARABLE_FAMILY, isIntravascular,
  computeAccumulationRatios, computeExposureMargin, computeSexFold,
  computeDoseProportionality, shortDoseProportionalityNote,
  shortSparseNote, validateSparseDesign,
  detectEventDrivenColumns, normalizeEventDrivenObservations,
  SUBJECT_COL as NM_SUBJECT_COL, TAD_COL as NM_TAD_COL, CONC_COL as NM_CONC_COL,
  DOSE_COL as NM_DOSE_COL, DERIVED_INFUSION_DURATION_COL as NM_INFUSION_DUR_COL,
  poolCompositeByCell, CELL_BW_COL, mixedCompoundMessage, denseProfilesViable, reportableResults,
  excludedCompositeRows, writeReviewColumns,
  sparseUnitKey,
  LZ_DEPENDENT_PARAMS, SPARSE_LZ_DEPENDENT_COLS, SPARSE_UNIT_PREFIX,
  ACCEPT_FLAGGED, isLzDependentParam, compositeGroupKey,
  computeRecipeFingerprint, deriveStage,
  defaultSexAssignment, uniqueSexValues, applySexMap, sexMapChangesValues, SEX_SYNONYMS,
  toPlatformNull, log10OrNull, finiteMinMax, numAt, f64Column,
  formatOverrideTimestamp, AuditLog, countDoseClusters,
  findColByPattern, computeSparseForRecipe, sparseCompositeUnits, sparseDoseSource, groupByCandidates,
} from '../adapter';
import type {
  RouteColumnResult, BioavailabilityResult, ProfileIndex, AccumulationResult,
  DoseProportionalityResult, SparseNcaResult, AucBasis,
  EventDrivenDetection, EventContract, ParameterWarning, SexCode, AcceptanceThresholds,
  ObservationsCtx, SparseComputeMapping,
} from '../adapter';

import {
  configureResultsDataFrame, configureResultsGrid, configureAggregatedGrid, configureTkGrid,
  applyResultsColumnProfile, applyResultsDefaultSort, buildFacetSubjectOrder,
  REVIEW_COL, REVIEW_COMPACT_TAG, PROBLEM_RANK_COL, DEFAULT_RESULTS_THRESHOLDS,
  excludedDecisionAt, excludedCellTipSpec, EXCLUDED_TEXT, fitExcludedColumns, REVIEW_EXCLUDED,
  buildFacetVerdictDescriptions, buildFacetLambdaZParams, buildFacetVerdictColors,
  refreshAcceptanceColumnsInPlace, writeSubjectOverrideFlags, setSubjectOverrideFlag, reasonsVisibilityStale,
  OVERRIDE_UNKNOWN_USER,
  readAcceptance, gridBandTextColour, boldGridCell,
  buildHeaderTooltips, buildAcceptanceTooltipSpec, paramValueTooltipSpec,
  RESULTS_EXTENDED_TAG,
  configureBioavailabilityGrid,
  configureDoseProportionalityGrid, configureDoseProportionalityDnGrid,
  configureSparseNcaGrid, configureSparseCompositeGrid,
  PIPELINE_COL, writePipelineColumns, pipelineReasonsLine,
  REASON_PILLS_TAG, applySparseReasonPills,
} from './results-view';
import type {
  ResultsColorThresholds, FacetLambdaZParams, SparseCohort, PipelineReasonsView,
} from './results-view';
import {detectPowerGridSupport, PILLS_MIN_POWERGRID, PILLS_MIN_POWERGRID_DEBUG} from './powergrid-support';
import {richTip, bindInputTip} from './rich-tooltip';
import type {RichTipSpec} from './rich-tooltip';
import {STATUS_NEUTRAL} from './status-palette';
import {startNcaTour} from '../tour/nca-tour';
import {
  resetNotices, collectNotice, flushNotices, refreshNotices, attachBannerHost, collectNoticesFrom,
} from './compute-notices';
import type {ComputeNotice, NoticeAnchor} from './compute-notices';

/** Closed-set Route strings; sourced from the adapter's lookup so the form
 *  cannot drift out of sync with sci-comp's RouteCode mapping. */
const ROUTE_CHOICES: ReadonlyArray<string> = Object.keys(ROUTE_STRING_TO_CODE);

/** Unit dropdowns (sci-comp public types). Default order matches spec §7.1
 *  units row; the first element is the default. */
const TIME_UNIT_CHOICES: ReadonlyArray<TimeUnit> = ['h', 'min', 's', 'd'];
const CONC_UNIT_CHOICES: ReadonlyArray<ConcentrationUnit> =
  ['ng/mL', 'µg/mL', 'µg/L', 'mg/L', 'µmol/L', 'nmol/L'];
const DOSE_UNIT_CHOICES: ReadonlyArray<DoseUnit> =
  ['mg/kg', 'mg', 'µg', 'ng', 'µg/kg'];
const BODY_WEIGHT_UNIT_CHOICES: ReadonlyArray<BodyWeightUnit> = ['kg', 'g'];

/** Closed-set enums for the Method-settings section. Sourced from sci-comp
 *  `nca.AucMethod` / `nca.BlqRule` / `nca.LambdaZStrategy.mode`. Listing them
 *  inline keeps the form construction simple and avoids importing a deeper
 *  set of helper types just for `items:` arrays. The contract is enforced by
 *  the cast at `buildRules()` time. */
const AUC_METHOD_CHOICES: ReadonlyArray<string> =
  ['linear', 'log-linear', 'linear-up-log-down'];
const BLQ_RULE_CHOICES: ReadonlyArray<string> =
  ['set-zero', 'set-half-lloq', 'exclude', 'missing'];
/** Only `auto-best-fit` is wired in v1 — `manual-points` and
 *  `manual-time-range` require the interactive λz override UI (click points
 *  on the chart) which is a separate cycle (PRD VR-10..VR-12). */
const LAMBDAZ_MODE_CHOICES: ReadonlyArray<string> = ['auto-best-fit'];

/** Picker filters. `isNumeric` accepts FLOAT and INT — Phase 3 extractor
 *  handles both with `INT_NULL → NaN` translation (NOTE-04 closure). The BLQ
 *  picker uses the shared `isBoolLikeNumeric` (BOOL or INT-with-{0,1} guard)
 *  imported from the adapter — symmetric with the autodetect BLQ filter so a
 *  user whose CSV has an INT 0/1 BLQ column can ALSO override via the
 *  dropdown (would otherwise be filtered out by BOOL-only). */
const isStringOrInt = (c: DG.Column): boolean =>
  c.type === DGNS.COLUMN_TYPE.STRING || c.type === DGNS.COLUMN_TYPE.INT;
const isNumeric = (c: DG.Column): boolean =>
  c.type === DGNS.COLUMN_TYPE.FLOAT || c.type === DGNS.COLUMN_TYPE.INT;

/** view.temp key holding the persistent {@link NcaController} for this view — one
 *  instance across Computes so the interactive λz editor can recompute a single
 *  subject in place (VR-11). */
const NCA_CONTROLLER_TAG = 'ncaController';

// ── In-place "Update analysis" re-save (analysis-update-in-place, F2–F4) ──────
/** Set on `view.temp` when the open analysis is STUDY-BACKED (reopened from a study, or
 *  promoted right after Add-to-study, F4/D7) — the "Update analysis" widget's save target.
 *  Absent ⇒ the widget stays disabled. Exported so the reopen + promote paths can set it. */
export const STUDY_BINDING_TAG = 'ncaStudyBinding';
/** `true` on `view.temp` once the analysis has UNSAVED edits (any input/unit/rule/threshold
 *  change or λz point override) since the last save — drives the widget's enabled state. */
const ANALYSIS_DIRTY_TAG = 'ncaAnalysisDirty';
/** The mount closure stashes its `refreshUpdateEnabled` here so module-level {@link markAnalysisDirty}
 *  (called from the λz override sites, which are not in the closure) can repaint the widget. */
const UPDATE_REFRESH_TAG = 'ncaUpdateBtnRefresh';
/** `true` on `view.temp` while a reopened analysis is REPLAYING its saved recipe into the inputs —
 *  the replayed `input.value` setters fire `onChanged` (possibly async, hence the reopen's deferred
 *  compute), which must NOT count as a user edit. The reopen path sets this before the replay and
 *  clears it once the deferred compute has settled. Exported so nca-analysis-open can gate it. */
export const ANALYSIS_REPLAYING_TAG = 'ncaAnalysisReplaying';

/** Identity of the saved analysis the open view is bound to — the SAVE target (nca-workflow-stages
 *  F-S2.4 generalises it from "a study analysis" to "a study OR a standalone analysis"). The type
 *  name is kept (every reader already imports it); `studyId` is the CONTAINER id either way. */
export interface StudyBinding {
  /** Which store the container lives in — a study, or standalone Analyses (D4). */
  readonly container: 'study' | 'standalone';
  /** The container id (a study id, or a standalone analysis container id). */
  readonly studyId: string;
  readonly analysisId: string;
  readonly label: string;
  /** The stored snapshot the view opened (F11: a Save compares it with the stored row). */
  readonly source: SourcePointer;
  /** `computeSnapshotHash` of that snapshot — a compute on a frame with another hash ran on
   *  unsaved data (PR-2 / PR-13). Null when unknown. */
  readonly sourceHash: string | null;
  /** The persisted dialog-time passthrough (SAR/annotation) selection, re-injected into the
   *  summary cache on update so the study grid keeps its structure columns. */
  readonly passthroughCols: readonly PassthroughColumnSpec[];
  /** The container's display label (the STUDY label for a study binding) — shown in the stage hint
   *  and on the "✓ In study" button. Absent on a binding made before it was carried. */
  readonly containerLabel?: string;
  /** The caller may only view this study (SD-13): Save, Add to study and every persisting write are
   *  disabled. */
  readonly readOnly?: boolean;
}

/** Mark the open analysis dirty (has unsaved edits) and repaint the SAVE widget. Unbound views
 *  are marked too (nca-workflow-stages F-S2.4 / D2): a never-saved analysis always has something
 *  to save. Idempotent. Module-level so the λz override / review-decision sites (outside the
 *  mount closure) can call it; the closure repaints via the stashed {@link UPDATE_REFRESH_TAG}. */
export function markAnalysisDirty(view: DG.TableView): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  // No-op when already dirty, or while the view is still REPLAYING its saved recipe on reopen
  // (those onChanged are not user edits).
  if (bag[ANALYSIS_DIRTY_TAG] === true || bag[ANALYSIS_REPLAYING_TAG] === true) return;
  bag[ANALYSIS_DIRTY_TAG] = true;
  (bag[UPDATE_REFRESH_TAG] as (() => void) | undefined)?.();
}

/** The view's binding when it belongs to a STUDY (not standalone / unsaved), else null. */
function studyBindingOf(view: DG.TableView): StudyBinding | null {
  const b = (view.temp as Record<string, unknown>)[STUDY_BINDING_TAG] as StudyBinding | null;
  return b?.container === 'study' ? b : null;
}

/** True when the view holds unsaved changes (nca-workflow-stages AD-8). */
export function isAnalysisDirty(view: DG.TableView): boolean {
  return (view.temp as Record<string, unknown>)[ANALYSIS_DIRTY_TAG] === true;
}

// ── Persisted workflow state handed to the pane on reopen (nca-workflow-stages S2) ──────
/** view.temp key: the {@link PendingWorkflowState} a reopen stashes BEFORE the first compute. */
const PENDING_WORKFLOW_TAG = 'ncaPendingWorkflowState';

/** What a reopened analysis carries into the live controller: its keyed λz overrides (D22)
 *  and its review-decision ledgers (D13). Seeded once, before the first compute, so the very
 *  first compute already applies the overrides and reconciles the decisions. */
export interface PendingWorkflowState {
  readonly lzOverrides: ReadonlyArray<LzOverride>;
  readonly reviewDecisions: ReadonlyArray<ReviewLedger>;
}

/** Stash persisted workflow state for the pane's next compute (reopen / Revert). If a
 *  controller already exists (Revert after a reopen) it is applied immediately. */
export function setPendingWorkflowState(view: DG.TableView, state: PendingWorkflowState): void {
  const bag = view.temp as Record<string, unknown>;
  const controller = bag[NCA_CONTROLLER_TAG] as NcaController | null;
  if (controller != null) {
    controller.setLzOverrides(state.lzOverrides);
    controller.setReviewLedgers(state.reviewDecisions);
    return;
  }
  bag[PENDING_WORKFLOW_TAG] = state;
}

/** Seed a freshly created controller from the stashed state (once), then clear it. */
function applyPendingWorkflowState(view: DG.TableView, controller: NcaController): void {
  const bag = view.temp as Record<string, unknown>;
  const pending = bag[PENDING_WORKFLOW_TAG] as PendingWorkflowState | null;
  if (pending == null) return;
  controller.setLzOverrides(pending.lzOverrides);
  controller.setReviewLedgers(pending.reviewDecisions);
  bag[PENDING_WORKFLOW_TAG] = null;
}

/** view.temp key: the inputs of the last compute's review state ({@link ReviewCtx}). */
export const REVIEW_CTX_TAG = 'ncaReviewCtx';

/** view.temp key: the pipeline run's own reasons for this analysis ({@link PipelineReasonsView}),
 *  set when a pipeline-saved entry is opened (nca-pipeline-save-rework F10); null otherwise. */
const PIPELINE_REASONS_TAG = 'ncaPipelineReasons';

/** Set (or clear) the pipeline reasons this view shows — the `Pipeline` column of the results and
 *  composite grids, and the Review-pane line. */
export function setPipelineReasons(view: DG.TableView, pv: PipelineReasonsView | null): void {
  (view.temp as Record<string, unknown>)[PIPELINE_REASONS_TAG] = pv;
}

function pipelineReasonsOf(view: DG.TableView): PipelineReasonsView | null {
  return ((view.temp as Record<string, unknown>)[PIPELINE_REASONS_TAG] as PipelineReasonsView | null) ?? null;
}

/** The Review-pane pipeline line of one unit, or undefined. */
function pipelineLineFor(view: DG.TableView, unitKey: string): string | undefined {
  const pv = pipelineReasonsOf(view);
  const reasons = pv?.byKey.get(unitKey);
  return pv !== null && reasons !== undefined && reasons.length > 0 ? pipelineReasonsLine(pv, reasons) : undefined;
}

/** The per-compute inputs {@link NcaController.reviewSnapshot} needs beyond the controller's own
 *  state — captured once per Compute so every later reader (Review column, stage strip, save,
 *  Confirm dry run) derives the SAME units and fingerprints. */
export interface ReviewCtx {
  readonly thresholds: ResultsColorThresholds;
  readonly viable: boolean;
  readonly composites: ReadonlyArray<CompositeUnitInput & {fingerprint: string}>;
  readonly memberHashByProfile: ReadonlyMap<string, string> | undefined;
}

// ── Confirmed recipe + live stage (nca-workflow-stages F-S1.5 / F-S3.3 / AD-5) ──────────
/** view.temp key: the {@link ConfirmedRecipe} of the last successful Confirm (compute). */
const CONFIRMED_TAG = 'ncaConfirmed';
/** view.temp key: the mount closure's `getCurrentRecipe`, for module-level stage derivation. */
const RECIPE_FN_TAG = 'ncaRecipeFn';

/** The last confirmed recipe of this view, or null (never confirmed → Setup). */
export function getConfirmed(view: DG.TableView): ConfirmedRecipe | null {
  return ((view.temp as Record<string, unknown>)[CONFIRMED_TAG] as ConfirmedRecipe | null) ?? null;
}

/** Set / clear the confirmed recipe (a Confirm, a reopen of a confirmed entry). */
export function setConfirmed(view: DG.TableView, c: ConfirmedRecipe | null): void {
  (view.temp as Record<string, unknown>)[CONFIRMED_TAG] = c;
}

/** The live recipe fingerprint of this view's Setup inputs (null before the pane mounted). */
export function currentRecipeFingerprint(view: DG.TableView): string | null {
  const fn = (view.temp as Record<string, unknown>)[RECIPE_FN_TAG] as (() => AnalysisRecipe) | null;
  return fn == null ? null : computeRecipeFingerprint(fn());
}

/** The review state + derived stage of this view RIGHT NOW (D1: derived, never stored as
 *  truth). The single function the stage strip, the hint, the Review column, SAVE and the
 *  D14 dialog read — no view re-derives a stage or a count. */
export function liveWorkflowState(view: DG.TableView): {
  stage: StageResult; units: ReadonlyArray<ReviewUnit>; fingerprints: ReadonlyMap<string, string>;
} {
  const bag = view.temp as Record<string, unknown>;
  const controller = bag[NCA_CONTROLLER_TAG] as NcaController | null;
  const ctx = bag[REVIEW_CTX_TAG] as ReviewCtx | null;
  const snap = controller != null && ctx != null ?
    controller.reviewSnapshot(ctx.thresholds, {
      viable: ctx.viable, composites: ctx.composites, memberHashByProfile: ctx.memberHashByProfile,
    }) : {units: [] as ReviewUnit[], fingerprints: new Map<string, string>()};
  const confirmed = getConfirmed(view);
  const stage = deriveStage({
    confirmedFingerprint: confirmed?.fingerprint ?? null,
    currentFingerprint: currentRecipeFingerprint(view) ?? '',
    hasResults: controller?.results != null,
    units: snap.units.map((u) => ({need: u.need, decided: u.decision !== null})),
  });
  return {stage, units: snap.units, fingerprints: snap.fingerprints};
}

/** Everything the analysis entry persists about the analyst's work (v15) — read at save time
 *  and by the clean-view writeback (AD-8: written together with the state it derives from). */
export interface WorkflowPersist {
  readonly confirmed: ConfirmedRecipe | null;
  readonly lzOverrides: ReadonlyArray<LzOverride>;
  readonly reviewDecisions: ReadonlyArray<ReviewLedger>;
  readonly stageCache: StageCache | null;
}

/** The v15 workflow state of this view, for persistence. */
export function workflowPersistOf(view: DG.TableView): WorkflowPersist {
  const controller = (view.temp as Record<string, unknown>)[NCA_CONTROLLER_TAG] as NcaController | null;
  const pending = (view.temp as Record<string, unknown>)[PENDING_WORKFLOW_TAG] as PendingWorkflowState | null;
  const {stage} = liveWorkflowState(view);
  return {
    confirmed: getConfirmed(view),
    // Before the first compute the state still sits in the pending stash (a Setup-stage save of
    // a reopened analysis must not drop it).
    lzOverrides: controller?.getLzOverrides() ?? pending?.lzOverrides ?? [],
    reviewDecisions: controller?.getReviewLedgers() ?? pending?.reviewDecisions ?? [],
    stageCache: {stage: stage.stage, pending: stage.pending, notApplicable: stage.notApplicable},
  };
}

/** The view controller's per-row λz action lookup (AD-2), for the facet markers; undefined
 *  before the first compute. */
function lzActionOfView(view: DG.TableView): ((rowIdx: number) => LzAction) | undefined {
  const controller = (view.temp as Record<string, unknown>)[NCA_CONTROLLER_TAG] as NcaController | null;
  return controller == null ? undefined : (r: number) => controller.currentLzAction(r);
}

// ── Stage navigation (nca-workflow-stages S3) ────────────────────────────────────────────
/** view.temp key: the last compute's presentation closure — re-docks every result panel from
 *  the HELD controller results without recomputing (a no-change return to Review/Report, P6). */
const PRESENT_TAG = 'ncaPresentFn';
/** view.temp key: the mount closure's stage refresh (strip + ribbon + Review column), so
 *  module-level edit sites (a λz point toggle, a review decision) keep the strip current. */
const STAGE_REFRESH_TAG = 'ncaStageRefresh';
/** view.temp key: the stage currently SHOWN (ViewStage). */
export const VIEW_STAGE_TAG = 'ncaViewStage';

/**
 * The `Review` column is stage-scoped (D15 / approved design): visible on Review λz, and in the
 * Full profile; hidden in Compact on Report and on Setup. Grid-side visibility only — the frame
 * keeps the column. No-op before the column exists.
 */
function applyReviewColumnVisibility(view: DG.TableView): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const grid = bag[RESULTS_GRID_TAG] as DG.Grid | null;
  const df = grid?.dataFrame ?? null;
  if (grid == null || df == null || df.col(REVIEW_COL) === null) return;
  // The profile reads this frame tag, so a later Compact/Full toggle keeps the stage scoping.
  df.setTag(REVIEW_COMPACT_TAG, viewStageOf(view) === 'review' ? 'true' : 'false');
  applyResultsColumnProfile(grid, resultsExtended(view));
}

/**
 * Write / refresh the Review column (+ hidden detail) from the controller's review snapshot —
 * after a compute's reconcile, after a decision, after a λz point edit. The ONE writer; the grid,
 * the strip and the Review pane all read the same snapshot (AC-S3.1.1).
 */
export function refreshReviewColumns(view: DG.TableView): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  refreshSparseExclusions(view); // the composite table reads the same review state (SF-1 (b))
  const controller = bag[NCA_CONTROLLER_TAG] as NcaController | null;
  const results = controller?.results ?? null;
  if (results == null) return;
  const added = writeReviewColumns(results, liveWorkflowState(view).units);
  // The pipeline run's own reasons (F10): written once per results frame (a compute brings a new
  // frame), keyed by SUBJECT like the Review column.
  const pv = pipelineReasonsOf(view);
  const subj = results.col('SUBJECT');
  const pipelineChanged = subj !== null && (results.col(PIPELINE_COL) === null) !== (pv === null) &&
    writePipelineColumns(results, pv, (i) => subj.getString(i));
  const grid = bag[RESULTS_GRID_TAG] as DG.Grid | null;
  if (grid == null) return;
  if (added || pipelineChanged) applyReviewColumnVisibility(view);
  // An excluded value must read whole as `{value} (excl.)` — widen those columns if needed.
  const rc = results.col(REVIEW_COL);
  fitExcludedColumns(grid, LZ_DEPENDENT_PARAMS,
    (i) => rc !== null && !rc.isNone(i) && rc.getString(i) === REVIEW_EXCLUDED);
  grid.invalidate();
}

/** Refresh this view's stage UI (no-op before the pane mounted). */
export function refreshStageUi(view: DG.TableView): void {
  ((view.temp as Record<string, unknown>)[STAGE_REFRESH_TAG] as (() => void) | null)?.();
}

/** The stage currently shown in this view. */
export function viewStageOf(view: DG.TableView): ViewStage {
  return ((view.temp as Record<string, unknown>)[VIEW_STAGE_TAG] as ViewStage | null) ?? 'setup';
}

/** Unit-key prefix of a UC-04 sparse composite (`sparse:<group>`, AD-1). */
const SPARSE_KEY_PREFIX = SPARSE_UNIT_PREFIX;

/** The reportable (Exclude-λz-blanked) clone of a results frame for THIS view's current
 *  decisions (AD-6) — the one frame every aggregate consumer reads. */
export function reportableFrameOf(view: DG.TableView, results: DG.DataFrame): DG.DataFrame {
  const controller = (view.temp as Record<string, unknown>)[NCA_CONTROLLER_TAG] as NcaController | null;
  if (controller == null) return results;
  return reportableResults(results, controller.excludedKeys()).df;
}


/** UI sentinel for an unassigned per-value sex choice (the choice-input
 *  placeholder). Distinct from any `SexCode` so the Compute-disable gate detects
 *  "not yet mapped". Pane-only — not part of the adapter map contract. */
const SEX_UNASSIGNED = '—';

/** Max distinct values for a plausible sex axis. Above this the column is not a
 *  sex column (Compute is disabled rather than building N inputs). Pane-only. */
const SEX_MAX_DISTINCT = 8;

/** Display label of each value/column source mode. The radio shows the label; the recipe keeps
 *  the lowercase `SourceMode`, so every read and write goes through `readSourceMode` /
 *  `writeSourceMode` (mapping-source-radio, user-directed 2026-10-19). */
const SOURCE_MODE_LABEL = {value: 'Value', column: 'Column'} as const;

/** A `<X> type` value/column switch: a horizontal radio, laid out like the Add-to-study
 *  `Study` radio (`.nca-radio-row`). */
function sourceModeInput(label: string): DG.InputBase<string | null> {
  const input = ui.input.radio(label, {
    items: [SOURCE_MODE_LABEL.value, SOURCE_MODE_LABEL.column], value: SOURCE_MODE_LABEL.value,
  });
  input.root.classList.add('nca-radio-row');
  return input;
}

function readSourceMode(input: DG.InputBase<string | null>): Exclude<SourceMode, 'none'> {
  return input.value === SOURCE_MODE_LABEL.column ? 'column' : 'value';
}

/** `'none'` (a recipe with no source set) shows as `Value`, the form's default. */
function writeSourceMode(input: DG.InputBase<string | null>, mode: SourceMode): void {
  input.value = mode === 'column' ? SOURCE_MODE_LABEL.column : SOURCE_MODE_LABEL.value;
}

/**
 * Build the canonical `ColumnContract` over the event-driven NORMALIZED frame
 * (AD-1). The transform emits fixed column names (Subject / TAD / Conc / Dose_mg
 * / ~RouteCode) plus the hoisted body-weight column (under its source name).
 * Concentration carries the user-selected source unit; AMT is absolute mg (AC-A6).
 */
function buildNonmemContract(
  srcDf: DG.DataFrame, concUnit: ConcentrationUnit, bodyWeightName: string | null,
): ColumnContract {
  const obs: {-readonly [K in keyof ColumnContract['observations']]: ColumnContract['observations'][K]} = {
    subjectId: {name: NM_SUBJECT_COL, required: true, rawType: 'string'},
    time: {name: NM_TAD_COL, required: true, rawType: 'float64', units: 'h'},
    concentration: {name: NM_CONC_COL, required: true, rawType: 'float64', units: concUnit},
    dose: {name: NM_DOSE_COL, required: false, rawType: 'float64', units: 'mg'},
    route: {name: DERIVED_ROUTE_CODE_COL, required: false, rawType: 'int32'},
  };
  if (bodyWeightName !== null && srcDf.col(bodyWeightName) !== null)
    obs.bodyWeight = {name: bodyWeightName, required: false, rawType: 'float64', units: 'kg'};
  // Per-subject infusion duration (Slice A2 / BD-2): the transform emits
  // `~InfusionDuration` (NaN for RATE=0 subjects → extractor falls back to the
  // scalar). Reference it so RATE>0 IV-infusion subjects get their own duration.
  if (srcDf.col(NM_INFUSION_DUR_COL) !== null)
    obs.infusionDuration = {name: NM_INFUSION_DUR_COL, required: false, rawType: 'float64', units: 'h'};
  return {observations: obs};
}

/**
 * AD-5 mandatory route-confirm dialog. When a NONMEM file has no explicit ROUTE
 * column, route is inferred from RATE/CMT to {IV-bolus, IV-infusion, extravascular}
 * — but PO vs SC is not recoverable, and a wrong IV/extravascular label mislabels
 * CL vs CL/F (~3-fold error at F≈30%). Block Compute until the scientist confirms.
 */
function showRouteConfirmDialog(onConfirm: (code: RouteCode) => void): void {
  const choice = ui.input.choice<string>('Route', {value: 'IV-bolus', items: [...ROUTE_CHOICES]});
  const msg = ui.divText(
    'No explicit ROUTE column was found — route was inferred from RATE / CMT and ' +
    'cannot distinguish PO from SC. Confirm the administration route; it sets ' +
    'whether CL or CL/F is reported.');
  msg.style.maxWidth = '420px';
  msg.style.whiteSpace = 'normal';
  ui.dialog('Confirm route — NONMEM ingestion')
    .add(msg).add(choice.root)
    .onOK(() => onConfirm(ROUTE_STRING_TO_CODE[choice.value as string]))
    .show();
}

/** `view.temp` keys for the docked results grid handle. Re-Compute uses these
 *  to swap the DataFrame in place rather than re-docking.
 *
 *  Since the 2026-08-13 layout the results grid FILL-docks onto the SOURCE-data
 *  grid, so `RESULTS_NODE_TAG` names the shared TABLES group
 *  (`study.csv` | `NCA results` | `By group` | `Bioavailability` | …). Every
 *  companion table already anchors here and needed no change — the charts moved,
 *  the tables did not. */
export const RESULTS_GRID_TAG = 'ncaResultsGrid';
const RESULTS_NODE_TAG = 'ncaResultsNode';

/**
 * The results grid's column profile, resolved from `view.temp` — **the DEFAULT
 * lives here and nowhere else.** Four call sites read this state, and a restated
 * `=== true` in each of them is exactly how a default drifts apart.
 *
 * Unset → `false` (**Compact**): the grid opens on the scan profile.
 */
function resultsExtended(view: DG.TableView): boolean {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  // `== null` catches BOTH the dart-side null of a never-set key (BUG-06) and
  // undefined; only an explicit `true` selects Full.
  const stored = bag[RESULTS_EXTENDED_TAG] as boolean | null | undefined;
  return stored === true;
}

/** view.temp key: the scientist's `Violated` column switch (Settings), once they flipped it. */
const REASON_PILLS_VIEW_TAG = 'ncaReasonPills';

/**
 * Is the `Violated` pills column shown? The scientist's switch when set; otherwise the app-start
 * PowerGrid check ({@link detectPowerGridSupport}) — off on a PowerGrid that cannot draw the pills
 * readably. **The default lives here and nowhere else** (the {@link resultsExtended} rule).
 */
function reasonPillsShown(view: DG.TableView): boolean {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const stored = bag[REASON_PILLS_VIEW_TAG] as boolean | null | undefined;
  return stored == null ? detectPowerGridSupport().pills : stored;
}

/** Apply the `Violated` switch to the results grid (frame tag read by the column profile) and to
 *  the Sparse NCA grid. Safe before either grid exists. */
function applyReasonPillsVisibility(view: DG.TableView): void {
  const show = reasonPillsShown(view);
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const grid = bag[RESULTS_GRID_TAG] as DG.Grid | null;
  const df = grid?.dataFrame ?? null;
  if (grid != null && df != null) {
    df.setTag(REASON_PILLS_TAG, show ? 'true' : 'false');
    applyResultsColumnProfile(grid, resultsExtended(view));
  }
  const sparse = bag[SPARSE_PARAM_GRID_TAG] as DG.Grid | null;
  if (sparse != null) applySparseReasonPills(sparse, show);
}

/**
 * `view.temp` key for the CHARTS group: the list of dock nodes currently living
 * in the top half (`Mean ± SD` | `Combined` | `Facet` | `Composite`).
 *
 * Post-Compute layout (user-directed 2026-08-13):
 *
 * ```
 * ┌ Column mapping ┬ Mean ± SD │ Combined │ Facet │ Composite ┐
 * │      0.25      │        charts group (up, 0.5)            │
 * │                ├ study.csv │ ●NCA results │ By group │ …  ┤
 * │                │        tables group (RESULTS_NODE_TAG)   │
 * └────────────────┴──────────────────────────────────────────┘
 * ```
 *
 * **Why a LIST and not a single anchor node.** A tab group has no handle of its
 * own in this API — you address it through any member. Before this layout the
 * charts group was addressed through `view.grid`'s node, and that was safe for a
 * reason worth naming: the source-data tab is PERMANENT, so the anchor could
 * never die under us. This layout moves that permanent member down to the tables
 * group, leaving a charts group in which **every** member is conditionally
 * closed — `Mean ± SD` when the dense frame loses its time/conc column or the
 * profile comes back empty, `Combined` when a mapping input is cleared, `Facet`
 * on a subject-set change or >20 subjects, `Composite` when the design stops
 * being sparse. Pinning one node as "the anchor" would therefore hand the next
 * dock a dead handle, and `dock(FILL, <dead node>)` does not reliably throw — it
 * can land the panel at the view root instead, which reads as a layout the user
 * never asked for rather than as an error.
 *
 * So the zone is the set of its members, each dock re-picks the first one the dock
 * manager still resolves ({@link isDockedPanel} — NOT "attached to the document",
 * which is a different and wrong test; see that function), and the group is
 * re-created only when genuinely nothing survives.
 */
const CHARTS_ZONE_TAG = 'ncaChartsPanels';

/** Height the charts group takes when it is created, as a fraction of the tables
 *  group it splits off. Written once; the user drags from here. */
const CHARTS_ZONE_RATIO = 0.5;

/**
 * Is this panel still docked anywhere?
 *
 * **`dockManager.findNode` is the authority, and the two obvious alternatives are
 * both WRONG — measured on the real host, not reasoned about** (see
 * `dock-layout-tests` "PROBE"):
 *
 * | probe | single panel | after a 2nd tab joins its group |
 * |---|---|---|
 * | `container.containerElement.isConnected` | `true` | **`false` for the inactive tab** |
 * | `findNode(container.containerElement)` | **`false`** | `false` |
 * | `findNode(<element passed to dock()>)` | `true` | `true` |
 *
 * So `isConnected` does not mean "alive" — it means "is the ACTIVE tab": dock-spawn
 * detaches an inactive tab's content from the document. A liveness check built on
 * it prunes every background chart, and once it prunes them all, the next dock
 * silently opens a SECOND charts row beside the first. That is exactly the defect
 * the first draft of this layout would have shipped.
 *
 * And the node must be looked up by the element ORIGINALLY PASSED TO `dock()` —
 * `container.containerElement` is a different element that `findNode` does not
 * resolve even for a healthy single panel. Hence this module tracks the charts
 * group as the LIST OF ELEMENTS IT DOCKED, never as `DockNode` handles.
 *
 * Exported for `dock-layout-tests`, which re-measures the table above on every
 * run so a platform change cannot quietly invalidate it.
 */
export function isDockedPanel(view: DG.TableView, el: HTMLElement): boolean {
  try {
    return (view.dockManager.findNode(el) ?? null) !== null;
  } catch {
    return false;
  }
}

/**
 * Bring `NCA results` to the front of the tables group.
 *
 * **Why this is needed even though results dock FIRST.** They have to: the
 * results panel is what establishes the tables group every companion table then
 * FILLs into, and what the charts row splits off. So dock order cannot be
 * reversed — but the consequence of docking first is that the LAST decorator to
 * arrive owns the visible tab, and on a dose-ranging study that is
 * `Dose-Normalized`. The scientist's Compute finished with the main result
 * buried three tabs deep (reported 2026-08-13).
 *
 * The API path is measured, not assumed (`dock-layout-tests`, "a background tab
 * can be brought back to the front"): `setActiveChild` lives on the PARENT's
 * container and takes the target node's container. Backgrounding is not closing —
 * the companion tables stay docked and one click away.
 */
/** Bring the docked panel whose node is stored under `nodeTag` in `view.temp` to the
 *  front of its tab group. Shared by {@link activateResultsTab} and
 *  {@link activateMeanSdTab}: the tab-focus mechanism is identical (only the tag
 *  differs), so it lives once rather than being copied per tab. A null node (nothing
 *  docked under the tag yet / tab closed) or an unexpected dock shape is a cosmetic
 *  no-op, never a thrown error. The focus call itself
 *  (`node.parent.container.setActiveChild(node.container)`) is pinned by
 *  dock-layout-tests ("a background tab can be brought back to the front"). */
function activateDockedTab(view: DG.TableView, nodeTag: string): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const node = bag[nodeTag] as DG.DockNode | null;
  if (node == null) return; // nothing docked under this tag yet
  try {
    node.parent?.container?.setActiveChild(node.container);
  } catch {
    // Dock tree not shaped as expected on this build. Cosmetic: the tab stays where
    // the platform left it rather than the Compute / selection change failing.
  }
}

/** Front the "NCA results" tab — the main result must not stay buried when a later
 *  decorator docks last into the same group (see dock-layout-tests). */
function activateResultsTab(view: DG.TableView): void {
  activateDockedTab(view, RESULTS_NODE_TAG);
}

/** Front the "Mean ± SD" chart tab — used when a Group-by selection changes, so the
 *  scientist lands on the visualization of the split they just chose (user-directed
 *  2026-08-28). Exported for dock-layout-tests (same precedent as {@link isDockedPanel}
 *  being exported for that suite). */
export function activateMeanSdTab(view: DG.TableView): void {
  activateDockedTab(view, MEAN_SD_CHART_NODE_TAG);
}

/**
 * Entering Report (user-directed 2026-09-30): the `Mean ± SD` tab comes to the front, and the context
 * panel shows the `Combined` concentration–time chart in an OPEN `Profiles` group — the population
 * overview the report is about, beside the per-group summary. The panel chart is a second native
 * line chart on the Combined chart's own frame + options (one viewer cannot sit in two places).
 * BUILD-ONCE, re-bound on every entry: the chart is never docked in a view, so `viewer.close()` on
 * it throws (no dock/view parent) — reuse, like the λz readout panel, never close-and-recreate.
 */
function showReportPanel(view: DG.TableView, activateTab = true): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  if (activateTab) activateMeanSdTab(view);
  const combined = bag[PROFILE_CHART_TAG] as DG.Viewer | null;
  const opts = bag[PROFILE_CHART_OPTS_TAG] as Partial<DG.ILineChartSettings> | null;
  if (combined == null || opts == null || combined.dataFrame == null) return;
  let panel = bag[REPORT_PANEL_TAG] as {
    root: HTMLElement; chart: DG.Viewer; pane: DG.AccordionPane; title: HTMLElement; host: HTMLElement;
    /** The Combined viewer whose resizes the panel follows, and that subscription. */
    observed: DG.Viewer | null; sizeSub: Subscription | null;
  } | null;
  // Legend at the TOP RIGHT (user-directed 2026-10-13, after Bottom and Right). `Always`, not
  // the default `Auto`: under Auto the panel's legend was squeezed into a space too small for its
  // labels (clipped, with an inner scrollbar — QA 2026-10-13). The study overview's drill-in chart
  // pairs its legend with `Always` the same way (nca-study-overview.ts).
  const panelOpts: Partial<DG.ILineChartSettings> =
    {...opts, legendVisibility: 'Always', legendPosition: 'RightTop'};
  if (panel == null) {
    const chart = DGNS.Viewer.lineChart(combined.dataFrame, panelOpts);
    // A definite height: a line chart has no intrinsic size and collapses in an auto-height pane
    // (the same recipe as the λz readout's profile pane).
    const chartHost = ui.div([chart.root]);
    chartHost.style.width = '100%';
    // The SAME height as the docked Combined chart (user-directed 2026-10-13), kept in step when
    // the dock is resized; see `matchDockedHeight` + the re-subscription below.
    chartHost.style.height = `${REPORT_PANEL_FALLBACK_HEIGHT}px`;
    chart.root.style.width = '100%';
    chart.root.style.height = '100%';
    const acc = ui.accordion();
    // The table's name as the panel title, styled exactly like the λz readout's subject title
    // (user-directed 2026-09-30) — the panel says WHAT it shows.
    const title = panelTitleEl();
    acc.addTitle(title);
    const pane = acc.addPane('Profiles', () => chartHost, true);
    acc.root.style.width = '100%';
    panel = {root: acc.root, chart, pane, title, host: chartHost, observed: null, sizeSub: null};
    bag[REPORT_PANEL_TAG] = panel;
  } else {
    // Re-bind to the current frame + options (a re-Compute or a semi-log switch changed them).
    panel.chart.dataFrame = combined.dataFrame;
    panel.chart.setOptions(panelOpts);
  }
  panel.title.textContent = view.dataFrame?.name ?? view.name;
  panel.pane.expanded = true;
  // Follow THIS Combined viewer's resizes: a re-Compute may have re-docked a new one, so the
  // subscription moves with it (architect review 2026-10-13 — a once-only subscription stayed on
  // the first viewer). Re-measured on every entry as well.
  if (panel.observed !== combined) {
    panel.sizeSub?.unsubscribe();
    const host = panel.host;
    panel.sizeSub = ui.onSizeChanged(combined.root).subscribe(() => matchDockedHeight(combined, host));
    panel.observed = combined;
  }
  matchDockedHeight(combined, panel.host);
  const root = panel.root;
  setTimeout(() => {
    try {
      grok.shell.windows.showContextPanel = true;
      grok.shell.o = root;
    } catch {/* no shell */}
  }, 0);
}

/** The Report panel chart's height before the docked Combined chart has been laid out. */
const REPORT_PANEL_FALLBACK_HEIGHT = 450;

/** Size the Report panel chart to the docked Combined chart's height (user-directed 2026-10-13).
 *  On Report the Combined tab usually sits behind `Mean ± SD` (height 0), so the nearest laid-out
 *  ancestor is measured instead — the tab area both charts share. Unchanged while nothing is laid
 *  out yet (the fallback stands). */
function matchDockedHeight(combined: DG.Viewer, host: HTMLElement): void {
  let el: HTMLElement | null = combined.root;
  while (el !== null && el.offsetHeight === 0) el = el.parentElement;
  const h = el?.offsetHeight ?? 0;
  if (h > 0) host.style.height = `${h}px`;
}

/** A context-panel title element: 18px / --grey-4 / single-line ellipsis — shared by the λz readout
 *  (subject) and the Report panel (table name), so the two titles cannot drift. */
function panelTitleEl(): HTMLElement {
  const title = ui.div([]);
  Object.assign(title.style, {
    fontSize: '18px', height: '18px', color: 'var(--grey-4, #9497A0)', display: 'block',
    whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis',
    cursor: 'default', userSelect: 'all',
  });
  return title;
}

/** The TABLES group — the node every companion table FILLs into and the charts
 *  group splits off. Prefers the handle `swapOrDockResults` stored over
 *  `findNode`, which this file's own history records as unreliable once the dock
 *  tree has been rearranged; `findNode` stays as the last-resort fallback.
 *
 *  `explicit` is the node threaded down from `runCompute` — the same value, passed
 *  rather than looked up, by decorators that already hold it. */
function resultsAnchor(
  view: DG.TableView, explicit: DG.DockNode | null = null,
): DG.DockNode | null {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  // `!= null`: `view.temp` returns dart-side null for a never-set key (BUG-06).
  return (explicit ?? bag[RESULTS_NODE_TAG] ?? view.dockManager.findNode(view.grid.root)) ?? null;
}

/**
 * Dock a chart as a tab of the charts group, creating that group on first use.
 *
 * Tab order is dock order, so the call sequence in `runCompute` IS the tab
 * sequence the user sees: `Mean ± SD`, `Combined`, `Facet` (and `Composite`,
 * which docks later from the sparse path).
 */
function dockChartTab(view: DG.TableView, el: HTMLElement, title: string): DG.DockNode {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  // Re-resolve every remembered panel through the dock manager on each call. A
  // chart closed since the last dock drops out here; a chart that is merely
  // sitting in the background stays, because {@link isDockedPanel} asks the dock
  // manager rather than the DOM. `!= null`: `view.temp` returns dart-side null
  // for a never-set key (BUG-06).
  const known = (bag[CHARTS_ZONE_TAG] as HTMLElement[] | null) ?? [];
  const live: HTMLElement[] = [];
  let host: DG.DockNode | null = null;
  for (const panel of known) {
    const node = view.dockManager.findNode(panel) ?? null;
    if (node === null) continue;
    live.push(panel);
    host ??= node;
  }
  const node = host !== null ?
    view.dockManager.dock(el, DGNS.DOCK_TYPE.FILL, host, title) :
    view.dockManager.dock(
      el, DGNS.DOCK_TYPE.TOP, resultsAnchor(view), title, CHARTS_ZONE_RATIO);
  live.push(el);
  bag[CHARTS_ZONE_TAG] = live;
  return node;
}
/** `view.temp` keys for the per-group aggregated grid (#11). Docked as a
 *  FILL-tab on the results node ("By group" next to "NCA results"); swapped in
 *  place on re-Compute / Group-by change, closed when Group-by → "(none)". */
const AGG_GRID_TAG = 'ncaAggGrid';
const AGG_NODE_TAG = 'ncaAggNode';
/** `view.temp` keys for the stratum-level F% (bioavailability) grid (per-row-
 *  route Slice D). Docked as a FILL-tab on the results node ("Bioavailability")
 *  when the frame carries both IV and extravascular routes + a dose + a Sex
 *  metadata column; swapped in place on re-Compute, closed when not applicable. */
const FPCT_GRID_TAG = 'ncaFpctGrid';
const FPCT_NODE_TAG = 'ncaFpctNode';
/** `view.temp` keys for the toxicokinetics aggregate grids (tk-aggregates AD-6).
 *  Docked as FILL-tabs on the results node ("Accumulation" / "Exposure Margin" /
 *  "Sex Comparison") when the frame is repeat-dose (≥ 2 StudyDay values) + a Sex
 *  metadata column; swapped in place on re-Compute, closed when not applicable. */
const ACC_GRID_TAG = 'ncaAccGrid';
const ACC_NODE_TAG = 'ncaAccNode';
const MARGIN_GRID_TAG = 'ncaMarginGrid';
const MARGIN_NODE_TAG = 'ncaMarginNode';
const SEXFOLD_GRID_TAG = 'ncaSexFoldGrid';
const SEXFOLD_NODE_TAG = 'ncaSexFoldNode';
/** `view.temp` keys for the dose-proportionality grids (Smith power model,
 *  FR-410..413 / UC-03). Docked as FILL-tabs on the results node ("Dose
 *  Proportionality" = the per-Route×Sex β / 90% CI / verdict table; "Dose-
 *  Normalized" = the per-subject dn-AUC / dn-Cmax table) when the frame carries
 *  a dose column + a Sex metadata column AND spans ≥ 2 dose levels (a dose-
 *  ranging design); swapped in place on re-Compute, closed when not applicable.
 *  The P4 log-log scatter + CI-band and dn-AUC bar charts are a separate
 *  step. */
const DP_GRID_TAG = 'ncaDosePropGrid';
const DP_NODE_TAG = 'ncaDosePropNode';
const DPDN_GRID_TAG = 'ncaDosePropDnGrid';
const DPDN_NODE_TAG = 'ncaDosePropDnNode';
/** `view.temp` keys for the sparse / destructive-sampling panels (UC-04). The
 *  honest-uncertainty parameter grid, the composite-profile grid, and the VR-07
 *  composite scatter (raw individuals + mean overlay). */
export const SPARSE_PARAM_GRID_TAG = 'ncaSparseParamGrid';
const SPARSE_PARAM_NODE_TAG = 'ncaSparseParamNode';
const SPARSE_COMP_GRID_TAG = 'ncaSparseCompGrid';
const SPARSE_COMP_NODE_TAG = 'ncaSparseCompNode';
const SPARSE_CHART_TAG = 'ncaSparseChart';
const SPARSE_CHART_NODE_TAG = 'ncaSparseChartNode';
/** `view.temp` keys for the VR-03 group mean ± SD profile chart — a FILL-tab
 *  sibling to the Combined / Facet / Composite tabs.
 *  - CHART: the `DG.Viewer` (scatterPlot with per-group Mean + ±SD lines).
 *  - NODE: its `DG.DockNode`, kept so {@link closeMeanSdChart} can close it.
 *  - CTX: the post-Compute render context ({@link MeanSdCtx}) — the normalized
 *    observations frame + resolved time/conc/dose names — stored so the Group-by
 *    control can re-split the chart WITHOUT a recompute. The chart reflects the
 *    LAST RUN, like the other profile tabs (see {@link PROFILE_CONC_NAME_TAG}). */
export const MEAN_SD_CHART_TAG = 'ncaMeanSdChart';
export const MEAN_SD_CHART_NODE_TAG = 'ncaMeanSdChartNode';
const MEAN_SD_CTX_TAG = 'ncaMeanSdCtx';
/** `view.temp` key for the full (group × time) profile frame backing the chart —
 *  kept so the semi-log toggle can rebuild the per-scale chart frame (the log
 *  frame drops all-BLQ rows; see {@link meanSdChartFrame}) without a recompute. */
const MEAN_SD_PROFILE_TAG = 'ncaMeanSdProfile';
/** `view.temp` key for the rxjs subscription bridging the Accumulation grid's
 *  current (Dose×Sex) cell → the Combined profile-overlay filter (AD-5). */
const TK_OVERLAY_SUB_TAG = 'ncaTkOverlaySub';
/** `view.temp` key for the concentration-time profile chart. The chart lives
 *  in a `FILL`-docked tab next to the main data grid (same central area). */
export const PROFILE_CHART_TAG = 'ncaProfileChart';
/** `view.temp` key: the options the `Combined` chart was last rendered with — so the Report context
 *  panel's copy ({@link showReportPanel}) is the SAME chart, never a re-derived look-alike. */
const PROFILE_CHART_OPTS_TAG = 'ncaProfileChartOpts';
/** `view.temp` key: the Report context panel, built once and reused ({@link showReportPanel}). */
export const REPORT_PANEL_TAG = 'ncaReportPanel';
/** `view.temp` key for the concentration column name the profile chart + facet
 *  grid were last BUILT with (i.e. at the last Compute). The charts reflect the
 *  last run, not the live mapping — they are not rebuilt on input change — so
 *  the semi-log toggle must re-bind against THIS column, not `concInput.value`,
 *  which the user may have edited since (the stale chart frame would not contain
 *  the new column). */
const PROFILE_CONC_NAME_TAG = 'ncaProfileConcName';
/** `view.temp` key for the rxjs subscription that bridges
 *  `results.selection` → `chartDf.filter` (see
 *  {@link wireResultsSelectionToProfileFilter}). Stored so each Compute can
 *  dispose the previous wiring before binding a new one — the `results`
 *  DataFrame instance is replaced on every Compute, and the old subscription
 *  points at a stale source. */
const SELECTION_SUB_TAG = 'ncaResultsSelectionSub';

/** `view.temp` keys for the per-subject facet grid (a stack of scatterPlots,
 *  one per subject, packed into an `ui.splitH/splitV` grid). The grid lives
 *  as a FILL-tab next to the overlay {@link PROFILE_CHART_TAG}.
 *  - DIV: the outer `HTMLDivElement` returned by {@link buildFacetGrid}.
 *  - NODE: its `DG.DockNode`, kept so {@link ensureFacetGrid} can close it.
 *  - PLOTS: the `DG.Viewer[]` (one per subject), for `dataFrame =`-swap when
 *    the subject set is unchanged between Computes / mapping refreshes, and
 *    for the dimming bridge.
 *  - SUBJECTS: the ordered subject-key list the current grid was built from.
 *    Compared on subsequent ensureFacetGrid calls to decide
 *    swap-in-place vs. full rebuild. */
export const FACET_GRID_DIV_TAG = 'ncaFacetGridDiv';
const FACET_GRID_NODE_TAG = 'ncaFacetGridNode';
const FACET_PLOTS_TAG = 'ncaFacetPlots';
const FACET_SUBJECTS_TAG = 'ncaFacetSubjects';
/** Per-tile point-role column (M1): 'included' / 'excluded' in the λz fit. Bound
 *  as the scatterPlot `markersColumnName`; marker shapes are THROWAWAY (final
 *  treatment is a Rule-3 design decision). */
const FACET_ROLE_COL = '_lz_role';
/** Per-tile hidden source-df row index per point (VR-11 C0): maps a clicked facet
 *  point back to its observations-frame row for the λz override. */
const FACET_SRCROW_COL = '_srcRow';
/** `view.temp` key holding the closure that rebuilds the facet from the
 *  controller's current results frame after an interactive λz override (C0). */
const FACET_REFRESH_TAG = 'ncaFacetRefresh';
/** `view.temp` key for the AC-7.2 dependent-aggregate refresh closure: re-runs the
 *  By-group / %F / dose-prop / TK decorators against the patched results frame after
 *  an interactive λz override, so the aggregates never lag the per-subject numbers.
 *  Overwritten each Compute (captures that Compute's decorator inputs). */
const AGG_REFRESH_TAG = 'ncaAggRefresh';
/** `view.temp` key for the per-tile click→override subscriptions (C0). Disposed +
 *  re-established whenever the facet frames are (re)bound. */
const FACET_CLICK_SUBS_TAG = 'ncaFacetClickSubs';
/** `view.temp` re-entrancy guard: true while an override-triggered facet refresh
 *  is in flight, so the frame swap's currentRow reset does not recurse. */
const FACET_CLICK_BUSY_TAG = 'ncaFacetClickBusy';
/** `view.temp` keys for the single-subject λz editor (F4, docked tab — a WALKING
 *  SKELETON for the design look). The editor is one enlarged facet-style plot bound
 *  to the results-grid current row, reusing the whole C0/C1 override + refresh path. */
const EDITOR_NODE_TAG = 'ncaLzEditorNode';
export const EDITOR_PLOT_TAG = 'ncaLzEditorPlot';
const EDITOR_HEADER_TAG = 'ncaLzEditorHeader';
const EDITOR_SUBJECT_TAG = 'ncaLzEditorSubject';
/** The docked editor CONTENT element (`divV([header, plot.root])`). Kept so a re-Compute
 *  can ask the dock manager whether the editor is still live (`findNode`) and REUSE it in
 *  place instead of close()+re-dock — the duplicate-tab footgun ensureFacetGrid documents. */
const EDITOR_CONTENT_TAG = 'ncaLzEditorContent';
/** `view.temp` keys for the curated λz context-panel readout: the build-once accordion
 *  shell (root + subject-name title + content host) and the `(subject) => void` closure
 *  that refreshes it. */
export const LZ_READOUT_PANEL_TAG = 'ncaLzReadoutPanel';
export const LZ_READOUT_TAG = 'ncaLzReadoutFn';
/** The subject the context panel currently shows — read by the panel's Reset-all control so it
 *  acts on the displayed subject. Set on every `showLzReadout`. */
const LZ_READOUT_SUBJECT_TAG = 'ncaLzReadoutSubject';
/** The F5 per-point decisions GRID viewer + its cell-click subscription — disposed + rebuilt per
 *  gesture (a build-once context-panel viewer goes blank after close→reopen, like the chart). */
const LZ_F5_GRID_TAG = 'ncaLzF5Grid';
const LZ_F5_CLICK_TAG = 'ncaLzF5ClickSub';
/** The subject the F5 grid last rendered — so a SAME-SUBJECT gesture (a toggle) REUSES the grid in
 *  place (no flicker, scroll kept) instead of rebuilding it. */
const LZ_F5_SUBJECT_TAG = 'ncaLzF5Subject';
/** The F5 State-glyph verdict hue, held in the bag so the (reused) onCellPrepare re-reads the
 *  current subject's hue after a verdict-changing toggle without re-registering the handler. */
const LZ_F5_HUE_TAG = 'ncaLzF5Hue';
/** F5 grid row height (px). ONE source for both the pinned `rowHeight` option and the host-height
 *  math (`rows * LZ_F5_ROW_HEIGHT + header`), so the two cannot drift (architect nit, 2026-09-03). */
const LZ_F5_ROW_HEIGHT = 24;
/** The context-panel's own single-subject profile chart (F#2). Held so it is `.close()`d and
 *  rebuilt per gesture — a build-once viewer in the context panel goes blank after the panel is
 *  closed & reopened (the same platform behaviour the study-overview panel works around). */
const LZ_READOUT_CHART_TAG = 'ncaLzReadoutChart';
/** The context-panel profile chart's own point-click subscription — a DISTINCT tag from the
 *  editor's `EDITOR_CLICK_SUB_TAG` so the two plots' subscriptions never clobber each other. */
const LZ_READOUT_CLICK_SUB_TAG = 'ncaLzReadoutClickSub';
/** LocalStorage key remembering whether the λz context-panel's profile-chart pane is EXPANDED:
 *  collapsed by default, but the user's toggle is remembered — within the session by the
 *  build-once accordion, and across sessions by this pref (user-directed 2026-09-02). Read/write
 *  are try/catch-guarded (localStorage can throw / be unavailable in some contexts). */
const LZ_CHART_EXPANDED_KEY = 'nca.lz.chartExpanded';
function readChartExpandedPref(): boolean {
  try {return localStorage.getItem(LZ_CHART_EXPANDED_KEY) === 'true';} catch {return false;}
}
function writeChartExpandedPref(expanded: boolean): void {
  try {localStorage.setItem(LZ_CHART_EXPANDED_KEY, expanded ? 'true' : 'false');} catch {/* unavailable */}
}
/** The chart context threaded into {@link showLzReadout} so its profile pane can rebuild the
 *  single-subject plot (the readout closure captures these from {@link ensureLambdaZEditor}). */
interface LzChartCtx {
  denseDf: DG.DataFrame; subjectName: string; timeName: string; concName: string; useLogY: boolean;
}
const EDITOR_REFRESH_TAG = 'ncaLzEditorRefresh';
const EDITOR_ROWSUB_TAG = 'ncaLzEditorRowSub';
const EDITOR_CLICK_SUB_TAG = 'ncaLzEditorClickSub';
/** `view.temp` key for the rxjs subscription that bridges
 *  `results.selection` → per-facet-plot opacity (see
 *  {@link wireResultsSelectionToFacetDim}). Replaced on every Compute. */
const FACET_SELECTION_SUB_TAG = 'ncaFacetSelectionSub';
/** `view.temp` key for the rxjs subscription that bridges the By-group aggregate
 *  grid's `selection` → the Mean±SD chart's curve filter (see
 *  {@link wireByGroupSelectionToMeanSdFilter}). Re-wired on every applyGrouping
 *  (the agg frame is rebuilt per Group-by change / Compute). */
const BYGROUP_SELECTION_SUB_TAG = 'ncaByGroupSelectionSub';
/** `view.temp` key for the current By-group→Mean±SD filter `apply` closure. Stored
 *  so a Mean±SD DataFrame SWAP that bypasses the selection subscription (the Semi-log
 *  toggle rebuilds a per-scale frame whose filter defaults to all-true) can re-apply
 *  the live selection filter to the new frame — otherwise the chart silently reverts
 *  to all curves while the grid stays selected. */
const BYGROUP_APPLY_TAG = 'ncaByGroupApply';

/** Maps a compute-notice's logical anchor to the `view.temp` grid/viewer tag of
 *  its home panel. The notice banner overlays that panel's `root`; when the panel
 *  is not docked, {@link resolveNoticeRoot} falls back to the results grid (B3). */
const NOTICE_GRID_TAG_BY_ANCHOR: Record<NoticeAnchor, string> = {
  byGroup: AGG_GRID_TAG,
  bioavailability: FPCT_GRID_TAG,
  doseProp: DP_GRID_TAG,
  accumulation: ACC_GRID_TAG,
  exposureMargin: MARGIN_GRID_TAG,
  sexComparison: SEXFOLD_GRID_TAG,
  sparse: SPARSE_PARAM_GRID_TAG,
  composite: SPARSE_COMP_GRID_TAG,
  meanSd: MEAN_SD_CHART_TAG,
  results: RESULTS_GRID_TAG,
};

/** Resolve a notice anchor to the DOM root to overlay: the home panel's
 *  grid/viewer root when docked, else the always-present results grid root (B3),
 *  else `null` (nothing docked yet). Reads dart-side `view.temp` with `!= null`
 *  (never-set keys return null, not undefined — BUG-06 footgun). */
function resolveNoticeRoot(view: DG.TableView, anchor: NoticeAnchor): HTMLElement | null {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const home = bag[NOTICE_GRID_TAG_BY_ANCHOR[anchor]] as DG.Viewer | null;
  if (home != null && home.root != null) return home.root;
  const results = bag[RESULTS_GRID_TAG] as DG.Viewer | null;
  return (results != null && results.root != null) ? results.root : null;
}

/**
 * AC-7.2 — live-replace the notice banner of each dependent-aggregate panel after an
 * interactive λz override, grouped by anchor → the panel's OWN docked grid root (the
 * SAME `refreshNotices` mechanism By-group uses on every Group-by change). `collected`
 * is that override's freshly-evaluated %F/dose-prop/TK notices ({@link collectNoticesFrom}).
 *
 * Keyed on each panel's OWN grid root (via {@link NOTICE_GRID_TAG_BY_ANCHOR}), NOT the
 * B3 results-fallback root, so a panel with a resolved warning is CLEARED (empty filter →
 * `refreshNotices` removes its banner) and no fallback-to-results banner is ever clobbered
 * (`refreshNotices` only touches banners parented to the root it is given). A panel that
 * is not docked has a null root and is skipped. By-group is NOT here — `applyGrouping`
 * already refreshes its own banner.
 */
function refreshDependentPanelBanners(
  view: DG.TableView, collected: ReadonlyArray<ComputeNotice>,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  // The dense results-derived decorators the AC-7.2 refresh re-runs (bioavailability /
  // dose-proportionality / toxicokinetics); each contributes notices under one anchor.
  const anchors: NoticeAnchor[] = ['bioavailability', 'doseProp', 'accumulation'];
  for (const anchor of anchors) {
    const home = bag[NOTICE_GRID_TAG_BY_ANCHOR[anchor]] as DG.Viewer | null;
    const root = (home != null && home.root != null) ? home.root : null;
    if (root === null) continue; // panel not docked → nothing to refresh
    refreshNotices(view, collected.filter((n) => n.anchor === anchor), root);
  }
}

/**
 * The facet viewport: how many plots stand side by side, and how many rows fit
 * on screen before the grid scrolls.
 *
 * **The facet grid is ALWAYS built** (user-directed 2026-08-13). It previously
 * refused above 20 subjects, on the DiffStudio rationale that a smaller tile is
 * unreadable on a 1920×1080 laptop — but refusing to draw is not the only answer
 * to "too small to read": keeping the tile at its 5×4 size and letting the grid
 * scroll gives the same readable tile AND the remaining subjects. The suppression
 * toast is gone with it.
 *
 * The cap's SECOND justification was cost, and dropping it needed a number rather
 * than a shrug. Measured on this host (`grok test --host local`, one
 * `Viewer.lineChart` per subject): **20 → 62 ms · 60 → 137 ms · 120 → 292 ms** —
 * linear, ≈ 2.4 ms per viewer, so a 200-subject study costs about half a second
 * once. That is a pause, not a hang, which is what makes "always build it"
 * affordable at any N a preclinical study reaches.
 */
const FACET_VIEWPORT_COLS = 5;
const FACET_VIEWPORT_ROWS = 4;

/**
 * Floor on a facet tile's height. **`FACET_VIEWPORT_ROWS` is a target, not a
 * mandate** — a row is `max(paneHeight / ROWS, this)`.
 *
 * Measured, after shipping the mandate and being shown the result: on a 686 px
 * facet pane, four rows meant 171 px per tile and DG's line chart dropped its
 * axes entirely — a strip of curve with no scale, which is worse than scrolling.
 * A plot left to size itself takes **300 px** (probed on this host), so that is
 * the floor: the height the viewer itself considers enough.
 *
 * The consequence is deliberate. On a typical ~700–1000 px pane the grid shows
 * two to three rows and scrolls for the rest, rather than four unreadable ones;
 * on a tall pane rows grow past 300 px and four fit exactly, which is the 5×4
 * base as asked. Readability wins the tie because an unreadable tile answers no
 * question — the same reasoning that retired the old hard cut-off, applied to
 * the other axis.
 */
export const FACET_TILE_MIN_PX = 300;

/** Above this many plots the hand-tuned split layout gives way to the scrolling
 *  `FACET_VIEWPORT_COLS`-wide grid. **User-directed = 10** (this cycle): split
 *  covers N ≤ 10 (see {@link getMaxFacetCols} for the per-N shapes), everything
 *  larger scrolls. Decoupled from the 5×4 scrolling viewport on purpose — that
 *  viewport still governs the scrolling tile HEIGHT, not this split/scroll cutoff. */
const FACET_SCROLL_ABOVE = 10;

/** Category-column name carried on every per-subject facet frame so a scatterPlot
 *  tile can wear its verdict hue (or the neutral fallback). scatterPlot has no
 *  `lineColor`/`markerColor` (unlike the retired lineChart) — colour rides on a
 *  column via `colorColumnName`/`linesByColumnName`, exactly like the Mean±SD
 *  scatterPlot. The per-tile colour depends on this column being CATEGORICAL, not
 *  numeric — see {@link paintFacetSeries} for the `subj '<subject>'` token that
 *  guarantees it (bare numeric-looking values made DG infer a numeric column, and
 *  categorical colour-coding silently did not apply — every tile went default
 *  blue; manual QA 2026-08-26). */
const FACET_COLOR_COL = 'color';

/** Structural type for an rxjs `Subscription` handle. Defined locally so the
 *  package does not pull in a direct rxjs import for one type. */
type Unsubscribable = {unsubscribe(): void};

export interface MountMappingPaneOptions {
  /** Overrides `DEFAULT_RULES` for tests / future presets. */
  readonly rules?: NcaRules;
  /**
   * AD-5 replay: pre-fill EVERY parameter-affecting mapping input from a saved
   * recipe so a subsequent Compute reproduces byte-identical parameters. Applied
   * after autodetect (overrides it) and seeds the `rules` closure from
   * `recipe.ncaRules` (R1 SY-F5 — pass `{recipe, rules: recipe.ncaRules}`).
   */
  readonly recipe?: AnalysisRecipe;
  /**
   * Lazy-hydration summary hook (AC-D1): fired after each successful Compute with
   * the per-profile results + the live recipe, so the coordinator can write the
   * summary cache (decision h). Never fired on validation-failed. The math is
   * NEVER read back from this cache (thesis #1) — recipe in, results recomputed.
   */
  readonly onComputed?: (ctx: ComputeCallbackCtx) => void;
  /**
   * Provenance of `df` (threaded from the launch path) — used by "Save to project"
   * to record where the analysis data lives (AC-C4 populate flow). Null/absent →
   * a best-effort in-session tableId, which the share dialog flags as owner-only.
   */
  readonly source?: SourcePointer;
  /**
   * When the analysis was launched from a study's "+ add analysis" entry, the
   * "Add to study" dialog pre-selects this study (the user can still switch).
   */
  readonly targetStudyId?: string;
}


/** Payload of {@link MountMappingPaneOptions.onComputed}. */
export interface ComputeCallbackCtx {
  readonly perProfile: ReadonlyArray<PerProfileResult>;
  readonly recipe: AnalysisRecipe;
  /** The post-`appendGroupingColumns` results DataFrame (carries the grouping
   *  metadata columns) — the per-(compound × route × dose-norm) substrate
   *  `buildSummaryCache` aggregates over (v2, study-unified-grid AD-1). */
  readonly resultsDf: DG.DataFrame;
  /** The observations substrate for the inline curve thumbnail (v5,
   *  chemist-sar-grid S4). `null` when no dense frame is available (defensive —
   *  the compute path always has one, so this is populated in practice). */
  readonly observations: ObservationsCtx | null;
  /** True when the view holds unsaved changes at compute time (nca-workflow-stages AD-8): the
   *  study-backed writeback then skips the cache, so the hub never shows numbers / a stage for
   *  a state that was not saved (D1). Filled by the pane, so the callback factory never reads
   *  the dirty tag across the views/ import cycle. */
  readonly dirty: boolean;
  /** The v15 workflow state this compute leaves (confirmed recipe, keyed overrides, reconciled
   *  review ledgers, derived stage) — persisted together with the cache it describes (AD-8). */
  readonly workflow: WorkflowPersist;
}

/**
 * Imperative handle returned by {@link mountMappingPane} (AD-5). `getCurrentRecipe`
 * serializes the live form into an exhaustive {@link AnalysisRecipe};
 * `triggerCompute` fires the SAME code path as the play icon — used by the
 * cold-deep-link hydration to compute without a user click (AC-D1: click and
 * deep-link share one path).
 */
export interface MappingPaneHandle {
  getCurrentRecipe(): AnalysisRecipe;
  /** Run the pipeline. `opts.silent` (set by the replay-on-open hydration path)
   *  suppresses the interactive confirm dialogs (per-kg; route-confirm via the a1
   *  auto-resolve) and the transient success-path toasts (composite-detected /
   *  sparse-design) — the user already saw them when the analysis was first
   *  computed and saved. A user-initiated Compute (play icon, or a manual re-run
   *  after open) passes no opts → everything shows. Audit records (WARNINGS_JSON)
   *  and persistent in-pane notice banners are unaffected. */
  triggerCompute(opts?: {silent?: boolean}): void;
}

/**
 * Mount the column-mapping form as a left-docked pane on the given table
 * view. Compute → `runNcaPipeline` → render validation errors *or* dock /
 * swap the results grid below the main grid in the right 70% area.
 */
export function mountMappingPane(
  view: DG.TableView,
  df: DG.DataFrame,
  opts: MountMappingPaneOptions = {},
): MappingPaneHandle {
  // R1 SY-F5: seed the `rules` closure from the replay recipe when hydrating, so
  // `buildRules()`'s closure-fallbacks (adjRSquaredFactor / extrapWarnPctAumc)
  // resolve to the replayed NcaRules — not DEFAULT_RULES.
  const rules = opts.rules ?? opts.recipe?.ncaRules ?? DEFAULT_RULES;

  const subjectInput = ui.input.column('Subject', {
    table: df, filter: isStringOrInt,
  });
  bindInputTip(subjectInput, richTip({
    title: 'Subject',
    summary: 'Column that identifies individual subjects.',
    detail: 'Each unique value defines a profile for NCA computation.',
  }));
  const timeInput = ui.input.column('Time', {
    table: df, filter: isNumeric,
  });
  bindInputTip(timeInput, richTip({
    title: 'Time',
    summary: 'Column with sampling time points.',
  }));
  const concInput = ui.input.column('Concentration', {
    table: df, filter: isNumeric,
  });
  bindInputTip(concInput, richTip({
    title: 'Concentration',
    summary: 'Column with measured concentrations.',
  }));

  // ── Event-driven (NONMEM) source mode (AD-1/AD-2) ─────────────────────────
  // Detect the NONMEM layout (EVID + AMT + a numeric EVID∈{0..4}); when present,
  // PROPOSE the mode via a toggle the scientist confirms (never silent, AC-A5).
  // The pane is re-mounted per upload, so detection + toggle state reset on a
  // dataset change (AC-A8). `confirmedRoute` is filled by the route-confirm
  // dialog (AD-5) when route must be inferred from RATE/CMT.
  const eventDetection: EventDrivenDetection = detectEventDrivenColumns(df);
  let confirmedRoute: RouteCode | null = null;
  const nonmemToggle = ui.input.toggle('Event-driven (NONMEM)', {value: false});
  bindInputTip(nonmemToggle, richTip({
    title: 'Event-driven (NONMEM)',
    summary: 'Treat the table as a NONMEM event record (EVID / AMT / CMT / RATE / ' +
      'MDV [+ ADDL / II]).',
    detail: 'Dose rows (EVID=1) are filtered out, observations are re-timed to ' +
      'time-after-dose, ADDL/II repeat schedules are expanded, and AMT becomes the ' +
      'absolute-mg dose of record. Proposed automatically when the layout is detected.',
  }));
  // Hint shown when the mode is proposed but not yet confirmed/applied.
  const nonmemHint = ui.divText(
    'NONMEM event layout detected — leave on to ingest EVID/AMT/ADDL/II, or turn off ' +
    'for standard column mapping.', 'nca-mapping-pane-sparse-hint');
  nonmemHint.style.display = 'none';

  // Source-mode radios gate dose / route / infusion / body weight / LLOQ between a
  // value and a column. The corresponding scalar / column inputs are shown or
  // hidden based on the choice (see `linkSourceVisibility` below).
  const doseSourceInput = sourceModeInput('Dose type');
  bindInputTip(doseSourceInput, richTip({
    title: 'Dose type',
    summary: 'How dose is provided: a single value applied to all subjects, or a ' +
      'per-row column from the table.',
  }));
  const doseColInput = ui.input.column('Dose', {
    table: df, filter: isNumeric, nullable: true,
  });
  bindInputTip(doseColInput, richTip({
    title: 'Dose (column)',
    summary: 'Per-row dose column.',
  }));
  const doseScalarInput = ui.input.float('Dose', {
    value: 1,
  });
  bindInputTip(doseScalarInput, richTip({
    title: 'Dose',
    summary: 'Scalar dose applied to all subjects.',
  }));
  const routeSourceInput = sourceModeInput('Route type');
  bindInputTip(routeSourceInput, richTip({
    title: 'Route type',
    summary: 'How route is provided: a single route applied to all subjects, or a ' +
      'per-row column for mixed-route analysis.',
  }));
  const routeScalarInput = ui.input.choice<string>('Route', {
    value: 'PO', items: [...ROUTE_CHOICES],
  });
  bindInputTip(routeScalarInput, richTip({
    title: 'Route',
    summary: 'Administration route applied to all subjects.',
    detail: 'Determines CL vs CL/F and Vz vs Vz/F labeling in the results.',
    def: routeScalarInput.value as string, // single-source: the input's own default
  }));
  // Per-row route COLUMN picker (per-row-route AD-1). String or int; nullable.
  // When the Route type toggle is 'column', its values are tolerantly
  // normalized to RouteCode in the Compute handler. (Pattern A: the toggle
  // makes scalar/column mutually exclusive, so ROUTE_SCALAR_COLUMN_CONFLICT is
  // unreachable from this form — that code stays live only for programmatic
  // callers.)
  const routeColInput = ui.input.column('Route', {
    table: df, filter: isStringOrInt, nullable: true,
  });
  bindInputTip(routeColInput, richTip({
    title: 'Route (column)',
    summary: 'Per-row administration-route column (e.g. IV / PO / SC).',
    detail: 'Enables mixed-route analysis: each subject is labeled CL or CL/F by its own route.',
  }));
  const infusionSourceInput = sourceModeInput('Infusion duration type');
  bindInputTip(infusionSourceInput, richTip({
    title: 'Infusion duration type',
    summary: 'How infusion duration is provided: a single value or a per-row column.',
    detail: 'Only relevant for the IV-infusion route.',
  }));
  // Optional per-row infusion-duration column — used ONLY to disambiguate bare
  // "IV" into IV-bolus vs IV-infusion per row (AD-3 / AC-A8). sci-comp v0.6 does
  // not otherwise consume per-row infusion duration (it is inert), so this
  // column affects route resolution, not the AUC math.
  const infusionDurationColInput = ui.input.column('Infusion duration', {
    table: df, filter: isNumeric, nullable: true,
  });
  bindInputTip(infusionDurationColInput, richTip({
    title: 'Infusion duration (column)',
    summary: 'Optional per-row infusion duration.',
    detail: 'A positive value resolves a bare "IV" route to IV-infusion (else IV-bolus). Units = the Time column.',
  }));
  // Scalar counterpart of the infusion-duration column (kept contiguous with its
  // type toggle + column — Pattern A). Visibility is governed by
  // `refreshInfusionVisibility` (route relevance × source toggle), installed below.
  const infusionDurationInput = ui.input.float('Infusion duration', {
    value: 0,
  });
  bindInputTip(infusionDurationInput, richTip({
    title: 'Infusion duration',
    summary: 'Duration of the IV infusion, in the units of the Time column.',
    detail: 'Required for IV-infusion route; ignored otherwise.',
    def: `${infusionDurationInput.value} (unspecified)`,
  }));
  const bodyWeightSourceInput = sourceModeInput('Body weight type');
  bindInputTip(bodyWeightSourceInput, richTip({
    title: 'Body weight type',
    summary: 'How body weight is provided: a single value applied to all subjects, or a ' +
      'per-subject column.',
  }));
  const bodyWeightColInput = ui.input.column('Body weight', {
    table: df, filter: isNumeric, nullable: true,
  });
  bindInputTip(bodyWeightColInput, richTip({
    title: 'Body weight (column)',
    summary: 'Per-subject body weight column.',
    detail: 'When provided, CL and Vz are weight-normalized.',
  }));
  const blqFlagColInput = ui.input.column('BLQ flag column', {
    table: df, filter: isBoolLikeNumeric, nullable: true,
  });
  bindInputTip(blqFlagColInput, richTip({
    title: 'BLQ flag column',
    summary: 'Below-Limit-of-Quantification (BLQ) flag column.',
    detail: 'Rows marked TRUE are treated as below-limit measurements.',
  }));
  const lloqSourceInput = sourceModeInput('LLOQ type');
  bindInputTip(lloqSourceInput, richTip({
    title: 'LLOQ type',
    summary: 'How the Lower Limit of Quantification (LLOQ) is provided: a single value for ' +
      'all measurements, or a per-row column from the table.',
  }));
  const lloqColInput = ui.input.column('LLOQ', {
    table: df, filter: isNumeric, nullable: true,
  });
  bindInputTip(lloqColInput, richTip({
    title: 'LLOQ (column)',
    summary: 'Column with Lower Limit of Quantification (LLOQ) values per row.',
    // LLOQ has no unit picker of its own — it is expressed in the Concentration
    // unit and rescaled by the same factor (`CONC_TO_MG_L`, unit-conversion.ts).
    // Stated HERE, on the field, because the field is where the scientist types
    // the number: entering 0.5 against a ng/mL (or µg/L) concentration column
    // while reading it as mg/L is a 1000x error in which samples classify as BLQ.
    // (mg/L and µg/mL are NOT such a pair — both are factor 1, numerically
    // identical; the real 1000x steps are mg/L vs µg/L and mg/L vs ng/mL.)
    // Since the inline-unit rows landed, LLOQ is the only numeric field in
    // Context with no unit beside it, so its own tooltip has to carry what the
    // row no longer shows.
    detail: 'Expressed in the Concentration unit — rescaled by the same factor.',
  }));
  const lloqScalarInput = ui.input.float('LLOQ', {
    value: 0,
  });
  bindInputTip(lloqScalarInput, richTip({
    title: 'LLOQ',
    summary: 'Lower Limit of Quantification (LLOQ) applied uniformly to all measurements.',
    // See the LLOQ (column) tooltip above — same reason, same wording.
    detail: 'Expressed in the Concentration unit — rescaled by the same factor.',
  }));

  // ---- Per-profile scalars previously hardcoded to null in ProfileScalars.
  // Body weight is EITHER this uniform scalar OR a per-subject column, selected
  // by `bodyWeightSourceInput` (Pattern A — exactly one is visible). The
  // extractor-level column precedence (profile-inputs.ts:160) remains as a
  // defensive fallback but is no longer reachable from this form. Value 0 /
  // null → treated as unspecified (sci-comp gets null, weight-normalised
  // params are skipped).
  const bodyWeightScalarInput = ui.input.float('Body weight', {
    value: 0,
  });
  bindInputTip(bodyWeightScalarInput, richTip({
    title: 'Body weight',
    summary: 'Body weight value applied when no Body weight column is mapped.',
    detail: 'Used for weight-normalized CL/Vz parameters.',
    def: `${bodyWeightScalarInput.value} (unspecified)`,
  }));
  // ---- Subject Sex (sex-stratifier designation, Slice A/B). The Sex role was
  // previously bound to one hardcoded `/^sex$/i` regex at four sites with no
  // override, so a non-`sex`-named or coded-`1/2` column silently lost F% /
  // dose-proportionality / TK stratification. This picker DESIGNATES the sex
  // column; a synonym-matching column is AUTO-BOUND on mount (see the auto-bind
  // block after the section is assembled — user-directed, supersedes the earlier
  // suggest-only A2 design), and the value→M/F/U mapping inputs then surface for
  // review. Every reference NCA tool (Phoenix sort-vars, PKNCA grouping, Pumas
  // `group=`) lets the user pick the stratifier; this mirrors that. The picker is
  // NOT added to ColumnContract — the designated column stays carried metadata
  // (appendGroupingColumns) the decorators read by name.
  const sexColInput = ui.input.column('Sex column', {
    table: df, filter: isStringOrInt, nullable: true,
  });
  bindInputTip(sexColInput, richTip({
    title: 'Sex column',
    summary: 'Pick the per-subject column that encodes sex.',
    detail: 'Leave empty if your data has no sex axis — sex-stratified parameters ' +
      '(%F, dose-proportionality, TK) will be skipped. Map each value to M/F/U below.',
  }));
  // `>8`-distinct guard surfaced as an INLINE DG validator (red under the picker),
  // not a panel text line — keeps the section terse (user-directed, supersedes the
  // earlier architect "no DG validator" finding). The `sexTooMany` Compute-disable
  // flag is computed separately in `refreshSexMapping` (the validator is the visual;
  // the flag is the gate). The validator reads picker state, so its `value` arg is
  // unused.
  sexColInput.addValidator(() => {
    const c = sexColInput.value;
    if (c === null) return null;
    const n = uniqueSexValues(c).length;
    return n > SEX_MAX_DISTINCT ?
      `${n} distinct values — too many for a sex axis; pick a column with ≤${SEX_MAX_DISTINCT} categories` :
      null;
  });
  /** Single source of truth for the designated sex column (A3). Reads picker
   *  state — NOT a name heuristic. Replaces the four inline `/^sex$/i` reads. */
  const getSexColumnName = (): string | null => sexColInput.value?.name ?? null;
  // Container for the dynamic per-value → M/F/U mapping inputs (Slice B fills it
  // on column designation; empty in the picker-only state).
  const sexMapContainer = ui.div([], 'nca-sex-map-container');
  // Slice B state: one choice input per distinct raw sex value, plus the two
  // Compute-disable flags `refreshComputeEnabled` reads. `sexTooMany` (column has
  // > SEX_MAX_DISTINCT values — implausible as a sex axis) and `sexUnassigned`
  // (a value still mapped to the "—" sentinel) each BLOCK Compute, so F% is never
  // computed on an unresolved / mis-picked stratification (honest uncertainty).
  const sexValueInputs: Array<{raw: string; input: DG.InputBase<string | null>}> = [];
  // onChanged subscriptions on the per-value inputs — torn down before each
  // rebuild so re-selecting a column does not orphan stale subscriptions.
  const sexValueSubs: Array<{unsubscribe(): void}> = [];
  let sexTooMany = false;
  let sexUnassigned = false;
  // AD-5 replay: a recipe's raw→M/F/U map, staged here so `refreshSexMapping`
  // (fired by setting `sexColInput.value`) pre-fills each per-value input from the
  // recipe instead of the name heuristic. Null in the normal (non-replay) path.
  let pendingSexMap: Readonly<Record<string, string>> | null = null;
  // Last successful compute's per-profile results — kept so "Save to project" can
  // build the summary cache without re-running the math (thesis #1). Null until
  // the first successful Compute.
  let lastPerProfile: ReadonlyArray<PerProfileResult> | null = null;
  // The matching post-appendGroupingColumns results frame for the last successful
  // Compute — kept alongside `lastPerProfile` so "Add to study" hands the grouping
  // substrate to `buildSummaryCache` (v2, study-unified-grid AD-1). Set/cleared in
  // lockstep with `lastPerProfile`.
  let lastResults: DG.DataFrame | null = null;
  // The matching observations substrate for the inline curve thumbnail (v5,
  // chemist-sar-grid S4) — the same normalized dense frame + role columns the
  // mean±SD chart pools from. Set/cleared in lockstep with `lastResults` so
  // "Add to study" hands it to `buildSummaryCache` (the study overview is
  // cache-only, so the thumbnail curve must be built here, at compute time).
  let lastObservations: ObservationsCtx | null = null;
  // The post-compute summary-cache hook. Starts as opts.onComputed (set by the reopen path), and is
  // ALSO installed by promote-after-save (F4/D7) so a just-added analysis caches on recompute like a
  // reopened one — via the shared makeStudyBackedOnComputed factory. Held in a `let` (not read from
  // the readonly opts) precisely so promote can install it after mount.
  let studyOnComputed = opts.onComputed;
  // Replay-open: the persisted group-by column(s) to restore once the first Compute
  // populates the Group-by candidates (applyRecipe runs before that). Consumed
  // (set null) on first restore; the live user selection wins thereafter.
  let pendingGroupByCols: readonly string[] | null = null;
  // The live multi-column Group-by selection as RESULTS-column NAMES (the source of
  // truth for both consumers; the picker's `DG.Column[]` value binds to a per-Compute
  // frame, so we track names, v8 multi-column-group-by). Empty = no grouping.
  let groupBySelectedNames: string[] = [];
  // Guards the picker's onValueChanged while we set its value PROGRAMMATICALLY (rebind
  // on Compute / replay restore) so the compute path's own applyGrouping/chart refresh
  // isn't double-fired mid-rebuild.
  let suppressGroupByChange = false;
  /** Build the resolved raw→M/F/U map from the per-value inputs (skips the
   *  unassigned sentinel — the Compute-disable gate guarantees none remain). */
  const buildSexMap = (): Map<string, SexCode> => {
    const m = new Map<string, SexCode>();
    for (const {raw, input} of sexValueInputs) {
      const v = input.value;
      if (v === 'M' || v === 'F' || v === 'U') m.set(raw, v);
    }
    return m;
  };
  // ---- Sparse / destructive design (UC-04). Three role mappings the sparse
  // pipeline previously resolved ONLY by name heuristic (`findColByPattern` in
  // `applySparseNca`): a non-conventionally-named sparse dataset (e.g.
  // `ScheduledTime` / `AnimalNo`) got no Sparse panels and no explanation. These
  // explicit controls make the roles mappable + the topology declarable. The
  // toggle gates them: OFF ⇒ today's pure-heuristic path (zero regression); ON ⇒
  // the picked values override the heuristic. Pickers PRE-FILL from the heuristic
  // and the toggle AUTO-CHECKS when a nominal-time column is detected (see the
  // autodetect block), so conventionally-named data still "just works".
  const sparseToggleInput = ui.input.toggle('Sparse design', {
    value: false,
  });
  bindInputTip(sparseToggleInput, richTip({
    title: 'Sparse design',
    summary: 'Enable explicit sparse / destructive-sampling design mapping.',
    detail: 'Auto-enabled when a nominal-time column is detected. When off, sparse ' +
      'roles are resolved by column-name heuristic (the prior behaviour).',
  }));
  const nominalTimeColInput = ui.input.column('Nominal time', {
    table: df, filter: isNumeric, nullable: true,
  });
  bindInputTip(nominalTimeColInput, richTip({
    title: 'Nominal time',
    summary: 'The nominal (scheduled) sampling-time column the sparse composite mean profile is pooled on.',
    detail: 'This is the trigger for sparse NCA; distinct from the per-animal observed Time ' +
      'used by the dense pipeline.',
  }));
  const animalIdColInput = ui.input.column('Animal ID', {
    table: df, filter: isStringOrInt, nullable: true,
  });
  bindInputTip(animalIdColInput, richTip({
    title: 'Animal ID',
    summary: 'Animal / subject-pool identifier used to derive the destructive-vs-batch ' +
      'sampling topology.',
    def: 'Subject column',
  }));
  // Declared sampling topology. `(auto)` (default) keeps today's behaviour: the
  // topology is derived from the data (and cross-checked against a detected
  // SamplingType column). A non-auto choice DECLARES the topology study-wide; it
  // is routed to sci-comp `sparseAuc`'s AD-2 cross-check, so a declaration that
  // disagrees with the animal × nominal-time data surfaces SPARSE_TOPOLOGY_MISMATCH.
  const SAMPLING_TYPE_AUTO = '(auto)';
  const samplingTypeInput = ui.input.choice<string>('Sampling type', {
    value: SAMPLING_TYPE_AUTO,
    items: [SAMPLING_TYPE_AUTO, 'destructive', 'batch', 'serial'],
  });
  bindInputTip(samplingTypeInput, richTip({
    title: 'Sampling type',
    summary: 'Declared sampling topology.',
    detail: '(auto) derives it from the data. destructive / batch / serial declares it and is ' +
      'cross-checked against the data — a mismatch is flagged, not silently overridden.',
    def: SAMPLING_TYPE_AUTO, // single-source: the same constant the input is initialized with
  }));
  // Guidance hint (AC3): shown ONLY when the sparse toggle is on but no
  // nominal-time column is mapped — so a sparse intent never renders silently.
  // Conservative default for an inline form hint: neutral muted text, smallest size.
  const sparseHint = ui.divText(
    'Map a nominal-time column to run sparse NCA.', 'nca-mapping-pane-sparse-hint');
  sparseHint.style.display = 'none';

  // ---- Method settings (previously hardcoded in DEFAULT_RULES). Initial
  // values come from the supplied `rules` (which is `opts.rules ??
  // DEFAULT_RULES`); the click handler reads the live input values via
  // `buildRules()` before each Compute, so user edits flow through.
  const aucMethodInput = ui.input.choice<string>('AUC method', {
    value: rules.aucMethod,
    items: [...AUC_METHOD_CHOICES],
  });
  bindInputTip(aucMethodInput, richTip({
    title: 'AUC method',
    summary: 'Trapezoidal AUC integration scheme.',
    detail: 'linear-up-log-down uses linear interpolation between samples where concentration ' +
      'is rising or constant (absorption / distribution) and log-linear between samples where ' +
      'it is falling (any declining interval).',
  }));
  const compensatedInput = ui.input.bool('Compensated summation', {
    value: rules.compensatedSummation,
  });
  bindInputTip(compensatedInput, richTip({
    title: 'Compensated summation',
    summary: 'Use compensated summation in AUC accumulators — reduces floating-point drift ' +
      'on long profiles with values spanning many orders of magnitude.',
    detail: 'Slower; rarely changes the answer.',
  }));

  const blqPreInput = ui.input.choice<string>('BLQ before first', {
    value: rules.blq.preFirstMeasurable, items: [...BLQ_RULE_CHOICES],
  });
  bindInputTip(blqPreInput, richTip({
    title: 'BLQ before first',
    summary: 'Action on BLQ points before the first measurable observation.',
  }));
  const blqEmbeddedInput = ui.input.choice<string>('BLQ embedded', {
    value: rules.blq.embedded, items: [...BLQ_RULE_CHOICES],
  });
  bindInputTip(blqEmbeddedInput, richTip({
    title: 'BLQ embedded',
    summary: 'Action on BLQ points between two measurable observations.',
  }));
  const blqAfterInput = ui.input.choice<string>('BLQ first after last', {
    value: rules.blq.afterLast, items: [...BLQ_RULE_CHOICES],
  });
  bindInputTip(blqAfterInput, richTip({
    title: 'BLQ first after last',
    summary: 'Action on the first BLQ point after the last measurable observation.',
  }));
  const blqConsecutiveInput = ui.input.choice<string>('BLQ trailing', {
    value: rules.blq.consecutiveAfterLast, items: [...BLQ_RULE_CHOICES],
  });
  bindInputTip(blqConsecutiveInput, richTip({
    title: 'BLQ trailing',
    summary: 'Action on the second and subsequent consecutive BLQ points after the last ' +
      'measurable observation.',
  }));

  // λz selection mode — only auto-best-fit is reachable in v1; the
  // manual-points / manual-time-range modes require an interactive
  // click-to-include UI on the chart (separate cycle). The input is kept
  // declared so {@link buildRules} can read a value, but it is NOT
  // appended to the form: a single-option dropdown is a dead control
  // (changing nothing). When manual modes ship, this hide will be
  // reverted and the input appended back to the form.
  const lambdaModeInput = ui.input.choice<string>('λz mode', {
    value: rules.lambdaZ.mode, items: [...LAMBDAZ_MODE_CHOICES],
  });
  const lambdaMinPointsInput = ui.input.int('λz min points', {
    value: rules.lambdaZ.minPoints, min: 2,
  });
  bindInputTip(lambdaMinPointsInput, richTip({
    title: 'λz min points',
    summary: 'Minimum number of terminal points required to accept a candidate λz fit.',
    detail: 'Lower values fit more profiles but make individual fits less defensible.',
  }));
  const lambdaMinRSquaredInput = ui.input.float('λz min adj-R²', {
    value: rules.lambdaZ.minRSquared, min: 0, max: 1,
  });
  bindInputTip(lambdaMinRSquaredInput, richTip({
    title: 'λz min adj-R²',
    summary: 'Minimum adjusted R² for a candidate window during λz POINT SELECTION.',
    detail: 'The sci-comp technical gate: fits below this are not considered by the ' +
      'core. This is NOT the acceptance threshold — that is a separate field below. ' +
      'Keep it soft; the acceptance judgment lives in "adj-R² acceptance threshold".',
  }));
  // adj-R² ACCEPTANCE threshold (R_acc) — lambdaz-adjr2-threshold-spec-v2 §3.2. A
  // SEPARATE field from the point-selection gate above: the value entered here is applied
  // VERBATIM as the acceptance criterion (adjR² < R_acc → Flagged), with no derived
  // recomputation. Not a member of NcaRules (SAP/presentation, thesis #1); it feeds the
  // results-view thresholds directly (see the thresholds projection at Compute).
  const adjRSqAcceptanceInput = ui.input.float('adj-R² acceptance threshold', {
    // Restored from the replayed recipe on open (persistence); else the default.
    value: opts.recipe?.acceptanceRSquared ?? DEFAULT_ACCEPTANCE_RSQUARED,
    min: ACCEPTANCE_RSQUARED_MIN, max: ACCEPTANCE_RSQUARED_MAX,
  });
  bindInputTip(adjRSqAcceptanceInput, richTip({
    title: 'adj-R² acceptance threshold',
    summary: 'Acceptance bar for the finished λz fit.',
    detail: 'adj. R² < this → profile Flagged. Applied as entered (not derived). ' +
      'Must be ≥ λz min adj-R².',
  }));
  // Consistency constraint (§3.2): R_acc ≥ minR — the gate cannot be softer than the
  // acceptance bar. Uses DG's own validator surface (inline red error); re-validated when
  // the gate changes so the message tracks either edit.
  adjRSqAcceptanceInput.addValidator((): string | null => {
    const rAcc = adjRSqAcceptanceInput.value;
    const minR = lambdaMinRSquaredInput.value;
    return (rAcc != null && minR != null && rAcc < minR) ?
      `Must be ≥ λz min adj-R² (${minR}).` : null;
  });
  lambdaMinRSquaredInput.onChanged.subscribe(() => adjRSqAcceptanceInput.validate());
  const lambdaExcludeCmaxInput = ui.input.bool('Exclude Cmax from λz', {
    value: rules.lambdaZ.excludeCmax,
  });
  bindInputTip(lambdaExcludeCmaxInput, richTip({
    title: 'Exclude Cmax from λz',
    summary: 'When ON, the Cmax point is excluded from the λz candidate window.',
    detail: 'Standard practice for orally administered profiles where the absorption phase ' +
      'distorts the early terminal slope.',
  }));
  const lambdaAdjRSqFactorInput = ui.input.float('λz tie-break factor', {
    value: rules.lambdaZ.adjRSquaredFactor ?? 1e-4, min: 0,
  });
  bindInputTip(lambdaAdjRSqFactorInput, richTip({
    title: 'λz tie-break factor',
    summary: 'When two candidate fits have adj-R² within this distance, the one with more points wins.',
    formula: 'score = adj-R² + factor · n',
    detail: 'Set 0 to disable tie-breaking entirely. Changing it shifts λz selection on ' +
      'profiles whose candidate windows have comparable adj-R².',
    // single-source: same expression that initializes the input value
    def: `${rules.lambdaZ.adjRSquaredFactor ?? 1e-4} (reference-suite-validated)`,
  }));

  const extrapWarnPctInput = ui.input.int('Extrapolation warn %', {
    value: rules.extrapWarnPct, min: 0, max: 100,
  });
  bindInputTip(extrapWarnPctInput, richTip({
    title: 'Extrapolation warn %',
    summary: 'Soft warning threshold for % AUC extrapolated beyond t_last.',
    detail: 'Above this %, the AUCinf value gets a warning flag in the parameter warnings stream.',
  }));
  const extrapErrorPctInput = ui.input.int('Extrapolation error %', {
    value: rules.extrapErrorPct, min: 0, max: 100,
  });
  bindInputTip(extrapErrorPctInput, richTip({
    title: 'Extrapolation error %',
    summary: 'Hard error threshold for % AUC extrapolated beyond t_last.',
    detail: 'Above this %, AUCinf is considered unreliable.',
  }));

  // Unit inputs — load-bearing after `nca-unit-semantics` ships
  // (NOTE-07 closure). Selected unit strings drive adapter-side rescaling
  // to canonical h / mg/L / mg / kg before sci-comp's unit-agnostic
  // computeNca runs. Conversion lives in `adapter/unit-conversion.ts`
  // (AD-1, single rescaling point at the boundary).
  //
  // Default canonical units (mg/L / h / mg) — flow through `scaleProfileInputs`
  // as the identity factor when the user does not override.
  //
  // CAPTION IS DELIBERATELY EMPTY on all four unit pickers (mapping-pane-inline-units,
  // 2026-07-28). Each rides as the second control on its value's own form row —
  // `unitRow(timeInput, timeUnitInput)` below — so the value and the unit that scales
  // it are read together. A caption here would restate the value's own caption one
  // control to the left ("Dose  [1] Dose [mg]"), the redundancy the move exists to remove.
  // The control stays identifiable two ways:
  // the unit tokens are self-describing (h / mg/L / mg / kg), and each rich tooltip
  // carries its OWN title ("Time unit", …) rather than inheriting the caption.
  //
  // Measured consequence (live probe, grok test 2026-07-28): with caption '', DG still
  // builds a `captionLabel` but lays it out at 0px. `bindInputTip` binds BOTH
  // captionLabel and the editor (views/rich-tooltip.ts:76), so the surviving hover
  // target is the dropdown itself — the whole visible control. The rescale explanation
  // (the only in-form statement of the canonical h / mg/L / mg / kg conversion, and the
  // direct cover for the silent mis-scale behind GAP-11 / GAP-13) stays reachable.
  // ONE rich-tooltip factory per unit, defined once and shared by BOTH the horizontal
  // and the vertical (dual-mode) picker of that unit — `richTip(...)` returns a reusable
  // `() => HTMLElement` renderer, so binding the same factory to two inputs is safe and
  // the two pickers' tooltips can never drift (architect SIMPLIFY).
  const timeUnitInput = ui.input.choice<string>('', {
    value: 'h', items: [...TIME_UNIT_CHOICES],
  });
  const timeUnitTip = richTip({
    title: 'Time unit',
    summary: 'Unit for the Time column; values are rescaled to canonical h before computation.',
    detail: 'Also applies to any scalar Infusion duration value — ensure both are in the same unit.',
    def: timeUnitInput.value as string, // single-source: the input's own default
  });
  bindInputTip(timeUnitInput, timeUnitTip);
  const concUnitInput = ui.input.choice<string>('', {
    value: 'mg/L', items: [...CONC_UNIT_CHOICES],
  });
  const concUnitTip = richTip({
    title: 'Concentration unit',
    summary: 'Unit for the Concentration column; values are rescaled to canonical mg/L before computation.',
    detail: 'The same factor also rescales LLOQ — it is expressed in the concentration unit.',
    def: concUnitInput.value as string,
  });
  bindInputTip(concUnitInput, concUnitTip);
  const doseUnitInput = ui.input.choice<string>('', {
    value: 'mg', items: [...DOSE_UNIT_CHOICES],
  });
  const doseUnitTip = richTip({
    title: 'Dose unit',
    summary: 'Unit for the Dose column; values are rescaled to canonical mg before computation.',
    detail: 'Per-kg variants (mg/kg, µg/kg) are multiplied by body weight (kg).',
    def: doseUnitInput.value as string,
  });
  bindInputTip(doseUnitInput, doseUnitTip);

  // Body weight unit picker (AD-6 / E7). Defaults to kg; pick g for rodent
  // studies reporting grams. A species-absent heuristic warning fires when
  // bodyWeight > 15 with bwUnit='kg' (catches mouse-as-kg silent 1000×
  // error).
  const bodyWeightUnitInput = ui.input.choice<string>('', {
    value: 'kg', items: [...BODY_WEIGHT_UNIT_CHOICES],
  });
  const bodyWeightUnitTip = richTip({
    title: 'Body weight unit',
    summary: 'Unit for the Body weight column or scalar; values are rescaled to kg before computation.',
    detail: 'Pick g for rodent studies reporting grams.',
    def: bodyWeightUnitInput.value as string,
  });
  bindInputTip(bodyWeightUnitInput, bodyWeightUnitTip);
  // (The non-canonical concentration-unit rescale — and that LLOQ shares the same
  // factor — is surfaced by the Concentration unit input's own rich tooltip
  // above, so no separate on-form notice is needed.)

  // --- Dual-mode unit pickers (mapping-pane-unit-dual-mode) --------------------
  // The docked `ui.form` flips between a "wide" label-left layout and a "narrow"
  // label-on-top layout (the platform toggles the `ui-form-condensed` class on the
  // form by width). The inline `unitRow` (divH) picker reads correctly only in the
  // wide layout; in narrow it orphans past the value's right edge (the reported bug).
  // So each unit carries TWO synchronised pickers: the HORIZONTAL one (above, via
  // `unitRow`, on the value's own row) shown in WIDE mode, and a VERTICAL one (here,
  // its own captioned row placed under the value) shown in NARROW mode.
  // `applyUnitPickerMode` (after dock) toggles which is visible off the
  // `ui-form-condensed` class; `linkUnitPair` keeps the pair equal so every existing
  // consumer keeps reading the ORIGINAL (horizontal) input unchanged. A deliberate
  // crutch — the platform exposes no clean single-control seam for this.
  const timeUnitInputV = ui.input.choice<string>('Time unit',
    {value: timeUnitInput.value as string, items: [...TIME_UNIT_CHOICES]});
  bindInputTip(timeUnitInputV, timeUnitTip); // shared factory — cannot drift from horizontal
  const concUnitInputV = ui.input.choice<string>('Concentration unit',
    {value: concUnitInput.value as string, items: [...CONC_UNIT_CHOICES]});
  bindInputTip(concUnitInputV, concUnitTip);
  const doseUnitInputV = ui.input.choice<string>('Dose unit',
    {value: doseUnitInput.value as string, items: [...DOSE_UNIT_CHOICES]});
  bindInputTip(doseUnitInputV, doseUnitTip);
  const bodyWeightUnitInputV = ui.input.choice<string>('Body weight unit',
    {value: bodyWeightUnitInput.value as string, items: [...BODY_WEIGHT_UNIT_CHOICES]});
  bindInputTip(bodyWeightUnitInputV, bodyWeightUnitTip);
  // Keep each horizontal/vertical pair equal. The value-equality guard (set only when
  // different) is loop-safe whether DG fires `onChanged` synchronously or async — the
  // echo write finds the values already equal and stops. Consumers read the HORIZONTAL
  // (original) input; editing the vertical one in narrow mode mirrors back to it (and
  // so to `onAnyChange` + the recipe seam) via the horizontal input's own subscriptions.
  const linkUnitPair = (
    h: DG.InputBase<string | null>, v: DG.InputBase<string | null>,
  ): void => {
    h.onChanged.subscribe(() => {if (v.value !== h.value) v.value = h.value;});
    v.onChanged.subscribe(() => {if (h.value !== v.value) h.value = v.value;});
  };
  linkUnitPair(timeUnitInput, timeUnitInputV);
  linkUnitPair(concUnitInput, concUnitInputV);
  linkUnitPair(doseUnitInput, doseUnitInputV);
  linkUnitPair(bodyWeightUnitInput, bodyWeightUnitInputV);

  const errorPanel = ui.divV([], 'nca-mapping-pane-error-panel');
  errorPanel.style.display = 'none';

  // Compute trigger as a solid FA play icon (pattern: compute-utils
  // function-views/src/fitting-view.ts:602-613). Gating works through
  // (1) a click-time guard, (2) icon color, (3) tooltip. Two reasons the
  // icon may be inactive:
  //   - required columns not picked yet (Subject / Time / Concentration)
  //   - last Compute already ran successfully and no input has changed
  //     since ("consumed"; re-activates on any onChanged from any input).
  let consumed = false;
  /** Confirm is enabled (as last computed by `refreshComputeEnabled`) — the stage strip reads it to
   *  make Review a Confirm link on Setup (user-directed 2026-10-13). */
  let confirmReady = false;

  // Y-axis mode for both the Combined overlay chart and every Facet plot.
  // Default 'log' because terminal-phase reading is the PK convention. The
  // Semi-log toggle in the Settings dialog (opened by the ribbon gear) flips
  // this; pre-Compute changes only update the state (no charts exist yet),
  // post-Compute changes call applyYAxisModeToCharts to re-bind yColumnNames
  // on the live viewers.
  let useLogY = true;
  // Latching flag: the Settings gear (which opens the Semi-log + Compact dialog)
  // is hidden until the first successful Compute, then revealed by re-calling
  // setRibbonPanels. Both toggles only affect post-Compute surfaces (charts /
  // results grid), so the gear appears with them — same timing as Group by.
  // Subsequent Computes do NOT re-call (would replace the entire ribbon).
  let settingsRevealed = false;

  // --- Per-group aggregation (#11) -----------------------------------------
  // "No grouping" is the nullable choice's blank option (value === null) — the
  // per-subject results view only. No explicit sentinel item.
  // Latching flag mirroring `settingsRevealed`: the Group-by ribbon group is
  // hidden until the first Compute that yields ≥ 1 grouping candidate.
  let groupByRevealed = false;

  // Friendly display label for a results-column name (the same header the results
  // grid shows, set by `swapOrDockResults`): native keys read "Route" / "CL basis" /
  // "Status"; carried-metadata columns have no friendlyName and fall back to the name.
  const groupByFriendly = (resultsDf: DG.DataFrame, name: string): string =>
    resultsDf.col(name)?.getTag('friendlyName') ?? name;

  /** The LIVE results frame (the docked results grid's), or null before the first Compute.
   *  Same lookup `applyGrouping` does — the Group-by picker's value must always be
   *  re-resolved against THIS frame (see {@link dedupeGroupBySelection}). */
  const liveResultsDf = (): DG.DataFrame | null => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const tempBag = view.temp as Record<string, any>;
    // `== null` not `=== undefined`: view.temp returns null for an unset key (BUG-06).
    const grid = tempBag[RESULTS_GRID_TAG] as DG.Grid | null;
    return grid == null ? null : grid.dataFrame;
  };

  // Rebind the multi-column picker to the fresh Compute's results frame, scoping its
  // candidates to `candidates` (F1: the carried denseDf-present columns — every
  // offerable key is honoured by BOTH the aggregate grid AND the mean±SD chart, which
  // resolves the frame asymmetry; results-only virtual keys like STATUS are excluded).
  // Restores the selection by NAME (replay-pending first, else the live selection),
  // dropping any name no longer offered. Setting `.value` is done under a suppress
  // guard so the compute path's own applyGrouping/chart refresh isn't double-fired.
  const updateGroupByColumns = (
    candidates: ReadonlyArray<string>,
    resultsDf: DG.DataFrame,
  ): void => {
    const candidateSet = new Set(candidates);
    // The rebind is INSIDE the suppress guard (BUG-11): `ChangeTable` re-resolves the
    // checked columns against the new frame and fires `onChanged` with a transient value
    // (the old frame's columns merged with the new frame's). Left unguarded, the handler
    // below wrote that transient into `groupBySelectedNames` — and on a later Compute, with
    // `pendingGroupByCols` already consumed, `wanted` read it back as the user's selection,
    // laundering the transient into the persisted recipe.
    suppressGroupByChange = true;
    try {
      ui.input.setColumnsInputTable(groupByInput, resultsDf, (col: DG.Column) => candidateSet.has(col.name));
      // Restore selection by name — replay-pending is one-shot (consumed here so the
      // user's live selection wins on later Computes), else keep the live selection.
      const wanted = pendingGroupByCols ?? groupBySelectedNames;
      if (pendingGroupByCols !== null) pendingGroupByCols = null;
      const restored = dedupeGroupBySelection(
        wanted.filter((nm) => candidateSet.has(nm)), resultsDf);
      groupBySelectedNames = restored.map((c) => c.name);
      groupByInput.value = restored;
    } finally {
      suppressGroupByChange = false;
    }
  };

  // Recompute (or close) the aggregated grid for the current Group-by selection,
  // reading the live results DataFrame from the docked results grid so it always
  // reflects the latest Compute. Empty selection = no grouping (v8: the selection is
  // now a composite tuple of column names, not a single label).
  const applyGrouping = (): void => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const tempBag = view.temp as Record<string, any>;
    // `== null` not `=== undefined`: view.temp returns null for an unset key (BUG-06).
    const grid = tempBag[RESULTS_GRID_TAG] as DG.Grid | null;
    // The By-group mixed-basis/dose-norm notice tracks the LIVE Group-by selection
    // (re-evaluated here on every Group-by change, no recompute), so it renders its
    // own banner via `refreshNotices` rather than the per-Compute collector. The
    // aggregate grid uniquely owns its root, so the banner is captured against that
    // root specifically (never the results fallback) — cleared when grouping is off.
    const priorAggRoot = (tempBag[AGG_GRID_TAG] as DG.Viewer | null)?.root ?? null;
    if (grid == null || groupBySelectedNames.length === 0) {
      refreshNotices(view, [], priorAggRoot);
      closeAggregated(view);
      return;
    }
    // Aggregate the REPORTABLE frame (AD-6 / U1): Exclude λz values do not enter the grouped
    // statistics, so N_<param> of a λz-dependent parameter is the per-parameter n (D11).
    const resultsDf = reportableFrameOf(view, grid.dataFrame);
    const present = groupBySelectedNames.filter((nm) => resultsDf.col(nm) !== null);
    if (present.length === 0) {
      // Every selected column vanished (mapping changed) — clear to no grouping.
      refreshNotices(view, [], priorAggRoot);
      closeAggregated(view);
      return;
    }
    const label = present.map((nm) => groupByFriendly(resultsDf, nm)).join(' · ');
    // group-by-cl-basis-explicit-control F1/F3: the aggregate now honours the user's keys
    // EXACTLY (no silent CL_BASIS auto-append). Where a pooled cell mixes CL and CL/F, the
    // corrupted family cells are WITHHELD (NE) by the blank post-pass below, and the banner
    // discloses it — banner + blank driven by the SAME `mixedBasisKeys` set (never desync).
    // `resolveAggregateGroupKeys` buckets by the full tuple, so a composite CELL decides
    // mixing (not merely one dimension). CL / Vz / Vss are per kg on every row
    // (nca-cl-bw-normalised-units AD-4), so mixed engine dose bases need no withholding.
    const {keys, mixedBasis, mixedBasisKeys} = resolveAggregateGroupKeys(resultsDf, present);
    const byGroupNotices: ComputeNotice[] = [];
    if (mixedBasis) {
      byGroupNotices.push({
        key: 'MIXED_CL_BASIS_AGGREGATE', anchor: 'byGroup', severity: 'warning',
        text:
          `"${label}" pools true CL (IV) and apparent CL/F (EV): ` +
          `CL, Vz, Vss, MRT and AUMC are not comparable across bases and are withheld (NE). ` +
          `Add CL basis or Route to "Group by" to separate them. ` +
          `t½, λz and %extrap are shown but assume an elimination-limited terminal phase — ` +
          `not valid under flip-flop absorption.`,
      });
    }
    // The aggregate grid bold-fonts fragile N cells (N < FRAGILE_N) via an
    // onCellPrepare hook attached at grid creation (F4, swapOrDockAggregated).
    const agg = aggregateNcaResults(resultsDf, keys);
    // F3 Decision A: withhold the basis-corrupted family cells in the pooled mixed rows.
    // The blank is a POST-PASS on the returned frame — it is not an aggregation change —
    // keyed on the SAME `mixedBasisKeys` that fired the banner above. (The parenthetical here
    // used to say `aggregateNcaResults` was "byte-untouched"; that stopped being true when
    // nca-dg-native-missing-values F3 deleted its NaN→FLOAT_NULL shadow clone. The claim
    // this comment is actually making — post-pass, not aggregator edit — still holds, and
    // the PKNCA 60/60 parity pin is unaffected: no aggregated VALUE changed, only which
    // layer writes the null.)
    blankMixedBasisFamily(agg, mixedBasisKeys, keys, IV_EV_INCOMPARABLE_FAMILY);
    // F-S5.2 / D11: how many of each group's profiles an Exclude λz decision withheld — so a
    // λz-dependent MEAN whose EVERY contributor was excluded reads "not reportable", distinct from
    // a blank not-estimated cell (thesis #4). Counted from the RAW frame by the same composite key
    // the aggregate groups on; a hidden helper column read by the grid's single onCellPrepare.
    appendExcludedCounts(view, grid.dataFrame, agg, keys);
    swapOrDockAggregated(view, agg);
    // VR-21 linked selection: selecting group rows in the aggregate grid filters the
    // Mean±SD chart to those groups' curves (analogous to results→profile at group grain).
    // Re-wired here because `agg` is a fresh frame on every Group-by change / Compute.
    wireByGroupSelectionToMeanSdFilter(view, agg, keys);
    // Aggregate grid is now (re)docked — render the banner over its fresh root.
    const aggRoot = (tempBag[AGG_GRID_TAG] as DG.Viewer | null)?.root ?? priorAggRoot;
    refreshNotices(view, byGroupNotices, aggRoot);
  };

  // VR-03 — re-split the docked "Mean ± SD" chart by the current Group-by
  // selection WITHOUT a recompute (the chart reflects the LAST RUN, like the
  // other profile tabs — see MEAN_SD_CTX_TAG). Group-key resolution:
  // Group-by selection (when that column is present in the observations frame)
  // → dose column → single pooled "All" curve. Results-only keys (STATUS, and
  // ROUTE when route is a scalar) are absent from the observations frame, so
  // they fall back to dose / All rather than throwing. Called with
  // `createIfMissing = true` from Compute (initial dock) and `false` from the
  // Group-by change (refresh-only).
  const refreshMeanSdChart = (createIfMissing: boolean): void => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const tempBag = view.temp as Record<string, any>;
    const ctx = tempBag[MEAN_SD_CTX_TAG] as MeanSdCtx | null;
    if (ctx == null) return;
    const existing = tempBag[MEAN_SD_CHART_TAG] as DG.Viewer | null;
    if (existing == null && !createIfMissing) return;

    // Group-key resolution (v8, multi-column-group-by): the selected columns that are
    // present in the observations frame. 0 → single pooled "All" curve; 1 → that
    // column; ≥2 → a synthesized composite key column on a clone (compositeGroupKey).
    // F2: there is NO dose fallback — an empty selection is ONE pooled curve, always;
    // the tool never silently auto-stratifies by dose (the disclosure notice discloses
    // the cross-dose pool instead).
    // Per-row VALUE getter for a selected key on the CHART's frame. For the ROUTE key
    // this maps the observations' per-row route CODE (DERIVED_ROUTE_CODE_COL) to the
    // CANONICAL route label (ROUTE_CODE_TO_STRING) — the SAME label the results ROUTE
    // column (and the aggregate grid) carries — so the chart legend reads `IV-bolus/PO/SC`,
    // NOT the raw source `IV/PO/SC` that `denseDf.col('ROUTE')` would resolve to via
    // case-insensitive lookup (grid↔chart label consistency). Every other key uses its
    // own denseDf column verbatim.
    const routeCodeCol = ctx.denseDf.col(DERIVED_ROUTE_CODE_COL);
    // The selected keys splittable on the CHART frame. A key is chart-splittable when it
    // has its own observations column, OR (CL_BASIS, F4) it can be SYNTHESIZED from the
    // per-row route code — CL_BASIS is a results-side per-row tag with no denseDf column,
    // so without this it would silently drop from the chart while still splitting the grid
    // (the frame asymmetry this topic closes). The offer gate already withholds CL_BASIS
    // unless routeCodeCol is present, so the two frames stay in lock-step.
    const selNames = groupBySelectedNames.filter((nm) =>
      ctx.denseDf.col(nm) !== null || (nm === 'CL_BASIS' && routeCodeCol !== null));
    const chartValueOf = (nm: string): (i: number) => string => {
      if (nm === 'ROUTE' && routeCodeCol !== null) {
        return (i) => routeCodeCol.isNone(i) ? '' :
          (ROUTE_CODE_TO_STRING[routeCodeCol.get(i) as number] ?? '');
      }
      // F4 (group-by-cl-basis-explicit-control): CL basis split on the CHART frame,
      // synthesized from the per-row route code — route 1:1-determines basis. TWO-step:
      // routeCode → canonical route string (ROUTE_CODE_TO_STRING) → isIntravascular
      // predicate → 'CL' (IV) / 'CL/F' (EV). These are the SAME labels the grid CL_BASIS
      // column carries (grid == chart legend). isIntravascular is STRING-typed (no int
      // overload), so the code must be resolved to a route string first.
      if (nm === 'CL_BASIS' && routeCodeCol !== null) {
        return (i) => {
          if (routeCodeCol.isNone(i)) return '';
          const routeStr = ROUTE_CODE_TO_STRING[routeCodeCol.get(i) as number] ?? '';
          return routeStr === '' ? '' : (isIntravascular(routeStr) ? 'CL' : 'CL/F');
        };
      }
      const col = ctx.denseDf.col(nm)!;
      return (i) => col.isNone(i) ? '' : col.getString(i);
    };
    // A single key needs a synthesized column only when its chart values differ from the
    // raw denseDf column: ROUTE (canonical relabel) or CL_BASIS (route-code → basis).
    const needsSynthesis = (nm: string): boolean =>
      (nm === 'ROUTE' || nm === 'CL_BASIS') && routeCodeCol !== null;

    let chartCtx = ctx;
    let groupName: string | null = null;
    if (selNames.length === 1 && !needsSynthesis(selNames[0]))
      groupName = selNames[0];
    else if (selNames.length >= 1) {
      const getters = selNames.map(chartValueOf);
      const dd = ctx.denseDf.clone();
      // Chart group is DISPLAY-ONLY (colour + line split + legend) — it never joins
      // with the drill-in / curve / warning keys, so it uses the human-readable ` · `
      // separator, NOT the byte-identical U+0001 `compositeGroupKey` (whose control
      // char leaked into the legend as a `␁` box). Natural sort still orders numerics
      // within the join (`M · 10` < `M · 100`).
      dd.columns.add(DGNS.Column.fromStrings(MEAN_SD_COMPOSITE_COL,
        Array.from({length: dd.rowCount},
          (_v, i) => getters.map((g) => g(i)).join(GROUP_LABEL_SEP))));
      chartCtx = {...ctx, denseDf: dd};
      groupName = MEAN_SD_COMPOSITE_COL;
    }

    const {groupCount, minCellN, fragileCells} =
      ensureMeanSdChart(view, chartCtx, groupName, useLogY);
    // Notices are collected only on the Compute pass (createIfMissing) — see B4 note:
    // the chart is created once then swapped in place, so gating on the fresh-DOM flag
    // let a banner show once and never re-show; gating on createIfMissing re-shows it
    // each Compute (runs after resetNotices / before flushNotices). The Group-by
    // re-split path (createIfMissing=false) deliberately does NOT collect (no flush /
    // no reset → a notice there would leak stale into the next Compute).
    if (createIfMissing) {
      // F4: PER-CELL fragility (multi-column-group-by) — a heterogeneous composite
      // (most cells N=6, one N=1) keeps the global median above the old `medianN<=1`
      // trigger while a fragile cell is plotted silently (thesis §4). Fire when ANY
      // cell is fragile (minCellN < FRAGILE_N) and NAME the fragile cell(s); collapse
      // to a summary past FRAGILE_LIST_MAX (avoids the gap-18 flood at high composite
      // cardinality). Retains the observed-time serial-drift advice when applicable.
      if (groupCount > 0 && minCellN < FRAGILE_N) {
        collectNotice(view, {
          key: 'MEAN_SD_FRAGMENTED', anchor: 'meanSd', severity: 'info',
          text: fragilityNoticeText(fragileCells, ctx.nominalUsed),
        });
      }
      // F2: cross-dose pooling disclosure (multi-column-group-by BUG-08 fix). When the
      // effective grouping omits Dose on multi-dose data (empty selection, or any
      // selected key that doesn't include the dose column), the pooled curve mixes dose
      // levels — concentration scales with dose, so the pooled mean is misleading. This
      // is off the correctness path (the curve is drawn either way); the notice makes
      // the pool honest and points at the fix. Gated by a log-space dose clusterer so
      // it does NOT fire on a single-nominal-level covariate (theoph: 10 per-kg values,
      // 1 nominal level).
      const doseName = ctx.doseName;
      const doseSelected = doseName !== null && selNames.includes(doseName);
      if (!doseSelected && doseName !== null && ctx.denseDf.col(doseName) !== null) {
        const levels = countDoseClusters(ctx.denseDf.col(doseName)!);
        if (levels >= 2) {
          collectNotice(view, {
            key: 'MEAN_SD_CROSS_DOSE_POOL', anchor: 'meanSd', severity: 'info',
            text:
              `Mean ± SD pools across ${levels} dose levels; concentration scales with ` +
              `dose — select Dose in "Group by" to stratify.`,
          });
        }
      }
    }
  };

  // Detect a single source column mapped to MORE THAN ONE input (e.g. Time and
  // Concentration both pointing at the same column). That configuration computes
  // but then fails / produces garbage downstream, so it must block Compute.
  // Inactive paired columns are `null` (linkSourceVisibility clears them when the
  // source toggle is 'value'), so only ACTIVE column mappings are checked.
  // Returns `[columnName, [input labels…]]` per conflicting column.
  const duplicateColumnMappings = (): Array<[string, string[]]> => {
    const picks: ReadonlyArray<{label: string; col: DG.Column | null}> = [
      {label: 'Subject', col: subjectInput.value},
      {label: 'Time', col: timeInput.value},
      {label: 'Concentration', col: concInput.value},
      {label: 'Dose', col: doseColInput.value},
      {label: 'Route', col: routeColInput.value},
      {label: 'Infusion duration', col: infusionDurationColInput.value},
      {label: 'Body weight', col: bodyWeightColInput.value},
      {label: 'BLQ flag', col: blqFlagColInput.value},
      {label: 'LLOQ', col: lloqColInput.value},
    ];
    const byName = new Map<string, string[]>();
    for (const {label, col} of picks) {
      if (col === null) continue;
      const labels = byName.get(col.name);
      if (labels === undefined) byName.set(col.name, [label]);
      else labels.push(label);
    }
    return [...byName.entries()].filter(([, labels]) => labels.length > 1);
  };

  const refreshComputeEnabled = (): void => {
    const missing: string[] = [];
    if (subjectInput.value === null) missing.push('Subject');
    if (timeInput.value === null) missing.push('Time');
    if (concInput.value === null) missing.push('Concentration');
    const required = missing.length === 0;
    const dupes = duplicateColumnMappings();
    // §3.2 consistency (lambdaz-adjr2-threshold-spec-v2): the acceptance bar R_acc must
    // be ≥ the point-selection gate minR. A softer bar is a DEAD criterion — every
    // reported fit already has adjR² ≥ minR, so adjR² < R_acc can never fire. HARD gate,
    // like the sex constraints below: block Compute, not just tint the field (the inline
    // validator alone did not stop a Compute — user QA 2026-08-25).
    const rAcc = adjRSqAcceptanceInput.value;
    const minR = lambdaMinRSquaredInput.value;
    const rAccBelowGate = rAcc != null && minR != null && rAcc < minR;
    // Subject-sex mapping gates (Slice B): an implausible sex column (too many
    // values) or an unresolved coded value blocks Compute — F% must never run on
    // a stratification the scientist has not fully resolved.
    // C-6: the tooltip LISTS every failing check (the same terms runCompute hard-guards on),
    // not just the first one.
    const problems: string[] = [];
    if (!required) problems.push(`Map these required columns: ${missing.join(', ')}.`);
    if (dupes.length > 0) {
      problems.push(`Each column may map to only one input. Conflict: ` +
        `${dupes.map(([n, l]) => `"${n}" → ${l.join(' & ')}`).join('; ')}.`);
    }
    if (sexTooMany) {
      problems.push(`The Sex column has too many distinct values for a sex axis — pick a ` +
        `column with ≤${SEX_MAX_DISTINCT} categories, or clear it.`);
    }
    if (sexUnassigned) {
      problems.push('Assign a sex (M / F / U) to every value in the Subject Sex section — ' +
        '%F cannot be computed until all values are mapped.');
    }
    if (rAccBelowGate) {
      problems.push(`The adj-R² acceptance threshold (${rAcc}) is below λz min adj-R² (${minR}) — ` +
        'raise it to at least the point-selection gate.');
    }
    // `consumed` (set after a successful compute, cleared on any input change — as the old ▶ icon
    // behaved, user-directed 2026-09-29): a confirmed Setup with nothing changed has nothing to
    // recompute, so Confirm dims until an input changes. It still isn't `disabled` — the tooltip
    // must stay hoverable to explain WHY.
    const enabled = problems.length === 0 && !consumed;
    confirmReady = enabled;
    confirmBtn.style.opacity = enabled ? '' : '0.5';
    confirmBtn.style.cursor = enabled ? 'pointer' : 'not-allowed';
    // Read at hover: after a Setup edit of a confirmed analysis the button is `Apply` (refreshStage).
    ui.tooltip.bind(confirmBtn, () => problems.length > 0 ?
      ui.divV(problems.map((t) => ui.divText(t))) :
      consumed ? 'Already confirmed — change any input to re-run' :
        liveWorkflowState(view).stage.setupChanged ?
          'Apply the changes and recompute — then review the λz fits' :
          'Confirm the setup and compute — then review the λz fits');
    // Via the view hook, not the closure const: this runs during mount, before the strip exists.
    refreshStageUi(view);
  };

  // Reads the live Method-settings inputs and assembles a fresh NcaRules
  // object for each Compute. Called inside the click handler so user edits
  // between Computes flow through. Casts are safe by construction: every
  // input's `items:` array is a const literal of the corresponding sci-comp
  // enum, so the runtime value is always a valid enum member.
  const buildRules = (): NcaRules => ({
    aucMethod: aucMethodInput.value as NcaRules['aucMethod'],
    compensatedSummation: compensatedInput.value ?? false,
    blq: {
      preFirstMeasurable: blqPreInput.value as NcaRules['blq']['preFirstMeasurable'],
      embedded: blqEmbeddedInput.value as NcaRules['blq']['embedded'],
      afterLast: blqAfterInput.value as NcaRules['blq']['afterLast'],
      consecutiveAfterLast: blqConsecutiveInput.value as NcaRules['blq']['consecutiveAfterLast'],
    },
    lambdaZ: {
      mode: lambdaModeInput.value as NcaRules['lambdaZ']['mode'],
      minPoints: lambdaMinPointsInput.value ?? rules.lambdaZ.minPoints,
      minRSquared: lambdaMinRSquaredInput.value ?? rules.lambdaZ.minRSquared,
      excludeCmax: lambdaExcludeCmaxInput.value ?? rules.lambdaZ.excludeCmax,
      // Fall back to the initial `rules.lambdaZ.adjRSquaredFactor` (which
      // comes from DEFAULT_RULES = 1e-4 unless `opts.rules` overrides) so
      // that clearing the field preserves the reference-suite-validated
      // value. Earlier `?? undefined` silently mapped null → undefined →
      // sci-comp default 0, breaking PKNCA parity on Theoph subject 6.
      adjRSquaredFactor: lambdaAdjRSqFactorInput.value ??
        rules.lambdaZ.adjRSquaredFactor,
    },
    extrapWarnPct: extrapWarnPctInput.value ?? rules.extrapWarnPct,
    extrapErrorPct: extrapErrorPctInput.value ?? rules.extrapErrorPct,
    // FR-200 %AUMCextrap threshold. No dedicated input yet (PP-2c — DEFAULT_RULES
    // suffices); project from the base rules so the value threads into the
    // rule-hash and the results-view band without re-reading inputs.
    extrapWarnPctAumc: rules.extrapWarnPctAumc,
  });

  // ── AD-5 recipe seam (lossless round-trip; SCIENCE gate AC-D6) ───────────────
  // `getCurrentRecipe()` serializes EVERY parameter-affecting input — columns
  // (each scalar/column triple), the NONMEM/sparse mode toggles + confirmedRoute,
  // the 4 units, and the full NcaRules (= buildRules() output, the canonical
  // serializer — NOT a parallel one), plus the adj-R² ACCEPTANCE threshold R_acc as
  // its OWN recipe field (it is not part of NcaRules — SAP/presentation, thesis #1 —
  // yet must lock with the analysis, §7). `applyRecipe()` is its inverse, assigning
  // each input back so a Compute reproduces byte-identical parameters AND the same
  // acceptance verdicts. A dropped field would silently change the compute PATH or the
  // verdicts, so the enumeration here MUST stay in lockstep with the
  // schema's AnalysisRecipe.
  const getCurrentRecipe = (): AnalysisRecipe => {
    const sexValueMap: Record<string, string> = {};
    for (const {raw, input} of sexValueInputs) {
      const v = input.value;
      if (v === 'M' || v === 'F' || v === 'U') sexValueMap[raw] = v;
    }
    const columns: RecipeColumns = {
      subject: subjectInput.value?.name ?? null,
      time: timeInput.value?.name ?? null,
      concentration: concInput.value?.name ?? null,
      doseSource: readSourceMode(doseSourceInput),
      doseCol: doseColInput.value?.name ?? null,
      doseScalar: doseScalarInput.value ?? 0,
      routeSource: readSourceMode(routeSourceInput),
      routeScalar: routeScalarInput.value as string,
      routeCol: routeColInput.value?.name ?? null,
      infusionSource: readSourceMode(infusionSourceInput),
      infusionDurationCol: infusionDurationColInput.value?.name ?? null,
      infusionDurationScalar: infusionDurationInput.value ?? 0,
      bodyWeightSource: readSourceMode(bodyWeightSourceInput),
      bodyWeightCol: bodyWeightColInput.value?.name ?? null,
      bodyWeightScalar: bodyWeightScalarInput.value ?? 0,
      blqFlagCol: blqFlagColInput.value?.name ?? null,
      lloqSource: readSourceMode(lloqSourceInput),
      lloqCol: lloqColInput.value?.name ?? null,
      lloqScalar: lloqScalarInput.value ?? 0,
      sexCol: sexColInput.value?.name ?? null,
      sexValueMap,
    };
    const modes: RecipeModes = {
      nonmem: nonmemToggle.value ?? false,
      // R1 SY-F1: the confirmed NONMEM route closure var — round-tripped so replay
      // does not re-trigger the route-confirm dialog (R2-F1 deserialization below).
      confirmedRoute: confirmedRoute,
      sparse: {
        on: sparseToggleInput.value ?? false,
        nominalTimeCol: nominalTimeColInput.value?.name ?? null,
        animalIdCol: animalIdColInput.value?.name ?? null,
        // R1 SY-F3: serialize the '(auto)' sentinel VERBATIM, never null.
        samplingType: samplingTypeInput.value as SamplingType,
      },
    };
    const units: RecipeUnits = {
      time: timeUnitInput.value as string,
      concentration: concUnitInput.value as string,
      dose: doseUnitInput.value as string,
      bodyWeight: bodyWeightUnitInput.value as string,
    };
    // Persist the "Group by" ribbon selection as the underlying RESULTS-COLUMN names
    // (v8: a set, in pick order) so the cache-only study grid can regroup. Empty = no
    // grouping. Sourced from the live selection names (kept in sync by onChanged).
    const groupByCols = [...groupBySelectedNames];
    // v10 (study-grid-passthrough-columns, PP-1): passthrough columns are a DIALOG-time pick
    // (openAddToStudyDialog), NOT a mapping-pane input — so the pane's recipe emits [] and the
    // add-to-study dialog OVERRIDES it. Emitting [] (not undefined) keeps the recipe round-trip
    // deep-compare lossless (recipe-round-trip-tests.ts: no []/undefined drift).
    return {
      columns, modes, ncaRules: buildRules(), units, groupByCols, passthroughCols: [],
      // R_acc rides the recipe as its own field (NOT part of NcaRules — SAP/presentation,
      // thesis #1) so save/load locks the acceptance threshold with the analysis (§7).
      acceptanceRSquared: adjRSqAcceptanceInput.value ?? DEFAULT_ACCEPTANCE_RSQUARED,
    };
  };

  /** Resolve a recipe-stored column name against the live df (a name that no
   *  longer exists drops to null, exactly as a stale autodetect would). */
  const colByName = (name: string | null): DG.Column | null =>
    (name !== null && df.col(name) !== null) ? df.col(name) : null;

  const applyRecipe = (recipe: AnalysisRecipe): void => {
    const c = recipe.columns;
    subjectInput.value = colByName(c.subject);
    timeInput.value = colByName(c.time);
    concInput.value = colByName(c.concentration);
    writeSourceMode(doseSourceInput, c.doseSource);
    doseColInput.value = colByName(c.doseCol);
    doseScalarInput.value = c.doseScalar;
    writeSourceMode(routeSourceInput, c.routeSource);
    routeScalarInput.value = c.routeScalar;
    routeColInput.value = colByName(c.routeCol);
    writeSourceMode(infusionSourceInput, c.infusionSource);
    infusionDurationColInput.value = colByName(c.infusionDurationCol);
    infusionDurationInput.value = c.infusionDurationScalar;
    writeSourceMode(bodyWeightSourceInput, c.bodyWeightSource);
    bodyWeightColInput.value = colByName(c.bodyWeightCol);
    bodyWeightScalarInput.value = c.bodyWeightScalar;
    blqFlagColInput.value = colByName(c.blqFlagCol);
    writeSourceMode(lloqSourceInput, c.lloqSource);
    lloqColInput.value = colByName(c.lloqCol);
    lloqScalarInput.value = c.lloqScalar;
    // Modes (compute-PATH toggles: they change the analytical output).
    nonmemToggle.value = recipe.modes.nonmem;
    // R2-F1: pre-fill the confirmedRoute closure so the route dialog never fires
    // on replay (the recipe already carries the resolved route).
    confirmedRoute = recipe.modes.confirmedRoute;
    sparseToggleInput.value = recipe.modes.sparse.on;
    nominalTimeColInput.value = colByName(recipe.modes.sparse.nominalTimeCol);
    animalIdColInput.value = colByName(recipe.modes.sparse.animalIdCol);
    samplingTypeInput.value = recipe.modes.sparse.samplingType;
    // Units.
    timeUnitInput.value = recipe.units.time;
    concUnitInput.value = recipe.units.concentration;
    doseUnitInput.value = recipe.units.dose;
    bodyWeightUnitInput.value = recipe.units.bodyWeight;
    // Method rules.
    const r = recipe.ncaRules;
    aucMethodInput.value = r.aucMethod;
    compensatedInput.value = r.compensatedSummation;
    blqPreInput.value = r.blq.preFirstMeasurable;
    blqEmbeddedInput.value = r.blq.embedded;
    blqAfterInput.value = r.blq.afterLast;
    blqConsecutiveInput.value = r.blq.consecutiveAfterLast;
    lambdaModeInput.value = r.lambdaZ.mode;
    lambdaMinPointsInput.value = r.lambdaZ.minPoints;
    lambdaMinRSquaredInput.value = r.lambdaZ.minRSquared;
    // R_acc lives on the recipe itself (not in NcaRules). `?? default` keeps a pre-persistence
    // recipe (or a hand-edited file) loading at the nca-studio default rather than blank.
    adjRSqAcceptanceInput.value = recipe.acceptanceRSquared ?? DEFAULT_ACCEPTANCE_RSQUARED;
    lambdaExcludeCmaxInput.value = r.lambdaZ.excludeCmax;
    if (r.lambdaZ.adjRSquaredFactor !== undefined && r.lambdaZ.adjRSquaredFactor !== null)
      lambdaAdjRSqFactorInput.value = r.lambdaZ.adjRSquaredFactor;
    extrapWarnPctInput.value = r.extrapWarnPct;
    extrapErrorPctInput.value = r.extrapErrorPct;
    // Sex column + value map LAST: stage the pending map, then set the column so
    // refreshSexMapping (fired by the change) pre-fills the per-value inputs from
    // the recipe; also call it directly in case the setter's onChanged is deferred.
    pendingSexMap = c.sexValueMap;
    sexColInput.value = colByName(c.sexCol);
    refreshSexMapping();
    // Group-by: the candidates are only known after the first Compute populates
    // them (appendGroupingColumns runs there), so stage the persisted column(s) and
    // let `updateGroupByColumns` restore the selection once they are available.
    pendingGroupByCols = recipe.groupByCols;
  };

  // Treat 0 (or null/undefined) as "unspecified" for the optional scalars
  // that the form exposes — there is no domain-meaningful bodyWeight or
  // infusionDuration of 0, so we collapse it to null which sci-comp
  // interprets as "skip the weight-normalised parameters / no infusion".
  const scalarOrNull = (v: number | null): number | null =>
    (v != null && v > 0) ? v : null;

  const runCompute = (runOpts?: {
    silent?: boolean;
    /** Confirm pre-warn (F-S3.5 / AD-7): run the pipeline on the edited Setup with a THROW-AWAY
     *  controller + AuditLog and report how many review decisions a real Confirm would reset.
     *  Nothing on screen, in the controller, or in the session audit log changes. */
    dryRun?: (resets: number) => void;
  }): void => {
    // `silent` is set ONLY by the replay-on-open hydration path (AC-D1): suppress
    // the interactive confirm dialogs + transient success toasts the user already
    // saw when first computing this analysis. A user click (play icon / manual
    // re-run after open) passes no opts → silent=false → everything shows again.
    // (Named `runOpts` to avoid shadowing the outer `opts` = MountMappingPaneOptions
    // that runPipeline reads for `onComputed`.)
    const silent = runOpts?.silent ?? false;
    if (consumed)
      return;
    if (subjectInput.value === null || timeInput.value === null || concInput.value === null)
      return;
    // The greyed icon does not block clicks — hard-guard against a column mapped
    // to more than one input (would compute then fail / mislabel downstream).
    if (duplicateColumnMappings().length > 0)
      return;
    // Hard-guard the sex-mapping gates too (the greyed icon is only a visual): an
    // implausible sex column (> 8 distinct) or an unresolved value must never
    // compute — F% would stratify on a mis-picked / half-mapped axis.
    if (sexTooMany || sexUnassigned)
      return;
    // Hard-guard §3.2: the adj-R² acceptance bar must be ≥ the point-selection gate.
    // The inline field validator alone did not stop a Compute (user QA 2026-08-25); the
    // greyed icon is only a visual, so this is the real block — mirrors the sex guards.
    {
      const rAcc = adjRSqAcceptanceInput.value;
      const minR = lambdaMinRSquaredInput.value;
      if (rAcc != null && minR != null && rAcc < minR)
        return;
    }

    // ── Event-driven (NONMEM) normalization (AD-1) ───────────────────────────
    // When the mode is active, transform the uploaded df into the canonical
    // observations frame; EVERY consumer below runs on `srcDf` with the
    // normalized column names (AC-B10 four-consumer threading). Standard mode
    // leaves srcDf === df and the names equal to the picked inputs (unchanged).
    const nmActive = (nonmemToggle.value ?? false) && eventDetection.isEventDriven;
    let srcDf = df;
    let nmSubjectName = subjectInput.value.name;
    let nmTimeName = timeInput.value.name;
    let nmConcName = concInput.value.name;
    let nmContract: ColumnContract | null = null;
    let nmRouteWarn: ReadonlyMap<string, ReadonlyArray<ParameterWarning>> | undefined;
    if (nmActive) {
      const dayName = findColByPattern(df, /^day$|study.?day/i);
      const ec: EventContract = {
        id: subjectInput.value.name, time: timeInput.value.name, dv: concInput.value.name,
        amt: eventDetection.amt!.name, evid: eventDetection.evid!.name,
        cmt: eventDetection.cmt?.name, rate: eventDetection.rate?.name, mdv: eventDetection.mdv?.name,
        addl: eventDetection.addl?.name, ii: eventDetection.ii?.name, ss: eventDetection.ss?.name,
        route: routeColInput.value?.name, bodyWeight: bodyWeightColInput.value?.name,
        day: dayName ?? undefined,
        // Pass every column as a covariate candidate; `collectCovariateNames`
        // (event-driven.ts) strips the event-mechanics columns (EVID/AMT/TIME/DV/
        // CMT/RATE/MDV/ADDL/II/SS + id/time/dv) before hoisting, so only study
        // metadata (Sex / DOSE_MGKG / NOAEL / …) survives onto the obs rows.
        covariates: df.columns.names(),
      };
      const norm = normalizeEventDrivenObservations(df, ec, {confirmedRoute});
      if (norm.errors.length > 0) {
        renderValidationErrors(errorPanel, norm.errors as ValidationError[]);
        return;
      }
      // AD-5 mandatory route-confirm gate (AC-C2): inferred route → block Compute.
      if (norm.routeConfirmRequired) {
        if (!silent) {
          showRouteConfirmDialog((code) => {confirmedRoute = code; onAnyChange();});
          return;
        }
        // a1 (silent replay): a saved analysis reached the route gate with no
        // confirmedRoute. INVARIANT (verified): this is unreachable from a normal
        // save — the Add-to-study action no-ops until `lastPerProfile !== null`
        // (saveBtn handler), which is set only by a SUCCESSFUL compute, and a
        // successful compute on a NONMEM-inferred dataset requires passing this very
        // gate (which sets confirmedRoute). So a saved NONMEM-inferred recipe always
        // carries confirmedRoute ≠ null; null here means an out-of-band recipe
        // (hand-edited JSON, or a pre-`confirmedRoute` schema migration). Don't block
        // a passive open with an undismissable dialog: proceed with the auto-inferred
        // route already encoded in `norm.observations` (the route the original
        // compute used; the inference is deterministic on the same data). The
        // per-subject inferred-route warnings in WARNINGS_JSON are the audit record.
        // Surface the anomaly to the dev console (thesis #4) without a UI element.
        // Falls through — no early return.
        console.warn('NCA: silent replay reached the NONMEM route-confirm gate with ' +
          'confirmedRoute=null (out-of-band recipe) — proceeding with the auto-inferred ' +
          'route; verify the add-to-study save gate. See WARNINGS_JSON for per-subject routes.');
      }
      clearValidationErrors(errorPanel);
      srcDf = norm.observations;
      nmSubjectName = NM_SUBJECT_COL; nmTimeName = NM_TAD_COL; nmConcName = NM_CONC_COL;
      nmRouteWarn = norm.warningsBySubject;
      nmContract = buildNonmemContract(srcDf, concUnitInput.value as ConcentrationUnit,
        bodyWeightColInput.value?.name ?? null);
    }

    // ── Subject-sex value→M/F/U map ──────────────────────────────────────────
    // Apply the designated value→sex map to the SOURCE frame here — BEFORE the
    // composite pooling AND the carriage below — so BOTH the NONMEM composite
    // cell key (`poolCompositeByCell`) and the dense carried stratum read the
    // resolved M/F/U value, never a raw `1`/`2` that would silently merge/split
    // strata downstream; keep this the single application point. Standard mode aliases `srcDf === df`, so clone
    // first to keep the user's displayed grid on its raw tokens; the NONMEM frame
    // is already a private working copy. An all-identity map (clean M/F data) is
    // skipped — the Slice-A fast path, no clone, decorators read the column as-is.
    const sexMapName = getSexColumnName();
    const sexMap = buildSexMap();
    if (sexMapName !== null && srcDf.col(sexMapName) !== null && sexMapChangesValues(sexMap)) {
      if (srcDf === df) srcDf = df.clone();
      applySexMap(srcDf, sexMapName, sexMap);
    }

    // ── Destructive-composite pooling (IMPL-GAP-2, auto-detect) ──────────────
    // For a destructive-composite TK design (e.g. 06_rat_tk28d_nonmem — animals
    // sampled once each), the per-pseudo-subject profiles are degenerate (≈1 obs
    // each), so the dense per-subject pipeline produces hundreds of single-point
    // profiles. The meaningful analysis is the composite-mean per
    // dose×sex×occasion CELL (the PKNCA per-cell card). Auto-detected via
    // per-pseudo-subject non-viability: the dense 03/serial NONMEM path stays
    // byte-unchanged (it IS dense-viable → never pools). Pooling feeds ONLY the
    // dense path (pipeline + charts + grouping + TK); the per-animal Sparse panel
    // keeps the un-pooled normalized frame. Surfaced via a toast + self-evident
    // cell-keyed subject names (never silent).
    let denseDf = srcDf;
    let nmComposite = false;
    let compositeDoseGroupName: string | null = null;
    // NONMEM composite cells: the member-observation hash per cell (AC-S1.4.3) — a review
    // decision on a cell resets when an animal moves between cells, even if the cell id does not.
    let compositeMemberHash: ReadonlyMap<string, string> | undefined;
    if (nmActive && !denseProfilesViable(srcDf, NM_SUBJECT_COL)) {
      compositeDoseGroupName = findColByPattern(srcDf, /dose.?mg.?kg|dose.?level/i);
      // Designated sex column (A3) — replaces the old `/^sex$/i` heuristic. The
      // value→M/F/U map (Slice B) is already applied to `srcDf` BEFORE this point,
      // so the composite cell key (`D{dose}_{sex}_Day{occ}`) keys on the resolved
      // M/F/U value, not a raw `1`/`2`.
      // poolCompositeByCell tolerates a name absent from srcDf (key omits sex).
      const sexName = getSexColumnName();
      const bwName = bodyWeightColInput.value?.name ?? null;
      const {pooled, cellCount, memberHashByCell, mixedCompound, warningsByCell} = poolCompositeByCell(srcDf, {
        doseGroupCol: compositeDoseGroupName, sexCol: sexName, bodyWeightCol: bwName,
      });
      // BUG-07 pattern: a cell holding animals of ≥ 2 compounds would average them into one
      // composite — refuse the compute instead (same code + wording as the cassette block).
      if (mixedCompound.length > 0) {
        renderValidationErrors(errorPanel, mixedCompound.map((m) => ({
          code: 'MULTI_ANALYTE_NOT_SUPPORTED' as const,
          message: mixedCompoundMessage(m, `Composite cell '${m.group}'`),
          column: m.column,
          details: {kind: 'cassette-mixed-group', column: m.column, group: m.group, values: [...m.values]},
        })));
        return;
      }
      compositeMemberHash = memberHashByCell;
      denseDf = pooled;
      nmComposite = true;
      // Rebuild the contract over the pooled frame. The cell's body weight is the pooler's cell
      // mean over the members with a dose and a weight (nca-cl-bw-normalised-units AD-5b), never
      // one anchor animal's: per-kg CL = mean dose / mean BW / AUC, exact for weight-based dosing.
      nmContract = buildNonmemContract(denseDf, concUnitInput.value as ConcentrationUnit,
        denseDf.col(CELL_BW_COL) !== null ? CELL_BW_COL : null);
      // The cell's per-kg notes (partial BW, not weight-based) join that cell's profile warnings.
      if (warningsByCell.size > 0) {
        const merged = new Map<string, ParameterWarning[]>();
        for (const [k, v] of nmRouteWarn ?? new Map<string, ReadonlyArray<ParameterWarning>>()) merged.set(k, [...v]);
        for (const [cell, w] of warningsByCell) merged.set(cell, [...(merged.get(cell) ?? []), ...w]);
        nmRouteWarn = merged;
      }
      if (!silent) {
        grok.shell.info(
          `Destructive-composite NONMEM detected — pooled ${srcDf.rowCount} observations into ` +
          `${cellCount} dose×sex×day composite-mean profiles.`);
      }
    }

    // Per-row-route Option A: when a Route column is mapped, tolerantly encode
    // its strings to a transient int32 RouteCode column (DERIVED_ROUTE_CODE_COL)
    // BEFORE buildContract. Blocking errors (unmappable column / IV-adjacent
    // unmappable rows) render in the error panel and abort; per-subject warnings
    // (bare-IV default, within-subject inconsistency, …) flow into the per-
    // profile WARNINGS_JSON via runNcaPipeline options. Skipped in NONMEM mode —
    // the transform already encoded the derived ~RouteCode column.
    let routeCol: DG.Column | null = nmActive ? denseDf.col(DERIVED_ROUTE_CODE_COL) : null;
    let routeWarn: RouteColumnResult['warningsBySubject'] | undefined;
    if (!nmActive && routeColInput.value !== null) {
      const norm = normalizeRouteColumn(routeColInput.value, subjectInput.value, {
        infusionDurationCol: infusionDurationColInput.value,
        infusionScalarPositive: (infusionDurationInput.value ?? 0) > 0,
      });
      if (norm.errors.length > 0) {
        renderValidationErrors(errorPanel, norm.errors);
        return;
      }
      // Add the transient RouteCode column to `srcDf` — the frame the pipeline +
      // carriage actually read. When a sex map cloned `srcDf` (standard mode),
      // this lands on the clone, not the user's grid `df`; absent a clone
      // `srcDf === df`, so behaviour is unchanged (the column is added to `df` and
      // hidden from the grid below).
      const existing = srcDf.col(DERIVED_ROUTE_CODE_COL);
      if (existing !== null) srcDf.columns.remove(DERIVED_ROUTE_CODE_COL);
      srcDf.columns.add(DGNS.Column.fromInt32Array(DERIVED_ROUTE_CODE_COL, norm.codes));
      routeCol = srcDf.col(DERIVED_ROUTE_CODE_COL);
      routeWarn = norm.warningsBySubject;
      // Keep the internal encoded column out of the user's main grid (no-op when
      // it was added to a sex-map clone rather than the grid frame).
      const gc = view.grid.col(DERIVED_ROUTE_CODE_COL);
      if (gc !== null) gc.visible = false;
    }
    if (nmActive) routeWarn = nmRouteWarn as RouteColumnResult['warningsBySubject'] | undefined;

    // Resolve the sparse-design role columns (UC-04) so they enter the contract
    // when the Sparse-design toggle is on. Mirrors the heuristic fallback in
    // `applySparseNca`: explicit picker > Subject (animal) / name heuristic. The
    // SamplingType role is the auto-detected column (the choice declares a scalar
    // topology, routed separately via `declaredTopology`); when a SamplingType
    // column exists, populating it here mounts `validateSparseDesign` below.
    const sparseOn = sparseToggleInput.value ?? false;
    let sparseNominalCol: DG.Column | null = null;
    let sparseAnimalCol: DG.Column | null = null;
    let sparseTypeCol: DG.Column | null = null;
    if (sparseOn) {
      // In NONMEM mode the pooling key is TAD (the normalized time) and the
      // animal is the pseudo-subject — both on `srcDf`.
      const nomName = nmActive ? NM_TAD_COL :
        (nominalTimeColInput.value?.name ?? findColByPattern(df, /nominal.*time/i));
      sparseNominalCol = nomName !== null ? srcDf.col(nomName) : null;
      sparseAnimalCol = nmActive ? srcDf.col(NM_SUBJECT_COL) : (animalIdColInput.value ?? subjectInput.value);
      const stName = findColByPattern(srcDf, /sampling.?type|^stype$/i);
      sparseTypeCol = stName !== null ? srcDf.col(stName) : null;
    }

    const contract = nmContract ?? buildContract({
      subject: subjectInput.value, time: timeInput.value, conc: concInput.value,
      doseCol: doseColInput.value, bodyWeightCol: bodyWeightColInput.value,
      blqFlagCol: blqFlagColInput.value, lloqCol: lloqColInput.value,
      routeCol,
      nominalTimeCol: sparseNominalCol, animalIdCol: sparseAnimalCol,
      samplingTypeCol: sparseTypeCol,
    });
    // Mount the AD-2 pre-compute topology cross-check (UC-04). `validateSparseDesign`
    // no-ops unless all three sparse roles are on the contract; when a SamplingType
    // column declares 'destructive' but animals span ≥ 2 nominal times, it surfaces
    // SPARSE_TOPOLOGY_MISMATCH as early mapping-time feedback (before the bootstrap
    // runs) — distinct from sparseAuc's post-compute note in the Sparse NOTES column.
    // Transient design-warning toasts — suppress on a silent replay-open (the user
    // saw them on the first compute); the SPARSE NOTES column + WARNINGS_JSON keep
    // the persistent record (thesis #4).
    if (sparseOn && !nmActive && !silent) {
      for (const w of validateSparseDesign(srcDf, contract))
        grok.shell.warning(`${w.code}: ${w.message}`);
    }
    // Scalar values are read UNCONDITIONALLY here, regardless of each input's
    // Pattern-A source toggle. This is intentional: when a source is 'column'
    // the column path takes over downstream (buildContract adds the column to the
    // contract; the extractor reads it per row / per subject), so the (hidden)
    // scalar is simply ignored — it never needs to be zeroed.
    const scalars: ProfileScalars = {
      dose: doseScalarInput.value ?? 0,
      route: ROUTE_STRING_TO_CODE[routeScalarInput.value as string],
      doseUnits: doseUnitInput.value as DoseUnit,
      concentrationUnits: concUnitInput.value as ConcentrationUnit,
      timeUnits: timeUnitInput.value as TimeUnit,
      bodyWeight: scalarOrNull(bodyWeightScalarInput.value),
      bodyWeightUnit: bodyWeightUnitInput.value as BodyWeightUnit,
      infusionDuration: scalarOrNull(infusionDurationInput.value),
      lloq: lloqScalarInput.value ?? 0,
    };

    // Pipeline body as a closure so the AD-5 confirmation dialog (below) can
    // invoke it on OK without duplicating ~100 lines. Captures `df`, `contract`,
    // `scalars`, and the view/input closures.
    // The sparse panel's mapping for THIS compute — shared by the real sparse panel and the
    // Confirm dry run so both compute the SAME composites.
    const sparseMappingNow = (acceptance: AcceptanceThresholds, bootstrapIterations?: number): SparseComputeMapping => {
      const sparseEnabled = sparseToggleInput.value ?? false;
      const stChoice = samplingTypeInput.value;
      const declaredTopology: 'destructive' | 'batch' | 'serial' | null =
        (sparseEnabled && (stChoice === 'destructive' || stChoice === 'batch' || stChoice === 'serial')) ?
          stChoice : null;
      // BUG-19 AD-6: the NONMEM dose basis (absolute `Dose_mg`, never AMT × BW) — one pure helper.
      const dose = sparseDoseSource({
        mode: nmComposite ? 'nonmem-composite' : (nmActive ? 'nonmem' : 'standard'),
        doseColName: doseColInput.value?.name ?? null, doseUnit: doseUnitInput.value as DoseUnit,
        compositeDoseGroupCol: compositeDoseGroupName,
      });
      return {
        subjectName: nmSubjectName, timeName: nmTimeName, concName: nmConcName, doseColName: dose.doseColName,
        bwColName: bodyWeightColInput.value?.name ?? null,
        units: {
          time: timeUnitInput.value as TimeUnit, concentration: concUnitInput.value as ConcentrationUnit,
          dose: dose.doseUnit, bodyWeight: bodyWeightUnitInput.value as BodyWeightUnit,
        },
        // The route codes the dense run reads (AD-5 / S-1): standard mode's derived column on `srcDf`
        // — only when THIS compute wrote it (`routeCol`; a stale one from an unmapped Route column
        // stays on the frame) — NONMEM's normalized `~RouteCode`; else the route scalar.
        route: {
          codesCol: (nmActive ? srcDf.col(DERIVED_ROUTE_CODE_COL) !== null : routeCol !== null) ?
            DERIVED_ROUTE_CODE_COL : null,
          scalar: routeScalarInput.value as string,
        },
        lloq: {colName: lloqColInput.value?.name ?? null, scalar: lloqScalarInput.value ?? 0},
        blqColName: nmActive ? null : (blqFlagColInput.value?.name ?? null),
        nominalTimeName: nmActive ? (sparseEnabled ? NM_TAD_COL : null) :
          (sparseEnabled ? (nominalTimeColInput.value?.name ?? null) : null),
        animalIdName: nmActive ? (sparseEnabled ? NM_SUBJECT_COL : null) :
          (sparseEnabled ? (animalIdColInput.value?.name ?? null) : null),
        declaredTopology, acceptance, bootstrapIterations,
      };
    };
    // The Confirm pre-warn dry run (F-S3.5 / AD-7) runs the SAME pipeline on the edited Setup, on a
    // throw-away controller seeded with the live overrides + decisions and a THROW-AWAY AuditLog
    // (runNcaPipeline otherwise appends to the session singleton — phantom compute entries on a
    // Cancel; architect finding 2). Counts the decisions the real Confirm would reset.
    const runDryRun = (ncaRules: NcaRules, report: (resets: number) => void): void => {
      const live = (view.temp as Record<string, unknown>)[NCA_CONTROLLER_TAG] as NcaController | null;
      const tmp = new NcaController(denseDf, contract, scalars, ncaRules,
        {log: new AuditLog(), ...(routeWarn !== undefined ? {routeWarnings: routeWarn} : {})});
      tmp.setLzOverrides(live?.getLzOverrides() ?? []);
      tmp.setReviewLedgers(live?.getReviewLedgers() ?? []);
      if (tmp.computeAll().status !== 'completed') {report(0); return;}
      const t: ResultsColorThresholds = {
        extrapWarnPct: ncaRules.extrapWarnPct, extrapErrorPct: ncaRules.extrapErrorPct,
        lambdaMinRSquared: ncaRules.lambdaZ.minRSquared,
        lambdaAcceptanceRSquared: adjRSqAcceptanceInput.value ?? DEFAULT_ACCEPTANCE_RSQUARED,
        lambdaMinPoints: ncaRules.lambdaZ.minPoints, extrapWarnPctAumc: ncaRules.extrapWarnPctAumc,
      };
      const viable = nmComposite ? denseProfilesViable(denseDf, NM_SUBJECT_COL) :
        denseProfilesViable(df, subjectInput.value?.name ?? null);
      // Composites without docking and with a short bootstrap — no bootstrap-derived field is in
      // any decision fingerprint (AD-7 determinism premise), so the count is exact.
      const sparse = computeSparseForRecipe(srcDf, sparseMappingNow(t, 20)).result;
      const composites = sparseCompositeUnits(sparse);
      const snap = tmp.reviewSnapshot(t, {viable, composites, memberHashByProfile: compositeMemberHash});
      report(tmp.reconcileDecisions(snap.fingerprints));
    };
    const runPipeline = (): void => {
    // Build the rules ONCE — both `runNcaPipeline` and the grid's
    // color-coding thresholds project from this same object. Deriving
    // thresholds from the assembled NcaRules (rather than re-reading
    // inputs) guarantees that the color legend reflects exactly what
    // sci-comp computed with, even if `buildRules()` ever grows
    // cross-field validation / clamping / coercion.
      const ncaRules = buildRules();
      // Route the whole-cohort Compute through a persistent NcaController (one per
      // view, via view.temp) so the interactive λz editor can later recompute a
      // SINGLE subject in place (VR-11). computeAll() wraps runNcaPipeline and is
      // byte-identical to it on an override-free frame (AC-3.2); inputs are
      // refreshed from live UI each Compute (AC-1.5). `== null` catches the DG-bag
      // never-set null (guardrails).
      if (runOpts?.dryRun !== undefined) {
        runDryRun(ncaRules, runOpts.dryRun);
        return;
      }
      const ctrlBag = view.temp as Record<string, unknown>;
      let controller = ctrlBag[NCA_CONTROLLER_TAG] as NcaController | null;
      if (controller == null) {
        controller = new NcaController(denseDf, contract, scalars, ncaRules,
          routeWarn !== undefined ? {routeWarnings: routeWarn} : {});
        ctrlBag[NCA_CONTROLLER_TAG] = controller;
      } else
        controller.setInputs(denseDf, contract, scalars, ncaRules, routeWarn);
      // A reopened analysis hands its persisted λz overrides + review decisions to the pane
      // BEFORE the first compute (nca-workflow-stages F-S2.1/F-S2.2): they are seeded into the
      // controller once, so the very first compute already applies the analyst's overrides.
      applyPendingWorkflowState(view, controller);

      const outcome = controller.computeAll();
      if (outcome.status === 'validation-failed') {
        renderValidationErrors(errorPanel, outcome.errors);
        // Don't mark consumed — user should be able to click again after fixing
        // the error (the click on its own without input changes is a re-attempt).
        return;
      }
      clearValidationErrors(errorPanel);
      // AC-S2.1.3: overrides whose point vanished or whose concentration changed were ignored by
      // this compute and dropped from the store — never silently re-applied to a different value.
      {
        const dropped = controller.consumeReconciled();
        if (dropped > 0) {
          grok.shell.warning(
            `${dropped} manual λz override(s) reverted to auto — the source ` +
            `(time, conc) changed since they were recorded.`);
        }
      }
      // A successful compute of the live Setup IS the Confirm (D17): record the confirmed recipe
      // snapshot + fingerprint (AD-5), unless this exact recipe is already the confirmed one (a
      // silent reopen replay, or a recompute of an unchanged Setup — the timestamp then keeps
      // the moment the analyst actually confirmed). A migrated entry carries an empty
      // fingerprint and is re-stamped here.
      {
        const recipeNow = getCurrentRecipe();
        const fpNow = computeRecipeFingerprint(recipeNow);
        if (getConfirmed(view)?.fingerprint !== fpNow)
          setConfirmed(view, {recipe: recipeNow, fingerprint: fpNow, at: new Date().toISOString()});
      }
      // Everything below docks / refreshes the result panels from `outcome` (this compute's
      // held results). Kept as ONE closure so a no-change return to Review / Report re-docks it
      // without recomputing (P6 / AC-S3.4.2).
      const present = (): void => {
      // B4: start a fresh per-Compute notice collection and tear down the prior
      // run's banners. A notice the user dismissed last run reappears now if it
      // still applies; one that no longer applies disappears. Decorators below
      // contribute via `collectNotice`; `flushNotices` renders them at the end.
      // DELIBERATELY placed AFTER the validation-failed early return (1278): a
      // failed Compute is a no-op re-attempt (see 1280) that does NOT re-dock the
      // results grid, so the prior run's results AND their banners both stay —
      // resetting here would strip the banners while leaving the stale results
      // displayed (e.g. a pooled MEAN_CL that should be withheld, with no explaining banner), an
      // incoherent, less-honest state. Reset tracks the results
      // it describes: cleared+rebuilt exactly when new results replace the old.
        resetNotices(view);
        // DOCK ORDER IS THE LAYOUT, not just crash-avoidance. Results FIRST: it
        // FILLs onto `view.grid`, forming the tables group that the charts row
        // then splits off. Doing it the other way round would leave
        // `findNode(view.grid.root)` pointing into a re-parented tree — the
        // documented dock-manager footgun. Keep this order.
        const thresholds: ResultsColorThresholds = {
          extrapWarnPct: ncaRules.extrapWarnPct,
          extrapErrorPct: ncaRules.extrapErrorPct,
          lambdaMinRSquared: ncaRules.lambdaZ.minRSquared,
          // R_acc is the ONE threshold NOT projected from NcaRules (it has no sci-comp
          // counterpart — SAP/presentation, thesis #1). Read straight from its own input;
          // fall back to the default if the field was cleared. See the field's contract in
          // ResultsColorThresholds and the AD-4 note.
          lambdaAcceptanceRSquared: adjRSqAcceptanceInput.value ?? DEFAULT_ACCEPTANCE_RSQUARED,
          lambdaMinPoints: ncaRules.lambdaZ.minPoints,
          // Project from the just-built NcaRules (not a re-read input) per the
          // Rules→presentation projection invariant (FR-200 PP-2b).
          extrapWarnPctAumc: ncaRules.extrapWarnPctAumc,
        };
        const resultsNode = swapOrDockResults(view, outcome.results, thresholds);
        // #11: propagate constant-per-subject metadata (compound / dose / sex / ...)
        // into the results frame so the user can group by it, then refresh the
        // Group-by control's options. `outcome.perProfile` is in results-row order,
        // so `appendGroupingColumns` aligns by index. Re-running on every Compute
        // keeps the candidate list and any live aggregation in sync with a changed
        // mapping. See adapter/grouping-columns.ts + adapter/results-aggregator.ts.
        const profiles = outcome.perProfile.map((p) => p.index);
        const carried = appendGroupingColumns(outcome.results, denseDf, profiles, contract);
        // F7 (VR-11): project the controller's per-subject λz-override state onto the
        // results frame as the hidden `_lz_overridden` flag, so the ACCEPTANCE cell shows
        // the ✏️ indicator + names the analyst. Written BEFORE the profile re-apply below
        // so the service column is hidden in the same pass (an `_lz*` column). On a fresh
        // study this is empty; on a reopened study with sticky overrides it repopulates
        // the ✏️ immediately after Compute.
        writeSubjectOverrideFlags(
          outcome.results, controller.overriddenSubjects(), controller.overriddenSubjectTimes());
        // Re-apply the column profile now that the grouping columns EXIST. Ordering
        // is load-bearing: `swapOrDockResults` above configured the grid while the
        // frame still held only the fixed schema, so without this the carried
        // compound / dose / sex columns land wherever DG appends them (the far
        // right) instead of beside the subject — where a Discovery lead scanning a
        // cassette needs them. Cheap and idempotent.
        {
        // eslint-disable-next-line @typescript-eslint/no-explicit-any
          const bag = view.temp as Record<string, any>;
          const grid = bag[RESULTS_GRID_TAG] as DG.Grid | null;
          if (grid != null) {
            applyResultsColumnProfile(grid, resultsExtended(view));
            // Re-assert the problem-first sort AFTER the grouping-column mutation +
            // profile re-apply above, so it is the last word on the row order on every
            // Compute — robust to the frame mutation clearing a prior sort (Q6).
            applyResultsDefaultSort(grid);
          }
        }
        // F1 (multi-column-group-by): the Group-by CANDIDATES are the carried metadata
        // columns PRESENT IN THE OBSERVATIONS FRAME (`denseDf`) — so every offerable key
        // is honoured by BOTH consumers (aggregate grid reads results; mean±SD chart
        // reads denseDf), resolving the frame asymmetry.
        const groupCandidates = groupByCandidates(carried, outcome.results, denseDf, DERIVED_ROUTE_CODE_COL);
        // ROUTE is re-offered as an EXPLICIT key (user decision 2026-07-09, "Variant A"):
        // a DMPK chemist compares IV vs PO, and grouping by Route also separates the
        // CL/CL·F basis (route 1:1-determines it), so a mixed-basis group's withheld family
        // becomes populated. It is offered ONLY when BOTH frames can split by it (the rule lives
        // in `groupByCandidates`, shared with the headless demo seed). Carry the user's OWN route
        // column name onto the internal ALL-CAPS `ROUTE` key via friendlyName, so grid headers /
        // aggregate labels read "Route" (their column) instead of the internal "ROUTE". (The
        // stored key + the results column stay canonical `ROUTE` — only the display label changes.)
        if (readSourceMode(routeSourceInput) === 'column' && routeColInput.value !== null)
          outcome.results.col('ROUTE')?.setTag('friendlyName', routeColInput.value.name);
        // group-by-cl-basis-explicit-control F4: CL basis re-offered as an EXPLICIT key,
        // symmetric to ROUTE. CL basis is COARSER than Route (CL/F pools PO+SC), so it is a
        // distinct analysis axis, not a redundant one. Offered ONLY when the results frame
        // actually spans >1 basis AND the chart can synthesize the split from the per-row
        // route code (route 1:1-determines basis) — the same lock-step that keeps grid==chart
        // (see chartValueOf CL_BASIS branch). Its friendlyName "CL basis" heads the grid group
        // column; the internal key stays canonical `CL_BASIS`. DOSE_NORM stays IMPLICIT (a
        // data artifact, not an axis a chemist picks — user decision).
        if (groupCandidates.includes('CL_BASIS')) outcome.results.col('CL_BASIS')!.setTag('friendlyName', 'CL basis');
        updateGroupByColumns(groupCandidates, outcome.results);

        // Per-subject profiles are dense enough for the dense decorators AND for per-profile
        // review units (a destructive / batch design is not: its per-animal rows are Not
        // applicable and the group composite is reviewed — SF-3). Composite mode judges the POOLED
        // frame. (Hoisted from the decorator block below: the review reconcile needs it first.)
        const denseViable = nmComposite ?
          denseProfilesViable(denseDf, NM_SUBJECT_COL) :
          denseProfilesViable(df, subjectInput.value?.name ?? null);
        // D13 auto-reset for the per-profile units, BEFORE any aggregate reads the exclusion set:
        // a decision the fit no longer supports moves to history and must not exclude anything.
        // Sparse-composite ledgers are reconciled after the sparse panel computes (below).
        const denseReview = controller.reviewSnapshot(thresholds,
          {viable: denseViable, memberHashByProfile: compositeMemberHash});
        const denseResets = controller.reconcileDecisions(denseReview.fingerprints,
          (k: string) => !k.startsWith(SPARSE_KEY_PREFIX));
        // The ONE exclusion feed (AD-6 / U1): every aggregate consumer below reads this blanked
        // clone, never the raw results frame — an Exclude λz profile is treated exactly like a
        // not-estimated one. The per-subject grid keeps the raw frame (value visible + marked).
        const reportable = (): DG.DataFrame => reportableFrameOf(view, outcome.results);

        // ── Charts row ────────────────────────────────────────────────────────
        // TAB ORDER IS DOCK ORDER, and the user-directed order (2026-08-13) is
        // Mean ± SD → Combined → Facet. So the Mean±SD block runs HERE, ahead of
        // the other two, instead of below with the table decorators its
        // `nmDoseName` sibling serves. This block also CREATES the charts group
        // (see {@link CHARTS_ZONE_TAG}), so it must stay above `ensureProfileChart`
        // for a second, independent reason. `Facet` needs no ordering care — it
        // docks from a 100 ms timeout and therefore always lands last.
        //
        // Composite mode groups the decorators (TK accumulation, dose-prop, F%) on
        // the NOMINAL per-kg dose (e.g. DOSE_MGKG) — stable across days within a
        // dose level — not the per-cell mean absolute Dose_mg (which drifts with
        // each day's animal body weights and would split Day-1 from Day-N). In
        // NONMEM mode the dose lives in the normalized `Dose_mg` column. Declared
        // here rather than beside its F% / dose-prop / TK consumers because the
        // Mean±SD render context needs it too — ONE declaration, read by both.
        // The SAME dose column the sparse compute reads (AD-6: one shared expression).
        const nmDoseName = sparseDoseSource({
          mode: nmComposite ? 'nonmem-composite' : (nmActive ? 'nonmem' : 'standard'),
          doseColName: doseColInput.value?.name ?? null, doseUnit: doseUnitInput.value as DoseUnit,
          compositeDoseGroupCol: compositeDoseGroupName,
        }).doseColName;
        // VR-03 group mean ± SD profile chart. Time key = a mapped nominal-time
        // column when present (or TAD in NONMEM mode — already nominal), else the
        // rendered observed time; a serial-drift toast fires from refreshMeanSdChart
        // if observed-time pooling fragments to ~1 obs/cell. Stored in temp so the
        // Group-by control re-splits the chart without a recompute (chart reflects
        // the LAST RUN, like Combined/Facet). Consumes the SAME normalized
        // observations frame `denseDf` the dense profile tabs render from.
        const meanSdNominalName = nmActive ? NM_TAD_COL :
          ((nominalTimeColInput.value?.name != null &&
          denseDf.col(nominalTimeColInput.value.name) !== null) ?
            nominalTimeColInput.value.name : null);
        (view.temp as Record<string, unknown>)[MEAN_SD_CTX_TAG] = {
          denseDf, timeName: meanSdNominalName ?? nmTimeName, concName: nmConcName,
          doseName: nmDoseName, nominalUsed: meanSdNominalName !== null,
        } satisfies MeanSdCtx;
        refreshMeanSdChart(/* createIfMissing */ true);
        // NONMEM mode (AC-B10 / NOTE-12): the chart consumers read the NORMALIZED
        // frame (TAD as time, EVID=1 dose rows already dropped → no FLOAT_NULL conc).
        ensureProfileChart(
          view, denseDf, nmSubjectName, nmTimeName, nmConcName,
          /* createIfMissing */ true,
          useLogY,
        );
        // Order the facets to MATCH the results-grid problem-first sort (item 2). The
        // rank rides on `outcome.results` (written by configureResultsDataFrame in the
        // swap above), keyed by subject exactly as the facet builder keys its subjects.
        const facetSubjectOrder = buildFacetSubjectOrder(outcome.results);
        // Per-subject λz acceptance verdict → each tile's bottom caption. Built from
        // the SAME results frame + thresholds as the grid's ACCEPTANCE cell, so a
        // tile's caption cannot disagree with its row's verdict.
        const facetDescriptions = buildFacetVerdictDescriptions(
          outcome.results, thresholds, controller.nonContiguousSubjects());
        // Spike `facet-lambdaz-line-points` (VR-10): per-subject λz line params, read
        // from the SAME results frame keyed by SUBJECT, for the facet terminal-line
        // overlay + point-role marking.
        const facetLambdaZ = buildFacetLambdaZParams(outcome.results);
        // Experiment: colour the facet tiles by ACCEPTANCE verdict (3-way, matches the
        // grid cell) instead of the positional palette.
        const facetColors = buildFacetVerdictColors(outcome.results, thresholds);
        // A destructive / batch design (not dense-viable) has one or two samples per animal: a
        // per-animal facet tile or λz editor has nothing to fit and only shows platform errors /
        // a lone point. The composite is the reviewed unit there (SF-3) — so neither tab is
        // docked (an open one closes: null names = the close path) (user-directed 2026-10-01).
        ensureFacetGrid(
          view, denseDf, denseViable ? nmSubjectName : null, nmTimeName, nmConcName,
          /* createIfMissing */ true,
          useLogY,
          facetSubjectOrder,
          facetDescriptions,
          facetLambdaZ,
          facetColors,
          // F6: honest used-point markers from the real λz fit membership (not the
          // time-window proxy) — the single source shared with the click-toggle.
          controller.fitRowsBySubject(),
        );
        // C0 (VR-11): after an interactive λz override, rebuild the facet from the
        // controller's IN-PLACE-patched results frame (the numeric λz columns +
        // verdict hue/caption re-derive from readAcceptance, which reads the patched
        // numerics). Captures the current chart context; overwritten each Compute.
        (view.temp as Record<string, unknown>)[FACET_REFRESH_TAG] = (subject: string): void => {
          const r = controller.results;
          if (r == null) return;
          // AC-7.3: re-derive the ACCEPTANCE verdict / flags / rank cells from the patched
          // numerics so the grid verdict cannot lag the facet caption (O(N) cell writes,
          // no viewer work — cheap).
          refreshAcceptanceColumnsInPlace(r, thresholds);
          // F7 (AC-7.1): add / clear this subject's ✏️ override flag, then force the
          // results grid to repaint — the ✏️ is drawn in the ACCEPTANCE `onCellPrepare`,
          // which a hidden-sibling change alone would not re-trigger.
          setSubjectOverrideFlag(
            r, subject, controller.overrideUserFor(subject), controller.overrideTimeFor(subject));
          {
          // eslint-disable-next-line @typescript-eslint/no-explicit-any
            const grid = (view.temp as Record<string, any>)[RESULTS_GRID_TAG] as DG.Grid | null;
            grid?.invalidate();
            // The first / last Flagged profile shows / hides the compact `Violated` column.
            if (grid != null && reasonsVisibilityStale(grid, resultsExtended(view)))
              applyResultsColumnProfile(grid, resultsExtended(view));
          }
          // Perf (VR-11): re-render ONLY the edited tile's scatterPlot — the click no
          // longer forces an N-scatterPlot re-render (the measured bottleneck, slow past
          // ~20 subjects). Data-prep still walks all rows/subjects (see refreshFacetSubject),
          // cheap at N=64. The full facet (order + all tiles) is rebuilt on the next full
          // Compute; an override never re-sorts (D10), so the tile stays put.
          refreshFacetSubject(
            view, subject, denseDf, nmSubjectName, nmTimeName, nmConcName, useLogY, r, thresholds,
            controller.fitRowsBySubject(), controller.nonContiguousSubjects());
        };
        // Bridge results-grid selection → profile-chart filter: selecting one or
        // more rows in the results grid hides the curves of unselected subjects
        // in the chart. Re-wired on every Compute because `outcome.results` is a
        // fresh DataFrame instance — the prior subscription targets a stale one.
        wireResultsSelectionToProfileFilter(view, outcome.results, nmSubjectName);
        // Same wiring for the facet grid, but applied as per-plot opacity rather
        // than a DataFrame filter (each facet plot owns its own per-subject DF).
        wireResultsSelectionToFacetDim(view, outcome.results);
        // F4 (WALKING SKELETON): dock the single-subject λz editor tab bound to the
        // results-grid current row. Reuses the facet plot machinery + the whole override
        // path; a design look — placement/arrangement is being reviewed live.
        ensureLambdaZEditor(
          view, denseDf, nmSubjectName, nmTimeName, nmConcName, useLogY, outcome.results, thresholds,
          /* dock */ denseViable);
        // Stage-scoped ribbon (D3 / D20): Settings shows on Review + Report, Group by on Report only
        // and only when the input carried grouping candidates (since v8 the candidates are the
        // carried denseDf metadata). `applyStageRibbon` rebuilds the ribbon on each stage change.
        settingsRevealed = true;
        groupByRevealed = groupCandidates.length > 0;
        // Refresh the aggregated grid against the fresh results frame (re-aggregates
        // the current Group-by selection, or closes the tab if "(none)").
        applyGrouping();
        // The dense per-subject decorators (F%, dose-proportionality, TK) all need
        // a per-subject AUCinf, which needs a fittable terminal phase (≥ 3 points).
        // A destructive / batch SPARSE design has ~1-2 samples/subject, so every
        // per-subject profile is degenerate (AUCinf NA) and those decorators are
        // analytically inapplicable — the composite Sparse panels are the correct
        // surface, so we SKIP them (and close any prior tabs) rather than spam
        // inapplicable-stratum warnings. Serial designs (many samples/subject) stay
        // viable and are unaffected (this is NOT the `isSparseDesign` pooling test,
        // which is true for serial too — it is a per-subject degeneracy test).
        // Composite mode: viability is judged on the POOLED frame (≥ minPoints per
        // cell), not the raw df (where each destructive animal has ~1-2 rows) — the
        // pooled composites are dense, so TK / F% / dose-prop decorators apply.
        // Slice D: stratum-level F% grid. Docks "Bioavailability" next to results
        // when the frame is mixed-route AND carries a per-row dose + a Sex metadata
        // column; otherwise closes any prior F% tab. Consumes the post-
        // appendGroupingColumns frame (which carries Sex / Dose).
        // In NONMEM mode the route is always present (the derived ~RouteCode
        // column). `nmDoseName` is declared with the charts row above — it feeds
        // the Mean±SD render context as well as the decorators here.
        const hasRoute = nmActive ? (routeCol !== null) : (routeColInput.value !== null);
        // Resolve the designated sex column for the dense decorators (F% / dose-
        // proportionality / TK). The decorators stratify on a CARRIED metadata
        // column, so the picked name only applies when `appendGroupingColumns`
        // actually carried it (constant-per-subject). A name that was NOT carried —
        // because the user designated a non-constant-per-subject column (AC-A5) —
        // resolves to null → F%/TK stay absent + the needs-Sex banner fires (safe
        // degradation, never a garbled partial stratification). Replaces the old
        // per-decorator inline `carried.find(/^sex$/i)`.
        const designatedSex = getSexColumnName();
        const sexForDecorators =
        (designatedSex !== null && carried.includes(designatedSex)) ? designatedSex : null;
        if (denseViable)
          applyBioavailability(view, reportable(), sexForDecorators, nmDoseName, hasRoute);
        else closeBioavailability(view);
        // Dose proportionality (Smith power model, FR-410..413 / UC-03). Docks
        // "Dose Proportionality" (per-Route×Sex β / 90% CI / verdict) + "Dose-
        // Normalized" next to results when the frame carries a per-subject dose +
        // a Sex column AND spans ≥ 2 dose levels; otherwise closes any prior tabs.
        // Consumes the post-appendGroupingColumns frame (Sex / Dose / ROUTE). The
        // P4 charts (log-log scatter + CI band, dn bar) are a design-gated follow-up.
        if (denseViable) {
          applyDoseProportionality(view, reportable(), sexForDecorators, nmDoseName, resultsNode,
            nmComposite ? 'auclast' : 'aucinf');
        } else closeDoseProportionality(view);
        let sparseResult: SparseNcaResult | null = null;
        // UC-04 sparse / destructive-sampling NCA. Docks "Sparse NCA" (closed-form
        // AUClast + Holder/Bailer CI + topology + derived + exploratory bootstrap),
        // "Composite Profile" (mean/SD/%CV/n/%BLQ per nominal time), and a VR-07
        // composite scatter when the SOURCE frame is a sparse design (a distinct
        // nominal-time column + pooled animals); otherwise closes any prior tabs.
        // Consumes the RAW observations `df` (sparse needs the per-animal samples).
        // The sparse mapping is built by `sparseMappingNow` — the SAME one the Confirm dry run uses,
        // so both compute identical composites (role overrides / declared topology / live R_acc).
        sparseResult = applySparseNca(view, srcDf, sparseMappingNow(thresholds), resultsNode, useLogY);
        // tk-aggregates AD-5/AD-6: repeat-dose TK aggregates. Docks "Accumulation"
        // / "Exposure Margin" / "Sex Comparison" next to results when the frame is
        // repeat-dose (carries a StudyDay column with ≥ 2 distinct days) + Dose +
        // Sex; otherwise closes any prior TK tabs. Consumes the post-
        // appendGroupingColumns frame (Sex / Dose / StudyDay / NOAEL / human dose).
        // The TK panels FILL-dock onto the results node returned by
        // `swapOrDockResults` (threaded through as the explicit anchor — no reliance
        // on findNode after the dock tree changed; see swapOrDockToxGrid's
        // FIRST-COMPUTE note for the null-vs-undefined guard fix).
        if (denseViable) {
          applyToxicokinetics(
            view, reportable(), carried, sexForDecorators, nmDoseName,
            denseDf, profiles, nmSubjectName, nmTimeName, resultsNode);
        } else closeToxicokinetics(view);
        // Sparse-composite review units (UC-04, SF-3): their decisions reconcile against the
        // composite's own fingerprint (members + fit summary + verdict), once the panel computed.
        const compositeUnits = sparseCompositeUnits(sparseResult);
        const sparseResets = controller.reconcileDecisions(
          new Map(compositeUnits.map((c) => [sparseUnitKey(c.groupKey), c.fingerprint])),
          (k: string) => k.startsWith(SPARSE_KEY_PREFIX));
        // The review state of this compute — read by the Review column, the stage strip, the
        // save-time stageCache and the Confirm dry run (never re-derived by a view).
        (view.temp as Record<string, unknown>)[REVIEW_CTX_TAG] = {
          thresholds, viable: denseViable, composites: compositeUnits,
          memberHashByProfile: compositeMemberHash,
        } satisfies ReviewCtx;
        refreshReviewColumns(view);
        applyReviewColumnVisibility(view);
        // F-S3.6: decisions the fit no longer supports were reset (D13) — say so, as a compute
        // NOTICE (Info), never a new WARNINGS_JSON code (the gap-18 revert lesson).
        if (denseResets + sparseResets > 0) {
          const n = denseResets + sparseResets;
          collectNotice(view, {
            key: 'REVIEW_DECISIONS_RESET', anchor: 'results', severity: 'info',
            text: `${n} review decision${n === 1 ? '' : 's'} reset — the fit changed since ` +
            `${n === 1 ? 'it was' : 'they were'} recorded.`,
          });
        }
        // AC-7.2 dependent-aggregate refresh. After an interactive λz override patches
        // a subject's derived params, the aggregates that read the results frame (By-
        // group Mean±SD/%CV, %F, dose-proportionality, TK) are stale. This closure
        // re-runs those decorators against the SAME (now-patched) `outcome.results`
        // frame, so their DATA tracks the per-subject numbers. Captured per Compute
        // (its decorator inputs are this run's), stored for the facet click flow.
        //
        // NOT refreshed: the concentration-time Mean±SD CHART — it pools the raw
        // `denseDf`, which an λz override never touches (only the derived params move).
        //
        // BANNERS refresh LIVE too, not just the numbers: the %F/dose-prop/TK decorators
        // re-`collectNotice`; instead of discarding those (which would leave a stale
        // banner — e.g. TK's Formula-B RAC_DIVERGENCE is computed FROM Day-1 λz/AUCinf, so
        // an override can newly raise or resolve it), `collectNoticesFrom` captures each
        // decorator's freshly-evaluated notices and `refreshDependentPanelBanners`
        // live-replaces each panel's banner in place — the SAME `refreshNotices` mechanism
        // By-group uses. In-place throughout: each decorator's swap-or-dock updates the
        // existing grid's frame (no re-dock / focus steal), and an override leaves every
        // structural precondition (route / dose / sex / dose-levels) intact, so they stay
        // on the swap path.
        (view.temp as Record<string, unknown>)[AGG_REFRESH_TAG] = (): void => {
          const r = controller.results;
          if (r == null) return;
          applyGrouping(); // By-group aggregate — in-place swap + its own refreshNotices banner
          if (denseViable) {
          // Same exclusion feed as the Compute pass (AD-6): a decision made since then counts.
            const rep = reportableFrameOf(view, r);
            const collected = collectNoticesFrom(view, () => {
              applyBioavailability(view, rep, sexForDecorators, nmDoseName, hasRoute);
              applyDoseProportionality(view, rep, sexForDecorators, nmDoseName, resultsNode,
                nmComposite ? 'auclast' : 'aucinf');
              applyToxicokinetics(view, rep, carried, sexForDecorators, nmDoseName,
                denseDf, profiles, nmSubjectName, nmTimeName, resultsNode);
            });
            refreshDependentPanelBanners(view, collected);
          }
        };
        // All decorators have docked (or closed) their grids — render the collected
        // notices as banners over their home panels now that the anchor roots exist.
        // (By-group renders its own banner live in `applyGrouping` — it also tracks
        // the Group-by control, which re-splits without a recompute.)
        flushNotices(view, (anchor) => resolveNoticeRoot(view, anchor));
        // LAST, after every companion table has docked into the same group: results
        // dock first by necessity, so without this the final decorator to arrive
        // owns the visible tab and the main result finishes a Compute buried.
        // A destructive design reports through the composite table — front it (user-directed
        // 2026-10-01); otherwise the per-subject results, as before.
        if (!denseViable && sparseResult !== null) activateDockedTab(view, SPARSE_PARAM_NODE_TAG);
        else activateResultsTab(view);
        consumed = true;
        refreshComputeEnabled();
        // Capture for "Add to study" (build the summary cache without recompute) +
        // enable the now-clickable Add-to-study ribbon button. `outcome.results` is
        // the post-appendGroupingColumns frame (grouping metadata carried above) —
        // the substrate buildSummaryCache aggregates over (study-unified-grid AD-1).
        lastPerProfile = outcome.perProfile;
        lastResults = outcome.results;
        // v5 (chemist-sar-grid S4): capture the SAME observations substrate the
        // mean±SD chart pools from (the MEAN_SD_CTX_TAG values set just above), so
        // the cached inline thumbnail IS the analysis's own geomean curve, just small.
        lastObservations = {
          df: denseDf, timeName: meanSdNominalName ?? nmTimeName,
          concName: nmConcName, doseName: nmDoseName,
          // Route-split key (chemist-sar-grid S4, review fix): the pipeline's own
          // per-row derived RouteCode column — present on `denseDf` whenever a route is
          // mapped, in BOTH standard mode (encoded from the mapped Route column above)
          // and NONMEM mode (encoded by the event-driven transform). The curve builder
          // reads these canonical codes directly (no re-normalization), so a compound's
          // IV and PO arms get SEPARATE curves for either ingestion mode. Null for a
          // scalar-route analysis (one route per compound → compound-level key already
          // route-correct).
          routeCodeName: denseDf.col(DERIVED_ROUTE_CODE_COL) !== null ? DERIVED_ROUTE_CODE_COL : null,
        };
        refreshSaveEnabled();
      };
      (view.temp as Record<string, unknown>)[PRESENT_TAG] = present;
      present();
      // D17: a successful Confirm moves to Review λz; a silent reopen lands on the stage where the
      // analysis stopped (F-S3.8 — Complete lands on Report).
      if (silent) {
        const st = liveWorkflowState(view).stage;
        goToStage(st.stage === 'complete' ? 'report' : (st.stage === 'review' ? 'review' : 'setup'));
        // A loaded analysis keeps the side panel on the Browse tree (user-directed 2026-10-12).
        // Hiding the toolbox at open is not enough on a large table: it came back when the
        // results were docked here (QA 2026-10-12). Re-applied now and once more on the next
        // tick, after the docked viewers' first render — the moment QA saw it return.
        keepBrowseInFront();
        setTimeout(keepBrowseInFront, 0);
      } else goToStage('review');
      // The docked results grid takes keyboard focus. A Browse-tree preview listens for this to
      // hand focus back to the tree, so arrow keys keep walking the tree (nca-tree-browser.ts).
      grok.events.fireCustomEvent(NCA_RESULTS_PRESENTED_EVENT, view);
      // AD-5 / AC-D1: hand the per-profile results + the live recipe + the grouping
      // substrate to the coordinator-driven hydration path so it can write the
      // summary cache (decision h). Math is never read back from the cache (thesis #1).
      if (studyOnComputed !== undefined) {
        studyOnComputed({
          perProfile: outcome.perProfile, recipe: getCurrentRecipe(),
          // The summary cache aggregates the REPORTABLE frame (AD-6): excluded λz values never
          // reach a study-grid statistic.
          resultsDf: reportableFrameOf(view, outcome.results),
          observations: lastObservations,
          dirty: isAnalysisDirty(view),
          workflow: workflowPersistOf(view),
        });
      }
    };

    // No per-kg confirmation gate: CL / Vz /
    // Vss are per kg on every analysis now, so a `/kg` dose without body weight is no longer the
    // exceptional basis — what it lacks (the absolute companions) is stated by PER_KG_NO_BW.
    runPipeline();
  };
  // D17: `Confirm` at the bottom of the mapping pane replaces the ribbon play icon. It runs the
  // SAME runCompute; when Setup is unchanged and results are held it just moves forward (P6 —
  // nothing recomputes). DG's raised `bigButton` — the pane's ONE primary action (user-directed
  // 2026-10-13) — dimmed (not disabled) while any check fails, so the tooltip can list what to fix
  // (C-6). Reads `Apply` after a Setup edit of a confirmed analysis (refreshStage).
  const confirmBtn = ui.bigButton('Confirm', () => onConfirm());
  confirmBtn.name = 'button-Confirm';

  // Any onChanged from any input clears `consumed` and refreshes the icon.
  // Required inputs (Subject / Time / Conc) flip ready / not-ready; optional
  // inputs always flip "consumed" back to false so the icon re-activates.
  const onAnyChange = (): void => {
    consumed = false;
    refreshComputeEnabled(); // also refreshes the stage strip (P6 lock / Revert link)
    // F3: any input / unit / rule / threshold change is an unsaved edit for a study-backed
    // analysis (no-op otherwise — a fresh analysis uses Add-to-study, not Update).
    markAnalysisDirty(view);
  };
  subjectInput.onChanged.subscribe(onAnyChange);
  timeInput.onChanged.subscribe(onAnyChange);
  concInput.onChanged.subscribe(onAnyChange);
  doseSourceInput.onChanged.subscribe(onAnyChange);
  doseColInput.onChanged.subscribe(onAnyChange);
  doseScalarInput.onChanged.subscribe(onAnyChange);
  routeSourceInput.onChanged.subscribe(onAnyChange);
  routeScalarInput.onChanged.subscribe(onAnyChange);
  routeColInput.onChanged.subscribe(onAnyChange);
  infusionSourceInput.onChanged.subscribe(onAnyChange);
  infusionDurationColInput.onChanged.subscribe(onAnyChange);
  bodyWeightSourceInput.onChanged.subscribe(onAnyChange);
  bodyWeightColInput.onChanged.subscribe(onAnyChange);
  blqFlagColInput.onChanged.subscribe(onAnyChange);
  lloqSourceInput.onChanged.subscribe(onAnyChange);
  lloqColInput.onChanged.subscribe(onAnyChange);
  lloqScalarInput.onChanged.subscribe(onAnyChange);
  bodyWeightScalarInput.onChanged.subscribe(onAnyChange);
  infusionDurationInput.onChanged.subscribe(onAnyChange);
  aucMethodInput.onChanged.subscribe(onAnyChange);
  compensatedInput.onChanged.subscribe(onAnyChange);
  blqPreInput.onChanged.subscribe(onAnyChange);
  blqEmbeddedInput.onChanged.subscribe(onAnyChange);
  blqAfterInput.onChanged.subscribe(onAnyChange);
  blqConsecutiveInput.onChanged.subscribe(onAnyChange);
  // `lambdaModeInput` is intentionally NOT subscribed — it is never
  // appended to the form (see its declaration above), so its onChanged
  // can never fire. Keeping it out of this list prevents the misleading
  // impression that toggling λz mode would re-enable Compute.
  lambdaMinPointsInput.onChanged.subscribe(onAnyChange);
  lambdaMinRSquaredInput.onChanged.subscribe(onAnyChange);
  adjRSqAcceptanceInput.onChanged.subscribe(onAnyChange);
  lambdaExcludeCmaxInput.onChanged.subscribe(onAnyChange);
  lambdaAdjRSqFactorInput.onChanged.subscribe(onAnyChange);
  extrapWarnPctInput.onChanged.subscribe(onAnyChange);
  extrapErrorPctInput.onChanged.subscribe(onAnyChange);
  // Unit inputs — load-bearing after `nca-unit-semantics`. Any change
  // must clear the "consumed" gate so Compute re-activates and the user
  // sees the rescaled output.
  timeUnitInput.onChanged.subscribe(onAnyChange);
  concUnitInput.onChanged.subscribe(onAnyChange);
  doseUnitInput.onChanged.subscribe(onAnyChange);
  bodyWeightUnitInput.onChanged.subscribe(onAnyChange);
  // Sparse-design controls (UC-04). Any change re-arms Compute so the sparse
  // panels recompute with the new mapping / declared topology.
  sparseToggleInput.onChanged.subscribe(onAnyChange);
  nominalTimeColInput.onChanged.subscribe(onAnyChange);
  animalIdColInput.onChanged.subscribe(onAnyChange);
  samplingTypeInput.onChanged.subscribe(onAnyChange);
  refreshComputeEnabled();

  // Set a control's source-toggle visibility. Writes BOTH the inline display
  // (immediate effect) AND a `data-src-hidden` flag that `appendCollapsibleSection`'s
  // expand/collapse logic reads, so a section expand does NOT clobber a
  // source-hidden control (the two mechanisms otherwise fight over inline
  // `display` on the same element — see updateWgts). `dataset.srcHidden` maps to
  // the `data-src-hidden` attribute. Declared before its callers
  // (`refreshInfusionVisibility`, `linkSourceVisibility`) to avoid a const TDZ.
  //
  // Assumption (holds by construction): every source toggle and the inputs that
  // drive `refreshInfusionVisibility` live in the SAME collapsible section as
  // the controls they show/hide, so a toggle can only fire while that section is
  // expanded. We therefore write `display` directly without consulting the
  // section's collapsed state. If a source toggle is ever moved to a different
  // section (or driven programmatically while collapsed), this would need to
  // also respect the collapse state.
  const setControlVisible = (el: HTMLElement, visible: boolean): void => {
    el.dataset.srcHidden = visible ? '0' : '1';
    el.style.display = visible ? '' : 'none';
  };

  // Sparse-section reveal (UC-04). The toggle gates the three role pickers via the
  // same `data-src-hidden` mechanism the source toggles use (so a section expand
  // does not clobber a toggled-off control). The guidance hint (AC3) shows only in
  // the actionable state: toggle on, but no nominal-time column mapped.
  const refreshSparseVisibility = (): void => {
    const on = sparseToggleInput.value ?? false;
    setControlVisible(nominalTimeColInput.root, on);
    setControlVisible(animalIdColInput.root, on);
    setControlVisible(samplingTypeInput.root, on);
    const needHint = on && nominalTimeColInput.value === null;
    setControlVisible(sparseHint, needHint);
  };
  sparseToggleInput.onChanged.subscribe(() => refreshSparseVisibility());
  nominalTimeColInput.onChanged.subscribe(() => refreshSparseVisibility());

  // Rebuilds the per-value → M/F/U mapping inputs whenever the designated sex
  // column changes (Slice B). Scans the WHOLE column (all rows, IV and EV) so no
  // token is missed and the IV↔EV value-equality join cannot break. > 8 distinct
  // values ⇒ implausible sex axis (disable Compute); ≤ 8 ⇒ one choice per value,
  // pre-filled from the honest-uncertainty heuristic (text tokens resolved,
  // ambiguous codes left unassigned). Inputs are rebuilt from scratch each time
  // (no per-column cache — architect SIMPLIFY).
  const refreshSexMapping = (): void => {
    for (const s of sexValueSubs) s.unsubscribe();
    sexValueSubs.length = 0;
    sexValueInputs.length = 0;
    sexMapContainer.innerHTML = '';
    sexTooMany = false;
    sexUnassigned = false;
    const col = sexColInput.value;
    if (col === null) {setControlVisible(sexMapContainer, false); refreshComputeEnabled(); return;}

    const uniques = uniqueSexValues(col);
    if (uniques.length > SEX_MAX_DISTINCT) {
      // Too many categories — not a sex axis. Disable Compute (the `sexTooMany`
      // flag is read by refreshComputeEnabled → greys + tooltips the icon). The
      // user-facing explanation is the INLINE picker validator (red), not a panel
      // text line — build no per-value inputs and keep the container empty.
      sexTooMany = true;
      sexColInput.validate(); // surface the >8 validator red immediately
      setControlVisible(sexMapContainer, false);
      refreshComputeEnabled();
      return;
    }

    const syncUnassigned = (): void => {
      sexUnassigned = sexValueInputs.some(({input}) => input.value === SEX_UNASSIGNED);
    };
    for (const raw of uniques) {
      const def = defaultSexAssignment(raw);
      const label = raw === '' ? '(blank)' : raw;
      // AD-5 replay precedence: a pending recipe map wins over the name heuristic
      // so a saved Sex designation round-trips exactly (AC-D6).
      const replayed = pendingSexMap?.[raw];
      const initial = (replayed === 'M' || replayed === 'F' || replayed === 'U') ?
        replayed : (def ?? SEX_UNASSIGNED);
      const input = ui.input.choice<string>(label, {
        items: [SEX_UNASSIGNED, 'M', 'F', 'U'],
        value: initial,
      });
      bindInputTip(input, richTip({
        title: `Sex value "${label}"`,
        summary: 'Map this raw value to a sex stratum (M / F / U).',
        detail: def === null ?
          'Ambiguous code — assign it from the study’s documented convention. %F ' +
          'cannot be computed until every value is mapped.' :
          'Pre-filled from an unambiguous token; change it if your data differs.',
      }));
      // Unassigned value → INLINE red validator on this input (no panel text). The
      // `sexUnassigned` Compute-disable flag is the gate; this is the per-value
      // visual marking exactly which value still needs an M/F/U.
      input.addValidator((v) => v === SEX_UNASSIGNED ? 'Assign M / F / U' : null);
      sexValueSubs.push(input.onChanged.subscribe(() => {syncUnassigned(); input.validate(); onAnyChange();}));
      sexValueInputs.push({raw, input});
      sexMapContainer.appendChild(input.root);
      input.validate(); // surface the unassigned red immediately on build (coded 1/2)
    }
    syncUnassigned();
    setControlVisible(sexMapContainer, true);
    refreshComputeEnabled();
  };
  // On column designation (manual pick OR the autodetect auto-bind below) rebuild
  // the per-value → M/F/U inputs and re-arm Compute.
  sexColInput.onChanged.subscribe(() => {refreshSexMapping(); onAnyChange();});

  // Infusion-duration controls have TWO visibility gates that AND together
  // (route relevance × Pattern-A source). `infusionDuration` is consumed only
  // by the IV-infusion route, so the whole infusion group is hidden otherwise;
  // within the relevant group, `infusionSourceInput` selects scalar XOR column.
  // Bespoke (not `linkSourceVisibility`) because that helper unconditionally
  // shows one control, whereas here the group can be hidden entirely.
  // (Replaces the former scalar-only `refreshInfusionVisibility` and the
  // route-scalar-disable `refreshRouteSource`; route scalar/column visibility
  // is now Pattern A via `linkSourceVisibility(routeSourceInput, …)` below.)
  // `forcedRelevant` lets the autodetect block assert "a route column is mapped"
  // without re-reading `routeSourceInput.value` — which, right after a
  // programmatic `writeSourceMode(routeSourceInput, 'column')`, still reads the stale
  // pre-autodetect 'value' at mount (same getter-staleness as UX-01). Without
  // it the initial call took the scalar branch ('PO' → not relevant) and hid the
  // infusion controls for an autodetected route column (PATTERN-GAP-UX01).
  const refreshInfusionVisibility = (forcedRelevant?: boolean): void => {
    // Layer 1 — relevance. 'column' route: relevant whenever a route column is
    // mapped (per-row route is unknown at form time; preserves the bare-IV
    // disambiguation affordance). Scalar route: only the IV-infusion value.
    const relevant = forcedRelevant ?? (readSourceMode(routeSourceInput) === 'column' ?
      (routeColInput.value !== null) :
      (routeScalarInput.value === 'IV-infusion'));
    if (!relevant) {
      setControlVisible(infusionSourceInput.root, false);
      setControlVisible(infusionDurationInput.root, false);
      setControlVisible(infusionDurationColInput.root, false);
      return;
    }
    // Layer 2 — Pattern-A source toggle within the relevant group. Clear the
    // column BEFORE the visibility writes (mirrors `linkSourceVisibility`) so no
    // observer ever sees a hidden-but-populated column.
    const uniform = readSourceMode(infusionSourceInput) === 'value';
    if (uniform && infusionDurationColInput.value !== null)
      infusionDurationColInput.value = null; // keep buildContract unambiguous
    setControlVisible(infusionSourceInput.root, true);
    setControlVisible(infusionDurationInput.root, uniform);
    setControlVisible(infusionDurationColInput.root, !uniform);
  };
  // Wrapped in `() => …` so the onChanged-emitted value is NOT forwarded as
  // `forcedRelevant` (a truthy payload would wrongly force the relevant branch).
  // User-driven changes settle the getters, so the unforced read is correct.
  routeScalarInput.onChanged.subscribe(() => refreshInfusionVisibility());
  routeSourceInput.onChanged.subscribe(() => refreshInfusionVisibility());
  routeColInput.onChanged.subscribe(() => refreshInfusionVisibility());
  infusionSourceInput.onChanged.subscribe(() => refreshInfusionVisibility());
  // Initial call is wired BEFORE autodetect (below, next to the
  // `linkSourceVisibility` calls) so it runs on the genuine default route state
  // (scalar 'PO' → not relevant → infusion hidden) and the autodetect route
  // block force-shows via `refreshInfusionVisibility(true)`.

  // Hide the scalar input when its companion column input is set — the scalar
  // exists only to satisfy `ProfileScalars` when the column path is unused.
  // Generic over the scalar-input value type `S` so it serves both the numeric
  // scalars (Dose / LLOQ / Body weight / Infusion) and the Route scalar, which
  // is a `choice<string>`. The body only toggles visibility and clears the
  // column on switch-back to value — it never reads the scalar value — so
  // widening `S` is safe.
  // Returns `apply` so the autodetect block can FORCE column-mode visibility
  // right after it programmatically flips the source to `'column'`, without
  // re-reading the choice getter — which does not reflect a just-assigned value
  // synchronously at mount (UX-01: the stale read left the scalar showing for
  // autodetected-to-column groups). `apply(forcedUniform)` bypasses the getter.
  // `onColumnRestored` (optional): fired when a stashed column is re-applied on
  // switch-back to column mode. The restore is a programmatic set, which does NOT
  // fire `colInput.onChanged` (DG model, see UX-01), so any listener that must
  // react to the column reappearing is notified explicitly. Used by the route
  // group so the infusion controls re-show on a value→column round-trip
  // (PATTERN-GAP-UX01: refreshInfusionVisibility fires before this restore and
  // would otherwise read the not-yet-restored null column).
  const linkSourceVisibility = <S>(
    sourceInput: DG.InputBase<string | null>,
    scalarInput: DG.InputBase<S>,
    colInput: DG.InputBase<DG.Column | null>,
    onColumnRestored?: () => void,
  ): ((forcedUniform?: boolean) => void) => {
    // Remembers the last picked column so a value→column round-trip restores it
    // instead of leaving the picker empty. The invariant `colInput.value === null
    // ⇔ source is 'value'` is preserved (every Run-time consumer keys
    // dose/route/lloq/BW column-vs-scalar off `<x>ColInput.value !== null`, e.g.
    // buildContract + the route normalizer), so in value mode the column is still
    // cleared — but the selection is stashed and re-applied on switch-back.
    let lastCol: DG.Column | null = null;
    const apply = (forcedUniform?: boolean): void => {
      const uniform = forcedUniform ?? (readSourceMode(sourceInput) === 'value');
      if (uniform) {
        // Switching to value: stash + clear the column so the contract builder
        // unambiguously falls through to the scalar path.
        if (colInput.value !== null) {
          lastCol = colInput.value;
          colInput.value = null;
        }
      } else if (colInput.value === null && lastCol !== null) {
        // Switching back to column: restore the previously picked column.
        // (If `df` were swapped out between the clear and this restore, `lastCol`
        // would be a stale reference — but this pane is rebuilt per DataFrame, so
        // a swap re-runs the whole mount and resets `lastCol`. Same stale-mapping
        // failure mode as any other column input; not worsened here.)
        colInput.value = lastCol;
        onColumnRestored?.();
      }
      setControlVisible(scalarInput.root, uniform);
      setControlVisible(colInput.root, !uniform);
    };
    sourceInput.onChanged.subscribe(() => apply());
    apply();
    return apply;
  };
  // Column auto-fill (spec §3.1, decisions 1A/2b/3a/4a). Runs once at initial
  // mount; user overrides persist for the view's lifetime. Silent UX — same
  // affordance as a manual pick.
  //
  // Ordering: autodetect runs AFTER all `ui.input.column(...)` constructors
  // and AFTER the onAnyChange subscriptions (so the assignments below trigger
  // refreshComputeEnabled via the existing onChange wiring).
  //
  // Source-toggle visibility is wired BEFORE autodetect (UX-01 fix). Wiring it
  // after used to leave each group's initial `apply()` reading the choice getter
  // to discover the post-autodetect source — but a programmatically-assigned
  // `writeSourceMode(<group>SourceInput, 'column')` is NOT reflected by that getter
  // synchronously at mount, so `apply()` read the stale default `'value'` and
  // left the scalar showing next to the active column picker (Ed U2/U3, dose +
  // route). Now: wire first (establishes the value-mode default + subscribes
  // onChanged for later manual toggles), then the autodetect block calls the
  // captured `apply<X>Vis(false)` after each `'column'` flip to FORCE column-mode
  // visibility without re-reading the unsettled getter.
  const applyDoseVis = linkSourceVisibility(doseSourceInput, doseScalarInput, doseColInput);
  const applyLloqVis = linkSourceVisibility(lloqSourceInput, lloqScalarInput, lloqColInput);
  const applyBwVis = linkSourceVisibility(bodyWeightSourceInput, bodyWeightScalarInput, bodyWeightColInput);
  // Route scalar is a choice<string> (handled by the generic linkSourceVisibility).
  // A restored route column is always infusion-relevant, so force-show the
  // infusion group on restore (refreshInfusionVisibility runs before the restore
  // on the same onChanged and would otherwise read the still-null column).
  const applyRouteVis = linkSourceVisibility(
    routeSourceInput, routeScalarInput, routeColInput,
    () => refreshInfusionVisibility(true));
  // Initial infusion visibility on the genuine default route state (before any
  // autodetect flip) — hides the group for the default scalar 'PO' route.
  refreshInfusionVisibility();

  const detected = detectMapping(df);
  if (detected.subject !== undefined) subjectInput.value = detected.subject;
  if (detected.time !== undefined) timeInput.value = detected.time;
  if (detected.concentration !== undefined) concInput.value = detected.concentration;
  // AC-A7: in Event-driven mode the transform owns AMT (`'amt'` ∈ DOSE_SYNONYMS),
  // so suppress the standard AMT→dose autofill — the dose slot is populated only
  // from the transform output (the emitted Dose_mg column).
  if (detected.dose !== undefined && !eventDetection.isEventDriven) {
    doseColInput.value = detected.dose;
    writeSourceMode(doseSourceInput, 'column');
    applyDoseVis(false); // force column-mode visibility (getter not yet settled)
  }
  if (detected.bodyWeight !== undefined) {
    bodyWeightColInput.value = detected.bodyWeight;
    writeSourceMode(bodyWeightSourceInput, 'column');
    applyBwVis(false);
  }
  if (detected.blqFlag !== undefined) blqFlagColInput.value = detected.blqFlag;
  if (detected.route !== undefined) {
    routeColInput.value = detected.route;
    writeSourceMode(routeSourceInput, 'column');
    applyRouteVis(false);
    // Force-show the infusion group: a mapped route column is always
    // infusion-relevant ('column' branch), but the stale `routeSourceInput`
    // getter would mis-read 'value' here (PATTERN-GAP-UX01).
    refreshInfusionVisibility(true);
  }
  // (No infusion autodetect — DetectedMapping has no infusion field, so
  // infusionSourceInput correctly stays 'value'.)
  if (detected.lloq !== undefined) {
    lloqColInput.value = detected.lloq;
    writeSourceMode(lloqSourceInput, 'column');
    applyLloqVis(false);
  }

  // Unit autofill (GAP-11): when the detector parsed an unambiguous unit token
  // from a unit-suffixed column name (e.g. `Conc_ng_mL` → 'ng/mL'), pre-select
  // the matching dropdown so a scientist who maps the column but forgets to
  // switch the unit does not silently mis-scale CL/Vz/AUC by 1000×. Guarded by
  // the live CHOICES arrays so a detector value outside the current dropdown
  // set is ignored rather than assigned. Manual override is preserved — these
  // inputs stay editable.
  if (detected.concentrationUnit !== undefined &&
      CONC_UNIT_CHOICES.includes(detected.concentrationUnit as ConcentrationUnit))
    concUnitInput.value = detected.concentrationUnit;
  if (detected.doseUnit !== undefined &&
      DOSE_UNIT_CHOICES.includes(detected.doseUnit as DoseUnit))
    doseUnitInput.value = detected.doseUnit;
  if (detected.timeUnit !== undefined &&
      TIME_UNIT_CHOICES.includes(detected.timeUnit as TimeUnit))
    timeUnitInput.value = detected.timeUnit;
  if (detected.bodyWeightUnit !== undefined &&
      BODY_WEIGHT_UNIT_CHOICES.includes(detected.bodyWeightUnit as BodyWeightUnit))
    bodyWeightUnitInput.value = detected.bodyWeightUnit;

  // Sparse-design pre-fill (UC-04, zero-regression anchor). Pre-fill the role
  // pickers from the SAME name heuristic `applySparseNca` uses, and AUTO-CHECK the
  // toggle when a nominal-time column is present — so a conventionally-named sparse
  // dataset (e.g. `04_mouse_sparse_destructive.csv`: NominalTime_hr / MouseID /
  // SamplingType) reaches the picked-value path with byte-identical output, no
  // clicks. Animal ID defaults to the mapped Subject (the heuristic's own default),
  // falling back to the animal/mouse/subject name pattern. Sampling type stays
  // `(auto)` = derive-from-data. User edits override all of this for the view's life.
  const nominalGuess = findColByPattern(df, /nominal.*time/i);
  if (nominalGuess !== null) {
    nominalTimeColInput.value = df.col(nominalGuess);
    sparseToggleInput.value = true;
  }
  animalIdColInput.value = subjectInput.value ??
    (() => {const a = findColByPattern(df, /animal|mouse|subject/i); return a !== null ? df.col(a) : null;})();
  refreshSparseVisibility();

  // Event-driven (NONMEM) mode (AD-2). Lock the dose unit to absolute mg (AC-A6
  // — NONMEM AMT is absolute mass; leaving the pane at mg/kg would re-multiply by
  // BW in scaleProfileInputs → ~4× CL inflation) and surface the proposal hint.
  const refreshNonmemMode = (): void => {
    const on = nonmemToggle.value ?? false;
    // The mode is only ACTIVE when the layout is actually NONMEM (the detection
    // gate) — the SAME condition Compute uses for `nmActive`. Gating the dose-unit
    // lock on `on` alone would lock mg on a non-NONMEM table the user manually
    // toggled on, while Compute silently runs standard mode → wrong CL if the real
    // dose is mg/kg. Keep the UI side effect consistent with what Compute does.
    const active = on && eventDetection.isEventDriven;
    // Route hint visibility through setControlVisible (sets BOTH style.display AND
    // data-src-hidden), so a Profile-data section expand does not clobber it back
    // to visible — the same mechanism the sparse hint uses. A raw style.display
    // here is overwritten by appendCollapsibleSection's updateWgts on expand,
    // surfacing the "NONMEM detected" hint on non-NONMEM tables (e.g. theoph).
    setControlVisible(nonmemHint, eventDetection.isEventDriven && !on);
    if (active) doseUnitInput.value = 'mg'; // syncs to doseUnitInputV via linkUnitPair
    doseUnitInput.enabled = !active; // locked read-only only when NONMEM mode is active
    doseUnitInputV.enabled = !active; // dual-mode mirror carries the same lock
  };
  // Propose the mode when the NONMEM layout is detected (never silent, AC-A5):
  // pre-check the toggle; the scientist confirms by leaving it on through Compute.
  if (eventDetection.isEventDriven) nonmemToggle.value = true;
  nonmemToggle.onChanged.subscribe(() => {confirmedRoute = null; refreshNonmemMode(); onAnyChange();});
  refreshNonmemMode();

  // Concentration-time profile chart + facet grid are built ONLY by the
  // Compute handler (`createIfMissing = true`); they are NOT refreshed
  // reactively on mapping changes. The docked viewers reflect the LAST RUN —
  // editing Subject/Time/Concentration after a Compute leaves the charts in
  // place until the user clicks Confirm again. (Reactive rebuilds both
  // surprised the user — the chart should be a result of RUN, not a live
  // preview — and crashed on transient invalid mappings, e.g. Time momentarily
  // pointed at the Concentration column; see `chartColumnsDistinct`.) The
  // Y-axis log/linear toggle is a DISPLAY change, not a mapping change, so it
  // still updates live via `applyYAxisModeToCharts` (setOptions, no rebuild).
  // R1 C2: re-evaluate the Compute button after the autodetect block in case
  // Datagrok's `input.value` setter is asynchronous (the onAnyChange listeners
  // above would not have fired by now). If the setter is synchronous, this
  // call is idempotent.
  refreshComputeEnabled();

  // (Source-toggle visibility AND infusion visibility were wired + force-applied
  // above, around autodetect — no post-autodetect getter read-back, which would
  // be stale at mount. See UX-01 / PATTERN-GAP-UX01.)

  // Single `ui.form` with section headers + inputs interleaved (DiffStudio
  // pattern — see DiffStudio/src/app.ts:1750-1769 `getInputsForm`). The form
  // is created empty, then headers (`ui.h2`) and `input.root` elements are
  // appended in order. CSS rule `.nca-mapping-pane h2 { grid-column: 1 / -1 }`
  // (app-styles.css) makes section headers span the full form grid width.
  const form = ui.form([]);
  // Three top-level collapsible sections. Data + Optional default open;
  // Settings default collapsed (low-frequency control vs the always-needed
  // mapping + scalars). Sub-sections inside Settings (AUC / BLQ / λz /
  // Extrapolation) are not individually collapsible — collapsing Settings
  // hides all of them at once. `lambdaModeInput` is intentionally absent
  // from the form — see its declaration for the no-effect-yet rationale.
  // Time and Concentration carry their unit picker inline (`unitRow`) — the
  // standalone `Units` section was retired here (mapping-pane-inline-units,
  // 2026-07-28), because a unit read five rows away from the value it scales is
  // a unit the scientist does not read at all.
  const reapplyProfileData = appendCollapsibleSection(form, 'Profile data',
    [subjectInput.root,
      unitRow([timeInput], timeUnitInput),
      timeUnitInputV.root, // narrow-mode row under Time (dual-mode; hidden in wide)
      unitRow([concInput], concUnitInput),
      concUnitInputV.root, // narrow-mode row under Concentration
      nonmemToggle.root, nonmemHint],
    /* initialExpanded */ true,
    richTip({
      title: 'Profile data',
      summary: 'Required column mappings: subject identifier, sampling time, and measured concentration.',
      detail: 'Compute stays disabled until all three are set. Time and Concentration ' +
        'carry their source unit alongside; values are rescaled to canonical h / mg/L ' +
        'before computation.',
    }));
  const reapplyContext = appendCollapsibleSection(form, 'Context', [
    // Each group is `type toggle, scalar, column` (Pattern A — exactly one of
    // scalar/column is visible per the toggle). Visibility is display-driven,
    // so all roots are appended unconditionally.
    //
    // Dose and Body weight additionally fold their unit picker onto the value's
    // own row via `unitRow`. Both Pattern-A controls go INSIDE the wrapper, so
    // `linkSourceVisibility` keeps toggling the inner roots exactly as before and
    // the unit rides next to whichever of the two is currently shown. `Dose  [1]`
    // read on its own could not distinguish 1 mg from 1 mg/kg — a factor of body
    // weight straight into CL and Vz.
    doseSourceInput.root,
    unitRow([doseScalarInput, doseColInput], doseUnitInput),
    doseUnitInputV.root, // narrow-mode row under Dose (dual-mode; hidden in wide)
    routeSourceInput.root, routeScalarInput.root, routeColInput.root,
    infusionSourceInput.root, infusionDurationInput.root, infusionDurationColInput.root,
    bodyWeightSourceInput.root,
    unitRow([bodyWeightScalarInput, bodyWeightColInput], bodyWeightUnitInput),
    bodyWeightUnitInputV.root, // narrow-mode row under Body weight
    blqFlagColInput.root,
    lloqSourceInput.root, lloqScalarInput.root, lloqColInput.root,
  ], /* initialExpanded */ true,
  richTip({
    title: 'Context',
    summary: 'Dosing context and per-row qualifiers: dose, route, infusion duration ' +
      '(IV-infusion only), body weight, LLOQ — plus the BLQ flag.',
    detail: 'Each is set as a uniform value OR a per-row column via its type toggle. ' +
      'Drives derived parameters (CL, Vz, weight normalisation) and BLQ handling. ' +
      'Dose and Body weight carry their source unit alongside the value; LLOQ is ' +
      'expressed in the Concentration unit and is rescaled by the same factor.',
  }));
  // Subject Sex section (A1) — collapsible, DEFAULT EXPANDED, directly under
  // Context. Always visible (not hidden/collapsed) because sex stratification is
  // science-critical: a silently-unstratified F% would otherwise go unnoticed.
  // The picker designates the sex column; `sexMapContainer` holds the per-value →
  // M/F/U inputs (Slice B), populated on designation.
  appendCollapsibleSection(form, 'Subject Sex', [
    sexColInput.root, sexMapContainer,
  ], /* initialExpanded */ true,
  richTip({
    title: 'Subject Sex',
    summary: 'Designate which column holds subject sex, and map its values to M/F/U.',
    detail: 'Used to stratify bioavailability (%F), dose-proportionality and TK by sex. ' +
      'Leave the picker empty if your data has no sex axis.',
  }));
  refreshSexMapping(); // initial empty render (no candidate / pre-auto-bind)
  // Autodetect → AUTO-BIND (user-directed; supersedes the suggest-only A2 design,
  // peer-review finding #2). When a synonym-matching column exists, set it as the
  // Sex column so the value→M/F/U mapping inputs surface immediately for review.
  // Deferred a tick so the picker is fully initialized and its `onChanged` fires
  // (which builds the dependent per-value inputs). No candidate → picker stays
  // empty. Residual risk (a synonym-named non-sex column auto-binding) is bounded:
  // the mapping inputs are shown for review and coded values block Compute until
  // assigned (`sexUnassigned`).
  // AD-5: skip the heuristic auto-bind when replaying a recipe — `applyRecipe`
  // (below) sets the Sex column explicitly, and this deferred setter would
  // otherwise clobber it ~25ms later.
  if (opts.recipe === undefined) {
    const sexCandidate = findColByPattern(df, SEX_SYNONYMS);
    if (sexCandidate !== null)
      setTimeout(() => {sexColInput.value = df.col(sexCandidate);}, 25);
  }
  // Sparse / destructive design section (UC-04) — collapsible, default COLLAPSED
  // (only relevant to sparse studies; the dense majority should not pay viewport
  // for it). Sits AFTER Units, BEFORE Rules (per user request 2026-06-16). The
  // toggle is always visible when the section is expanded; the three role pickers
  // + the guidance hint are gated by it (data-src-hidden, set in
  // refreshSparseVisibility). Pre-fill + auto-on already ran in the autodetect
  // block, so on a conventionally-named sparse dataset the toggle is already
  // checked here.
  appendCollapsibleSection(form, 'Sparse / destructive design', [
    sparseToggleInput.root,
    nominalTimeColInput.root, animalIdColInput.root, samplingTypeInput.root,
    sparseHint,
  ], /* initialExpanded */ false,
  richTip({
    title: 'Sparse / destructive design',
    summary: 'Map the nominal-time and animal-ID roles and optionally declare the ' +
      'sampling topology (destructive / batch / serial).',
    detail: 'Enable the toggle to override the column-name heuristic; auto-enabled when a ' +
      'nominal-time column is detected.',
    note: 'Drives the composite-profile sparse NCA (closed-form AUC + Holder/Bailer CI) ' +
      'shown in the Sparse panels.',
  }));
  appendCollapsibleSection(form, 'Rules', [
    ...buildSubsection('AUC', [aucMethodInput, compensatedInput],
      richTip({
        title: 'AUC',
        summary: 'Area-under-curve integration.',
        detail: 'Which trapezoidal scheme integrates concentration over time, and whether ' +
          'the running sum uses compensated arithmetic to reduce floating-point drift.',
      })),
    ...buildSubsection('BLQ', [
      blqPreInput, blqEmbeddedInput, blqAfterInput, blqConsecutiveInput,
    ],
    richTip({
      title: 'BLQ',
      summary: 'Below-LOQ point handling.',
      detail: 'Action applied to BLQ samples in each phase of the profile: before the first ' +
        'measurable observation, between measurables, and after the last measurable.',
    })),
    ...buildSubsection('λz', [
      lambdaMinPointsInput, lambdaMinRSquaredInput, adjRSqAcceptanceInput,
      lambdaExcludeCmaxInput, lambdaAdjRSqFactorInput,
    ],
    richTip({
      title: 'λz',
      summary: 'Terminal-phase slope (λz) auto-best-fit criteria.',
      detail: 'Minimum points in the candidate window, minimum adjusted R² to accept the fit, ' +
        'whether Cmax is excluded, and the tie-break factor between candidate windows of ' +
        'comparable quality.',
    })),
    ...buildSubsection('Extrapolation', [
      extrapWarnPctInput, extrapErrorPctInput,
    ],
    richTip({
      title: 'Extrapolation',
      summary: 'Thresholds on the % of AUCinf that comes from extrapolation beyond the last ' +
        'measured time point.',
      detail: 'Above the warning threshold AUCinf is flagged; above the error threshold it is ' +
        'considered unreliable.',
    })),
  ], /* initialExpanded */ false,
  richTip({
    title: 'Rules',
    summary: 'NCA method rules — how the analysis runs, not what data goes in.',
    detail: 'AUC integration scheme, BLQ-point handling, λz selection criteria, and ' +
      'extrapolation thresholds.',
    note: 'Defaults are sensible for most preclinical PK profiles; expand to override.',
  }));

  // `Discard` (P6 / F-S3.4; `Revert` until 2026-10-13, user-directed): re-applies the confirmed
  // recipe after a Setup edit, so an abandoned edit never loses the review work it would have
  // reset. Shown only then, on any stage. A regular button LEFT of Apply (user-directed 2026-10-01).
  const revertBtn = ui.button('Discard', () => {
    const c = getConfirmed(view);
    if (c === null) return;
    applyRecipe(c.recipe);
    // DG input setters can settle asynchronously — re-read the stage once they have.
    setTimeout(() => refreshStage(), 50);
  }, 'Discard the changes and restore the setup you last confirmed');
  revertBtn.style.display = 'none'; // shown only after a Setup edit of a confirmed analysis
  // Discard + Confirm/Apply live in their OWN `ui.form` (its `addButtons` row, right-aligned like
  // the Review pane's reason form), PINNED under the scrolling mapping form (user-directed
  // 2026-10-13). A sticky row inside the mapping form vanished in its narrow (condensed) layout;
  // a separate form outside the scroll area is always on screen, whatever the form's layout.
  const actionForm = ui.form([]);
  actionForm.classList.add('nca-mapping-actions');
  ui.forms.addButtons(actionForm, [revertBtn, confirmBtn]);
  const btnRow = actionForm.querySelector<HTMLElement>('.ui-input-buttons > .ui-input-editor');
  if (btnRow !== null) btnRow.style.justifyContent = 'flex-end';
  const paneRoot = ui.divV([
    ui.divV([form, errorPanel], 'nca-mapping-pane-scroll'),
    actionForm,
  ], 'nca-mapping-pane');

  // Default dock ratio (mapping-pane-unit-dual-mode, user-directed): 0.25 ≈ 480px on
  // a 1920 screen / 341px on a 1366 laptop. On a full-size screen this opens in the
  // NON-condensed (label-left) layout, so the unit pickers show INLINE on their value
  // rows from the start; on a 1366 laptop 341px sits right at the ~346px condensed
  // threshold, so it may open in the narrow (label-on-top) layout — either way the
  // dual-mode pickers render correctly (inline when wide, their own vertical rows when
  // narrow). The divider stays draggable.
  //
  // History: 0.25 is the ORIGINAL default, from before UX-02 narrowed it to 0.18 (Ed U1
  // — "pane too wide, dragged narrower on every dataset"). The dual-mode pickers changed
  // that calculus (a narrow pane is no longer broken), and the user restored 0.25 as the
  // balance between revealing inline units and not eating grid width. (0.33 was tried
  // this session and dialled back to 0.25.)
  view.dockManager.dock(paneRoot, 'left', null, 'Column mapping', 0.25);

  // ── Stage strip + stage navigation (nca-workflow-stages F-S3.1 / F-S3.2 / F-S3.4) ──────────
  const strip = createWorkflowStrip((st) => goToStage(st), () => startReview(view), () => {
    const b = studyBindingOf(view); // study bindings only (standalone → null)
    if (b === null) return;
    void openStudyOverview(b.studyId);
  });
  view.root.prepend(strip.root);
  // Setup / Review λz / Report (D7). Starts at Setup; a compute or a reopen moves it.
  let viewStage: ViewStage = 'setup';
  // Assigned once the Add-to-study button exists (declared much later in this closure — a `let`
  // hook avoids a temporal-dead-zone read when refreshStage runs during mount).
  let addBtnRefresh: (() => void) | null = null;
  const refreshStage = (): void => {
    const st = liveWorkflowState(view).stage;
    const b = (view.temp as Record<string, unknown>)[STUDY_BINDING_TAG] as StudyBinding | null;
    const studyLabel = b?.container === 'study' ? (b.containerLabel ?? null) : null;
    strip.update({stage: st, viewStage, studyLabel, confirmable: confirmReady});
    addBtnRefresh?.();
    // Discard pairs with Apply on EVERY stage (user-directed 2026-10-13): the mapping pane stays
    // docked across stages, so an edit made while on Review / Report must be discardable there too.
    revertBtn.style.display = st.setupChanged ? '' : 'none';
    // Back on Setup, the first edit of a confirmed analysis turns Confirm into Apply (user-directed
    // 2026-10-13); Discard (or a fresh compute) turns it back.
    confirmBtn.textContent = st.setupChanged ? 'Apply' : 'Confirm';
    // A Setup edit locks Review/Report (P6); if the analyst is looking at one of them when the
    // lock arrives (e.g. an edit from the tour), they stay put — the strip shows the lock.
  };
  // Stage-scoped widgets (D3 / D20): show exactly the shown stage's widgets; hide the rest.
  // Ribbon layout (user-directed 2026-09-29): the PERSISTENT controls — SAVE (leftmost, D2) and the
  // guide — form the first group and never move. The view tools follow as ONE group: Settings (once
  // there is something to present — after the first compute, off Setup), then Group by on Report.
  // Report's labelled Add to study closes the ribbon. Appending at the end means a stage change never
  // shifts a control that was already on screen (stable positions), while D3's hide-not-dim holds.
  const applyStageRibbon = (): void => {
    const tools: HTMLElement[] = [];
    if (viewStage !== 'setup' && settingsRevealed) tools.push(settingsBtn);
    if (viewStage === 'report' && groupByRevealed) tools.push(groupByInput.root);
    const groups: HTMLElement[][] = [[updateBtn, infoBtn]];
    if (tools.length > 0) groups.push(tools);
    if (viewStage === 'report') groups.push([saveBtn]); // ⊞ Add to study — Report only (D3)
    view.setRibbonPanels(groups);
  };
  // Stage navigation is INDICATOR-ONLY (user-directed 2026-09-29): it never closes or re-docks a
  // panel. The mapping pane and the results/charts stay docked across every stage — closing and
  // re-presenting them dropped input state and left duplicate chart tabs (QA 2026-09-29). A stage
  // change only updates the strip, the stage-scoped ribbon actions, and the Review-column
  // visibility. The results are presented ONCE per compute (`present()` in runPipeline).
  const goToStage = (target: ViewStage): void => {
    const st = liveWorkflowState(view).stage;
    if (target !== 'setup' && !reachableStages(st).has(target)) {
      // Review clicked on Setup = Confirm (user-directed 2026-10-13): the same action as the button,
      // pre-warning included; a successful compute moves to Review itself.
      if (target === 'review' && confirmReady) onConfirm();
      return;
    }
    viewStage = target;
    (view.temp as Record<string, unknown>)[VIEW_STAGE_TAG] = target;
    applyStageRibbon();
    refreshStage();
    applyReviewColumnVisibility(view);
    if (target === 'report') showReportPanel(view);
  };
  (view.temp as Record<string, unknown>)[GOTO_STAGE_TAG] = goToStage;
  (view.temp as Record<string, unknown>)[STAGE_REFRESH_TAG] = (): void => {
    refreshReviewColumns(view);
    refreshStage();
  };
  // Confirm (D17): compute the current Setup. `consumed` (nothing changed since the last compute)
  // makes it a no-op — the button is dimmed, and moving to Review/Report is done by clicking the
  // stage strip, not by re-confirming (user-directed 2026-09-29). A dry-run pre-warning fires
  // first when the recipe changed AND review decisions would be reset (F-S3.5 / D13).
  const onConfirm = (): void => {
    if (consumed) return; // nothing to recompute; navigate via the strip
    const st = liveWorkflowState(view).stage;
    const controller = (view.temp as Record<string, unknown>)[NCA_CONTROLLER_TAG] as NcaController | null;
    const hasDecisions = (controller?.getReviewLedgers() ?? []).some((l) => l.current !== null);
    if (!hasDecisions || !st.setupChanged) {runCompute(); return;}
    runCompute({dryRun: (n) => {
      if (n === 0) {runCompute(); return;}
      // F-S3.5: "N review decisions will be reset" → Continue / Cancel.
      const dlg = ui.dialog(`${n} review decision${n === 1 ? '' : 's'} will be reset`)
        .add(ui.divText(`The fit changes for ${n === 1 ? 'that profile' : 'those profiles'}. The reset ` +
          'decisions stay in the history.', {style: {maxWidth: '420px', whiteSpace: 'normal'}}))
        .onOK(() => runCompute());
      dlg.show();
      const ok = dlg.root.querySelector('button[name="button-OK"]') as HTMLElement | null;
      if (ok !== null) ok.textContent = 'Continue';
    }});
  };

  // --- Dual-mode unit-picker visibility (mapping-pane-unit-dual-mode) ----------
  // Toggle the horizontal picker (on the value row, inside its `unitRow` divH) vs the
  // vertical picker (its own row) off the form's `ui-form-condensed` class — measured
  // live: present in the narrow label-on-top layout, absent in wide. The horizontal
  // picker reads correctly only wide; the vertical one is a normal stacked row that
  // reads correctly narrow.
  //  - HORIZONTAL pickers hide via a plain `style.display` write on the picker root:
  //    they live INSIDE the `unitRow` divH (NOT a direct section child), so a section
  //    collapse hides the wrapper for them and `appendCollapsibleSection`'s `updateWgts`
  //    loop never touches the picker itself — a direct write is collapse-safe here.
  //  - VERTICAL pickers ARE direct section children. `setControlVisible` records the mode
  //    intent as `data-src-hidden`, but its own direct `display` write would pop the row
  //    back under a COLLAPSED header when this fires from the resize-driven observer
  //    (decision-audit finding — the failure mode `setControlVisible`'s own doc comment
  //    warns about). So AFTER setting the flags we call each affected section's `reapply`
  //    handle, which recomputes `display` as (expanded AND mode-shown) — composing the
  //    mode axis with the collapse axis instead of overriding it. The asymmetry with the
  //    horizontal write is load-bearing: the horizontal picker's collapse is handled by
  //    its wrapper, the vertical picker's must be composed here.
  const unitPairs: ReadonlyArray<readonly [
    DG.InputBase<string | null>, DG.InputBase<string | null>,
  ]> = [
    [timeUnitInput, timeUnitInputV],
    [concUnitInput, concUnitInputV],
    [doseUnitInput, doseUnitInputV],
    [bodyWeightUnitInput, bodyWeightUnitInputV],
  ];
  const applyUnitPickerMode = (condensed: boolean): void => {
    for (const [h, v] of unitPairs) {
      h.root.style.display = condensed ? 'none' : ''; // horizontal: wide only (collapse-safe)
      setControlVisible(v.root, condensed); // vertical: record mode intent (data-src-hidden)
    }
    // Compose the mode intent with each section's collapse state (the vertical rows live
    // in Profile data + Context), so a mode flip while a section is collapsed leaves it hidden.
    reapplyProfileData();
    reapplyContext();
  };
  // Initial state + react to every wide/narrow flip. The platform sets the class after
  // layout, so: read it now, observe the form for later toggles, and re-read once after
  // layout settles (covers the class being added between this line and the observe call).
  applyUnitPickerMode(form.classList.contains('ui-form-condensed'));
  const unitModeObserver = new MutationObserver(() => {
    applyUnitPickerMode(form.classList.contains('ui-form-condensed'));
  });
  unitModeObserver.observe(form, {attributes: true, attributeFilter: ['class']});
  setTimeout(() => applyUnitPickerMode(form.classList.contains('ui-form-condensed')), 100);
  // Tear the observer down when the view closes (mirrors the `closeSub` idiom below), so a
  // closed mapping pane does not leak a live observer holding the form + the 8 unit inputs.
  const unitModeCloseSub = grok.events.onViewRemoved.subscribe((v: DG.ViewBase) => {
    if (v === view) {unitModeObserver.disconnect(); unitModeCloseSub.unsubscribe();}
  });

  // Y-axis toggle. Re-binds yColumnNames on the overlay chart and on every
  // facet plot (via the pre-derived `log(<conc>)` column or the raw `<conc>`
  // column from `buildChartDataFrame` / `buildFacetDataFrames`). No-op
  // before Compute (no live viewers yet); the state still flips so the
  // first Compute renders in the chosen mode.
  const applyYAxisModeToCharts = (): void => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const tempBag = view.temp as Record<string, any>;
    // UC-04 sparse Composite scatter — its own fixed log column ('Conc' /
    // 'log(Conc)'), independent of the mapped conc name, so it is toggled BEFORE
    // the concName guard below (which gates only the conc-name-derived charts).
    const sparseChart = tempBag[SPARSE_CHART_TAG] as DG.Viewer | null;
    if (sparseChart != null)
      sparseChart.setOptions({y: useLogY ? 'log(Conc)' : 'Conc'});
    // VR-03 Mean ± SD scatterPlot — its own fixed pre-derived columns
    // (Mean/Lower/Upper ↔ logMean/logLower/logUpper), independent of the mapped
    // conc name, so toggled here (before the concName guard below, which gates only
    // the conc-name-derived charts). Swaps the mean y AND the ±SD whisker columns.
    const meanSdChart = tempBag[MEAN_SD_CHART_TAG] as DG.Viewer | null;
    const meanSdProfile = tempBag[MEAN_SD_PROFILE_TAG] as DG.DataFrame | null;
    if (meanSdChart != null && meanSdProfile != null) {
      // Rebuild the per-scale frame (log drops all-BLQ rows) and swap the bound
      // statistic: arithmetic Mean±SD (linear) ↔ geometric mean ×/÷ gSD (semi-log)
      // + native log axis. See ensureMeanSdChart / meanSdChartFrame.
      const chartDf = meanSdChartFrame(meanSdProfile, useLogY);
      const opts: Record<string, unknown> = {
        y: useLogY ? 'geoMean' : 'Mean',
        yWhiskerMinColumnName: useLogY ? 'geoLower' : 'Lower',
        yWhiskerMaxColumnName: useLogY ? 'geoUpper' : 'Upper',
        yAxisType: useLogY ? 'logarithmic' : 'linear',
        // Description text tracks the scale; showSizeSelector / showTitle /
        // descriptionPosition are set once at create and persist (setOptions merges).
        description: useLogY ? 'Group geometric mean ×/÷ gSD (semi-log)' : 'Group mean ± SD',
      };
      // Re-pin BOTH axes for the NEW scale: y to the whisker extent (caps stay in
      // frame), x to the finite Time extent — the log frame drops all-BLQ rows
      // (often the latest times), so the X extent can shrink vs the linear frame.
      const range = meanSdYRange(chartDf, useLogY);
      if (range !== null) {opts.yMin = range.min; opts.yMax = range.max;}
      const xRange = meanSdXRange(chartDf);
      if (xRange !== null) {opts.xMin = xRange.min; opts.xMax = xRange.max;}
      // TEMP: guard ONLY the frame swap — a different-row-count dataFrame leaves
      // DG's scatterPlot pointHitTest BitSet stale (a platform bug); setOptions is left
      // unguarded so real config bugs surface.
      try {meanSdChart.dataFrame = chartDf;} catch {/* stale hit-test BitSet */}
      meanSdChart.setOptions(opts);
      // The swapped-in per-scale frame starts with an all-true filter — re-apply the
      // live By-group selection filter (if any) so a Semi-log toggle does not silently
      // desync the chart from the grid selection. No-op when grouping is off.
      (tempBag[BYGROUP_APPLY_TAG] as (() => void) | null)?.();
    }
    // Use the conc column the charts were BUILT with (last Compute), NOT the
    // live mapping input — charts are not rebuilt on input change, so the live
    // value may differ and would not exist in the stale chart frame.
    const concName = tempBag[PROFILE_CONC_NAME_TAG] as string | null;
    if (concName == null) return;
    const yColName = useLogY ? `log(${concName})` : concName;
    const overlayChart = tempBag[PROFILE_CHART_TAG] as DG.Viewer | null;
    const scalePatch: Partial<DG.ILineChartSettings> = {
      yColumnNames: [yColName],
      description: useLogY ? 'Log concentration–time profile' : 'Concentration–time profile',
    };
    if (overlayChart != null) overlayChart.setOptions(scalePatch);
    // Keep the recorded Combined options current, and the Report panel's copy of that chart in step
    // (it is the SAME chart — it must follow the scale switch too).
    const opts = tempBag[PROFILE_CHART_OPTS_TAG] as Partial<DG.ILineChartSettings> | null;
    if (opts != null) tempBag[PROFILE_CHART_OPTS_TAG] = {...opts, ...scalePatch};
    (tempBag[REPORT_PANEL_TAG] as {chart: DG.Viewer} | null)?.chart.setOptions(scalePatch);
    const facetPlots = tempBag[FACET_PLOTS_TAG] as DG.Viewer[] | null;
    if (facetPlots != null) {
      // Spike `facet-lambdaz-line-points`: the λz line is inherently SEMI-LOG (a
      // log-linear fit is a straight line only in log space). Its content lives on
      // the frame and never changes with the toggle — we only GATE its visibility
      // (`showDataframeFormulaLines: useLogY`), so on the linear axis it is hidden
      // rather than mis-drawn in log-values near zero. No per-tile clear/add here
      // (that neither refreshed reliably nor was cheap — it flooded the console).
      for (const p of facetPlots)
        p.setOptions({yColumnName: yColName, showDataframeFormulaLines: useLogY});
    }
  };

  // Y-axis mode toggle — standard DG `ui.input.toggle` with its built-in label.
  // Default `useLogY = true` (declared earlier) → switch reads ON → chart renders
  // in semi-log mode on first Compute. Lives in the Settings dialog (opened by the
  // ribbon gear); the caption "Semi-log Y scale" is self-describing, so the dialog
  // carries no category headings. The dialog re-appends this input's `.root` on every
  // open. The input OBJECT persists in this closure, so its value and `onChanged`
  // reactivity survive dialog close/reopen — no separate ribbon alias is kept.
  //
  // **State synchronisation.** Single source of truth = the input's `onChanged`
  // event; the `if (newValue !== useLogY)` short-circuit prevents redundant
  // `applyYAxisModeToCharts()` invocations on echo / no-op events.
  //
  // Tooltip is rich-content, bound via {@link bindInputTip} (caption + editor)
  // for the same reason as every other input on this panel — see its docstring.
  const semiLogToggleInput = ui.input.toggle('Semi-log Y scale', {value: useLogY});

  // Results-grid column profile. The CONTROL reads `Compact`, default ON
  // (user-directed 2026-08-13); the STATE it writes is still
  // `RESULTS_EXTENDED_TAG`, i.e. `extended === !compact`.
  //
  // The label was inverted rather than the stored flag: `extended` is the
  // parameter name of `resultsGridVisibleOrder` / `applyResultsColumnProfile` and
  // the key `swapOrDockResults` replays on every re-Compute, so flipping the
  // storage would have rippled through both modules and their tests to rename one
  // caption. The inversion therefore lives at this single UI edge — the ONLY
  // place the two conventions meet, and it is spelled out on both lines below so
  // no reader has to infer it.
  //
  // Default is unchanged in BEHAVIOUR: the grid already opened compact (the old
  // `Full` toggle was OFF by default). Only the caption's polarity changed, so no
  // stored recipe or `view.temp` value needs migrating.
  //
  // Lives in the Settings dialog (opened by the ribbon gear), captioned "Compact NCA
  // results" — self-describing, so the dialog carries no category headings
  // (user-directed 2026-09-20, superseding compact-toggle-banner-strip B1, which had it
  // on the results-grid banner strip). The input OBJECT persists in this closure — the
  // dialog re-appends its `.root` on every open — so its value and reactivity survive
  // close/reopen; no separate banner alias is kept.
  //
  // Reads the grid from `view.temp` AT CHANGE TIME rather than capturing it.
  // `swapOrDockResults` keeps one grid instance alive across Computes and only
  // rebinds its DataFrame, so the lookup is stable — and unlike a captured
  // reference it cannot go stale if that ever changes. A null grid (toggled
  // before the first Compute) is a no-op: the tag still records the choice, and
  // `swapOrDockResults` applies it when the grid appears.
  //
  // Plain `tooltipText`, not a `richTip`: two short clauses (one per toggle position)
  // have nothing for the title/summary/detail structure to carry, and the rich panel's
  // own chrome would dwarf them. Both positions appear ("ON —" / "OFF —"), matching
  // Semi-log's tooltip and the rule that a toggle tip answers "what happens if I flip
  // this", not "what is it doing now". OFF names what Full ADDS (study metadata, route)
  // so the scientist knows what compacting hides; ON keeps "subject", which says
  // identity survives it.
  const compactToggleInput = ui.input.toggle('Compact NCA results', {
    value: !resultsExtended(view),
    tooltipText: 'ON — subject + NCA parameters only. OFF — all columns (study metadata, route).',
  });
  compactToggleInput.onChanged.subscribe(() => {
    const compact = compactToggleInput.value ?? true;
    const extended = !compact;
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const bag = view.temp as Record<string, any>;
    bag[RESULTS_EXTENDED_TAG] = extended;
    const grid = bag[RESULTS_GRID_TAG] as DG.Grid | null;
    if (grid != null) applyResultsColumnProfile(grid, extended);
  });

  // `Violated` pills column (acceptance-reason-pills, user-directed 2026-10-19). Default = the
  // app-start PowerGrid check: a PowerGrid that cannot draw the pills readably starts it OFF, and
  // the tooltip says why. Plain `tooltipText` with both positions, like `Compact NCA results`.
  // Applies live to the results AND Sparse NCA grids; the flag rides `view.temp`, so a re-Compute
  // keeps the choice.
  const pgSupport = detectPowerGridSupport();
  const violatedToggleInput = ui.input.toggle('Violated column', {
    value: reasonPillsShown(view),
    tooltipText: 'ON — a Violated column beside Acceptance shows each violated λz criterion as a ' +
      'pill. OFF — hidden; the reasons stay in the Acceptance tooltip.' + (pgSupport.pills ? '' :
      ` Off by default: ${pgSupport.version === null ? 'PowerGrid is not installed' :
        `PowerGrid ${pgSupport.version} draws the pills unreadably`} — needs PowerGrid ` +
        `${PILLS_MIN_POWERGRID} or later, or a ${PILLS_MIN_POWERGRID_DEBUG}.X debug build.`),
  });
  violatedToggleInput.onChanged.subscribe(() => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    (view.temp as Record<string, any>)[REASON_PILLS_VIEW_TAG] = violatedToggleInput.value ?? false;
    applyReasonPillsVisibility(view);
  });

  bindInputTip(semiLogToggleInput, richTip({
    title: 'Y-axis scale',
    summary: 'ON — log10 concentration axis. PK convention: terminal-phase slope reads ' +
      'as a straight line.',
    detail: 'OFF — linear concentration axis.',
  }));

  semiLogToggleInput.onChanged.subscribe((): void => {
    const newValue = semiLogToggleInput.value ?? false;
    if (newValue !== useLogY) {
      useLogY = newValue;
      applyYAxisModeToCharts();
    }
  });

  // "Group by" — per-group aggregation control (#11), now a MULTI-column picker
  // (F1, multi-column-group-by): 0 columns = ungrouped; N columns = composite
  // grouping keyed by the tuple. Hidden until the first Compute that yields grouping
  // candidates (revealed on the ribbon like the Y-axis toggle). Candidates are
  // rebound each Compute (`updateGroupByColumns`) and scoped to the observations
  // frame so both consumers (aggregate grid + mean±SD chart) honour every key.
  // Stock `ui.input.columns` (adopt the platform control).
  // Bound to `df` initially as a placeholder; the first Compute rebinds it to the
  // results frame with the candidate filter.
  //
  // `showSelectedColsOnTop: false` is LOAD-BEARING, not cosmetic (BUG-11). The option
  // defaults to TRUE, and the platform then renders the selected columns in a block on top
  // WITHOUT removing them from the main list — so every selected column is drawn twice, the
  // dialog's counter reads double (3 selected → "6 checked"), and OK commits ALL ticked rows,
  // i.e. the doubled set. That is the GENERATOR of the unbounded growth this bug was reported
  // for: 3 → 6 → 12 …, laundered through `getCurrentRecipe()` into the persisted recipe.
  // Measured directly against the platform widget (throwaway probe, 2026-07-17): default →
  // 6 checked for 3 selected; `false` → 3. The `dedupeGroupBySelection` guard makes the
  // selection CONVERGE even with the default on, but this removes the generator rather than
  // absorbing it. Cost: the picker's selected rows no longer float to the top of the list —
  // a fair trade for a list that shows each column exactly once.
  //
  // NO TEST GUARDS THIS LINE (deliberate — the DOM probe that measured it was deleted rather
  // than shipped: it depended on DG's private CSS + a timing race, and would have gone red on
  // a future DG fix). If you drop the option, the dedup guard keeps the DATA clean and only
  // the rendered list doubles — a visual regression, caught by eye, not by CI. When DG stops
  // double-listing, this option can simply go. See TODO.md BUG-11.
  const groupByInput = ui.input.columns('Group by',
    {table: df, value: [], showSelectedColsOnTop: false});
  bindInputTip(groupByInput, richTip({
    title: 'Group by',
    summary: 'Summarise NCA parameters per group (mean / SD / %CV / N). Pick one or more.',
    detail: 'Pick one or more metadata columns carried from the input table (e.g. compound, ' +
      'dose, sex); several columns group by the composite tuple (e.g. Sex × Dose). Leave empty ' +
      'for per-subject results only.',
    note: 'A group that pools true CL (IV) and apparent CL/F (EV) withholds CL, Vz, Vss, ' +
      'MRT and AUMC (shown as NE) — add CL basis or Route to separate them.',
  }));
  // Group-by drives BOTH the aggregated parameter grid and the VR-03 Mean ± SD
  // chart's split — one control, two consumers (re-split together, no recompute).
  // The suppress guard skips this while we set the value programmatically on rebind.
  groupByInput.onChanged.subscribe(() => {
    if (suppressGroupByChange) return;
    const raw = groupByInput.value ?? [];
    const frame = liveResultsDf();
    // Pre-Compute there is no results frame to resolve against and no consumer to feed —
    // the first Compute rebinds the picker and restores from `pendingGroupByCols` anyway.
    if (frame === null) return;
    // BUG-11: the picker's OK commits `value ++ tableColumns` ticked BY NAME, so a plain
    // OK re-adds every already-selected name (measured 5 → 8 on a no-op OK). Collapse to
    // the distinct names as LIVE columns of this frame, and — only when that actually
    // changed something — write the collapsed value back so the widget converges instead of
    // growing. Without the write-back the ribbon caption and the "Select columns..." list
    // would keep showing (and re-committing) duplicates even though every consumer
    // downstream reads the deduped names. Guarded, and the guard is what stops it looping:
    // our own assignment re-enters this handler and would otherwise recurse.
    const cols = dedupeGroupBySelection(raw.map((c: DG.Column) => c.name), frame);
    groupBySelectedNames = cols.map((c) => c.name);
    if (cols.length !== raw.length) {
      suppressGroupByChange = true;
      try {
        groupByInput.value = cols;
      } finally {
        suppressGroupByChange = false;
      }
    }
    applyGrouping();
    refreshMeanSdChart(false);
    // F3/D3: Group-by is a persisted recipe field (`groupByCols`), so a grouping change is an
    // unsaved edit — mark dirty so "Update analysis" can persist it. `groupByInput` is wired
    // separately from `onAnyChange`, so this is its own hook. No-op unless study-backed / replaying.
    markAnalysisDirty(view);
    // Land the scientist on the Mean ± SD chart — the visualization of the split
    // they just picked — but ONLY when a group is actually selected. Clearing the
    // selection (cols empty) removes the grouping and pools to one curve, which is
    // not "a choice to look at", so we leave the active tab where it is
    // (user-directed 2026-08-28).
    if (cols.length > 0) activateMeanSdTab(view);
  });

  // Interactive guide trigger on the ribbon. Resolves the dynamic refs
  // (results grid + profile chart) at click time so the tour reflects the
  // current post-Compute state without a separate refresh path. Pure
  // read of `view.temp`; never mutates view state.
  //
  // View-close cleanup: the tour's dim overlay (`#nca-tour-overlay`) is
  // appended to `document.body`, OUTSIDE the analysis view, so if the user
  // closes the view mid-tour the overlay would orphan and persist on top
  // of the Datagrok shell. Subscribing `onViewRemoved` and synthesising a
  // Done click guarantees the overlay is removed alongside the view.
  const infoBtn = ui.iconFA('info-circle', () => {
    // eslint-disable-next-line @typescript-eslint/no-explicit-any
    const tempBag = view.temp as Record<string, any>;
    // The tour gates every post-Compute page on `isOnScreen`, which starts with
    // `el.isConnected` — and since the 2026-08-13 layout the results grid is a TAB
    // of the source data, so its root is detached whenever another table tab is in
    // front. Its docstring's premise ("dock manager can close any subset") held
    // while results had their own split; now a perfectly-docked results grid reads
    // as closed and the tour silently drops its first and most important page.
    // Front the tab before the anchors are measured — the same courtesy the
    // pre-Compute pages already get from the describer's `ensureVisible` hook.
    activateResultsTab(view);
    const doneBtn = startNcaTour({
      mainGrid: view.grid.root,
      mappingPane: paneRoot,
      subjectInput, timeInput, concInput, computeBtn: confirmBtn, settingsBtn,
      addToStudyBtn: saveBtn,
      // Group-by control: only surfaced once revealed (a Compute produced
      // grouping candidates); null otherwise so the tour skips its page.
      groupByControl: groupByRevealed ? groupByInput.root : null,
      resultsGrid: (tempBag[RESULTS_GRID_TAG] ?? null) as DG.Grid | null,
      facetGrid: (tempBag[FACET_GRID_DIV_TAG] ?? null) as HTMLDivElement | null,
      profileChart: (tempBag[PROFILE_CHART_TAG] ?? null) as DG.Viewer | null,
    });
    // Null when startNcaTour surfaced a warning instead of starting a
    // tour (zero on-screen anchors) — nothing to clean up on view-close.
    if (doneBtn == null) return;
    const closeSub = grok.events.onViewRemoved.subscribe((v: DG.ViewBase) => {
      if (v === view) {
        doneBtn.click();
        closeSub.unsubscribe();
      }
    });
  }, 'Run interactive guide');
  infoBtn.classList.add('nca-ribbon-guide-icon');

  // Settings gear — opens a modal with the two view-presentation toggles. Each caption
  // is self-describing ("Semi-log Y scale", "Compact NCA results"), so the dialog carries
  // NO category headings (user-directed 2026-09-20). The toggle INPUT objects live in this
  // closure (single source of truth, handlers already wired above); each open re-appends
  // their roots into a fresh form, so their value + reactivity persist across close/reopen.
  //
  // Both changes apply LIVE via each toggle's `onChanged`, so OK merely confirms and
  // closes; Cancel / × close the same way (the choice is already applied either way).
  // `onOK` is what puts the OK button on the dialog.
  //
  // Single-instance gate mirrors the Add-to-study icon (nca-add-to-study-dialog-dedup):
  // a closure boolean set SYNCHRONOUSLY on click and lifted on `onClose`, so rapid
  // clicks can never stack a second dialog. No `await` here (nothing to load), so the
  // add-to-study `.catch` re-enable path is unnecessary — the flag is set immediately
  // before `.show()` and cleared on any dismissal (OK / Cancel / ×).
  let settingsDialogOpen = false;
  const settingsBtn = ui.iconFA('cog', () => {
    if (settingsDialogOpen) return; // a dialog is already open — do not stack a second
    settingsDialogOpen = true;
    // Fresh form each open; the persistent toggle roots are re-appended (no headings,
    // no custom form styling — DG sizes the dialog to its content).
    const form = ui.form([]);
    form.append(semiLogToggleInput.root, compactToggleInput.root, violatedToggleInput.root);
    // onOK adds the OK button (changes already applied live, so it just closes).
    // Subscribe onClose BEFORE show() so no dismissal is missed (nca-add-to-study.ts).
    const dlg = ui.dialog('Settings').add(form).onOK(() => {});
    dlg.onClose.subscribe(() => {settingsDialogOpen = false;});
    dlg.show();
    // Fit the dialog to its toggles (user-directed 2026-09-30): the default dialog box left a
    // tall empty band under them.
    dlg.root.style.height = 'auto';
    dlg.root.style.minHeight = '0';
  }, 'Settings — chart scale and results columns');

  // "Add to study" (AC-C4 populate flow). Active only after a successful Compute —
  // the entry's summary cache is built from the last `perProfile` (no recompute,
  // thesis #1) and the table is written into the study folder as the analysis data
  // file. Inactive state is just LIGHTER (still clickable, no action) — not
  // pointer-events:none — per user-directed 2026-06-29.
  // Single-instance gate: repeated clicks on the icon used to STACK a new "Add analysis to study"
  // dialog each time (openAddToStudyDialog unconditionally `.show()`s a fresh one). Because the
  // dialog opens after an `await` (study list + user context), the gate is set SYNCHRONOUSLY here —
  // before the await — so clicks during that window are also suppressed; it is lifted when the
  // dialog closes (onClosed, any dismissal) or if it fails to open (the `.catch` below), so it can
  // never latch permanently (nca-add-to-study-dialog-dedup).
  let addToStudyDialogOpen = false;
  // Add to study — a text-only button (user-directed 2026-09-29): the completing action of the
  // Report stage. It STAYS ACTIVE after the analysis is in a study (D5 "one copy" reversed,
  // user-directed 2026-10-05): each Add creates a NEW analysis from the current confirmed results,
  // so one source table can back several study analyses. No duplicate guard — the scientist decides
  // when to click. A cancelled dialog leaves it active.
  const saveBtn = ui.button('Add to study', () => {
    // Inactive until a successful Compute: both are set in lockstep at the fire
    // site, so the null guard covers the grouping substrate too.
    if (lastPerProfile === null || lastResults === null) return;
    if (addToStudyDialogOpen) return; // a dialog is already open — do not stack a second
    addToStudyDialogOpen = true;
    refreshSaveEnabled(); // dim + not-allowed while open
    const reenable = (): void => {addToStudyDialogOpen = false; refreshSaveEnabled();};
    const perProfile = lastPerProfile;
    const results = lastResults;
    void (async () => {
      // D5: a view bound to a STANDALONE analysis MOVES it. Unsaved edits are saved first, so the
      // move carries exactly what is on screen (one move implementation reads the stored entry).
      // D14 soft gate (F-S6.5): the dialog itself opens on the not-reviewed warning when profiles are
      // still to review (Cancel / Add anyway), then reveals its inputs.
      const pending = liveWorkflowState(view).stage.pending;
      const binding = paneTemp[STUDY_BINDING_TAG] as StudyBinding | null;
      const moveFrom = binding?.container === 'standalone' ? binding : undefined;
      if (moveFrom !== undefined && paneTemp[ANALYSIS_DIRTY_TAG] === true) await runUpdate();
      await openAddToStudyDialog({
        df, recipe: getCurrentRecipe(),
        // The summary cache of a study entry aggregates the REPORTABLE frame (AD-6).
        perProfile, resultsDf: reportableFrameOf(view, results), observations: lastObservations,
        presetStudyId: opts.targetStudyId,
        workflow: workflowPersistOf(view), moveFrom, pendingReview: pending,
        onSaved: (saved) => promoteToStudyBacked(saved), // F4/D7: SAVE targets the study entry now
        onClosed: reenable, // dialog dismissed (OK/Cancel/×) or never opened → lift the gate
      });
    })().catch((e) => {
      // A throw BEFORE `.show()` (e.g. user-context load) rejects here without an onClose — lift
      // the gate so the icon is not stuck disabled. The normal path resolves without error, so
      // this never double-fires with `onClosed`.
      reenable();
      grok.shell.error(`Could not open dialog: ${e instanceof Error ? e.message : String(e)}`);
    });
  });
  // Standard DG tooltip (lazy → reflects the current active/inactive state); the
  // cursor is the default DG icon cursor (same as the Compute play icon).
  ui.tooltip.bind(saveBtn, () => {
    if (lastPerProfile === null) return 'Compute the analysis first, then add it to a study';
    // Already in a study: the button stays active and a click adds the CURRENT results as a new
    // analysis (one table can back several study analyses — user-directed 2026-10-05).
    return studyBindingOf(view) !== null ?
      'Add the current results to a study (creates a new analysis)' :
      'Add this analysis to a study';
  });
  const refreshSaveEnabled = (): void => {
    // Inactive (lighter) when: (1) pre-Compute (still hoverable — the click is a no-op,
    // user-directed 2026-06-29); (2) a dialog is already open (also not-allowed cursor). Being in a
    // study NO LONGER dims it (D5 reversed, user-directed 2026-10-05): a re-add creates a new
    // analysis, so the button must stay active.
    saveBtn.style.opacity = (lastPerProfile !== null && !addToStudyDialogOpen) ? '' : '0.4';
    saveBtn.style.cursor = addToStudyDialogOpen ? 'not-allowed' : '';
  };
  // The stage refresh (which runs on every binding change) re-evaluates the button.
  addBtnRefresh = refreshSaveEnabled;
  refreshSaveEnabled(); // initial: lighter (inactive) until first Compute

  // ── "Update analysis" in-place re-save (analysis-update-in-place, F2–F4) ─────────────
  // In-place re-save of a STUDY-BACKED analysis: re-writes the source CSV + the study.json
  // entry so post-save edits (λz point overrides + any input/unit/rule/threshold change)
  // persist — distinct from "Add to study", which creates a NEW analysis. Enabled only while
  // study-backed AND dirty; the lit (green `ui-btn-ok`) button IS the unsaved-changes indicator
  // (D6). Copies saveBtn's dim + lazy-tooltip pattern (Rule 0). The control is now the stock DG
  // save button (cloud-upload + "Save"), replacing the earlier floppy-disk `ui.iconFA('save')`
  // icon (user-directed 2026-09-14). English-only (D5).
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const paneTemp = view.temp as Record<string, any>;
  // Native DG button (was a floppy `ui.iconFA('save')` icon) — the standard raised/green save
  // button (`ui-btn-raised ui-btn-ok`) with a cloud-upload glyph + "Save" label, matching the
  // platform's own save button (the stock DG save button is the reference). CSS uppercases the label to
  // "SAVE". `name` mirrors the platform button for parity.
  const updateBtn = ui.button([ui.iconFA('cloud-upload'), 'Save'], () => void runUpdate());
  updateBtn.classList.add('ui-btn-raised', 'ui-btn-ok');
  updateBtn.name = 'button-Save';
  /** SD-13: a View-only study — the server would refuse every write. */
  const viewOnly = (): boolean => (paneTemp[STUDY_BINDING_TAG] as StudyBinding | null)?.readOnly === true;
  ui.tooltip.bind(updateBtn, () => {
    if (viewOnly()) return VIEW_ONLY_TIP;
    if (paneTemp[ANALYSIS_DIRTY_TAG] !== true) return 'No unsaved changes';
    return paneTemp[STUDY_BINDING_TAG] == null ?
      `Save this analysis to ${UNGROUPED_SECTION}` :
      'Save changes to this analysis (data, settings, λz edits and review decisions)';
  });
  const refreshUpdateEnabled = (): void => {
    // SD-13: a View-only study shows SAVE disabled the platform way (grey `d4-disabled`, the tooltip
    // on DG's hover overlay) — the "nothing to save" dim below reads as an active button. Applied
    // once: the binding's access does not change for the life of the view.
    if (viewOnly()) {
      if (!updateBtn.classList.contains('d4-disabled')) ui.setDisabled(updateBtn, true, VIEW_ONLY_TIP);
      return;
    }
    // D2: SAVE is lit whenever there is something to save — bound or not (a never-saved analysis
    // is dirty from its first edit; its first save asks for a name, then goes to Ungrouped).
    const enabled = paneTemp[ANALYSIS_DIRTY_TAG] === true;
    // Visual dimming — NOT native `disabled` — so the state-dependent tooltip still explains WHY
    // there is nothing to save (a `disabled` button suppresses hover in some browsers). The button
    // stays clickable; the `runUpdate` guard makes an inactive click a no-op (user-directed
    // 2026-09-14). Inactive drops the green `ui-btn-ok` to the neutral raised button + dim.
    updateBtn.classList.toggle('ui-btn-ok', enabled);
    updateBtn.style.opacity = enabled ? '' : '0.5';
    updateBtn.style.cursor = enabled ? 'pointer' : 'not-allowed';
  };
  const runUpdate = async (): Promise<void> => {
    // Guard: dimmed when clean or View-only — a stray click is a no-op.
    if (paneTemp[ANALYSIS_DIRTY_TAG] !== true || viewOnly()) return;
    const binding = paneTemp[STUDY_BINDING_TAG] as StudyBinding | null;
    // D1: Save works at ANY stage — with or without results. The cache is built from the
    // REPORTABLE frame of the last confirmed compute (AD-6), and the workflow state rides along.
    const input = {
      df, recipe: getCurrentRecipe(),
      computed: lastPerProfile !== null && lastResults !== null ?
        {perProfile: lastPerProfile, resultsDf: reportableFrameOf(view, lastResults),
          observations: lastObservations} : null,
      workflow: workflowPersistOf(view),
    };
    if (binding == null) {
      // The FIRST save asks for the analysis name, pre-filled with the table name (user-directed
      // 2026-09-30, superseding D2's silent first save). A name already taken is made unique by
      // saveNewStandaloneAnalysis (the toast names the final label). Cancel saves nothing.
      const nameInput = ui.input.string('Name', {value: view.dataFrame?.name || view.name});
      const dlg = ui.dialog('Save analysis').add(ui.form([nameInput])).onOK(() => {
        void (async () => {
          const bound = await saveNewStandaloneAnalysis(input, (nameInput.value ?? '').trim() || view.name);
          if (bound === null) return;
          bindSavedAnalysis(bound, bound.label);
        })();
      });
      dlg.show();
      nameInput.input.focus();
      return;
    }
    const res = await updateAnalysisInPlace({...input, binding});
    if (res.rebound !== undefined) {
      // F11: the source was replaced while the view was open and the work was saved beside as a
      // NEW analysis — bind to it, so later saves and the summary writeback target it.
      bindSavedAnalysis(res.rebound, res.rebound.containerLabel ?? binding.containerLabel ?? '');
    } else if (res.ok) {
      paneTemp[ANALYSIS_DIRTY_TAG] = false;
      refreshUpdateEnabled();
    }
  };
  // Stash so module-level markAnalysisDirty (called from the λz override sites, outside this
  // closure) can repaint the widget after flipping the dirty flag.
  paneTemp[UPDATE_REFRESH_TAG] = refreshUpdateEnabled;
  refreshUpdateEnabled(); // initial: disabled (not yet study-backed, or clean)
  // The reopen path binds the view AFTER this mount, then refreshes the stage: repaint SAVE there
  // too, so a View-only binding greys it out before any edit.
  addBtnRefresh = (): void => {
    refreshSaveEnabled();
    refreshUpdateEnabled();
  };
  // Promote the view to study-backed right after Add-to-study (F4/D7) so "Update" targets THIS
  // analysis without a reopen; the reopen path sets the same binding after mounting.
  // Bind the view to a just-saved analysis (first silent save, Add-to-study, or a move) so SAVE
  // targets it without a reopen (F4/D7). `containerLabel` is the study label (study) or the
  // analysis label (standalone) for the deep-link URL.
  const bindSavedAnalysis = (saved: StudyBinding, containerLabel: string): void => {
    paneTemp[STUDY_BINDING_TAG] = {...saved, containerLabel} satisfies StudyBinding;
    paneTemp[ANALYSIS_DIRTY_TAG] = false;
    // F4/D7: from now on, every recompute of a CLEAN view refreshes the container's summary cache
    // (like a reopened analysis) via the SAME factory the reopen path uses (AD-8 skips it while
    // dirty). Writes only the cache + workflow state, never the source CSV (SAVE's job).
    studyOnComputed = makeStudyBackedOnComputed(saved.analysisId, saved.label, saved.passthroughCols, df,
      saved.source, saved.sourceHash, () => markAnalysisDirty(view));
    if (saved.container === 'study') updateNcaStudyPath(view, containerLabel, saved.label);
    else updateNcaAnalysisPath(view, saved.label);
    refreshUpdateEnabled();
    // The hint now names the study; the Add-to-study button STAYS active (a re-add creates another
    // analysis — user-directed 2026-10-05), so refreshStage only re-evaluates the Compute-gating.
    refreshStage();
  };
  const promoteToStudyBacked = (saved: SavedAnalysis): void => {
    bindSavedAnalysis({
      container: 'study', studyId: saved.studyId, analysisId: saved.analysisId, label: saved.label,
      source: saved.source, sourceHash: saved.sourceHash, passthroughCols: saved.passthroughCols,
    }, saved.studyLabel);
  };

  // Unify the gaps around ribbon-group dividers (user-directed 2026-06-30). Bare
  // icons sit tighter to a divider than inputs (whose caption carries padding), so
  // add an explicit edge margin to the icons + the trailing choice input. One knob
  // (`RIBBON_EDGE`) drives every icon↔divider gap — to match the toggler-label→divider reference.
  const RIBBON_EDGE = '8px';
  // Default flat ui.button — visually secondary to the primary raised-green SAVE (updateBtn).
  saveBtn.style.marginLeft = RIBBON_EDGE; // divider → Add to study
  // SAVE → guide (in-group), guide → the divider before the stage's controls.
  infoBtn.style.marginLeft = RIBBON_EDGE;
  infoBtn.style.marginRight = RIBBON_EDGE;
  // Settings gear: icon-edge margins on both sides so the middle group reads balanced
  // whether it holds just the gear or [Settings, Group by] (Report).
  settingsBtn.style.marginLeft = RIBBON_EDGE; // divider → Settings
  settingsBtn.style.marginRight = RIBBON_EDGE; // Settings → Group by (in-group), or → divider
  // Group-by → Settings is now an in-group gap (both live in the middle group).
  groupByInput.root.style.marginRight = RIBBON_EDGE;

  // Ribbon groups (stable positions, user-directed 2026-09-29): [SAVE, ⓘ] | [Settings, Group by]
  // (stage-scoped: off Setup; Group by on Report only) | [Add to study] (Report only). Compute
  // is gone from the ribbon — the pane's Confirm button computes (D17). See applyStageRibbon.
  applyStageRibbon();

  // AD-5 replay LAST — after autodetect + form assembly + all wiring — so the
  // saved recipe overrides every autodetected default. The hydration path then
  // calls `triggerCompute()` to compute on the replayed mapping (AC-D1).
  if (opts.recipe !== undefined)
    applyRecipe(opts.recipe);

  (view.temp as Record<string, unknown>)[RECIPE_FN_TAG] = getCurrentRecipe;
  // First paint of the stage strip + Revert link (QA 2026-09-29: without it the strip stayed EMPTY
  // — zero height — until the first Confirm/edit, and the Revert link showed un-hidden).
  refreshStage();
  return {getCurrentRecipe, triggerCompute: runCompute};
}

/** Dock the results grid as a TAB of the source-data grid, or swap the DataFrame
 *  in place if a previous results grid is already docked.
 *
 *  `FILL` onto the platform's own `view.grid` node (user-directed 2026-08-13,
 *  superseding the `'down', 0.5` split): both are tables the scientist reads, so
 *  they share one group — and whichever is in front gets the full height of the
 *  tables row instead of half of it. The charts row above is what the freed split
 *  now pays for (see {@link CHARTS_ZONE_TAG}).
 *
 *  This FILL is also the FIRST dock operation of the post-Compute layout, and
 *  that is what makes it safe: the historical footgun — a `FILL` next to
 *  `view.grid` raising inside the dock manager — only ever reproduced AFTER a
 *  split had already re-parented that node.
 *
 *  Re-Compute branch (`existingGrid.dataFrame = results`): Datagrok keeps the
 *  grid instance alive but rebinds its data. Grid-side settings such as
 *  `gridCol.isTextColorCoded`, JSON-column visibility, and the registered
 *  `onCellTooltip` handler may be reset when the underlying dataframe changes;
 *  re-applying `configureResultsGrid` is idempotent and cheap, so we do it
 *  unconditionally on every swap to keep presentation behaviour stable. */
function swapOrDockResults(
  view: DG.TableView,
  results: DG.DataFrame,
  thresholds: ResultsColorThresholds,
): DG.DockNode | null {
  configureResultsDataFrame(results, thresholds);

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  // `!= null` not `!== undefined`: view.temp returns null for a never-set key
  // (BUG-06 root cause); harmless here since results docks first, but kept
  // consistent so the latent `null.dataFrame` crash can't surface.
  const existingGrid = tempBag[RESULTS_GRID_TAG] as DG.Grid | null;

  // The scientist's compact/Full choice lives on the VIEW, not the frame: this
  // function replaces the DataFrame on every Compute, so a `df.temp` tag (the
  // study grid's mechanism, where the frame persists) would silently reset their
  // choice each run. `!= null` because view.temp returns null for an unset key.
  const extended = resultsExtended(view);
  // The `Violated` switch rides the frame, which is new each Compute — re-stamp it before the
  // profile below reads it.
  results.setTag(REASON_PILLS_TAG, reasonPillsShown(view) ? 'true' : 'false');

  if (existingGrid != null) {
    existingGrid.dataFrame = results;
    configureResultsGrid(existingGrid, thresholds, extended);
    // Re-Compute reuses the existing results node — return it so callers can
    // anchor sibling tabs (F% / TK panels) to the same node.
    return (tempBag[RESULTS_NODE_TAG] as DG.DockNode | undefined) ?? null;
  }

  const grid = DGNS.Viewer.grid(results);
  configureResultsGrid(grid, thresholds, extended);
  // A Review cell opens the review of its profile (user-directed 2026-09-30): "To review" is styled
  // as a link, and a decided cell is how the analyst revisits a decision once nothing is left
  // (the Complete hint names it). Subscribed ONCE per grid — the grid survives a re-Compute (only
  // its frame is swapped), so the frame is read at click time. Deferred a tick so the grid's own
  // click handling (current row / context panel) settles before the Review pane is shown.
  grid.onCellClick.subscribe((gc: DG.GridCell) => {
    if (!gc.isTableCell || gc.tableColumn?.name !== REVIEW_COL) return;
    const v = gc.cell.value;
    if (v === null || v === '') return;
    const subj = grid.dataFrame.col('SUBJECT');
    const row = gc.tableRowIndex;
    if (subj === null || row === null || row < 0) return;
    const subject = subj.getString(row);
    setTimeout(() => startReview(view, subject), 0);
  });
  const mainGridNode = view.dockManager.findNode(view.grid.root) ?? null;
  // Bare `attachBannerHost(grid.root)` — the same wrapper every other docked panel
  // uses, carrying this panel's compute-notice banners. The Compact toggle no longer
  // rides here as a leading control: it moved into the ribbon Settings dialog together
  // with Semi-log (user-directed 2026-09-20, superseding compact-toggle-banner-strip B1).
  const node = view.dockManager.dock(
    attachBannerHost(grid.root), DGNS.DOCK_TYPE.FILL, mainGridNode, 'NCA results');
  tempBag[RESULTS_GRID_TAG] = grid;
  tempBag[RESULTS_NODE_TAG] = node;
  return node;
}

// The results panel used to be wrapped in an extra `divV([controlStrip,
// bannerHost])`. That wrapper is GONE (user-directed 2026-08-13), and its removal
// is the fix for a real layout defect, not a tidy-up:
//
//   wrap  = divV([strip, body])   height:100%  minHeight:0
//   ├ strip                       flex:0 0 auto  padding:2px 8px   ≈26px
//   └ body = attachBannerHost(…)  height:100%   ← claims the WHOLE column
//
// `body` asked for 100% of `wrap` while `strip` had already taken ~26px of the
// same flex column, so the children summed to 100% + 26px and the panel
// overflowed by exactly the strip's height — sliding the grid's last row under
// the dock tab bar. (`minHeight: 0` was set on the CONTAINER, not on `body`, the
// flex item that actually had to shrink, so nothing reined it in.)
//
// Every sibling panel — By group, Bioavailability, Dose Proportionality,
// Dose-Normalized, Sparse NCA, Composite Profile, the TK grids — docks a bare
// `attachBannerHost(grid.root)`, and not one of them has ever shown the overlap.
// The results panel was the only one carrying an extra flex layer and the only
// one with the bug, so dropping the wrapper returns it to the reference pattern
// rather than working around the symptom.
//
// The `Compact` toggle it once hosted rode the ribbon beside Semi-log / Group by from
// 2026-08-13, then the results grid's banner strip from 2026-08-28, and as of
// 2026-09-20 lives in the ribbon Settings dialog together with Semi-log (see the gear
// in `mountMappingPane`). `swapOrDockResults` therefore docks a bare
// `attachBannerHost(grid.root)` — the reference pattern every sibling panel uses,
// carrying compute-notice banners only. The download button is REMOVED — see the note
// at `RESULTS_EXTENDED_TAG`'s toggle for what that costs.

/** Dock (or swap in place) the per-group aggregated grid as a FILL-tab on the
 *  results node — "By group" sits next to "NCA results" in the bottom area
 *  (no extra viewport budget consumed). Mirrors {@link swapOrDockResults}'s
 *  swap-on-recompute pattern via `view.temp` handles. */
/** Hidden helper column on the By-group aggregate: profiles in the group whose λz-dependent
 *  values an Exclude λz decision withheld (F-S5.2). Never shown; read by the cell hook. */
const EXCLUDED_N_COL = '_EXCLUDED_N';

function hideHelperColumn(grid: DG.Grid, name: string): void {
  const gc = grid.columns.byName(name);
  if (gc !== null) gc.visible = false;
}

/** Write {@link EXCLUDED_N_COL} onto `agg`: per aggregate row, the number of excluded profiles in
 *  that group, keyed by the SAME composite group key the aggregator used. */
function appendExcludedCounts(
  view: DG.TableView, raw: DG.DataFrame, agg: DG.DataFrame, keys: ReadonlyArray<string>,
): void {
  const controller = (view.temp as Record<string, unknown>)[NCA_CONTROLLER_TAG] as NcaController | null;
  const byGroup = reportableResults(raw, controller?.excludedKeys() ?? new Set()).excludedCountBy(keys);
  const keyCols = keys.map((k) => agg.col(k)).filter((c): c is DG.Column => c !== null);
  const out = new Int32Array(agg.rowCount);
  if (keyCols.length === keys.length)
    for (let i = 0; i < agg.rowCount; i++) out[i] = byGroup.get(compositeGroupKey(keyCols, i)) ?? 0;
  if (agg.col(EXCLUDED_N_COL) !== null) agg.columns.remove(EXCLUDED_N_COL);
  agg.columns.add(DGNS.Column.fromInt32Array(EXCLUDED_N_COL, out));
}

/** True when an aggregate cell is a λz-dependent MEAN / SD / %CV that is empty BECAUSE every
 *  contributor was excluded (N_<P> = 0 and the group has excluded profiles) — "not reportable",
 *  not "not estimated". */
function notReportableCell(agg: DG.DataFrame | null, colName: string, row: number): boolean {
  if (agg == null || row < 0) return false;
  const m = /^(MEAN|SD|CV_PCT)_(.+)$/.exec(colName);
  if (m === null || !isLzDependentParam(m[2])) return false;
  const cell = agg.col(colName);
  const n = agg.col(`N_${m[2]}`);
  const ex = agg.col(EXCLUDED_N_COL);
  if (cell === null || n === null || ex === null || !cell.isNone(row)) return false;
  return Number(n.get(row)) === 0 && Number(ex.get(row)) > 0;
}

function swapOrDockAggregated(view: DG.TableView, agg: DG.DataFrame): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  // `!= null` not `!== undefined`: view.temp returns null for a never-set key, so
  // the first dock would otherwise run `null.dataFrame = agg` (BUG-06 sibling).
  const existingGrid = tempBag[AGG_GRID_TAG] as DG.Grid | null;
  if (existingGrid != null) {
    existingGrid.dataFrame = agg;
    // Re-apply presentation: the swapped-in frame is freshly built, so its stat
    // columns carry no friendlyName / format tags until re-tagged.
    configureAggregatedGrid(existingGrid, FRAGILE_N);
    hideHelperColumn(existingGrid, EXCLUDED_N_COL);
    return;
  }
  const grid = DGNS.Viewer.grid(agg);
  // FRAGILE_N is passed so the header tooltips state the SAME threshold the
  // bold-font hook below applies — one constant, one sentence.
  configureAggregatedGrid(grid, FRAGILE_N);
  hideHelperColumn(grid, EXCLUDED_N_COL);
  // F4 (multi-column-group-by): bold-font the per-param N cells where N < FRAGILE_N.
  // Attached ONCE on grid creation (the handler reads the live cell value, so it
  // survives dataFrame swaps — re-attaching per swap would stack handlers). NOT
  // setConditional (a value-heatmap is forbidden for fragility flags); bold survives
  // grayscale and scales per-row. Mirrors the study-overview grid's n<4 bold precedent.
  grid.onCellPrepare((gc) => {
    if (!gc.isTableCell || gc.tableColumn === null) return;
    const nm = gc.tableColumn.name;
    // F-S5.2: a λz-dependent statistic withheld because EVERY contributor was excluded.
    if (notReportableCell(grid.dataFrame, nm, gc.cell.rowIndex)) {
      gc.customText = 'not reportable';
      gc.style.textColor = 0xFF9497A0; // muted grey-4, like the per-subject excluded cells
      return;
    }
    if (nm !== 'N' && !nm.startsWith('N_')) return; // per-param contributor counts + group N
    const v = gc.cell.value;
    if (typeof v !== 'number' || !(v < FRAGILE_N)) return;
    boldGridCell(gc); // ONE shared bold idiom (results-view.ts) — no re-typed twin (rule 5)
  });
  // Anchor to the results node so the aggregated grid is a sibling TAB of the
  // per-subject results grid; fall back to the main grid node if results aren't
  // docked yet (defensive — Compute always docks results first).
  const anchor = resultsAnchor(view);
  const node = view.dockManager.dock(attachBannerHost(grid.root), DGNS.DOCK_TYPE.FILL, anchor, 'By group');
  tempBag[AGG_GRID_TAG] = grid;
  tempBag[AGG_NODE_TAG] = node;
}

/** Compute and dock the stratum-level F% grid when the results frame supports
 *  it (mixed route + per-row dose + a Sex metadata column), else close it. The
 *  trigger is a successful Compute (mirrors the auto-docked profile chart); no
 *  separate ribbon control. `carried` is the list of metadata column names
 *  appended by `appendGroupingColumns` — the Sex key is matched case-insensitively
 *  within it. */
function applyBioavailability(
  view: DG.TableView,
  results: DG.DataFrame,
  sexName: string | null,
  doseColName: string | null,
  hasRouteColumn: boolean,
): void {
  if (!hasRouteColumn || doseColName === null || sexName === null) {
    // Distinguish the suppression reasons so a scientist with a mixed-route frame
    // that has route + dose but no Sex column knows WHY no F% appears (sex
    // stratification is mandatory — AD-6 — so the absence is correct, but silent
    // suppression is opaque). The no-route / no-dose cases are self-evidently F%-
    // inapplicable, so they stay silent; only the route+dose-present-but-no-Sex
    // case gets an explanatory toast.
    if (hasRouteColumn && doseColName !== null && sexName === null) {
      // F% tab is not docked here (no Sex → no F%), so this banner falls back to
      // the results grid (B3 — the "not applicable" signal must not be dropped).
      collectNotice(view, {
        key: 'F_PCT_NEEDS_SEX', anchor: 'bioavailability', severity: 'info',
        // Action path is the Subject Sex picker (A4): the frame has route + dose
        // but no DESIGNATED, constant-per-subject sex column. Point the scientist
        // at the picker (or at adding a sex axis if the data has none) — the
        // candidate suggestion itself is the separate static hint in the section.
        text:
          'Bioavailability (%F) needs a per-subject Sex column to stratify the IV ' +
          'reference (pooling sexes biases %F by up to ~36%). Pick the Sex column ' +
          'in the Subject Sex section, or add one if your data has none.',
      });
    }
    closeBioavailability(view);
    return;
  }
  let result: BioavailabilityResult;
  try {
    result = computeBioavailability(results, {doseCol: doseColName, stratifyCols: [sexName]});
  } catch (e) {
    // Surface WHY F% didn't appear instead of failing silently (Ed BUG-06
    // secondary defect) — the route+dose+Sex preconditions held, so a failure here
    // is unexpected and worth showing.
    grok.shell.warning(`Bioavailability (%F) could not be computed: ${(e as Error)?.message ?? e}`);
    closeBioavailability(view);
    return;
  }
  if (result.table.rowCount === 0) {
    // No EV strata → no F% tab → banner over the results grid (B3 fallback).
    collectNotice(view, {
      key: 'F_PCT_NO_ROWS', anchor: 'bioavailability', severity: 'info',
      text:
        'Bioavailability (%F) produced no rows — the frame has no extravascular ' +
        '(non-IV) strata to divide by an IV reference. No %F table to show.',
    });
    closeBioavailability(view);
    return;
  }
  if (result.warnings.length > 0) {
    collectNotice(view, {
      key: 'F_PCT_WARNINGS', anchor: 'bioavailability', severity: 'warning',
      text: result.warnings.map((w) => w.message).join(' '),
    });
  }
  swapOrDockBioavailability(view, result.table);
}

/** Dock (or swap in place) the F% grid as a FILL-tab on the results node. */
function swapOrDockBioavailability(view: DG.TableView, fpct: DG.DataFrame): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  // `view.temp[key]` returns dart-side **null** (not undefined) for a never-set
  // key, so guard with `!= null` — the old `!== undefined` was true for that null
  // and ran `null.dataFrame = fpct` on the first Compute → "Cannot set properties
  // of null (setting 'dataFrame')", swallowed by applyBioavailability's catch, so
  // the tab never docked (Ed re-validation 2026-06-10, BUG-06; same root cause as
  // the TK-panel dock fix in swapOrDockToxGrid).
  const existingGrid = tempBag[FPCT_GRID_TAG] as DG.Grid | null;
  if (existingGrid != null) {
    existingGrid.dataFrame = fpct;
    configureBioavailabilityGrid(existingGrid); // fresh frame → re-apply headers/format/tooltips
    return;
  }
  const grid = DGNS.Viewer.grid(fpct);
  configureBioavailabilityGrid(grid);
  const anchor = resultsAnchor(view);
  const node = view.dockManager.dock(attachBannerHost(grid.root), DGNS.DOCK_TYPE.FILL, anchor, 'Bioavailability');
  tempBag[FPCT_GRID_TAG] = grid;
  tempBag[FPCT_NODE_TAG] = node;
}

/** Close and forget the F% grid tab. No-op when none is docked. */
function closeBioavailability(view: DG.TableView): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  const node = tempBag[FPCT_NODE_TAG] as DG.DockNode | null;
  if (node != null) {
    try {
      view.dockManager.close(node);
    } catch {
      // Node already gone — ignore.
    }
  }
  tempBag[FPCT_GRID_TAG] = undefined;
  tempBag[FPCT_NODE_TAG] = undefined;
}

// =====================================================================
// Dose proportionality (Smith power model, FR-410..413 / UC-03).
// Faithful copy of the bioavailability dock lifecycle above: compute on the
// post-appendGroupingColumns frame, dock the verdict + dn grids as FILL-tabs,
// surface the honest-uncertainty warning stream, swap in place on re-Compute.
// =====================================================================

/**
 * Dock the dose-proportionality grids when the frame is a dose-ranging design
 * (a per-subject dose column + a Sex metadata column + ≥ 2 distinct doses).
 * The per-Route×Sex Smith verdict (β, 90% CI, mode) is the headline; the dn-AUC
 * table is the secondary surface (AD-9). The P4 charts are design-gated.
 */
function applyDoseProportionality(
  view: DG.TableView,
  results: DG.DataFrame,
  sexName: string | null,
  doseColName: string | null,
  resultsNode: DG.DockNode | null,
  aucBasis: AucBasis = 'aucinf',
): void {
  const doseCol = doseColName !== null ? results.col(doseColName) : null;
  // Dose proportionality needs a per-subject dose, a Sex stratifier (the headline
  // is per Route×Sex), and an actual dose range. The no-dose / no-sex cases are
  // self-evidently inapplicable, so they close silently.
  if (doseColName === null || sexName === null || doseCol === null) {
    closeDoseProportionality(view);
    return;
  }
  const distinctDoses = new Set<number>();
  for (let i = 0; i < results.rowCount; i++) {
    if (doseCol.isNone(i)) continue;
    const d = Number(doseCol.get(i));
    if (Number.isFinite(d) && d > 0) distinctDoses.add(d);
  }
  if (distinctDoses.size < 2) {
    // Single dose level — not a dose-ranging study; nothing to analyze.
    closeDoseProportionality(view);
    return;
  }

  let result: DoseProportionalityResult;
  try {
    result = computeDoseProportionality(results, {
      doseCol: doseColName,
      sexCol: sexName,
      // ROUTE is always written by the results-writer; a single-route frame
      // collapses to one route stratum (no harm).
      routeCol: 'ROUTE',
      // Composite/sparse frames have a uniformly-NaN AUCinf (no fittable terminal
      // phase), so the default 'aucinf' basis self-reports "inapplicable" even
      // though valid AUClast exists. Caller passes 'auclast' there — the TK
      // convention (AUCτ basis) for repeat-dose dose-proportionality.
      basis: aucBasis,
    });
  } catch (e) {
    grok.shell.warning(`Dose proportionality could not be computed: ${(e as Error)?.message ?? e}`);
    closeDoseProportionality(view);
    return;
  }
  if (result.table.rowCount === 0) {
    closeDoseProportionality(view);
    return;
  }
  // Surface the honest-uncertainty stream as a SHORT summary toast: the same
  // code repeats across strata (e.g. SMITH_INAPPLICABLE per Route×Sex×param), so
  // a joined paragraph is unreadable and the transient balloon scrolls past. We
  // dedupe by code with a count; the full per-stratum sentences stay in
  // result.warnings, and the verdict grid carries per-stratum VERDICT/MODE.
  if (result.warnings.length > 0) {
    const counts = new Map<string, number>();
    for (const w of result.warnings) counts.set(w.code, (counts.get(w.code) ?? 0) + 1);
    const parts = [...counts.entries()].map(([code, n]) =>
      `${shortDoseProportionalityNote(code)}${n > 1 ? ` ×${n}` : ''}`);
    collectNotice(view, {
      key: 'DOSE_PROP_NOTES', anchor: 'doseProp', severity: 'warning',
      // No "see the … tab" trailer: the banner already sits over the Dose
      // Proportionality grid (and on the B3 fallback there is no such tab).
      text: `Dose proportionality · ${parts.join(', ')}`,
    });
  }
  swapOrDockDoseProportionality(view, result.table, result.dnTable, resultsNode);
}

/** Dock (or swap in place) the dose-proportionality verdict + dn grids as
 *  FILL-tabs on the results node. Mirrors {@link swapOrDockBioavailability}. */
function swapOrDockDoseProportionality(
  view: DG.TableView,
  verdict: DG.DataFrame,
  dn: DG.DataFrame,
  resultsNode: DG.DockNode | null,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  // `!= null` (not `!== undefined`): `view.temp[key]` returns dart-side null for
  // a never-set key (BUG-06); the old `!== undefined` ran `null.dataFrame = df`
  // on the first Compute and the tab never docked.
  const existingVerdict = tempBag[DP_GRID_TAG] as DG.Grid | null;
  const existingDn = tempBag[DPDN_GRID_TAG] as DG.Grid | null;
  if (existingVerdict != null) {
    existingVerdict.dataFrame = verdict;
    configureDoseProportionalityGrid(existingVerdict);
    if (existingDn != null) {existingDn.dataFrame = dn; configureDoseProportionalityDnGrid(existingDn);}
    return;
  }
  const anchor = resultsAnchor(view, resultsNode);
  const verdictGrid = DGNS.Viewer.grid(verdict);
  configureDoseProportionalityGrid(verdictGrid);
  const verdictNode = view.dockManager.dock(
    attachBannerHost(verdictGrid.root), DGNS.DOCK_TYPE.FILL, anchor, 'Dose Proportionality');
  tempBag[DP_GRID_TAG] = verdictGrid;
  tempBag[DP_NODE_TAG] = verdictNode;

  const dnGrid = DGNS.Viewer.grid(dn);
  configureDoseProportionalityDnGrid(dnGrid);
  const dnNode = view.dockManager.dock(
    dnGrid.root, DGNS.DOCK_TYPE.FILL, verdictNode ?? anchor, 'Dose-Normalized');
  tempBag[DPDN_GRID_TAG] = dnGrid;
  tempBag[DPDN_NODE_TAG] = dnNode;
}

/** Close and forget the dose-proportionality grid tabs. No-op when none docked. */
function closeDoseProportionality(view: DG.TableView): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  for (const nodeTag of [DP_NODE_TAG, DPDN_NODE_TAG]) {
    const node = tempBag[nodeTag] as DG.DockNode | null;
    if (node != null) {
      try {
        view.dockManager.close(node);
      } catch {
        // Node already gone — ignore.
      }
    }
  }
  tempBag[DP_GRID_TAG] = undefined;
  tempBag[DP_NODE_TAG] = undefined;
  tempBag[DPDN_GRID_TAG] = undefined;
  tempBag[DPDN_NODE_TAG] = undefined;
}

// =====================================================================
// UC-04 sparse / destructive-sampling NCA panels.
// Faithful copy of the dose-proportionality dock lifecycle: detect applicability
// on a successful Compute, run `computeSparseNca` over the RAW observations
// (sparse is a cross-animal composite — it needs the per-animal samples, not the
// per-subject results frame), dock the VR-07 composite scatter + the honest-
// uncertainty parameter grid + the composite-profile grid next to results.
// =====================================================================

// `denseProfilesViable` moved to adapter/review-units.ts (nca-workflow-stages F-S1.3) so the
// review-unit classifier and the hub-side stage read the SAME predicate; imported above.

// `findColByPattern`, `isSparseDesign`, the sparse compute (`computeSparseForRecipe`) and
// `sparseCompositeUnits` moved to the adapter (nca-pipeline-save-rework F1) so the pipeline
// save computes the same composites; imported above.

/**
 * Detect a sparse / destructive design on the source frame and, when present,
 * compute + dock the sparse panels. Mirrors {@link applyDoseProportionality}.
 * Closes any prior sparse tabs when the current mapping is not a sparse design.
 *
 * @param mapping subject/time/conc names + units from the mapping form.
 */
function applySparseNca(
  view: DG.TableView,
  df: DG.DataFrame,
  mapping: SparseComputeMapping,
  resultsNode: DG.DockNode | null,
  useLogY: boolean,
): SparseNcaResult | null {
  const {result, blocked, error, elapsedMs} = computeSparseForRecipe(df, mapping);
  if (error !== null) {
    // Sparse tab is not docked on failure → banner over the results grid (B3).
    collectNotice(view, {
      key: 'SPARSE_COMPUTE_FAILED', anchor: 'sparse', severity: 'warning',
      text: `Sparse NCA could not be computed: ${error}`,
    });
  }
  // BUG-19 AD-2: a stratum that was not computed (a sparse IV-infusion route; a dose group with no
  // sample left) is a notice, never a silent bolus or a vanished composite. One line per stratum.
  for (const b of blocked) {
    collectNotice(view, {
      key: 'SPARSE_ROUTE_BLOCKED', anchor: 'sparse', severity: 'warning',
      text: `Sparse NCA: ${b.route}${b.doseKey !== '' ? `, dose ${b.doseKey}` : ''} not computed — ${b.reason}.`,
    });
  }
  // AC-C3 / AD-12: measure the synchronous bootstrap wall-clock against NFR-03
  // (< 10 s). Over budget ⇒ advise wiring the worker pool (a measured trigger,
  // not a time-based deferral). Sparse N is small by design, so v1 runs sync.
  if (result !== null && elapsedMs > 10000) {
    collectNotice(view, {
      key: 'SPARSE_BOOTSTRAP_SLOW', anchor: 'sparse', severity: 'warning',
      text: `Sparse NCA + bootstrap took ${(elapsedMs / 1000).toFixed(1)} s ` +
        `(> NFR-03 10 s budget) — consider the worker pool for this dataset size.`,
    });
  }
  if (result === null || result.groups.length === 0) {
    closeSparseNca(view);
    return null;
  }

  // Honest-uncertainty stream — a SHORT toast summary (the transient balloon
  // can't hold the full sentences); the persistent NOTES column of the Sparse
  // NCA grid carries the per-group advisories, never silent (thesis #4).
  const distinctCodes = [...new Set(result.warnings.map((w) => w.code))];
  if (distinctCodes.length > 0) {
    const labels = distinctCodes.map((c) => shortSparseNote(c));
    collectNotice(view, {
      key: 'SPARSE_NOTES', anchor: 'sparse', severity: 'warning',
      // No "see the NOTES column" trailer: the banner sits over the Sparse grid,
      // whose NOTES column carries the per-group detail right below it.
      text: `Sparse NCA · ${distinctCodes.length} note(s): ${labels.join(', ')}`,
    });
  }

  ensureSparseChart(view, result, useLogY);
  swapOrDockSparse(view, result, resultsNode, mapping.acceptance);
  return result;
}

/** view.temp key: the sparse parameter grid's current-row subscription (composite review). */
const SPARSE_ROW_SUB_TAG = 'ncaSparseRowSub';
/** view.temp key: navigate the Sparse NCA table to a composite by its review key (`sparse:<g>`). */
const SPARSE_GOTO_TAG = 'ncaSparseGoto';
/** view.temp keys: the composite group keys aligned with the sparse parameter-table rows, and
 *  each row's Exclude λz decision (or null) — recomputed on every review-state refresh. */
const SPARSE_GROUP_KEYS_TAG = 'ncaSparseGroupKeys';
const SPARSE_EXCLUDED_TAG = 'ncaSparseExcluded';

/**
 * Recompute which sparse composites are EXCLUDED and repaint their table: Exclude λz on a composite
 * must be as real as on a dense profile, so its CL / t½ / Vz are marked too. The
 * mask comes from the controller's `excludedKeys()` — the SAME set the dense reportable feed
 * reads — via the adapter's {@link excludedCompositeRows}; the decision record rides along for
 * the tooltip. No-op before a sparse compute.
 */
function refreshSparseExclusions(view: DG.TableView): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const keys = bag[SPARSE_GROUP_KEYS_TAG] as string[] | null;
  const controller = bag[NCA_CONTROLLER_TAG] as NcaController | null;
  if (keys == null || controller == null) return;
  const mask = excludedCompositeRows(keys, controller.excludedKeys());
  const byKey = new Map(liveWorkflowState(view).units.map((u) => [u.key, u]));
  bag[SPARSE_EXCLUDED_TAG] = keys.map((g, i) => mask[i] ? (byKey.get(sparseUnitKey(g))?.decision ?? null) : null);
  const grid = bag[SPARSE_PARAM_GRID_TAG] as DG.Grid | null;
  if (grid == null) return;
  fitExcludedColumns(grid, SPARSE_LZ_DEPENDENT_COLS, (i) => mask[i] === true);
  grid.invalidate();
}

/**
 * The context panel for ONE UC-04 group composite (F-S4.8 / U3 / SF-3): the composite's verdict,
 * led by the pooled-cohort cue (AC-S4.8.0), and the SAME Review pane a dense profile gets —
 * Accept as is / Exclude λz / Acknowledge with a reason, never point editing.
 */
function showCompositeReview(
  view: DG.TableView, g: SparseNcaResult['groups'][number], cohort: SparseCohort,
  thresholds: ResultsColorThresholds,
  /** Move to another composite (its key) — `Next to review` walks the pending composites. */
  goTo: (key: string) => void,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const controller = bag[NCA_CONTROLLER_TAG] as NcaController | null;
  if (controller == null) return;
  const render = (): void => {
    const {stage, units, fingerprints} = liveWorkflowState(view);
    const key = sparseUnitKey(g.groupKey);
    const unit = units.find((u) => u.key === key);
    const title = ui.divText(g.groupKey, {style: {fontSize: '18px', color: 'var(--grey-4, #9497A0)'}});
    const verdict = ui.divText(g.acceptance.verdict === ACCEPT_FLAGGED ?
      g.acceptance.reasons.join(' · ') : g.acceptance.verdict, {style: {fontWeight: 'bold', whiteSpace: 'normal'}});
    const cue = ui.divText(pooledCohortCue(cohort), {style: {whiteSpace: 'normal', color: 'var(--grey-5, #7D8086)'}});
    // Full-width host + accordion, exactly as the dense λz readout (`showLzReadout`): without a
    // definite width the Review pane's `ui.form` collapses to its inline layout (label and
    // checkbox run together — manual QA 2026-09-30).
    const reviewHost = ui.div([]);
    reviewHost.style.width = '100%';
    const acc = ui.accordion('nca-composite-review');
    acc.root.style.width = '100%';
    // AC-S4.8.0: the pooled-cohort cue LEADS the verdict — never read as one animal's profile.
    acc.addTitle(ui.divV([title, cue, verdict], {style: {gap: '6px'}}));
    acc.addPane('Review', () => reviewHost, true);
    if (unit !== undefined) {
      const pendingComposites = units.filter((u) => u.kind === 'composite' && u.need === 'review' &&
        u.decision === null).map((u) => u.key);
      const after = (): void => {
        markAnalysisDirty(view);
        refreshStageUi(view);
        render();
      };
      renderReviewPane(reviewHost, {
        unit, selected: [], pendingOrder: pendingComposites, thresholds,
        editable: viewStageOf(view) === 'review' && !stage.setupChanged,
        pipelineLine: pipelineLineFor(view, unit.key),
        onDecide: (keys, action, reasons) => {
          for (const k of keys) controller.decide(k, action, reasons, fingerprints.get(k) ?? '');
          after();
        },
        onReset: (k) => {controller.resetDecision(k); after();},
        onNext: () => {
          if (pendingComposites.length === 0) return;
          const at = pendingComposites.indexOf(key);
          goTo(pendingComposites[(at + 1) % pendingComposites.length] ?? pendingComposites[0]);
        },
      });
    }
    setTimeout(() => {
      try {grok.shell.windows.showContextPanel = true; grok.shell.o = acc.root;} catch {/* no shell */}
    }, 0);
  };
  render();
}

/** The pooled-cohort shape of each composite, aligned with the parameter-table rows (AC-S4.8.0). */
function sparseCohortsOf(result: SparseNcaResult): SparseCohort[] {
  return result.groups.map((g) => ({
    nAnimals: g.nAnimals, minPerTimepoint: g.minTimepointN, maxPerTimepoint: g.maxTimepointN,
    reasons: g.acceptance.reasons,
    // The C0_EXTRAP tooltip: the method of the shown anchor (`sparseAuc`), the share from the dense run.
    c0: g.c0 !== null ? {method: g.c0.method, pctAucBackExtrap: g.c0Provenance?.pctAucBackExtrap ?? NaN} : null,
  }));
}

/** Dock (or swap) the sparse parameter grid + composite grid as FILL-tabs on the
 *  results node. Mirrors {@link swapOrDockDoseProportionality}. */
function swapOrDockSparse(
  view: DG.TableView, result: SparseNcaResult, resultsNode: DG.DockNode | null,
  /** The thresholds the composites were judged with — the grid's adj-R² band sits at the SAME
   *  `R_acc` the verdict used. */
  acceptance: AcceptanceThresholds,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  const existingParam = tempBag[SPARSE_PARAM_GRID_TAG] as DG.Grid | null;
  const existingComp = tempBag[SPARSE_COMP_GRID_TAG] as DG.Grid | null;
  // UC-04 composite review (F-S4.8 / U3): the sparse grid's current row opens the SAME Review pane
  // for that group composite — no point editing. Re-bound per compute (a fresh parameter frame).
  (tempBag[SPARSE_ROW_SUB_TAG] as {unsubscribe(): void} | null)?.unsubscribe?.();
  const thresholds = {...DEFAULT_RESULTS_THRESHOLDS, ...acceptance};
  const cohorts = sparseCohortsOf(result);
  tempBag[SPARSE_GROUP_KEYS_TAG] = result.groups.map((g) => g.groupKey);
  tempBag[SPARSE_EXCLUDED_TAG] = null;
  // The pipeline run's own reasons per composite (F10), keyed like the app's composite units.
  writePipelineColumns(result.parameterTable, pipelineReasonsOf(view),
    (i) => sparseUnitKey(result.groups[i].groupKey));
  const excludedAt = (i: number): ReviewDecisionRecord | null =>
    (tempBag[SPARSE_EXCLUDED_TAG] as (ReviewDecisionRecord | null)[] | null)?.[i] ?? null;
  // `Next to review` on a composite (F-S4.3): move the parameter table's current row to that
  // composite — the row subscription below then opens its Review pane.
  const goToComposite = (key: string): void => {
    const idx = result.groups.findIndex((g) => sparseUnitKey(g.groupKey) === key);
    if (idx < 0) return;
    // Already the current row → no row-change event; open its pane directly.
    if (result.parameterTable.currentRowIdx === idx)
      showCompositeReview(view, result.groups[idx], cohorts[idx], thresholds, goToComposite);
    else result.parameterTable.currentRowIdx = idx;
  };
  tempBag[SPARSE_GOTO_TAG] = goToComposite; // startReview's entry into a composite
  tempBag[SPARSE_ROW_SUB_TAG] = result.parameterTable.onCurrentRowChanged.subscribe(() => {
    const i = result.parameterTable.currentRowIdx;
    const g = i >= 0 ? result.groups[i] : undefined;
    if (g !== undefined) showCompositeReview(view, g, cohorts[i], thresholds, goToComposite);
  });
  if (existingParam != null) {
    existingParam.dataFrame = result.parameterTable;
    configureSparseNcaGrid(existingParam, thresholds, cohorts, excludedAt);
    applySparseReasonPills(existingParam, reasonPillsShown(view));
    refreshSparseExclusions(view);
    if (existingComp != null) {
      existingComp.dataFrame = result.compositeTable;
      configureSparseCompositeGrid(existingComp);
    }
    return;
  }
  const anchor = resultsAnchor(view, resultsNode);
  const paramGrid = DGNS.Viewer.grid(result.parameterTable);
  configureSparseNcaGrid(paramGrid, thresholds, cohorts, excludedAt);
  applySparseReasonPills(paramGrid, reasonPillsShown(view));
  const paramNode = view.dockManager.dock(
    attachBannerHost(paramGrid.root), DGNS.DOCK_TYPE.FILL, anchor, 'Sparse NCA');
  tempBag[SPARSE_PARAM_GRID_TAG] = paramGrid;
  tempBag[SPARSE_PARAM_NODE_TAG] = paramNode;
  refreshSparseExclusions(view);

  const compGrid = DGNS.Viewer.grid(result.compositeTable);
  configureSparseCompositeGrid(compGrid);
  const compNode = view.dockManager.dock(
    attachBannerHost(compGrid.root), DGNS.DOCK_TYPE.FILL, paramNode ?? anchor, 'Composite Profile');
  tempBag[SPARSE_COMP_GRID_TAG] = compGrid;
  tempBag[SPARSE_COMP_NODE_TAG] = compNode;
}

/**
 * VR-07 — composite scatter of individuals + a connected composite-mean line.
 * Reuses the DG scatterPlot primitive: raw per-animal
 * points are one colour series per dose group; the composite mean per nominal
 * time is overlaid as a distinct series (`D{dose} mean`) AND connected into a
 * per-dose line via the viewer's native `linesByColumnName` / `linesOrderColumnName`.
 *
 * The line is drawn only through the means: the `LineGroup` column gives every
 * mean point of a dose the SAME group (`D{dose}` → one polyline) and every
 * raw point a UNIQUE group (`~raw_{k}` → a singleton, so no line) — raw points
 * stay isolated markers (connecting them would imply a single animal's profile,
 * which is false in a destructive design). Lines are ordered by nominal time.
 *
 * Redraw-safe: `chartDf` is rebuilt and `options` (incl. the line settings) is
 * re-applied on BOTH paths — viewer creation and the `existing != null` re-bind
 * (`setOptions`) — so the line survives a re-Compute.
 */
/** The composite chart's x column: nominal time in the composite engine's unit (BUG-19 AC-6.1). */
const SPARSE_CHART_X = 'Nominal time (h)';

function ensureSparseChart(view: DG.TableView, result: SparseNcaResult, useLogY: boolean): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  const nomTime: number[] = []; const conc: number[] = [];
  const dose: string[] = []; const kind: string[] = []; const lineGroup: string[] = [];
  let rawSeq = 0;
  // The dose label is the dose level as the source spells it (BUG-19 AD-3), the route leading it
  // when the analysis holds more than one route.
  const multiRoute = new Set(result.groups.map((g) => g.route)).size > 1;
  for (const g of result.groups) {
    const label = `${multiRoute ? `${g.route} ` : ''}D${g.doseLevel.key}`;
    for (const p of g.rawPoints) {
      if (!Number.isFinite(p.conc)) continue;
      nomTime.push(p.nominalTime); conc.push(p.conc); dose.push(label); kind.push('raw');
      lineGroup.push(`~raw_${rawSeq++}`); // '~' sorts after 'D' ⇒ keeps D10/D50 at palette 0/1; unique ⇒ no line
    }
    for (const tp of g.composite.timepoints) {
      if (!(tp.mean > 0)) continue;
      nomTime.push(tp.time); conc.push(tp.mean); dose.push(label); kind.push('mean');
      lineGroup.push(label); // 'D10'/'D50' — same value as `dose` ⇒ the mean line gets the dose colour
    }
  }
  // Pre-derived semi-log Y (PK convention): log10(conc) for conc > 0, else
  // FLOAT_NULL — so the renderer never sees log(0)/log(neg) and BLQ/zero points
  // simply drop out (the mean line honestly breaks at a below-LLOQ timepoint).
  // Mirrors the dense profile chart's `log(<conc>)` column approach.
  const logConc = new Float32Array(conc.length);
  for (let i = 0; i < conc.length; i++) {
    // FLOAT_NULL is a tiny POSITIVE sentinel — guard it explicitly so it can't
    // pass `> 0` and plot a spurious log10(sentinel) point (matches the dense
    // chart's buildChartDataFrame guard).
    logConc[i] = log10OrNull(conc[i]);
  }
  // The composite engine reads nominal time in hours whatever the input unit (BUG-19 AC-6.1), so
  // the x axis names its unit.
  const chartDf = DGNS.DataFrame.fromColumns([
    DGNS.Column.fromFloat64Array(SPARSE_CHART_X, Float64Array.from(nomTime)),
    DGNS.Column.fromFloat64Array('Conc', Float64Array.from(conc)),
    DGNS.Column.fromFloat32Array('log(Conc)', logConc),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'Dose', dose),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'Kind', kind),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'LineGroup', lineGroup),
  ]);
  // Freeze a NATURAL category order on the colour column so the legend reads dose
  // labels in magnitude order (D10, D30, D100, D300), not the string-lexical
  // default (D10, D100, D30, D300). Platform's own "Sort → Natural" mechanism —
  // `setCategoryOrder` persists the `.category-order` tag viewers honour. Metadata
  // only: no plotted value or palette-consistency change (markers + mean lines still
  // key off the same 'D{dose}' categories). See buildMeanSdProfile for the rationale.
  const doseCol = chartDf.col('Dose')!;
  doseCol.setCategoryOrder(
    doseCol.categories.slice().sort((a, b) => a.localeCompare(b, undefined, {numeric: true})));
  // Colour by Dose (2 categories) and connect the means by LineGroup, whose mean
  // categories are the SAME 'D10'/'D50' strings → they land on the SAME alphabetical
  // palette indices as Dose, so each mean line matches its dose's marker colour.
  // Raw vs mean is shown by marker shape (Kind), not colour. (DG's categorical
  // palette is assigned in alphabetical category order — verified empirically.)
  const options = {
    x: SPARSE_CHART_X, y: useLogY ? 'log(Conc)' : 'Conc',
    colorColumnName: 'Dose',
    markersColumnName: 'Kind',
    showRegressionLine: false,
    linesByColumnName: 'LineGroup',
    linesOrderColumnName: SPARSE_CHART_X,
    linesWidth: 2,
    title: useLogY ?
      'Composite profile (semi-log) — individuals + mean line' :
      'Composite profile — individuals + mean line',
  };
  const existing = tempBag[SPARSE_CHART_TAG] as DG.Viewer | null;
  if (existing != null) {
    existing.dataFrame = chartDf;
    existing.setOptions(options);
    return;
  }
  const chart = DGNS.Viewer.scatterPlot(chartDf, options);
  // Charts group, beside Combined / Facet (user decision 2026-08-13): it plots the
  // same concentration–time axes they do. Its two companion TABLES (`Sparse NCA`,
  // `Composite Profile`) stay in the tables group — the split follows the
  // surface, not the feature.
  const node = dockChartTab(view, chart.root, 'Composite');
  tempBag[SPARSE_CHART_TAG] = chart;
  tempBag[SPARSE_CHART_NODE_TAG] = node;
}

/** Close and forget every sparse tab. No-op when none docked. */
function closeSparseNca(view: DG.TableView): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  for (const nodeTag of [SPARSE_PARAM_NODE_TAG, SPARSE_COMP_NODE_TAG, SPARSE_CHART_NODE_TAG]) {
    const node = tempBag[nodeTag] as DG.DockNode | null;
    if (node != null) {
      try {
        view.dockManager.close(node);
      } catch {
        // Node already gone — ignore.
      }
    }
  }
  const chart = tempBag[SPARSE_CHART_TAG] as DG.Viewer | null;
  if (chart != null) {
    try {
      chart.close();
    } catch {
      // Already closed — ignore.
    }
  }
  // The composite-review row subscription holds the dead parameter table — release it too.
  (tempBag[SPARSE_ROW_SUB_TAG] as {unsubscribe(): void} | null)?.unsubscribe?.();
  tempBag[SPARSE_ROW_SUB_TAG] = undefined;
  tempBag[SPARSE_GROUP_KEYS_TAG] = undefined;
  tempBag[SPARSE_EXCLUDED_TAG] = undefined;
  tempBag[SPARSE_GOTO_TAG] = undefined;
  tempBag[SPARSE_PARAM_GRID_TAG] = undefined;
  tempBag[SPARSE_PARAM_NODE_TAG] = undefined;
  tempBag[SPARSE_COMP_GRID_TAG] = undefined;
  tempBag[SPARSE_COMP_NODE_TAG] = undefined;
  tempBag[SPARSE_CHART_TAG] = undefined;
  tempBag[SPARSE_CHART_NODE_TAG] = undefined;
}

// =====================================================================
// VR-03 — group mean ± SD concentration-time profile chart.
// Pools the normalized observations frame into a (group × time) mean ± SD
// profile (`buildMeanSdProfile`) and renders Mean + ±SD lines per group, coloured
// by group, in a FILL-tab. The data primitive is the adapter helper; this is the
// wiring + long-frame render (mirrors `ensureSparseChart`'s scatterPlot idiom).
// =====================================================================

/** Render context captured at Compute so the Group-by control can re-split the
 *  "Mean ± SD" chart without a recompute. `denseDf` is the same normalized
 *  observations frame the Combined/Facet tabs render from. `nominalUsed` records
 *  whether `timeName` is a nominal-time key (a mapped nominal-time column, or TAD
 *  in NONMEM mode) — when true the serial-drift fragmentation toast is suppressed
 *  (binning is already on a common grid). */
interface MeanSdCtx {
  readonly denseDf: DG.DataFrame;
  readonly timeName: string;
  readonly concName: string;
  readonly doseName: string | null;
  readonly nominalUsed: boolean;
}

/** Synthetic constant-group column name — a single pooled "All" curve when no
 *  group key resolves (correct for single-dose data like theoph: one mean ± SD
 *  across all subjects IS the canonical VR-03 view). */
const MEAN_SD_ALL_GROUP = 'NcaGroup';

/** Composite Group-by key column synthesized on a denseDf clone for the mean±SD
 *  chart when ≥2 grouping columns are selected (v8, multi-column-group-by). */
const MEAN_SD_COMPOSITE_COL = '__nca_meansd_group';

/** Human-readable separator for DISPLAY composite labels (legend, notice text) — the
 *  readable counterpart to the internal U+0001 keying separator. */
const GROUP_LABEL_SEP = ' · ';

/**
 * Collapse a Group-by picker value to the SET of distinct column names it means, as columns
 * of `frame` (first-seen order = the user's pick order; a name absent from `frame` is dropped).
 *
 * **Why both halves are load-bearing** (BUG-11, 2026-07-17). The selection is semantically a
 * SET: grouping by [A, B, A] partitions identically to [A, B], so a repeated name adds no
 * information — it only doubles the key in every LABEL (ribbon caption, composite chart legend,
 * aggregate group label) and, once `getCurrentRecipe` persists it into `recipe.groupByCols`,
 * replays on every reopen.
 *  - **Dedupe by name:** the stock `ui.input.columns` dialog lists its rows as
 *    `value ++ tableColumns` but ticks the checkboxes BY NAME, so every already-selected name
 *    is ticked TWICE (once in its `value` row, once in its table row) and a plain OK commits
 *    the union. Measured on the user's study: a no-op open→OK grew the selection 5 → 8 (+3 =
 *    the 3 distinct names) and the legend to 8 segments. Unbounded on repeat.
 *  - **Re-resolve against `frame`:** a value column carried over from a PREVIOUS Compute's
 *    frame is not identity-equal to this frame's same-named column, so the dialog lists it as
 *    an EXTRA row on top of the table's own — which is what seeds the duplicate in the first
 *    place. Handing the widget only live columns keeps its list one-row-per-name.
 *
 * `buildSummaryCache` already reduces its own key list this way (`Array.from(new Set(…))`),
 * which is why the aggregate NUMBERS were never wrong — a duplicated key is a no-op partition —
 * and the damage stayed in the labels + the persisted recipe. This is that same invariant,
 * enforced at the pane's two write sites instead of downstream of them.
 *
 * Module-level + exported (not a `mountMappingPane` closure) because it depends on nothing but
 * its arguments and carries the invariant the pane is judged on — the tests pin it against a
 * real frame.
 */
export function dedupeGroupBySelection(
  names: ReadonlyArray<string>, frame: DG.DataFrame,
): DG.Column[] {
  const out: DG.Column[] = [];
  const seen = new Set<string>();
  for (const nm of names) {
    if (seen.has(nm)) continue;
    const live = frame.col(nm);
    if (live === null) continue;
    seen.add(nm);
    out.push(live);
  }
  return out;
}

/** A (group × time) cell with N below this is "fragile" — its mean is a small-sample
 *  estimate whose limits must be surfaced, not hidden (thesis §4). Set to 4 to MATCH
 *  the shipped study-overview grid precedent (`nca-study-overview.ts` bold-fonts `n < 4`,
 *  flagging N=1/2/3), so the same "fragile estimate" concept reads identically across
 *  the study grid, the aggregate grid, and the Mean±SD banner (cross-view consistency —
 *  the spec's F4 bullet names that precedent as the reuse target). */
export const FRAGILE_N = 4;
/** Above this many fragile cells the fragility banner COLLAPSES to a count + worst-cell
 *  summary instead of enumerating — avoids the gap-18 noisy-warning flood at high
 *  composite cardinality (each added key multiplies buckets toward N=1). Absolute count
 *  ONLY — no fraction trigger (R2 NEW-3: a fraction with no floor fires at trivial
 *  cardinality). The enumerate/collapse/no-collapse branches are unit-tested via the
 *  exported {@link fragilityNoticeText} (multi-column-group-by-tests.ts). */
export const FRAGILE_LIST_MAX = 6;

/** One fragile (group × time) cell for the mean±SD fragility banner. `group` is the
 *  (possibly composite, U+0001-separated) group key. */
export interface FragileCell {readonly group: string; readonly time: number; readonly n: number;}

/** Humanize a (possibly composite) group key (U+0001 → ` · `) for a notice / label. */
function humanizeGroupKey(group: string): string {
  return group.split(KEY_SEP).join(GROUP_LABEL_SEP);
}

/** The MEAN_SD_FRAGMENTED banner text (F4, multi-column-group-by). Enumerates the
 *  fragile (group × time) cells worst-first, or collapses to a count + worst cell past
 *  {@link FRAGILE_LIST_MAX}. Retains the actionable "map a nominal-time column" advice
 *  when pooling on OBSERVED time (the observed-time serial-drift case). */
export function fragilityNoticeText(cells: ReadonlyArray<FragileCell>, nominalUsed: boolean): string {
  const sorted = cells.slice().sort((a, b) => a.n - b.n); // worst (smallest N) first
  const fmt = (c: FragileCell): string => `N=${c.n} for ${humanizeGroupKey(c.group)} at ${c.time} h`;
  const body = sorted.length > FRAGILE_LIST_MAX ?
    `${sorted.length} (group × time) cells have N < ${FRAGILE_N}; worst ${fmt(sorted[0])} — ` +
      `small-N group means are fragile.` :
    `${sorted.map(fmt).join('; ')} — small-N group means are fragile.`;
  const advice = nominalUsed ? '' :
    ' Map a nominal-time column (Sparse / destructive design ▸ Nominal time) so serial ' +
    'samples align onto a common time grid.';
  return `Mean ± SD: ${body}${advice}`;
}


/** Y-axis [min, max] (in DATA units) that fully contains the ±SD whiskers so the
 *  caps don't clip against the plot edge — the viewer otherwise auto-scales to the
 *  mean only. Reads the band-edge columns for the CURRENT scale: `Lower`/`Upper`
 *  (linear, arithmetic) or `geoLower`/`geoUpper` (semi-log, geometric — always
 *  positive, so a log axis is well-defined). Linear uses additive padding; log
 *  uses multiplicative padding (symmetric in log space). Returns null when no
 *  finite band exists (caller defers to viewer auto-range). */
function meanSdYRange(df: DG.DataFrame, useLogY: boolean): {min: number; max: number} | null {
  const loE = finiteMinMax(df.col(useLogY ? 'geoLower' : 'Lower'));
  const upE = finiteMinMax(df.col(useLogY ? 'geoUpper' : 'Upper'));
  if (loE === null || upE === null) return null;
  const mn = loE.min; const mx = upE.max;
  if (mn > mx) return null;
  if (useLogY) return {min: mn * 0.9, max: mx * 1.1}; // geo edges > 0 → log-safe
  const pad = (mx - mn) * 0.05 || Math.abs(mx) * 0.05 || 1;
  return {min: mn - pad, max: mx + pad};
}

/** X-axis [min, max] pinned to the finite Time extent (+2% additive pad) so the
 *  viewer doesn't pad into negative time / leave a dead zone where stray markers
 *  render. Re-evaluated per scale because the log frame drops all-BLQ rows (often
 *  the latest timepoints), shrinking the X extent. Null → defer to auto-range. */
function meanSdXRange(df: DG.DataFrame): {min: number; max: number} | null {
  const e = finiteMinMax(df.col('Time'));
  if (e === null || e.max <= e.min) return null;
  const pad = (e.max - e.min) * 0.02;
  return {min: e.min - pad, max: e.max + pad};
}

/** Per-scale chart frame. Linear → the full profile (BLQ/zero cells plot legibly
 *  at y=0). Semi-log → a clone with all-BLQ cells (`geoMean` = FLOAT_NULL,
 *  nPos=0) dropped: a log scatterPlot otherwise renders those null-y points at the
 *  axis floor and the line connects through them (bottom artefact). Rows with a
 *  geometric value (≥1 measurable sample) are kept — their geo band is positive. */
function meanSdChartFrame(profile: DG.DataFrame, useLogY: boolean): DG.DataFrame {
  if (!useLogY) return profile;
  const gm = profile.col('geoMean')!;
  const mask = DGNS.BitSet.create(profile.rowCount);
  for (let i = 0; i < profile.rowCount; i++) if (!gm.isNone(i)) mask.set(i, true);
  return profile.clone(mask);
}

/** Create / update / tear down the VR-03 group mean ± SD chart. `groupName` null
 *  (or absent from `denseDf` — e.g. a results-only key like STATUS) → a single
 *  pooled "All" curve. Returns whether this call CREATED the tab and the median N
 *  per (group × time) cell, so the caller can surface the serial-drift toast on
 *  first render. */
function ensureMeanSdChart(
  view: DG.TableView, ctx: MeanSdCtx, groupName: string | null, useLogY: boolean,
): {created: boolean; medianN: number; groupCount: number; minCellN: number;
    fragileCells: FragileCell[]} {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  const {denseDf, timeName, concName} = ctx;
  if (denseDf.col(timeName) === null || denseDf.col(concName) === null) {
    closeMeanSdChart(view);
    return {created: false, medianN: 0, groupCount: 0, minCellN: 0, fragileCells: []};
  }

  // Resolve the group column; synthesize a constant "All" group when none applies
  // (never silently pool across dose levels — the caller resolves dose first).
  let workDf = denseDf; let gName = groupName;
  if (gName === null || denseDf.col(gName) === null) {
    workDf = denseDf.clone();
    workDf.columns.add(DGNS.Column.fromList(
      DGNS.COLUMN_TYPE.STRING, MEAN_SD_ALL_GROUP,
      new Array<string>(workDf.rowCount).fill('All')));
    gName = MEAN_SD_ALL_GROUP;
  }

  const profile = buildMeanSdProfile(workDf, {timeCol: timeName, concCol: concName, groupCol: gName});
  if (profile.rowCount === 0) {
    closeMeanSdChart(view);
    return {created: false, medianN: 0, groupCount: 0, minCellN: 0, fragileCells: []};
  }

  // Render the wide (group × time) profile directly as a scatterPlot:
  //  - mean point per (group × time), connected into one line per group
  //    (`linesByColumnName='Group'` + order by Time);
  //  - colour by Group (`colorColumnName='Group'`) — markers AND the mean line
  //    share the group colour (both keyed off the same 'Group' categories →
  //    identical alphabetical palette indices);
  //  - **native ±SD whiskers** via `yWhiskerMin/MaxColumnName`, with a NATIVE
  //    logarithmic y-axis (`yAxisType`) for semi-log.
  // Scale-appropriate statistic (the science fix): the LINEAR axis shows the
  // ARITHMETIC mean ± SD (`Mean` / `Lower` / `Upper`); the SEMI-LOG axis shows the
  // GEOMETRIC mean ×/÷ gSD (`geoMean` / `geoLower` / `geoUpper`). Arithmetic
  // mean−SD goes ≤ 0 at high CV / late timepoints → undefined on a log axis; the
  // geometric band is always positive and is the PK-standard for a semi-log
  // profile (the correct measure per scale, not a render hack).
  // Per-scale chart frame (log drops all-BLQ rows so the viewer can't render them
  // at the axis floor). Store the full profile so the semi-log toggle can rebuild.
  tempBag[MEAN_SD_PROFILE_TAG] = profile;
  const chartDf = meanSdChartFrame(profile, useLogY);
  const yCol = useLogY ? 'geoMean' : 'Mean';
  const loCol = useLogY ? 'geoLower' : 'Lower';
  const upCol = useLogY ? 'geoUpper' : 'Upper';
  const range = meanSdYRange(chartDf, useLogY);
  const options: Record<string, unknown> = {
    x: 'Time',
    y: yCol,
    colorColumnName: 'Group',
    linesByColumnName: 'Group',
    linesOrderColumnName: 'Time',
    linesWidth: 2,
    yWhiskerMinColumnName: loCol,
    yWhiskerMaxColumnName: upCol,
    yAxisType: useLogY ? 'logarithmic' : 'linear',
    showRegressionLine: false,
    // Suppress scatterPlot edge decorations that render as stray strips/points at
    // the plot bottom/right (axis histograms, filtered-out ghosts, drop lines).
    showXHistogram: false,
    showYHistogram: false,
    showFilteredOutPoints: false,
    showDropLines: false,
    // User request 2026-06-23: hide the Size selector; show the caption as the
    // viewer DESCRIPTION (top) instead of a title bar.
    showSizeSelector: false,
    showTitle: false,
    description: useLogY ? 'Group geometric mean ×/÷ gSD (semi-log)' : 'Group mean ± SD',
    descriptionVisibilityMode: 'Always',
    descriptionPosition: 'Top',
  };
  // Pin the y-axis to the whisker extent so the ±SD caps don't clip against the
  // plot edge (the viewer otherwise auto-scales to the mean only).
  if (range !== null) {options.yMin = range.min; options.yMax = range.max;}
  const xRange = meanSdXRange(chartDf);
  if (xRange !== null) {options.xMin = xRange.min; options.xMax = xRange.max;}

  const existing = tempBag[MEAN_SD_CHART_TAG] as DG.Viewer | null;
  if (existing != null) {
    // TEMP: guard ONLY the frame swap — a different-row-count dataFrame leaves
    // DG's scatterPlot pointHitTest BitSet stale (a platform bug); setOptions stays
    // unguarded so real config bugs surface.
    try {existing.dataFrame = chartDf;} catch {/* stale hit-test BitSet */}
    existing.setOptions(options);
  } else {
    const chart = DGNS.Viewer.scatterPlot(chartDf, options);
    // FIRST tab of the charts group, so this call is normally the one that splits
    // the charts row off the tables group. `runCompute` calls it immediately
    // before `ensureProfileChart` for exactly that reason — tab order is dock
    // order, and the user-directed order is Mean±SD → Combined → Facet.
    const node = dockChartTab(view, attachBannerHost(chart.root), 'Mean ± SD');
    tempBag[MEAN_SD_CHART_TAG] = chart;
    tempBag[MEAN_SD_CHART_NODE_TAG] = node;
  }

  const nCol = profile.col('N')!;
  const tCol = profile.col('Time')!;
  const gCol = profile.col('Group')!;
  const ns: number[] = [];
  const groups = new Set<string>();
  const fragileCells: FragileCell[] = [];
  let minCellN = Number.POSITIVE_INFINITY;
  for (let i = 0; i < profile.rowCount; i++) {
    const nv = nCol.get(i) as number;
    ns.push(nv);
    if (nv < minCellN) minCellN = nv;
    const g = String(gCol.get(i));
    groups.add(g);
    // F4: capture EVERY fragile (group × time) cell (N < FRAGILE_N), not just the
    // global-median outcome — a heterogeneous composite hides a fragile cell behind a
    // healthy median otherwise (thesis §4).
    if (nv < FRAGILE_N) fragileCells.push({group: g, time: tCol.get(i) as number, n: nv});
  }
  // Median N per (group × time) cell (lower-median for even N) — retained for the
  // observed-time serial-drift signal; F4 uses minCellN + fragileCells for the banner.
  ns.sort((a, b) => a - b);
  const medianN = ns.length === 0 ? 0 : ns[Math.floor((ns.length - 1) / 2)];
  if (!Number.isFinite(minCellN)) minCellN = 0;
  return {created: existing == null, medianN, groupCount: groups.size, minCellN, fragileCells};
}

/** Close and forget the "Mean ± SD" tab. No-op when none is docked. */
function closeMeanSdChart(view: DG.TableView): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  const node = tempBag[MEAN_SD_CHART_NODE_TAG] as DG.DockNode | null;
  if (node != null) {
    try {
      view.dockManager.close(node);
    } catch {
      // Node already detached — ignore.
    }
  }
  const chart = tempBag[MEAN_SD_CHART_TAG] as DG.Viewer | null;
  if (chart != null) {
    try {
      chart.close();
    } catch {
      // Already closed — ignore.
    }
  }
  tempBag[MEAN_SD_CHART_TAG] = undefined;
  tempBag[MEAN_SD_CHART_NODE_TAG] = undefined;
  // Clear the render context + cached profile too (completeness, matching
  // closeSparseNca) so a stale post-close state can't be mistaken for "computed".
  tempBag[MEAN_SD_CTX_TAG] = undefined;
  tempBag[MEAN_SD_PROFILE_TAG] = undefined;
}

// =====================================================================
// tk-aggregates AD-5 / AD-6 — toxicokinetic aggregate panels.
// Faithful copies of the bioavailability dock lifecycle above: detect
// applicability on a successful Compute, compute the adapter primitives over the
// post-appendGroupingColumns frame, dock FILL-tabs next to the results node.
// =====================================================================

/** Read a study-level constant (NOAEL / human dose / human BW) from a carried
 *  results column — the first non-null row. Null when absent / all-null. These
 *  ride the source CSV as constant-per-subject metadata (06 schema:
 *  `NOAEL_mg_kg` / `HumanDose_mg` / `HumanBodyWeight_kg`); a dedicated mapping UI
 *  for them (+ an empirical human-AUC channel) is the tracked DG-1 / RG-3 follow-up. */
function readStudyConstant(results: DG.DataFrame, name: string): number | null {
  const col = results.col(name);
  if (col === null) return null;
  for (let i = 0; i < results.rowCount; i++)
    if (!col.isNone(i)) return Number(col.get(i));
  return null;
}

/** Max finite value of a source column (τ default = max NominalTime — AD-1 (b)). */
function maxFiniteInColumn(df: DG.DataFrame, name: string): number {
  const col = df.col(name);
  if (col === null) return NaN;
  let mx = NaN;
  const n = df.rowCount;
  for (let i = 0; i < n; i++) {
    const v = numAt(col, i);
    if (Number.isFinite(v) && (!Number.isFinite(mx) || v > mx)) mx = v;
  }
  return mx;
}

/** Append a per-subject `TLAST` (last sampling time) column to the results frame
 *  for the `tlast == τ` guard (D4). Row i ↔ `profiles[i]` (appendGroupingColumns
 *  alignment). No-op if already present or the time column is missing. */
function appendTlastColumn(
  results: DG.DataFrame, df: DG.DataFrame,
  profiles: ReadonlyArray<ProfileIndex>, timeName: string,
): void {
  if (results.col('TLAST') !== null) return;
  const timeCol = df.col(timeName);
  if (timeCol === null || results.rowCount !== profiles.length) return;
  const out = new Float64Array(profiles.length);
  for (let i = 0; i < profiles.length; i++) {
    let mx = NaN;
    for (const r of profiles[i].observationRows) {
      const t = numAt(timeCol, r);
      if (Number.isFinite(t) && (!Number.isFinite(mx) || t > mx)) mx = t;
    }
    out[i] = mx;
  }
  // `f64Column`, NOT the bare constructor: `TLAST` is added to the RESULTS frame
  // (and read back as a result column by `computeAccumulationRatios` via
  // `options.tlastCol`), so a profile with no finite sampling time must land as
  // the platform null. A raw NaN here would be invisible to `isNone`, counted as
  // present by every DG-native consumer, and would poison this column's
  // statistics. Found by the Arm A inter-rater pass, 2026-08-03 — this writer sits
  // in `views/`, which F1's scope table excluded wholesale on the assumption that
  // `views/` builds only input frames; that assumption was wrong for this one
  // function.
  results.columns.add(f64Column('TLAST', out));
}

/** Generic swap-or-dock for a TK aggregate grid (FILL-tab on the results node).
 *  Mirrors {@link swapOrDockBioavailability}. */
function swapOrDockToxGrid(
  view: DG.TableView, table: DG.DataFrame, gridTag: string, nodeTag: string, title: string,
  anchorNode: DG.DockNode | null,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  // FIRST-COMPUTE dock fix (root cause). `view.temp[key]` returns **null** for a
  // never-set key, NOT undefined (DG bag convention — see `ensureProfileChart`'s
  // `!= null` check). The old guard `existingGrid !== undefined` was true for that
  // null, so the first Compute wrongly took the swap branch and ran
  // `null.dataFrame = table` → "Cannot set properties of null (setting
  // 'dataFrame')". It was caught and the panel only appeared on a re-Compute
  // (when the key had been set to undefined by closeToxGrid). Guard with `!= null`.
  const existingGrid = tempBag[gridTag] as DG.Grid | null;
  if (existingGrid != null) {
    existingGrid.dataFrame = table;
    configureTkGrid(existingGrid); // fresh frame → re-apply friendly headers/format/tooltips
    return;
  }
  const grid = DGNS.Viewer.grid(table);
  configureTkGrid(grid);
  // Anchor to the results node returned by swapOrDockResults (explicit, no
  // findNode-after-topology-change); fall back to temp / findNode defensively.
  const anchor = resultsAnchor(view, anchorNode);
  tempBag[nodeTag] = view.dockManager.dock(attachBannerHost(grid.root), DGNS.DOCK_TYPE.FILL, anchor, title);
  tempBag[gridTag] = grid;
}

/** Close and forget a single TK aggregate grid tab. No-op when none is docked. */
function closeToxGrid(view: DG.TableView, gridTag: string, nodeTag: string): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  const node = tempBag[nodeTag] as DG.DockNode | null;
  if (node != null) {
    try {
      view.dockManager.close(node);
    } catch {
      // Node already gone — ignore.
    }
  }
  tempBag[gridTag] = undefined;
  tempBag[nodeTag] = undefined;
}

/** Close all TK aggregate tabs + dispose the overlay-selection bridge. */
function closeToxicokinetics(view: DG.TableView): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  closeToxGrid(view, ACC_GRID_TAG, ACC_NODE_TAG);
  closeToxGrid(view, MARGIN_GRID_TAG, MARGIN_NODE_TAG);
  closeToxGrid(view, SEXFOLD_GRID_TAG, SEXFOLD_NODE_TAG);
  const sub = tempBag[TK_OVERLAY_SUB_TAG] as Unsubscribable | null;
  if (sub != null) sub.unsubscribe();
  tempBag[TK_OVERLAY_SUB_TAG] = undefined;
}

/** Compute and dock the TK aggregates when the results frame is repeat-dose
 *  (a carried StudyDay column with ≥ 2 distinct days) + Dose, else close them.
 *  Sex is mandatory (sex-stratified accumulation / fold); a repeat-dose frame
 *  with Dose but no Sex degrades with an explanatory toast (AC-14), mirroring
 *  {@link applyBioavailability}. */
function applyToxicokinetics(
  view: DG.TableView,
  results: DG.DataFrame,
  carried: ReadonlyArray<string>,
  sexName: string | null,
  doseColName: string | null,
  df: DG.DataFrame,
  profiles: ReadonlyArray<ProfileIndex>,
  subjectName: string,
  timeName: string | null,
  resultsNode: DG.DockNode | null,
): void {
  const studyDayName = carried.find((c) => /^study_?day$/i.test(c)) ?? null;
  // Not a repeat-dose frame → silently inapplicable (like F% on a single route).
  if (studyDayName === null || doseColName === null) {
    closeToxicokinetics(view);
    return;
  }
  const dayCol = results.col(studyDayName);
  if (dayCol === null) {closeToxicokinetics(view); return;}
  const days = new Set<string>();
  for (let i = 0; i < results.rowCount; i++) if (!dayCol.isNone(i)) days.add(dayCol.getString(i));
  if (days.size < 2) {closeToxicokinetics(view); return;} // single-day study: no accumulation

  // Repeat-dose frame confirmed; Sex is mandatory for the aggregates (AC-14).
  // `sexName` is the designated (and carried) sex column, resolved by the caller.
  if (sexName === null) {
    // TK tabs are not docked here (no Sex → no aggregates) → banner over the
    // results grid (B3). Anchored to 'accumulation', the primary TK panel.
    collectNotice(view, {
      key: 'TK_NEEDS_SEX', anchor: 'accumulation', severity: 'info',
      // Action path is the Subject Sex picker (A4) — point at it, or at adding a
      // sex axis if the data has none. See F_PCT_NEEDS_SEX above.
      text:
        'TK aggregates (accumulation ratio, sex comparison) need a per-subject Sex ' +
        'column to stratify. Pick the Sex column in the Subject Sex section, or add ' +
        'one if your data has none.',
    });
    closeToxicokinetics(view);
    return;
  }

  const tau = timeName !== null ? maxFiniteInColumn(df, timeName) : NaN;
  if (timeName !== null) appendTlastColumn(results, df, profiles, timeName);
  const warnMsgs: string[] = [];

  // --- AD-1 / AD-2: accumulation ratio + theoretical ceiling.
  let acc: AccumulationResult | null = null;
  try {
    acc = computeAccumulationRatios(results, {
      tau, doseCol: doseColName, studyDayCol: studyDayName, stratifyCols: [sexName],
      tlastCol: results.col('TLAST') !== null ? 'TLAST' : undefined,
    });
    for (const w of acc.warnings) warnMsgs.push(w.message);
    if (acc.table.rowCount > 0)
      swapOrDockToxGrid(view, acc.table, ACC_GRID_TAG, ACC_NODE_TAG, 'Accumulation', resultsNode);
    else closeToxGrid(view, ACC_GRID_TAG, ACC_NODE_TAG);
  } catch {
    closeToxGrid(view, ACC_GRID_TAG, ACC_NODE_TAG);
  }

  // --- AD-3: NOAEL exposure/safety margin (when NOAEL + human dose are carried).
  const noael = readStudyConstant(results, 'NOAEL_mg_kg');
  const humanDoseMg = readStudyConstant(results, 'HumanDose_mg');
  const humanBw = readStudyConstant(results, 'HumanBodyWeight_kg');
  if (noael !== null && humanDoseMg !== null && humanBw !== null) {
    try {
      const margin = computeExposureMargin(results, {
        noaelMgKg: noael, humanDoseMg, humanBodyWeightKg: humanBw,
        doseCol: doseColName, studyDayCol: studyDayName, stratifyCols: [sexName],
      });
      for (const w of margin.warnings) warnMsgs.push(w.message);
      if (margin.table !== null && margin.table.rowCount > 0)
        swapOrDockToxGrid(view, margin.table, MARGIN_GRID_TAG, MARGIN_NODE_TAG, 'Exposure Margin', resultsNode);
      else closeToxGrid(view, MARGIN_GRID_TAG, MARGIN_NODE_TAG);
    } catch {
      closeToxGrid(view, MARGIN_GRID_TAG, MARGIN_NODE_TAG);
    }
  } else
    closeToxGrid(view, MARGIN_GRID_TAG, MARGIN_NODE_TAG);


  // --- AD-4: sex fold-difference over a (Dose × Sex) aggregate.
  try {
    const agg = aggregateNcaResults(results, [doseColName, sexName]);
    const fold = computeSexFold(agg, {sexCol: sexName});
    if (fold.table.rowCount > 0)
      swapOrDockToxGrid(view, fold.table, SEXFOLD_GRID_TAG, SEXFOLD_NODE_TAG, 'Sex Comparison', resultsNode);
    else closeToxGrid(view, SEXFOLD_GRID_TAG, SEXFOLD_NODE_TAG);
  } catch {
    closeToxGrid(view, SEXFOLD_GRID_TAG, SEXFOLD_NODE_TAG);
  }

  // Accumulation + exposure-margin warnings are merged into one notice (matching
  // the prior joined balloon — channel change only, not content) anchored to the
  // primary TK panel (Accumulation); falls back to the results grid if that grid
  // produced no rows (B3).
  if (warnMsgs.length > 0) {
    collectNotice(view, {
      key: 'TK_WARNINGS', anchor: 'accumulation', severity: 'warning',
      text: warnMsgs.join(' '),
    });
  }

  // --- AD-5: selecting an Accumulation (Dose×Sex) cell overlays that cell's
  // Day 1 vs Day N profiles in the Combined chart (reuses the chart + semi-log
  // toggle verbatim; only the cell → two-series binding is new).
  if (acc !== null)
    wireAccumulationToProfileOverlay(view, acc.table, results, sexName, doseColName, subjectName);
}

/** Bridge the Accumulation grid's current (Dose×Sex) cell → the Combined
 *  profile-overlay filter: the chart shows only that cell's pseudo-subjects
 *  (its Day 1 + Day N profiles). Mirrors
 *  {@link wireResultsSelectionToProfileFilter}; re-wired on every Compute. */
function wireAccumulationToProfileOverlay(
  view: DG.TableView,
  accTable: DG.DataFrame,
  results: DG.DataFrame,
  sexName: string,
  doseColName: string,
  subjectName: string,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  const prev = tempBag[TK_OVERLAY_SUB_TAG] as Unsubscribable | null;
  if (prev != null) prev.unsubscribe();

  const resSubjectCol = results.col('SUBJECT');
  const resSexCol = results.col(sexName);
  const resDoseCol = results.col(doseColName);
  const accSexCol = accTable.col(sexName);
  const accDoseCol = accTable.col('DOSE');
  if (resSubjectCol === null || resSexCol === null || resDoseCol === null ||
      accSexCol === null || accDoseCol === null)
    return;

  const apply = (): void => {
    const chart = tempBag[PROFILE_CHART_TAG] as DG.Viewer | null;
    if (chart == null) return;
    const chartDf = chart.dataFrame;
    if (chartDf == null) return;
    const chartSubjectCol = chartDf.col(subjectName);
    if (chartSubjectCol === null) {chartDf.filter.setAll(true); return;}
    const idx = accTable.currentRowIdx;
    if (idx < 0) {chartDf.filter.setAll(true); return;}
    const sex = accSexCol.isNone(idx) ? '' : accSexCol.getString(idx);
    const dose = accDoseCol.get(idx) as number;
    // Pseudo-subjects of this (Dose×Sex) cell (its Day 1 + Day N profiles).
    const subjects = new Set<string>();
    for (let i = 0; i < results.rowCount; i++) {
      const s = resSexCol.isNone(i) ? '' : resSexCol.getString(i);
      const d = numAt(resDoseCol, i);
      if (s === sex && d === dose) subjects.add(String(resSubjectCol.get(i)));
    }
    chartDf.filter.init((row) => subjects.has(String(chartSubjectCol.get(row))));
  };

  tempBag[TK_OVERLAY_SUB_TAG] = accTable.onCurrentRowChanged.subscribe(apply);
}

/** Close and forget the aggregated grid tab (Group-by → "(none)"). No-op when
 *  none is docked. */
function closeAggregated(view: DG.TableView): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  const node = tempBag[AGG_NODE_TAG] as DG.DockNode | null;
  if (node != null) {
    try {
      node.container.destroy();
      view.dockManager.close(node);
    } catch {
      // Node already gone (e.g. user closed the tab manually) — ignore.
    }
  }
  tempBag[AGG_GRID_TAG] = undefined;
  tempBag[AGG_NODE_TAG] = undefined;
  // Drop the linked-selection subscription and clear any curve filter it left on the
  // Mean±SD chart — grouping is off, so the chart shows every curve again (VR-21 teardown).
  const sub = tempBag[BYGROUP_SELECTION_SUB_TAG] as Unsubscribable | null;
  if (sub != null) sub.unsubscribe();
  tempBag[BYGROUP_SELECTION_SUB_TAG] = undefined;
  tempBag[BYGROUP_APPLY_TAG] = undefined;
  const meanSdChart = tempBag[MEAN_SD_CHART_TAG] as DG.Viewer | null;
  if (meanSdChart?.dataFrame != null) meanSdChart.dataFrame.filter.setAll(true);
}

/** Build a CHART-ONLY DataFrame from the user-picked Subject / Time /
 *  Concentration columns of the main DataFrame. Column names in the chart
 *  DataFrame match the source names verbatim so the chart's
 *  `xColumnName` / `yColumnNames` / `splitColumnName` keep the user-meaningful
 *  identifiers.
 *
 *  Transforms applied:
 *  - **Subject → STRING** (`String(value)` per row). DG line chart's split
 *    parameter requires a string column, but the adapter accepts
 *    string-or-int subjects (e.g. `theoph.csv` quoted "1"…"12" parses as
 *    INT). A STRING-typed copy satisfies the chart contract without
 *    constraining the adapter.
 *  - **Time → Float64 copy** with DG `FLOAT_NULL` translated to NaN.
 *  - **Conc → Float32 copy** of raw values, no substitution.
 *  - **log(Conc) → Float32**, `log10(v)` only when `v` is a valid
 *    positive concentration (`v !== FLOAT_NULL && v > 0`), otherwise
 *    `FLOAT_NULL`. The explicit FLOAT_NULL guard matters: DG's
 *    `FLOAT_NULL` is a tiny *positive* sentinel (~2.68e-34) — it passes
 *    a bare `v > 0` check and would produce a bogus `log10 ≈ -33.6`
 *    point on the chart for every null observation. Zeros, negatives,
 *    NaN, INT_NULL also correctly map to null via `v > 0`.
 *    Pre-derived so the chart can switch its y-binding to a log-scale
 *    series without recomputing per render.
 *
 *  The main DataFrame is NOT mutated — column-mapping pickers and the
 *  compute pipeline keep their unmodified view of the input. */
function buildChartDataFrame(
  mainDf: DG.DataFrame,
  subjectName: string,
  timeName: string,
  concName: string,
): DG.DataFrame {
  const subjectCol = mainDf.col(subjectName);
  const timeCol = mainDf.col(timeName);
  const concCol = mainDf.col(concName);
  if (subjectCol === null || timeCol === null || concCol === null)
    throw new Error('buildChartDataFrame: required column missing from main DataFrame');

  const n = mainDf.rowCount;

  const subjectValues = new Array<string>(n);
  for (let i = 0; i < n; i++) subjectValues[i] = String(subjectCol.get(i));
  const newSubjectCol = DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, subjectName, subjectValues);

  const timeRaw = timeCol.getRawData();
  const timeOut = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const v = timeRaw[i];
    timeOut[i] = (v === DGNS.FLOAT_NULL) ? NaN : v;
  }
  const newTimeCol = DGNS.Column.fromFloat64Array(timeName, timeOut);

  const concRaw = concCol.getRawData();
  const concOut = new Float32Array(n);
  const logOut = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    const v = concRaw[i];
    // Re-null guard (AC-F1) for any null-conc path (BLQ-as-null, sparse predose).
    // In Event-driven mode the EVID=1 dose rows (null DV) are already dropped
    // upstream, so this is residual.
    //
    // The `Number.isNaN` half is the load-bearing one: `isNone()` returns FALSE for
    // a raw NaN, so an unguarded NaN renders as the literal text "NaN" and collapses
    // the column's stats. The `=== FLOAT_NULL` half is belt-and-braces —
    // `Math.fround(FLOAT_NULL) === FLOAT_NULL`, so narrowing to Float32 preserves the
    // sentinel and its `isNone`. CORRECTION: this comment previously said the guard
    // stopped FLOAT_NULL from "copying a near-zero value"; that mechanism was
    // measured false on 2026-07-30 (see buildFacetDataFrames for the full note and
    // the probe reference). Do not re-derive the narrowing story.
    concOut[i] = toPlatformNull(v);
    logOut[i] = log10OrNull(v);
  }
  const newConcCol = DGNS.Column.fromFloat32Array(concName, concOut);
  const newLogCol = DGNS.Column.fromFloat32Array(`log(${concName})`, logOut);

  return DGNS.DataFrame.fromColumns([newSubjectCol, newTimeCol, newConcCol, newLogCol]);
}

/** True when Subject / Time / Concentration map to three DISTINCT source
 *  columns. The chart builders assemble a DataFrame whose column names mirror
 *  the source verbatim, so a name collision (e.g. Time mapped to the same
 *  column as Concentration) makes `DG.DataFrame.fromColumns` throw
 *  "Column named 'X' already exists". The Compute dup-guard
 *  (`duplicateColumnMappings`) blocks a RUN with colliding columns, but these
 *  builders take raw `string` arguments and do not structurally require that
 *  gate to have run first — so they self-guard at their own boundary as
 *  defense-in-depth. On a collision the builder is a no-op, leaving any existing
 *  chart intact rather than crashing. */
function chartColumnsDistinct(
  subjectName: string,
  timeName: string,
  concName: string,
): boolean {
  return subjectName !== timeName &&
    subjectName !== concName &&
    timeName !== concName;
}

/** Create or update the concentration-time profile chart for the current
 *  mapping selection. Lifecycle is driven by the three required mapping
 *  fields (Subject / Time / Concentration):
 *
 *  - All three present, chart does not exist, `createIfMissing` = true
 *    → create the chart (right-split next to the main+results column).
 *  - All three present, chart exists → update its column bindings via
 *    `setOptions`.
 *  - All three present, chart does not exist, `createIfMissing` = false
 *    → no-op (caller is the mapping `onChanged` subscription, which only
 *    refreshes an already-created chart).
 *  - Any of the three null → close the chart if it exists.
 *
 *  **Lifecycle ordering (important).** The first chart creation must happen
 *  AFTER the results grid is docked, not at mount time:
 *    1. User clicks Compute → `swapOrDockResults` FILLs results onto `view.grid`,
 *       forming the tables group.
 *    2. `refreshMeanSdChart(true)` splits the charts group off the top of it.
 *    3. THEN this function FILLs `Combined` in as the charts group's second tab.
 *  Subsequent mapping changes call `ensureProfileChart(..., false)` which
 *  only updates the already-docked chart.
 *
 *  The historical footgun this ordering also avoids: a `FILL`-dock next to
 *  `view.grid` AFTER the results grid had split the area raised inside DG's dock
 *  manager, because the main-grid leaf was no longer the right anchor. The
 *  current layout FILLs onto that leaf as its very first operation, before any
 *  split re-parents it, so the hazard is now structurally out of reach rather
 *  than merely avoided by convention.
 *
 *  Chart binds to derived columns (string subject, positive-only conc)
 *  rather than the user-picked columns directly. The derived columns are
 *  hidden in the main data grid so the user sees only the original CSV
 *  columns. See {@link ensureStringSubject} and {@link ensurePositiveConc}.
 *
 *  Y axis binds to the pre-derived `log(<conc>)` column from
 *  {@link buildChartDataFrame} — PK convention reads terminal-phase off a
 *  log Y, and pre-deriving the values (`v > 0 ? log10(v) : FLOAT_NULL`)
 *  side-steps the renderer ever seeing log(0) / log(negative).
 *  `yAxisType` stays `'linear'` because the *values* are already
 *  log-transformed; the facet grid uses the same column for the same
 *  reason. */
function ensureProfileChart(
  view: DG.TableView,
  df: DG.DataFrame,
  subjectName: string | null,
  timeName: string | null,
  concName: string | null,
  createIfMissing: boolean,
  useLogY: boolean,
): void {
  // `view.temp` is a DG-side bag, NOT a plain JS object. Missing keys come
  // back as `null` (not `undefined`), so the existence check uses `== null`
  // to catch both. The annotation reads `| null` for the same reason.
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  const existing = tempBag[PROFILE_CHART_TAG] as DG.Viewer | null;

  if (subjectName === null || timeName === null || concName === null) {
    if (existing != null) {
      existing.close();
      tempBag[PROFILE_CHART_TAG] = null;
    }
    return;
  }

  // Skip the rebuild on a Subject/Time/Concentration column collision — see
  // `chartColumnsDistinct`. Leaves any existing chart intact rather than
  // crashing in `buildChartDataFrame`.
  if (!chartColumnsDistinct(subjectName, timeName, concName)) return;

  if (existing == null && !createIfMissing) return;

  // Build a fresh chart-only DataFrame from the current mapping. Column
  // names match the source verbatim so chart options stay user-meaningful.
  // The main DataFrame is NEVER mutated — mapping inputs keep their original
  // column listing, no `__chart_*` ghosts in the picker dropdowns.
  const chartDf = buildChartDataFrame(df, subjectName, timeName, concName);
  const yColName = useLogY ? `log(${concName})` : concName;
  // Record the conc column this render is bound to so the semi-log toggle
  // (applyYAxisModeToCharts) re-binds against the rendered column, not a
  // since-edited mapping input. Written on both the create and re-bind paths.
  tempBag[PROFILE_CONC_NAME_TAG] = concName;

  const options: Partial<DG.ILineChartSettings> = {
    xColumnName: timeName,
    yColumnNames: [yColName],
    splitColumnName: subjectName,
    yAxisType: 'linear' as const,
    showMarkers: 'Always' as const,
    // Hide the split picker — subject is fixed for this chart and the
    // picker would let the user accidentally re-bind the split column.
    showSplitSelector: false,
    // Domain-meaningful header pinned to the top of the chart. Description
    // is a viewer-level title with markup support; visibility forced so it
    // shows regardless of the viewer's auto-hide threshold.
    description: useLogY ? 'Log concentration–time profile' : 'Concentration–time profile',
    descriptionPosition: 'Top' as const,
    descriptionVisibilityMode: 'Always' as const,
  };
  tempBag[PROFILE_CHART_OPTS_TAG] = options;

  if (existing != null) {
    // Re-bind the existing viewer to the new chart DataFrame, then apply
    // the (possibly changed) column options.
    existing.dataFrame = chartDf;
    existing.setOptions(options);
    return;
  }

  const chart = DGNS.Viewer.lineChart(chartDf, options);
  // Second tab of the charts group — `refreshMeanSdChart` runs just above this in
  // `runCompute` and has already created the group. (Earlier FILL attempts
  // crashed when the chart bound to the same DataFrame as the main grid — the
  // chart uses its own chart-only DataFrame, removing that aliasing.)
  dockChartTab(view, chart.root, 'Combined');
  tempBag[PROFILE_CHART_TAG] = chart;
}

/** Bridge results-grid selection → profile-chart filter.
 *
 *  When the user selects one or more rows in the results grid, the profile
 *  chart shows only the curves for the selected subjects (rows in the chart's
 *  underlying DataFrame whose subject column matches a selected `SUBJECT`).
 *  Empty selection → all curves visible.
 *
 *  Lifecycle: re-wired on every Compute. Each Compute installs a brand-new
 *  `results` DataFrame; the prior subscription points at a stale source and
 *  must be disposed before subscribing to the new one. Empty initial
 *  selection on the new `results` is fine — `apply()` runs once on wire-up
 *  and resets the chart filter to "all curves" via the trueCount=0 branch.
 *
 *  Selection lives on `results`; the filter is applied to whichever DataFrame
 *  the chart is currently bound to — read lazily inside `apply()` so a
 *  subsequent mapping change (which rebuilds chartDf in
 *  {@link ensureProfileChart}) does not strand the handler on an orphaned
 *  DataFrame.
 *
 *  Edge case: if mapping changes between Computes and chartDf is rebuilt with
 *  a different Subject column name, the captured `subjectName` no longer
 *  resolves in the new chartDf. We treat that as "no filter to apply" and
 *  reset to setAll(true) — matches the product decision to NOT preserve the
 *  filter across mapping changes (no curves staying hidden by stale selection
 *  state). */
function wireResultsSelectionToProfileFilter(
  view: DG.TableView,
  results: DG.DataFrame,
  subjectName: string,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;

  const prevSub = tempBag[SELECTION_SUB_TAG] as Unsubscribable | null;
  if (prevSub != null) prevSub.unsubscribe();

  const resultsSubjectCol = results.col('SUBJECT');
  if (resultsSubjectCol === null) return;

  const apply = (): void => {
    const chart = tempBag[PROFILE_CHART_TAG] as DG.Viewer | null;
    if (chart == null) return;
    const chartDf = chart.dataFrame;
    if (chartDf == null) return;
    const chartSubjectCol = chartDf.col(subjectName);
    if (chartSubjectCol === null) {
      chartDf.filter.setAll(true);
      return;
    }
    const selection = results.selection;
    if (selection.trueCount === 0) {
      chartDf.filter.setAll(true);
      return;
    }
    const selected = new Set<string>();
    const indexes = selection.getSelectedIndexes();
    for (let i = 0; i < indexes.length; i++)
      selected.add(String(resultsSubjectCol.get(indexes[i])));
    chartDf.filter.init((row) => selected.has(String(chartSubjectCol.get(row))));
  };

  tempBag[SELECTION_SUB_TAG] = results.selection.onChanged.subscribe(apply);
  apply();
}

/** Group rows of `mainDf` by `subjectName` and emit one chart-shaped
 *  DataFrame per subject. Each per-subject DataFrame has 3 columns:
 *  `<timeName>` (Float64, FLOAT_NULL→NaN), `<concName>` (Float32 raw copy),
 *  and `log(<concName>)` (Float32, `v > 0 ? log10(v) : FLOAT_NULL`). Mirrors
 *  the transform rules of {@link buildChartDataFrame}.
 *
 *  Subjects appear in **first-occurrence order** in the result array. That
 *  order is what the facet-grid color cycle keys off — keeping it stable
 *  across Computes (provided the source DataFrame is stable) means a given
 *  subject keeps its color on every re-render. */
export function buildFacetDataFrames(
  mainDf: DG.DataFrame,
  subjectName: string,
  timeName: string,
  concName: string,
  subjectLz?: ReadonlyMap<string, FacetLambdaZParams>,
  fitRowsBySubject?: ReadonlyMap<string, ReadonlySet<number>>,
  lzActionOf?: (rowIdx: number) => LzAction,
): ReadonlyArray<{subject: string, df: DG.DataFrame, lz?: FacetLambdaZParams | null}> {
  const subjectCol = mainDf.col(subjectName);
  const timeCol = mainDf.col(timeName);
  const concCol = mainDf.col(concName);
  if (subjectCol === null || timeCol === null || concCol === null)
    throw new Error('buildFacetDataFrames: required column missing from main DataFrame');
  // Interactive λz override state (VR-11, F6): drives the marker so a deliberately
  // excluded point is visible even when it still sits inside the fitted window. Read from the
  // controller's keyed store (nca-workflow-stages AD-2) — the SAME source the click-toggle reads.

  const n = mainDf.rowCount;

  const subjectKeys: string[] = [];
  const rowsBySubject = new Map<string, number[]>();
  for (let i = 0; i < n; i++) {
    const key = String(subjectCol.get(i));
    let rows = rowsBySubject.get(key);
    if (rows === undefined) {
      rows = [];
      rowsBySubject.set(key, rows);
      subjectKeys.push(key);
    }
    rows.push(i);
  }

  const timeRaw = timeCol.getRawData();
  const concRaw = concCol.getRawData();
  const logColName = `log(${concName})`;

  return subjectKeys.map((subject) => {
    const rows = rowsBySubject.get(subject)!;
    const lz = subjectLz?.get(subject) ?? null;
    // F6 (VR-11) — REAL λz fit membership, the honest marker source: the source rows
    // the reported fit actually used (translated from `pointsUsed` by the controller).
    // `null` ⇒ no controller context (pre-Compute / non-interactive) → fall back to the
    // time-window reconstruction below. When present it WINS: the window rule mislabels
    // a point that sits inside [tFirst,tLast] but is absent from a gappy manual set (the
    // bug the user hit — a reset point stayed 'included' though it had left the fit).
    const fitRows = fitRowsBySubject?.get(subject) ?? null;
    // M1 point-role window (fallback only): a point is 'included' iff its time falls in
    // the reported terminal window [tFirst, tLast] AND it is measurable. NaN bounds (a
    // subject with no fit) make every comparison false → all points 'excluded'.
    const winLo = lz !== null && Number.isFinite(lz.tFirst) ? lz.tFirst : NaN;
    const winHi = lz !== null && Number.isFinite(lz.tLast) ? lz.tLast : NaN;
    const m = rows.length;
    const timeOut = new Float64Array(m);
    const concOut = new Float32Array(m);
    const logOut = new Float32Array(m);
    const roleOut = new Array<string>(m);
    // Source-df row index per point — lets the interactive λz editor map a clicked
    // facet point back to its observations-frame row (VR-11, C0). Never displayed.
    const srcRowOut = new Int32Array(m);
    for (let j = 0; j < m; j++) {
      const r = rows[j];
      srcRowOut[j] = r;
      const t = timeRaw[r];
      timeOut[j] = (t === DGNS.FLOAT_NULL) ? NaN : t;
      const c = concRaw[r];
      // Re-null guard on BOTH columns, mirroring {@link buildChartDataFrame}.
      // (Symbol reference, not a line pin — a line number here goes stale on the
      // sibling's next edit, and this comment's whole point is that its citations
      // are checkable.)
      // The linear column was a bare copy until 2026-07-30 (user QA): `8849b44`
      // guarded the sibling builder and never touched this function, while the
      // comment here already claimed to mirror it.
      //
      // WHAT EACH BRANCH IS FOR — measured on the live platform:
      //
      //  - `Number.isNaN(c)` is the LOAD-BEARING half. `col.isNone()` returns FALSE
      //    for a raw IEEE NaN — DG does not consider it missing. An unguarded NaN
      //    therefore reaches the grid as the literal text "NaN", and a single one
      //    collapses min/max/avg/sum/stdev on the whole column to NaN.
      //
      //  - `c === FLOAT_NULL` is belt-and-braces, NOT a fix. FLOAT_NULL
      //    (2.6789344063684636e-34) is exactly representable in Float32 —
      //    `Math.fround(FLOAT_NULL) === FLOAT_NULL` — and the narrowed value still
      //    reads as `isNone`. Narrowing does NOT drop the sentinel. Kept only so
      //    this line stays character-identical to the sibling; do not cite it as
      //    the reason the guard exists.
      //
      // CORRECTION TO NOTE-12 LORE: that note attributed a spurious ~0 chart point
      // to FLOAT_NULL "not surviving" the Float32 narrowing. The measurement above
      // contradicts that mechanism. `8849b44` was a two-part fix — structural (drop
      // EVID=1 dose rows upstream) plus this residual guard — and the evidence says
      // the STRUCTURAL half did the work. Do not re-derive the narrowing story.
      //
      // The LOG line's guard is separately load-bearing for a different reason:
      // FLOAT_NULL passes `> 0`, so without the explicit check it would store a
      // REAL log10(2.68e-34) ≈ -33.6 instead of a gap.
      concOut[j] = toPlatformNull(c);
      const logVal = log10OrNull(c);
      logOut[j] = logVal;
      // Measurable ⇔ the log guard kept the value (positive, finite, not the
      // FLOAT_NULL sentinel) — the SAME eligibility sci-comp's λz fit uses, so an
      // embedded BLQ / ≤0 point inside the window is NOT mislabelled 'included'
      // (thesis #4). Window test is inclusive of the reported tFirst/tLast bounds.
      const measurable = logVal !== DGNS.FLOAT_NULL;
      const inWindow = timeOut[j] >= winLo && timeOut[j] <= winHi; // false when winLo/winHi NaN
      const act = lzActionOf !== undefined ? lzActionOf(r) : 'auto';
      if (fitRows !== null) {
        // Honest membership (F6): the reported fit's real rows drive the marker, so it
        // cannot disagree with the fit. In fit → filled (auto-in OR force-in). Not in
        // fit: a deliberate `force-out` renders as the OUTLIER (×); anything else
        // (auto-excluded, or a force-in the min-points floor dropped) renders hollow.
        roleOut[j] = fitRows.has(r) ? 'included' :
          act === 'force-out' ? 'force-out' : 'excluded';
      } else {
        // Fallback (no controller context): the time-window reconstruction. A
        // `force-out` still wins as OUTLIER; `force-in` renders included; an `auto`
        // point uses the window rule.
        roleOut[j] = act === 'force-out' ? 'force-out' :
          act === 'force-in' ? 'included' :
            (measurable && inWindow) ? 'included' : 'excluded';
      }
    }
    const df = DGNS.DataFrame.fromColumns([
      DGNS.Column.fromFloat64Array(timeName, timeOut),
      DGNS.Column.fromFloat32Array(concName, concOut),
      DGNS.Column.fromFloat32Array(logColName, logOut),
      DGNS.Column.fromStrings(FACET_ROLE_COL, roleOut),
      DGNS.Column.fromInt32Array(FACET_SRCROW_COL, srcRowOut),
    ]);
    return {subject, df, lz};
  });
}

/**
 * Number of **columns** in the split facet grid for `n` total plots — how many
 * plots stand side by side (= the row width). {@link buildFacetGrid} allocates
 * exactly this many columns and deals plots into them round-robin, so the column
 * count decides the row shapes: e.g. N=5 over 3 columns → rows of 3 then 2.
 *
 * Governs the SPLIT layout only, i.e. `n <= FACET_SCROLL_ABOVE` (10). Past that the
 * grid switches to the fixed {@link FACET_VIEWPORT_COLS}-wide scrolling layout and
 * this table is not consulted — which is why it only covers 1..10.
 *
 * **User-directed table** (2026, this cycle) — the requested per-N shapes:
 *   1→1, 2→2, 3→3 (single row); 4→2 (2×2); 5→3 (3+2); 6→3 (3+3); 7→4 (4+3);
 *   8→4 (4+4); 9→5 (5+4); 10→5 (5+5).
 */
export function getMaxFacetCols(n: number): number {
  switch (n) {
  case 1:
    return 1;
  case 2:
    return 2;
  case 3:
    return 3;
  case 4:
    return 2;
  case 5:
  case 6:
    return 3;
  case 7:
  case 8:
    return 4;
  case 9:
  case 10:
    return 5;
  default:
    // N>10 uses the scrolling layout (this is never consulted there); safe fallback.
    return FACET_VIEWPORT_COLS;
  }
}

/** Variant-A per-tile colour for the scatterPlot facet: paint the tile's
 *  verdict hue ({@link buildFacetVerdictColors}, or `STATUS_NEUTRAL` when absent)
 *  onto a constant category column the plot colours by.
 *  scatterPlot has no `lineColor`/`markerColor`, so the hue must live ON the
 *  DataFrame (in {@link FACET_COLOR_COL}'s categorical colour tags) rather than
 *  as a viewer option — which is why BOTH the create path ({@link buildFacetGrid})
 *  and the in-place `dataFrame`-swap reuse path ({@link ensureFacetGrid}) call this
 *  after (re)building a frame: a fresh frame carries no colour tags, and
 *  `colorColumnName`/`linesByColumnName` on the viewer would resolve to a colourless
 *  column without it. Adds the column on first sight (a facet frame otherwise holds
 *  only Time/Conc/log(Conc)); idempotent thereafter. */
function paintFacetSeries(df: DG.DataFrame, subject: string, color: number): void {
  // Category token `subj '<subject>'` (user-directed 2026-08-26): an explicit,
  // unambiguously-string label. CONFIRMED root cause of the all-tiles-blue bug —
  // a bare numeric-looking subject ("1") made DG infer a NUMERIC column, so the
  // categorical colour-coding set below silently did not apply. The prefix keeps
  // the column categorical.
  const category = `subj '${subject}'`;
  let cat = df.col(FACET_COLOR_COL);
  if (cat === null) {
    const arr = new Array<string>(df.rowCount).fill(category);
    cat = df.columns.add(DGNS.Column.fromStrings(FACET_COLOR_COL, arr));
  }
  cat.colors.setCategorical({[category]: color});
}

/** Spike `facet-lambdaz-line-points` (M3): per-category marker SHAPE on the
 *  {@link FACET_ROLE_COL} column so 'included' vs 'excluded' λz points read
 *  differently. Mirrors {@link paintFacetSeries}: the coding lives on the COLUMN
 *  (a `.marker-coding` tag), so a fresh frame from {@link buildFacetDataFrames}
 *  carries none and BOTH the create and reuse paths must re-apply it after a
 *  frame swap. Shapes are THROWAWAY (final treatment is a Rule-3 design call). */
function applyFacetPointMarkers(df: DG.DataFrame): void {
  const col = df.col(FACET_ROLE_COL);
  if (col === null) return;
  // included = filled circle, auto-excluded = hollow circle (CIRCLE_BORDER),
  // deliberately force-out = OUTLIER (D9). Crosses stay reserved (user-directed).
  col.markers.assign('included', DGNS.MARKER_TYPE.CIRCLE);
  col.markers.assign('excluded', DGNS.MARKER_TYPE.CIRCLE_BORDER);
  col.markers.assign('force-out', DGNS.MARKER_TYPE.OUTLIER);
}

/** Spike `facet-lambdaz-line-points` (M4): draw the subject's λz terminal line on
 *  the tile via native Datagrok `formulaLines` — `ln(C)=intercept−λz·t`, confined
 *  to `[tFirst,tLast]`. Cleared+re-added (idempotent) so the create path, the
 *  reuse-path frame swap, and (indirectly) the semi-log toggle can all rely on it.
 *
 *  **The line content is inherently SEMI-LOG** — a log-linear fit is a straight line
 *  only in log space. It is ALWAYS written to the frame here (against the pre-derived
 *  `log(<conc>)` column); the LINEAR axis HIDES it via the viewer's
 *  `showDataframeFormulaLines: useLogY` gate rather than re-drawing it (an `exp(…)`
 *  linear form did not evaluate in DG formula-lines).
 *
 *  **Vertical anchor is DISPLAY-derived, not sci-comp's intercept (option A).** The
 *  fit's `intercept` is in the concentration units computeNca was run on, which can
 *  DIFFER from the facet's displayed conc column (e.g. computeNca on µg/mL vs a
 *  `Conc_ng_mL` display column — a ×1000 gap sinks the line ~log10(1000)=3 decades).
 *  So we keep sci-comp's SLOPE (λz, unit-independent, validated) but re-anchor the
 *  line vertically to the facet's OWN plotted points: `a = mean over the window's
 *  measurable points of [log10(C_display) + (λz/ln10)·t]`, giving
 *  `log10(C) = a − (λz/ln10)·t`. This overlays the terminal points in whatever units
 *  the tile shows. (The underlying display-vs-compute unit gap is a separate upstream
 *  concern — logged as a gap; here we only make the OVERLAY faithful.) The
 *  include/exclude MARKERS are scale-independent and stay on both axes. No line when
 *  the fit failed (NaN/≤0 λz) or the window holds no measurable point. Style/colour
 *  are THROWAWAY (Rule-3 design call). */
export function applyFacetLambdaZLine(
  df: DG.DataFrame,
  lz: FacetLambdaZParams | null,
  concName: string,
  timeName: string,
): void {
  const fl = df.meta.formulaLines;
  fl.clear();
  if (lz === null || !Number.isFinite(lz.lambdaZ) || lz.lambdaZ <= 0) return;
  const timeCol = df.col(timeName);
  const logCol = df.col(`log(${concName})`);
  if (timeCol === null || logCol === null) return;
  const lo = Number.isFinite(lz.tFirst) ? lz.tFirst : -Infinity;
  const hi = Number.isFinite(lz.tLast) ? lz.tLast : Infinity;
  const slopeLog = lz.lambdaZ / Math.LN10; // slope of log10(C) vs t (line descends by −slopeLog)
  // Anchor sci-comp's slope to the facet's own display-unit window points:
  // a = mean[ log10(C_display) + slopeLog·t ] over measurable points in [tFirst,tLast].
  const tRaw = timeCol.getRawData();
  const logRaw = logCol.getRawData();
  const n = df.rowCount;
  let sum = 0;
  let cnt = 0;
  for (let k = 0; k < n; k++) {
    const lv = logRaw[k];
    if (lv === DGNS.FLOAT_NULL || !Number.isFinite(lv)) continue; // BLQ / ≤0 / gap
    const t = tRaw[k];
    if (t === DGNS.FLOAT_NULL || !(t >= lo && t <= hi)) continue;
    sum += lv + slopeLog * t;
    cnt++;
  }
  if (cnt === 0) return; // no plotted point in the window to anchor to
  const a = sum / cnt;
  const formula = `\${log(${concName})} = ${a} - ${slopeLog} * \${${timeName}}`;
  // Black + POSITIVE zIndex so the line draws ON TOP of the coloured data line (chart
  // data sits at depth 0; a negative zIndex hid the λz line behind the points it
  // overlaps in the terminal window). `title: 'λz'` labels the line AND suppresses the
  // raw-formula fallback DG prints when a line has no title; `description` is
  // deliberately OMITTED — adding it triggered a DG on-hover null-error
  // (`html_utils.dart htmlMouseOverItemVisibility`) on every mousemove over the tile.
  // Colour/width/style are a Rule-3 design call.
  const item: DG.FormulaLine = {
    title: 'λz', formula, color: '#000000', width: 2, style: 'dashed',
    visible: true, opacity: 100, zIndex: 20,
  };
  if (Number.isFinite(lz.tFirst)) item.min = lz.tFirst;
  if (Number.isFinite(lz.tLast)) item.max = lz.tLast;
  fl.addLine(item);
}

/** Tag each tile so CSS drops the separator on the grid's OUTER right/bottom edge
 *  (`.nca-facet-last-col`/`.nca-facet-last-row` — a separator divides two tiles, so
 *  the perimeter carries none). Positions are row-major with `cols` columns — the
 *  same fill order for BOTH layouts (the split grid deals tiles round-robin into
 *  `cols` columns, which visits row 0 left-to-right, then row 1, …; the scrolling
 *  CSS grid fills row-major directly). A tile keeps its right border only if a tile
 *  sits to its right in the SAME row, and its bottom border only if a tile sits
 *  directly below — so a partial last row trims correctly (a tile with an empty cell
 *  beneath it is a bottom edge). Position-only, so it survives a same-set reorder
 *  (count and slots unchanged) with no re-tag. */
function tagFacetEdgeTiles(plots: ReadonlyArray<DG.Viewer>, cols: number): void {
  const n = plots.length;
  const lastRow = Math.floor((n - 1) / cols);
  for (let i = 0; i < n; i++) {
    // Keep the right separator for every column except the rightmost — not just where
    // a tile sits DIRECTLY to the right. In a partial last row a tile that is not in
    // the last column (e.g. col 3 of a 5-wide grid) still needs its right border so
    // the vertical separator runs down the full height instead of stopping at the row
    // above. Symmetric to `hasBottom` below.
    const hasRight = (i % cols) < (cols - 1);
    // Keep the bottom separator for EVERY row except the visually last one — not just
    // where a tile sits DIRECTLY below. In a partial last row the columns above the
    // empty slots (`i + cols >= n`) still need their bottom border so the horizontal
    // separator spans the full width instead of stopping short under the last column.
    const hasBottom = Math.floor(i / cols) < lastRow;
    if (!hasRight) plots[i].root.classList.add('nca-facet-last-col');
    if (!hasBottom) plots[i].root.classList.add('nca-facet-last-row');
  }
}

/** Build the per-subject facet grid as a `ui.splitH(ui.splitV(...))` of
 *  independent scatter plots — same layout idiom as DiffStudio's
 *  `getFacetPlot` (`app.ts:3370-3435`). Each plot binds its own per-subject
 *  DataFrame to `xColumnName=time` and `yColumnName=yColName`, connecting the
 *  points into a line via `linesOrderColumnName=time` (scatterPlot draws no line
 *  until an order column is set); colour = the subject's verdict hue
 *  ({@link buildFacetVerdictColors}, `STATUS_NEUTRAL` when absent), and
 *  carried on the frame by {@link paintFacetSeries} (scatterPlot has no
 *  `lineColor`/`markerColor`); the subject is the always-on plot TITLE and
 *  its λz acceptance verdict (`subjectDescriptions`) the always-on bottom caption.
 *  Every tile root carries `nca-facet-tile` for the dashed inter-tile separators.
 *
 *  Mirrors the Mean±SD scatterPlot (`ensureMeanSdChart`): same
 *  `linesOrderColumnName`/`linesWidth`/`colorColumnName` line-drawing recipe, and
 *  the same suppression of scatterPlot edge decorations (axis histograms,
 *  filtered-out ghosts, drop lines) that otherwise render as stray strips.
 *
 *  Axes are shown but their column-selectors are hidden (`showXSelector`,
 *  `showYSelector` both false) — each facet plot's X/Y binding is fixed at build
 *  time and the selectors would only invite accidental re-binding.
 *  `autoLayout: false` keeps the chart
 *  drawing area stable across the grid (DG's auto-layout heuristic
 *  reshuffles axes / paddings based on per-plot data extent, which made
 *  neighbouring tiles render at different effective sizes).
 *  Empty trailing cells in shorter columns get a placeholder `ui.div('')` so
 *  `ui.splitV` makes all columns equal-height. */
export function buildFacetGrid(
  perSubjectDfs: ReadonlyArray<{subject: string, df: DG.DataFrame, lz?: FacetLambdaZParams | null}>,
  timeName: string,
  yColName: string,
  subjectDescriptions?: ReadonlyMap<string, string>,
  concName?: string,
  useLogY?: boolean,
  subjectColors?: ReadonlyMap<string, number>,
): {root: HTMLDivElement, plots: DG.Viewer[]} {
  const n = perSubjectDfs.length;
  const scrolling = n > FACET_SCROLL_ABOVE;
  // Scrolling layout is a flat grid: one column bucket is unused, and plots are
  // appended in subject order rather than dealt round-robin, because a CSS grid
  // fills row-major and dealing would scramble the reading order.
  // Never allocate more side-by-side columns than there are tiles — otherwise the
  // split layout leaves an empty trailing column and the grid does not fill the pane
  // width (getMaxFacetCols returns 4 for N=1..3, so 3 subjects left a 4th blank slot).
  const maxInRow = scrolling ? 1 : Math.min(getMaxFacetCols(n), n);
  const cols: DG.Viewer[][] = new Array(maxInRow);
  for (let i = 0; i < maxInRow; i++) cols[i] = [];

  const plots: DG.Viewer[] = [];
  let idx = 0;
  for (let i = 0; i < n; i++) {
    const {subject, df} = perSubjectDfs[i];
    // Tile hue = ACCEPTANCE-verdict colour when supplied (3-way, matches the grid cell),
    // else STATUS_NEUTRAL — a no-verdict tile must NOT wear a colour that reads as a verdict.
    const color = subjectColors?.get(subject) ?? STATUS_NEUTRAL;
    // Variant-A colour: paint the tile hue onto the frame BEFORE the viewer reads
    // `colorColumnName` (scatterPlot has no lineColor/markerColor).
    paintFacetSeries(df, subject, color);
    // Spike `facet-lambdaz-line-points`: per-point marker shape (M3) + λz line
    // (M4) — both live on the frame, so applied here at build (and re-applied on
    // the reuse-path frame swap, exactly like paintFacetSeries).
    applyFacetPointMarkers(df);
    if (concName !== undefined)
      applyFacetLambdaZLine(df, perSubjectDfs[i].lz ?? null, concName, timeName);
    const plot = DGNS.Viewer.scatterPlot(df, {
      xColumnName: timeName,
      yColumnName: yColName,
      autoLayout: false,
      showXAxis: true,
      showYAxis: true,
      // Gridlines on both axes (was missing — the tiles read as floating points).
      showVerticalGridLines: true,
      showHorizontalGridLines: true,
      // Hide EVERY interactive selector chip (X/Y axis, colour, size/marker) — a
      // facet tile's bindings are fixed at build time and the chips only invite
      // accidental re-binding (user-directed 2026-08-26).
      showXSelector: false,
      showYSelector: false,
      showSizeSelector: false,
      yAxisType: 'linear' as const,
      markerType: DGNS.MARKER_TYPE.CIRCLE,
      // Connect the points into a profile: scatterPlot draws no line until an
      // order column is set, so order by time and group the single series by the
      // colour column (one category per tile → one coloured polyline).
      linesOrderColumnName: timeName,
      linesByColumnName: FACET_COLOR_COL,
      linesWidth: 2,
      colorColumnName: FACET_COLOR_COL,
      showColorSelector: false,
      // Spike `facet-lambdaz-line-points`: marker shape encodes λz include/exclude
      // (M3); native formula lines carry the λz terminal line (M4). Both persist
      // across the reuse-path setOptions merge — only the frame-level coding/line
      // is re-applied there. markersColumnName is bound only when the role column
      // is present (a real facet frame always has it; a degenerate/mirror frame
      // without it stays unbound rather than pointing markers at a missing column).
      ...(df.col(FACET_ROLE_COL) !== null ? {markersColumnName: FACET_ROLE_COL} : {}),
      // λz line is semi-log-only: shown when this tile is on the log axis, hidden
      // on linear (the line content is always present on the frame; this GATES it).
      showDataframeFormulaLines: useLogY === true,
      // Colour is a per-tile decoration (one category), not a scale the reader
      // decodes — the tile TITLE already names the subject. Hide the colour legend
      // the viewer adds once `colorColumnName` is set (user-directed 2026-08-26).
      legendVisibility: 'Never' as const,
      // Suppress scatterPlot edge decorations that render as stray strips/points
      // at the plot bottom/right (Mean±SD lesson — see ensureMeanSdChart).
      showXHistogram: false,
      showYHistogram: false,
      showFilteredOutPoints: false,
      showDropLines: false,
      // A facet click is an OVERRIDE action, not row navigation — DG's current-row
      // highlight (the black ring on the clicked point, and its axis-projected
      // crosshair strip) is unwanted noise that persists across the frame swap. The
      // override has its OWN feedback: the marker changes to filled / hollow / ×.
      // Hover feedback (showMouseOverPoint) stays on so the target point is legible.
      showCurrentPoint: false,
      // Subject + violated-criteria on ONE always-on line at the TOP (user-directed
      // 2026-09-02): the viewer TITLE is OFF (its text is not colourable — canvas), so
      // the whole `subject: reasons` label rides the DESCRIPTION, which IS a stylable
      // DOM element the `--nca-tile-color` CSS rule tints to the tile's verdict hue.
      showTitle: false,
      description: `${subject}: ${subjectDescriptions?.get(subject) ?? ''}`,
      descriptionPosition: 'Top' as const,
      descriptionVisibilityMode: 'Always' as const,
    });
    // Tag the root so the selection→dim bridge can look up plots by subject.
    plot.root.dataset.ncaFacetSubject = subject;
    // Dashed inter-tile separators (see `.nca-facet-tile` in app-styles.css).
    plot.root.classList.add('nca-facet-tile');
    // Verdict hue → the tile's title + caption (via `--nca-tile-color` + the CSS rule
    // in app-styles.css), so the subject title and the violated-criteria caption read in
    // the SAME colour as the tile's line/markers (user-directed 2026-09-02).
    plot.root.style.setProperty('--nca-tile-color', DGNS.Color.toHtml(color));
    // Smooth opacity transition for the dim animation.
    plot.root.style.transition = 'opacity 120ms';
    plots.push(plot);
    cols[idx].push(plot);
    idx = (idx + 1) % maxInRow;
  }

  if (scrolling) return {root: buildScrollingFacetGrid(plots), plots};

  // Trim the outer right/bottom perimeter. Row-major width here is `maxInRow`
  // (the split grid's column count); the scrolling path tags with its own width.
  tagFacetEdgeTiles(plots, maxInRow);
  const rowsCount = cols[0].length;
  const col2Roots = (col: DG.Viewer[]): HTMLElement[] => {
    const roots: HTMLElement[] = col.map((p) => p.root);
    while (roots.length < rowsCount) roots.push(ui.div(''));
    return roots;
  };
  const splitCols = cols.map((col) => ui.splitV(col2Roots(col)));
  const root = ui.splitH(splitCols);
  return {root, plots};
}

/**
 * Lay `plots` out as a scrolling `FACET_VIEWPORT_COLS`-wide grid whose tiles keep
 * the height they would have in a `COLS × ROWS` viewport.
 *
 * **Why not `ui.splitH`/`ui.splitV` here.** Splitters DISTRIBUTE the available
 * space: 60 plots in five split columns of twelve would each get a twelfth of the
 * pane height and be unreadable — the exact outcome the old hard cut-off existed
 * to prevent. A scroll container instead keeps the tile at its viewport size and
 * moves the surplus below the fold, which is what "always show the facet" has to
 * mean if it is to mean anything.
 *
 * The row height cannot be a CSS percentage: in a scrolling grid the content box
 * grows with the content, so a percentage row would resolve against a height it
 * is itself defining. It is therefore measured from the host and re-measured on
 * resize. Writing `gridAutoRows` changes CONTENT height, never the container's
 * own (it is `height: 100%` of its dock pane), so the observer cannot feed
 * itself — and the write is skipped when the value is unchanged anyway.
 */
export function buildScrollingFacetGrid(plots: ReadonlyArray<DG.Viewer>): HTMLDivElement {
  // This grid is exactly `FACET_VIEWPORT_COLS` wide — trim the outer right/bottom
  // perimeter at that width (the split path uses its own column count).
  tagFacetEdgeTiles(plots, FACET_VIEWPORT_COLS);
  const root = ui.div(plots.map((p) => p.root), 'nca-facet-scroll') as HTMLDivElement;
  // A grid item's default `min-width: auto` is its min-CONTENT width, and a DG
  // viewer's is 400 px — measured, not guessed. With five 229 px tracks in a
  // 1145 px pane that overflowed to 1316 px, and `overflow-x: hidden` sliced the
  // last two columns off with no way to reach them (reported 2026-08-13).
  for (const p of plots) {
    p.root.style.minWidth = '0';
    p.root.style.width = '100%';
  }
  root.style.display = 'grid';
  root.style.gridTemplateColumns = `repeat(${FACET_VIEWPORT_COLS}, minmax(0, 1fr))`;
  root.style.overflowY = 'auto';
  root.style.overflowX = 'hidden';
  root.style.height = '100%';
  root.style.width = '100%';

  let applied = -1;
  const sizeRows = (): void => {
    // `clientHeight` is the VISIBLE box — the viewport we are dividing — not the
    // scrolled content height.
    const h = root.clientHeight;
    if (h <= 0) return; // not laid out yet; the observer fires again when it is
    // Fill the pane when there is room; never crush below what the chart needs
    // to draw its axes ({@link FACET_TILE_MIN_PX}).
    const rowPx = Math.max(FACET_TILE_MIN_PX, Math.floor(h / FACET_VIEWPORT_ROWS));
    if (rowPx === applied) return;
    applied = rowPx;
    root.style.gridAutoRows = `${rowPx}px`;
  };
  sizeRows();
  // The facet docks from a 100 ms timeout and its pane is user-resizable, so the
  // first useful measurement usually arrives after this function returns.
  new ResizeObserver(sizeRows).observe(root);
  return root;
}

/** Create / update / tear down the per-subject facet grid alongside the
 *  overlay {@link PROFILE_CHART_TAG}. Mirrors the lifecycle shape of
 *  {@link ensureProfileChart}; the differences are:
 *
 *  - **No cap.** The grid is built for any non-empty subject set (user-directed
 *    2026-08-13). Past one viewport ({@link FACET_SCROLL_ABOVE}) the layout
 *    switches to the scrolling {@link FACET_VIEWPORT_COLS}-wide grid rather than
 *    refusing to draw. Only an EMPTY set returns early.
 *  - **Subject-set diff.** Re-Compute and mapping changes can keep the same
 *    subject set (cheap path: swap each plot's DataFrame in place) or change
 *    it (full rebuild). Compared by ordered key array against
 *    {@link FACET_SUBJECTS_TAG}.
 *  - **Deferred dock.** First dock happens inside `setTimeout(100ms)` —
 *    same anti-race pattern DiffStudio uses (`app.ts:1696`) to escape the
 *    in-flight FILL dock op for the overlay chart that runs in the same
 *    tick. */
// ── Single-subject λz editor (F4, docked tab — WALKING SKELETON for the design look) ──
//
// One enlarged facet-style plot bound to the results-grid current row, with a live
// readout. It reuses the WHOLE shipped override path: the plot frame carries the same
// FACET_SRCROW_COL / FACET_ROLE_COL / dataset.ncaFacetSubject a facet tile does, so its
// click routes through `handleFacetPointClick` (override → grid + facet + aggregates all
// refresh), and `handleFacetPointClick` additionally fires EDITOR_REFRESH so the editor
// plot re-renders too. On refresh the plot's dataFrame is SWAPPED (not the viewer), so
// the click subscription survives — exactly like `refreshFacetSubject`.

/** Dock the λz-editor tab (once) and bind it to the results-grid current row. */
/** Reset EVERY point of the context panel's current subject to `auto` (VR-11 Reset-all) — clears
 *  all manual force-in/force-out for that subject in one action, then drives the SAME refresh chain
 *  a point-click does (facet tile + aggregates + editor + context readout). Source rows come off
 *  the panel's own chart frame (built for this subject; falls back to the editor plot); resetting an
 *  already-auto point is a no-op. */
function resetAllLzOverrides(view: DG.TableView): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  // AC-S3.9: the SAME guard as a point click — a reset is a λz edit too (review F1, 2026-09-30).
  if (!lzEditingAllowed(view)) return;
  const controller = bag[NCA_CONTROLLER_TAG] as NcaController | null;
  const subject = bag[LZ_READOUT_SUBJECT_TAG] as string | null;
  const frame = (bag[LZ_READOUT_CHART_TAG] as DG.Viewer | null)?.dataFrame ??
    (bag[EDITOR_PLOT_TAG] as DG.Viewer | null)?.dataFrame ?? null;
  const srcCol = frame?.col(FACET_SRCROW_COL) ?? null;
  if (controller == null || subject == null || subject === '' || frame == null || srcCol == null)
    return;
  const rows = new Set<number>();
  for (let r = 0; r < frame.rowCount; r++)
    if (!srcCol.isNone(r)) rows.add(srcCol.get(r) as number);
  const actions: OverrideAction[] = [...rows].map((rowIdx) => ({rowIdx, action: 'auto'}));
  if (actions.length === 0) return;
  try {
    controller.overrideLambdaZPoints(subject, actions);
    refreshStageUi(view); // a reset can resolve (D11(1)) or replace (D13) a review decision — AC-S3.1.1
    (bag[FACET_REFRESH_TAG] as ((s: string) => void) | null)?.(subject);
    (bag[AGG_REFRESH_TAG] as (() => void) | null)?.();
    (bag[EDITOR_REFRESH_TAG] as (() => void) | null)?.();
    (bag[LZ_READOUT_TAG] as ((s: string) => void) | null)?.(subject);
    markAnalysisDirty(view); // F3: a λz reset is an unsaved edit (study-backed only)
  } catch (e) {
    grok.shell.error(`λz reset failed: ${e instanceof Error ? e.message : String(e)}`);
  }
}

/** The λz marker legend (AC-6.1): the 3-state point vocabulary the facet tiles / editor / context
 *  chart all draw — filled ● IN the fit, hollow ○ excluded (auto), ✕ forced out (a deliberate
 *  force-out). Static; neutral grey so it does not compete with the verdict-coloured plot. Lives on
 *  the context panel's "λz points" pane (user-directed 2026-09-02). */
function buildLzMarkerLegend(colour?: number): HTMLElement {
  // All three GLYPHS take the subject's VERDICT hue (matches the plot markers); ● / ○ enlarged;
  // the LABELS stay neutral grey (user-directed 2026-09-02). Grey fallback when no verdict colour.
  const hue = colour !== undefined ? DGNS.Color.toHtml(colour) : 'var(--grey-6)';
  const item = (glyph: string, label: string, glyphPx: number, tip: RichTipSpec): HTMLElement => {
    const g = ui.divText(glyph);
    g.style.cssText = `font-size: ${glyphPx}px; line-height: 1; color: ${hue}; margin-right: 4px;`;
    const row = ui.divH([g, ui.divText(label)]);
    row.style.cssText = 'align-items: center; gap: 3px;';
    ui.tooltip.bind(row, richTip(tip)); // each marker explains itself (user-directed 2026-09-02)
    return row;
  };
  const legend = ui.divH([
    item('●', 'in λz fit', 16, {title: 'In the λz fit',
      summary: 'This point IS used in the reported terminal-slope (λz) regression — auto-selected ' +
        'by the best-fit search, or force-included by you.'}),
    item('○', 'excluded', 16, {title: 'Excluded from λz',
      summary: 'This point is NOT in the λz fit — auto-excluded by the best-fit window (e.g. Cmax / ' +
        'early distribution points), or a force-in the minimum-points floor could not keep.'}),
    item('✕', 'forced out', 13, {title: 'Forced out',
      summary: 'You deliberately excluded this point from λz (force-out). Click it again on the ' +
        'plot (or here) to restore it to auto.'}),
  ]);
  // padding-top gives the enlarged ●/○ glyphs (16px, line-height 1) vertical room — without it their
  // tops are clipped by the pane boundary (user-found fix, 2026-09-03).
  legend.style.cssText =
    'gap: 14px; align-items: center; font-size: 12px; color: var(--grey-6); padding-top: 15px;';
  return legend;
}

/** F5: the per-point λz DECISIONS grid for `subject` — a NATIVE `DG.Viewer.grid` with one row per
 *  real observation point, `State (●/○/✕, verdict hue) · Time · Conc`, time-ordered. Clicking a row
 *  toggles that point via the SAME {@link applyLzPointToggle} the plot click uses (auto →
 *  force-out/in → auto), so the grid and the plot cannot disagree. Reads the subject's frame from
 *  {@link buildEditorFrame}; synthetic points (`_srcRow < 0` — the prepended C0 / augmented point)
 *  are not real rows and are skipped. On a SAME-SUBJECT gesture (a toggle) the grid is REUSED in
 *  place — only State / _role are updated + `invalidate()` — so it neither flickers nor loses its
 *  scroll; only a subject / row-count change rebuilds the viewer (a build-once context-panel viewer
 *  goes blank after close→reopen, so the rebuild path disposes + recreates). `_srcRow` / `_role`
 *  ride hidden on the frame so the click can resolve the toggle without a scan.
 *
 *  SCOPE vs AC-5.1: the synthesis F5 was a READ-ONLY grid of the OVERRIDDEN points only
 *  (`_lz_action ≠ auto`), carrying `ANALYST` (`_lz_user`) + `RECORDED` (`_lz_time`) provenance.
 *  Directed live (2026-09-03), this became an ALL-POINTS INTERACTIVE grid (every real point,
 *  click-to-toggle — reversing D8's "editing on the plot only" for v1). The who/when `Analyst`/
 *  `Recorded` columns SHIP here (2026-09-08, spec `nca-lambdaz-override-whowhen-surfacing`): read
 *  from `_lz_user`/`_lz_time` via the controller's `lzProvenanceForRow`, VISIBLE only when the
 *  subject has ≥1 override (`setLzProvenanceColsVisible`, req#3 — a point toggle shows them, a full
 *  reset hides them), with English header tooltips. Presentation only — the backing `_lz_user`/
 *  `_lz_time` were already written on every override; who is also surfaced via the results-grid ✎
 *  tooltip. No edit history (current-state per point, RUO/AD-04). */
/** Show/hide the λz-points grid's who/when columns (`Analyst` / `Recorded`) — visible
 *  only when the subject has at least one manual override (req#3, D14 / AC-5.1). Shared
 *  by the full-rebuild and the in-place toggle-refresh paths so the two cannot disagree. */
function setLzProvenanceColsVisible(grid: DG.Grid, visible: boolean): void {
  for (const n of ['Analyst', 'Recorded']) {
    const gc = grid.col(n);
    if (gc !== null) gc.visible = visible;
  }
}

function renderLzPointsGrid(
  view: DG.TableView, host: HTMLElement, subject: string, ctx: LzChartCtx,
  resultsDf: DG.DataFrame, thresholds: ResultsColorThresholds,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const disposePrior = (): void => {
    const prior = bag[LZ_F5_GRID_TAG] as DG.Grid | null;
    if (prior !== null)
      try {prior.close();} catch {/* stale viewer — already detached */}
    bag[LZ_F5_GRID_TAG] = null;
    (bag[LZ_F5_CLICK_TAG] as {unsubscribe(): void} | null)?.unsubscribe?.();
  };

  const srcDf = buildEditorFrame(
    view, subject, ctx.denseDf, ctx.subjectName, ctx.timeName, ctx.concName, resultsDf, thresholds)?.df ?? null;
  const tCol = srcDf?.col(ctx.timeName) ?? null;
  const cCol = srcDf?.col(ctx.concName) ?? null;
  const roleCol = srcDf?.col(FACET_ROLE_COL) ?? null;
  const srcCol = srcDf?.col(FACET_SRCROW_COL) ?? null;
  if (srcDf == null || tCol == null || cCol == null || roleCol == null || srcCol == null) {
    disposePrior(); ui.empty(host); host.style.height = ''; return;
  }

  const controller = bag[NCA_CONTROLLER_TAG] as NcaController | null;
  const glyphOf = (role: string): string => role === 'included' ? '●' : role === 'force-out' ? '✕' : '○';
  const states: string[] = []; const times: number[] = []; const concs: number[] = [];
  const srcRows: number[] = []; const roles: string[] = [];
  // Who/when provenance per row (D14 / AC-5.1): read the SAME `_lz_user`/`_lz_time` the
  // override wrote, via the controller's row-keyed accessor. `analysts`/`recorded` are blank
  // for auto rows; `anyOverride` gates the two columns' visibility (req#3 — a force-in/out
  // shows them, a full reset back to auto hides them again).
  const analysts: string[] = []; const recorded: string[] = []; let anyOverride = false;
  for (let r = 0; r < srcDf.rowCount; r++) {
    if (srcCol.isNone(r)) continue;
    const sr = srcCol.get(r) as number;
    if (sr < 0) continue; // synthetic prepended point — not a togglable observation
    const role = roleCol.getString(r);
    states.push(glyphOf(role)); times.push(numAt(tCol, r)); concs.push(numAt(cCol, r));
    srcRows.push(sr); roles.push(role);
    const prov = controller?.lzProvenanceForRow(sr) ?? null;
    if (prov !== null) {
      anyOverride = true;
      analysts.push(prov.user !== '' ? prov.user : OVERRIDE_UNKNOWN_USER);
      recorded.push(formatOverrideTimestamp(prov.time));
    } else {
      analysts.push(''); recorded.push('');
    }
  }
  if (states.length === 0) {
    disposePrior(); ui.empty(host); host.style.height = ''; return;
  }

  // State-glyph hue in the bag so the (reused) onCellPrepare re-reads THIS subject's verdict.
  bag[LZ_F5_HUE_TAG] = buildFacetVerdictColors(resultsDf, thresholds).get(subject) ?? STATUS_NEUTRAL;

  // REUSE the existing grid IN PLACE on a SAME-SUBJECT gesture — update only State / _role and
  // repaint (NO dispose → no flicker, and the vertical scroll is kept natively). A toggle never
  // changes the point count, so the row set matches; only a subject / row-count change (or no prior
  // grid) falls through to a full rebuild.
  const prior = bag[LZ_F5_GRID_TAG] as DG.Grid | null;
  const priorDf = prior?.dataFrame ?? null;
  if (prior !== null && priorDf !== null && bag[LZ_F5_SUBJECT_TAG] === subject &&
      priorDf.rowCount === states.length) {
    const stateCol = priorDf.col('State');
    const roleColD = priorDf.col('_role');
    const analystColD = priorDf.col('Analyst');
    const recordedColD = priorDf.col('Recorded');
    if (stateCol !== null && roleColD !== null) {
      for (let i = 0; i < states.length; i++) {
        stateCol.set(i, states[i]);
        roleColD.set(i, roles[i]);
        analystColD?.set(i, analysts[i]);
        recordedColD?.set(i, recorded[i]);
      }
      // req#3: the toggle just changed override state → re-evaluate who/when column visibility.
      setLzProvenanceColsVisible(prior, anyOverride);
      prior.invalidate(); // repaint — onCellPrepare re-reads the (possibly changed) verdict hue
      return;
    }
  }

  // Full rebuild (subject change / first render).
  disposePrior();
  ui.empty(host);
  host.style.height = '';
  // Time / Concentration headers = the ACTUAL source-column names the user picked in the mapping
  // inputs (ctx.timeName / ctx.concName), so the panel grid carries the same column identities as the
  // study data and the column-role pickers (user-directed 2026-09-03).
  const timeHeader = ctx.timeName; const concHeader = ctx.concName;
  const frame = DGNS.DataFrame.fromColumns([
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'State', states),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.FLOAT, timeHeader, times),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.FLOAT, concHeader, concs),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'Analyst', analysts),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, 'Recorded', recorded),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.INT, '_srcRow', srcRows),
    DGNS.Column.fromList(DGNS.COLUMN_TYPE.STRING, '_role', roles),
  ]);
  const timeC = frame.col(timeHeader); if (timeC !== null) timeC.meta.format = '0.###';
  const concC = frame.col(concHeader); if (concC !== null) concC.meta.format = '0.###';

  const grid = DGNS.Viewer.grid(frame);
  // Unknown option keys are ignored by DG — safe to request them; guard anyway.
  try {
    grid.setOptions(
      {allowEdit: false, showRowHeader: false, showCurrentRowIndicator: false, rowHeight: LZ_F5_ROW_HEIGHT});
  } catch {/* keep going with defaults */}
  for (const n of ['_srcRow', '_role']) {
    const gc = grid.col(n);
    if (gc !== null) gc.visible = false;
  }
  const stateGc = grid.col('State');
  if (stateGc !== null) stateGc.width = 46;
  const analystGc = grid.col('Analyst');
  if (analystGc !== null) analystGc.width = 110;
  const recordedGc = grid.col('Recorded');
  if (recordedGc !== null) recordedGc.width = 122;
  setLzProvenanceColsVisible(grid, anyOverride); // req#3: hidden until the subject has an override

  // State glyph in the subject's verdict hue (read from the bag so a reuse-path repaint tracks it).
  grid.onCellPrepare((gc) => {
    const h = (bag[LZ_F5_HUE_TAG] as number | undefined) ?? STATUS_NEUTRAL;
    if (gc.isTableCell && gc.gridColumn?.name === 'State') gc.style.textColor = h;
  });
  // Call-to-action tooltip on the grid's OWN cells (user-directed 2026-09-02 — on the grid
  // elements, not a separate icon).
  grid.onCellTooltip((gc, x, y) => {
    if (gc.isColHeader) {
      const n = gc.tableColumn?.name;
      if (n === 'Analyst') {
        ui.tooltip.show('Who forced this point into or out of the λz fit.', x, y); return true;
      }
      if (n === 'Recorded') {
        ui.tooltip.show('When the override was recorded (local time).', x, y); return true;
      }
      return false;
    }
    if (!gc.isTableCell) return false;
    ui.tooltip.show('Click a point — here or on the plot — to force it into or out of the λz fit; ' +
      'click again to reset to auto.', x, y);
    return true;
  });
  // Click a row → toggle that point (the SAME handler as a plot click; onCellClick is a real click
  // event, so a programmatic current-row change cannot re-fire it).
  bag[LZ_F5_CLICK_TAG] = grid.onCellClick.subscribe((gc: DG.GridCell) => {
    if (!gc.isTableCell) return;
    const i = gc.cell.rowIndex;
    const srCol = frame.col('_srcRow'); const rlCol = frame.col('_role');
    if (i < 0 || srCol === null || rlCol === null) return;
    applyLzPointToggle(view, subject, srCol.get(i) as number, rlCol.getString(i));
  });

  // Full width; height sized to show EVERY row (no internal scrollbar). A capped grid scrolled
  // internally and its scroll captured the mouse wheel, so the legend below stayed hidden; showing
  // all rows lets the context panel scroll as one and the legend stays reachable. rowHeight pinned
  // to 24 above, header ~34.
  host.style.width = '100%';
  host.style.height = `${states.length * LZ_F5_ROW_HEIGHT + 34}px`;
  grid.root.style.width = '100%';
  grid.root.style.height = '100%';
  grid.root.style.boxSizing = 'border-box';
  host.appendChild(grid.root);
  bag[LZ_F5_GRID_TAG] = grid;
  bag[LZ_F5_SUBJECT_TAG] = subject;
}

function ensureLambdaZEditor(
  view: DG.TableView, denseDf: DG.DataFrame,
  subjectName: string, timeName: string, concName: string, useLogY: boolean,
  resultsDf: DG.DataFrame, thresholds: ResultsColorThresholds,
  /** False in a non-viable (destructive) design: the readout closures are still bound (a row
   *  click shows that animal's panel) but no editor tab is docked, and an open one closes. */
  dock = true,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  // Current subject = the results grid's current row (fallback: first row).
  const subjCol = resultsDf.col('SUBJECT');
  if (subjCol === null || resultsDf.rowCount === 0) return;
  const rowToSubject = (r: number): string =>
    (r >= 0 && r < resultsDf.rowCount && !subjCol.isNone(r)) ? subjCol.getString(r) : subjCol.getString(0);

  // Curated context-panel readout for an ARBITRARY subject — the point-click handler
  // (facet OR editor) and the row-select both drive it (user-directed subject-relevant
  // context). Captures this Compute's resultsDf + thresholds.
  bag[LZ_READOUT_TAG] = (subject: string): void => {
    // On Report the context panel stays on the population overview (user-directed 2026-09-30): a
    // row click / point click does not swap in the per-profile readout — the Report panel is
    // re-shown instead (no tab switch), so the grid's own click cannot leave it blank.
    if (viewStageOf(view) === 'report') {
      showReportPanel(view, false);
      return;
    }
    showLzReadout(view, resultsDf, subject, thresholds,
      {denseDf, subjectName, timeName, concName, useLogY});
  };

  // (Re)bind the current-row subscription to THIS results frame (fresh per Compute).
  (bag[EDITOR_ROWSUB_TAG] as {unsubscribe(): void} | null)?.unsubscribe?.();
  bag[EDITOR_ROWSUB_TAG] = resultsDf.onCurrentRowChanged.subscribe(() => {
    const s = rowToSubject(resultsDf.currentRowIdx);
    refreshLambdaZEditor(view, s, denseDf, subjectName, timeName, concName, useLogY, resultsDf, thresholds);
  });
  // EDITOR_REFRESH: re-render whatever subject the editor currently shows (called by
  // handleFacetPointClick after an override, from EITHER the editor or a facet tile).
  bag[EDITOR_REFRESH_TAG] = (): void => {
    const s = (bag[EDITOR_SUBJECT_TAG] as string | null) ?? rowToSubject(resultsDf.currentRowIdx);
    refreshLambdaZEditor(view, s, denseDf, subjectName, timeName, concName, useLogY, resultsDf, thresholds);
  };

  const first = rowToSubject(resultsDf.currentRowIdx >= 0 ? resultsDf.currentRowIdx : 0);

  if (!dock) {
    const node = bag[EDITOR_NODE_TAG] as DG.DockNode | null;
    if (node != null) try {view.dockManager.close(node);} catch {/* already gone */}
    bag[EDITOR_NODE_TAG] = null;
    bag[EDITOR_PLOT_TAG] = null;
    bag[EDITOR_CONTENT_TAG] = null;
    return;
  }

  // REUSE the docked editor IN PLACE across Computes. This function used to close() the
  // prior node and immediately re-dock synchronously on EVERY Compute — the exact
  // close-then-redock footgun ensureFacetGrid documents (dockManager.close does not
  // reliably detach before the re-dock, so a unit/dose-change recompute left a SECOND
  // "λz editor" plot). The click subscription is viewer-level (survives a frame swap —
  // see refreshFacetSubject), so keeping the viewer + dock node and swapping only the
  // frame via refreshLambdaZEditor is the same in-place strategy the facet uses. The
  // fresh-Compute closures were just rebound to THIS resultsDf above; pass showReadout
  // = false so the recompute does not hijack the context panel. Build+dock only when
  // there is no live editor (findNode on the docked content, as dockChartTab itself does).
  const livePlot = bag[EDITOR_PLOT_TAG] as DG.Viewer | null;
  const liveContent = bag[EDITOR_CONTENT_TAG] as HTMLElement | null;
  if (livePlot != null && liveContent != null && view.dockManager.findNode(liveContent) != null) {
    refreshLambdaZEditor(
      view, first, denseDf, subjectName, timeName, concName, useLogY, resultsDf, thresholds, false);
    return;
  }

  const frame = buildEditorFrame(view, first, denseDf, subjectName, timeName, concName, resultsDf, thresholds);
  if (frame === null) return;
  bag[EDITOR_SUBJECT_TAG] = first;

  // The PLOT fills the whole tab. Subject data lives in the DG context panel — driven by
  // selecting the subject's results-grid row (handleFacetPointClick sets `currentRowIdx`),
  // exactly as clicking that row does (user-directed 2026-09-02), not a hand-rolled panel.
  const logColName = `log(${concName})`;
  const plot = DGNS.Viewer.scatterPlot(frame.df, {
    xColumnName: timeName,
    yColumnName: logColName, // semi-log (the λz line only reads in log)
    showXAxis: true, showYAxis: true,
    showVerticalGridLines: true, showHorizontalGridLines: true,
    showXSelector: false, showYSelector: false, showSizeSelector: false, showColorSelector: false,
    yAxisType: 'linear' as const,
    markerType: DGNS.MARKER_TYPE.CIRCLE, markerDefaultSize: 10,
    linesOrderColumnName: timeName, linesByColumnName: FACET_COLOR_COL, linesWidth: 2,
    colorColumnName: FACET_COLOR_COL,
    ...(frame.df.col(FACET_ROLE_COL) !== null ? {markersColumnName: FACET_ROLE_COL} : {}),
    showDataframeFormulaLines: true,
    legendVisibility: 'Never' as const,
    showCurrentPoint: false,
    // Copy the facet tile's edge-decoration suppression (Rule 0 — same viewer options):
    // histograms / filtered-out points / drop-lines render as stray strips/points at the
    // plot bottom/right (the Mean±SD lesson) and autoLayout fights the fixed binding.
    autoLayout: false,
    showXHistogram: false, showYHistogram: false,
    showFilteredOutPoints: false, showDropLines: false,
    // No VIEWER title — DG binds it to the dock caption ('λz editor'), which is not the
    // subject. The SUBJECT is shown by our own header strip above the plot instead.
    showTitle: false,
  });
  plot.root.dataset.ncaFacetSubject = first;
  bag[EDITOR_PLOT_TAG] = plot;
  wireFacetPlotClick(view, plot, EDITOR_CLICK_SUB_TAG);

  // Subject header ABOVE the plot (the "title" the user reads); the dock caption stays
  // the tab name. Updated in place on retarget / override. The plot goes DIRECTLY into
  // the flex column (no ui.box wrapper — that squeezed the plot to ~40% width) with
  // `flex:1 1 auto` + `minHeight:0` so it fills the remaining width AND height.
  const header = ui.div();
  header.style.textAlign = 'center';
  header.style.fontSize = '16px';
  header.style.margin = '2px 0';
  header.style.flex = '0 0 auto';
  bag[EDITOR_HEADER_TAG] = header;
  renderEditorHeader(header, resultsDf, first, thresholds,
    (bag[NCA_CONTROLLER_TAG] as NcaController | null)?.nonContiguousSubjects());
  plot.root.style.flex = '1 1 auto';
  plot.root.style.minHeight = '0';
  plot.root.style.width = '100%';
  // The λz-point controls (Reset-all + legend, and the F5 decisions grid) live on the CONTEXT
  // PANEL (native accordion, user-directed 2026-09-02), not floated over the canvas — the editor
  // tab stays just header + plot.
  const content = ui.divV([header, plot.root]);
  content.style.height = '100%';
  content.style.width = '100%';
  bag[EDITOR_CONTENT_TAG] = content; // liveness anchor for the in-place reuse path above
  const existing = bag[EDITOR_NODE_TAG] as DG.DockNode | null;
  if (existing != null)
    try {view.dockManager.close(existing);} catch {/* already gone */}
  // Dock FILL into the SAME charts group as the Facet (user-directed 2026-09-02): the
  // editor is a tab OVER the facet plot, sharing its top region. `dockChartTab` FILLs
  // into that group (self-creating it if it is not up yet) — the exact call the Facet
  // itself uses. Close the prior editor node FIRST so a recompute never leaves a
  // duplicate tab. The facet joins the same group later (its own deferred dockChartTab).
  bag[EDITOR_NODE_TAG] = dockChartTab(view, content, 'λz editor');
}

/** Swap the editor plot's frame for `subject` (the viewer is kept, so the click
 *  subscription survives) and refresh the context-panel readout. No-op if not docked. */
function refreshLambdaZEditor(
  view: DG.TableView, subject: string, denseDf: DG.DataFrame,
  subjectName: string, timeName: string, concName: string, _useLogY: boolean,
  resultsDf: DG.DataFrame, thresholds: ResultsColorThresholds, showReadout = true,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const plot = bag[EDITOR_PLOT_TAG] as DG.Viewer | null;
  if (plot == null) return;
  const frame = buildEditorFrame(view, subject, denseDf, subjectName, timeName, concName, resultsDf, thresholds);
  if (frame === null) return;
  bag[EDITOR_SUBJECT_TAG] = subject;
  try {plot.dataFrame = frame.df;} catch {/* stale hit-test BitSet, as in refreshFacetSubject */}
  plot.setOptions({xColumnName: timeName, yColumnName: `log(${concName})`, showDataframeFormulaLines: true});
  plot.root.dataset.ncaFacetSubject = subject;
  const header = bag[EDITOR_HEADER_TAG] as HTMLElement | null;
  if (header != null) {
    renderEditorHeader(header, resultsDf, subject, thresholds,
      (bag[NCA_CONTROLLER_TAG] as NcaController | null)?.nonContiguousSubjects());
  }
  // A row-select / click retarget refreshes the curated context readout; a Compute-time
  // in-place reuse passes showReadout = false so it does not hijack the context panel.
  if (showReadout)
    (bag[LZ_READOUT_TAG] as ((s: string) => void) | null)?.(subject);
}

/** Build a fresh single-subject facet-style FRAME (series/marker/line coding applied)
 *  + the DG-native readout for `subject`. Returns null for a degenerate frame / unknown
 *  subject. The frame carries the facet service columns so `handleFacetPointClick` works
 *  verbatim. NO plot is created here (the caller keeps one viewer and swaps the frame). */
function buildEditorFrame(
  view: DG.TableView, subject: string, denseDf: DG.DataFrame,
  subjectName: string, timeName: string, concName: string,
  resultsDf: DG.DataFrame, thresholds: ResultsColorThresholds,
): {df: DG.DataFrame} | null {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const controller = bag[NCA_CONTROLLER_TAG] as NcaController | null;
  const lzMap = buildFacetLambdaZParams(resultsDf);
  const entries = buildFacetDataFrames(
    denseDf, subjectName, timeName, concName, lzMap, controller?.fitRowsBySubject(), lzActionOfView(view));
  const found = entries.find((e) => e.subject === subject) ?? entries[0];
  if (found === undefined) return null;
  const df = found.df;
  const colour = buildFacetVerdictColors(resultsDf, thresholds).get(found.subject) ?? STATUS_NEUTRAL;
  paintFacetSeries(df, found.subject, colour);
  applyFacetPointMarkers(df);
  applyFacetLambdaZLine(df, found.lz ?? null, concName, timeName);
  return {df};
}

/**
 * Render the editor's subject header IN PLACE: **bold subject** + `: {violated criteria}`
 * (the SAME short-form phrases as the ACCEPTANCE cell tooltip / facet caption, via
 * `readAcceptance`; an Accepted / Not-estimated subject shows that verdict word). Both
 * the subject and the reasons take the subject's VERDICT COLOUR — the same hue the plot
 * line/markers use (`buildFacetVerdictColors`), so the header matches its chart.
 */
function renderEditorHeader(
  header: HTMLElement, resultsDf: DG.DataFrame, subject: string, thresholds: ResultsColorThresholds,
  nonContiguous?: ReadonlySet<string>,
): void {
  const subjEl = ui.divText(subject);
  subjEl.style.fontWeight = 'bold';
  subjEl.style.display = 'inline';
  // Reuse the SAME caption + hue the facet tiles use (no re-derivation — one source, so
  // the header cannot disagree with the grid ACCEPTANCE / facet caption): the short-form
  // verdict/reasons from `buildFacetVerdictDescriptions` (threaded the SAME `nonContiguous`
  // set the facet caption uses, so the "· non-contiguous λz" note is preserved — this
  // editor is exactly where such selections are made) and the verdict colour from
  // `buildFacetVerdictColors`, both keyed by subject.
  const caption = buildFacetVerdictDescriptions(resultsDf, thresholds, nonContiguous).get(subject) ?? '';
  const colour = buildFacetVerdictColors(resultsDf, thresholds).get(subject);
  if (caption === '' || colour === undefined) {
    header.style.color = '';
    header.replaceChildren(subjEl);
    return;
  }
  const rest = ui.divText(`: ${caption}`);
  rest.style.display = 'inline';
  header.style.color = DGNS.Color.toHtml(colour); // inherited by both inline children
  header.replaceChildren(subjEl, rest);
}

/** The ACCEPTANCE status + violated-criteria reasons on ONE verdict-coloured line — the SAME text
 *  the facet caption / editor header show (`buildFacetVerdictDescriptions`, threaded the same
 *  `nonContiguous` set, so the "· non-contiguous λz" note carries). For Flagged the caption IS the
 *  reasons → "Flagged: {reasons}"; for Accepted / Not estimated the caption equals the verdict
 *  word → shown alone. Sits directly UNDER the subject title in the context panel (user-directed
 *  2026-09-02), not inside the parameters pane. */
function buildStatusLine(
  resultsDf: DG.DataFrame, subject: string, thresholds: ResultsColorThresholds,
  nonContiguous?: ReadonlySet<string>,
): HTMLElement {
  const el = ui.divText('');
  el.style.fontWeight = 'bold';
  el.style.whiteSpace = 'normal'; // reasons can be long — allow wrap
  const i = findSubjectRow(resultsDf, subject);
  if (i < 0) return el;
  const a = readAcceptance(
    resultsDf.col('SPAN_RATIO'), resultsDf.col('LAMBDAZ_ADJ_R2'),
    resultsDf.col('LAMBDAZ_NPOINTS'), resultsDf.col('PCT_EXTRAP'), i, thresholds);
  // Just the caption — NO leading verdict word (user-directed 2026-09-02): for Flagged this is
  // the violated-criteria reasons (the colour already conveys "flagged"); Accepted / Not estimated
  // keep their word because the caption IS that word. Identical text to the facet tile caption.
  const caption = buildFacetVerdictDescriptions(resultsDf, thresholds, nonContiguous).get(subject) ?? a.verdict;
  el.textContent = caption;
  const colour = buildFacetVerdictColors(resultsDf, thresholds).get(subject);
  if (colour !== undefined) el.style.color = DGNS.Color.toHtml(colour);
  // Same tooltip as the grid's ACCEPTANCE cell — ONE source (`buildAcceptanceTooltipSpec`), so the
  // status line and the grid cell cannot disagree (user-directed 2026-09-02).
  ui.tooltip.bind(el, richTip(buildAcceptanceTooltipSpec(resultsDf, i, thresholds)));
  return el;
}

/**
 * CURATED "NCA parameters" pane content: ALL main NCA parameters as label-over-value cards —
 * "{name}, {unit}" on top (unit from the column's canonical `units` tag, which the writer sets —
 * mL/min/kg and L/kg for CL / Vz) + the value below, coloured by the SAME grid band where the
 * column is banded (`gridBandTextColour`, so a panel value cannot disagree with its grid cell).
 * The ACCEPTANCE status rides the title area now ({@link buildStatusLine}), not this pane; the
 * manual-override note is NOT shown here (user-directed 2026-09-02) — it stays in the ACCEPTANCE
 * tooltip (`buildAcceptanceTooltipSpec`, on the status line + the grid cell). Curated on purpose —
 * a native `DG.Row` view leaks the hidden `_PROBLEM_RANK` helper + raw WARNINGS/PROVENANCE JSON.
 */
function buildLambdaZReadout(
  resultsDf: DG.DataFrame, subject: string, thresholds: ResultsColorThresholds,
): HTMLElement {
  const i = findSubjectRow(resultsDf, subject);
  if (i < 0) return ui.divText('No data for this subject.');
  const headerTips = buildHeaderTooltips(thresholds); // one map → each label's name tooltip
  // AC-S5.1.4: an EXCLUDED profile's λz-dependent values read here exactly as in the grid — kept,
  // muted, ` (excl.)`, with the same tooltip (one decision read, one tooltip spec).
  const excluded = excludedDecisionAt(resultsDf, i);
  // TWO-COLUMN grid, mirroring the study-overview Forms viewer's LOOK (screen-2, user-directed
  // 2026-09-14): "{name}, {unit}" in a left column that auto-sizes to the widest label
  // (`max-content`), the value LEFT-aligned in a shared column right beside it (`1fr`) — so every
  // value starts at the SAME x and the two read as a tidy table, NOT the old flush-right
  // `space-between` gap. NO dividers at all (user-directed 2026-09-14): no vertical rule between the
  // columns (a `column-gap` separates them) AND no per-row horizontal line — the earlier bottom
  // border was dropped on the user's request; row rhythm comes from the grid's `row-gap`. This
  // stays a curated HTML build, not the native "Forms" viewer (which would leak `_PROBLEM_RANK` +
  // raw WARNINGS/PROVENANCE JSON — see the doc comment above). Unit from the column's `units` tag
  // (canonical — the writer set it: CL in mL/min/kg, Vz in L/kg). The NAME tooltip is the
  // grid's HEADER tooltip; the VALUE tooltip is the per-value band/fragility spec — BOTH shared with
  // the results grid (user-directed 2026-09-02), so the panel's explanations cannot disagree.
  const CELL_CSS = 'margin: 0; align-self: baseline;';
  const paramCells = (label: string, col: string, dp: number): [HTMLElement, HTMLElement] => {
    const c = resultsDf.col(col);
    const raw = (c !== null && !c.isNone(i)) ? numAt(c, i) : NaN;
    const unit = c?.getTag('units') ?? '';
    const l = ui.h2(unit !== '' ? `${label}, ${unit}` : label); // param name as ui.h2 (user-directed)
    // Darker label colour to match the study-overview Forms-viewer field labels (user-directed
    // 2026-09-02) — DG's standard dark label grey, not the pale --grey-4. `nowrap` keeps a long
    // label ("AUC extrapolated, %") on one line so `max-content` sizes the column to it.
    l.style.cssText = CELL_CSS + 'color: var(--grey-6); white-space: nowrap;';
    const hdr = headerTips[col];
    if (hdr !== undefined) ui.tooltip.bind(l, richTip(hdr));
    const isExcl = excluded !== null && isLzDependentParam(col) && Number.isFinite(raw);
    const v = ui.divText(Number.isFinite(raw) ? `${raw.toFixed(dp)}${isExcl ? ' (excl.)' : ''}` : '—');
    v.style.cssText = CELL_CSS; // left-aligned in its own column (screen-2 alignment)
    if (isExcl) {
      v.style.color = DGNS.Color.toHtml(EXCLUDED_TEXT); // muted, never a band hue (not reported)
      ui.tooltip.bind(v, richTip(excludedCellTipSpec(excluded)));
      return [l, v];
    }
    const colr = gridBandTextColour(col, raw, thresholds);
    if (colr !== null) v.style.color = DGNS.Color.toHtml(colr);
    ui.tooltip.bind(v, richTip(paramValueTooltipSpec(resultsDf, col, i, thresholds)));
    return [l, v];
  };
  // All main params, in the screenshot order. Decimals MIRROR the grid's NUMBER_FORMATS so the
  // panel value matches the grid cell for the SAME subject — including LAMBDA_Z, which IS a grid
  // column (Full profile, NUMBER_FORMATS.LAMBDA_Z = '0.000' = 3 dp). A 4-dp panel would disagree
  // with the grid's "0.069" for the identical value — the "cannot disagree" invariant this panel
  // exists to hold (decision-auditor, rule 16).
  const specs: ReadonlyArray<readonly [string, string, number]> = [
    ['λz', 'LAMBDA_Z', 3], ['t½', 'HALF_LIFE', 3], ['adj. R²', 'LAMBDAZ_ADJ_R2', 3],
    ['points in λz', 'LAMBDAZ_NPOINTS', 0], ['span', 'SPAN_RATIO', 2],
    ['AUC extrapolated', 'PCT_EXTRAP', 1], ['AUC∞', 'AUC_INF', 3],
    ['CL', 'CL', 3], ['Vz', 'VZ', 3], ['MRT', 'MRT', 3],
  ];
  const cells: HTMLElement[] = [];
  for (const [label, col, dp] of specs)
    cells.push(...paramCells(label, col, dp));
  const grid = ui.div(cells);
  // `max-content` label column (sized to the widest label) + `1fr` value column, values sharing one
  // left edge; `column-gap` separates them WITHOUT a divider rule; `row-gap` gives the vertical
  // rhythm now that the per-cell padding/border is gone; `align-items: baseline` lines up each row's
  // label and value on their text baseline.
  grid.style.cssText =
    'display: grid; grid-template-columns: max-content 1fr; column-gap: 24px; row-gap: 8px; ' +
    'align-items: baseline;';
  return grid;
}

/** Results-frame row index for `subject`, or -1. One helper for every subject→row lookup
 *  (the editor click + the readout), replacing the hand-rolled scans (architect). */
function findSubjectRow(resultsDf: DG.DataFrame, subject: string): number {
  const subjCol = resultsDf.col('SUBJECT');
  if (subjCol === null) return -1;
  for (let r = 0; r < resultsDf.rowCount; r++)
    if (!subjCol.isNone(r) && subjCol.getString(r) === subject) return r;
  return -1;
}

/**
 * (Re)build the context panel's OWN single-subject profile chart into `host` (F#2, user-directed
 * 2026-09-02). Reuses {@link buildEditorFrame} (same series / markers / λz-line coding as the
 * editor), but as a DISPLAY chart: NO click-to-override wiring and NO verdict-coloured caption
 * (the panel header already carries the subject + verdict). Disposed +
 * rebuilt per gesture (a build-once context-panel viewer goes blank after close→reopen), tracked
 * on `LZ_READOUT_CHART_TAG`.
 */
function renderLzProfileChart(
  view: DG.TableView, host: HTMLElement, subject: string, ctx: LzChartCtx,
  resultsDf: DG.DataFrame, thresholds: ResultsColorThresholds,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  // Dispose the prior chart viewer in try/catch — `Viewer.close()` on a stale/detached viewer
  // THROWS in `grok_Viewer_Close` (optional chaining does NOT catch that), and disposing a stale
  // viewer must not break the panel (the SAME guard renderParamsForm uses in the study-overview).
  const prior = bag[LZ_READOUT_CHART_TAG] as DG.Viewer | null;
  if (prior !== null)
    try {prior.close();} catch {/* stale viewer — already detached */}
  bag[LZ_READOUT_CHART_TAG] = null;
  ui.empty(host);
  const frame = buildEditorFrame(
    view, subject, ctx.denseDf, ctx.subjectName, ctx.timeName, ctx.concName, resultsDf, thresholds);
  if (frame === null) {
    host.appendChild(ui.divText('No profile to plot.'));
    return;
  }
  const yCol = ctx.useLogY ? `log(${ctx.concName})` : ctx.concName;
  const plot = DGNS.Viewer.scatterPlot(frame.df, {
    xColumnName: ctx.timeName, yColumnName: yCol,
    showXAxis: true, showYAxis: true,
    showVerticalGridLines: true, showHorizontalGridLines: true,
    showXSelector: false, showYSelector: false, showSizeSelector: false, showColorSelector: false,
    yAxisType: 'linear' as const,
    markerType: DGNS.MARKER_TYPE.CIRCLE, markerDefaultSize: 8,
    linesOrderColumnName: ctx.timeName, linesByColumnName: FACET_COLOR_COL, linesWidth: 2,
    colorColumnName: FACET_COLOR_COL,
    ...(frame.df.col(FACET_ROLE_COL) !== null ? {markersColumnName: FACET_ROLE_COL} : {}),
    showDataframeFormulaLines: ctx.useLogY, // λz line reads only in log space
    legendVisibility: 'Never' as const,
    showCurrentPoint: false,
    autoLayout: false,
    showXHistogram: false, showYHistogram: false,
    showFilteredOutPoints: false, showDropLines: false,
    showTitle: false, // no verdict-coloured caption — the panel header above carries identity
  });
  plot.root.style.width = '100%';
  plot.root.style.height = '100%';
  plot.root.dataset.ncaFacetSubject = subject; // handleFacetPointClick resolves the subject from this
  host.appendChild(plot.root);
  bag[LZ_READOUT_CHART_TAG] = plot;
  // Click-to-override on the context chart too (SAME handler as the facet/editor; distinct
  // subscription tag). handleFacetPointClick → override → showReadout rebuilds this chart, so
  // the plot's own click event is fully read before the rebuild closes it.
  wireFacetPlotClick(view, plot, LZ_READOUT_CLICK_SUB_TAG);
}

/**
 * Show `subject`'s curated readout in the DG context panel, using the proven
 * `updateCompoundPanel` recipe (nca-study-overview.ts): a BUILD-ONCE accordion shell
 * (a bare freshly-built element set to `grok.shell.o` does not reliably render; the
 * reused `divV([acc.root])` root DOES), refreshed by emptying + re-filling its content
 * host, then `showContextPanel = true` + `grok.shell.o = root`. Deferred to the next tick
 * so it wins over the viewer object DG selects synchronously on a point-click. The second
 * pane holds the subject's own profile chart ({@link renderLzProfileChart}, rebuilt per gesture).
 */
function showLzReadout(
  view: DG.TableView, resultsDf: DG.DataFrame, subject: string, thresholds: ResultsColorThresholds,
  ctx: LzChartCtx,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  let panel = bag[LZ_READOUT_PANEL_TAG] as
    {root: HTMLElement; host: HTMLElement; title: HTMLElement; statusHost: HTMLElement;
      legendHost: HTMLElement; f5Host: HTMLElement; chartHost: HTMLElement;
      chartPane: DG.AccordionPane; reviewHost: HTMLElement; reviewPane: DG.AccordionPane;
      resetBtn: HTMLElement} | null;
  if (panel == null) {
    // Subject-name TITLE + the ACCEPTANCE status line BOTH ride the accordion's NATIVE title area
    // (`addTitle`), EXACTLY as the study-overview panel does — so they INHERIT Datagrok's own
    // panel-header padding and line up with the pane header below AND the native context panel
    // (user-directed "same indents", 2026-09-02). The status sits DIRECTLY under the title (its own
    // request 2026-09-02), not inside the parameters pane. `nca-ct-title` text style on the title
    // (18px / --grey-4 / ellipsis); both are set per refresh below.
    const title = panelTitleEl();
    const statusHost = ui.div([]);
    // Even vertical rhythm (user-directed 2026-09-02): the title→status gap and the status→params
    // gap match the gap between param rows. Symmetric top+bottom margin.
    statusHost.style.margin = '8px 0';
    const host = ui.div([]);
    host.style.width = '100%';
    // Second pane: the subject's own C–t profile chart. A definite height is required — the
    // scatterPlot has no intrinsic size and would collapse to 0 in an auto-height accordion pane.
    const chartHost = ui.div([]);
    chartHost.style.width = '100%';
    chartHost.style.height = '300px';
    // λz-point controls pane (native accordion, user-directed 2026-09-02): the marker legend +
    // Reset-all; the F5 per-point decisions grid joins here next. Reset-all reads the panel's
    // CURRENT subject at click time (LZ_READOUT_SUBJECT_TAG), so the build-once button serves every
    // subject.
    // Reset-all as a text button (user-directed 2026-09-02 — back to the button).
    const resetBtn = ui.button('Reset all', () => resetAllLzOverrides(view),
      'Reset every manual λz point decision for this subject to auto');
    const legendHost = ui.div([]); // re-filled per subject (verdict-hue markers) in the refresh below
    const f5Host = ui.div([]); // F5 per-point decisions grid — rebuilt per subject in the refresh below
    // Design B (user-directed 2026-09-02): an action toolbar row (Reset, right-aligned) ABOVE the
    // grid; the legend BELOW the grid as a key. Actions on top, data in the middle, key at the foot.
    const actionRow = ui.divH([resetBtn]);
    actionRow.style.cssText = 'align-items: center; justify-content: flex-start;';
    const lzPointsHost = ui.divV([actionRow, f5Host, legendHost]);
    lzPointsHost.style.cssText = 'gap: 6px; padding: 2px 0; width: 100%;';
    const acc = ui.accordion('nca-lz-readout');
    acc.addTitle(ui.divV([title, statusHost]));
    // F-S4.2 (approved design): the Review pane leads — the decision and `Next to review` sit
    // together, directly under the verdict line they act on (D15).
    const reviewHost = ui.div([]);
    reviewHost.style.width = '100%';
    const reviewPane = acc.addPane('Review', () => reviewHost, true);
    acc.addPane('NCA parameters', () => host, true);
    acc.addPane('λz points', () => lzPointsHost, true);
    // Profile chart COLLAPSED by default; opens to the remembered pref if the user expanded it in a
    // prior session (within a session the build-once accordion keeps the toggle on its own).
    const chartPane = acc.addPane('Concentration–time profile', () => chartHost, readChartExpandedPref());
    acc.root.style.width = '100%';
    const root = ui.divV([acc.root]);
    root.style.width = '100%';
    // Fill the context-panel viewport and scroll INTERNALLY (same recipe as the facet-scroll root
    // above): without this, when the expanded panes (esp. the full-height λz-points grid) overrun
    // the panel, the bottom — the marker legend — is CLIPPED instead of scrolled into view.
    root.style.height = '100%';
    root.style.overflowY = 'auto';
    root.style.overflowX = 'hidden';
    panel = {root, host, title, statusHost, legendHost, f5Host, chartHost, chartPane, reviewHost, reviewPane,
      resetBtn};
    bag[LZ_READOUT_PANEL_TAG] = panel;
  }
  panel.title.textContent = subject;
  // AC-S3.9: λz point editing — Reset all included — exists only where it is allowed (Review λz,
  // Setup unchanged since Confirm); elsewhere the action is not offered at all.
  panel.resetBtn.style.display = lzEditingAllowed(view) ? '' : 'none';
  bag[LZ_READOUT_SUBJECT_TAG] = subject; // the "λz points" pane's Reset-all acts on this subject
  // Recolour the "λz points" legend to THIS subject's verdict hue (matches the plot markers).
  ui.empty(panel.legendHost);
  panel.legendHost.appendChild(
    buildLzMarkerLegend(buildFacetVerdictColors(resultsDf, thresholds).get(subject)));
  // F5: (re)build the native per-point decisions grid for THIS subject (states reflect the current
  // fit; the viewer is disposed + rebuilt inside).
  renderLzPointsGrid(view, panel.f5Host, subject, ctx, resultsDf, thresholds);
  // Status + reason directly under the title (rebuilt per refresh; the controller's non-contiguous
  // set carries the "· non-contiguous λz" note, as on the facet caption / editor header).
  const controller = bag[NCA_CONTROLLER_TAG] as NcaController | null;
  ui.empty(panel.statusHost);
  panel.statusHost.appendChild(
    buildStatusLine(resultsDf, subject, thresholds, controller?.nonContiguousSubjects()));
  renderReviewForSubject(view, panel.reviewHost, resultsDf, subject, thresholds);
  // On the Review λz stage the Review pane is always shown OPEN (user-directed 2026-09-30): it is
  // where the decision is made, so a collapse from an earlier visit (the build-once accordion keeps
  // it) must not hide it.
  if (viewStageOf(view) === 'review') panel.reviewPane.expanded = true;
  ui.empty(panel.host);
  panel.host.appendChild(buildLambdaZReadout(resultsDf, subject, thresholds));
  // F#2: (re)build the subject's profile chart in the second pane (per-gesture, disposed inside).
  renderLzProfileChart(view, panel.chartHost, subject, ctx, resultsDf, thresholds);
  // Persist the user's expand/collapse choice for the chart pane (captured on each gesture, so it
  // survives a page reload / new session; within the session the build-once accordion already keeps it).
  writeChartExpandedPref(panel.chartPane.expanded);
  const root = panel.root;
  setTimeout(() => {
    try {
      grok.shell.windows.showContextPanel = true;
      grok.shell.o = root;
      // The F5 grid was built while the panel was still unmounted (0-width); now it is on screen,
      // nudge a re-layout so its canvas measures the FULL pane width.
      (bag[LZ_F5_GRID_TAG] as DG.Grid | null)?.invalidate?.();
    } catch {/* no shell */}
  }, 0);
}

/**
 * AC-S3.9 (closes GAP-LZ-STALE-RULES / lambdaz-manual-point-selection AC-4.6): λz points are
 * editable only on Review λz and only while the held results match the confirmed Setup — an
 * override can never be committed against rules the analyst edited but did not confirm. The ONE
 * predicate every λz edit path consults (point click, Reset all).
 */
export function lzEditingAllowed(view: DG.TableView): boolean {
  return viewStageOf(view) === 'review' && !liveWorkflowState(view).stage.setupChanged;
}

/** view.temp key: the mount closure's stage navigation (so module-level {@link startReview} can move
 *  to the Review stage). */
export const GOTO_STAGE_TAG = 'ncaGotoStage';

/**
 * Start (or continue) the λz review at `subject` — or, without one, at the first profile still to
 * review in the problem-first order (user-directed 2026-09-30: the strip hint and the "To review"
 * cells are the way in; review is NOT auto-started after Confirm — the analyst first sees the
 * results as a whole). Moves to the Review λz stage, selects the profile's row, brings the
 * "λz editor" tab to the front, and opens the context panel on its Review pane.
 */
export function startReview(view: DG.TableView, subject?: string): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const resultsDf = (bag[RESULTS_GRID_TAG] as DG.Grid | null)?.dataFrame ?? null;
  let target = subject;
  if (target === undefined) {
    const pending = liveWorkflowState(view).units.filter((u) => u.need === 'review' && u.decision === null);
    // Dense profiles first, in problem-first order; then the UC-04 composites — a destructive
    // design has only composites to review (before 2026-10-01 the hint did nothing there).
    const dense = new Set(pending.filter((u) => u.kind === 'profile').map((u) => u.key));
    target = (resultsDf != null ? pendingInProblemOrder(resultsDf, dense)[0] : undefined) ??
      pending.find((u) => u.kind === 'composite')?.key;
  }
  if (target === undefined) return;
  if (viewStageOf(view) !== 'review') (bag[GOTO_STAGE_TAG] as ((s: ViewStage) => void) | null)?.('review');
  // A composite is reviewed from the Sparse NCA table: front it and select that composite's row
  // (its row subscription opens the composite Review pane).
  if (target.startsWith(SPARSE_KEY_PREFIX)) {
    activateDockedTab(view, SPARSE_PARAM_NODE_TAG);
    (bag[SPARSE_GOTO_TAG] as ((key: string) => void) | null)?.(target);
    return;
  }
  if (resultsDf == null) return;
  const r = rowOfSubject(resultsDf, target);
  if (r >= 0) resultsDf.currentRowIdx = r; // drives the editor + facet (existing subscription)
  activateDockedTab(view, EDITOR_NODE_TAG);
  (bag[LZ_READOUT_TAG] as ((s: string) => void) | null)?.(target);
}

/** The results-frame row of `subject`, or -1. */
function rowOfSubject(resultsDf: DG.DataFrame, subject: string): number {
  const c = resultsDf.col('SUBJECT');
  if (c === null) return -1;
  for (let i = 0; i < resultsDf.rowCount; i++) if (c.getString(i) === subject) return i;
  return -1;
}

/** Pending profile keys in the shipped problem-first order (`_PROBLEM_RANK`, then SUBJECT) —
 *  the order `Next to review` walks, independent of any manual grid re-sort (F-S4.3). */
function pendingInProblemOrder(resultsDf: DG.DataFrame, pending: ReadonlySet<string>): string[] {
  const subj = resultsDf.col('SUBJECT');
  const rank = resultsDf.col(PROBLEM_RANK_COL);
  if (subj === null) return [];
  const rows: {key: string; r: number}[] = [];
  for (let i = 0; i < resultsDf.rowCount; i++) {
    const key = subj.getString(i);
    if (!pending.has(key)) continue;
    rows.push({key, r: rank !== null && !rank.isNone(i) ? Number(rank.get(i)) : 99});
  }
  rows.sort((a, b) => a.r - b.r || a.key.localeCompare(b.key, undefined, {numeric: true}));
  return rows.map((x) => x.key);
}

/**
 * Fill the Review pane for one results-row subject from the LIVE workflow state (never
 * re-derived here): its unit, the other selected pending units (bulk), and the problem-first
 * pending order. Decisions go through the controller (audited there); afterwards the Review
 * column, the strip, the dependent aggregates (an Exclude λz changes what they read — AD-6) and
 * the dirty flag refresh.
 */
function renderReviewForSubject(
  view: DG.TableView, host: HTMLElement, resultsDf: DG.DataFrame, subject: string,
  thresholds: ResultsColorThresholds,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const controller = bag[NCA_CONTROLLER_TAG] as NcaController | null;
  const {stage, units, fingerprints} = liveWorkflowState(view);
  const unit = units.find((u) => u.key === subject);
  if (controller == null || unit === undefined) {
    ui.empty(host);
    return;
  }
  const pending = new Set(units.filter((u) => u.need === 'review' && u.decision === null).map((u) => u.key));
  const order = pendingInProblemOrder(resultsDf, pending);
  const selectedKeys = new Set<string>();
  const subj = resultsDf.col('SUBJECT');
  if (subj !== null) for (const i of resultsDf.selection.getSelectedIndexes()) selectedKeys.add(subj.getString(i));
  const afterChange = (): void => {
    markAnalysisDirty(view);
    refreshStageUi(view);
    (bag[AGG_REFRESH_TAG] as (() => void) | null)?.();
    (bag[LZ_READOUT_TAG] as ((s: string) => void) | null)?.(subject);
  };
  renderReviewPane(host, {
    unit,
    // F-S4.4: a bulk block needs ≥ 2 SELECTED units — the current row joins the batch only when it
    // is itself selected (review F8: an unselected current row was silently added to the batch).
    selected: selectedKeys.has(subject) ?
      units.filter((u) => u.key !== subject && selectedKeys.has(u.key)) : [],
    pendingOrder: order,
    thresholds,
    editable: viewStageOf(view) === 'review' && !stage.setupChanged,
    pipelineLine: pipelineLineFor(view, unit.key),
    onDecide: (keys, action, reasons) => {
      for (const k of keys) controller.decide(k, action, reasons, fingerprints.get(k) ?? '');
      afterChange();
    },
    onReset: (k) => {
      controller.resetDecision(k);
      afterChange();
    },
    onNext: () => {
      if (order.length === 0) return;
      const at = order.indexOf(subject);
      const next = order[(at + 1) % order.length] ?? order[0];
      const r = rowOfSubject(resultsDf, next);
      if (r >= 0) resultsDf.currentRowIdx = r; // drives the editor + facet (existing subscription)
      (bag[LZ_READOUT_TAG] as ((s: string) => void) | null)?.(next);
    },
  });
}

/** Wire a facet-style plot's point-click through the shared λz override handler, storing the
 *  subscription under `subTag`. Shared by the docked λz editor (`EDITOR_CLICK_SUB_TAG`) and the
 *  context-panel profile chart (`LZ_READOUT_CLICK_SUB_TAG`) — DISTINCT tags so the two plots'
 *  subscriptions never clobber each other. The plot MUST carry `dataset.ncaFacetSubject`
 *  (handleFacetPointClick reads it to resolve the subject). */
function wireFacetPlotClick(view: DG.TableView, plot: DG.Viewer, subTag: string): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  (bag[subTag] as {unsubscribe(): void} | null)?.unsubscribe?.();
  try {
    const sp = plot as unknown as DG.ScatterPlotViewer;
    const clicked = sp.onPointClicked;
    if (clicked != null) {
      bag[subTag] = clicked.subscribe((e: DG.EventData<DG.RowDataArgs>) =>
        handleFacetPointClick(view, plot, e?.args?.rowId ?? -1));
    }
  } catch {/* not a scatterPlot (degenerate frame) */}
}

/**
 * Wire each facet tile's point-click to a λz override (VR-11, C0 walking skeleton).
 * Uses the scatterPlot's own {@link DG.ScatterPlotViewer.onPointClicked} — DG's
 * precise per-point hit-test (`{rowId}`), NOT `onCurrentRowChanged`, which is a
 * navigation/selection signal that fires unreliably for clustered inner points and
 * line hits (only isolated terminal points registered). The clicked row maps to the
 * source observation row ({@link FACET_SRCROW_COL}); the subject recomputes via the
 * {@link NcaController} and the facet rebuilds ({@link FACET_REFRESH_TAG}).
 * Subscriptions are disposed + re-established on every (re)bind.
 */
function wireFacetClickToOverride(view: DG.TableView, plots: DG.Viewer[]): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const prior = bag[FACET_CLICK_SUBS_TAG] as Array<{unsubscribe(): void}> | null;
  if (prior != null) for (const s of prior) try {s.unsubscribe();} catch {/* already gone */}
  const subs: Array<{unsubscribe(): void}> = [];
  for (const plot of plots) {
    try {
      const sp = plot as unknown as DG.ScatterPlotViewer;
      const clicked = sp.onPointClicked;
      if (clicked == null) continue;
      subs.push(clicked.subscribe((e: DG.EventData<DG.RowDataArgs>) =>
        handleFacetPointClick(view, plot, e?.args?.rowId ?? -1)));
    } catch {/* not a scatterPlot tile (degenerate/mirror frame) */}
  }
  bag[FACET_CLICK_SUBS_TAG] = subs;
}

/**
 * Toggle ONE λz point for `subject` — the SINGLE source shared by a plot point-click
 * ({@link handleFacetPointClick}) and an F5 decisions-grid row-click ({@link buildLzPointsGrid}),
 * so the two control surfaces cannot diverge. `role` is the point's current `_lz_role`
 * ('included' / 'excluded' / 'force-out') — REAL λz fit membership (F6), the SAME source the
 * marker draws from. A manually-forced point resets to `auto` (undo); an `auto` point is forced
 * OUT if it is currently IN the fit, else forced IN. Then drives the full refresh chain: the facet
 * tile, the dependent aggregates (AC-7.2 — By-group / %F / dose-prop / TK), the λz editor, and the
 * context readout (LAST, so the acted subject wins). Re-entrancy-guarded via FACET_CLICK_BUSY_TAG.
 */
function applyLzPointToggle(view: DG.TableView, subject: string, srcRow: number, role: string): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const controller = bag[NCA_CONTROLLER_TAG] as NcaController | null;
  if (controller == null || subject === '' || bag[FACET_CLICK_BUSY_TAG] === true) return;
  if (!lzEditingAllowed(view)) return; // AC-S3.9
  const stored = controller.currentLzAction(srcRow);
  const action: 'auto' | 'force-in' | 'force-out' =
    stored !== 'auto' ? 'auto' : (role === 'included' ? 'force-out' : 'force-in');
  bag[FACET_CLICK_BUSY_TAG] = true;
  try {
    controller.overrideLambdaZPoints(subject, [{rowIdx: srcRow, action}]);
    refreshStageUi(view); // an edit can resolve (D11(1)) or replace (D13) a review decision
    (bag[FACET_REFRESH_TAG] as ((s: string) => void) | null)?.(subject);
    (bag[AGG_REFRESH_TAG] as (() => void) | null)?.();
    (bag[EDITOR_REFRESH_TAG] as (() => void) | null)?.();
    (bag[LZ_READOUT_TAG] as ((s: string) => void) | null)?.(subject);
    markAnalysisDirty(view); // F3: a λz point toggle is an unsaved edit (study-backed only)
  } catch (e) {
    grok.shell.error(`λz override failed: ${e instanceof Error ? e.message : String(e)}`);
  } finally {
    bag[FACET_CLICK_BUSY_TAG] = false;
  }
}

function handleFacetPointClick(view: DG.TableView, plot: DG.Viewer, rowId: number): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  if (bag[FACET_CLICK_BUSY_TAG] === true) return; // belt-and-braces re-entrancy guard
  const frame = plot.dataFrame;
  if (frame == null || rowId < 0 || rowId >= frame.rowCount) return;
  const subject = plot.root.dataset.ncaFacetSubject;
  const srcCol = frame.col(FACET_SRCROW_COL);
  const roleCol = frame.col(FACET_ROLE_COL);
  if (subject == null || subject === '' || srcCol == null || roleCol == null) return;
  applyLzPointToggle(view, subject, srcCol.get(rowId) as number, roleCol.getString(rowId));
}

/**
 * Update ONE facet tile after a single-subject λz override — the fast interactive
 * path (VR-11). Swaps only that subject's single scatterPlot dataFrame + re-derives
 * its verdict hue/caption, instead of re-rendering all N tiles. The fix is in the
 * VIEWER dimension — **one** `scatterPlot.dataFrame` swap + `setOptions` instead of N
 * (the N-scatterPlot re-render was the measured bottleneck, slow past ~20 subjects).
 * Data-prep is NOT O(1): `buildFacetDataFrames` + `buildFacetVerdictColors/Descriptions`
 * still walk all rows / all N subjects and discard all but one — cheap in-memory
 * typed-array passes at N=64, but a per-subject builder would be the next optimization
 * if the ratified reactive epic (D16) raises click frequency. `onPointClicked` is a
 * viewer-level event, so swapping the frame does NOT break the click subscription — no
 * re-wire needed. The whole facet (order + all tiles) is rebuilt on a full Compute via
 * {@link ensureFacetGrid}; an override never re-sorts (D10), so the tile keeps its spot.
 */
function refreshFacetSubject(
  view: DG.TableView, subject: string, denseDf: DG.DataFrame,
  subjectName: string, timeName: string, concName: string, useLogY: boolean,
  results: DG.DataFrame, thresholds: ResultsColorThresholds,
  fitRowsBySubject?: ReadonlyMap<string, ReadonlySet<number>>,
  nonContiguous?: ReadonlySet<string>,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = view.temp as Record<string, any>;
  const plots = bag[FACET_PLOTS_TAG] as DG.Viewer[] | null;
  if (plots == null) return;
  const plot = plots.find((p) => p.root.dataset.ncaFacetSubject === subject);
  if (plot == null) return;
  const lzMap = buildFacetLambdaZParams(results);
  const entry = buildFacetDataFrames(
    denseDf, subjectName, timeName, concName, lzMap, fitRowsBySubject, lzActionOfView(view))
    .find((f) => f.subject === subject);
  if (entry == null) return;
  const df = entry.df;
  const colour = buildFacetVerdictColors(results, thresholds).get(subject);
  paintFacetSeries(df, subject, colour ?? STATUS_NEUTRAL);
  // Keep the title/caption hue in sync with the (possibly changed) verdict after an
  // override — the same `--nca-tile-color` the tile's create path set.
  plot.root.style.setProperty('--nca-tile-color', DGNS.Color.toHtml(colour ?? STATUS_NEUTRAL));
  applyFacetPointMarkers(df);
  applyFacetLambdaZLine(df, entry.lz ?? null, concName, timeName);
  const yColName = useLogY ? `log(${concName})` : concName;
  try {plot.dataFrame = df;} catch {/* stale scatterPlot hit-test BitSet, as elsewhere */}
  plot.setOptions({
    xColumnName: timeName,
    yColumnName: yColName,
    showDataframeFormulaLines: useLogY,
    showTitle: false,
    description: `${subject}: ${buildFacetVerdictDescriptions(results, thresholds, nonContiguous).get(subject) ?? ''}`,
    descriptionPosition: 'Top' as const,
  });
}

function ensureFacetGrid(
  view: DG.TableView,
  df: DG.DataFrame,
  subjectName: string | null,
  timeName: string | null,
  concName: string | null,
  createIfMissing: boolean,
  useLogY: boolean,
  subjectOrder?: ReadonlyMap<string, number>,
  subjectDescriptions?: ReadonlyMap<string, string>,
  subjectLz?: ReadonlyMap<string, FacetLambdaZParams>,
  subjectColors?: ReadonlyMap<string, number>,
  fitRowsBySubject?: ReadonlyMap<string, ReadonlySet<number>>,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  /** Subject → its bottom-caption verdict text; '' when absent (a degenerate
   *  results frame), which renders as no caption rather than a wrong one. */
  const descOf = (subject: string): string => subjectDescriptions?.get(subject) ?? '';
  const existingDiv = tempBag[FACET_GRID_DIV_TAG] as HTMLDivElement | null;
  const existingNode = tempBag[FACET_GRID_NODE_TAG] as DG.DockNode | null;
  const existingPlots = tempBag[FACET_PLOTS_TAG] as DG.Viewer[] | null;
  const existingSubjects = tempBag[FACET_SUBJECTS_TAG] as ReadonlyArray<string> | null;

  const closeExisting = (): void => {
    if (existingNode != null) {
      try {
        view.dockManager.close(existingNode);
      } catch {/* node may already be detached */}
    } else if (existingDiv != null)
      existingDiv.remove();

    tempBag[FACET_GRID_DIV_TAG] = null;
    tempBag[FACET_GRID_NODE_TAG] = null;
    tempBag[FACET_PLOTS_TAG] = null;
    tempBag[FACET_SUBJECTS_TAG] = null;
  };

  if (subjectName === null || timeName === null || concName === null) {
    if (existingDiv != null) closeExisting();
    return;
  }

  // Skip the rebuild on a Subject/Time/Concentration column collision — see
  // `chartColumnsDistinct`. Leaves any existing grid intact rather than
  // crashing in `buildFacetDataFrames`.
  if (!chartColumnsDistinct(subjectName, timeName, concName)) return;

  if (existingDiv == null && !createIfMissing) return;

  const perSubjectDfs = buildFacetDataFrames(
    df, subjectName, timeName, concName, subjectLz, fitRowsBySubject, lzActionOfView(view));

  // NO upper bound (user-directed 2026-08-13). The facet grid used to refuse
  // above `MAX_FACET_PROFILES` subjects and explain itself with an info toast;
  // both are gone. Past one viewport `buildFacetGrid` switches to the scrolling
  // layout, which keeps the tile readable instead of refusing to draw it.
  // An EMPTY set still returns early — there is nothing to plot, and that is not
  // a suppression to disclose.
  if (perSubjectDfs.length === 0) {
    if (existingDiv != null) closeExisting();
    return;
  }

  // Reorder the tiles to the caller's problem-first key (item 2) so the facet reads
  // in the SAME order as the results grid. Stable sort: subjects absent from the map
  // (a null-subject edge, or an undefined order) fall to the tail keeping their
  // natural order.
  const ordered = subjectOrder === undefined ? perSubjectDfs :
    [...perSubjectDfs].sort((a, b) =>
      (subjectOrder.get(a.subject) ?? Number.POSITIVE_INFINITY) -
      (subjectOrder.get(b.subject) ?? Number.POSITIVE_INFINITY));

  const newSubjects = ordered.map((p) => p.subject);
  const yColName = useLogY ? `log(${concName})` : concName;

  // SAME SUBJECT SET (any order) → reuse the already-docked tiles in place. The
  // check is SET-based, not sequence-based, on purpose: a threshold-recompute keeps
  // the same subjects but the problem-first ORDER can change, and taking the
  // close+re-dock rebuild path on a mere reorder duplicated the panel (the dock is
  // deferred by setTimeout, and closing the prior node did not reliably detach it).
  // The tiles are positional, so rebinding position `i` to `ordered[i]` — data,
  // label, and subject tag — reorders the facet visually with no rebuild.
  const existingSet = existingSubjects != null ? new Set(existingSubjects) : null;
  const sameSubjectSet =
    existingPlots != null && existingSet != null &&
    existingSubjects!.length === newSubjects.length &&
    newSubjects.every((s) => existingSet.has(s));

  if (sameSubjectSet) {
    for (let i = 0; i < existingPlots!.length; i++) {
      const {subject, df} = ordered[i];
      // Colour: ACCEPTANCE-verdict hue when supplied (follows the subject, not the
      // position), else STATUS_NEUTRAL. Paint the fresh frame before it is
      // bound — the colour tags live on the column, and a frame straight from
      // buildFacetDataFrames carries none.
      paintFacetSeries(df, subject,
        subjectColors?.get(subject) ?? STATUS_NEUTRAL);
      // Spike `facet-lambdaz-line-points`: marker coding (M3) + λz line (M4) live on
      // the frame → re-apply on the fresh swapped-in frame (as with the colour tag).
      applyFacetPointMarkers(df);
      applyFacetLambdaZLine(df, ordered[i].lz ?? null, concName, timeName);
      // TEMP: guard ONLY the frame swap — a different-row-count dataFrame leaves
      // DG's scatterPlot pointHitTest BitSet stale (a platform bug); setOptions stays
      // unguarded so real config bugs surface.
      try {existingPlots![i].dataFrame = df;} catch {/* stale hit-test BitSet */}
      existingPlots![i].setOptions({
        xColumnName: timeName,
        yColumnName: yColName,
        // λz line visibility follows the current scale (semi-log only); the line
        // content was just (re)written onto the fresh frame above.
        showDataframeFormulaLines: useLogY,
        // The `subject: reasons` label (now the DESCRIPTION, title off — see create)
        // must follow the reorder, not stay put (setOptions merges, so the position /
        // visibility options set at create persist and need no restating here).
        showTitle: false,
        description: `${subject}: ${descOf(subject)}`,
      });
      // Keep the subject→tile lookup (the selection-dim bridge) correct post-reorder.
      existingPlots![i].root.dataset.ncaFacetSubject = subject;
    }
    tempBag[FACET_SUBJECTS_TAG] = newSubjects;
    wireFacetClickToOverride(view, existingPlots!); // frames were swapped → re-wire
    return;
  }

  if (existingDiv != null) closeExisting();

  const {root: facetRoot, plots} =
    buildFacetGrid(ordered, timeName, yColName, subjectDescriptions, concName, useLogY, subjectColors);
  tempBag[FACET_GRID_DIV_TAG] = facetRoot;
  tempBag[FACET_PLOTS_TAG] = plots;
  tempBag[FACET_SUBJECTS_TAG] = newSubjects;
  wireFacetClickToOverride(view, plots);

  window.setTimeout(() => {
    try {
      // The charts group, not `findNode(profileChart.root)`: this dock runs 100 ms
      // late, after every other panel has rearranged the tree around it, and
      // `dockChartTab` re-picks a member that is still attached at THAT moment.
      tempBag[FACET_GRID_NODE_TAG] = dockChartTab(view, facetRoot, 'Facet');
    } catch {/* dock-manager race; user can still see overlay tab */}
  }, 100);
}

/** Bridge results-grid selection → per-facet-plot opacity. Empty selection
 *  resets all plots to fully opaque; otherwise plots whose subject is not in
 *  the selected set are dimmed to opacity 0.12 (animated via the
 *  `transition: opacity 120ms` set in {@link buildFacetGrid}).
 *
 *  Re-wired on every Compute because `results` is a fresh DataFrame each
 *  time; the prior subscription targets a stale source. */
/**
 * Scroll a facet tile into view inside its own scroll container — and nowhere
 * else.
 *
 * **Deliberately not `Element.scrollIntoView()`.** That walks every scrollable
 * ancestor, so in a docked app it can scroll the pane, the dock area and the
 * page to reveal one tile — moving surfaces the scientist never asked to move.
 * The arithmetic below adjusts exactly one `scrollTop`: the facet grid's.
 *
 * The container is located by CLASS, not by walking up to the first scrollable
 * ancestor. That walk was the first draft and it is wrong in the case that
 * matters most: when the facet grid happens to fit, the walk sails past it and
 * finds the dock pane or the page instead — scrolling the whole app to reveal a
 * tile that was already visible. `closest` can only ever find the facet grid, so
 * in the split layout (no such element) this is a no-op, which is correct: a
 * split layout has nothing off-screen to scroll to.
 */
export function scrollFacetTileIntoView(tile: HTMLElement): void {
  const scroller = tile.closest<HTMLElement>('.nca-facet-scroll');
  if (scroller === null) return;
  const c = scroller.getBoundingClientRect();
  const t = tile.getBoundingClientRect();
  // `block: 'nearest'` semantics: move only if the tile is actually outside, and
  // only far enough to bring it in. A tile already on screen must not jump.
  if (t.top < c.top) scroller.scrollTop += t.top - c.top;
  else if (t.bottom > c.bottom) scroller.scrollTop += t.bottom - c.bottom;
}

function wireResultsSelectionToFacetDim(
  view: DG.TableView,
  results: DG.DataFrame,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  const prevSub = tempBag[FACET_SELECTION_SUB_TAG] as Unsubscribable | null;
  if (prevSub != null) prevSub.unsubscribe();

  const resultsSubjectCol = results.col('SUBJECT');
  if (resultsSubjectCol === null) return;

  const apply = (): void => {
    const plots = tempBag[FACET_PLOTS_TAG] as DG.Viewer[] | null;
    if (plots == null) return;
    const selection = results.selection;
    if (selection.trueCount === 0) {
      for (const p of plots) p.root.style.opacity = '1';
      return;
    }
    const selected = new Set<string>();
    const indexes = selection.getSelectedIndexes();
    for (let i = 0; i < indexes.length; i++)
      selected.add(String(resultsSubjectCol.get(indexes[i])));
    let firstSelected: HTMLElement | null = null;
    for (const p of plots) {
      const subj = p.root.dataset.ncaFacetSubject ?? '';
      const hit = selected.has(subj);
      p.root.style.opacity = hit ? '1' : '0.12';
      if (hit && firstSelected === null) firstSelected = p.root;
    }
    // Dimming alone is a half-answer once the grid scrolls: with 40 subjects the
    // highlighted tile is usually below the fold, so clicking a results row lit up
    // something the scientist could not see (reported 2026-08-13). Bring the first
    // selected tile into view — first in DOM order, i.e. the topmost, so a
    // multi-row selection scrolls to its start rather than an arbitrary member.
    if (firstSelected !== null) scrollFacetTileIntoView(firstSelected);
  };

  tempBag[FACET_SELECTION_SUB_TAG] = results.selection.onChanged.subscribe(apply);
  apply();
}

/** Linked selection: By-group aggregate grid → Mean±SD chart, at GROUP grain
 *  (VR-21). Selecting group rows in the aggregate grid filters the Mean±SD chart
 *  to those groups' curves; an empty selection shows all. The group analogue of
 *  {@link wireResultsSelectionToProfileFilter} (which links per-SUBJECT rows to the
 *  profile chart) — same shape: subscribe to `selection.onChanged`, read the chart
 *  frame at event time (never a stale capture), `trueCount === 0` → show all.
 *
 *  Group identity is reconstructed the SAME way the chart labels its curves: the
 *  per-key group columns carried on `agg` (the user's `keys`, in order) joined by
 *  {@link GROUP_LABEL_SEP}. This is byte-identical to the chart's `Group` values
 *  because the grid and chart share the same per-key labels
 *  (group-by-cl-basis-explicit-control F4/F5: canonical ROUTE / CL·F, numeric-key
 *  `getString` identity). A single key joins to the value itself (no separator). If
 *  the labels ever fail to match (e.g. a key present in results but not the chart
 *  frame), the `Set` simply misses and nothing is filtered — safe degradation.
 *
 *  Re-wired on every {@link applyGrouping} because `agg` is a fresh frame per
 *  Group-by change / Compute; torn down in {@link closeAggregated} when grouping is
 *  cleared. */
function wireByGroupSelectionToMeanSdFilter(
  view: DG.TableView,
  agg: DG.DataFrame,
  keys: ReadonlyArray<string>,
): void {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const tempBag = view.temp as Record<string, any>;
  const prevSub = tempBag[BYGROUP_SELECTION_SUB_TAG] as Unsubscribable | null;
  if (prevSub != null) prevSub.unsubscribe();

  const keyCols = keys
    .map((k) => agg.col(k))
    .filter((c): c is DG.Column => c !== null);
  if (keyCols.length === 0) return; // group columns not carried → cannot match → no link

  const apply = (): void => {
    const chart = tempBag[MEAN_SD_CHART_TAG] as DG.Viewer | null;
    if (chart == null) return;
    const chartDf = chart.dataFrame;
    if (chartDf == null) return;
    const groupCol = chartDf.col('Group');
    if (groupCol === null) {
      chartDf.filter.setAll(true); // no group axis on the chart → nothing to narrow to
      return;
    }
    const selection = agg.selection;
    if (selection.trueCount === 0) {
      chartDf.filter.setAll(true);
      return;
    }
    const selected = new Set<string>();
    const indexes = selection.getSelectedIndexes();
    for (let i = 0; i < indexes.length; i++) {
      const row = indexes[i];
      selected.add(keyCols.map((c) => c.isNone(row) ? '' : c.getString(row)).join(GROUP_LABEL_SEP));
    }
    // Safe degradation (not a blank chart): if NONE of the selected group labels exist
    // among the chart's curve labels, the two frames key on a different dimension set —
    // e.g. a results-only key (STATUS, scalar-ROUTE) that the aggregate carries but the
    // chart's observation frame drops. Filtering would hide every curve; instead show all.
    if (!groupCol.categories.some((cat) => selected.has(cat))) {
      chartDf.filter.setAll(true);
      return;
    }
    chartDf.filter.init((row) => selected.has(String(groupCol.get(row))));
  };

  tempBag[BYGROUP_SELECTION_SUB_TAG] = agg.selection.onChanged.subscribe(apply);
  tempBag[BYGROUP_APPLY_TAG] = apply; // re-appliable after a Mean±SD frame swap (Semi-log)
  apply();
}

/**
 * Compose ONE form row carrying a value control and the unit picker that scales
 * it, side by side (mapping-pane-inline-units, 2026-07-28).
 *
 * The unit pickers used to live in a standalone `Units` section, several rows of
 * scroll away from the values they qualify. Reading `Dose  [1]` there told you
 * nothing about whether that was 1 mg or 1 mg/kg — and that difference is a
 * multiplication by body weight, i.e. a silent CL/Vz mis-scale. Putting the two
 * on one row makes the value unreadable WITHOUT its unit, which is the point.
 *
 * `values` is a list because the Pattern-A groups (Dose, Body weight) carry TWO
 * mutually-exclusive controls — a scalar and a column picker — of which exactly
 * one is visible at a time. Both are handed in and the hidden one is simply a
 * `display: none` flex child, so `linkSourceVisibility` / `setControlVisible`
 * keep operating on the inner roots and need no change at all. Single-control
 * rows (Time, Concentration) pass a one-element list.
 *
 * Measured on the live platform (probe, grok test 2026-07-28) before adopting
 * this shape:
 *  - `ui.divH` takes its flex from the `d4-flex-row` CSS CLASS and sets NO inline
 *    style, so the `style.display = 'none'` → `style.display = ''` round-trip that
 *    {@link setControlVisible} and {@link appendCollapsibleSection}'s `updateWgts`
 *    both perform restores `flex` rather than collapsing the row to `block`.
 *  - `.ui-input-root` does NOT use `display: contents`, and a wrapped row's caption
 *    lands at the same x as an unwrapped row's (delta 0px) — the wrapper does not
 *    break the form's label column.
 */
function unitRow(values: ReadonlyArray<DG.InputBase>, unit: DG.InputBase): HTMLElement {
  // Tag the picker so the width contract in `.nca-unit-row` can single it out:
  // it shows `h` / `mg/L` / `µmol/L`, never a column name, and DG's default
  // 200px editor would eat the row inside a ~346px pane (see app-styles.css).
  unit.root.classList.add('nca-unit-picker');
  return ui.divH([...values.map((v) => v.root), unit.root], 'nca-unit-row');
}

/** Build (but do NOT append) the elements of a sub-section: an `h3` header
 *  followed by the inputs' `.root` HTMLElements. Returned as a flat array
 *  for spreading into a parent collapsible section's `children`. h3 spans
 *  the full grid row via the `.nca-mapping-pane h3` CSS rule. Optional
 *  `tooltip` binds a hover tooltip to the h3 — a plain string or a `richTip()`
 *  factory (same as the parent section's tooltip in
 *  {@link appendCollapsibleSection}). */
function buildSubsection(
  title: string,
  inputs: ReadonlyArray<DG.InputBase>,
  tooltip?: string | (() => HTMLElement),
): HTMLElement[] {
  const h3 = ui.h3(title);
  // Indented to the parent section's TITLE text (16px chevron slot + 8px margin, see
  // appendCollapsibleSection), so AUC / BLQ / λz / … read as children of Rules (user-directed
  // 2026-10-13). Inline, like the section chevrons: the `form h3` CSS rule does not reach a ui.form.
  h3.style.paddingLeft = '24px';
  if (tooltip !== undefined) ui.tooltip.bind(h3, tooltip);
  return [h3, ...inputs.map((i) => i.root)];
}

/** Append a collapsible top-level section. Click anywhere on the header
 *  toggles the section open/closed. Mirrors DiffStudio's
 *  `getCategoryWidget` (utils.ts:665) idiom — chevron-right when closed,
 *  chevron-down when open — but adapted to our flat `ui.form` grid: the
 *  chevrons live INSIDE the `h2` (via `prepend`) so the header keeps the
 *  full-row span from `.nca-mapping-pane form h2`. `children` includes
 *  every element to be hidden when collapsed: sub-section headers, inputs'
 *  roots, AND row wrappers built by {@link unitRow} (four of the Profile-data
 *  / Context children are wrappers, not raw roots). Hide is via
 *  `.hidden = true` (i.e. `display: none`). */
function appendCollapsibleSection(
  form: HTMLElement,
  title: string,
  children: ReadonlyArray<HTMLElement>,
  initialExpanded: boolean,
  tooltip?: string | (() => HTMLElement),
): () => void {
  // Two-chevron pattern from DiffStudio `getCategoryWidget`
  // (DiffStudio/src/utils.ts:665). One chevron is visible at a time;
  // the other is `display: none`. Anti-jump trick: BOTH chevrons are
  // forced to a fixed 16px-wide inline-block slot with `text-align:
  // center`, applied as **inline styles** so they bypass the FA cascade
  // and any class-based specificity battles — the title text always
  // lands at `16 (slot) + 8 (margin) = 24px` from the h2 start,
  // independent of which chevron is currently shown. Inline styles
  // (specificity 1,0,0,0) trump every CSS rule on the element.
  //
  // Children are hidden via `.style.display = 'none'` (NOT `.hidden`)
  // because DG input wrappers carry their own `display:` rule that
  // out-specifies the user-agent `[hidden]` selector — `.hidden = true`
  // is a no-op on those.
  const stylize = (el: HTMLElement): void => {
    el.style.display = 'inline-block';
    el.style.boxSizing = 'border-box';
    el.style.width = '16px';
    el.style.marginRight = '8px';
    el.style.fontSize = '11px';
    el.style.lineHeight = '1';
    el.style.textAlign = 'center';
    el.style.verticalAlign = 'middle';
    el.style.color = 'var(--grey-4)';
    el.style.cursor = 'pointer';
  };

  let expanded = initialExpanded;
  const updateWgts = (e: boolean): void => {
    expanded = e;
    // Naming convention from DiffStudio: `chevronToOpen` = the icon you
    // CLICK to open (shown only while collapsed); `chevronToClose` = the
    // icon you CLICK to close (shown only while expanded). Naming maps
    // to user-intent, not to current-state visibility.
    chevronToOpen.style.display = e ? 'none' : 'inline-block';
    chevronToClose.style.display = e ? 'inline-block' : 'none';
    // On expand, restore each child to visible UNLESS a source-toggle has marked
    // it hidden (`data-src-hidden="1"`, set by `setControlVisible`). Without this
    // guard the blanket `display=''` would clobber per-control Pattern-A
    // visibility (scalar/column hidden by its 'value' / 'column' toggle).
    //
    // FOOTGUN for anyone extending `unitRow`: this loop does NOT recurse into
    // descendants, so for a wrapped row the guard protects the WRAPPER, which
    // nothing ever marks `srcHidden` — the inner Pattern-A visibility survives
    // because it is simply never touched here, not because the guard defends it.
    // A group that needs the WHOLE row hidden (the infusion-duration group is the
    // live candidate — `refreshInfusionVisibility` hides its three controls
    // individually today) must call `setControlVisible` on the wrapper itself, or
    // an expand will restore an empty flex row.
    for (const c of children)
      c.style.display = e ? (c.dataset.srcHidden === '1' ? 'none' : '') : 'none';
  };

  const chevronToOpen = ui.iconFA('chevron-right', () => updateWgts(true), 'Expand');
  const chevronToClose = ui.iconFA('chevron-down', () => updateWgts(false), 'Collapse');
  stylize(chevronToOpen);
  stylize(chevronToClose);

  const header = ui.h2(title);
  header.prepend(chevronToOpen, chevronToClose);
  header.style.cursor = 'pointer';
  if (tooltip !== undefined) ui.tooltip.bind(header, tooltip);
  // Header text click also toggles — only when the click target IS the
  // h2 itself (so the chevrons' own click handlers continue to handle
  // chevron-icon clicks without double-firing).
  header.addEventListener('click', (ev) => {
    if (ev.target === header) updateWgts(!expanded);
  });

  form.append(header);
  for (const c of children) form.append(c);
  updateWgts(initialExpanded);
  // Reapply handle: recomputes every child's display against the CURRENT expand state
  // + each child's `data-src-hidden` flag, without toggling expand/collapse. Callers
  // that flip a child's `data-src-hidden` outside a header click (e.g. the dual-mode
  // unit-picker `applyUnitPickerMode`, which fires from a resize-driven MutationObserver)
  // call this so the flip composes with collapse instead of writing `display` blindly —
  // a raw write would pop a single row back under a collapsed header (decision-audit fix).
  return () => updateWgts(expanded);
}

function buildContract(picked: {
  subject: DG.Column | null; time: DG.Column | null; conc: DG.Column | null;
  doseCol: DG.Column | null; bodyWeightCol: DG.Column | null;
  blqFlagCol: DG.Column | null; lloqCol: DG.Column | null;
  /** Encoded int32 RouteCode column (per-row-route Option A). When present the
   *  contract gains `observations.route`; the extractor reads it per profile. */
  routeCol?: DG.Column | null;
  /** Sparse-design roles (UC-04, AD-11). All optional; populated only when the
   *  Sparse-design toggle is on. The dense pipeline ignores them; their presence
   *  on the contract is what mounts `validateSparseDesign` (the AD-2 pre-compute
   *  topology cross-check). Absent on dense studies ⇒ byte-unchanged. */
  nominalTimeCol?: DG.Column | null;
  animalIdCol?: DG.Column | null;
  samplingTypeCol?: DG.Column | null;
}): ColumnContract {
  if (picked.subject === null || picked.time === null || picked.conc === null)
    throw new Error('buildContract: required column not picked (Compute should have been disabled).');

  const subjectRawType: ColumnRawType =
    picked.subject.type === DGNS.COLUMN_TYPE.STRING ? 'string' : 'int32';
  const animalRawType: ColumnRawType | null = picked.animalIdCol != null ?
    (picked.animalIdCol.type === DGNS.COLUMN_TYPE.STRING ? 'string' : 'int32') : null;

  const contract: ColumnContract = {
    observations: {
      subjectId: {name: picked.subject.name, required: true, rawType: subjectRawType},
      time: {name: picked.time.name, required: true, rawType: 'float64'},
      concentration: {name: picked.conc.name, required: true, rawType: 'float64'},
      ...(picked.doseCol !== null ?
        {dose: {name: picked.doseCol.name, required: false, rawType: 'float64' as const}} :
        {}),
      ...(picked.routeCol !== null && picked.routeCol !== undefined ?
        {route: {name: picked.routeCol.name, required: false, rawType: 'int32' as const}} :
        {}),
      ...(picked.bodyWeightCol !== null ?
        {bodyWeight: {name: picked.bodyWeightCol.name, required: false, rawType: 'float64' as const}} :
        {}),
      ...(picked.blqFlagCol !== null ?
        {blqFlag: {name: picked.blqFlagCol.name, required: false, rawType: 'uint8' as const}} :
        {}),
      ...(picked.lloqCol !== null ?
        {lloq: {name: picked.lloqCol.name, required: false, rawType: 'float64' as const}} :
        {}),
      ...(picked.nominalTimeCol != null ?
        {nominalTime: {name: picked.nominalTimeCol.name, required: false, rawType: 'float64' as const}} :
        {}),
      ...(picked.animalIdCol != null && animalRawType !== null ?
        {animalId: {name: picked.animalIdCol.name, required: false, rawType: animalRawType}} :
        {}),
      ...(picked.samplingTypeCol != null ?
        {samplingType: {name: picked.samplingTypeCol.name, required: false, rawType: 'string' as const}} :
        {}),
    },
  };
  return contract;
}

function renderValidationErrors(panel: HTMLElement, errors: ReadonlyArray<ValidationError>): void {
  panel.innerHTML = '';
  panel.style.display = 'block';
  panel.appendChild(ui.divText(
    // Confirm is the compute trigger since nca-workflow-stages (D17 — there is no Compute button).
    errors.length === 1 ?
      'Validation failed — fix this error and click Confirm again:' :
      `Validation failed — fix these ${errors.length} errors and click Confirm again:`,
    'nca-mapping-pane-error-panel-title',
  ));
  for (const e of errors) {
    const where = e.column !== undefined ?
      (e.row !== undefined ? ` on '${e.column}' row ${e.row}` : ` on '${e.column}'`) :
      '';
    panel.appendChild(ui.divText(`• ${e.code}${where}: ${e.message}`,
      'nca-mapping-pane-error-item'));
  }
  grok.shell.warning(`NCA mapping: ${errors.length} validation error${errors.length === 1 ? '' : 's'}.`);
}

function clearValidationErrors(panel: HTMLElement): void {
  panel.innerHTML = '';
  panel.style.display = 'none';
}
