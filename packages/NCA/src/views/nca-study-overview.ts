/**
 * Project overview — the cached two-tier overview rendered when a
 * project URL (`/apps/Nca/<study label>`) or a list row opens.
 *
 * **Cache-only render.** This view reads ONLY each analysis's
 * `summaryCache` — it issues ZERO `runNcaPipeline` calls. Opening a project is a
 * read; full results recompute only when the user opens an analysis (lazy, AD-4).
 *
 * **Drift mark (AC-D3, decision g).** After the instant cache render, a
 * background source re-probe marks analyses whose source has drifted with a
 * "stale — refreshes on open" indicator. The actual recompute is the lazy
 * on-open path (the spec-mandated single compute path, AD-4/AD-5); a fully
 * headless background recompute is a follow-up (it requires an off-screen view
 * lifecycle — see TODO / audit).
 */

import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {
  buildStudyGrid, buildCurveFrame, buildOverlayCurveFrame,
  PARAM_TAG, STAT_TAG, PASSTHROUGH_TAG, PASSTHROUGH_STRUCT, ANALYSIS_ID_COL, COMPOUND_ID_COL, FLAGS_COL,
  ANALYSIS_COL, ROUTE_COL, CURVE_COL, studyGridRowTiers, normalizeGroupKey, dedupGroupValuesByCanonical,
  finiteMinMax, studyCellExclusion, CROSS_DOSE_POOL_CODE, uniqueLabel, getNcaRepository, isEntryStale,
  studyCellBasisStale, STALE_ENGINE_FLAG, studyCurvePending, pendingCurveIds, withCurves,
} from '../adapter';
import type {NcaStudy, AnalysisEntry, CompoundTierRow, SourcePointer, LoadedStudy} from '../adapter';
import {AGG_PARAM, aggStatSentence, EXCLUDED_TEXT} from './results-view';
import {richTipEl, richTip} from './rich-tooltip';
import {updateNcaStudyPath, ncaDeepLinkPath} from './nca-studies-router';
import {openAnalysis, messageView} from './nca-analysis-open';
import {blueIcon, disabledIcon, td, shortDate, NO_EDIT_TIP, CANNOT_UNDO} from './nca-studies-list';
import {launchAnalysis} from './nca-launch';

/**
 * Render the cached two-tier overview for `studyId`. Returns a tombstone /
 * not-authorized message view on failure rather than throwing (AD-8).
 */
export async function renderStudyOverview(studyId: string): Promise<DG.View> {
  let project: LoadedStudy;
  try {
    // The server shows the study only to whoever may view it.
    project = await getNcaRepository().loadStudy(studyId);
  } catch (e) {
    return messageView('Study unavailable',
      `Study '${studyId}' could not be opened (${e instanceof Error ? e.message : String(e)}). ` +
      'It may have been deleted, or you may not have access to it.');
  }

  // The view's content re-renders in place; `refresh` re-reads the study from
  // storage (picks up new analyses, updated counts/summaries) and re-probes drift.
  const root = ui.divV([], 'nca-study-overview');
  // Fill the view height so the grid below can flex-grow into the remaining space
  // (rather than sit at a fixed height with dead space beneath it).
  root.style.height = '100%';
  const view = DG.View.fromRoot(root);
  view.name = project.label;
  updateNcaStudyPath(view, project.label); // label-keyed URL (router resolves back)

  // Extended-columns profile state (chemist-sar-grid S3): OFF by default. Captured
  // here so it persists across in-view refresh() re-renders (like the analysis
  // view's ribbon toggles persist across Compute).
  let extendedView = false;

  /** Bumped by every render: a curve load answers only the render that started it. */
  let generation = 0;
  const render = (s: LoadedStudy): void => {
    const mine = ++generation;
    ui.empty(root);
    const tierOne = buildTierOne(s, () => {void refresh();});

    const title = ui.h1(s.label);
    title.style.margin = '0';

    // "Analyses" section header carries the Refresh + "+ add analysis" actions.
    const analysesTitle = ui.h2('Analyses');
    analysesTitle.style.margin = '0';
    const refreshIcon = blueIcon('sync', () => {void refresh();}, 'Refresh this study');
    const addIcon = s.access.canEdit ? blueIcon('plus', () => addAnalysisDialog(s), 'Add an analysis to this study') :
      disabledIcon('plus', NO_EDIT_TIP);
    const icons = ui.divH([refreshIcon, addIcon]);
    icons.style.alignItems = 'center';
    icons.style.gap = '12px';
    const analysesHeader = ui.divH([analysesTitle, icons]);
    analysesHeader.style.alignItems = 'center';
    analysesHeader.style.gap = '8px';
    // Extra space between the study title (h1) and the Analyses section (h2).
    analysesHeader.style.marginTop = '20px';

    // "Per-compound NCA results" — the L2 unified grid (study-unified-grid AD-4):
    // one row per (compound × route × dose-norm) case across every analysis,
    // rendered cache-only (computes nothing). Header only when there is a study.
    const gridTitle = ui.h2('Per-compound NCA results');
    gridTitle.style.margin = '0';
    ui.tooltip.bind(gridTitle, richTip({
      title: 'Per-compound NCA results',
      summary:
        'Cached Mean / %CV / n per compound (or per group) across all analyses. ' +
        'Cells flag fragility: n < 4, %AUCextrap > 20%. Computes nothing — read from cache.',
    }));

    const {root: gridEl, grid} = buildUnifiedGrid(s, extendedView);
    // "Full" ribbon-style toggle (chemist-sar-grid S3), OFF by default:
    // compact = the SAR lead trio (AUC∞ · MRT · %F) + CL · t½, mean only (Discovery
    // DMPK fast scan); ON reveals the full Mean|%CV|n set + secondary params
    // (ADME/PKPD Modeler). Only offered when there is a grid to profile.
    // One short native tooltip via the input option (shows on the label + the
    // switch). User-directed 2026-07-07 (rename Extended→Full, terse tip).
    const extendedToggle = ui.input.toggle('Full', {
      value: extendedView,
      tooltipText: 'Show the full parameter table.',
    });
    extendedToggle.onChanged.subscribe(() => {
      extendedView = extendedToggle.value ?? false;
      if (grid === null) return;
      applyStudyColumnProfile(grid, extendedView); // also re-stamps STUDY_EXTENDED_TAG
      // If a compound-detail panel is already open, rebuild its "Forms" viewer so its fields
      // track the new profile too. Guarded on an EXISTING panel + active rows so flipping the
      // toggle never force-opens the context panel (renderParamsForm touches no `grok.shell.o`).
      const df = grid.dataFrame;
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      const panel = (df.temp as Record<string, any>)[CTPROFILE_PANEL_TAG] as CompoundPanel | null;
      if (panel != null && activeRowRefs(df).length > 0) renderParamsForm(panel, df);
    });
    // Title left, the "Full" view-toggle pushed to the right edge of the grid
    // (space-between) so it reads as a table view-control, not part of the title
    // (user-directed layout, 2026-07-07).
    const gridHeader = ui.divH(grid !== null ? [gridTitle, extendedToggle.root] : [gridTitle]);
    gridHeader.style.alignItems = 'center';
    gridHeader.style.justifyContent = 'space-between';
    gridHeader.style.marginTop = '24px';

    // Content column fills the view height; the fixed-height rows (title / Analyses
    // / grid header) take their natural size and `gridEl` flex-grows into the rest.
    const content = ui.divV([
      title,
      analysesHeader,
      tierOne,
      gridHeader,
      gridEl,
    ]);
    content.style.height = '100%';
    content.style.minHeight = '0'; // allow the flex child (grid) to shrink + scroll
    root.appendChild(content);
    // AC-D3 / F-N.6 stale mark: a stored snapshot cannot drift (SD-4), so only an engine change
    // can make cached numbers stale — a synchronous check, no source probe. The recompute is the
    // lazy on-open path.
    markStale(s, tierOne);
    // ems-modification step 5: the result rows came without their curves; load the pending ones
    // (one request) and re-render the thumbnails and the curve panel through this same path.
    const ids = pendingCurveIds(s);
    if (ids.length > 0) {
      getNcaRepository().loadCurves(ids)
        .then((curves) => {if (mine === generation) render(withCurves(s, curves));})
        .catch((e) => grok.shell.error('Could not load the curves: ' + (e instanceof Error ? e.message : String(e))));
    }
  };

  const refresh = async (): Promise<void> => {
    try {
      const fresh = await getNcaRepository().loadStudy(studyId);
      view.name = fresh.label;
      render(fresh);
    } catch (e) {
      grok.shell.error(
        `Could not refresh study: ${e instanceof Error ? e.message : String(e)}`);
    }
  };

  render(project);
  return view;
}

/**
 * "+ add analysis" from a study: a Hub-style dialog (table only) that opens the
 * analysis (mapping pane) on the chosen table, PRE-TARGETED to `study` — so when
 * the user maps, computes, and clicks "Add to study", this study is pre-selected.
 * (Mapping + compute must happen in the analysis view, so this launches it rather
 * than persisting directly.) No name input — the view is titled with the table name
 * and Add to study proposes the analysis name (user-directed 2026-10-06).
 */
function addAnalysisDialog(study: NcaStudy): void {
  const tableInput = ui.input.table('Table', {
    tooltipText: 'Table with subject IDs, time points, and concentration values.',
  });
  ui.dialog('New analysis')
    .add(tableInput.root)
    .onOK(() => {
      const df = tableInput.value;
      if (df === null) {
        grok.shell.warning('Pick a table first.');
        return;
      }
      // DG.DataFrame's runtime `.id` is absent from the platform typings; an
      // unsaved in-session table has an empty id, so we fall back to df.name as
      // the source ref below (the resolver matches on id OR name).
      // eslint-disable-next-line @typescript-eslint/no-explicit-any
      const tableId = (df as any).id as string | undefined;
      const source: SourcePointer = {
        kind: 'tableId', ref: (tableId != null && tableId.length > 0) ? tableId : df.name,
      };
      launchAnalysis(df, df.name, source, study.id);
    })
    .show();
}

// ── Unified per-compound NCA grid (study-unified-grid AD-4) ─────────────────

/** Mean-column number formats by display-parameter. Value params use `'0.000'`
 *  (3 decimals — the codebase's documented proxy for "3 significant digits",
 *  matching the per-subject results grid's `NUMBER_FORMATS` in results-view.ts,
 *  so a parameter formats IDENTICALLY in both grids). Percentages (%CV,
 *  %AUCextrap) use `'0.0'` (1 decimal — 3 would be false precision); n is an
 *  integer (no format). */
const MEAN_FORMAT: Readonly<Record<string, string>> = {
  'Cmax': '0.000', 'AUC∞': '0.000', 't½': '0.000', 'Tmax': '0.000',
  'CL': '0.000', 'CL/F': '0.000', 'Vz': '0.000', 'Vz/F': '0.000',
  'MRT': '0.000', 'Vss': '0.000', '%AUCextrap': '0.0',
  // %F is a percentage (like %AUCextrap) → 1 decimal; 3 would be false precision.
  '%F': '0.0',
};

// ── Header-tooltip glossary bridge (study-grid-header-tooltips spike) ─────────
// The study grid tags parameters by DISPLAY label (`Cmax`, `AUC∞`, …); the shared
// {@link AGG_PARAM} gloss (results-view.ts) is keyed by INTERNAL name (`CMAX`,
// `AUC_INF`, …). This map bridges the two so a `Mean Cmax` tooltip here reuses the
// SAME gloss DEFINITION as the analysis-view aggregated grid, without a parallel
// gloss table (Rule 5). The friendly token still renders in THIS grid's own label
// convention (`Cmax`, not the aggregated grid's `C_max`) so the tooltip matches its
// own header. Covers the route-uniform params BASIS_GLOSS does NOT own — CL/CL·F,
// Vz/Vz·F and Vss carry their per-column basis meaning via BASIS_GLOSS and never
// reach here. MRT is route-uniform (one column for all routes) → bridged here.
const STUDY_PARAM_TO_AGG: Readonly<Record<string, string>> = {
  'Cmax': 'CMAX', 'AUC∞': 'AUC_INF', 't½': 'HALF_LIFE',
  'MRT': 'MRT', 'Tmax': 'TMAX', '%AUCextrap': 'PCT_EXTRAP',
};

// Per-column basis gloss (design B3): the basis-ROUTED pairs CL vs CL/F and Vz vs
// Vz/F are already SEPARATE columns, so each header names its own basis. `Vss` is
// included here too (v3) — it has NO `/F` counterpart (an apparent steady-state
// volume is not defined; it is IV-only, `study-grid.ts`), but it reuses this table
// for its intravascular-basis gloss rather than the generic AGG bridge. The analysis
// view has no per-column equivalent (it carries basis via the aggregator `label` tag
// + the CL_BASIS stratification note), so this true-vs-apparent text is study-grid-
// local — essential domain semantics, not a duplicated gloss (CL_BASIS meaning, TODO note b).
const BASIS_GLOSS: Readonly<Record<string, string>> = {
  'CL': 'true clearance (Dose/AUC∞), intravascular basis (F = 1)',
  'CL/F': 'apparent clearance (Dose/AUC∞), extravascular basis (divided by unknown F)',
  'Vz': 'terminal-phase volume of distribution, intravascular basis (F = 1)',
  'Vz/F': 'apparent terminal-phase volume, extravascular basis (divided by unknown F)',
  'Vss': 'volume of distribution at steady state (MRT·CL), intravascular basis (F = 1)',
  // %F (v4) — cross-profile, no AGG_PARAM entry (F_PCT is not a within-group
  // aggregate), so like Vss it carries its study-grid-local gloss here rather than
  // via the generic AGG bridge. Sex-stratified IV reference; uncapped; blank/NE
  // where a compound has no paired IV arm.
  '%F': 'absolute bioavailability, mean of the EV dose-normalized AUC∞ ÷ the ' +
    'sex-matched IV reference; uncapped, blank where a compound has no IV arm',
};

// Leading (non-parameter) columns — study-grid-specific meaning with no analysis-
// view equivalent (`Analysis` provenance, `Route`, `C–t`, `Warnings`). The per-key
// grouping columns are NOT here — they get their own header (the key's friendly label)
// and carry no fixed tooltip. Mirrors AGG_GROUP_TOOLTIPS in shape (results-view.ts).
const LEADING_TOOLTIPS: Readonly<Record<string, {title: string; desc: string; note?: string}>> = {
  [ANALYSIS_COL]: {title: 'Analysis',
    desc: 'The analysis this cached row came from (provenance) — one study spans many analyses.'},
  [ROUTE_COL]: {title: 'Route', desc: 'Route of administration for this case.'},
  [CURVE_COL]: {title: 'C–t',
    desc: 'Inline semi-log geometric-mean concentration-time thumbnail — one pooled line ' +
      'per case (cross-dose pools are disclosed). Click the row for the full semi-log profile ' +
      'with the ×/÷ gSD band; select several rows to overlay their curves.'},
  [FLAGS_COL]: {title: 'Warnings',
    desc: 'Per-case warnings and not-estimable flags — shown in the cell (and carried ' +
      'into a saved table); hover a cell for the full text.',
    note: 'NE = not estimable.'},
};

/** The gloss definition for a study-grid parameter (display label). BASIS_GLOSS for
 *  the basis-routed columns (CL/CL·F, Vz/Vz·F) plus Vss, else the shared
 *  {@link AGG_PARAM} gloss via {@link STUDY_PARAM_TO_AGG}; `null` when the label is
 *  unknown to both. Exported +
 *  pure so a unit test can pin that EVERY `buildStudyGrid`-emitted `PARAM_TAG`
 *  resolves — the silent-map-drift guard (mirrors the exported `aggDisplayNames`
 *  contract, GAP-10; study-grid-header-tooltips review). */
export function studyParamGloss(param: string): string | null {
  return BASIS_GLOSS[param] ?? AGG_PARAM[STUDY_PARAM_TO_AGG[param] ?? '']?.gloss ?? null;
}

/** The stat-sentence for a study-grid stat tag. `n` uses a study-grid-local sentence
 *  (this grid has no cohort-N column, so the aggregated grid's "≤ the group's N"
 *  wording would cite a nonexistent column); `Mean` / `%CV` reuse the aggregated
 *  grid's shared wording. `null` for an unknown stat. Exported + pure for the same
 *  drift guard as {@link studyParamGloss}. */
export function studyStatSentence(stat: string, param: string): string | null {
  if (stat === 'n')
    return 'Number of subjects contributing a non-null value to this parameter.';
  return aggStatSentence(stat, param);
}

// ── Chemist column profile (chemist-sar-grid S3) ─────────────────────────────
// The grid leads with the DMPK optimization signal. COMPACT (default, Discovery
// DMPK) shows the SAR lead trio + CL/t½ as MEAN-ONLY columns for a fast analog-
// series scan; EXTENDED (ADME/PKPD Modeler) reveals the full Mean|%CV|n set + the
// secondary params. Column VISIBILITY only — the values are byte-identical; a
// hidden column still sorts once revealed (design B4 / thesis #2 DMPK-first).

/** Compact param order (left→right): the SAR lead trio (AUC∞ · MRT · %F) then
 *  CL · t½. Only each param's Mean column is shown in compact. CL/CL·F and (later)
 *  Vz basis-routed pairs share a slot — whichever basis is present on the screen. */
const COMPACT_PARAMS: ReadonlyArray<string> = ['AUC∞', 'MRT', '%F', 'CL', 'CL/F', 't½'];
/** Extended param order: lead trio + CL/t½ first (so toggling keeps the eye in
 *  place), then the secondary zone in the spec's order `Cmax · Vss · Vz ·
 *  %AUCextrap` (Vz/F rides with Vz — the basis-routed EV counterpart). */
const EXTENDED_PARAMS: ReadonlyArray<string> = [
  'AUC∞', 'MRT', '%F', 'CL', 'CL/F', 't½', 'Cmax', 'Vss', 'Vz', 'Vz/F', '%AUCextrap',
];
/** Params projected as a single standalone column (no Mean|%CV|n triple, design
 *  B3) — their column IS the bare display label. */
const STANDALONE_PARAMS: ReadonlySet<string> = new Set(['%AUCextrap']);

/** The fixed (non-param, non-per-key) study-grid columns — used to isolate the
 *  per-key grouping columns (anything present that is neither one of these nor a
 *  `PARAM_TAG`-carrying param column). */
const KNOWN_STUDY_COLS: ReadonlySet<string> = new Set([
  ANALYSIS_COL, ROUTE_COL, CURVE_COL, FLAGS_COL, ANALYSIS_ID_COL, COMPOUND_ID_COL]);

/** Ordered list of VISIBLE column names for the chemist profile. Leading context
 *  (Analysis / per-key grouping columns / Route / C–t) always shows; then the param
 *  columns per the active profile; then Warnings (Full profile only). Only names actually
 *  present in `df` are included (a param group the study dropped as all-empty — e.g. %F on
 *  an all-IV screen — is simply absent; a study where no analysis groups → no key columns). */
/** Temp-bag tag: the study grid's active "Full" profile (`true` = extended). Stamped on the grid
 *  frame by {@link applyStudyColumnProfile} so the context panel's frame-bound "Forms" viewer can
 *  read the SAME profile the grid shows without threading the render-closure local. */
export const STUDY_EXTENDED_TAG = 'ncaStudyExtended';

/** The column names the "NCA parameters" Forms viewer should render — the grid's currently-visible
 *  columns ({@link studyGridVisibleOrder} for the active {@link STUDY_EXTENDED_TAG} profile) MINUS
 *  the C–t thumbnail (a `DATA_FRAME` column that has no form field) — so the pane mirrors the grid:
 *  compact when the "Full" toggle is off, the full Mean|%CV|n set when on. The back-pointer columns
 *  are already absent from the visible order. */
export function studyParamsFormFields(df: DG.DataFrame): string[] {
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const extended = ((df.temp as Record<string, any>)[STUDY_EXTENDED_TAG] as boolean | null) === true;
  return studyGridVisibleOrder(df, extended).filter((n) => n !== CURVE_COL);
}

export function studyGridVisibleOrder(df: DG.DataFrame, extended: boolean): string[] {
  const present = new Set(df.columns.names());
  const order: string[] = [];
  const push = (name: string): void => {if (present.has(name)) order.push(name);};
  // Leading context (Analysis / per-key grouping columns / Route) then the inline C–t
  // thumbnail — per the ratified column order `Group(s) · C–t · AUC∞ · …` (S4). The
  // thumbnail shows in BOTH compact and extended profiles (it is the "C–t" lead, not a param).
  push(ANALYSIS_COL);
  // Per-key grouping columns (study-grid-grouping-columns): always one named column per
  // grouping key (the generic `Group` column is retired). They are the present columns that
  // are neither a known leading/trailing column nor a param column (a param carries
  // PARAM_TAG; a per-key column does not), in build (first-appearance) order — each a
  // first-class filterable/sortable column. Zero when no analysis groups.
  for (const name of df.columns.names()) {
    if (KNOWN_STUDY_COLS.has(name)) continue;
    if ((df.col(name)?.getTag(PARAM_TAG) ?? '') !== '') continue;
    // v10: passthrough (annotation) columns are NOT grouping keys — they go in their own zone
    // below (PP-3: without this discriminator they would leak into the key zone).
    if ((df.col(name)?.getTag(PASSTHROUGH_TAG) ?? '') !== '') continue;
    order.push(name);
  }
  // v10 passthrough zone (study-grid-passthrough-columns): the carried annotation columns, AFTER
  // the grouping keys and BEFORE Route (structure-first, as buildStudyGrid emitted them — user-
  // approved 2026-07-14: "structure beside the number"). Discriminated by PASSTHROUGH_TAG.
  for (const name of df.columns.names())
    if ((df.col(name)?.getTag(PASSTHROUGH_TAG) ?? '') !== '') order.push(name);
  push(ROUTE_COL); push(CURVE_COL);
  const params = extended ? EXTENDED_PARAMS : COMPACT_PARAMS;
  for (const p of params) {
    if (STANDALONE_PARAMS.has(p)) {push(p); continue;}
    push(`${p} Mean`);
    if (extended) {push(`${p} %CV`); push(`${p} n`);}
  }
  // Warnings only in the Full profile — the compact scan leads with the SAR signal.
  // EXPLICIT user decision 2026-07-07 that OVERRIDES ratified AC-3c/R1-F4 ("an unfit
  // case shows NE, never an ambiguous blank"): the compact/default view no longer
  // surfaces the NE / not-estimable flag, so a NaN AUC∞ reads as a blank cell there.
  // Accepted tradeoff: a decluttered default view over the NE signal there.
  if (extended) push(FLAGS_COL);
  return order;
}

/**
 * Apply the compact/extended column profile to the grid: order the visible
 * columns per {@link studyGridVisibleOrder}, hide everything else (including the
 * always-hidden back-pointer + warnings source). The frame is unchanged, so
 * sorting a hidden-then-revealed column still works (chemist-sar-grid S3 DoD).
 */
function applyStudyColumnProfile(grid: DG.Grid, extended: boolean): void {
  const df = grid.dataFrame;
  // Stamp the active profile on the grid frame so the context panel's "Forms" viewer (bound to
  // this SAME frame) can mirror the grid's visible columns — the single write-point, since this
  // function is called on both the initial build and every "Full" toggle flip.
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  (df.temp as Record<string, any>)[STUDY_EXTENDED_TAG] = extended;
  const visible = studyGridVisibleOrder(df, extended);
  const visibleSet = new Set(visible);
  // Full order = visible (in profile order) then the rest (hidden), so setOrder
  // gets every column name (DG orders the listed columns; trailing hidden ones
  // keep the frame complete).
  const hidden = df.columns.names().filter((n) => !visibleSet.has(n));
  grid.columns.setOrder([...visible, ...hidden]);
  for (const col of df.columns) {
    const gc = grid.col(col.name);
    if (gc === null) continue;
    // The back-pointer + full-warnings source are NEVER user-visible.
    if (col.name === ANALYSIS_ID_COL) {gc.visible = false; continue;}
    gc.visible = visibleSet.has(col.name);
  }
}

/**
 * The cache-only unified grid section. One native `DG.Grid` over
 * {@link buildStudyGrid}; empty-state caption when no analysis has computed
 * per-compound rows yet. Sized to fit the overview viewport (15" laptop budget)
 * — a fixed-height box the grid fills. Returns the section element plus
 * the grid (null on the empty state) so the caller can wire the Extended toggle
 * to {@link applyStudyColumnProfile} (chemist-sar-grid S3).
 */
function buildUnifiedGrid(study: LoadedStudy, extended: boolean): {root: HTMLElement; grid: DG.Grid | null} {
  // Engine-stale analyses (v14 F-N.6) are marked through ONE predicate (`isEntryStale`, SD-9),
  // shared by BOTH readers of the cache: the grid row (Warnings cell) and the context panel's
  // cached curve (disclosure bar). Two readers, one definition — they cannot disagree.
  const isStale: StalePredicate = (a) => isEntryStale(a);
  const df = buildStudyGrid(study, isStale);
  if (df.rowCount === 0) {
    const empty = ui.divText('No computed results yet — open an analysis to compute it.');
    empty.style.color = 'var(--grey-4, #9497A0)';
    empty.style.marginTop = '12px';
    return {root: empty, grid: null};
  }
  const grid = DG.Viewer.grid(df);
  configureStudyGrid(grid, study, isStale);
  applyStudyColumnProfile(grid, extended);
  const box = ui.div([grid.root], 'nca-study-grid');
  // Flex-grow into the remaining view height (user-directed 2026-07-07 — was a fixed
  // 360px that left dead space below). `min-height` floors it on a short viewport;
  // `min-height: 0` on the flex box lets the grid scroll internally instead of
  // overflowing the page.
  box.style.flex = '1';
  box.style.minHeight = '240px';
  // Tight gap between the section title and the grid (user-directed 2026-07-07 —
  // match the grid's own header→row spacing, not the looser 12px).
  box.style.marginTop = '4px';
  grid.root.style.width = '100%';
  grid.root.style.height = '100%';
  return {root: box, grid};
}

/**
 * Grid-side presentation for the unified grid: read-only, hidden back-pointer +
 * warnings source, parameter-qualified headers (`"Mean Cmax"` etc. — the app's
 * own aggregated-grid convention; native spanning column groups do NOT render
 * in this embedded standalone grid, see below), per-column number formats, and
 * **format-not-color** fragility flags: `n < 4` and `%AUCextrap > 20%` render in
 * **bold** (text emphasis, grayscale-surviving — design B4 forbids a value
 * heatmap / value color; audit-checklist C-05/C-27). Flags cell = NE only; the
 * full warning stream is on its tooltip (B4: no Warnings column in v1). A
 * structure column, if present, renders as a molecule (AC-4b — no such column
 * ships in v1; the guard is forward-compatible for the SMILES slice).
 */
/**
 * Set the parameter-qualified headers + number formats on a study-grid frame.
 * Exported + grid-free so the header contract is unit-testable (the "undefined
 * undefined" regression: `getTag` returns **`undefined`** — NOT null / '' — for
 * an unset tag, so the untagged leading columns MUST be skipped, else
 * `${stat} ${param}` writes the literal "undefined undefined" as their header).
 *
 * Headers carry BOTH stat and parameter (`"Mean Cmax"` / `"%CV Cmax"` /
 * `"n Cmax"`) — the app's own aggregated "By group" grid convention
 * (results-view.ts `configureAggregatedGrid`, which renders correctly). Each
 * parameter's three columns are emitted adjacently by buildStudyGrid, so they
 * read as a visual group. Native spanning column groups (`df.temp['.columnGroups']`)
 * do NOT render in the embedded standalone `DG.Viewer.grid` used here (verified
 * against the live render; Admetica's spanning groups work only on a docked
 * TableView grid) — so the header text must be self-describing.
 */
export function applyStudyGridHeaders(df: DG.DataFrame): void {
  for (const col of df.columns) {
    const param = col.getTag(PARAM_TAG);
    const stat = col.getTag(STAT_TAG);
    if (param === null || param === undefined || param === '') continue;
    if (stat !== null && stat !== undefined && stat !== '') {
      col.setTag('friendlyName', `${stat} ${param}`); // e.g. "Mean Cmax"
      if (stat === 'Mean') col.meta.format = MEAN_FORMAT[param] ?? '0.000';
      else if (stat === '%CV') col.meta.format = '0.0';
    } else
      col.meta.format = MEAN_FORMAT[param] ?? '0.0'; // standalone (%AUCextrap)
  }
}

/** In-cell curve thumbnail options. Native `plot.line`
 *  (standard DG viewers only, NO custom SVG/canvas renderer): one line per dose `Series`, semi-log Y, all chrome
 *  off (it is a sparkline — the drill-in big graph carries the axes/band/legend). */
const THUMB_OPTS: Record<string, unknown> = {
  xColumnName: 'Time',
  yColumnNames: ['Conc'],
  splitColumnName: 'Series',
  yAxisType: 'logarithmic',
  showXAxis: false, showYAxis: false,
  showVerticalGridLines: false, showHorizontalGridLines: false,
  showXSelector: false, showYSelectors: false,
  showTopPanel: false, showMouseOverRowLine: false,
  showAggrSelectors: false, legendVisibility: 'Never',
};

/** Per-cell axis ranges pinned to THIS curve's own extent so the sparkline fills
 *  the cell (auto-range otherwise leaves it flat + top-pinned, and jams the end
 *  points against the left/right edges). Y is semi-log → LOG-space padding, GENEROUS
 *  at the bottom so the lowest marker clears the edge; X is linear time → additive
 *  padding both sides so the first/last markers aren't clipped. Omits an axis whose
 *  range is degenerate (single point) → auto-range there. User-directed 2026-07-07. */
function thumbRanges(frame: DG.DataFrame): Record<string, number> {
  const out: Record<string, number> = {};
  // `finiteMinMax` (adapter-utils) is the shared isNone-FIRST scan — reused here
  // rather than respelled. The bare `conc.get(i)` + `Number.isFinite(v) && v > 0`
  // this replaced was a live trap once `Conc` began carrying the platform null:
  // `Column.get()` returns the raw sentinel even on a null cell (probe PR-2), and
  // `DG.FLOAT_NULL` = 2.68e-34 is FINITE and POSITIVE, so it passed BOTH guards,
  // became `mn`, and `Math.log10(mn) ≈ −33.57` stretched this semi-log axis about
  // 40 decades — flattening every thumbnail into a line at the top of its cell.
  // `positiveOnly` keeps the log axis's own requirement (a real 0 has no
  // position there) without re-admitting the sentinel.
  const cRange = finiteMinMax(frame.col('Conc'), true);
  if (cRange !== null && cRange.max > cRange.min) {
    const lo = Math.log10(cRange.min); const hi = Math.log10(cRange.max); const range = hi - lo;
    out.yMin = Math.pow(10, lo - Math.max(range * 0.4, 0.3)); // generous bottom → no clip
    out.yMax = Math.pow(10, hi + Math.max(range * 0.12, 0.1));
  }
  const tRange = finiteMinMax(frame.col('Time'));
  if (tRange !== null && tRange.max > tRange.min) {
    const pad = (tRange.max - tRange.min) * 0.08; // breathing room so end points aren't jammed
    out.xMin = tRange.min - pad; out.xMax = tRange.max + pad;
  }
  return out;
}

/** Grid row height (px) that gives the inline thumbnail enough vertical room to
 *  read as a curve at the 15" viewport (the grid sits in
 *  a fixed 360px box that scrolls, so taller rows never overflow the page). */
const THUMB_ROW_HEIGHT = 48;
/** Thumbnail column width (px) — a compact sparkline slot. */
const THUMB_COL_WIDTH = 120;

/** Engine-staleness predicate threaded from the coordinator (F-N.6) so the grid row and
 *  the context panel mark the SAME cached analyses. Defaults to never-stale for callers
 *  (and tests) that build a grid without a coordinator. */
type StalePredicate = (entry: AnalysisEntry) => boolean;

function configureStudyGrid(
  grid: DG.Grid, study: NcaStudy, isStale: StalePredicate = () => false,
): void {
  grid.setOptions({allowEdit: false});
  const df = grid.dataFrame;

  // Hide the back-pointer id columns (analysis id + composite-identity compoundId).
  const idGc = grid.col(ANALYSIS_ID_COL);
  if (idGc !== null) idGc.visible = false;
  const cidGc = grid.col(COMPOUND_ID_COL);
  if (cidGc !== null) cidGc.visible = false;

  // v5 inline thumbnail column: HTML cell + widened row/column so the native
  // `plot.line` in `onCellPrepare` has room. Present only when a case carries a
  // curve (buildStudyGrid drops the column otherwise, hide-all-empty parity).
  const curveGc = grid.col(CURVE_COL);
  if (curveGc !== null) {
    curveGc.cellType = 'html';
    curveGc.width = THUMB_COL_WIDTH;
    grid.props.rowHeight = THUMB_ROW_HEIGHT;
  }

  applyStudyGridHeaders(df);

  // Fragility = bold text (format/glyph, NOT value color — design B4). Bold
  // survives grayscale and is co-located with the value. n<4 → fragile fit;
  // %AUCextrap>20 → unreliable extrapolation. Native cell-style hook (mirrors
  // results-view onCellPrepare); no custom cell renderer.
  grid.onCellPrepare((gc) => {
    if (!gc.isTableCell || gc.tableColumn === null) return;
    // v5 inline thumbnail: render the compound's semi-log geomean curve in-cell
    // with the native `plot.line` (no custom renderer). The cell value is the
    // per-compound DataFrame buildStudyGrid stashed (or null on a case with no curve).
    if (gc.tableColumn.name === CURVE_COL) {
      // A curve still being loaded is not "no curve": a muted ellipsis until it arrives.
      // DESIGN-REVIEW: no table match for a transient loading cell — the most conservative
      // fallback (neutral grey text, no hue), the register of the overview's empty caption.
      if (studyCurvePending(df, gc.cell.rowIndex)) {
        const wait = ui.divText('…');
        wait.style.color = 'var(--grey-4, #9497A0)';
        gc.style.element = wait;
        return;
      }
      const frame = gc.cell.value as DG.DataFrame | null;
      if (frame !== null && frame !== undefined && frame.rowCount > 0) {
        // Per-cell axis ranges so the curve fills the cell (not flat + top-pinned)
        // with breathing room so end markers aren't clipped at the edges.
        const el = frame.plot.line({...THUMB_OPTS, ...thumbRanges(frame)}).root;
        el.style.width = '100%';
        el.style.height = '100%';
        gc.style.element = el;
      }
      return;
    }
    // F-S5.2: a λz-dependent Mean / %CV withheld because EVERY contributor was excluded by review
    // reads `not reportable` (muted) — never the blank of a not-estimated value. Same wording and
    // grey as the analysis view's By-group grid.
    if (studyCellExclusion(df, gc.tableColumn.name, gc.cell.rowIndex)?.notReportable === true) {
      gc.customText = 'not reportable';
      gc.style.textColor = EXCLUDED_TEXT;
      return;
    }
    const v = gc.cell.value;
    if (typeof v !== 'number') return;
    const param = gc.tableColumn.getTag(PARAM_TAG);
    const stat = gc.tableColumn.getTag(STAT_TAG);
    const fragile = (stat === 'n' && v < 4) || (param === '%AUCextrap' && v > 20);
    if (!fragile) return;
    const base = gc.style.font;
    gc.style.font = (base !== null && base !== undefined && base.length > 0) ?
      (/bold/i.test(base) ? base : `bold ${base}`) : 'bold 13px sans-serif';
  });

  // v10 (study-grid-passthrough-columns): a passthrough column persisted as a STRUCTURE renders
  // as a 2D molecule. Driven by PASSTHROUGH_TAG === 'structure' (the persisted isStructure
  // designation captured at the add-to-study dialog), SUPERSEDING the retired literal
  // 'Structure'/'SMILES' name-match — so a `Canonical_SMILES`/`Parent_Structure` column renders
  // too, and a text column literally named 'SMILES' does NOT (P3: designation, not name).
  for (const col of df.columns)
    if (col.getTag(PASSTHROUGH_TAG) === PASSTHROUGH_STRUCT) col.semType = DG.SEMTYPE.MOLECULE;

  // v10 (study-grid-passthrough-columns, RG-B): a many-structure passthrough cell shows a bare '—'
  // (a MOLECULE-semType cell can't host a text count); the distinct COUNT is disclosed here on hover
  // (user-approved 2026-07-14). The row's tier is resolved via the SHARED studyGridRowTiers order
  // (same function buildStudyGrid used → grid row index maps 1:1), and the entry by the column's
  // canonical key. `— (N)`-in-cell was rejected: it would deviate from the approved tooltip form.
  const rowTiers = studyGridRowTiers(study);

  // Tooltips: (a) header cells explain each parameter/fragility rule; (b) a
  // Flags data-cell surfaces the FULL warning stream (kept out of the visible
  // cell, which shows only NE — design B4: no Warnings column in v1).
  grid.onCellTooltip((cell, x, y) => {
    // (0) Many-structure passthrough cell → disclose the distinct structure count on hover.
    if (cell.isTableCell && cell.tableColumn?.getTag(PASSTHROUGH_TAG) === PASSTHROUGH_STRUCT) {
      const rowIdx = cell.tableRowIndex;
      if (rowIdx === null) return false;
      const tier = rowTiers[rowIdx]?.r ?? null;
      if (tier === null) return false;
      const canonical = normalizeGroupKey(cell.tableColumn.name);
      const entry = tier.passthrough.find((p) => normalizeGroupKey(p.key) === canonical);
      if (entry !== undefined && entry.values.length > 1) {
        ui.tooltip.show(
          ui.divText(`${entry.values.length} distinct structures in this row`,
            {style: {maxWidth: '360px', whiteSpace: 'normal'}}), x, y);
        return true;
      }
      return false;
    }
    // (a0) A CL / Vz / Vss cell blanked because the analysis was cached under an older reporting
    // basis (nca-cl-bw-normalised-units AD-8, architect B3): say why it is empty and when it fills.
    if (cell.isTableCell && cell.tableColumn !== null && cell.tableRowIndex !== null &&
        studyCellBasisStale(df, cell.tableColumn.name, cell.tableRowIndex)) {
      ui.tooltip.show(richTipEl({
        title: STALE_ENGINE_FLAG,
        summary: 'These results were saved before CL, Vz and Vss were reported per kg (mL/min/kg, L/kg), ' +
          'so the old values are not shown under the per-kg header. They are recomputed when someone ' +
          'with edit access opens the analysis.',
      }), x, y);
      return true;
    }
    // (a') A parameter cell of a row with Exclude λz decisions → `n = a of b` (F-S5.2).
    if (cell.isTableCell && cell.tableColumn !== null) {
      const rowIdx = cell.tableRowIndex;
      const ex = rowIdx === null ? null : studyCellExclusion(df, cell.tableColumn.name, rowIdx);
      if (ex !== null) {
        ui.tooltip.show(richTipEl({
          title: ex.notReportable ? 'not reportable' : `n = ${ex.n} of ${ex.total}`,
          summary: ex.notReportable ?
            `Every profile's λz was excluded by review (n = 0 of ${ex.total}).` :
            `${ex.total - ex.n} of ${ex.total} profiles excluded by review (Exclude λz) — not in this statistic.`,
        }), x, y);
        return true;
      }
    }
    // (b) Flags data cell → the cell's own text (no title), so a truncated
    // multi-warning cell is fully readable on hover. The text lives in the cell
    // (exportable); the tooltip just reflects it (user-directed 2026-07-07).
    if (cell.isTableCell && cell.tableColumn?.name === FLAGS_COL) {
      const rowIdx = cell.tableRowIndex;
      if (rowIdx === null) return false;
      const text = df.col(FLAGS_COL)?.getString(rowIdx) ?? '';
      if (text === '') return false;
      // Split on the per-message newline boundary (NOT '; ' — messages contain that).
      const parts = text.split('\n').filter((p) => p.length > 0);
      const tip = parts.length === 1 ?
        ui.divText(parts[0], {style: {maxWidth: '360px', whiteSpace: 'normal'}}) :
        ui.divV(parts.map((p) => ui.divText(`• ${p}`, {style: {maxWidth: '360px', whiteSpace: 'normal'}})));
      ui.tooltip.show(tip, x, y);
      return true;
    }
    if (!cell.isColHeader) return false;
    const name = cell.tableColumn?.name;
    if (name === undefined || name === null) return false;

    // (c) Leading (non-parameter) columns — study-grid-specific meaning.
    const leading = LEADING_TOOLTIPS[name];
    if (leading !== undefined) {
      ui.tooltip.show(richTipEl({title: leading.title, summary: leading.desc, note: leading.note}), x, y);
      return true;
    }

    // (a) Parameter columns — same structured `.nca-richtip` look as the inputs and
    // the aggregated grid: title → stat sentence (summary) → gloss (detail) →
    // fragility rule (note). The gloss DEFINITION is shared with the analysis view
    // (formatted to this grid's own labels — `Cmax`, not `C_max`).
    const col = df.col(name);
    const param = col?.getTag(PARAM_TAG) ?? null;
    const stat = col?.getTag(STAT_TAG) ?? null;
    if (param === null) return false;
    const hasStat = stat !== null && stat !== '';
    const glossText = studyParamGloss(param);
    const gloss = glossText !== null ? `${param} — ${glossText}.` : undefined;
    // Fragility rules this grid encodes as bold text (study-specific).
    const note = stat === 'n' ? 'Bold when n < 4 (fragile estimate).' :
      param === '%AUCextrap' ? 'Bold when > 20% (unreliable extrapolation).' : undefined;
    // With a stat, the stat sentence leads and the gloss is the muted detail; the
    // stat-less `%AUCextrap` has no sentence, so its gloss becomes the summary lede.
    const statSentence = hasStat ? studyStatSentence(stat as string, param) : null;
    ui.tooltip.show(statSentence !== null ?
      richTipEl({title: `${stat} ${param}`, summary: statSentence, detail: gloss, note}) :
      richTipEl({title: param, summary: gloss ?? '', note}), x, y);
    return true;
  });

  // v5 (chemist-sar-grid S4c) → context-panel-ct-profile: row selection / current-row → the
  // right context panel (master-detail). BUILD-ONCE, reactive (F4): the accordion + Forms viewer
  // + chart host are created once (stored on `df.temp`) and UPDATED per gesture, never rebuilt
  // per click (no viewer/subscription leak). The chart is CACHE-ONLY — it binds each row's
  // persisted `curve` (the SAME data the thumbnail draws), so it corresponds by construction (F1).
  // Gesture contract (BQ-4): a non-empty selection overlays every selected row's curve (F3); an
  // empty selection falls back to the single current row.
  df.onCurrentRowChanged.subscribe(() => {
    // A single click sets currentRow (selection stays empty) → single-curve path. When a
    // selection is active, the selection handler owns the panel — don't fight it.
    if (df.selection.trueCount === 0) updateCompoundPanel(study, grid, isStale);
  });
  wireSelectionOverlay(study, grid, isStale);
}

/** Temp-bag tag: the build-once compound-detail panel bundle (context-panel-ct-profile F4),
 *  stored on the grid frame's `df.temp` so it is created once per grid frame and reused across
 *  row clicks (no per-click rebuild → no Forms-viewer/subscription leak). */
const CTPROFILE_PANEL_TAG = 'ncaCtProfilePanel';
/** Temp-bag tag: the `df.selection.onChanged` overlay subscription (context-panel-ct-profile F3).
 *  Copies the {@link wireByGroupSelectionToMeanSdFilter} temp-bag teardown idiom
 *  (`BYGROUP_SELECTION_SUB_TAG`): clear the prior sub before subscribing so re-wiring is
 *  idempotent (AC-F3.4 — one live sub after open→close→reopen). */
const CTPROFILE_SELECTION_SUB_TAG = 'ncaCtProfileSelectionSub';

/** Drive the context panel off `df.selection.onChanged` (F3, copies
 *  {@link wireByGroupSelectionToMeanSdFilter}). Idempotent: unsubscribe the prior sub on this
 *  frame first, so open→close→reopen leaves exactly one live handler (AC-F3.4). Stored on the
 *  grid frame's `df.temp`; guarded with `!= null` (the `view.temp`/`df.temp` dart-boundary
 *  returns null — NOT undefined — for a never-set key; code-quality-guardrails). */
function wireSelectionOverlay(
  study: NcaStudy, grid: DG.Grid, isStale: StalePredicate = () => false,
): void {
  const df = grid.dataFrame;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = df.temp as Record<string, any>;
  const prev = bag[CTPROFILE_SELECTION_SUB_TAG] as {unsubscribe: () => void} | null;
  if (prev != null) prev.unsubscribe();
  bag[CTPROFILE_SELECTION_SUB_TAG] = df.selection.onChanged.subscribe(() => {
    updateCompoundPanel(study, grid, isStale);
  });
}

// ── Row-click drill-in: compound detail context panel (chemist-sar-grid S4c) ──

/** Resolve the clicked grid row back to its cached {@link CompoundTierRow}. Matches
 *  on the hidden composite-identity back-pointer `compoundId` (the tier's frozen
 *  `compositeGroupKey` tuple) + `route` — display-form-independent, so it works for
 *  the single `Group` column AND the Variant-2 per-key columns (v8). Returns the
 *  analysis entry + the matched row (first match; the rare multi-dose-norm split
 *  shares the same compound-level curve/scalars).
 *
 *  Label-INDEPENDENT: matches the compute-time-frozen `compoundId`, NEVER a
 *  reconstructed `compoundId === entry.label` (which broke on rename — the same root
 *  cause that leaked the stale label into the Group column). A whole-analysis row's
 *  `compoundId` is the analysis label at COMPUTE time; the back-pointer carries that
 *  exact frozen value, so the equality holds regardless of later renames. */
export function resolveClickedRow(
  study: NcaStudy, analysisId: string, compoundId: string, route: string,
): {entry: AnalysisEntry; row: CompoundTierRow} | null {
  const entry = study.analyses.find((a) => a.id === analysisId);
  if (entry === undefined) return null;
  // v15: a null cache (never confirmed) has no tier row to resolve (AC-S1.7.2 narrowing).
  const row = entry.summaryCache?.compoundTier.find(
    (r) => r.compoundId === compoundId && r.route === route);
  return row !== undefined ? {entry, row} : null;
}

/** A resolved (analysis entry, cached tier row) pair for the context panel. */
type ResolvedRow = {entry: AnalysisEntry; row: CompoundTierRow};

/** F3 legibility threshold (context-panel-ct-profile): overlays past this many series get a
 *  NON-TRUNCATING advisory (colors may be hard to distinguish). Design-step value (Rule-3):
 *  12 — the top of the ~10–12-class ceiling of categorical (qualitative) palette convention
 *  (e.g. ColorBrewer qualitative schemes), NOT a hard cap. Every selected row is STILL overlaid;
 *  the caption is advice, never an excuse to drop the scientist's selection (LBC-3, thesis #3). */
const CTPROFILE_LEGIBILITY_THRESHOLD = 12;

/** The event-driven (NONMEM) decline message — kept this cycle (F7, BQ-5): the raw source is an
 *  event log; the cached-curve TAD-correctness fixture check (AC-F7.1) gates removing this. */
const NONMEM_DECLINE_MSG =
  'Open the analysis view for the event-driven (NONMEM) concentration–time profile — ' +
  'the study overview shows it there with time-after-dose normalization.';

/** Semi-log concentration-time profile chart options (context-panel-ct-profile F1). Binds the
 *  shared {@link buildCurveFrame}/{@link buildOverlayCurveFrame} frame (`Time`/`Conc` = geoMean /
 *  `geoLower`/`geoUpper` band / `Series`), so the panel draws the SAME cached curve the grid
 *  thumbnail draws (correspondence by construction). `Series` colours/keys the line(s): the
 *  single pooled curve for one row, or one line per selected row (row-identity labels) on overlay.
 *  The ×/÷ gSD band rides `geoLower`/`geoUpper` (present on a fresh curve; absent → no whisker +
 *  the F2 "recompute" disclosure). */
const PROFILE_GRAPH_OPTS: Record<string, unknown> = {
  x: 'Time',
  y: 'Conc',
  colorColumnName: 'Series',
  linesByColumnName: 'Series',
  linesOrderColumnName: 'Time',
  linesWidth: 2,
  yWhiskerMinColumnName: 'geoLower',
  yWhiskerMaxColumnName: 'geoUpper',
  yAxisType: 'logarithmic',
  showRegressionLine: false,
  showXHistogram: false, showYHistogram: false,
  showFilteredOutPoints: false, showDropLines: false,
  showSizeSelector: false, showTitle: false,
  // Legend at the top (Series = the overlaid row-identity labels), color-by
  // selector hidden — the panel drives colour, not the viewer chrome.
  legendVisibility: 'Always', legendPosition: 'Top',
  showColorSelector: false,
};

/** The row-identity label for a resolved row — `Analysis · group-values · route` — used both as
 *  the overlay `Series` category (F3) and the single-row panel title. A whole-analysis row omits
 *  the group segment (its `groupValues` is empty). */
export function rowIdentityLabel({entry, row}: ResolvedRow): string {
  // Dedup group values by the SAME typographic canonical id the study grid merges on
  // (`normalizeGroupKey`, first-seen wins — buildStudyGrid's `keyCols`/`keyCellValue`).
  // Two `groupValues` entries can share a canonical key when two result columns carry the
  // same `friendlyName` (nca-analysis-open groupValues keying) — the grid collapses them
  // into ONE key column, so joining ALL raw values here printed the value twice
  // (`… · 1 · F · 1 · …`). Matching the grid's merge keeps the label == its key columns.
  const grp = row.isWholeAnalysis ? '' : (row.groupValues.length > 0 ?
    dedupGroupValuesByCanonical(row.groupValues).map((g) => g.value).join(' · ') : row.compoundId);
  return [entry.label, grp, row.route].filter((s) => s !== '').join(' · ');
}

/** Number of subjects behind a row's cross-dose pool disclosure, or `null` when the row does not
 *  pool across dose levels (F6). Reads the already-computed `CROSS_DOSE_POOL` warning on the tier
 *  row — zero schema change; restores the *disclosure* half of the E3 pool-and-disclose contract. */
function crossDosePoolN(row: CompoundTierRow): number | null {
  return row.warnings.some((w) => w.code === CROSS_DOSE_POOL_CODE) ? row.n : null;
}

/** True when a row carries a curve whose band was NOT computed (a legacy pre-v11 cached curve —
 *  every series' `geoLower`/`geoUpper` is absent). This is the F2 honesty distinguisher: an
 *  absent band (→ "recompute to show") is NOT the same as a present band with a `FLOAT_NULL`
 *  element at a genuine n=1 point (honest no-whisker, no disclosure). */
function hasLegacyAbsentBand(row: CompoundTierRow): boolean {
  const curve = row.curve;
  if (curve === null || curve.series.length === 0) return false;
  return curve.series.some((s) => s.geoLower === undefined || s.geoUpper === undefined);
}

/**
 * A NON-MUTED, icon-led disclosure bar (context-panel-ct-profile F6/F2, LBC-1). Distinct from
 * {@link emptyCaption}'s muted grey "nothing to show" register — it is a caveat about what you
 * ARE seeing, so it must NOT read as an empty state (an unnoticed disclosure discloses nothing,
 * thesis #4). Follows the compute-notice-banner idiom: a NEUTRAL severity
 * glyph (no background hue) + non-muted text.
 * Design-step decision (Rule-3): `warning` glyph for the cross-dose caveat (a correctness caveat),
 * `info` glyph for the band-not-computed notice.
 */
function disclosureBar(text: string, severity: 'info' | 'warning'): HTMLElement {
  const icon = ui.iconFA(severity === 'warning' ? 'exclamation-triangle' : 'info-circle');
  icon.style.color = 'var(--grey-6)'; // neutral glyph, no hue
  icon.style.marginRight = '6px';
  icon.style.flexShrink = '0';
  const label = ui.divText(text, {style: {whiteSpace: 'normal'}});
  label.style.color = 'var(--grey-7)'; // NON-muted (NOT emptyCaption's --grey-4)
  const bar = ui.divH([icon, label], 'nca-ct-disclosure');
  bar.style.alignItems = 'baseline';
  bar.style.margin = '4px 0';
  bar.style.fontSize = '12px';
  return bar;
}

/** The disclosures for the currently-shown row set (F6 cross-dose + F2 band-not-computed +
 *  F3 legibility). Aggregated across an overlay: one bar each, never repeated per row. Exported so
 *  the disclosure-appears behaviour (AC-F6.1, AC-F2.2/F2.3) is unit-testable without a live panel. */
export function buildDisclosures(
  resolved: ReadonlyArray<ResolvedRow>,
  isStale: (entry: AnalysisEntry) => boolean = () => false,
): HTMLElement[] {
  const bars: HTMLElement[] = [];
  // F-N.6 (GROK-20960): the curve drawn here comes from `summaryCache`, the same
  // cache the study grid marks `stale (refreshes on open)` when a DIFFERENT (or
  // unknown) sci-comp version computed it. The grid row carries that cue, but the
  // panel is a second reader of the same cached numbers and must not present them
  // as current on its own — a curve is read as data, not as a quotation of the
  // row. Same coordinator predicate, so the two cues cannot disagree.
  const staleN = resolved.filter((x) => isStale(x.entry)).length;
  if (staleN > 0) {
    bars.push(disclosureBar(
      resolved.length === 1 ?
        'Computed by a different NCA engine version than the one installed — the cached ' +
          'curve and parameters may predate a numbers-changing fix; reopen the analysis to refresh.' :
        `${staleN} of these were computed by a different NCA engine version than the one ` +
          'installed — those cached curves may predate a numbers-changing fix; reopen each ' +
          'analysis to refresh.', 'warning'));
  }
  // F3: overlay legibility advisory (non-truncating — every series is still drawn).
  const plottable = resolved.filter(
    (x) => !x.entry.modes.nonmem && x.row.curve !== null && x.row.curve.series.length > 0);
  if (plottable.length > CTPROFILE_LEGIBILITY_THRESHOLD) {
    bars.push(disclosureBar(
      `${plottable.length} compounds overlaid — colours may be hard to distinguish; ` +
      'narrow the selection or open the analysis view for facets.', 'info'));
  }
  // F6: a pooled-across-dose-levels curve is a correctness caveat.
  const crossN = resolved.map((x) => crossDosePoolN(x.row)).filter((n): n is number => n !== null);
  if (crossN.length > 0) {
    bars.push(disclosureBar(
      resolved.length === 1 ?
        `Curve pools across ≥2 dose levels (n = ${crossN[0]}); concentration scales with dose — ` +
          'group by Dose to stratify.' :
        `${crossN.length} of the overlaid curves pool across ≥2 dose levels; concentration scales ` +
          'with dose — group by Dose to stratify.', 'warning'));
  }
  // F2: a legacy cached curve has no computed variability band → disclose "recompute to show".
  if (resolved.some((x) => hasLegacyAbsentBand(x.row))) {
    bars.push(disclosureBar(
      'Variability band not computed for a cached curve here — recompute the analysis to show ' +
      'the ×/÷ gSD band.', 'info'));
  }
  return bars;
}

/**
 * Render the cache-only concentration-time profile for the active row set into `host`
 * (context-panel-ct-profile F1/F3). Cache-only + synchronous — binds each row's persisted
 * `curve` (the SAME data the thumbnail draws) via the shared {@link buildCurveFrame} /
 * {@link buildOverlayCurveFrame}. NO source re-resolve, NO re-derivation (the divergence this
 * cycle deletes). NONMEM rows are declined this cycle (F7, BQ-5 — kept until the AC-F7.1 fixture
 * check confirms cached-curve TAD correctness); a row with no cached curve renders an honest
 * caption (the params pane still stands). Every selected row is overlaid — never truncated (LBC-3).
 */
export function renderProfileChart(host: HTMLElement, resolved: ReadonlyArray<ResolvedRow>): void {
  ui.empty(host);
  const plottable = resolved.filter(
    (x) => !x.entry.modes.nonmem && x.row.curve !== null && x.row.curve.series.length > 0);
  if (plottable.length === 0) {
    const anyNonmem = resolved.some((x) => x.entry.modes.nonmem);
    host.appendChild(emptyCaption(anyNonmem ? NONMEM_DECLINE_MSG :
      'No measurable concentration–time points to plot.'));
    return;
  }
  const frame = plottable.length === 1 ?
    buildCurveFrame(plottable[0].row.curve!, rowIdentityLabel(plottable[0])) :
    buildOverlayCurveFrame(plottable.map((x) => ({curve: x.row.curve!, label: rowIdentityLabel(x)})));
  const chart = DG.Viewer.scatterPlot(frame, PROFILE_GRAPH_OPTS);
  chart.root.style.width = '100%';
  chart.root.style.height = '100%';
  host.appendChild(chart.root);
}

/** The native multi-record "Forms" viewer type (context-panel-ct-profile F5, BQ-3b). Plural
 *  "Forms" (NOT the singular single-record "Form"): it renders each row as its OWN form card, so an
 *  N-row multi-select shows N parameter forms — matching the chart's N overlaid curves. Not in the
 *  typed `VIEWER` enum (a core-but-untyped viewer); `ViewerType` accepts any string, so it is
 *  created by name via `fromType`. */
const FORMS_VIEWER_TYPE = 'Forms';

/**
 * Create the native "Forms" viewer for the "NCA parameters" pane (F5, BQ-3b; RG-CTP-1 mechanism
 * M1). Bound to the STUDY GRID's OWN frame (`grid.dataFrame`) so it is ASSOCIATED with the grid and
 * follows its selection/current row NATIVELY. The "Forms" viewer chooses which rows to render via
 * boolean row-source flags (NOT a `rowSource` enum). They are set DYNAMICALLY from the live
 * selection so a row is never rendered twice:
 *  - a non-empty selection → `showSelectedRows` only (one card per selected row — the multi-select
 *    case, matching the chart's overlaid curves; `showCurrentRow` OFF so the current row, if it is
 *    also selected, is not drawn a SECOND time — the duplicate-card bug);
 *  - an empty selection → `showCurrentRow` only (a plain row-click shows that one row — not blank).
 * `showMouseOverRow` off (hover must not swap the pane). `colorCode` OFF (user-directed 2026-07-30,
 * superseding the 2026-07-16 "on"): the platform's value colouring carries no legend and no hover
 * explanation anywhere in this pane, so a coloured field could not tell the scientist what the
 * colour meant — the same defect that made MRT / Vss unreadable in the results grid. This pane also
 * has no colour vocabulary of its own (design B4 keeps the study grid on format-not-colour), so
 * turning it off removes an encoding rather than hiding information. Fields ARE
 * restricted via the platform's `fieldsColumnNames` to the grid's currently-visible columns
 * ({@link studyParamsFormFields}) so the pane MIRRORS the grid's "Full" toggle: the compact
 * Mean-only SAR set when off, the full Mean|%CV|n + secondary set when on. This SUPERSEDES the
 * earlier "keep all columns" choice (2026-07-14) — the form showing %CV/n while the grid hid them
 * was an inconsistency (user-directed 2026-07-16). Returns `null` if the platform cannot construct
 * the viewer (→ caption).
 */
function createStudyParamsViewer(df: DG.DataFrame): DG.Viewer | null {
  const hasSelection = df.selection.trueCount > 0;
  try {
    return DG.Viewer.fromType(FORMS_VIEWER_TYPE, df, {
      showSelectedRows: hasSelection,
      showCurrentRow: !hasSelection,
      showMouseOverRow: false,
      colorCode: false,
      fieldsColumnNames: studyParamsFormFields(df),
    });
  } catch {
    return null;
  }
}

/** The build-once compound-detail panel bundle (F4), stored on the grid frame's `df.temp`. */
interface CompoundPanel {
  readonly root: HTMLElement;
  readonly titleHost: HTMLElement;
  readonly discHost: HTMLElement;
  readonly chartHost: HTMLElement;
  /** Host the "Forms" viewer mounts into. */
  readonly formHost: HTMLElement;
  /** The current "Forms" viewer bound to the study-grid frame (M1). REBUILT per gesture (a
   *  build-once viewer goes blank after the context panel is closed & reopened — its detached DOM
   *  stops repainting), so this holds the latest instance to `.close()` before the next rebuild.
   *  Mutable; `null` until the first render / if the platform could not construct it. */
  form: DG.Viewer | null;
}

/**
 * Build the two-accordion SHELL once (F4): "NCA parameters" (host for the native "Forms" viewer,
 * F5/M1) + "Concentration–time profile" (disclosure bars + chart host), BOTH expanded. Cached on
 * `df.temp` (guarded `!= null` — the dart-boundary returns null, not undefined, for a never-set
 * key) so the accordion container + its expand state survive across gestures. The Forms viewer and
 * the chart are NOT built here — they are (re)built per gesture (`renderParamsForm` /
 * `renderProfileChart`), because a build-once viewer's detached DOM stops repainting after the
 * context panel is closed & reopened (documented platform behaviour; F4/AC-F4.2 reconciled). All
 * hosts are full-width so the chart + form fill the context panel.
 */
function ensureCompoundPanel(grid: DG.Grid): CompoundPanel {
  const df = grid.dataFrame;
  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  const bag = df.temp as Record<string, any>;
  const existing = bag[CTPROFILE_PANEL_TAG] as CompoundPanel | null;
  if (existing != null) return existing;

  // Match Datagrok's native accordion-title style (`.d4-accordion-title label`, user-directed) —
  // set INLINE because that platform rule targets a `<label>` child, so a class on this plain div
  // does not inherit it. 18px / --grey-4 / single-line ellipsis, exactly the native header look.
  const titleHost = ui.div([], 'nca-ct-title');
  Object.assign(titleHost.style, {
    fontSize: '18px', height: '18px', color: 'var(--grey-4, #9497A0)', display: 'block',
    whiteSpace: 'nowrap', overflow: 'hidden', textOverflow: 'ellipsis',
    cursor: 'default', userSelect: 'all',
  });
  const discHost = ui.div([], 'nca-ct-disclosures');
  const chartHost = ui.div([], 'nca-compound-profile');
  chartHost.style.height = '300px';
  chartHost.style.width = '100%';
  const formHost = ui.div([], 'nca-ct-params');
  formHost.style.width = '100%';
  formHost.style.minHeight = '120px'; // the viewer needs a definite box to lay out into
  // The "Forms" viewer is (re)built per gesture in updateCompoundPanel (renderParamsForm) — NOT
  // here — because a build-once viewer goes blank after the panel is closed & reopened.

  const profilePane = ui.divV([discHost, chartHost]);
  profilePane.style.width = '100%';
  const acc = ui.accordion('nca-ct-profile');
  // The row title rides the accordion's NATIVE title area (`addTitle`) so it inherits Datagrok's
  // own panel-header styling — consistent with the "NCA parameters" / "Concentration–time profile"
  // pane headers below it (no bespoke h2 typography). Updated per gesture via `titleHost`.
  acc.addTitle(titleHost);
  acc.addPane('NCA parameters', () => formHost, true);
  acc.addPane('Concentration–time profile', () => profilePane, true);
  acc.root.style.width = '100%';

  const root = ui.divV([acc.root], 'nca-compound-detail');
  root.style.width = '100%';
  const panel: CompoundPanel = {root, titleHost, discHost, chartHost, formHost, form: null};
  bag[CTPROFILE_PANEL_TAG] = panel;
  return panel;
}

/**
 * (Re)build the grid-bound "Forms" viewer into the panel's form host (F5/M1). Disposes the prior
 * instance first (`.close()`) so viewers do not accumulate on the grid frame, then creates a fresh
 * one bound to `grid.dataFrame` — this is why it survives close→reopen (a build-once viewer's
 * detached DOM stops repainting; the chart is rebuilt the same way). The fresh viewer, with its
 * selection-driven `showSelectedRows`/`showCurrentRow` flags (see `createStudyParamsViewer`),
 * immediately reflects the live selection/current row.
 */
function renderParamsForm(panel: CompoundPanel, df: DG.DataFrame): void {
  if (panel.form !== null) {
    try {panel.form.close();} catch {/* disposing a stale viewer must not break the panel */}
    panel.form = null;
  }
  ui.empty(panel.formHost);
  const form = createStudyParamsViewer(df);
  panel.formHost.appendChild(form !== null ? form.root : emptyCaption('NCA parameters unavailable.'));
  if (form !== null) {
    form.root.style.width = '100%';
    panel.form = form;
  }
}

/**
 * Update the build-once context panel for the current gesture (F3/F4/F5/F6). Resolves the active
 * rows (non-empty selection → overlay of selected rows; else the single current row — BQ-4), then
 * updates the panel: the grid-bound Forms viewer (params, via {@link renderParamsForm}), the chart,
 * and the disclosures. The accordion shell is built once (`ensureCompoundPanel`). Sets
 * `grok.shell.o` so the context panel shows. Exported as the drill-in entry point — the live grid's
 * row-click / selection subscriptions (`configureStudyGrid`) call it, and it is the unit under test
 * for the drill-in (no bespoke `showCompoundDetail` wrapper — that path is exactly this one).
 */
export function updateCompoundPanel(
  study: NcaStudy, grid: DG.Grid, isStale: StalePredicate = () => false,
): void {
  const df = grid.dataFrame;
  const resolved = activeRowRefs(df)
    .map((r) => resolveClickedRow(study, r.analysisId, r.compoundId, r.route))
    .filter((x): x is ResolvedRow => x !== null);
  if (resolved.length === 0) return;

  const panel = ensureCompoundPanel(grid);

  // Attach the (build-once) panel to the context panel FIRST, so its hosts have a real width when
  // the viewers measure themselves — a scatter/form created in a detached, 0-width host renders at
  // a default small size and does not reflow (the "chart not full width" symptom).
  grok.shell.windows.showContextPanel = true;
  grok.shell.o = panel.root;

  // Title (native accordion-title styling via addTitle): the single row's identity, or an
  // N-selected header. Plain text so the accordion's own header style applies (no bespoke h2).
  panel.titleHost.textContent = resolved.length === 1 ? rowIdentityLabel(resolved[0]) :
    `${resolved.length} compounds selected`;

  // Params pane (M1): (re)build the grid-bound "Forms" viewer — it survives close→reopen (a
  // build-once viewer goes blank once its DOM is detached). Bound to grid.dataFrame; its
  // selection-driven show-flags (createStudyParamsViewer) reflect this exact selection/current row.
  renderParamsForm(panel, df);

  // Profile pane: disclosures + chart (into the now-attached, full-width host).
  ui.empty(panel.discHost);
  for (const bar of buildDisclosures(resolved, isStale)) panel.discHost.appendChild(bar);
  renderProfileChart(panel.chartHost, resolved);
}

/** The active grid rows driving the panel (F3, BQ-4): the SELECTED rows when the selection is
 *  non-empty (overlay), else the single current row. Each resolved to its back-pointer triple
 *  (`analysisId`/`compoundId`/`route` — the hidden composite-identity columns), display-form
 *  independent. */
function activeRowRefs(df: DG.DataFrame): Array<{analysisId: string; compoundId: string; route: string}> {
  const read = (i: number): {analysisId: string; compoundId: string; route: string} => ({
    analysisId: df.col(ANALYSIS_ID_COL)?.getString(i) ?? '',
    compoundId: df.col(COMPOUND_ID_COL)?.getString(i) ?? '',
    route: df.col(ROUTE_COL)?.getString(i) ?? '',
  });
  const sel = df.selection;
  if (sel.trueCount > 0) return Array.from(sel.getSelectedIndexes()).map(read);
  const idx = df.currentRowIdx;
  return idx >= 0 ? [read(idx)] : [];
}

/** Muted caption for an unavailable/empty profile (honest — the params pane still stands; the
 *  graph is best-effort). MUTED grey by design — reserved for "nothing to show" empty states,
 *  distinct from {@link disclosureBar}'s non-muted caveat register (F6/LBC-1). */
function emptyCaption(text: string): HTMLElement {
  const el = ui.divText(text);
  el.style.color = 'var(--grey-4, #9497A0)';
  return el;
}

/** Tier-1: analyses as a minimal table — Name (clickable → open & compute) and
 *  Last computed. Other summary columns + the per-compound drill-in were removed
 *  (user-directed 2026-06-29) because their single-valued projection misrepresents
 *  multi-compound / multi-route / multi-dose analyses. */
function buildTierOne(study: NcaStudy, refresh: () => void): HTMLElement {
  if (study.analyses.length === 0) {
    // Same empty-state style as the Hub's "No studies yet." caption.
    const empty = ui.divText('No analyses in this study yet.');
    empty.style.marginLeft = '2px';
    empty.style.marginTop = '12px';
    empty.style.color = 'var(--grey-4, #9497A0)';
    return empty;
  }

  const table = document.createElement('table');
  table.className = 'nca-study-overview-tier1';
  // Pack columns to content — the table sits in a flex column that would otherwise
  // stretch it full-width; this matches the (block-hosted) Hub studies table.
  table.style.width = 'max-content';
  const thead = document.createElement('thead');
  const htr = document.createElement('tr');
  for (const h of ['Name', 'Last computed', '']) {
    const th = document.createElement('th');
    th.style.textAlign = 'left';
    th.style.padding = '0 28px 2px 0';
    // Sticky header so column labels stay visible while a long analysis list
    // scrolls under them (chemist-sar-grid, user-directed 2026-07-07).
    th.style.position = 'sticky';
    th.style.top = '0';
    // No visible fill — match the per-compound grid's own column header, which reads
    // as fill-less (user-directed 2026-07-08). Kept opaque at the PAGE surface colour
    // (`--grey-0`, not a distinct `--grey-1` tint, and not `transparent`) so the
    // sticky header still hides rows scrolling underneath while showing no box.
    th.style.background = 'var(--grey-0, #fff)';
    const h3 = ui.h3(h);
    h3.style.margin = '0'; // h3 default block margins otherwise inflate the header→rows gap
    // Grid-header grey — match the per-compound grid's column-header text colour
    // (user-directed 2026-07-08: --grey-6).
    h3.style.color = 'var(--grey-6)';
    th.appendChild(h3);
    htr.appendChild(th);
  }
  thead.appendChild(htr);
  table.appendChild(thead);

  const tbody = document.createElement('tbody');
  for (const a of study.analyses)
    tbody.appendChild(analysisRow(a, study, refresh));
  table.appendChild(tbody);

  // Cap the list height with an internal scroll (Option A) so a large analysis list
  // never pushes the per-compound grid below the fold — the grid stays visible; the
  // list scrolls within itself. ~7 rows before scrolling (user-directed 2026-07-07).
  const scroll = ui.div([table], 'nca-study-overview-tier1-scroll');
  scroll.style.maxHeight = '220px';
  scroll.style.overflowY = 'auto';
  scroll.style.marginTop = '12px'; // gap below the "Analyses" header (was on the table)
  return scroll;
}

function analysisRow(a: AnalysisEntry, study: NcaStudy, refresh: () => void): HTMLTableRowElement {
  const tr = document.createElement('tr');
  // Name: 20-char-capped link (full name in the tooltip) → open & compute.
  const MAX = 20;
  const display = a.label.length > MAX ? `${a.label.slice(0, MAX)}…` : a.label;
  const name = ui.divText(display, 'nca-analysis-label');
  name.dataset.analysisId = a.id;
  name.style.color = 'var(--blue-1, #2083D5)';
  name.style.cursor = 'pointer';
  name.addEventListener('click', () => {void openAnalysis(a.id);});
  ui.tooltip.bind(name, richTip({title: a.label, summary: 'Open & compute this analysis'}));

  // Right-side actions (analogous to the Hub study row): copy-link / rename / delete — the write
  // actions per the caller's access (SD-13).
  const linkBtn = blueIcon('share', () => copyAnalysisLink(study.label, a.label),
    'Copy a shareable link to this analysis');
  const editBtn = study.access.canEdit ?
    blueIcon('pen', () => renameAnalysisFlow(study.id, a, refresh), 'Rename analysis') :
    disabledIcon('pen', NO_EDIT_TIP);
  const deleteBtn = study.access.canEdit ?
    blueIcon('trash-alt', () => confirmDeleteAnalysis(a, refresh), 'Delete analysis') :
    disabledIcon('trash-alt', NO_EDIT_TIP);
  const actions = ui.divH([linkBtn, editBtn, deleteBtn]);
  actions.style.gap = '12px';

  tr.appendChild(td(name));
  // When the stored results last changed (`results_computed_at`) — not the last reopen.
  tr.appendChild(td(a.summaryCache === null || a.resultsComputedAt == null ? '—' : shortDate(a.resultsComputedAt)));
  tr.appendChild(td(actions));
  return tr;
}

/** Copy a deep-link to a specific analysis (`/apps/Nca/<study label>/<analysis label>`).
 *  Keyed on labels via the shared `ncaDeepLinkPath` so the link matches what
 *  `routeNcaApp` resolves (was uuid-keyed — the router resolves by label). */
export function copyAnalysisLink(studyLabel: string, analysisLabel: string): void {
  const url = `${window.location.origin}${ncaDeepLinkPath(studyLabel, analysisLabel)}`;
  const clip = window.navigator?.clipboard;
  if (clip?.writeText) {
    clip.writeText(url).then(
      () => grok.shell.info('Analysis link copied to clipboard.'),
      () => grok.shell.info(url));
  } else
    grok.shell.info(url);
}

/** Rename an analysis (label only) within its study, then refresh. The label stays unique within
 *  the study (label-keyed routing; excludes the analysis's own label). */
export function renameAnalysisFlow(studyId: string, a: {id: string; label: string}, refresh: () => void): void {
  const input = ui.input.string('Analysis name', {value: a.label});
  ui.dialog('Rename analysis')
    .add(input.root)
    .onOK(() => {
      void (async () => {
        try {
          const repo = getNcaRepository();
          const study = (await repo.listStudies()).find((s) => s.id === studyId);
          const taken = new Set((study?.analyses ?? []).filter((x) => x.id !== a.id).map((x) => x.label));
          await repo.renameAnalysis(a.id, uniqueLabel((input.value ?? '').trim() || a.label, taken));
          refresh();
        } catch (e) {
          grok.shell.error(`Could not rename: ${e instanceof Error ? e.message : String(e)}`);
        }
      })();
    })
    .show();
}

/** Final delete of an analysis of its study behind a blocking confirm (user decision 2026-10-21:
 *  no restore). */
export function confirmDeleteAnalysis(a: {id: string; label: string}, refresh: () => void): void {
  const name = document.createElement('b');
  name.textContent = a.label;
  const body = document.createElement('div');
  body.style.maxWidth = '420px';
  body.style.whiteSpace = 'normal';
  body.append('Delete analysis "', name, `"? ${CANNOT_UNDO}`);
  ui.dialog('Delete analysis')
    .add(body)
    .onOK(() => {
      void (async () => {
        try {
          await getNcaRepository().deleteAnalysis(a.id);
          grok.shell.info(`Deleted analysis '${a.label}'.`);
          refresh();
        } catch (e) {
          grok.shell.error(`Could not delete: ${e instanceof Error ? e.message : String(e)}`);
        }
      })();
    })
    .show();
}

/** Mark the rows whose cached numbers predate the installed engine (F-N.6, the ONE wording). */
function markStale(project: NcaStudy, tierOne: HTMLElement): void {
  for (const a of project.analyses) {
    if (!isEntryStale(a)) continue;
    const el = tierOne.querySelector(`[data-analysis-id="${a.id}"]`);
    if (el !== null) el.appendChild(ui.divText(` · ${STALE_ENGINE_FLAG}`, 'nca-analysis-stale'));
  }
}
