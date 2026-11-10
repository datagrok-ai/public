/**
 * Interactive guide for the NCA analysis view. Two mutually-exclusive
 * scenarios, determined by whether `resultsGrid` is non-null:
 *
 * **Pre-Compute (no results yet):** 6 pages — main source-data grid →
 * mapping pane → Subject / Time / Concentration inputs → Compute icon.
 * Each input page scrolls its anchor into view before the spotlight
 * renders, so a collapsed/scrolled-away input is still tour-reachable.
 *
 * **Post-Compute (results docked):** up to 6 pages — results table →
 * chart (facet preferred over Combined overlay) → Semi-log Y-axis
 * toggle → Group by (only when grouping candidates were produced) →
 * row-selection focus → Add to study. Every page is gated on its anchor being
 * actually on screen at click time (dock manager may have closed any
 * subset); missing anchors are silently skipped. If everything was
 * closed, `startNcaTour` surfaces a `grok.shell.warning` and returns
 * `null` rather than starting an empty tour.
 *
 * Triggered exclusively by the info icon next to the Compute icon (see
 * {@link mountMappingPane}). No autostart.
 *
 * The mapping pane resolves the dynamic refs at click time and hands
 * them in via {@link TourRefs} — this module never reads `view.temp`
 * directly, staying decoupled from the view-state keys.
 */

import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';

import {runDescriber} from './ui-describer';
import type {DescriptionPage} from './ui-describer';

/** References the tour needs from the mapping pane's mounting closure.
 *  The caller (`mountMappingPane`) resolves the dynamic refs
 *  (`resultsGrid`, `facetGrid`, `profileChart`) at click time, so the
 *  tour always sees the current post-Compute state. */
export interface TourRefs {
  /** The platform's built-in source-data grid (`view.grid.root`). */
  mainGrid: HTMLElement;
  /** The docked mapping pane root — the whole inputs panel. */
  mappingPane: HTMLElement;
  subjectInput: DG.InputBase;
  timeInput: DG.InputBase;
  concInput: DG.InputBase;
  /** The pane's Confirm button (D17 — it replaced the ribbon Compute icon). The field keeps its
   *  historical name; the tour step reads "Confirm the setup". */
  computeBtn: HTMLElement;
  /** The "Add to study" ribbon text button (Report stage only). Only ACTIONABLE after a
   *  successful Confirm and while the analysis is not in a study. The tour spotlights it only in
   *  the post-Confirm scenario. */
  addToStudyBtn: HTMLElement;
  /** Ribbon Settings gear (`cog`). Revealed on the ribbon after the first
   *  successful Compute (like Group by), not before. Opens the Settings dialog
   *  holding the Semi-log (chart scale) and Compact (results columns) toggles.
   *  The tour spotlights it only in the post-Compute scenario. */
  settingsBtn: HTMLElement;
  /** Ribbon "Group by" control root. `null` until a Compute produces
   *  grouping candidates (it is revealed on the ribbon only then). The tour
   *  shows its page only when non-null AND on screen. */
  groupByControl: HTMLElement | null;
  /** `null` before the first Compute; the docked results grid otherwise. */
  resultsGrid: DG.Grid | null;
  /** `null` pre-Compute, or when the subject set is empty; otherwise the
   *  wrapping div of the per-subject facet grid. (There is no upper subject
   *  bound since 2026-08-13 — a large study scrolls rather than suppressing the
   *  grid.) **Preferred** over `profileChart` as the tour's chart-page anchor
   *  when present. */
  facetGrid: HTMLDivElement | null;
  /** `null` before the first Compute or when the chart was closed; the
   *  Combined overlay chart viewer otherwise. Used as the chart-page
   *  anchor only when `facetGrid` is absent (large-N case). */
  profileChart: DG.Viewer | null;
}

/** Whether `el` is attached to the document AND occupies a non-empty
 *  rect that intersects the viewport. Used as the post-Compute
 *  "is this dockable still visible?" gate — dock manager can close any
 *  subset of {results, facet, combined chart} between Computes.
 *
 *  **`isConnected` here means VISIBLE, not ALIVE, and that is deliberate.** A
 *  background tab is detached from the document (measured — see
 *  `dock-layout-tests`), so it reads false while still being docked. For a
 *  spotlight tour that is the right test: an element with no rect cannot be
 *  highlighted. It does mean a page is dropped for a panel that merely sits
 *  behind another tab, so a CALLER that wants a specific page must front its tab
 *  first — `nca-mapping-pane`'s tour trigger calls `activateResultsTab` for
 *  exactly that reason. Do not "fix" this by switching to a docked-ness test:
 *  that would hand the tour an anchor it cannot draw on. */
function isOnScreen(el: HTMLElement): boolean {
  if (!el.isConnected) return false;
  const r = el.getBoundingClientRect();
  if (r.width === 0 || r.height === 0) return false;
  return r.bottom > 0 && r.top < window.innerHeight &&
    r.right > 0 && r.left < window.innerWidth;
}

/** Scroll `el` into view if any part of it sits outside the viewport.
 *  Used by the pre-Compute input pages so a collapsed parent section or
 *  a scrolled mapping pane still lets the tour spotlight land on a
 *  visible rect. `block: 'center'` keeps the popup well within the
 *  viewport rather than at an edge. */
function ensureVisible(el: HTMLElement): void {
  if (!el.isConnected) return;
  const r = el.getBoundingClientRect();
  const offscreen =
    r.bottom < 0 || r.top > window.innerHeight ||
    r.right < 0 || r.left > window.innerWidth;
  if (offscreen) el.scrollIntoView({block: 'center', inline: 'nearest'});
}

/** Build the page sequence and start the tour. Returns the Done button
 *  so callers can dismiss the tour programmatically (e.g. on view-close
 *  — `nca-tour-overlay` is attached to `document.body`, outside any
 *  view, so an active tour would otherwise persist). Returns `null` when
 *  the post-Compute scenario has zero on-screen anchors (every results
 *  dockable closed via dock manager) — the caller has nothing to clean
 *  up in that case. A warning is surfaced so the silent no-op is
 *  explained. */
export function startNcaTour(refs: TourRefs): HTMLButtonElement | null {
  const pages: DescriptionPage[] = refs.resultsGrid != null ?
    buildResultsPages(refs) :
    buildInputsPages(refs);
  if (pages.length === 0) {
    grok.shell.warning(
      'Interactive guide unavailable: no results or chart on screen. ' +
      'Re-open a closed pane via the dock manager, or click Confirm again.',
    );
    return null;
  }

  return runDescriber({
    pages,
    btnsText: {prev: 'Back', next: 'Next', done: 'Done'},
  });
}

/** Pre-Compute tour: source grid → mapping pane → Subject → Time →
 *  Concentration → Compute icon. Each input page has a `beforeShow` that
 *  scrolls it into view; the mapping pane and source grid are
 *  fixed-layout containers that are always visible. */
function buildInputsPages(refs: TourRefs): DescriptionPage[] {
  return [
    {
      root: refs.mainGrid,
      position: 'left',
      elements: {major: refs.mainGrid},
      description:
        '# Source data\n\n' +
        'The observations table you opened — one row per concentration ' +
        'measurement.',
    },
    {
      root: refs.mappingPane,
      position: 'right',
      elements: {major: refs.mappingPane},
      description:
        '# Mapping panel\n\n' +
        'Map columns, set units, tune NCA rules, and run the analysis.',
    },
    {
      root: refs.subjectInput.root,
      position: 'right',
      elements: {major: refs.subjectInput.root},
      beforeShow: () => ensureVisible(refs.subjectInput.root),
      description:
        '# Subject column\n\n' +
        'Identifies individual subjects. Each unique value defines one ' +
        'NCA profile.',
    },
    {
      root: refs.timeInput.root,
      position: 'right',
      elements: {major: refs.timeInput.root},
      beforeShow: () => ensureVisible(refs.timeInput.root),
      description:
        '# Time column\n\n' +
        'Sampling time points. Drives the x-axis and AUC integration windows.',
    },
    {
      root: refs.concInput.root,
      position: 'right',
      elements: {major: refs.concInput.root},
      beforeShow: () => ensureVisible(refs.concInput.root),
      description:
        '# Concentration column\n\n' +
        'Measured plasma concentrations. Feeds Cmax, AUC, λz, t½, CL, Vz.',
    },
    {
      // nca-workflow-stages F-S3.10: the ribbon Compute icon is gone (D17) — the step anchors to
      // `Confirm` at the bottom of the mapping pane, which may sit below the fold.
      root: refs.computeBtn,
      position: 'top',
      elements: {major: refs.computeBtn},
      beforeShow: () => ensureVisible(refs.computeBtn),
      description:
        '# Confirm the setup\n\n' +
        'Compute NCA parameters for every subject, then review the λz fits that need a decision.',
    },
  ];
}

/** Post-Compute tour: results table → chart (facet preferred) → Y-axis
 *  toggle → selection focus. Each page is gated on its anchor being
 *  actually on screen — dock manager close, drag-detach, or a
 *  collapsed ribbon group will silently drop the matching page. */
function buildResultsPages(refs: TourRefs): DescriptionPage[] {
  const grid = refs.resultsGrid!;
  const gridVisible = isOnScreen(grid.root);
  const facetVisible = refs.facetGrid != null && isOnScreen(refs.facetGrid);
  const overlayVisible = refs.profileChart != null && isOnScreen(refs.profileChart.root);
  const settingsVisible = isOnScreen(refs.settingsBtn);
  const groupByVisible = refs.groupByControl != null && isOnScreen(refs.groupByControl);
  const addToStudyVisible = isOnScreen(refs.addToStudyBtn);
  const chartVisible = facetVisible || overlayVisible;

  const pages: DescriptionPage[] = [];

  if (gridVisible) {
    pages.push({
      root: grid.root,
      position: 'top',
      elements: {major: grid.root},
      description:
        '# Results table\n\n' +
        'One row per subject with Cmax, t_max, AUC, half-life, CL, Vz, λz ' +
        'and diagnostics. Anomaly columns are colour-coded — hover any ' +
        'header for the rule and parameter description.',
    });
  }

  if (facetVisible) {
    pages.push({
      root: refs.facetGrid!,
      position: 'left',
      elements: {major: refs.facetGrid!},
      description:
        '# Per-subject profile grid\n\n' +
        'One small plot per subject — concentration vs time. Switch to the ' +
        'Combined tab to overlay every subject on shared axes.',
    });
  } else if (overlayVisible) {
    pages.push({
      root: refs.profileChart!.root,
      position: 'left',
      elements: {major: refs.profileChart!.root},
      description:
        '# Concentration–time profile\n\n' +
        'Every subject overlaid on a single chart. Shown when there are ' +
        'too many subjects (>20) to render a readable facet grid.',
    });
  }

  if (settingsVisible) {
    pages.push({
      root: refs.settingsBtn,
      position: 'bottom',
      elements: {major: refs.settingsBtn},
      description:
        '# Settings\n\n' +
        'Open the gear to change display settings: **Semi-log Y scale** switches the ' +
        'concentration axis between log10 (PK convention) and linear on the ' +
        'Combined and Facet charts, **Compact NCA results** switches the results ' +
        'table between the parameter-only and full column sets, and **Violated column** ' +
        'shows or hides the violated λz criteria as pills beside Acceptance.',
    });
  }

  if (groupByVisible) {
    pages.push({
      root: refs.groupByControl!,
      position: 'bottom',
      elements: {major: refs.groupByControl!},
      description:
        '# Group by\n\n' +
        'Summarise the per-subject results into a per-group table — mean, SD, ' +
        '%CV and N for each NCA parameter. Pick a metadata column carried from ' +
        'your input (e.g. compound, dose, sex) to dock a "By group" tab; leave ' +
        'blank for per-subject results only. Partial/missing values are excluded ' +
        'from each group mean (N reflects the actual contributors).',
    });
  }

  if (gridVisible && chartVisible) {
    const focusText = facetVisible ?
      'Click a row in the results table — that subject\'s tile stays ' +
      'bright while the others dim. Ctrl/Shift-click to keep several. ' +
      'Clear the selection to restore all tiles.' :
      'Click a row in the results table — only that subject\'s curve ' +
      'stays on the chart. Ctrl/Shift-click to compare. Clear the ' +
      'selection to bring all curves back.';
    pages.push({
      root: grid.root,
      position: 'top',
      elements: {major: grid.root},
      description: '# Focus on a subject\n\n' + focusText,
    });
  }

  // Final call-to-action: persist the computed analysis into a study. The icon
  // lives on the ribbon (always in the viewport), but is only meaningful after a
  // successful Compute — which is exactly the scenario this page belongs to. Gated
  // on-screen like every other post-Compute page (a collapsed ribbon would drop it).
  if (addToStudyVisible) {
    pages.push({
      root: refs.addToStudyBtn,
      position: 'bottom',
      elements: {major: refs.addToStudyBtn},
      description:
        '# Add to study\n\n' +
        'Save this analysis into a study — its parameter table and summary are ' +
        'written to the study without recomputing. Collect several analyses in one ' +
        'study to compare compounds or dose arms side by side. Open a study from the ' +
        'NCA browser tree in the left sidebar.',
    });
  }

  return pages;
}
