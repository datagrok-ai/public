---
feature: charts
target_layer: playwright
coverage_type: edge
priority: p1
realizes_atlas: [timelines-legend-click-to-filter-stability]
realizes: []
realized_as:
  - timelines-spec.ts
produced_from: atlas-driven
pyramid_layer: bug-focused
related_bugs:
  - GROK-19033
date_created: 2026-05-07
authored_by: orchestrator-test-designer-charts-migrate-2026-05-07
---

# Timelines viewer — legend filtering regression (GROK-19033)

Bug-focused regression scenario for the Charts package Timelines
viewer, covering GROK-19033: clicking a legend category to filter
Timelines could cause a white screen or console errors. This scenario
locks in the fix — legend click-to-filter, visual stability across
legend-visibility transitions (Auto / Always / Never), and a
mid-session `splitByColumnName` re-bind. This is not a general smoke
test for the Timelines viewer (there isn't one yet in this section) —
it specifically targets the legend-filtering regression.

## Setup

A clean Datagrok session is the only shared setup. Single dataset
across all scenario blocks: `System:AppData/ApiSamples/ae.csv`
(SDTM Adverse Events shape — USUBJID, AESTDY, AEENDY, AETERM,
AESEV, AESOC). This is the same dataset the Charts demo
`Visualization | General | Timelines` loads (per
`Charts/src/demos/demo.ts` `VIEWER_TABLES_PATH['Timelines']`) and
the same dataset the ApiSamples timelines example uses
(`scripts/ui/viewers/js viewers/timelines-viewer.js`). The file
`Charts/files/ae.csv` is the package-local copy; `ApiSamples`
publishes the same content under `System:AppData/ApiSamples/ae.csv`
which is reachable via `grok.data.files.openTable` and via
**Browse > Files > App Data > ApiSamples > ae.csv** in the UI.

## Scenarios

### Scenario 1: Legend click-to-filter on Timelines (GROK-19033 reproduction)

Steps:
1. Open `System:AppData/ApiSamples/ae.csv` (e.g. via
   **Browse > Files > App Data > ApiSamples**, or **File > Open**).
2. On the Menu Ribbon, click **Add viewer** and select **Timelines**.
   **Expected:** the Timelines viewer renders without console errors;
   default property values are set per atlas
   (splitByColumnName=`USUBJID`, startColumnName=`AESTDY`,
   endColumnName=`AEENDY`, colorColumnName=`AESOC` or first matching
   string column, marker=`circle`, lineWidth=`3`,
   showZoomSliders=`true`, legendVisibility=`Auto`).
3. Click the **Gear** icon on the Timelines viewer.
   **Expected:** the Context Panel opens with the Timelines
   property surface visible.
4. In the Context Panel, set **Color > Color Column Name** to
   **AESOC** (System Organ Class — string column with multiple
   distinct categories, drives the legend categories).
   **Expected:** the legend appears (visibility=Auto resolves to
   visible because there are categories) and lists the AESOC
   categories from `ae.csv` (e.g. `SKIN AND SUBCUTANEOUS TISSUE
   DISORDERS`, `GENERAL DISORDERS AND ADMINISTRATION SITE
   CONDITIONS`, `RENAL AND URINARY DISORDERS`,
   `METABOLISM AND NUTRITION DISORDERS`).
5. Click one of the AESOC legend categories (e.g.
   `SKIN AND SUBCUTANEOUS TISSUE DISORDERS`).
   **Expected:** the Timelines viewer filters cleanly to that
   category's intervals — only rows whose AESOC matches the
   clicked category remain rendered.
   **Expected (GROK-19033 invariant):** NO white screen, NO
   glitchy re-render, NO console error during the legend click.
   The viewer transitions smoothly from full to filtered state.
6. Click the same legend category a second time (toggle off).
   **Expected:** the filter clears; all AESOC categories' intervals
   are rendered again, with the same visual stability invariant
   (no white screens, no console errors).
7. Click two different legend categories in sequence (e.g.
   `RENAL AND URINARY DISORDERS`, then
   `METABOLISM AND NUTRITION DISORDERS`).
   **Expected:** each click toggles its category's visibility on
   the viewer; intermediate states are visually stable; final
   state shows only the categories whose legend entries are
   currently active.

### Scenario 2: Legend filter persists across legendVisibility transitions

Steps:
1. Continuing from Scenario 1's configured Timelines viewer
   (or repeat Setup + Scenario 1 Steps 1-4 to reach the same
   state). At least one AESOC legend category is currently
   filtered out (Scenario 1 Step 5 state).
2. In the Context Panel, set **legendVisibility** to **Always**.
   **Expected:** the legend remains visible regardless of
   category count; the previously toggled-off category remains
   filtered out (filter state survives the visibility-mode
   transition); no white screen / re-render glitch.
3. In the Context Panel, set **legendVisibility** to **Never**.
   **Expected:** the legend is hidden in the viewer; the
   underlying filter state is preserved (the viewer continues to
   exclude the previously toggled-off category's intervals).
4. In the Context Panel, set **legendVisibility** back to
   **Auto**.
   **Expected:** the legend reappears (Auto resolves to visible
   for >0 categories on a coloring column); the previously
   toggled-off category is still toggled off in the legend
   (visual state matches pre-transition); no white screen.
5. Click the previously toggled-off legend category to re-enable
   it.
   **Expected:** the category's intervals re-appear on the
   viewer, completing the round-trip; no console errors.

### Scenario 3: splitByColumnName re-bind mid-session preserves legend integrity

Steps:
1. Continuing from Scenario 2 (or repeat Setup + Scenario 1
   Steps 1-4 to reach a fresh AESOC-colored Timelines).
2. In the Context Panel, change **splitByColumnName** from
   `USUBJID` to a different categorical column on `ae.csv`
   (e.g. `AESEV` — Adverse Event Severity, values
   {MILD, MODERATE, SEVERE}). If `AESEV` is not selectable as
   a split-by column on the local atlas mapping, fall back to
   any other categorical string column listed in the
   property's selector dropdown.
   **Expected:** the viewer's lanes (rows) re-organize to one
   lane per distinct value of the new split-by column; the
   legend remains coloring-driven (still bound to AESOC);
   no white screen / no console error during the rebind.
3. Click an AESOC legend category in the new split-by
   configuration.
   **Expected:** the legend filter applies cleanly to the
   re-laned viewer (GROK-19033 invariant holds across
   split-by transitions); no white screen, no console error.
4. Revert **splitByColumnName** back to `USUBJID`.
   **Expected:** the original lane layout is restored; the
   legend filter from Step 3 either persists (if filter state
   is split-by-independent) or cleanly resets (acceptable
   either way) — the bug-class invariant is the absence of
   visual glitches and console errors, not a specific
   filter-persistence semantic.

## Notes

- This runs at the playwright/UI-driving layer because reproducing
  GROK-19033 requires actually clicking legend items in the DOM and
  checking for visual glitches across re-renders — that's not
  something a pure JS-API check can validate.
- There's no dedicated Timelines smoke scenario (Add Viewer →
  Timelines, basic property-panel sweep) in this section yet — this
  bug-focused scenario doesn't cover that ground. A general smoke
  scenario for Timelines is a planned follow-up, not covered here.
- The original GROK-19033 bug report was a brief 3-step repro (open
  data with Timelines, click a legend category, observe
  glitches/white screens). This scenario expands that into a
  deterministic sequence with an explicit dataset (`ae.csv`), an
  explicit legend-driving column (AESOC), and explicit category
  labels, so it's reliably reproducible.
- GROK-17222 (Line chart legend filtering) is the same bug class on a
  different viewer. It is not covered here — Line chart is a separate
  feature with its own test file.
