---
feature: charts
sub_features_covered:
  - charts.timelines
  - charts.timelines.legend-visibility
  - charts.timelines.split-by-column
  - charts.timelines.color-column
  - charts.timelines.start-column
  - charts.timelines.end-column
target_layer: playwright
coverage_type: edge
produced_from: atlas-driven
pyramid_layer: bug-focused
related_bugs:
  - GROK-19033
date_created: 2026-05-07
authored_by: orchestrator-test-designer-charts-migrate-2026-05-07
---

# Timelines viewer — legend filtering regression (GROK-19033)

Bug-focused regression scenario for the Charts package Timelines viewer.
Closes the coverage gap recorded in `scenario-chains/charts.yaml` rev 1
(`bug_match_attempts_skipped` for GROK-19033, skip_category
`no_affecting_scenarios`) — there is no other timelines.md in
TestTrack/Charts/, and the Charts atlas declares
`charts.timelines.*` sub_features that no manual scenario walks.

`pyramid_layer: bug_focused` per orchestrator Edit 4 SR routing
dispatch table — this scenario exists primarily to lock in the
GROK-19033 reproduction (legend click-to-filter; visual stability
across legend-visibility transitions Auto / Always / Never;
splitByColumnName re-bind mid-session). It is NOT the smoke slot for
the Timelines viewer; the smoke slot remains open (see Notes).

`coverage_type: edge` per the canonical enum
{smoke, regression, edge, perf}: this is a negative-path / bug-class
coverage of the legend filtering interaction. The companion
sister-bug GROK-17222 (Line chart legend filtering, Bug Library
Legend) and the GROK-19033 ticket's fix-iteration comment indicate
this class of bug regresses; an `edge` scenario is the dominant
shape. (See Notes for the `regression` alternative considered.)

`related_bugs: [GROK-19033]` — canonical bug-link. The atlas
`edge_cases` entry for GROK-19033 (charts.yaml lines ~1026-1038)
describes this exact reproduction class; this scenario implements
the test coverage marked `test_coverage.exists: false` in
`known_issues` for GROK-19033 (charts.yaml lines ~1116-1124).

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

- **target_layer: playwright rationale.** Sibling scenarios in
  `TestTrack/Charts/` (radar.md, sunburst.md, tree.md) all use
  `target_layer: playwright`, and the Charts package's existing
  `*-spec.ts` files (radar-spec.ts, sunburst-spec.ts, tree-spec.ts
  per existing-test-index.yaml) are at the playwright layer. The
  GROK-19033 reproduction requires real ECharts legend DOM
  interaction (legend item click, visual-stability assertion across
  re-renders) which the playwright layer handles natively;
  uitests-package would underspecify the visual-glitch invariant.

- **coverage_type: edge vs regression — judgment call.** Chose
  `edge` because the dominant shape is bug-class regression of a
  single negative-path interaction (legend click filtering
  glitching). The `regression` alternative was considered because
  the scenario walks 6 sub_features (legend-visibility,
  split-by-column, color-column, start-column, end-column, plus
  charts.timelines parent), but the scenario blocks all share the
  same bug-class invariant ("no white screen / no console error
  during legend filter transitions") rather than independent
  feature-shape regressions. The atlas `edge_cases` entry for
  GROK-19033 explicitly frames this as edge coverage class
  ("focus on visual stability"). `edge` also satisfies Critic
  Gate A's A-STRUCT-02 (≥1 edge/perf scenario in section is
  required when atlas has edge_cases — the Charts section
  currently has only `smoke` (radar.md) and `regression`
  (sunburst.md) plus tree.md; this scenario adds the missing
  `edge` slot).

- **sub_features_covered scope rationale.** Listed the three
  bug-affected sub_features per `bug-library/charts.yaml`
  GROK-19033 `affects` (`charts.timelines`,
  `charts.timelines.legend-visibility`,
  `charts.timelines.split-by-column`) plus three additional
  sub_features actually exercised by Steps 2 and 4 of Scenario 1
  (`charts.timelines.color-column` — required to populate the
  legend; `charts.timelines.start-column` and
  `charts.timelines.end-column` — required to render any
  intervals on the SDTM AESTDY/AEENDY shape). The remaining 6
  Timelines sub_features (show-open-intervals, event-tooltip,
  marker, line-width, date-format, axis-pointer,
  show-zoom-sliders) are NOT exercised here — bug-focused not
  exhaustive, per the SR routing intent.

- **Density (Critic A's A-STRUCT-01 mechanical check).** Each
  scenario block exercises ≥3 sub_features in interaction:
  Scenario 1 = {timelines, color-column, start-column,
  end-column, legend-visibility} → 5; Scenario 2 = {timelines,
  legend-visibility, color-column} → 3; Scenario 3 = {timelines,
  split-by-column, legend-visibility, color-column} → 4. Average
  density per scenario = 4.0, well above the ≥2 threshold.

- **Smoke slot deferred.** The Charts section currently has no
  `timelines.md` smoke scenario at any pyramid layer. This
  bug-focused scenario does NOT cover the canonical
  `Add Viewer → Timelines` smoke surface for the Timelines viewer
  (the legend-filter bug class is the focus). A separate smoke
  scenario is a follow-up — deferred per Lattice Rule 13 /
  A-MERIT-02 with cited dependency: "deferred — orchestrator SR
  dispatch routes this fallback to `bug_focused_candidates` slot
  only; smoke-slot scenario authoring is a separate Phase 1
  routing decision, not a sub-task of this dispatch."

- **Helpers (registered, available for downstream Automator).**
  From `helpers-registry.yaml :: grok_test_layer` (revision 1):
  - `findViewer(viewerName, view)` — locate the Timelines viewer
    after Add viewer (Steps Scenario 1.2, all property panel
    steps).
  - `isViewerPresent(viewers, viewerName)` — assertion in
    Scenario 1 Step 2.
  - `readDataframe(tableName)` — alternative to UI **File > Open**
    if the spec-style harness prefers programmatic loading.
  Sibling spec helpers (used by `radar-spec.ts`, `sunburst-spec.ts`)
  from `public/packages/UsageAnalysis/files/TestTrack/spec-login.ts`:
  - `softStep`, `loginToDatagrok`, `specTestOptions` — standard
    spec-layer wrappers; expected to be used by the downstream
    `timelines-spec.ts` Automator output.

- **Helper candidates surfaced (NOT invented here, registry write
  is downstream).** Two patterns recur in this scenario that have
  no direct registry helper:
  - `clickLegendCategory(viewer, categoryLabel)` — locate a legend
    item by text inside an EChart-rendered `.echarts-legend` (or
    Datagrok legend wrapper) and dispatch the click that toggles
    the category's filter. Used by Scenario 1 Steps 5/6/7,
    Scenario 2 Step 5, Scenario 3 Step 3.
  - `assertNoWhiteScreen(viewer, action)` — wrap a UI action and
    assert the viewer's root DOM has non-empty rendered children
    + no `console.error` matching `/Cannot read|undefined/`
    during/after. Used by every legend-click and visibility-mode
    transition in all three scenarios. Closely mirrors
    `isErrorBallon` semantics in the registry (line 100) but
    targets ECharts-rendered viewers rather than balloon
    notifications. Surface to registry maintainers — do not
    invent the helper inside the spec.

- **GROK-19033 specific reproduction note.** The bug ticket's
  `reproduction` is 3 steps (Open data with Timelines / Click
  legend category / Observed: glitches + white screens). This
  scenario expands those 3 steps into a deterministic UI-driven
  sequence with explicit dataset (ae.csv), explicit
  legend-driving column (AESOC for color), and explicit
  category labels (e.g. `SKIN AND SUBCUTANEOUS TISSUE DISORDERS`)
  so the spec is reproducible. The
  `expected: "Clicking a legend category filters Timelines to
  that category's intervals with no visual glitches or white
  screens"` from the bug ticket is the assertion class used in
  every "Expected (GROK-19033 invariant)" line.

- **Sister-bug awareness.** Bug-library `edge_case_for_atlas`
  for GROK-19033 cites GROK-17222 (Line chart legend filtering)
  as the same bug class manifesting on a different viewer. This
  scenario does NOT cover Line chart — that is a separate file
  under a different feature. Cited for awareness only;
  cross-viewer regression coverage is a separate Phase scope.

- **Dataset metadata (for parity with sibling scenarios).**
  Single dataset: `System:AppData/ApiSamples/ae.csv`.
