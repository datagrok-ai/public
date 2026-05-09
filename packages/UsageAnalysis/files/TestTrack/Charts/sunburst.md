---
feature: charts
sub_features_covered:
  - charts.sunburst
  - charts.sunburst.title
  - charts.sunburst.on-click
  - charts.sunburst.inherit-from-grid
  - charts.sunburst.include-nulls
  - charts.echart-base
target_layer: playwright
coverage_type: regression
pyramid_layer: integration
ui_coverage_responsibility:
  - add-viewer-sunburst
  - viewer-property-panel-gear
  - select-columns-dialog
  - viewer-context-menu-reset-view
  - sunburst-multi-selection
  - viewer-save-layout
  - viewer-apply-layout
ui_coverage_delegated_to: null
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/sunburst.md
migration_date: 2026-05-07
migration_report: sunburst-migration-report.md
related_bugs:
  - github-2954
  - github-3412
---

# Sunburst viewer

Multi-subsystem Sunburst integration scenario: viewer creation, property
panel + Select Columns dialog, hierarchy rendering, inherit-from-grid
coloring, include-nulls toggle, view reset (double-click + context
menu), multi-selection gestures (Click / Ctrl+Click / Ctrl+Shift+Click),
empty-category select/filter, project save → close-all → reopen
restoration, layout save + apply, old-layout compatibility (issue #2979),
and collaborative filtering (internal viewer filter ∧ panel filter).

Pyramid layer is `integration` per chain (charts.yaml rev 1) — Sunburst
co-exercises the viewer + Layouts + Projects + Filters subsystems.
This is NOT a single-viewer ui-smoke; do not consolidate.

## Setup

1. Authenticate as test user.
2. Open `System:DemoFiles/demog.csv` (target table for inherit-from-grid
   and collaborative-filtering sub-flows).
3. Open `System:DemoFiles/SPGI_v2.csv` (primary source driving
   hierarchy / nulls / project-save sub-flows).
4. For each opened table view, add a **Sunburst** viewer via
   **Add viewer** → **Sunburst**.
5. Cleanup: delete any project saved during the run (Step in
   `### Project save and reopen`); restore default layouts.

## Scenarios

### Viewer creation smoke

1. Confirm the Sunburst viewer is present on **SPGI_v2.csv** without
   console errors.
2. Confirm the Sunburst viewer is present on **demog.csv** without
   console errors.

### Property panel opens

1. Click the **Gear** icon in the Sunburst viewer (SPGI_v2 view).
2. Verify the context panel with viewer properties opens and shows
   the Sunburst property set.

### Table switching

1. With the property panel open, switch the bound table between
   **SPGI_v2.csv** and **demog.csv** via the property panel's table
   selector.
2. Verify the viewer re-renders against the new bound table without
   errors.

### Hierarchy configuration via Select Columns dialog

1. On the SPGI_v2 Sunburst, open the **Select Columns** dialog from
   the property panel.
2. Choose 2–4 columns and click **OK**.
3. Verify the hierarchy updates to reflect the selected columns
   (rings ordered top-to-bottom by selection order).
4. Reopen the **Select Columns** dialog and use the search box to
   locate a column, then click **Cancel**.
5. Verify no hierarchy changes are applied.
6. Verify segment label text is fully visible when there is enough
   space, otherwise hidden with a tooltip showing the value on hover.
7. Verify structural columns render at the correct ring size relative
   to their nesting level.

### Inherit from grid (demog.csv)

1. On the **demog.csv** Sunburst, open the **Select Columns** dialog
   and choose the **SEX** column; click **OK**.
2. Enable the **Inherit from grid** property in the Sunburst property
   panel.
3. Apply categorical coloring to the **SEX** column in the grid.
4. Verify the Sunburst viewer reflects the grid colors for the SEX
   ring.
5. Change the categorical coloring on the grid SEX column.
6. Verify the Sunburst viewer updates to reflect the new colors.

### Include nulls (SPGI_v2.csv)

1. On the SPGI_v2 Sunburst, open **Select Columns** and choose the
   **Core** and **R101** columns; click **OK**.
2. Enable the **Include nulls** property.
3. Verify grey segments appear for null values in either column.
4. Disable the **Include nulls** property.
5. Verify the grey null segments disappear.

### View reset

1. Drill into the Sunburst by clicking a ring segment, then
   double-click on empty space inside the viewer.
2. Verify the view resets to its initial (top-level) state.
3. Drill in again, then open the viewer context menu and click
   **Reset View**.
4. Verify the view resets to its initial state.

### Multi-selection behaviour

1. Click a single segment.
2. Verify exactly that segment is selected and the corresponding grid
   rows are selected.
3. Hold **Ctrl** and click an additional segment.
4. Verify both segments are selected and the grid row selection
   reflects the union.
5. Hold **Ctrl + Shift** and click one of the selected segments.
6. Verify that segment is deselected and the grid row selection
   reflects the new union.

### Select / filter on empty category (SPGI_v2.csv)

1. On the SPGI_v2 Sunburst, configure a column known to contain nulls
   (e.g. **Sampling Time**) via **Select Columns**.
2. Click on the null (grey) segment.
3. Verify the segment behaves like any other category — the
   corresponding rows (those with null in that column) are selected
   or filtered (per `onClick` action).

### Project save and reopen

1. On the SPGI_v2 Sunburst, configure 3–4 hierarchy columns via
   **Select Columns**.
2. Save the project (Save Project dialog → keep defaults → **OK**).
   Cancel any auto-share dialog that opens.
3. Close all views.
4. Reopen the saved project from **Browse > Dashboards** (or the
   recent-projects list).
5. Verify the Sunburst viewer is restored with the same hierarchy
   columns, bound table, and visible state.

### Layout save and apply

1. With the SPGI_v2 Sunburst still configured, save the current
   layout (via the **Layouts** menu / **Save** action — not Ctrl+S
   which saves the project).
2. Reset or change the Sunburst configuration (e.g. clear hierarchy
   columns) so the layout's effect is observable.
3. Apply the saved layout.
4. Verify the Sunburst viewer is restored to the configuration
   captured by the layout (hierarchy columns, properties).

### Old layout compatibility (issue #2979)

1. Apply / open the layout from issue
   [#2979](https://github.com/datagrok-ai/public/issues/2979) on
   the corresponding source table.
2. Verify the Sunburst viewer shows the columns the layout encodes,
   AND the columns presented in the **Select Columns** dialog are
   in sync with what the viewer is rendering.

> NOTE: The layout fixture for issue #2979 is referenced by URL in
> the original scenario but is NOT stored in the repo or
> `System:DemoFiles`. Provisioning is unresolved — see the
> Unresolved ambiguities section of the migration report. Downstream
> Automator must either commit a fixture layout, provision via
> `grok.dapi.layouts`, or `test.skip` this scenario citing the
> missing fixture.

### Collaborative filtering (demog.csv)

1. On the **demog.csv** Sunburst, configure 2–3 hierarchy columns
   via **Select Columns**.
2. Click a Sunburst segment to apply an internal viewer filter
   (`onClick` action set to **Filter**, or use the segment's
   filter context-menu entry).
3. Open the panel filter (filter widget) on the same view and
   apply a filter on a different column.
4. Verify the filter set in effect is the **intersection** of the
   internal Sunburst filter and the panel filter — the grid shows
   only rows matching both, and the Sunburst re-renders against
   the intersection.

## Notes

- **Origin:** migrated from `sunburst.md` original (9 numbered steps,
  several with sub-steps).
- **`pyramid_layer: integration`** per chain (`charts.yaml` rev 1):
  multi-subsystem coupling — Sunburst viewer + Property panel + Select
  Columns dialog + Layouts + Projects (save/close-all/reopen) + Filters
  (collaborative). Do not consolidate to ui-smoke; the smoke slot in
  the chain is owned by `radar.md`.
- **`target_layer: playwright`** because Step "Project save and reopen"
  requires save → close-all → reopen across a navigation boundary, and
  Step "Layout save and apply" requires UI menu interaction. Sibling
  spec `sunburst-spec.ts` already exists at the playwright layer per
  `existing-test-index.yaml` (currently `test.skip`'d).
- **Helpers (registry-only references; no helper invented here):**
  - `findViewer(viewerName, view)` — locate the Sunburst viewer.
  - `isViewerPresent(viewers, viewerName)` — assertion in
    `### Viewer creation smoke`.
  - `uploadProject(projectName, tableInfo, view, df)` — registry helper
    that may simplify `### Project save and reopen` Step 2; downstream
    Automator decides whether to use it or drive the Save dialog UI
    explicitly. (See `helpers-registry.yaml :: grok_test_layer`.)
- **Bug coverage cited in `related_bugs`:**
  - `github-2954` (Sunburst date-column selector validation) —
    `### Hierarchy configuration via Select Columns dialog` and
    atlas edge_cases reference column-type defense-in-depth. The
    explicit date-column repro is NOT exercised here (chain analyzer
    flagged this as a coverage gap; cite for awareness, do not invent
    the date-column step).
  - `github-3412` (cross-viewer Sunburst × Scatterplot color editor
    pollution) — `### Inherit from grid (demog.csv)` exercises the
    inherit-from-grid path that the bug pollutes; the cross-viewer
    Scatterplot pairing is NOT added here (single-viewer scenario per
    chain).
- **No fixture extraction** — chain `fixtures_extracted: []`.
  `sunburst.md` is independent (`depends_on: []`).
- **Self-cleaning** — Setup Step 5 deletes any project saved during
  the run.
- **Remediation cycle scope decision (charts-remediate-2026-05-09):**
  Critic E canonical subagent surfaced `### Project save and reopen`
  heading (L146-156) as uncovered in predecessor cycle
  charts-automator-only-2026-05-08 (FAIL with E-TRACE-02 +
  E-LAYER-COMPLIANCE-01). Migrator Decision: **KEEP "Project save
  and reopen"** as a softStep in sunburst-spec.ts. Rationale:
  Sunburst-specialty UI (Select-Columns dialog, Reset View, Save
  Layout) is exercised on the saved project in a way Radar-specific
  tests don't witness. The save/reopen serialization regression class
  (github-3412 for Sunburst×Scatterplot is a related precedent) can
  manifest differently per viewer. Distinct from
  radar-save-reopen-bug-spec.ts (Radar-specific GROK-18085
  reproduction). Automator implements between current Step 7
  (Layout) and Step 8 (#2979): save project, closeAll, find by
  name, open, assert Sunburst viewer present + hierarchy preserved.
  Cleanup deletes saved project in `finally`.
- **E-LAYER-COMPLIANCE-01 follow-up (deferred):** spec body
  currently uses `page.evaluate` + `dispatchEvent` exclusively for
  DOM driving; canonical critic flagged FAIL on strict regex.
  Per legend-* failed_attempts precedent, refactor to
  `page.locator(...).click()` is the canonical fix. **Deferred to
  next charts cycle** (charts-evaluate-extract-2026-05-09 in
  decision-log) — pairing the refactor with tree-spec.ts (which
  is being refactored in this remediation cycle) is the symmetric
  approach. Sunburst defer rationale: this cycle's primary scope
  is Project save+reopen softStep insertion; bundling a body-wide
  refactor risks scope creep and lengthens Validator B re-run
  surface. Documented in scenario Notes for traceability.
