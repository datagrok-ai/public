---
feature: charts
target_layer: playwright
coverage_type: regression
priority: p0
realizes_atlas: [charts.cp.configure-via-property-panel, charts.cp.persist-via-project-save-reopen]
realizes: [charts.sankey, charts.chord, charts.surface-plot, charts.radar, charts.tree, charts.globe, charts.word-cloud, charts.sunburst]
realized_as:
  - sunburst-spec.ts
pyramid_layer: integration
ui_coverage_responsibility:
  - add-viewer-sunburst
  - viewer-property-panel-gear
  - select-columns-dialog
  - viewer-context-menu-reset-view
  - viewer-save-layout
  - viewer-apply-layout
ui_coverage_delegated_to:
  sunburst-multi-selection: charts-ui.md
  sunburst-empty-category-click: charts-ui.md
produced_from: migrated
original_path: public/packages/UsageAnalysis/files/TestTrack/Charts/sunburst.md
migration_date: 2026-05-07
source_text_fixes: []
candidate_helpers: []
unresolved_ambiguities:
  - step-8-old-layout-fixture-sourcing-issue-2979
  - old-layout-file-location-ticket-commentary
  - inherit-from-grid-toggle-visibility-on-spgi-v2-vs-demog
  - helper-candidate-selecthierarchycolumns-viewer-columns
scope_reductions: []
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

This is an integration-level scenario: Sunburst is exercised together
with the property panel, the Select Columns dialog, Layouts, Projects,
and Filters — not as an isolated single-viewer smoke test.

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

> **Moved to** `charts-ui.md` (ui-only, canvas-gesture). The
> Click / Ctrl+Click / Ctrl+Shift+Click on canvas Sunburst segments
> has no automatable JS-API equivalent that exercises the actual
> hit-test + modifier-key dispatch. Manual scenario maintained
> separately; this section intentionally left blank in the
> playwright-layer scenario.

### Select / filter on empty category (SPGI_v2.csv)

> **Moved to** `charts-ui.md` (ui-only, canvas-gesture). Clicking
> the grey null-segment is canvas hit-testing on the null-category
> bucket — `df.selection.set` proxy does not verify the canvas
> hit-test routing. Manual scenario maintained separately.

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
> `System:DemoFiles`. Provisioning is unresolved. Downstream
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

- Sibling scenarios: `radar.md` owns the Add-Viewer + property-panel
  smoke coverage for the Charts section; this scenario is the
  integration-level counterpart and isn't meant to be merged into the
  smoke test.
- This runs at the playwright/UI-driving layer because the "Project
  save and reopen" and "Layout save and apply" steps need real UI
  navigation and menus, not just JS-API calls.
- `github-2954` (Sunburst date-column hierarchy handling) isn't
  specifically tested here — see `sunburst-date-column-bug.md`.
  `github-3412` (Sunburst × Scatterplot color-state pollution) isn't
  tested here either — see `sunburst-scatterplot-color-pollution-bug.md`.
- Setup Step 5 is self-cleaning: it deletes any project saved during
  the run.
- Project save/reopen is tested per-viewer rather than once centrally
  — the same save/reopen serialization bug class can manifest
  differently per viewer type. See also `radar-save-reopen-bug.md` for
  the Radar-specific GROK-18085 reproduction.
- The spec currently drives the DOM via `page.evaluate` +
  `dispatchEvent` rather than `page.locator(...).click()`. Refactoring
  to locator-based clicks is deferred to a follow-up pass (planned
  alongside the same refactor for `tree-spec.ts`) to avoid scope creep
  in this change.
