---
feature: linechart
realizes_atlas:
  - linechart.cp.legend-color-and-persistence
realizes:
  - viewers.line-chart
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - legend-color-and-persistence-spec.ts
related_bugs:
  - id: github-1498
    status: fixed
  - id: GROK-17278
    status: fixed
  - id: GROK-19825
    status: fixed
expected_results:
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      df.filter.trueCount is strictly less than the full row count — filtering
      by a legend category reduces visible rows.
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      The color assigned to the remaining category's line is unchanged from its
      pre-filter value (does not turn blue or shift; canvas render is not
      directly assertable, so assert via the persisted color property value).
  - anchor: "Scenario 1 Step 8"
    expectation: >-
      After resetting the legend filter, df.filter.trueCount is restored to the
      full row count (all categories visible again).
  - anchor: "Scenario 2 Step 8"
    expectation: >-
      After saving the layout and reopening from it, the per-category color
      property reads back the same value that was set before saving —
      persistence round-trip for layout (GROK-17278 regression guard).
  - anchor: "Scenario 2 Step 12"
    expectation: >-
      After saving the project via the SAVE button, closing all, and reopening
      the project, the Line chart viewer is restored and the per-category color
      reads back the value set before saving (GROK-17278 — project
      serialization preserves per-category color).
  - anchor: "Scenario 2 Step 13"
    expectation: >-
      The markers legend is present in the reopened project's Line chart viewer
      (GROK-19825 — the legend must not disappear on project reopen).
---

# Line Chart — Legend filter-color and layout/project persistence

## Setup

1. Open the SPGI dataset (use `readDataframe` helper for `SPGI.csv`).
2. Add a Line Chart viewer to the table view (use `findViewer` helper or add via the Toolbox).
3. In the viewer property panel, set **X** to a numeric column (e.g. `Chemical Space X`)
   and add at least one **Y** column (e.g. `Chemical Space Y`).
4. In the property panel under **Data**, set **Split** to a categorical column that has
   3–5 distinct categories (e.g. `Stereo Category`).
5. Confirm the legend renders one colored entry per category.
6. Assign a known per-category baseline color to each category (e.g. R_ONE=red,
   S_ABS=green, S_ACHIR=blue, S_PART=yellow, S_UNKN=magenta) so the github-1498
   color guard in Scenario 1 can assert an exact value rather than a default.

## Scenarios

### Scenario 1: Legend click-to-filter preserves original line colors

Steps:
1. Record the initial `df.filter.trueCount` (full dataset, all categories visible).
2. Note the color assigned to one category in the legend (read from the viewer's
   `colorColumnName`-related property, or a property-panel color swatch; this is
   the baseline for the color guard).
3. Click a **different** legend category entry to filter — only rows for that one
   category should remain visible.
4. Wait for the chart to re-render (no page error expected).
5. Assert: `df.filter.trueCount` is strictly less than the value recorded in Step 1.
6. Assert: the color property for the **remaining** category line is unchanged from
   Step 2 — it has NOT changed to blue or any other default color (github-1498
   regression guard).
7. Click the active legend category again (or use Reset View / reset legend filter)
   to restore all categories.
8. Assert: `df.filter.trueCount` equals the value recorded in Step 1.
9. Assert: no console errors or page errors occurred during Steps 3–8.

Expected:
- Filtering via a legend click reduces df.filter.trueCount to only the rows
  belonging to the clicked category.
- The remaining visible line's color property is identical to its pre-filter value —
  legend click-to-filter must not silently overwrite the color assignment.
- Resetting the filter restores the full row count with no errors.

### Scenario 2: Color changes persist through layout save/reopen and project save/reopen

Steps:
1. With the split line chart from Setup still active, change the color of one
   category in the legend (right-click the legend entry → Change Color, or use the
   property panel color picker).
2. Record the new color value (hex or name) as the expected persisted value.
3. Open **Save Layout** (via the viewer or table view toolbar).
4. Save the layout with a distinct name (e.g. `linechart-color-persist-test`).
5. Close all open tables and views.
6. Apply the saved layout (File → Layouts → open saved layout).
7. Wait for the chart to render.
8. Assert: the per-category color property reads back the same value recorded in
   Step 2 — color is restored from the layout round-trip (GROK-17278 regression
   guard; save→close→reopen is a real signal, not a setter re-read).
9. With the split Line Chart still active and the category color set, click the
   **SAVE** ribbon button, name the project distinctly, and confirm (**OK**).
   A "Share …" dialog appears after the save — dismiss it (**Cancel**).
   (Only the SAVE button captures the viewer layout — a JS-API
   `Project.create().addChild(saveLayout())` throws "Unable to add entity to the
   project", and `addChild(df)` alone restores the table but not the viewer.)
10. Close all. Reopen the saved project.
11. Locate the Line Chart viewer in the reopened project layout.
12. Assert: the Line chart viewer is restored and the per-category color property
    reads back the same value from Step 2 — color is restored from the project
    round-trip (GROK-17278 regression guard).
13. Assert: the markers legend is present in the reopened viewer (GROK-19825
    regression guard — legend must not disappear on project reopen).
14. Assert: no console errors or page errors from the spec's own steps (the
    platform's post-save Share dialog may log a benign NullError, which is
    excluded from the spec's error gate).

Expected:
- After saving the layout and reopening, the color applied to the category in Step 1
  is unchanged — saved layouts must serialize per-category color assignments.
- After saving the project, closing all, and reopening, the same color is present
  and the markers legend is visible — project serialization preserves both color
  customizations and the legend itself.
