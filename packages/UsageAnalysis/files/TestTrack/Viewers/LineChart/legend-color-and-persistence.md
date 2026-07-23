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
  - anchor: "S1: legend click filters + preserves remaining line colors"
    expectation: >-
      lc.filter.trueCount (the Line chart's viewer-local filter) is strictly
      less than the full row count — the legend click narrows the rows the
      viewer displays.
  - anchor: "legend click filters + preserves remaining line colors"
    expectation: >-
      The color assigned to the remaining category's line is unchanged from its
      pre-filter value (does not turn blue or shift; canvas render is not
      directly assertable, so assert via the persisted color property value).
  - anchor: "S1: re-click legend category resets filter to full count"
    expectation: >-
      After resetting the legend filter, lc.filter.trueCount (the viewer-local
      filter) is restored to the full row count (all categories visible again).
  - anchor: "S2: category color persists through layout round-trip (GROK-17278)"
    expectation: >-
      After saving the layout and reopening from it, the per-category color
      property reads back the same value that was set before saving —
      persistence round-trip for layout (GROK-17278 regression guard).
  - anchor: "S2 Steps 9-13: color AND markers legend survive a project save/close/reopen via the SAVE button (GROK-17278, GROK-19825)"
    expectation: >-
      After saving the project via the SAVE button, closing all, and reopening
      the project, the Line chart viewer is restored and the per-category color
      reads back the value set before saving (GROK-17278 — project
      serialization preserves per-category color).
  - anchor: "markers legend survive a project save/close/reopen"
    expectation: >-
      The markers legend is present in the reopened project's Line chart viewer
      (GROK-19825 — the legend must not disappear on project reopen).
---

# Line Chart — Legend filter-color and layout/project persistence

## Setup

1. Open the spgi-100 dataset (use `readDataframe` helper for `spgi-100.csv`).
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
1. Record the full row count (full dataset, all categories visible — the
   baseline for the viewer-local filter checks below).
2. Note the color assigned to one category in the legend (read from the viewer's
   `colorColumnName`-related property, or a property-panel color swatch; this is
   the baseline for the color guard).
3. Click a **different** legend category entry to filter — only rows for that one
   category should remain visible.
4. Wait for the chart to re-render (no page error expected).
5. Verify the viewer's own displayed-row count (its viewer-local filter) is
   strictly less than the value recorded in Step 1 — the legend click narrows
   the viewer's own filter, not the table's global filter.
6. Verify the color property for the **remaining** category line is unchanged from
   Step 2 — it has NOT changed to blue or any other default color (github-1498
   regression guard).
7. Click the active legend category again (or use Reset View / reset legend filter)
   to restore all categories.
8. Verify the viewer-local filter count equals the value recorded in Step 1.
9. Verify no console errors or page errors occurred during Steps 3–8.

Expected:
- Filtering via a legend click narrows the viewer-local filter count
  to only the rows belonging to the clicked category.
- The remaining visible line's color property is identical to its pre-filter value —
  legend click-to-filter must not silently overwrite the color assignment.
- Resetting the filter restores the full row count with no errors.

### Scenario 2: Color changes persist through layout save/reopen and project save/reopen

Steps:
1. With the split line chart from Setup still active, change the color of one
   category in the legend (right-click the legend entry → Change Color).
2. Record the new color value (hex or name) as the expected persisted value.
3. Save the layout: top menu **View > Layout > Save to Gallery**.
4. Confirm the layout was saved (a balloon reports it; the layout appears on the
   context panel).
5. Close the Line Chart viewer (the table view stays open — a layout needs an
   open table view to apply to).
6. Apply the saved layout (**View > Layout > Open Gallery**, click the saved
   layout) — the layout restores the closed viewer with the settings it had
   before closing.
7. Wait for the chart to render.
8. Assert: the per-category color property reads back the same value recorded in
   Step 2 — color is restored from the layout round-trip (GROK-17278 regression
   guard; save→close→reopen is a real signal, not a setter re-read).
9. With the split Line Chart still active and the category color set, click the
   **SAVE** ribbon button, name the project distinctly, and confirm (**OK**).
   A "Share …" dialog appears after the save — dismiss it (**Cancel**).
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

## Automation notes

Scenario 1: the "viewer-local filter" counts are read as `lc.filter.trueCount`
on the Line chart's own filter (distinct from the dataframe's global
`df.filter`); the dataset is opened via the `readDataframe` helper and the
viewer handle obtained via the `findViewer` helper.

Scenario 2, steps 5-6: the spec does not close the viewer — it clears the
per-category color (`colors.setCategorical({})`) as the stand-in perturbation,
then re-applies the saved layout onto the open table view and asserts the color
is restored. Closing the viewer by hand (step 5) exercises the same
restore-from-layout signal.
