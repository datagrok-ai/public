---
feature: barchart
realizes_atlas:
  - barchart.int.selected-filtered-rows-need-cumulative-aggr
realizes:
  - viewers.bar-chart
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - barchart-selected-filtered-aggr-spec.ts
related_bugs: []
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      Under count (cumulative), selecting a category drives the Selected Rows
      overlay for that category; the filter is untouched.
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      Under count, filtering drives the Filtered Rows overlay; the selection is
      unchanged.
  - anchor: "Scenario 1 Step 8"
    expectation: >-
      Switching to min (non-cumulative) preserves selection and filter state
      with the overlays absent.
  - anchor: "Scenario 1 Step 10"
    expectation: >-
      Switching back to sum (cumulative) with a fresh selection and filter
      restores the overlays.
  - anchor: "Scenario 2 Step 4"
    expectation: >-
      Under value-count (cumulative), selecting two categories drives the
      Selected Rows overlay for their combined rows; the filter is unchanged.
  - anchor: "Scenario 2 Step 6"
    expectation: >-
      Under value-count (the `values` aggregation), filtering re-renders the
      Filtered Rows overlay (canvas color delta); the selection is unchanged.
  - anchor: "Scenario 2 Step 7"
    expectation: >-
      Switching to avg (non-cumulative) preserves selection and filter state
      with the overlays absent.
---

# Bar Chart — Selected / Filtered Rows Overlays with Cumulative Aggregations

## Setup

1. Close all open tables and viewers (use the platform's Close All action).
2. Open the demog dataset — a table view appears.
3. Add a Bar Chart viewer to the table view — the viewer appears in the layout.
4. In the Bar Chart property panel (Context Panel), set **Split** (category column)
   to `race`.
5. Set the **Value** column to `age`.

## Scenarios

### Scenario 1: Overlay presence toggles with aggregation type (count / min / sum)

Steps:
1. In the Bar Chart property panel, set **Value Aggr Type** to `count`
   (a cumulative aggregation).
2. Enable **Show Selected Rows** in the Bar Chart property panel.
3. Enable **Show Filtered Rows** in the Bar Chart property panel.
4. In the TableView grid, select a subset of rows belonging to one category
   (e.g. Ctrl+click rows where `race = 'Asian'` — Ctrl+click toggles a row's
   selection; Shift+drag selects a range).
   Verify: the Selected Rows bar segment (lighter overlay) appears inside the
   corresponding bar on the chart.
5. Open the Filter panel and add a filter on `sex` to restrict visible rows.
6. Verify: the Filtered Rows bar segment appears inside each bar, reflecting
   the filtered-row subset count.
7. In the Bar Chart property panel, change **Value Aggr Type** to `min`
   (a non-cumulative aggregation).
8. Verify: the Selected Rows and Filtered Rows overlays are absent from all bars
   — neither overlay segment is rendered for a non-cumulative aggregation.
9. Remove the `sex` filter and clear the row selection (press Esc or
   deselect all rows in the grid).
10. Change **Value Aggr Type** to `sum` (cumulative). Re-select a row subset
    and re-enable the filter.
    Verify: the Selected Rows and Filtered Rows overlays reappear inside the
    bars, confirming the overlay is restored for additive aggregations.
11. Disable **Show Selected Rows** and **Show Filtered Rows** in the property
    panel; clear selection and filters.

Expected:
- With count aggregation: Selected Rows and Filtered Rows overlays render inside
  each bar proportionally to the selected/filtered row count.
- Switching to min (non-cumulative) removes all overlays — bars show no
  selected/filtered segment.
- Switching back to sum (cumulative) restores both overlays.

### Scenario 2: Overlay absent for avg; present for value-count (boundary check)

Steps:
1. In the Bar Chart property panel, set **Value Aggr Type** to `value-count`
   (the `values` property literal; cumulative — counts non-null values).
2. Enable **Show Selected Rows** in the property panel.
3. Select several rows across different categories in the grid (Ctrl+click each row).
4. Verify: the Selected Rows overlay segment appears inside bars for categories
   that contain selected rows.
5. Enable **Show Filtered Rows**; in the Filter panel, add a filter to reduce visible rows to a subset.
6. Change the filter subset while **Value Aggr Type** is `values` (value-count).
   Verify: the Filtered Rows overlay re-renders inside the bars (the bar segments
   recolor in place) reflecting the new filtered-row counts under the `values`
   aggregation, before switching to avg.
7. Change **Value Aggr Type** to `avg` (non-cumulative). Keep the selection
   and filter active.
   Verify: the Selected Rows and Filtered Rows overlays are absent from all bars
   — no overlay segment is rendered for avg.
8. Disable **Show Selected Rows** and **Show Filtered Rows**; clear selection
   and filters to restore baseline state.

Expected:
- value-count is cumulative: overlays appear.
- avg is non-cumulative: overlays disappear, confirming the guard applies to all
  non-additive aggregation types, not only min.

## Automation notes

Setup: views are closed via `grok.shell.closeAll()`; the demog dataset is opened
via `readDataframe('demog.csv')`; the viewer handle is obtained via
`findViewer('Bar Chart', view)`.
