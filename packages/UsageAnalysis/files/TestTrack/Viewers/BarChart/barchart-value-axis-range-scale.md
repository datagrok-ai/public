---
feature: barchart
realizes_atlas:
  - barchart.cp.value-axis-range-scale-scroll
  - barchart.int.selected-filtered-rows-need-cumulative-aggr
realizes:
  - viewers.bar-chart
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - barchart-value-axis-range-scale-spec.ts
related_bugs: []
expected_results:
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      With Value Min above the shortest bar, bars below it are clipped with a
      clipped-bar indicator at the base (GROK-19346).
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      With Value Max below the tallest bar, bars above it are clipped at the top
      and the value axis stops at the maximum.
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      The value-axis scroll bar navigates the constrained range without error.
  - anchor: "Scenario 1 Step 9"
    expectation: >-
      With a logarithmic value axis, positive-count bars re-scale and render
      without error.
  - anchor: "Scenario 1 Step 11"
    expectation: >-
      Switching the value axis back to linear restores linear proportional bar
      heights.
  - anchor: "Scenario 1 Step 13"
    expectation: >-
      Clearing Value Min and Value Max restores the full-range axis with no
      clipping.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      With Show Selected Rows enabled under count, a selection drives a per-bar
      overlay segment.
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      With Show Filtered Rows enabled under count, a filter drives a per-bar
      overlay segment.
  - anchor: "Scenario 2 Step 8"
    expectation: >-
      Under a non-cumulative aggregation (avg), the Selected/Filtered Rows
      overlays are suppressed.
---

# Bar Chart — Value-Axis Range, Scale, and Selected/Filtered Rows Overlays

## Setup

1. Close all open tables and viewers.
2. Open spgi-100 — a 100-row SPGI sample (`System:AppData/Chem/tests/spgi-100.csv`).
3. Add a Bar Chart viewer to the current table view.
4. In the in-chart **Category** selector, set the Split column to **Primary Series names**.
5. In the Context Panel > Data section, set **Value** to **CAST Idea ID** and **Value Aggr Type** to **count**.
   (Count of CAST Idea ID per category — positive-valued aggregation with sufficient range variation.)

## Scenarios

### Scenario 1: Value Min / Value Max clipped-bar indicators, scroll bar, and log-scale round-trip

Steps:
1. Observe the Bar Chart viewer with all bars fully rendered at the full data range.
   Note the approximate count of the shortest and tallest bars visible.
2. In the Context Panel > Value section, set **Value Min** to a value above the shortest bar's count
   (e.g. a value that clips at least one bar on the low end).
3. Verify that bars whose aggregated count falls below the specified **Value Min** are clipped:
   a clipped-bar indicator (a small notch or arrow) is rendered at the base of each affected bar
   (GROK-19346 regression guard: value-axis min/max constraints with clipped-bar indicators).
4. In the Context Panel > Value section, set **Value Max** to a value below the tallest bar's count
   (e.g. a value that clips at least one bar on the high end).
5. Verify that bars whose aggregated count exceeds the specified **Value Max** are clipped:
   a clipped-bar indicator is rendered at the top of each affected bar; the value axis ends at the
   specified maximum.
6. Confirm the value axis now shows a constrained range (Value Min to Value Max).
7. Locate the value-axis scroll bar (appears when the data range extends beyond the displayed axis
   window). Drag it to navigate within the constrained range.
   Verify the scroll bar responds — the chart shifts to reveal different portions of the value axis
   range without error or visual corruption.
8. In the Context Panel > Value section, set **Value Scale** to **log** (non-default; default is linear).
9. Verify that the chart re-renders with a logarithmic value axis: bars with positive non-zero counts
   are displayed with logarithmically-spaced heights; axis labels show logarithmic intervals (e.g. 1,
   10, 100); no error occurs and no bars disappear.
10. Confirm bars that had clipped-bar indicators still appear clipped under the log scale if applicable.
11. In the Context Panel > Value section, set **Value Scale** back to **linear** (revert from log).
    Verify bars return to linear proportional heights and axis labels reflect linear spacing; the chart
    is visually equivalent to the state before Step 8.
12. Clear **Value Min** (set to empty / remove) and clear **Value Max** (set to empty / remove).
13. Verify that all bars are now fully rendered without clipped-bar indicators; the value axis spans
    the full aggregated data range; the chart matches the initial baseline state from Step 1.

Expected:
- After Step 3: at least one bar shows a clipped-bar indicator at its base; no error is thrown.
- After Step 5: at least one bar shows a clipped-bar indicator at its top; the axis stops at Value Max.
- After Step 7: the scroll bar is present and dragging it shifts the visible axis range without error.
- After Step 9: logarithmic axis labels are visible; all positive-count bars are rendered correctly.
- After Step 11: linear axis with the same bars as before the log-scale switch; round-trip is clean.
- After Step 13: no clipped-bar indicators; full-range axis; chart is back to the baseline layout.

### Scenario 2: Show Selected Rows and Show Filtered Rows overlays require cumulative aggregation

Steps:
1. Confirm the current **Value Aggr Type** is **count** (cumulative aggregation from Setup Step 5).
2. In the Context Panel > Selection section, enable **Show Selected Rows**.
3. Select a subset of rows in the grid (e.g. click a row to select it; Ctrl+click to add more).
   Verify that an overlay bar segment representing the selected rows appears within each category bar
   that contains selected rows — the overlay is visible on top of the base bar.
4. In the Context Panel > Selection section, enable **Show Filtered Rows**.
5. Open the Filter panel and apply a filter to reduce the visible row count
   (e.g. filter the dataset so that at least one category has fewer filtered rows than total rows).
   Verify that an overlay bar segment representing the filtered rows appears within the relevant
   category bars in the Bar Chart.
6. Disable the filter (clear or remove it) and deselect all rows (press Esc).
7. Disable **Show Selected Rows** and **Show Filtered Rows** to restore baseline.
8. In the Context Panel > Data section, change **Value Aggr Type** to **avg** (non-cumulative aggregation).
   Enable **Show Selected Rows** and **Show Filtered Rows** again.
   Select a subset of rows; apply a filter.
9. Verify that no overlay bar segments appear for either selected or filtered rows when **avg** aggregation
   is active — the overlays are suppressed for non-cumulative aggregations
   (overlays only computed for count / sum / value-count).
10. Change **Value Aggr Type** back to **count** and disable both overlay options to restore baseline.

Expected:
- After Step 3: overlay segments for selected rows are visible within the relevant category bars.
- After Step 5: overlay segments for filtered rows are visible within the relevant category bars.
- After Step 9: no overlay segments appear for either selected rows or filtered rows under avg aggregation;
  bars render as standard single-segment bars.
