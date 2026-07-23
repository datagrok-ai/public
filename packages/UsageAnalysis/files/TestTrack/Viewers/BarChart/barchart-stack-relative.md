---
feature: barchart
realizes_atlas:
  - barchart.cp.stack-relative-negatives
realizes:
  - viewers.bar-chart
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-19480
    status: fixed
  - id: github-2659
    status: fixed
realized_as:
  - barchart-stack-relative-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      With Relative Values enabled and a Stack column set, each bar normalizes
      to equal width and splits into stacked segments — the width normalization
      produces a large canvas color delta and the chart is not blank
      (github-2659).
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      At least two stacked segments per bar are visible.
  - anchor: "Scenario 1 Step 7"
    expectation: >-
      With a sum aggregation whose totals include negative values, the stacked
      bars render without error (GROK-19480).
  - anchor: "Scenario 1 Step 9"
    expectation: >-
      Disabling Relative Values reverts bars to absolute widths — the
      de-normalization produces a large canvas color delta.
  - anchor: "Scenario 1 Step 11"
    expectation: >-
      Removing the Stack column collapses bars to single-segment with no legend
      column.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      With no Stack column, enabling Relative Values has no effect — bars keep
      absolute widths (not normalized), no legend, no error.
  - anchor: "Scenario 2 Step 5-6"
    expectation: >-
      Setting a Stack column while Relative Values is already on activates the
      normalized stacked view immediately — bars normalize to equal width and
      split into stacked segments (large canvas color delta), legend renders.
  - anchor: "Scenario 2 Step 7-8"
    expectation: >-
      Removing the Stack column reverts bars to absolute single-segment widths
      (Relative Values inert again, legend gone); re-adding the Stack column
      re-activates the normalized stacked view — the toggle is repeatable.
  - anchor: "Scenario 2 Step 9-10"
    expectation: >-
      Disabling Relative Values with no Stack column restores the baseline
      single-segment absolute-width bars, no error or blank chart.
---

# Bar Chart — Stacking, Relative Values, and Negative Aggregates

## Setup

1. Close all open tables and viewers.
2. Open spgi-100 — a 100-row SPGI sample (`System:AppData/Chem/tests/spgi-100.csv`).
3. Add a Bar Chart viewer to the current table view.
4. In the in-chart **Category** selector, set the Split column to **Primary Series Name**.
5. In the Context Panel > Data section, set **Value** to **Chemical Space X** and **Value Aggr Type** to **sum**.

## Scenarios

### Scenario 1: Relative Values with stacking and negative aggregates

Steps:
1. In the in-chart **Stack** selector, set the Stack column to **Scaffold Names** (non-default; no stack by default).
2. In the Context Panel > Data section, enable **Relative Values** (non-default; off by default).
3. Observe the Bar Chart viewer.
4. Verify that each outer bar normalizes to equal width and splits into stacked segments
   reflecting the proportional share of each **Scaffold Names** category within that **Primary Series Name** group.
   The chart must not be blank or empty (github-2659 regression guard: Relative Values + Stack must not blank the chart).
5. Verify that at least two stacked segments per bar are visible and that the segment proportions
   vary across the **Primary Series Name** categories.
6. In the Context Panel > Data section, confirm **Value Aggr Type** is **sum** and that the
   underlying **Chemical Space X** column contains negative values (mixed-sign numeric column).
7. Verify that the stacked bar chart renders without error even with negative sum totals —
   negative-value segments are handled correctly and no JavaScript errors appear in the console
   (GROK-19480 regression guard: stacked bars must not break on negative aggregated values).
8. In the Context Panel > Data section, disable **Relative Values** (revert to off).
9. Verify that bars revert to absolute single-segment widths reflecting the raw summed value
   of **Chemical Space X** per category; bar heights now differ across categories proportionally
   to their absolute sums.
10. In the in-chart **Stack** selector, remove the Stack column (set to **None**).
11. Verify that bars collapse to single-segment (no stacked breakdown) with no legend column
    shown; the chart matches a standard single-bar histogram layout.

Expected:
- After Step 4: each outer bar is equal width and divided into stacked segments; the chart is not blank.
- After Step 5: at least two distinct stacked segments per bar are visible with varying proportions.
- After Step 7: the chart renders correctly with no errors despite negative Chemical Space X sums.
- After Step 9: bars show absolute heights with no normalization; Relative Values is off.
- After Step 11: single-segment bars, no legend column, standard layout.

### Scenario 2: Relative Values requires a Stack column — inert alone, activates when Stack is added

Steps:
1. Confirm the Stack column is **None** (from the end state of Scenario 1).
2. In the Context Panel > Data section, enable **Relative Values**.
3. Verify that bars remain at absolute widths — Relative Values has no visual effect when
   no Stack column is set: no stacked segments, no legend, no error. The option is accepted
   but resolves to a no-op because the stacking precondition is unmet.
4. Confirm **Stack** is still **None**.
5. With **Relative Values** still enabled, set the Stack column to **Scaffold Names**.
6. Verify the chart immediately switches to the normalized stacked view — each outer bar
   normalizes to equal width and splits into stacked segments, and the legend renders.
7. Remove the Stack column (set to **None**). Verify bars revert to single-segment
   absolute-width rendering — Relative Values is still on but inert again, and the
   legend disappears.
8. Re-add the Stack column (**Scaffold Names**). Verify bars again normalize to equal
   width with a legend — the activation is repeatable, not a one-time event.
9. Remove the Stack column and disable **Relative Values** (back to defaults).
10. Verify the chart shows single-segment absolute-width bars — the baseline state,
    with no error or blank chart at any point of the toggle sequence.

Expected:
- After Step 3: no stacked segments appear; bars are indistinguishable from the
  non-Relative-Values state — the viewer is not blank and does not error.
- After Step 6: bars normalize to equal width and split into stacked segments
  (a large repaint delta on the canvas); the legend lists the Stack categories.
- After Steps 7-8: removing the Stack column deactivates the normalized view and
  re-adding it re-activates it — the behavior round-trips repeatably.
- After Step 10: baseline single-segment bars at absolute widths, no error.
