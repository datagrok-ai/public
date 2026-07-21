---
feature: barchart
realizes_atlas:
  - barchart.cp.sort-and-orient
realizes:
  - viewers.bar-chart
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-19480
    status: fixed
  - id: github-3417
    status: regression-risk
realized_as:
  - barchart-sort-orient-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      With Orientation Vertical and Bar Sort Descending, bars render vertically
      with heights decreasing left to right by aggregated value.
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      With stacking enabled on the negative-sum value column, the stacked bars
      keep rendering and the stack legend appears as a result of the stack
      action — absent before it, visible after (GROK-19480).
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      With negative-sum categories present, the vertical bars render without
      error, negatives below the zero baseline (GROK-19480).
  - anchor: "Scenario 1 Step 8"
    expectation: >-
      Reverting Orientation to Horizontal and Bar Sort to Ascending renders
      horizontal ascending bars matching the baseline.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      Switching to Vertical on small-magnitude counts does not leave excess
      whitespace above the tallest bar (github-3417).
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      Resetting Orientation to Auto returns the chart to the default layout with
      no residual whitespace artifact.
---

# Bar Chart — Sorting and Orientation

## Setup

1. Close all open tables and viewers.
2. Open spgi-100 — a 100-row SPGI sample (`System:AppData/Chem/tests/spgi-100.csv`).
3. Add a Bar Chart viewer to the current table view.
4. In the in-chart **Category** selector, set the Split column to **Primary Series names**.
5. In the Context Panel > Data section, set **Value** to **Chemical Space X** and **Value Aggr Type** to **sum**.

## Scenarios

### Scenario 1: Vertical orientation with Descending sort, stacking on negative sums, and negative values

Steps:
1. In the Context Panel > Style section, set **Orientation** to **Vertical** (non-default; default is Horizontal).
2. In the Context Panel > Style section, set **Bar Sort** to **Descending** (non-default; default is Ascending).
3. Verify that bars render vertically and that bar heights decrease from left to right, confirming
   the descending sort order by aggregated **Chemical Space X** sum is applied under vertical orientation.
4. In the Context Panel > Style section, set **Legend Visibility** to **Always** and enable stacking
   (set **Stack** to **Stereo Category**). Verify that stacking holds on the negative-sum
   **Chemical Space X** aggregation: the stack legend is visible in the viewer area (GROK-19480
   regression guard — stacking used to break for negative aggregated values). Then clear the Stack
   column to restore the un-stacked baseline for the following steps.
5. Confirm that the **Chemical Space X** column contains categories with negative sum totals.
   Verify that the chart renders correctly with negative-value bars displayed below the zero baseline
   without breaking the layout or producing a JavaScript error (GROK-19480 regression guard).
6. In the Context Panel > Style section, set **Bar Sort** to **Ascending** (revert sort).
7. Verify that bars are now in ascending order from left to right under vertical orientation.
8. In the Context Panel > Style section, set **Orientation** to **Horizontal** (revert orientation).
   Verify that bars render horizontally in ascending order; the layout returns to the expected
   horizontal baseline.

Expected:
- After Step 3: bars are vertical with decreasing heights left-to-right; values and colors are
  rendered correctly for all displayed categories.
- After Step 4: enabling the Stack column on the negative-sum Chemical Space X aggregation keeps
  stacking functional — the stack legend appears as a result of the stack action (it is absent
  before it) and no errors are raised (GROK-19480); clearing the Stack column restores the
  un-stacked baseline.
- After Step 5: chart remains functional with no errors despite negative Chemical Space X sums;
  negative bars appear below the baseline.
- After Step 8: chart reverts to horizontal bars in ascending order; no residual style artifact
  from vertical orientation remains.

### Scenario 2: Orientation round-trip and whitespace regression guard (github-3417)

Steps:
1. In the Context Panel > Data section, change **Value** to **CAST Idea ID** and **Value Aggr Type** to **count**
   (small-magnitude values, suitable for triggering the whitespace regression).
2. In the Context Panel > Style section, set **Orientation** to **Vertical**.
3. Observe the bar chart under vertical orientation with small-magnitude count values.
   Verify that available vertical space is used efficiently — the max bar height fills a reasonable
   proportion of the viewer area without a large empty whitespace gap above the tallest bar
   (github-3417 regression guard: excess whitespace with vertical orientation and small values).
4. In the Context Panel > Style section, set **Orientation** to **Horizontal** (intermediate revert).
   Verify bars render horizontally with no layout artifact.
5. In the Context Panel > Style section, set **Orientation** to **Auto** (full revert to default).
   Verify the chart returns to the auto-orientation layout with no residual whitespace artifact
   from the previous vertical-orientation state.

Expected:
- After Step 3: vertical bars occupy a reasonable fraction of the viewer's height; no large blank
  area appears above the tallest bar for small-count data.
- After Step 4: horizontal bars render correctly; layout is clean.
- After Step 5: orientation is back to the Auto (default) state; chart layout is indistinguishable
  from the initial state.
