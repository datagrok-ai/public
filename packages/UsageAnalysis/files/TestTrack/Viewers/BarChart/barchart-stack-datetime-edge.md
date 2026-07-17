---
feature: barchart
realizes_atlas:
  - barchart.int.stack-needs-additive-aggr
  - barchart.int.datetime-split-enables-map
realizes:
  - viewers.bar-chart
priority: p1
target_layer: playwright
coverage_type: edge
related_bugs: []
expected_results:
  - anchor: "Step 4"
    expectation: With a non-additive aggregation (avg) and a Stack column, no
      stacked segments appear and the legend column is suppressed.
  - anchor: "Step 5"
    expectation: Switching to sum (additive) with the same Stack column produces
      stacked bars with a visible legend.
  - anchor: "Step 7"
    expectation: Under count (additive), removing the Stack column collapses bars
      to single-segment with no legend column.
  - anchor: "Step 10"
    expectation: With a non-DateTime Split column, the Split Map selector is absent
      or disabled.
  - anchor: "Step 13"
    expectation: With a DateTime Split column (STARTED), the Split Map selector
      becomes available.
  - anchor: "Step 14"
    expectation: Changing Split Map from Year to Month re-categorizes bars by month.
  - anchor: "Step 15"
    expectation: Reverting Split Map to Year re-categorizes bars by year; the chart
      is not broken or blank.
realized_as:
  - barchart-stack-datetime-edge-spec.ts
---

# Bar Chart — Stack Aggregation Precondition and DateTime Split Map

## Setup

1. Close all open views.
2. Open `demog` (System:DemoFiles/demog.csv).
3. Add a Bar Chart viewer from the toolbox.
4. In the in-chart **Category** selector, set the Split column to **RACE**.

## Scenarios

### Scenario 1: Stacking requires additive aggregation

Steps:
1. In the in-chart **Stack** selector, set the Stack column to **SEX**.
2. In the **Value** selector, set the Value column to **AGE** with **avg** aggregation.
3. Observe the chart — verify that despite a Stack column being set, bars are rendered as single-segment (no stacked breakdown) and no stack legend column appears.
4. Verify no legend entries for the Stack column (SEX) are shown — the legend column is suppressed when the aggregation is non-additive (avg is not additive).
5. In the Context Panel under **Value**, change **Value Aggr Type** to **sum**. Verify the chart now shows stacked bars with visible stack-legend entries for each SEX category — confirming the stack table is only built for additive aggregations.
6. Change **Value Aggr Type** back to **count** (also additive). Verify stacked bars continue to appear.
7. Remove the Stack column (set to **None**). Verify bars collapse back to single-segment with no legend column.

Expected:
- With a non-additive aggregation (avg), no stacked segments appear and the legend column is suppressed — bars remain single-segment even though a Stack column is set.
- Switching to a SUM (additive) aggregation with the same Stack column produces stacked bars with a visible legend — confirming stacking requires an additive aggregation.
- With an additive aggregation (count) and a Stack column, removing the Stack column collapses bars back to single-segment — no legend column is shown.

### Scenario 2: DateTime Split column enables the Split Map selector

Steps:
8. In the in-chart **Value** selector, set the Value column to **AGE** with **count** aggregation.
9. In the in-chart **Category** selector, confirm the Split column is currently **RACE** (a string column, not DateTime). Inspect the Context Panel under **Category** — verify the **Split Map** (split function) selector is absent or disabled.
10. Verify the chart renders RACE categories without any date-map option available.
11. In the **Category** selector, change the Split column to **STARTED** (a DateTime column in demog.csv).
12. Inspect the Context Panel under **Category** — verify the **Split Map** selector is now available (enabled), offering options such as Year, Quarter, Month, Week, Day.
13. Verify the effective Split column name in the chart header/axis reflects the selected Split Map (e.g. "STARTED (Year)" by default).
14. Change **Split Map** to **Month**. Verify the chart re-categorizes bars by month — each bar represents a distinct month value; the chart is not blank or broken.
15. Change **Split Map** back to **Year**. Verify bars re-categorize by year and the chart displays correctly with the year-level grouping.

Expected:
- With a String or numeric (non-DateTime) Split column the Split Map selector is absent or disabled — no year/quarter/month map option is offered.
- With a DateTime Split column (STARTED), the Split Map selector becomes available — the effective split column name reflects the chosen map (e.g. Month, Quarter, Year).
- Changing Split Map from Year to Month re-categorizes the chart bars by month — categories reflect monthly grouping, not yearly.
- Reverting Split Map back to Year re-categorizes bars by year and the chart is not broken or blank.
