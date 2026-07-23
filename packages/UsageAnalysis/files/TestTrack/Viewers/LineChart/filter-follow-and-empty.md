---
feature: linechart
realizes_atlas:
  - linechart.cp.filter-follow-and-empty
realizes:
  - viewers.line-chart
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - filter-follow-and-empty-spec.ts
related_bugs:
  - id: github-2574
    status: fixed
  - id: GROK-18375
    status: fixed
  - id: GROK-20185
    status: fixed
expected_results:
  - anchor: "S1 Step 4: filter to zero rows, filter.trueCount === 0"
    expectation: >-
      df.filter.trueCount equals 0 — filtering all rows out leaves zero visible
      rows, consistent with the applied filter.
  - anchor: "S1 Step 5-6: hover empty log-axis chart raises no error (github-2574)"
    expectation: >-
      Hovering the empty chart area produces no new console errors or page
      errors (github-2574 regression guard — empty log-axis chart must not throw
      on hover).
  - anchor: "S1 Step 7-8: remove filter restores baseline row count"
    expectation: >-
      After removing the zero-rows filter, the filtered row count
      (df.filter.trueCount) returns to the full-dataset baseline recorded in
      Setup — the empty state is cleanly reversible.
  - anchor: "S2 Step 1-2: set X to a date column with Year-quarter time-split, no error"
    expectation: >-
      Setting X to the date column with the Year-quarter time split re-renders
      the chart with no new console or page errors — the GROK-18375
      preconfiguration (time-split X) is in place.
  - anchor: "S2 Step 5: filter while time-split X keeps rows > 0 (GROK-18375)"
    expectation: >-
      df.filter.trueCount is strictly greater than 0 — applying a structure
      filter while the X axis is time-split must not reduce the count to zero
      (GROK-18375 regression guard).
  - anchor: "S2 Step 7: remove filter, restore numeric X for next scenario"
    expectation: >-
      After removing the categorical filter and restoring the numeric X column
      with the time split cleared, the filtered row count returns to the
      baseline — the Scenario 2 state does not leak into Scenario 3.
  - anchor: "S3 Step 4-5: narrow X range, filter.trueCount drops below baseline (GROK-20185)"
    expectation: >-
      df.filter.trueCount changes from its initial value as the narrowed X
      range is applied through the Filter Panel — the live filter update
      propagates to the dataframe state.
  - anchor: "S3 Step 4-5: narrow X range"
    expectation: >-
      No new console errors or page errors are raised while the X range is
      narrowed and applied through the Filter Panel (GROK-20185 regression
      guard).
  - anchor: "S3 Step 7-8: restore full range, filter.trueCount returns to baseline"
    expectation: >-
      After restoring the slider to full range, df.filter.trueCount returns to
      the original full-dataset value — the filter is cleanly reversible.
---

# Line Chart — Filter interaction, follow-filter, and empty-chart resilience

## Setup

1. Open the spgi-100 dataset.
2. Add a Line Chart viewer to the table view (via the Toolbox).
3. In the viewer property panel, set **X** to a numeric column (e.g.
   `Chemical Space X`) and add at least one **Y** column (e.g.
   `Chemical Space Y`).
4. Enable **Axes Follow Filter** in the viewer property panel (Data section).
5. Open the **Filter Panel** (via the toolbar).
6. Record the initial filtered row count (all rows visible — baseline).

## Scenarios

### Scenario 1: Filter to zero rows and hover empty log-axis chart

Steps:
1. In the viewer property panel, switch the **X axis scale** to **Logarithmic**.
2. In the Filter Panel, set a numeric filter on the X column with a range that
   matches NO rows (e.g. set the max below the column minimum, such that the
   predicate is impossible to satisfy).
3. Wait for the chart to re-render (no page error expected during the update).
4. Verify the filtered row count equals `0`.
5. Move the mouse pointer over the empty chart canvas and pause briefly (simulate
   a hover without clicking).
6. Verify no new console errors or page errors appeared since the filter was
   applied (github-2574 regression guard — hovering an empty log-axis chart
   must not throw).
7. Reset the filter in the Filter Panel (remove the constraint) to restore
   full visibility.
8. Verify the filtered row count returns to the baseline recorded in Setup Step 6.
9. Switch the X axis scale back to **Linear**.

Expected:
- Filtering to zero visible rows with a logarithmic X axis is handled without
  errors; the filtered row count reads exactly `0`.
- Hovering the empty chart area raises no console or page errors — the
  empty-state hit-test path is defensive (github-2574 guard).
- Removing the filter fully restores the row count.

### Scenario 2: Time-split X axis with structure filter retains rows

Steps:
1. In the viewer property panel, set **X** to a date/time column (e.g. `Date`)
   and configure the time split: set **X Map** (Context Panel > X — the selector
   appears when X is a datetime column) to `year quarter`.
2. Confirm the chart re-renders with time-bucketed X values and no error.
3. In the Filter Panel, add a **Structure** filter on a SMILES column
   (e.g. `Smiles`) with a non-empty substructure query (e.g. a simple ring
   fragment that is present in several rows of spgi-100).
4. Wait for the filter to propagate.
5. Verify the filtered row count is strictly greater than `0` — a structure
   filter applied while the X axis uses time-split must NOT reduce all rows to
   zero (GROK-18375 regression guard).
6. Verify no console errors or page errors occurred during Steps 2–5.
7. Remove the structure filter from the Filter Panel and reset the X column to the
   numeric column used in Setup (restore for the next scenario).

Expected:
- The structure filter coexists with the time-split X axis setting and leaves a
  non-empty row set; GROK-18375 manifested as every row being incorrectly
  filtered out in this combination.
- No errors are raised during the structure filter application.

### Scenario 3: Narrowing the X range in the Filter Panel updates the chart live

Steps:
1. Confirm the X column is set to the numeric column from Setup and the axis
   scale is **Linear**.
2. Locate the X-column filter in the Filter Panel.
3. Record the current filtered row count (baseline — all rows).
4. Narrow the X range to roughly the middle 50% of the column's value span by
   applying the range through the Filter Panel filter.
5. Verify the filtered row count is strictly less than the baseline from Step 3 —
   the live filter update propagates to the dataframe state (GROK-20185
   regression guard; the narrowing must not leave the count unchanged).
6. Verify no new console errors or page errors occurred while the range was
   narrowed and applied.
7. Restore the full range through the same Filter Panel filter.
8. Verify the filtered row count returns to the baseline value from Step 3.

Expected:
- Narrowing the X range through the Filter Panel changes the filtered row count
  — the chart and the dataframe filter state update live (GROK-20185 guard).
- No errors are raised while the range is narrowed, applied, or restored.
- Restoring the full range brings the row count back to the original baseline.

## Automation notes

Setup: the dataset is opened via the `readDataframe` helper; the viewer is
located via the `findViewer` helper; the Filter Panel is opened via
`grok.shell.tv.getFiltersGroup()`; "filtered row count" is read as
`df.filter.trueCount`.

Scenario 1: the empty-chart state (Step 2, "matches NO
rows") is driven through a CATEGORICAL Filter Panel filter on the `Series`
column with an empty selection, not through a numeric X-range filter. A numeric
Filter Panel filter is an embedded Histogram whose range strip clamps the
selection to at least one bucket (`df.filter.trueCount` never reaches 0), so it
cannot produce the exactly-zero-rows state github-2574 needs; a categorical
filter with `selected: []` yields `trueCount === 0` deterministically. The
observable — an empty log-axis chart that must not throw on hover — is exercised
identically. Reaching zero rows via a hand-dragged numeric range stays a
human-side variation.

Scenario 2: the Structure/SMILES filter (Steps 3–5) is driven as a
CATEGORICAL Filter Panel filter on the structure-family `Series` column — the
Chem substructure filter cannot be driven headless. The categorical stand-in
exercises the same GROK-18375 path (a filter applied while the X axis is
time-split must not remove every row); a real substructure query stays a
human-side variation of this scenario.

Scenario 3: the range-slider DRAG is not scriptable in automation — the
Filter Panel numeric filter is an embedded Histogram whose range
strip is canvas-drawn with no DOM handle elements, canvas-strip drags leave the
filter untouched, and the min/max text inputs stay invisible to automation. The
range narrowing is therefore applied programmatically through the Filter
Panel's filter (`getFiltersGroup().updateOrAdd`), which still exercises the
GROK-20185 observable: the live `df.filter.trueCount` update and its
reversibility. A manual drag check stays a human-side variation of this
scenario.
