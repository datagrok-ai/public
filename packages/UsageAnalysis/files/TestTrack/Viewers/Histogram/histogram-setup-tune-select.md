---
feature: histogram
realizes_atlas:
  - histogram.cp.setup-tune-select
realizes:
  - viewers.histogram
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-19759
    status: fixed
  - id: GROK-18223
    status: fixed
  - id: github-2296
    status: fixed
expected_results:
  - anchor: "Scenario 1 Step 3"
    expectation: "Per-bin value counts are rendered on top of each bin (Show Values on).
      Actuation: the count glyphs are canvas-drawn text, not headless-readable, so
      the spec asserts the driving prop (showValues === true) plus the no-error
      render floor rather than the label pixels."
  - anchor: "Scenario 1 Step 4"
    expectation: "Increasing bin count to 50 does not crash or produce an error."
  - anchor: "Scenario 1 Step 5"
    expectation: "Bin count shown in the properties panel matches the in-plot bin
      selector (both read 50 after the in-plot change)."
  - anchor: "Scenario 1 Step 6"
    expectation: "Bin count shown in the in-plot selector matches the properties
      panel after setting via the panel (both read 30)."
  - anchor: "Scenario 1 Step 7"
    expectation: "Value column shown in the in-plot column selector matches the
      properties panel Value field (both read HEIGHT after the in-plot change)."
  - anchor: "Scenario 1 Step 8"
    expectation: "Value column shown in the properties panel matches the in-plot
      selector (both read WEIGHT after the panel change)."
  - anchor: "Scenario 1 Step 9"
    expectation: "Clicking a bin selects the rows in that bin; the status bar shows
      a non-zero selected count matching the bin height."
  - anchor: "Scenario 1 Step 10"
    expectation: "The current-row indicator is functional and hovering a bin
      raises no error: df.currentRowIdx clears to -1 and round-trips to a set
      row whose value lies within the value column's range. Actuation: a
      headless canvas hover drives NEITHER df.mouseOverRowIdx (stays -1) NOR
      df.currentRowIdx (recon 2026-07-21, demog.csv), so the hover itself is
      covered by the no-error floor and the indicator is exercised through its
      own set/read path."
  - anchor: "Scenario 2 Step 2"
    expectation: "Show Values is off; per-bin count labels are no longer rendered.
      Actuation: the labels are canvas-drawn text, not headless-readable, so the
      spec asserts the driving prop (showValues === false) plus the no-error
      render floor rather than the absence of label pixels."
  - anchor: "Scenario 2 Step 3"
    expectation: "Bin count reverts to 20; the properties panel and in-plot selector
      both show 20."
  - anchor: "Scenario 2 Step 4"
    expectation: "Value column reverts to AGE in both the properties panel and the
      in-plot selector."
realized_as:
  - histogram-setup-tune-select-spec.ts
---

# Histogram — Core setup, tuning, and bin selection

## Setup

1. Close all open tables and viewers.
2. Open `System:DemoFiles/demog.csv`.
3. Add a Histogram viewer to the table view.

## Scenarios

### Scenario 1: Set value column, tune bin count with Show Values, verify sync, click bin to select

Steps:
1. In the properties panel (Context Panel > Value), set **Value** to `AGE`.
2. In the properties panel (Context Panel > Value), enable **Show Values**.
3. Verify that per-bin value counts are rendered on top of each bin.
4. Set **Bins** to `50` via the in-plot bin count selector (the numeric stepper shown inside the viewer). Verify no crash or error occurs (GROK-18223).
5. Verify that the **Bins** field in the properties panel (Context Panel > Value > Bins) also reads `50` — the in-plot change must propagate to the panel (github-2296).
6. Set **Bins** to `30` via the properties panel. Verify that the in-plot bin selector also reads `30` — the panel change must propagate in-plot (github-2296).
7. In the in-plot column selector, switch the value column to `HEIGHT`. Verify that the properties panel **Value** field updates to `HEIGHT` (github-2296).
8. In the properties panel, set **Value** back to `WEIGHT`. Verify that the in-plot column selector updates to `WEIGHT` (github-2296). Per-bin counts displayed via Show Values should update to reflect the WEIGHT distribution (GROK-19759).
9. Click a bin in the histogram. Verify that the rows belonging to that bin become selected and the status bar shows a non-zero selected count equal to the bin height.
10. Hover over the selected bin (no error expected) and verify the current-row indicator is functional: it clears and round-trips to a set row whose value is within the value column's range. (A headless hover drives neither the current-row nor the mouse-over indicator, so those visual cues are a human-side check.)

Expected:
- Per-bin value counts render on top of bins when Show Values is on (GROK-19759).
- Increasing the bin count to 50 does not crash (GROK-18223).
- Bin count and value column stay synchronized between the in-plot controls and the properties panel in both directions (github-2296).
- Clicking a bin selects its rows; the status bar selected count matches the bin height.

### Scenario 2: Revert Show Values, bin count, and value column — round-trip

Steps:
1. Set **Value** to `AGE` via the properties panel.
2. In the properties panel, disable **Show Values**. Verify per-bin count labels are no longer rendered.
3. Set **Bins** to `20` via the properties panel (the default). Verify both the properties panel and the in-plot selector show `20`.
4. Verify the value column reads `AGE` in both the properties panel and the in-plot column selector.

Expected:
- Show Values off removes per-bin count labels.
- Bin count and value column return to the baseline values (20 bins, AGE column) in both the properties panel and in-plot controls — confirming a clean round-trip.
