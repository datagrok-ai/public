---
feature: histogram
realizes_atlas:
  - histogram.cp.split-and-color
realizes:
  - viewers.histogram
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-19761
    status: fixed
  - id: GROK-18399
    status: fixed
expected_results:
  - anchor: "S1: Color property rows active (opacity 1.0) with no split"
    expectation: >-
      Color property-grid rows (tr[name='prop-color'],
      tr[name='prop-color-aggr-type'], tr[name='prop-invert-color-scheme']) are
      active — getComputedStyle(tr).opacity === '1' for all three.
  - anchor: "S1: split RACE disables Color rows (opacity 0.5)"
    expectation: >-
      After setting a Split Column, all three Color property-grid rows become
      disabled — getComputedStyle(tr).opacity === '0.5' (or < 1.0) for all
      three.
  - anchor: "S1: range narrow under split keeps filter valid, no error"
    expectation: >-
      While the split column is set and the range slider is adjusted,
      df.filter.trueCount stays in the valid range [0, df.rowCount] and no new
      console error is raised.
  - anchor: "S1: clear split re-activates Color rows (opacity 1.0 round-trip)"
    expectation: >-
      After clearing the Split Column, the three Color property-grid rows
      re-activate — getComputedStyle(tr).opacity returns to '1' for all three
      (round-trip).
realized_as:
  - histogram-split-and-color-spec.ts
---

# Histogram — Split vs Color-coding Transition

## Setup

1. Close all open tables and viewers.
2. Open the **demog** dataset (Files > Demo Files > demog.csv).
3. Add a Histogram viewer to the table view.
4. Click the histogram so it becomes the current object and the Context Panel shows its properties.
5. In the property panel (Context Panel), set the **Value Column** to `AGE`.

## Scenarios

### Scenario 1: Split column disables color-coding, range filter survives, color re-activates on split clear

**Purpose:** Assert the product-state transition when a split column is set and cleared — the three Color property-grid rows must disable (opacity 0.5) when a split column is active and re-enable (opacity 1.0) when it is cleared. The range slider while split is set must not error and must keep df.filter valid (GROK-18399). Color UI disabling is the observable product-state signal for GROK-19761; a splitColumnName set-then-read is NOT asserted (vacuous round-trip).

Steps:

1. With no split column set, open the Context Panel if not already shown and make sure the histogram is the current object.
2. In the property panel, set **Color Column** to `SEX`, **Color Aggr Type** to `avg`, **Invert Color Scheme** to `true` (non-default lead — the color configuration is applied first, before any split).
3. In the Context Panel, locate the three Color rows — **Color**, **Color Aggr Type**, and **Invert Color Scheme**. Verify all three appear active (fully opaque, not dimmed) — color-coding is available while no split is set (GROK-19761 baseline).
4. Note the current filtered row count as the baseline (all rows pass the filter).
5. Set the histogram's **Split Column** property to `RACE` (a categorical column — non-default, lead).
   Look at the same three Color rows again — each must now appear dimmed (disabled). The split disabling the Color rows is the product-state signal (GROK-19761). **Do NOT assert a Split Column value round-trip.**
6. While the split column is `RACE`, use the in-histogram range slider to narrow the displayed range to a sub-range of `AGE` (roughly 20 to 60).
   Verify the filtered row count stays valid — between zero and the total row count, and not above the baseline from step 4.
   Verify no new console error has been raised since step 5 (GROK-18399).
7. Reset the range to the full column extent (clear min/max or drag the slider to full width) so the filter is no longer narrowed.
8. Clear the **Split Column** (set to empty/none).
   Look at the three Color rows once more — all three are active (undimmed) again (round-trip: color-coding re-activates, GROK-19761 regression guard).

Expected:
- Step 3: all three Color rows appear active (fully opaque) — color available, no split.
- Step 5: all three Color rows appear dimmed (disabled) after the split is set (GROK-19761).
- Step 6: the filtered row count stays within the valid range and no new console errors appear (GROK-18399).
- Step 8: all three Color rows return to the active (undimmed) state (round-trip).

## Automation notes

Scenario 1:
- The three Color rows are `tr[name='prop-color']`, `tr[name='prop-color-aggr-type']`, `tr[name='prop-invert-color-scheme']` under `.property-grid`. Wait for `state: 'attached'` (NOT `state: 'visible'` — the rows may sit inside a collapsed accordion pane); do NOT try to expand accordion panes and do NOT gate on visibility. "Active" reads as `getComputedStyle(tr).opacity === '1'`; "dimmed" as opacity `< 1` (expected `'0.5'`).
- The histogram is made current via `grok.shell.o = h`; the baseline and filtered row counts are read from `df.filter.trueCount` (valid range: `[0, df.rowCount]`); the range narrow in step 6 is driven by setting the viewer's min/max range properties (e.g. `valueMin = 20`, `valueMax = 60`).
- Setup mechanics: close views via `grok.shell.closeAll()`; open demog from `System:DemoFiles/demog.csv`; add the viewer via `view.addViewer('Histogram')` and capture it as `h`.

**Out of scope (do NOT assert):**
- The per-category legend rendered on the canvas (canvas-only, no reliable DOM host).
- Multi-series visual modes: Split Stack, Normalize Values, distribution lines/splines (canvas-only, no product-state signal beyond the props themselves — asserting them would be a vacuous round-trip).
- splitColumnName set-then-read (vacuous round-trip).
