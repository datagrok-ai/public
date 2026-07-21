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
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      Color property-grid rows (tr[name='prop-color'],
      tr[name='prop-color-aggr-type'], tr[name='prop-invert-color-scheme']) are
      active — getComputedStyle(tr).opacity === '1' for all three.
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      After setting a Split Column, all three Color property-grid rows become
      disabled — getComputedStyle(tr).opacity === '0.5' (or < 1.0) for all
      three.
  - anchor: "Scenario 1 Step 6"
    expectation: >-
      While the split column is set and the range slider is adjusted,
      df.filter.trueCount stays in the valid range [0, df.rowCount] and no new
      console error is raised.
  - anchor: "Scenario 1 Step 8"
    expectation: >-
      After clearing the Split Column, the three Color property-grid rows
      re-activate — getComputedStyle(tr).opacity returns to '1' for all three
      (round-trip).
realized_as:
  - histogram-split-and-color-spec.ts
---

# Histogram — Split vs Color-coding Transition

## Setup

1. Close all open tables and viewers (`grok.shell.closeAll()`).
2. Open the **demog** dataset (use `readDataframe('demog')` or navigate Files > demog.csv, open the table).
3. Add a Histogram viewer to the table view (`view.addViewer('Histogram')`); capture it as `h`.
4. Make the histogram the current object so the Context Panel shows its properties: `grok.shell.o = h`.
5. In the property panel (Context Panel), set the **Value Column** to `AGE`.

## Scenarios

### Scenario 1: Split column disables color-coding, range filter survives, color re-activates on split clear

**Purpose:** Assert the product-state transition when a split column is set and cleared — the three Color property-grid rows must disable (opacity 0.5) when a split column is active and re-enable (opacity 1.0) when it is cleared. The range slider while split is set must not error and must keep df.filter valid (GROK-18399). Color UI disabling is the observable product-state signal for GROK-19761; a splitColumnName set-then-read is NOT asserted (vacuous round-trip).

Steps:

1. With no split column set, open the Context Panel if not already shown. Ensure the histogram is the current object (`grok.shell.o = h`).
2. In the property panel, set **Color Column** to `SEX`, **Color Aggr Type** to `avg`, **Invert Color Scheme** to `true` (non-default lead — the color configuration is applied first, before any split).
3. Query the three Color property-grid rows via their `name` attributes under `.property-grid`:
   `tr[name='prop-color']`, `tr[name='prop-color-aggr-type']`, `tr[name='prop-invert-color-scheme']`.
   Wait for `state: 'attached'` (NOT `state: 'visible'` — rows may be inside a collapsed accordion pane).
   Read `getComputedStyle(tr).opacity` for each — it must equal `'1'` for all three (color-coding active, GROK-19761 baseline). **Do NOT try to expand accordion panes; do NOT gate on visibility.**
4. Record `const fullCount = df.filter.trueCount` (baseline: all rows in).
5. Set the histogram's **Split Column** property to `RACE` (a categorical column — non-default, lead).
   Re-query the same three `tr` nodes (they remain attached) and read their computed opacity.
   Each must now be `< 1` (expected: `'0.5'`) — the split disabling the Color property rows is the product-state signal (GROK-19761). **Do NOT assert splitColumnName round-trip.**
6. While the split column is `RACE`, use the in-histogram range-slider (or set `h.min` / `h.max` via JS) to narrow the displayed range to a sub-range of `AGE` (e.g. set min to 20, max to 60).
   Assert: `df.filter.trueCount` is in `[0, df.rowCount]` (valid, non-negative) and `df.filter.trueCount <= fullCount`.
   Assert: no new console error has been raised since step 5 (GROK-18399).
7. Reset the range to the full column extent (clear min/max or drag slider to full width) so df.filter is no longer narrowed.
8. Clear the **Split Column** (set to empty/none).
   Re-query the three Color `tr` nodes and read computed opacity — each must be `'1'` again (round-trip: color-coding re-activates, GROK-19761 regression guard).

Expected:
- Step 3: all three Color property-grid rows have computed opacity `'1'` (color active, no split).
- Step 5: all three Color property-grid rows have computed opacity `< 1` (opacity `'0.5'`) after split is set (GROK-19761).
- Step 6: `df.filter.trueCount` ∈ [0, df.rowCount]; no new console errors (GROK-18399).
- Step 8: all three Color property-grid rows return to computed opacity `'1'` (round-trip).

**Out of scope (do NOT assert):**
- The per-category legend rendered on the canvas (canvas-only, no reliable DOM host).
- Multi-series visual modes: Split Stack, Normalize Values, distribution lines/splines (canvas-only, no product-state signal beyond the props themselves — asserting them would be a vacuous round-trip).
- splitColumnName set-then-read (vacuous round-trip).
