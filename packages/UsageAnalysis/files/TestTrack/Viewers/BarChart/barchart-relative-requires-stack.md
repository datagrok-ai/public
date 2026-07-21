---
feature: barchart
realizes_atlas:
  - barchart.int.relative-requires-stack
realizes:
  - viewers.bar-chart
priority: p1
target_layer: playwright
coverage_type: regression
realized_as:
  - barchart-relative-requires-stack-spec.ts
related_bugs: []
expected_results:
  - anchor: "Scenario 1 Step 3"
    expectation: >-
      With no Stack column, enabling Relative Values has no effect — bars keep
      absolute widths (not normalized).
  - anchor: "Scenario 2 Step 6"
    expectation: >-
      Setting a Stack column with Relative Values on normalizes bars to equal
      width, split into stacked segments.
  - anchor: "Scenario 2 Step 7"
    expectation: >-
      Removing the Stack column reverts bars to absolute widths — Relative Values
      inert again.
  - anchor: "Scenario 3 Step 10"
    expectation: >-
      Disabling Relative Values with no Stack column restores the baseline
      single-segment bars, no error.
---

# Bar Chart — Relative Values requires Stack column

## Setup

1. Close all open tables and viewers.
2. Open spgi-100 — a 100-row SPGI sample (`System:AppData/Chem/tests/spgi-100.csv`).
3. Add a Bar Chart viewer to the current table view.
4. In the Context Panel > Data section, set **Split** to **Primary Series names**.
5. Set **Value** to **CAST Idea ID** and **Value Aggr Type** to **Count**.
6. Confirm that **Stack** is set to **None** (no Stack column, baseline state).

## Scenarios

### Scenario 1: Relative Values is inert without a Stack column

Steps:
1. Verify the chart shows single-segment bars with absolute widths for each
   **Primary Series names** category — no stacking is active and bars vary
   in width according to their count.
2. In the Context Panel > Data section, locate **Relative Values** and enable it
   (toggle to on / checked — non-default).
3. Verify that bars still display with their original absolute widths — no normalization
   to equal width occurs. The chart must not become blank or throw an error.
   Relative Values is active but Stack is None, so the option resolves to a no-op.
4. In the Context Panel, confirm **Stack** is still **None**.

Expected:
- With no Stack column set, enabling Relative Values has no visible effect —
  bars keep their absolute widths and are NOT normalized to equal width.
  The Relative Values property is accepted but resolves to a no-op (relativeValues
  effectively false) because the stacking precondition is unmet.

### Scenario 2: Relative Values activates when Stack column is added

Steps:
5. With **Relative Values** still enabled, set **Stack** to **Stereo Category** in the Context Panel > Data section.
6. Verify the chart immediately switches to a normalized stacked view — each outer bar
   (per Primary Series names category) now occupies the same total width, with stacked
   segments proportional to the Stereo Category distribution within that category.
7. Remove the Stack column — set **Stack** back to **None**. Verify bars revert to
   single-segment absolute-width rendering (Relative Values is still on but inert again
   without a Stack column).
8. Re-add **Stack** to **Stereo Category**. Verify bars again normalize to equal width confirming
   the toggle is repeatable and not a one-time event.

Expected:
- After setting a Stack column with Relative Values still enabled, bars
  immediately normalize to equal width and split into stacked segments —
  each outer bar occupies the same total width representing 100% of
  its category distribution.
- After removing the Stack column while Relative Values remains enabled,
  bars revert to absolute widths — Relative Values is again inert with
  no Stack column present.

### Scenario 3: Round-trip — disable Relative Values, verify baseline restoration

Steps:
9. With **Stack** set to **None** and **Relative Values** on, disable **Relative Values**
   (toggle to off / unchecked — back to default).
10. Verify the chart shows single-segment absolute-width bars — the chart is in its
    baseline state with no stacking and no normalization. No error or blank chart.

Expected:
- After disabling Relative Values with no Stack column set, the chart
  state is the same as the baseline — single-segment bars at absolute
  widths. No error or blank chart occurs during any step of this
  toggle sequence.
