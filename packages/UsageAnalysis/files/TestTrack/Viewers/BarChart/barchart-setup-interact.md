---
feature: barchart
realizes_atlas:
  - barchart.cp.setup-and-interact
realizes:
  - viewers.bar-chart
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-17532
    status: fixed
  - id: github-3188
    status: fixed
  - id: github-2562
    status: fixed
realized_as:
  - barchart-setup-interact-spec.ts
expected_results:
  - anchor: "Scenario 1 Step 1"
    expectation: >-
      Per-category bars render for each distinct Split value, heights reflecting
      the aggregated Value.
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      With On Click set to Filter, clicking a bar filters the grid to that bar's
      category only.
  - anchor: "Scenario 1 Step 5"
    expectation: >-
      Clicking a blank canvas zone releases the active click-filter — the grid
      returns to the full row range (observed dev behavior: tc 64 → 100,
      distinct 1 → 5).
  - anchor: "Scenario 1 Step 9"
    expectation: >-
      Double-clicking the canvas invokes Reset View and leaves the chart intact
      and error-free.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      With On Click reverted to Select, clicking a bar selects that category's
      rows without filtering (filter count unchanged).
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      Pressing Esc clears the selection (no rows highlighted).
---

# Bar Chart — Setup and core interaction (On Click Filter vs Select)

## Setup

1. Close all open tables and viewers.
2. Open spgi-100 — a 100-row SPGI sample (`System:AppData/Chem/tests/spgi-100.csv`).
3. Add a Bar Chart viewer to the current table view.
4. In the Context Panel > Data section, set **Split** to **Primary Series Name**.
5. Set **Value** to **CAST Idea ID** and **Value Aggr Type** to **Count**.

> Actuation note: the bar clicks target the dominant **Triazoles** category by
> fixed canvas fractions (x ≈ 0.55·w, y ≈ 0.392·h). Only that dominant bar is
> hit-testable headless; other bar positions register no hit (recon
> 2026-07-21). Switching the filter to a second bar (md Scenario 1 Steps 6-7)
> and Alt+drag zoom (md Scenario 1 Step 8) are documented reductions —
> headless-unreachable — so Reset View (Step 9) is exercised directly without a
> preceding zoom.

## Scenarios

### Scenario 1: On Click = Filter — click bar filters table; double-click Reset View clears zoom

Steps:
1. Verify that per-category bars render for each distinct value in the **Primary Series Name**
   column, with bar heights reflecting the Count of CAST Idea ID per category.
2. In the Context Panel > Style section, set **On Click** to **Filter** (non-default).
3. Click a bar segment in the Bar Chart (e.g. the tallest bar).
4. Verify that the grid table is filtered to rows matching that bar's category only;
   all other rows are hidden from the grid (filter indicator appears in the grid header).
5. Click blank canvas space inside the Bar Chart. Verify the click-filter releases and the grid
   returns to the full row range (observed dev behavior 2026-07-21: clicking an empty canvas zone
   under On Click=Filter clears the active bar filter rather than leaving it intact).
6. Click a different bar segment to switch the active filter to that category.
7. Verify the grid now shows only rows for the newly clicked category.
8. Alt+drag horizontally across two bars to zoom into that range.
9. Double-click the canvas to invoke Reset View; verify zoom clears and the chart
   returns to the default full-range view showing all categories.

Expected:
- After Step 1: per-category bars are visible and non-empty for at least two distinct categories.
- After Step 4: the grid row count drops to only those rows belonging to the clicked category;
  no rows from other categories appear.
- After Step 7: the grid row count changes to match the second clicked category.
- After Step 9: all category bars are again visible at full scale; no zoom offset remains.

### Scenario 2: On Click = Select — click bar selects rows without filtering; Esc clears

Steps:
1. In the Context Panel > Style section, revert **On Click** to **Select** (default).
2. Click a bar segment in the Bar Chart (e.g. the tallest bar).
3. Verify that the rows belonging to that category are highlighted (selected) in the grid;
   the total row count in the grid does NOT change (no filter is applied).
4. Press **Esc**.
5. Verify that the selection is cleared; no rows remain highlighted in the grid.

Expected:
- After Step 3: selected rows in the grid are highlighted, but the total visible row count
  equals the full dataset row count (On Click = Select does not filter).
- After Step 5: no grid rows are highlighted; selection count shows 0.
