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
      returns to the full row range.
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
  - anchor: "Scenario 2 Step 6"
    expectation: >-
      Ctrl+clicking a bar selects that category's rows in full; Ctrl+clicking a
      different bar grows the selection to both categories (additive, not
      replacement), with the filter untouched.
  - anchor: "Scenario 2 Step 8"
    expectation: >-
      Shift+dragging a rubber-band rectangle over two bars selects the covered
      categories' rows in full, without filtering.
  - anchor: "Scenario 2 Step 9"
    expectation: >-
      Pressing Esc after the gesture-built selection clears it back to zero.
  - anchor: "Scenario 3 Step 1"
    expectation: >-
      Setting categorical color coding on the Split column through the grid
      recolors the bars (large canvas delta above the settle floor), and after a
      save → close → clear → reload layout round-trip the Split column's color
      scheme is restored and the bar chart reopens.
---

# Bar Chart — Setup and core interaction (On Click Filter vs Select)

## Setup

1. Close all open tables and viewers.
2. Open spgi-100 — a 100-row SPGI sample (`System:AppData/Chem/tests/spgi-100.csv`).
3. Add a Bar Chart viewer to the current table view.
4. In the Context Panel > Data section, set **Split** to **Primary Series Name**.
5. Set **Value** to **CAST Idea ID** and **Value Aggr Type** to **Count**.

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
   returns to the full row range (clicking an empty canvas zone under On Click=Filter clears the
   active bar filter rather than leaving it intact).
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
6. **Ctrl+click** a bar — its category's rows become selected. Then **Ctrl+click** a
   different bar — the selection grows to include the second category as well
   (additive selection, not replacement); the filter stays untouched.
7. Verify both clicked categories are fully selected (every row of each category is
   highlighted) and the grid still shows the full row range.
8. Clear the selection, then **Shift+drag** a rubber-band rectangle over two bars —
   the covered categories' rows become selected in full, without filtering.
9. Press **Esc** — the gesture-built selection clears back to zero.

Expected:
- After Step 3: selected rows in the grid are highlighted, but the total visible row count
  equals the full dataset row count (On Click = Select does not filter).
- After Step 5: no grid rows are highlighted; selection count shows 0.
- After Steps 6-7: two distinct categories are selected, the selection count equals the
  combined row count of those categories, and the filter count equals the full dataset.
- After Step 8: the rubber-band selects at least the two covered categories in full;
  no filter is applied.
- After Step 9: the selection count returns to 0.

### Scenario 3: grid color coding on the Split column drives bar colors and persists across a layout round-trip

Steps:
1. In the grid, apply categorical **Color coding** to the **Primary Series Name** (Split)
   column and edit the scheme. Verify the bars in the Bar Chart recolor to match the column's
   category colors. Save the layout, close the bar chart, clear the column color coding, then
   apply the saved layout again and verify the color scheme is restored and the bar chart reopens.

Expected:
- After the color scheme is applied: the bars adopt the Split column's category colors
  (a large canvas repaint versus the uncolored baseline).
- After the layout round-trip: the Split column's categorical color scheme is restored
  (the color-coding tag matches the applied scheme) and the bar chart reopens intact.

## Automation notes

Scenarios 1-2: the bar clicks target the dominant **Triazoles** category by
fixed canvas fractions (x ≈ 0.55·w, y ≈ 0.392·h). Only that dominant bar is
hit-testable headless; other bar positions register no hit. Switching the
filter to a second bar (md Scenario 1 Steps 6-7)
and Alt+drag zoom (md Scenario 1 Step 8) are documented reductions —
headless-unreachable — so Reset View (Step 9) is exercised directly without a
preceding zoom. The additive-gesture steps (Scenario 2 Steps 6-9) do not use
fixed fractions: they locate individual bars by scanning the canvas for the
default bar-fill color (contiguous bands of #96d794 pixels), so two distinct
bars are clickable there.
