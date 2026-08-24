---
feature: pivottable
realizes_atlas:
  - pivottable.cp.rowsource-filter-selection-links
  - pivottable.int.filtering-requires-row-source-all
  - pivottable.int.selection-link-gated-by-row-source
realizes:
  - viewers.pivot-table
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
  - id: GROK-17726
    status: fixed
  - id: GROK-15518
    status: fixed
realized_as:
  - pivottable-rowsource-filter-selection-spec.ts

expected_results:
  - anchor: "Scenario 1 Step 4"
    expectation: >-
      With a SEX = "M" source filter active, the pivot's per-DIS_POP aggregate
      values differ from the unfiltered baseline and match an independent
      groupBy avg(AGE) computed over the filtered rows.
  - anchor: "Scenario 2 Step 3"
    expectation: >-
      Row Source is All and filteringEnabled is true, and the aggregated row set
      does not widen (the aggregation still runs over the SEX = "M" subset).
  - anchor: "Scenario 2 Step 5"
    expectation: >-
      A real cell click narrows df.filter.trueCount from the pre-click value to
      a smaller positive integer.
  - anchor: "Scenario 2 Step 7"
    expectation: >-
      With a SEX = "M" panel filter active, a real cell click REPLACES df.filter
      with the clicked DIS_POP group's full source bitset (single surviving key,
      count equals that group's full source-row total, non-M rows survive) and
      the count never grows past the pre-click baseline (GROK-17726
      monotonicity).
  - anchor: "Scenario 3 Step 5"
    expectation: >-
      With filteringEnabled=false under All, a cell click leaves
      df.filter.trueCount unchanged (I1 negative).
  - anchor: "Scenario 3 Step 9"
    expectation: >-
      At rowSource=Filtered with filteringEnabled=true, a cell click leaves
      df.filter.trueCount unchanged (I1 negative).
  - anchor: "Scenario 4 Step 4"
    expectation: >-
      At rowSource=Filtered, a real Ctrl+click on an aggregate row of the
      configured pivot (Group by DIS_POP, Aggregate avg(AGE), asserted via
      durable caption text) selects that group's full source rows via the
      selection-to-selection link (df.selection.trueCount equals the clicked
      group's source-row count, exactly one key) with df.filter untouched and no
      warning.
  - anchor: "Scenario 4 Step 9"
    expectation: >-
      At rowSource=All, a real Ctrl+click on a different aggregate row
      re-selects the newly clicked group's full source rows
      (df.selection.trueCount equals that group's source-row count, exactly one
      key) with df.filter untouched and no warning, while the pivot config
      (Group by DIS_POP, Aggregate avg(AGE)) survives the rebuild (asserted via
      durable caption text).
gate_verdicts:
  b:
    verdict: PASS
    cycle_id: direct-gate-b-2026-08-23-pivottable-rowsource-filter-selection-proved
    timestamp: 2026-08-23T21:00:00Z
    spec_runs:
      - spec: pivottable-rowsource-filter-selection-spec.ts
        result: passed
        attempts: 3
        duration_seconds: 18
        run_mode: headless-cold
        failure_keys: []
---

# Pivot Table — Row Source, filter and selection links

## Setup

1. Log in to Datagrok.
2. Close all open views.
3. Open the demog dataset via the File Browser (Files > App Files > Demo Files > demog.csv)
   and wait for the table view to appear.
4. Add the Pivot Table viewer from the Toolbox (Viewers section) using the viewer icon.
5. Wait for the tag-editor header to appear with three rows labeled Group by, Aggregate,
   and Pivot.
6. Configure a minimal cross-tab: add DIS_POP as the Group by column using the + icon next
   to the Group by row; confirm the Aggregate row shows avg(AGE) from auto-setup.

## Scenarios

### Scenario 1: Filter consumption — pivot recomputes when source filter changes

Steps:
1. Note the aggregate values visible in the pivot's inner grid for the top two
   DIS_POP groups (read the values directly from the pivot grid cells).
2. Open the Filter Panel via the toolbar icon (or View > Filter Panel) and wait for it to appear.
3. In the Filter Panel, locate the SEX column filter and select only "M"; wait for the
   filter to apply (the table row count in the status bar decreases).
4. Read the aggregate values for the same DIS_POP groups from the pivot's inner grid again,
   and verify they differ from the unfiltered baseline recorded in step 1.
5. Note the current filtered row count shown in the status bar.
6. Close the Filter Panel and clear the SEX filter (click the reset / clear-all button in the
   Filter Panel) before proceeding.

### Scenario 2: Cell-driven filtering under rowSource=All + filteringEnabled=true

Steps:
1. Ensure the pivot shows at least two DIS_POP rows in the inner grid (use the minimal setup
   from the Setup section if needed; close and re-add the viewer if the state is uncertain).
2. Open the Filter Panel and apply the SEX = "M" filter so the source dataframe is pre-filtered.
3. In the Pivot Table viewer's Settings panel (gear icon on the viewer header), set the
   Row Source property to All and confirm filteringEnabled is true (the default).
4. Note the current filtered row count shown in the status bar (call it T_before).
5. Click one aggregate cell in the pivot's inner grid (pick a non-empty DIS_POP cell).
6. Note the updated filtered row count in the status bar (call it T_after) and verify the
   active filters shown in the Filter Panel reflect the clicked cell's DIS_POP key.
7. Verify T_after equals the clicked DIS_POP group's full source-row count — the click
   replaces the SEX filter, it does not intersect it.
8. Clear the source filter (reset via the Filter Panel) before the next scenario.

### Scenario 3: Negative path — filter does NOT fire when conditions are not met

Steps:
1. Ensure the pivot is present with DIS_POP in Group by; Row Source should be All from
   Scenario 2 (re-open or re-configure if the state is uncertain).
2. In the Settings panel, set filteringEnabled to false.
3. Note the current filtered row count shown in the status bar.
4. Click one aggregate cell in the pivot's inner grid.
5. Verify the filtered row count in the status bar has not changed and the Filter Panel
   shows no new active filters.
6. Now set filteringEnabled back to true, and then set Row Source back to Filtered (the default
   value) in the Settings panel.
7. Note the current filtered row count shown in the status bar.
8. Click one aggregate cell in the pivot's inner grid.
9. Verify the filtered row count in the status bar has not changed and the Filter Panel
   shows no new active filters.
10. Restore filteringEnabled=true and Row Source=Filtered before the next scenario.

### Scenario 4: Selection link — aggregate row click moves selection

Steps:
1. Ensure the pivot is present with DIS_POP in Group by and Row Source = Filtered (the default).
2. Note the current selection count shown in the status bar or grid selection indicator.
3. Ctrl+click the DIS_POP key cell of one row in the pivot's inner grid to select that
   aggregate row (a plain click does not select — it only moves the grid's current cell).
4. Verify the selection count in the status bar has increased to reflect the source rows
   belonging to the selected DIS_POP group.
5. Ctrl+click the key cell of a different row in the pivot's inner grid.
6. Verify the selection count updates to reflect the newly selected group.
7. In the Settings panel, set Row Source to All.
8. Ctrl+click the key cell of another aggregate row in the pivot's inner grid.
9. Verify the selection count updates to reflect the source rows of the newly selected group.
10. Close all open views to clean up.

## Automation notes

- df.filter.trueCount read: `await page.evaluate(() => grok.shell.t.filter.trueCount)` where
  `grok.shell.t` is the demog table.
- df.selection.trueCount read: `await page.evaluate(() => grok.shell.t.selection.trueCount)`.
- dataFrame.rows.filters read: confirmed live — the collection iterates empty (the filter label
  is not JS-readable), so assertions use df.filter.trueCount instead.
- Independent group count (Scenario 2 Step 7): the clicked DIS_POP group's full source-row
  count, computed by iterating the DIS_POP column over the source rows.
- Row Source and filteringEnabled are viewer properties (rowSource, filteringEnabled),
  confirmed live.
- The cell click in Scenarios 2-3 must be a real browser click (Playwright `page.mouse.click`)
  — a simulated MouseEvent from `page.evaluate` does not fire the Dart listener; confirmed live
  that the trusted click fires it reliably.
- target_layer rationale: all four channels are driven by real UI events — a cell click is a
  Dart-side MouseEvent that the browser must dispatch; the filter and selection signals are read
  via the JS API but the trigger is always a UI action. The Settings panel for Row Source and
  filteringEnabled is a Playwright-accessible DOM panel.
- Confirmed live: the Row Source values read "Filtered" and "All"; filteringEnabled is a
  boolean viewer property; the pivot row headers show the DIS_POP values as-is.
- I1 closure (pivottable.int.filtering-requires-row-source-all): Scenario 3, Steps 2-9 — both
  negative conditions tested (filteringEnabled=false under All; Filtered+filteringEnabled=true).
- I3 closure (pivottable.int.selection-link-gated-by-row-source): Scenario 4, Steps 3-9 — the
  attached branch exercised at both Filtered and All (the detached Selected-family branch is out
  of scope per atlas decision).
- Won't-fix guard (GROK-14996): no in-viewer undo is asserted; the scenario clears the source
  filter via the Filter Panel before continuing.
- Row Source scope (binding): the aggregation always runs over the source rows for the active
  Row Source setting; Row Source does not change what rows are aggregated independently of the filter.
- No persistence tail: channel-isolation scenarios do not duplicate the layout + project
  round-trip tail present in the configure-crosstab scenario.
