---
feature: grid
realizes_atlas:
- grid.cp.rows-select-filter-navigate
- grid.int.frozen-rows-selection-mapping
- grid.int.sort-filter-shared-order
realizes:
- viewers.grid
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
- id: GROK-17455
  status: fixed
- id: GROK-18104
  status: fixed
realized_as:
- grid-rows-select-filter-navigate-spec.ts
expected_results:
- anchor: Step 3
  expectation: After clicking row 5's data cell, df.currentRowIdx equals the table index of visual row
    5 (converted via grid.gridRowToTable); df.selection.trueCount is 0 (a single click sets current
    row, not a selection).
- anchor: Step 5
  expectation: Ctrl+clicking a row-number cell adds that row to the selection and makes it current;
    the selection count grows by one.
- anchor: Step 6
  expectation: After pressing Ctrl+A, df.selection.trueCount equals df.rowCount (all rows selected).
- anchor: Step 7
  expectation: 'After pressing Esc, all three channels reset simultaneously: df.selection.trueCount
    is 0, df.columns.selected contains no selected columns, and df.currentRowIdx is -1.'
- anchor: Step 10
  expectation: After rapid Ctrl+clicks across five or more column headers, the console-error delta is
    0 — no console error (GROK-17455 guard).
- anchor: Step 13
  expectation: Arrow key presses move df.currentRowIdx by 1 each time in the expected direction; df.currentCol
    changes on left/right arrow.
- anchor: Step 15
  expectation: Ctrl+Home sets df.currentRowIdx to 0; Ctrl+End sets df.currentRowIdx to df.rowCount -
    1.
- anchor: Step 16
  expectation: The first PageDown press advances df.currentRowIdx by one visible page; the second PageDown
    press advances it by another full page (GROK-18104 guard — paging must not stall after the first
    invocation).
- anchor: Step 17
  expectation: 'Space toggles the current row''s membership in df.selection: pressing Space on a currently-unselected
    row adds it (df.selection.trueCount increases by 1); pressing Space again removes it (df.selection.trueCount
    decreases by 1).'
- anchor: Step 18
  expectation: Shift+Enter selects all rows sharing the current cell's value; df.selection.trueCount
    equals the count of rows with that value.
- anchor: Step 20
  expectation: With props.allowRowSelection set to false, pressing Shift+Down Arrow moves df.currentRowIdx
    without changing df.selection.trueCount (navigation proceeds but no rows are selected).
- anchor: Step 22
  expectation: With 2 rows pinned and the grid scrolled down, Ctrl+clicking three data-area rows leaves
    the selection holding exactly those table rows — the pinned rows shift the visual-to-table mapping,
    and this is the check that the correction is applied.
- anchor: Step 24
  expectation: With both df.filter active and a sort applied, navigating rows via arrow keys and reading
    grid.gridRowToTable(i) for each visual position yields the correct table indices (the combinedFilter
    / _order mapping is respected; calling getRowOrder() is avoided).
- anchor: Step 25
  expectation: A keyboard range-select (Shift+Down repeated) under the active sort+filter produces a
    df.selection whose selected indices are exactly the intended table indices (not raw visual offsets).
- anchor: Step 26
  expectation: Switching rowSource from All to Filtered shows only rows where df.filter is true (the
    visible-row count matches df.filter.trueCount); switching to Selected shows only selected rows;
    switching back to All restores the full row count — while df.filter.trueCount remains unchanged
    throughout (by-design divergence, not a bug).
- anchor: Step 27
  expectation: After typing a search term, cells matching it (substring, case-insensitive) show a tint
    colour readable via grid.cell(col, row).color that differs from a plain cell. Clearing the search
    text does NOT restore df.filter on its own, so the filter is restored explicitly before the next
    step. That the Toolbox Search pane opens and takes focus is platform behaviour, not the grid's —
    it is only a precondition here.
---

# Grid — Row Selection, Filter, and Keyboard Navigation

## Setup

1. Close all open views.
2. Open the demog dataset (demog.csv) to produce a Table View with the main grid.
3. Record the console-error count before any interactions begin (the baseline for
   console-error delta checks).

## Scenarios

### Scenario 1: Row selection — mouse clicks and Esc compound reset

Steps:
1. Focus the grid by clicking inside the data area.
2. Ctrl+click the row-number cell of visual row 5 — it becomes the current row and joins the
   selection.
3. Verify that row 5 becomes the current row and no rows are selected.
4. Hold Ctrl and click the row-header of visual row 15.
5. Verify that row 15 is added to the selection and becomes the current row.
6. Press Ctrl+A to select all rows.
7. Verify that all rows are selected.
8. Press Esc.
9. Verify the compound reset: no rows are selected, no columns are selected, and the
   current row is cleared.

Expected:
- Step 3: the current row is the table row behind visual row 5; the number of selected rows
  is 0 (single click sets current, not selection).
- Step 5: the Ctrl+clicked row joins the selection and becomes the current row.
- Step 9: every row in the table is selected.
- Step 11: 0 rows are selected; no columns are selected; the current row
  is -1.

### Scenario 2: Rapid Ctrl+click on column headers raises no error (GROK-17455)

Steps:
1. Close all and reopen demog to start from a clean state.
2. Rapidly Ctrl+click five or more column headers in quick succession.
3. Verify that no error appeared in the console.

Expected:
- Step 3: the console-error count is unchanged from before the clicks — no
  any console error (GROK-17455 guard). Whether the columns
  end up visibly selected is checked by hand (see the manual checklist).

### Scenario 3: Keyboard navigation — arrows, Ctrl+Home/End, PageDown twice (GROK-18104)

Steps:
1. Close all and reopen demog to start from a clean state.
2. Click any data cell to focus the grid and set a known current cell.
3. Press Down Arrow twice and then Right Arrow twice.
4. Verify that the current row and column have moved accordingly.
5. Press Ctrl+Home.
6. Verify that the current position is the first row.
7. Press Ctrl+End.
8. Verify that the current position is the last row.
9. Press PageDown.
10. Note the resulting current row position after the first PageDown.
11. Press PageDown again.
12. Verify that the second PageDown advances the position by another full page (not stuck
    at the same position as after the first press).

Expected:
- Step 4: the current row and the current column have each advanced by two from their starting
  values.
- Step 6: the current row is the first row of the table.
- Step 8: the current row is the last row of the table.
- Step 16: the current row after the second Page Down is further down than it was
  after the first PageDown — paging advances both times (GROK-18104 guard).

### Scenario 4: Space toggle and Shift+Enter value-group selection

Steps:
1. Close all and reopen demog to start from a clean state.
2. Ctrl+click the row-number cell of row 3 to make it current (the selection should
   be 0 after this single click).
3. Press Space.
4. Verify that row 3 is now selected (selection count is 1).
5. Press Space again.
6. Verify that row 3 is no longer selected (selection count is 0).
7. Navigate to a cell in the SEX column (which has a small number of distinct values).
8. Press Shift+Enter.
9. Verify that all rows sharing the same SEX value as the current cell are selected.

Expected:
- Step 4: 1 rows are selected.
- Step 6: 0 rows are selected.
- Step 9: the selected rows are every row whose SEX value
  matches the current cell's value.

### Scenario 5: Negative path — allowRowSelection false

Steps:
1. Close all and reopen demog to start from a clean state.
2. Open grid properties (the gear icon in the grid's toolbar) and set the "Allow Row
   Selection" toggle to off. Close the properties panel.
3. Navigate to any data row using the Down Arrow key.
4. Hold Shift and press Down Arrow three times.
5. Verify that navigation proceeds but no rows are selected.
6. Open grid properties again and restore "Allow Row Selection" to on.
7. Press Esc to clear state, then put the current row on row 2 with the Down Arrow and press
   Space to toggle it into the selection.
8. Verify that selection is working again.

Expected:
- Step 5: the current row has moved three rows down; 0 rows are selected (Shift+arrow
  navigates without selecting when allowRowSelection is false).
- Step 8: rows are selected again once row selection is switched back on.

### Scenario 6: Pinned-row interplay — selection coordinate mapping (typed interaction grid.frozen-rows-selection-mapping)

Steps:
1. Close all and reopen demog to start from a clean state.
2. Right-click the row-header of row 1 and choose "Pin Row" from the context menu; repeat
   for row 2 so that 2 rows are pinned at the top of the grid.
3. Scroll the grid down so that pinned rows remain visible at the top but the data area
   shows later rows.
4. Ctrl+click the row-number cells of 3 consecutive data-area rows (not the pinned ones) so each
   joins the selection.
5. Verify that the selected row indices match exactly the intended table rows.

Expected:
- Step 5: the selected rows are exactly the 3
  rows that were dragged over — not indices shifted by the 2-pin offset.

### Scenario 7: Sort + filter shared order mapping (typed interaction grid.sort-filter-shared-order)

Steps:
1. Close all and reopen demog to start from a clean state.
2. Click the filter icon in the grid's toolbar to open the Filter pane; apply a filter to
   the AGE column to keep only rows with AGE greater than 30. Close the Filter pane.
3. Double-click the AGE column header twice (first double-click gives descending, second
   gives ascending) so that the visible rows are sorted ascending by AGE.
4. Using Down Arrow, navigate through 5 visual rows and read the current-row table index
   at each step by observing the highlighted row in the grid.
5. Press Shift+Down Arrow 4 times to range-select 4 rows.
6. Verify that the selected indices are the correct table rows.
7. Verify that the filtered row count stays the same throughout these navigation steps.

Expected:
- Step 6: 4 rows are selected; the selected rows are
  that correspond to the rows visible under the sort+filter combination — confirmed by
  reading grid.gridRowToTable(i) for each visual position (getRowOrder() is never called
  on the grid to avoid materializing _order as a side effect).
- Step 7: the filtered row count is the same as right after the filter was applied
  in Step 2 — sorting and keyboard navigation do not modify the filter's true count.

### Scenario 8: Row source switching — Filtered / Selected / All divergence

Steps:
1. Close all and reopen demog to start from a clean state.
2. Apply a filter to the SEX column (keep only one value, e.g. M) via the Filter pane.
   Note how many rows the filter keeps.
3. Select about 5 rows via Ctrl+click on their row-headers.
4. In the grid's toolbar or Context Panel, switch the Row Source to "Filtered".
5. Verify that the visible row count matches the number of rows the filter keeps.
6. Switch the Row Source to "Selected".
7. Verify that the visible row count matches the number of selected rows.
8. Switch the Row Source back to "All".
9. Verify that all rows are visible again.
10. Verify that the filter itself is untouched throughout (switching the row source must not change what the filter keeps).

Expected:
- Step 5: the grid shows exactly the filtered rows, and the filtered row count itself does not change.
- Step 7: The grid displays the number of selected rows (5 from Step 3); the filtered row count does not change.
- Step 9: the grid shows every row of the table.
- Step 10: the filtered row count is identical at steps 5, 7 and 9 — what the grid shows and what the filter holds are deliberately allowed to differ.

### Scenario 9: Ctrl+F search — tint signal and explicit filter restore

Steps:
1. Close all and reopen demog to start from a clean state.
2. Press Ctrl+F to open the search pane.
3. The Toolbox Search pane opens with a text input — that is platform behavior and is only
   needed here so there is somewhere to type.
4. Type a value that appears in the AGE column (e.g. "35").
5. Verify that cells matching "35" (substring, case-insensitive) show a distinct tint color.
6. Select a non-matching cell. Verify its color is the plain background (no tint).
7. Clear the search input field.
8. Explicitly clear the filter (clearing the search text alone
   does not do this).
9. Verify that the grid shows all rows normally.

Expected:
- Step 3: The Toolbox Search pane is open with focus in its input.
- Step 5: grid.cell(col, row).color for a matching cell differs from the plain background
  color (a tint is applied).
- Step 6: grid.cell(col, row).color for a non-matching cell equals the plain background
  color.
- Step 9: All rows are visible in the grid after the explicit filter restore.

## Automation notes
- ANCHOR CONVENTION: the "Step N" labels on the expectations in this file's machine zone are the
  SPEC's step numbers — the automated check resolves them against the spec's step titles, which is
  why they run past the highest number you see below. The body numbers its steps per scenario, for
  a human executing one scenario at a time. The two numberings are independent on purpose; each
  expectation states its own observable, so a reader never has to follow the number to know what is
  being claimed. Do not renumber the body to chase them.

- the Ctrl+click-onto-an-existing-range case is manual and lives in the section's single manual checklist, grid-ui.md — not in a
  companion file beside this scenario. One list per section: a tester opening one file must not
  miss manual steps kept somewhere else.

- RULED OUT for mouse-driven row selection: a trusted DRAG on the row-number strip, and a
  Shift+click range built by mouse. Both were driven live and neither selects. Ctrl+click on the
  strip is NOT in that list — it works and is automated here. The manual half lives in grid-ui.md;
  there is no per-scenario companion file.

- The Space toggle and Shift+arrow paths are keyboard, not mouse — if the mouse gesture stays
  unreachable, they still prove the selection channel, and the mouse half becomes an honest
  manual line in grid-ui.md rather than a JS-API substitution.


### target_layer rationale
`playwright` — this scenario is a multi-step UI flow requiring real mouse clicks on canvas
row-headers, keyboard events (Ctrl+A, Esc, Tab corners, PageDown), drag gestures over data
cells, and toolbar interactions for Filter, Row Source, and Search. All assertions are
product-state reads (df-level: `df.selection.trueCount`, `df.currentRowIdx`,
`df.columns.selected`, `df.filter.trueCount`; grid-local: `grid.gridRowToTable(i)`,
`grid.cell().color`) accessible without pixel inspection — but the actuation is
necessarily browser-driven.

### Signal discipline
All expected-results assertions name product-state signals:
- Row selection state: `df.selection.trueCount`, `df.selection.getSelectedIndexes()`
- Current row / column: `df.currentRowIdx`, `df.currentCol`
- Column selection: `df.columns.selected` (names), `Column.isSelected`
- Filter truth count: `df.filter.trueCount` (unchanged by sort/navigation — BY DESIGN)
- Row-order mapping: `grid.gridRowToTable(i)` (never `grid.getRowOrder()` on an
  unsorted+unfiltered grid — it materializes `_order` as a side effect)
- Search tint: `grid.cell(col, row).color` vs. plain background
- Console errors: delta from baseline (0 for GROK-17455, 0 for all pagination steps)

### Typed interactions realized
- `grid.frozen-rows-selection-mapping` — Scenario 7 (pinned-row coordinate correction)
- `grid.sort-filter-shared-order` — Scenario 8 (sort + filter shared _order mapping)

### Known signal gaps
- Stats rows are canvas-only with no product-readable signal (see atlas cp description) —
  not exercised in this cp.
- The search tint is verified via `grid.cell().color` (a product read); no pixel-count
  assertion is used (near-vacuous for a grid and rejected by the assertion-strength
  doctrine).

### Bug guards embedded
- GROK-17455: Scenario 2, Step 5 — rapid Ctrl+click column selection must produce
  console-error delta 0.
- GROK-18104: Scenario 4, Step 16 — second PageDown must advance position, not stall.

---
{
  "order": 12,
  "datasets": ["System:DemoFiles/demog.csv"]
}
