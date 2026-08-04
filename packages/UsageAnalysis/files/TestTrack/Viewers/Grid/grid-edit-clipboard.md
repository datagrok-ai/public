---
feature: grid
realizes_atlas:
  - grid.cp.edit-clipboard
realizes:
  - viewers.grid
priority: p1
target_layer: playwright
coverage_type: regression
related_bugs:
  - id: GROK-20010
    status: fixed
  - id: GROK-20266
    status: fixed
realized_as:
  - grid-edit-clipboard-spec.ts
expected_results:
  - anchor: Step 4
    expectation: >-
      After pressing Enter to commit the edit, df.col('AGE').get(tableIdx)
      equals the newly typed value, the onCellValueEdited event count increased
      by 1, and the DOM ValueEditor element is no longer present in the grid
      root.
  - anchor: Step 6
    expectation: >-
      After pressing Esc inside an open cell editor, the AGE cell value is
      unchanged and the DOM ValueEditor element is no longer present in the grid
      root.
  - anchor: Step 8
    expectation: >-
      After pressing Delete on the current cell with no modifier key held,
      df.col('AGE').get(tableIdx) is null and the DOM ValueEditor is not
      present.
  - anchor: Step 11
    expectation: >-
      With props.allowEdit false, double-clicking an AGE cell leaves the DOM
      ValueEditor absent and a read-only warning balloon element appears in the
      DOM.
  - anchor: Step 13
    expectation: >-
      After restoring props.allowEdit to true and typing a digit on the current
      cell, df.col('AGE').get(tableIdx) equals the typed digit.
  - anchor: Step 16
    expectation: >-
      After pressing Ctrl+Shift+C on the current cell, the clipboard text equals
      the single cell value with no tabs or newlines.
  - anchor: Step 18
    expectation: >-
      After selecting 5 rows and pressing Ctrl+C, the multi-row copy of exactly
      5 selected rows is driven without raising a console or page error.
  - anchor: Step 21
    expectation: >-
      After Ctrl+A then Ctrl+C then Ctrl+V, the paste applies with a
      console-error delta of 0 and df.rowCount unchanged (GROK-20010 guard).
  - anchor: Step 23
    expectation: >-
      After pressing Ctrl+V into a different cell than the one copied,
      df.col('AGE').get(targetTableIdx) equals the previously copied cell value.
  - anchor: Step 27
    expectation: >-
      After pressing Shift+Del on the selected rows, df.rowCount has decreased
      by the number of deleted rows.
  - anchor: Step 28
    expectation: >-
      After pressing Ctrl+Z to undo the Shift+Del, df.rowCount returns to its
      value before the deletion.
---

# Grid — Cell Editing and Clipboard

## Setup

1. Close all open views.
2. Open the demog dataset (demog.csv) to produce a Table View with the main grid.
3. Record the starting console-error count as the baseline for all error-delta checks.
4. Record the table's row count as the baseline row count for undo/delete checks.

## Scenarios

### Scenario 1: Basic cell edit — commit with Enter and cancel with Esc

Steps:
1. Click any data cell in the AGE column to set it as the current cell. Note
   the current row's AGE value and its table-row index.
2. Double-click that same AGE cell to open the inline editor.
3. Verify that a DOM ValueEditor element is present inside the grid root while
   the editor is open.
4. Clear the editor input and type a new numeric value (e.g. 99), then press Enter.
5. Verify that the cell now shows the new value and the editor is dismissed.
6. Double-click the same cell again to re-open the editor.
7. Press Esc without typing.
8. Verify that the cell value is unchanged (the Esc cancelled the edit without
   modifying the data).

Expected:
- Step 3: A DOM ValueEditor input element is visible inside the grid root.
- Step 5: After Enter, the AGE cell's value equals 99; the the cell-edited event
  event count increased by 1; the inline editor is gone.
- Step 8: the AGE cell's value is still 99 (the Esc left the value intact).

### Scenario 2: Delete key clears current cell (no modifier)

Steps:
1. Close all and reopen demog to start from a clean state. Record the baseline.
2. Click an AGE cell to make it the current cell; note its table-row index and value.
3. Press Delete (with no Shift or Ctrl held).
4. Verify that the cell is cleared and no editor has opened.

Expected:
- Step 4: the AGE cell is empty — the cell was cleared. The
  DOM ValueEditor is not present (Delete did not open an editor).

### Scenario 3: Negative path — read-only grid rejects edits

Steps:
1. Close all and reopen demog to start from a clean state.
2. Open the grid properties via the gear icon in the grid toolbar. Set "Allow Edit"
   to off. Close the properties panel.
3. Double-click an AGE cell.
4. Verify that no editor opens and a warning balloon appears (because
   read-only notifications are enabled is on by default).
5. Open grid properties again and set "Allow Edit" back to on. Close the properties
   panel.
6. Click an AGE cell to make it current.
7. Type a single digit (e.g. 7) directly on the focused cell.
8. Press Enter.
9. Verify that the cell accepted the typed digit.

Expected:
- Step 4: The DOM ValueEditor is absent after the double-click; a warning balloon
  element is present in the DOM (read-only rejection is signalled visually).
- Step 9: the AGE cell's value equals 7 — the keystroke-replaces-value path
  works once editing is re-enabled.

### Scenario 4: Clipboard — single-cell Ctrl+Shift+C and multi-row Ctrl+C

Steps:
1. Close all and reopen demog to start from a clean state.
2. Click an AGE cell to make it the current cell. Note its value.
3. Press Ctrl+Shift+C (copy current cell value).
4. Verify that the clipboard holds the single cell value.
5. Ctrl+click the row-headers of 5 distinct rows to select them.
6. Press Ctrl+C (copy selection as TSV block).
7. Verify that the clipboard holds a TSV block with exactly 5 data rows.

Expected:
- Step 4: The clipboard text equals the single AGE value with no tab or newline
  characters.
- Step 7: The clipboard is a TSV block with 5 rows; each row is tab-separated and
  the columns match the grid's column order for those rows.

### Scenario 5: Ctrl+A → Ctrl+C → Ctrl+V does not error (GROK-20010)

Steps:
1. Close all and reopen demog to start from a clean state. Record the baseline
   console-error count and the table's row count.
2. Press Ctrl+A to select all rows.
3. Press Ctrl+C to copy.
4. Click a cell to set the paste target.
5. Press Ctrl+V to paste.
6. Verify that no console errors occurred and the row count is unchanged.

Expected:
- Step 6: console-error delta is 0 from the baseline; the table's row count equals the
  baseline row count (no rows were added or removed) — GROK-20010 guard.

### Scenario 6: Ctrl+V pastes into a different cell

Steps:
1. Close all and reopen demog to start from a clean state.
2. Click an AGE cell in row 1 to copy it (Ctrl+Shift+C). Note the value.
3. Click an AGE cell in row 10 to make it the paste target. Note its current value.
4. Press Ctrl+V.
5. Verify that row 10's AGE now holds the copied value from row 1.

Expected:
- Step 5: df.col('AGE').get(row10TableIdx) equals the value that was in row 1's
  AGE cell before pasting.

### Scenario 7: Shift+Del deletes selected rows; Ctrl+Z restores them

Steps:
1. Close all and reopen demog to start from a clean state. Record the table's row count.
2. Ctrl+click the row-headers of rows 3, 5, and 7 to select exactly 3 rows.
3. Press Shift+Del.
4. Verify that 3 rows have been removed.
5. Press Ctrl+Z.
6. Verify that the deleted rows have been restored.

Expected:
- Step 4: the table's row count equals the baseline minus 3.
- Step 6: the table's row count equals the original baseline (undo restored all 3 rows,
  routing through Command.runUndoable).

### Scenario 8: addNewRowOnLastRowEdit — auto-append on/off

Steps:
1. Close all and reopen demog to start from a clean state. Record the table's row count.
2. Open grid properties via the gear icon and enable "Add New Row On Last Row
   Edit". Close the properties panel.
3. Double-click the AGE cell in the very last row. Type a value and press Enter.
4. Verify that a new row was auto-appended.
5. Open grid properties and disable "Add New Row On Last Row Edit". Close the
   properties panel.
6. Double-click the AGE cell of what is now the last row. Type a value and press Enter.
7. Verify that no new row was appended.

Expected:
- Step 4: the table's row count equals the baseline plus 1 (auto-append fired).
- Step 7: the table's row count is still baseline plus 1 (no further auto-append when the
  setting is off).

## Automation notes
- The Dart editor-open handlers (double-click, Enter, digit input) and the plain-Delete handler are
  all gated on grid.currentGridCell, which ONLY the overlay's trusted mousedown hit-test assigns.
  Its JS getter reads `undefined` rather than a null GridCell, so it cannot be polled; the only
  observable proxy for "the click landed" is df.currentRowIdx plus df.currentCol.name. Settle on
  that after every seeding click, or the follow-up gesture outraces the assignment on a cold render
  and the open is silently dropped. When a genuine double-click still misses, re-seed with a
  trusted click and open with Enter — that path reuses the established cell with no second
  hit-test to race.
- the auto-append-on-last-row-edit case is manual and lives in the section's single manual checklist, grid-ui.md — not in a
  companion file beside this scenario. One list per section: a tester opening one file must not
  miss manual steps kept somewhere else.

### Signal discipline
All expected-results assertions use product-state signals:
- Edit state: `df.col('AGE').get(tableIdx)` (null for cleared cells, new value after commit)
- Event channel: `onCellValueEdited` event count delta
- Editor presence: DOM ValueEditor element presence/absence inside the grid root
- Clipboard content: clipboard text string (TSV for multi-row, single value for Ctrl+Shift+C)
- Row count: `df.rowCount` (for delete/undo and auto-append checks)
- Console errors: delta from baseline (0 for GROK-20010 guard)

### GROK-20010 guard (Ctrl+A → copy → paste)
Scenario 5 exercises the full select-all → copy → paste path to confirm the
the regression is absent. The guard is that the console stays clean — any error fails it,
not only the one the ticket happened to name.

### GROK-20266 guard (copy after Refresh)
The GROK-20266 guard (Toolbox Files Refresh breaks copy) is included in the atlas
cp description conditionally: "include only if recon shows a cheap Refresh path, else
the bug stays library-only." Recon must confirm whether a one-click Files Refresh is
reachable during automation; if it is, add it as Scenario 9 with a post-Refresh
Ctrl+C → Ctrl+V confirming the fresh value is pasted. This scenario is deferred
pending recon — no reliable automation path to a Files Refresh is confirmed. The
deferred step requires a live-browser check of the Toolbox > Files panel's Refresh
control location and availability.

### Shift+Del modifier guard
The atlas edge_case for Delete/Backspace vs Shift+Del is exercised by the interplay of
Scenario 2 (plain Delete clears a cell) and Scenario 7 (Shift+Del deletes rows) — the
two behaviors must be exercised and their effects confirmed separately to prove the
modifier is not mis-handled.

---
{
  "order": 13,
  "datasets": ["System:DemoFiles/demog.csv"]
}
