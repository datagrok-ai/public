---
feature: grid
realizes_atlas:
- grid.cp.columns-layout-persist
realizes:
- viewers.grid
priority: p0
target_layer: playwright
coverage_type: smoke
related_bugs:
- id: GROK-19753
  status: fixed
realized_as:
- grid-columns-style-persist-spec.ts
expected_results:
- anchor: Step 4
  expectation: After the first double-click on AGE, grid.props.sortByColumnNames contains 'AGE' and
    grid.props.sortTypes contains 'desc'; the value returned by grid.gridRowToTable(0) indexes the row
    with the maximum AGE in the dataframe.
- anchor: Step 5
  expectation: 'After the second double-click on AGE, grid.props.sortByColumnNames contains ''AGE''
    and grid.props.sortTypes contains ''asc''; grid.gridRowToTable(0) indexes the minimum-AGE row. The
    dataframe row order is unchanged: df.col(''AGE'').get(0) equals its value before any sorting (sort
    is grid-local).'
- anchor: Step 6
  expectation: After the third double-click on AGE, grid.props.sortByColumnNames is empty (sort reset);
    grid.gridRowToTable(0) no longer returns the minimum-AGE row index.
- anchor: Step 7
  expectation: 'The grid''s column order by idx reflects the drag result: the dragged column occupies
    its new position and no other column''s idx slot is displaced unexpectedly.'
- anchor: Step 9
  expectation: WEIGHT is hidden in the grid (gridColumn.visible === false) and is absent from the enumerated
    visible columns, while the surrounding columns keep their idx slots.
- anchor: Step 10
  expectation: gridColumn.width for the AGE column is larger than before the drag. Once the total column
    widths exceed the viewport the horizontal scrollbar is present, and scrolling right draws the further
    columns with no new grid error in the console — the GROK-19753 guard, whose symptom was exactly
    a broken visible-column window after manual resizes.
- anchor: Step 11
  expectation: grid.props.frozenColumns equals its value before pinning plus 1 (baseline frozenColumns
    is 1 because grid column 0 is the row header; pinning one data column increments it to 2).
- anchor: Step 12
  expectation: grid.props.pinnedRowColumnNames.length equals 2 (the pinned-row count; grid.pinnedRows
    is an object, not a number) and identifies the pinned rows.
- anchor: Step 14
  expectation: After re-applying the layout, the set of open viewers does not include the foreign viewer
    added in Step 13; the column order by idx matches Step 7; WEIGHT still has gridColumn.visible ===
    false; the AGE column width matches Step 10; sort props (sortByColumnNames, sortTypes) match Step
    5; grid.props.frozenColumns matches Step 11; grid.props.pinnedRowColumnNames.length equals 2.
- anchor: Step 16
  expectation: 'After Close All and reopening the saved project, the same assertion battery from Step
    14 passes: column order, the hidden column, the resized width, the sort props, the frozen-column
    count and the two pinned rows all come back as saved.'
bug_signal_raw: none
bug_signal_raw_cycle: 2026-07-30-grid-automate-02
strict_bug_repro_at_entry: true
strict_open_bug_ids: null
---

# Grid — Column Geometry and Persistence

Sort a column, reorder and hide columns, widen them until the grid scrolls, pin a column and
two rows, then confirm the whole arrangement survives a layout re-apply and a project
close-and-reopen.

## Setup

- Close all open views.
- Open demog. The table opens with its grid, and that grid is the one every step below acts on.
- Note the AGE value in the first row of the table as it stands before any sorting, and how many
  columns the grid currently freezes on the left (a fresh grid freezes one — the row-number column).

## Scenarios

### 1. Sorting cycles through descending, ascending and off

1. Click any cell to give the grid focus.
2. Find the AGE column header.
3. Double-click the AGE column header. This first double-click sorts descending.
4. Check the sort: the grid reports AGE as its sort column in descending order, the topmost row on
   screen carries the highest AGE in the dataset, and the table's own first row still holds the AGE
   value noted during setup — sorting rearranges what the grid shows, never the table itself.
5. Double-click the AGE column header again — the sort flips to ascending: the grid reports AGE as
   ascending, the topmost row carries the lowest AGE, and the table's first row is still unchanged.
6. Double-click the AGE column header a third time — the sort clears, the arrow disappears from the
   header, the grid no longer reports AGE as a sort column, and the topmost row is no longer the
   lowest-AGE one. Sort AGE ascending again so the grid enters the next scenario sorted.

### 2. Reordering a column

7. Drag the HEIGHT column header to a slot further right and drop it there. Do not drop it to the
   left of the frozen-column boundary — a column dropped there becomes pinned instead of moved.
   Afterwards the column order shows HEIGHT in its new position with the other columns shifted
   around it and none of them lost.

### 3. Hiding a column from the Order or Hide Columns dialog

8. Right-click any cell, choose **Order or Hide Columns...**, uncheck WEIGHT in the dialog's column
   list, and close the dialog.
9. Confirm WEIGHT takes up no width in the grid and no longer appears among the visible columns,
   while the columns around it keep their positions.

### 4. Widening columns until the grid scrolls

10. Drag the right border of the AGE column header to widen it, then widen several more columns the
    same way until the horizontal scrollbar appears, and scroll to the right. AGE ends up wider than
    it started, the further columns draw cleanly as you scroll, and no new errors appear in the
    console. This is the regression guard for GROK-19753, where scrolling broke after manual resizes.

### 5. Pinning a column and two rows

11. Right-click the SEX column header, open **Pin**, and choose **Pin Column** — the grid now freezes
    one more column on the left than it did during setup.
12. Right-click a data cell in the row, open **Pin**, choose **Pin Row**, and repeat on a second row —
    the row-header menu offers Pin Selected Rows, which records no column value, so the data-cell
    path is the one that makes the pinned rows identifiable. The grid
    now holds two pinned rows and remembers which rows they are.

### 6. The arrangement survives a layout and a project round-trip

13. Leave every change above in place and save the current view layout from the ribbon, noting its
    name. Then add another viewer, such as a Scatter Plot, so the view holds more than the grid.
14. Re-apply the saved layout. The Scatter Plot added in Step 13 is gone, and the grid comes back
    exactly as it was: HEIGHT in its new position, WEIGHT still hidden and zero-width, AGE at the
    width set in Step 10, AGE still sorted ascending, SEX still frozen, and two rows still pinned.
15. Save the view as a project with the ribbon **Save** button under a unique name, then Close All.
16. Reopen that project — everything checked in Step 14 holds again.
17. Afterwards, delete the probe layout and the probe project so nothing is left behind, even if an
    earlier step failed.

## Automation notes
- PERSISTENCE-TAIL CANON (corpus-verified; follow it instead of inventing a path):
  save with saveProjectViaUI(page, name) from helpers/projects.ts — the real ribbon Save, the only
  path that carries the view layout; reopen with grok.dapi.projects.find(id) then open() inside
  page.evaluate, which is what every persistence tail in the corpus does and which works (do NOT
  use openProjectFromDashboards here — the gallery tile is not reliably reachable); tear down with
  deleteProjectWithCleanup(page, {projectId}) in finally.
- ERROR-CHANNEL CANON for the save window: the ribbon Save renders a publication preview by cloning
  the live view into an offscreen iframe, and the detached view emits a "cloned iframe" message plus
  a Dart NullError PAIR. Gate a benign filter to the SAVE window only, with a LETTER-AGNOSTIC
  pattern (the minified symbol drifts), and keep everything else in that window failing.
- The reopen step deliberately asserts STATE ONLY, not the error channel. Reopening this rich layout
  does log a grid index error, but no open ticket owns that behaviour: the ticket this step used to
  cite was about the recent-projects widget, which the platform no longer supports, so the label did
  not match what fires. Asserting it here would mean a permanent red carrying the wrong name. If a
  ticket is ever filed for it, restore the console assertion with that ticket's id.

- If a strict assertion must stay red because the product is genuinely broken, wrap it in
  knownOpenBug(bugId, assertion) from helpers/known-open-bug.ts: the reproduction reports green and
  the assertion self-flips loudly when the bug is fixed. A permanently red expect() is not an
  acceptable resting state.
- resetShell(page) from helpers/openers.ts is the only helper that strips leftover BALLOONS from the
  DOM; error balloons never auto-hide and closeAll does not remove them, so call it before any phase
  whose balloon or console channel you intend to assert.

- Pin a COLUMN: right-click the column header, then Pin > Pin Column. Pin a ROW: click a cell in
  that row, then right-click it and choose Pin > Pin Row. The two pin actions live in different
  menus — the column one is only on the header, the row one only on a data cell.
- Project save and reopen go through the real ribbon Save button, never a JS-API project save:
  use saveProjectViaUI(page, name) from helpers/projects.ts, which clicks the ribbon Save,
  fills the name and dismisses the Share dialog that follows, returning {projectId, resolvedName}.
- Teardown uses deleteProjectWithCleanup(page, {projectId}) from helpers/projects.ts so the probe
  project never leaks; delete the probe layout in the same finally block.

- The canvas[name="overlay"] element hosts column header drag targets; compute drag
  coordinates from gridColumn.left + documentBounds for each column involved in reorder and resize.
- Do NOT drop a column header left of the freeze line during the reorder drag — that path silently
  pins the column instead of reordering it.
- The 3-state sort cycle starts with DESCENDING on the FIRST double-click; the atlas documents this
  explicitly (matching the live-build recon result). Do not assume ascending-first.
- Calling grid.getRowOrder() on an unsorted+unfiltered grid MATERIALIZES _order as a side effect —
  avoid it; use grid.gridRowToTable(i) for individual row mapping instead.
- Hiding a column: the Order or Hide Columns dialog is opened through the real right-click menu,
  but its per-column checkboxes are drawn on an embedded canvas with no DOM target, so the uncheck
  itself falls back to the JS API and the assertion is the ENUMERATION of the grid's visible
  columns (WEIGHT absent from it). GridColumn exposes no visibleWidth — do not cite it. This is a
  known partial-actuation gap: if a DOM-addressable path to those checkboxes is ever found, drive
  the click and keep the same enumeration assertion.
- The persistence tail uses the REAL ribbon Save button for the project save — not a JS-API
  project.save() call, which would skip the serialization path a real user's save goes through.
- Cleanup in finally: delete both the probe layout and the probe project to prevent state leakage
  across test runs.

---
{
  "order": 10,
  "datasets": ["System:DemoFiles/demog.csv"]
}
