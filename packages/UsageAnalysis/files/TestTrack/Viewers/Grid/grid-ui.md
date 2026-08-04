# Grid — Manual Test Checklist


All scenarios should start with the following sequence of events:
1. Close all open views
2. Open the dataset specified in the scenario
3. Add Grid if the scenario requires a Grid viewer
---

## Context Panel — Column Hamburger Menu Inline Filter


*Dataset: spgi-100*

1. Open spgi-100.csv (System:AppData/Chem/tests/spgi-100.csv), then hover the CAST Idea ID column header and click the hamburger icon that appears at its right
   edge — a popup opens with the column name and the Filter, Actions, Colors and Dev sections
2. In the Filter section, drag the left handle of the strip under the histogram towards the middle
   — the filtered row count in the status bar drops straight away
3. Press Add filter — the column's filter appears in the Filter Panel on the left, keeping the
   range just set
4. Hover the Stereo Category header and open its hamburger — the Filter section lists the
   categories with their row counts instead of a histogram
5. Click the R_ONE category — only that category stays ticked and the filtered row count follows;
   the Context Panel spells out both filters now in force
6. Move the mouse outside the popup — the popup closes and both filters stay

## Summary Columns — Form Designer 

*Dataset: demog*

1. Right-click any data cell and choose Add > Summary Columns > Design a Form... — a form designer
   opens as its own view named Summary, and a form column has already been added to the grid
2. Drag one of the field boxes to an empty spot in the designer — the field settles where it was
   dropped
3. Press Close and Apply on the ribbon — the designer view closes, the table view comes back, and
   the grid keeps the form column with the field in its new place
4. Click the form column header — the Context Panel switches to that column and offers its Actions
   and Renderer settings together with an Edit button
5. Press Edit — the form designer reopens with the layout exactly as it was left

## Context Panel — Permissions 

*Dataset: spgi-100*

1. Click the Chemist column header — the Context Panel shows that column
2. Expand Advanced > Permissions — it offers an "Editable by" field and a "Pin if editable" switch
3. Type your own user name into "Editable by" and press Enter
4. Click a Chemist cell and press Enter — the cell editor opens, so the column is yours to edit
5. Press Esc, then replace your name in "Editable by" with a colleague's user name and press Enter
6. Click a Chemist cell and press Enter — no editor opens
7. Press Delete on that cell — the edit is refused and a message states the column is only editable
   by that colleague

## Auto-Append on Last-Row Edit 

*Dataset: demog*

1. Note the table's row count
2. Open the grid settings and switch on "Add New Row On Last Row Edit", then close the panel
3. Double-click the AGE cell of the very last row, type a value and press Enter — a new empty row
   appears at the bottom and the row count grows by one; the cell holds what you typed
4. Switch the setting back off, then repeat the edit on the new last row — this time no row is
   appended and the row count stays as it was, while the cell still takes the value

## Ctrl+Click Adds a Row to an Existing Range 

*Dataset: demog*

1. Select rows 5 through 10 by dragging down the row-number strip — six rows highlight and row 10
   is the current row
2. Hold Ctrl and click the row-number cell of row 15 — it joins the selection, making seven
3. Confirm the current row stays on row 10: a Ctrl+click changes membership only, it does not move
   the current row

## Multi-Row Copy — Clipboard Content

*Dataset: demog*

1. Select 5 rows by dragging down the row-number strip, noting the values you see in them
2. Press Ctrl+C, then paste into a plain-text editor (Notepad or any text field outside the
   platform)
3. Confirm the pasted text holds 5 data lines, one per selected row, in the order the rows appear
   in the grid
4. Confirm each line is tab-separated and its values follow the grid's current column order,
   matching what the selected rows show — including any column you moved or hid before copying

## Selected Rows Color  

*Dataset: demog*

1. Open the grid settings and set Selected Rows Color to a clearly recognizable color
2. Select a few rows by dragging down the row-number strip — the selected rows are painted in
   that color, not the default one
3. Press Esc — the rows return to their normal appearance

## Row Selection by Mouse — Drag and Shift+Click  

*Dataset: demog*

1. Click the row-number cell of row 5, then Shift+click the row-number cell of row 10 — rows 5
   through 10 highlight in orange and row 10 becomes the current row
2. Click and drag down the row-number strip across several rows — the swept rows highlight as
   the pointer moves
3. Click and drag across data cells in the grid body, sweeping several rows and columns — the
   swept ROWS highlight; confirm whether the column headers highlight too (on the current build
   a block drag selects rows only)
4. Click a single row-number cell without modifiers — the selection collapses to that row

## Column Selection by Mouse 

*Dataset: demog*

1. Shift+click the AGE column header, then the HEIGHT header — both columns show as selected
2. Ctrl+click the AGE header again — its selection inverts while HEIGHT stays selected
3. Press Esc — the column selection clears

## Tab Navigation at the Row Edges 


*Dataset: demog*

1. Navigate to the last data cell of a row and press Tab — the current cell does not wrap past
   the right edge
2. Navigate to the first data cell of a row and press Shift+Tab — the current cell does not wrap
   past the left edge

## Block Selection by Mouse Drag 

1. Click a data cell in the AGE column and drag across to the WEIGHT column and several rows
   down — the dragged rows highlight in orange
2. Confirm visually whether the column headers highlight as well
3. Click elsewhere to clear the selection

## Column Stats — Visual Verification  

*Prerequisite: add the min and max stats rows through Add > Column Stats*

1. Verify the Min row shows the minimum AGE, HEIGHT, and WEIGHT values
2. Verify the Max row shows the maximum AGE, HEIGHT, and WEIGHT values

## Column Groups — Visual Behavior  

*Dataset: demog*

1. Select several columns (Shift+click their headers) and group them from the Context Panel
   Actions — a group band appears above the grouped headers
2. Expand and collapse the group from its band — the member columns hide and reappear
3. Change the group's name and color — the band reflects both
4. Ungroup from the header context menu — the band disappears and the columns stay in place

---
{
  "order": 102,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
