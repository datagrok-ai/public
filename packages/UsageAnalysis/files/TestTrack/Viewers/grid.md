# Grid tests (Playwright)

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog

## Sorting

1. Double-click the AGE column header — column sorts descending
2. Double-click the AGE column header again — column sorts ascending
3. Double-click the AGE column header a third time — sort resets, no arrow icon
4. Right-click the grid and select Sort... — the Sort dialog opens
5. In the Sort dialog, set first sort by SEX ascending, click the plus icon and add second sort by AGE descending, click OK
6. Right-click the RACE column header, expand Sort, select Natural — RACE values are sorted in natural string order

## Column Resizing

1. Double-click the right border of the AGE column header — column auto-sizes to fit content
2. Right-click the grid, expand Column Sizing, select Optimal — all columns resize to optimal width

## Column Reordering and Hiding

1. Right-click any cell and select Order or Hide Columns... — the dialog opens
2. Uncheck the WEIGHT column checkbox and click OK — WEIGHT column disappears from the grid
3. Double-click the bold separator where WEIGHT was hidden — WEIGHT column reappears

## Row Selection

1. Click row 5 — row 5 becomes the current row
2. Hold Shift and click row 10 — rows 5 through 10 are selected
3. Hold Ctrl and click row 15 — row 15 is added to the selection
4. Hold Ctrl+Shift and click row 8 — row 8 is removed from the selection
5. Press Ctrl+A — all rows are selected
6. Click any single row without modifiers — selection clears, clicked row becomes current
7. Hold Shift and press Down arrow several times — rows are selected one by one downward
8. Press ESC — selection is cleared

## Column Selection

1. Shift+click the SEX column header — SEX column is selected
2. Shift+drag across the RACE and DIS_POP column headers — all three columns are selected
3. Ctrl+click the SEX column header — SEX column selection is inverted

## Cell Editing

1. Double-click a cell in the AGE column — an inline editor appears with the current value
2. Type a new numeric value and press Enter — the value is updated in the grid
3. Press Ctrl+Z — the edit is undone, original value is restored
4. Double-click a cell in the SEX column — an inline editor appears
5. Change the value and press Escape — the edit is cancelled, original value remains
6. Double-click a cell in the STARTED column — a date editor appears, change the value and press Enter
7. Press Ctrl+Shift+Z — redo brings back the last edit

## Copy and Paste

1. Click a cell in the AGE column to make it current
2. Press Ctrl+C — the cell value is copied to the clipboard
3. Navigate to another cell in the AGE column and press Ctrl+V — the copied value is pasted
4. Select multiple rows using Shift+click, then press Ctrl+C — "Block copied" message appears
5. Press Ctrl+Z to undo the paste — the original value is restored
6. Select several rows and press Shift+Del — the selected rows are deleted
7. Press Ctrl+Z — the deleted rows are restored

## Context Menu — Data Cell

1. Right-click on a data cell — context menu appears
2. Expand Column Sizing submenu, select Optimal — all columns resize to optimal width
3. Right-click again, expand Grid Color Coding, select All — all numeric columns become color-coded
4. Right-click again, expand Grid Color Coding, select None — all color coding is removed
5. Right-click again, expand Add > Top, select Histogram — a histogram miniature appears in the column header area

## Column Header Context Menu

1. Right-click the AGE column header — context menu appears with column-specific options
2. Expand Sort submenu, select Ascending — column sorts ascending
3. Right-click the AGE column header, expand Color Coding, select Linear — AGE column cells become color-coded with a gradient
4. Right-click the SEX column header, expand Color Coding, select Categorical — each category gets a distinct color
5. Right-click the HEIGHT column header, expand Color Coding, select Pick Up Coloring
6. Right-click the WEIGHT column header, expand Color Coding, select Apply Coloring — WEIGHT column gets the same color scheme as HEIGHT
7. Right-click the AGE column header, expand Format submenu — verify format options are shown
8. Right-click the AGE column header, select Column Properties... — the column properties dialog opens

## Column Cell Style (Renderer)

1. Right-click the AGE column header, expand Style, select Percent Completed — each AGE cell renders as a progress bar proportional to its value
2. Right-click the AGE column header, expand Style, select Default — cells revert to plain numeric text
3. Right-click the SEX column header, expand Style — only Default option is available for a categorical column

## Keyboard Navigation

1. Click any cell to set focus
2. Press Down arrow — current cell moves one row down
3. Press Right arrow — current cell moves to the next column
4. Press Ctrl+Home — current cell jumps to the first row, first column
5. Press Ctrl+End — current cell jumps to the last row
6. Press Home — current cell jumps to the first column in the current row
7. Press End — current cell jumps to the last column
8. Press Page Down — grid scrolls down one page
9. Press Page Up — grid scrolls back up
10. Press Tab — current cell moves one cell to the right, wrapping to the next row at row end
11. Press Shift+Tab — current cell moves one cell to the left, wrapping to the previous row

## Pinned Rows and Columns

1. Right-click on a data cell in row 3, expand Pin, select Pin Row — the row appears pinned at the top of the grid
2. Select several rows, right-click one of them, expand Pin, select Pin Selected Rows — all selected rows are pinned
3. Right-click any pinned row, expand Pin, select Unpin All Rows — all rows are unpinned
4. Right-click the AGE column header, expand Pin, select Pin Column — AGE column moves to the frozen area on the left
5. Right-click the pinned AGE column header, expand Pin, select Unpin Column — column returns to its original position

## Frozen Columns Properties

1. Open Grid settings via the gear icon
2. Set Frozen Columns to 2 — the first two columns become frozen
3. Set Frozen Columns back to 1 — only the row header column is frozen
4. Toggle Show Column Labels off — column header text disappears
5. Toggle Show Column Labels back on — headers reappear
6. Change Col Labels Orientation to Vert — column header text rotates to vertical
7. Change Col Labels Orientation back to Auto

## Color Coding

1. Open Grid settings via the gear icon, set Color Coding to All — all numerical columns become color-coded
2. Set Color Coding to None — all color coding is removed
3. Set Color Coding back to Auto
4. Right-click the HEIGHT column header, expand Color Coding, select Linear — cells show a gradient based on values
5. Right-click the RACE column header, expand Color Coding, select Categorical — each category gets a distinct color
6. Right-click the HEIGHT column header, expand Color Coding, select Pick Up Coloring
7. Right-click the WEIGHT column header, expand Color Coding, select Apply Coloring — WEIGHT column gets the same color scheme as HEIGHT
8. Right-click the RACE column header, expand Color Coding, select Off — color coding is removed from RACE

## Summary Columns

1. Right-click any data cell, expand Add > Summary Columns, select Sparkline — a new summary column appears with sparkline charts for each row
2. Click the sparkline column header — the Context Panel shows renderer settings
3. Right-click the sparkline column header and select Remove — the summary column is deleted

## Column Stats

1. Right-click any data cell, expand Add > Column Stats, select Min — a stats row labeled "Min" appears below the data
2. Right-click any data cell, expand Add > Column Stats, select Max — a Max stats row is added
3. Right-click any data cell, expand Add > Column Stats, deselect Min — the Min stats row is removed

## Column Header Hamburger Menu

1. Hover over the AGE column header until the hamburger (≡) icon appears in the top-right of the header
2. Click the hamburger icon — a popup panel opens with column statistics and action items
3. Move the mouse away — the popup dismisses

## Search

1. Press Ctrl+F — a search box appears at the top of the grid
2. Type "Asian" — matching cells in the RACE column are highlighted
3. Press Down arrow in the search box — cursor moves to the next match
4. Clear the search text and type "AGE > 50" — rows matching the expression are highlighted
5. Close the search box — grid returns to normal display

## Row State Synchronization

1. Add a Scatter Plot viewer alongside the Grid (both should be visible)
2. Click row 10 in the Grid — the corresponding point is highlighted as current in the Scatter Plot
3. Click a different data point in the Scatter Plot — the current row in the Grid updates accordingly
4. Select several rows in the Grid using Shift+click — the selected points are highlighted in the Scatter Plot
5. Select several rows in the Scatter Plot — the corresponding rows are highlighted in the Grid


## Layout and Project Save / Restore — Full Cycle with Visual Check

1. Apply color coding, resize some columns, adjust row height
2. Save the layout
3. Add additional viewers (scatter plot, histogram)
4. Apply the saved layout — verify all formatting, styling, coloring, and column order are preserved
5. Save the project
6. Close All
7. Open the saved project — verify all settings are intact

## Table Switching (Grid Viewer)

1. Open spgi-100
2. Add a Grid viewer
3. With the added Grid viewer selected, open the Context Panel
4. Go to Data > Table and switch to spgi-100 — viewer re-binds to the new table

---
{
  "order": 15,
  "datasets": ["System:DemoFiles/demog.csv", "System:AppData/Chem/tests/spgi-100.csv"]
}
