# PC plot manual checklist

Manual checklist. Not included in PW automation.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add PC plot

## Selection (click on line)

1. Click on a line in the plot -- the corresponding row should be selected in the grid
2. Click on empty space -- selection should be cleared

## Area selection

1. On the viewer Shift-drag a rectangular area between two axes -- lines passing through the area should be selected
2. Shift-drag another area while holding Ctrl -- selection should be added to existing selection
3. Click on empty space -- selection should be cleared

## Current row and mouse-over

1. Click a row in the grid -- the corresponding line should be highlighted in the PC plot with the current line color
2. Hover over a line in the plot -- the line should be highlighted with the mouse-over color, and a tooltip should appear
3. Disable **Selection > Show Current Line** -- current row line no longer highlighted
4. Disable **Selection > Show Mouse Over Line** -- hover highlighting no longer shown
5. Enable **Selection > Show Mouse Over Row Group** -- hovering over a row highlights all rows sharing the same row group

## Show all lines with grid selection

1. Disable **Selection > Show All Lines**
2. Select some rows in the grid -- only selected lines should be shown on PC plot
3. Re-enable **Show All Lines** -- all lines visible again

## Viewer basics

1. Undock the viewer and move it around the screen
2. Dock the viewer in a different location
3. Open and close the Context Panel
4. Display the help window
5. Close the viewer and return it by **Edit > Undo** (Ctrl+Z)
6. Add a legend and change its position
7. Resize the viewer -- check mini-legend appearing and self-adjustable viewer layout

## Context menu: General and Tooltip

1. Right-click > **General** -- check all items
2. Right-click > **Tooltip** -- check all items
3. Open the Context Panel
4. Right-click > **Columns** -- verify that changes are consistent between the Context Panel and the context menu

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv"]
}
