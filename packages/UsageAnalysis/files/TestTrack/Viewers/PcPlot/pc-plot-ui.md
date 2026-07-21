# PC plot manual checklist

Human-only visual/gesture checks. The automatable steps (area selection, current-row line
highlight, the Show Current Line / Show Mouse Over Line / Show All Lines toggles, legend
add and positioning) are now covered automatically in the spec(s); what remains here needs
a human eye — hover colors and tooltips, "only selected lines visible" checks, window
operations, and context-menu consistency.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add PC plot

## Mouse-over and row group

1. Hover over a line in the plot -- the line should be highlighted with the mouse-over color, and a tooltip should appear
2. Enable **Selection > Show Mouse Over Row Group** -- hovering over a row highlights all rows sharing the same row group

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
6. Resize the viewer -- check mini-legend appearing and self-adjustable viewer layout

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
