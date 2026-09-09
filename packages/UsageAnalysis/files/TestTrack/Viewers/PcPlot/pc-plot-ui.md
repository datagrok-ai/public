---
feature: pcplot
target_layer: manual-only
coverage_type: smoke
manual_only_reason: |
  What remains here needs a human eye or trusted mouse input: hover colors and
  tooltips, window operations, context-menu consistency, and the grid-to-plot
  column drag.
---

# PC plot manual checklist

Human-only visual and gesture checks for the PC Plot.

All scenarios should start with the following sequence of events:
1. Close all
2. Open demog
3. Add PC plot

## Mouse-over and row group

1. Hover over a line in the plot -- the line should be highlighted with the mouse-over color, and a tooltip should appear
2. Enable **Selection > Show Mouse Over Row Group** -- hovering over a row highlights all rows sharing the same row group

## Viewer basics

1. Undock the viewer and move it around the screen
2. Dock the viewer in a different location
3. Open and close the Context Panel
4. Display the help window
5. Close the viewer and return it by **Edit > Undo** (Ctrl+Z)
6. Resize the viewer -- check mini-legend appearing and self-adjustable viewer layout

## Context menu: General and Tooltip

1. Right-click > **General** -- click through its items (Clone, Full Screen, Save to Gallery,
   Save as PNG, Embed..., and Close — leave Close for last, it removes the viewer) and verify
   each performs its named action; afterwards delete the gallery entry and discard the
   downloaded PNG
2. Right-click > **Tooltip** -- open **Edit...** and change the tooltip columns, then toggle
   the show/hide tooltip item and the remaining Tooltip items (e.g. the row-group
   tooltip) -- the hover tooltip follows each change
3. Open the Context Panel
4. Right-click > **Columns** -- verify the column set shown matches Context Panel > Value > Column Names

---
{
  "order": 101,
  "datasets": ["System:DemoFiles/demog.csv"]
}
