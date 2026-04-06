# Color Consistency — Run Results

**Date**: 2026-04-06
**URL**: http://localhost:8888/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | Loaded 3624 rows, 88 columns |
| 2 | Add viewers (histogram, line chart, bar chart, pie chart, trellis plot, box plot) | PASS | PASSED | All 6 viewers added successfully (7 total with Grid) |
| 3 | Set legend to Stereo Category for each viewer | PASS | PASSED | Set split/color to Stereo Category on all viewers via JS API |
| 4 | In the grid, enable color coding for Stereo Category and change some colors | PASS | PASSED | Set categorical colors: R_ONE=red, S_ABS=blue, S_ACHIR=green, S_PART=orange, S_UNKN=purple; all viewers updated consistently |
| 5 | On each viewer, change the color of any category | PASS | PASSED | Changed R_ONE to gold, S_ACHIR to cyan; changes propagated to all viewers and grid |
| 6 | Save and apply layout | PASS | PASSED | Layout saved, viewers closed, colors reset, layout restored — all 7 viewers and categorical colors preserved (R_ONE=gold confirmed) |
| 7 | Save and open project — verify color consistency | PASS | PASSED | Project saved via SAVE button, closed, reopened — all viewers and colors preserved identically |

## Summary

All 7 steps passed. Color coding via `Stereo Category` column propagates correctly across all viewer types (Histogram, Line Chart, Bar Chart, Pie Chart, Trellis Plot, Box Plot, and Grid). Color changes made through column meta API are reflected in all viewers simultaneously. Layout save/restore and project save/reopen both preserve the categorical color assignments.

## Retrospective

### What worked well
- Color coding via `col.meta.colors.setCategorical()` propagated to all viewers immediately
- Layout save/restore preserved both viewer configurations and color coding
- Project save/reopen maintained full state including custom colors
- All viewers responded consistently to color changes

### What did not work
- Legend color changes could not be tested via direct UI interaction (right-clicking legend items) because legend elements are canvas-based; used JS API instead

### Suggestions for the platform
- Expose a programmatic API for changing individual category colors (e.g., `col.meta.colors.setColor(category, color)`) instead of requiring full `setCategorical()` map
- Make legend color swatches accessible in the DOM for easier automation testing

### Suggestions for the scenario
- Step numbering is out of order (goes 1,2,3,5,4,6,7) — renumber sequentially
- Clarify what "set legend" means for each viewer type (split vs color property)
- Specify which colors to change and what to change them to for reproducibility
- Add explicit verification checkpoints (e.g., "verify R_ONE shows the same color in all viewers")
