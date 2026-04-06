# Structure Rendering — Run Results

**Date**: 2026-04-06
**URL**: http://localhost:8888/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | 3624 rows, 88 columns loaded |
| 2 | Add viewers (scatterplot, histogram, line chart, bar chart, pie chart, trellis plot, box plot) | PASS | PASSED | 8 viewers total (Grid + 7 added) |
| 3-4 | Add structure legend (Core) to each viewer | PASS | PASSED | Scatterplot: Marker=Core, Color=Core; other viewers: split/color by Core |
| 5 | Scatterplot: Marker=Core, Color=Series — check structure rendering | PASS | PASSED | Color legend shows Series categories; marker legend shows Core SMILES with structure icons (4 items with rendered molecule images) |
| 6 | Save and apply layout | PASS | PASSED | All 8 viewers restored with Color=Series, Marker=Core |
| 7 | Save and open project — verify legend and structure rendering | PASS | PASSED | All 8 viewers restored; scatterplot marker legend has 4 structure icons; structure rendering persists |

## Summary

All steps passed. The Core column (containing SMILES) renders molecule structure images in the scatterplot marker legend. Both color (Series) and marker (Core) legends display correctly together. Structure rendering persists through layout save/restore (project persistence). All 7 viewer types were successfully added and configured with Core-based legends.

## Retrospective

### What worked well
- Structure images render correctly in scatterplot marker legend items
- Combined color (Series) + marker (Core) legend displays both text and structure icons
- Trellis plot also renders structure images in its legend
- All 8 viewers restore correctly through layout save/restore
- Bar chart and histogram show Core SMILES in split legends with color differentiation

### What did not work
- Nothing — all steps passed

### Suggestions for the platform
- Consider rendering structure images in all viewer legends (histogram, line chart show only SMILES text, not images)

### Suggestions for the scenario
- Step 3 ("Add the structure legend to each viewer") could specify which property to set for each viewer type (Color, Split, Category, etc.) since different viewers use different property names
- Could add a step to verify structure image rendering specifically (check for canvas/img elements in legend items)
