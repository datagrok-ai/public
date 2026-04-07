# Structure Rendering — Run Results

**Date**: 2026-04-07
**URL**: http://localhost:8888/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open SPGI | PASS | PASSED | 3624 rows, 88 columns loaded |
| 2 | Add viewers (scatterplot, histogram, line chart, bar chart, pie chart, trellis plot, box plot) | PASS | PASSED | 8 viewers total (Grid + 7 added) |
| 3-4 | Add structure legend (Core) to each viewer | PASS | PASSED | Scatterplot: Marker=Core, Color=Core; other viewers: split/color/category by Core |
| 5 | Scatterplot: Marker=Core, Color=Series — check structure rendering | PASS | PASSED | Color legend shows Series categories; marker legend shows 4 Core SMILES with rendered structure icons |
| 6 | Save and apply layout | PASS | PASSED | All 8 viewers restored with Color=Series, Marker=Core |
| 7 | Save and open project — verify legend and structure rendering | PASS | PASSED | All 8 viewers restored; scatterplot marker legend has 4 structure icons; structure rendering persists |

## Timing

| Phase | Duration |
|-------|----------|
| Execute via grok-browser | ~35s |
| Spec file generation | ~3s |
| Spec script execution | 26s |

## Summary

All steps passed. The Core column (containing SMILES) renders molecule structure images in the scatterplot marker legend. Both color (Series) and marker (Core) legends display correctly together. Structure rendering persists through layout save/restore. All 7 viewer types were successfully added and configured with Core-based legends. Note: Chem package semType detection did not fire automatically for SPGI — had to set `col.semType = 'Molecule'` manually.

## Retrospective

### What worked well
- Structure images render correctly in scatterplot marker legend items (4 structure icons detected)
- Combined color (Series) + marker (Core) legend displays both text and structure icons
- All 8 viewers restore correctly through layout save/restore
- Bar chart, histogram, line chart show Core SMILES in split legends
- Updated spec to use CDP connection for consistency

### What did not work
- Chem package semType detection did not fire automatically for SPGI dataset — had to manually set `col.semType = 'Molecule'` on Core and Structure columns
- `v3d.props.color` style property access doesn't work for some viewer types — used `setOptions()` instead

### Suggestions for the platform
- Investigate why Molecule semType detection fails for SPGI Core column (valid SMILES like `CC1(N[*:1])CCN([*:2])C1`)
- Consider rendering structure images in all viewer legends (histogram, line chart show only SMILES text, not images)

### Suggestions for the scenario
- Step 3 could specify which property to set for each viewer type (Color, Split, Category)
- Add a step to verify structure image rendering specifically (check for canvas/img elements in legend items)
- Step numbering skips 3 (goes 1,2,4,5,6,7) — should be sequential
