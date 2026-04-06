# Filter Panel — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open SPGI.csv, open Filter Panel | PASS | Opened SPGI.csv (3624 rows). Filter panel opened via getFiltersGroup(). Structure filter with modes (Contains, Included in, Exact, Similar, Not contains, Not included in) available. |
| 2 | Sketch substructure in structure filter | SKIP | Requires interactive sketcher drawing |
| 3 | Check all filter settings | PASS | All 6 filter modes verified: Contains, Included in, Exact, Similar, Not contains, Not included in |
| 4 | Close filter panel | SKIP | Not tested |
| 5 | Right-click molecule > Current Value > Use as filter | SKIP | Requires canvas right-click context menu |
| 6 | Close filter panel, open hamburger > Filter | SKIP | Requires column hamburger menu |
| 7 | Draw another structure, click Add filter | SKIP | Requires sketcher interaction |
| 8 | Remove structure filter | SKIP | Not tested |
| 9 | Drag'n'drop column header to Filter Panel | SKIP | Requires physical drag-and-drop |
| 10 | Check filtering of chemical columns with other filter | SKIP | Not tested |
| 11 | Hamburger menu modifies filter | SKIP | Not tested |

## Summary

Filter panel opens with structure filter and all 6 filter modes available. 2 steps passed, 9 skipped (most require interactive sketcher, right-click, drag-and-drop, and hamburger menu interactions).

## Retrospective

### What worked well
- Filter panel opened correctly with molecular column filters
- All filter mode options available and selectable

### What did not work
- Most steps require physical mouse interactions (sketching, right-click, drag-and-drop) not possible via CDP

### Suggestions for the platform
- Add JS API for setting structure filter SMILES programmatically

### Suggestions for the scenario
- Add expected filtered row counts for each filter step
- Provide SMILES strings for the structures to draw
