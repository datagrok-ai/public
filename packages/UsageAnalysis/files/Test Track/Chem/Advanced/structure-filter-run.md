# Structure Filter — Run Results

**Date**: 2026-03-11
**URL**: https://public.datagrok.ai
**Status**: PARTIAL

## Steps

| # | Step | Result | Notes |
|---|------|--------|-------|
| 1 | Open SPGI.csv, open filter panel | PASS | Dataset opened (3624 rows, 88 cols). Filter panel shows structure filters for Structure, Core, R1-R3, R100, R101 columns with Sketch buttons. |
| 2 | Check filter settings (Contains, Included in, Exact, Similar) | PASS | All 6 filter modes available: Contains, Included in, Exact, Similar, Not contains, Not included in |
| 3 | Add structure filter, disable, close/reopen panel | SKIP | Requires sketcher interaction for drawing structures |
| 4 | Set filter, close panel, use Current value > Use as filter | SKIP | Requires right-click context menu on molecule cell |
| 5 | Set filter, close panel, hamburger menu > Filter, draw another | SKIP | Requires column hamburger menu interaction |
| 6 | Remove filter, Current value > Use as filter | SKIP | Requires right-click context menu |
| 7 | Open linked datasets, clone view, sync filters | SKIP | Complex multi-view scenario |

## Summary

Filter panel opens correctly with all structure columns detected and all 6 filter modes available. 2 steps passed, 5 skipped (require interactive sketcher drawing and right-click context menu interactions that are difficult to automate via CDP).

## Retrospective

### What worked well
- Filter panel opened automatically with all molecule columns detected
- All filter mode options (Contains, Included in, Exact, Similar, Not contains, Not included in) available
- Multiple structure columns (Structure, Core, R1-R3, R100, R101) all have individual filter controls

### What did not work
- Could not test sketcher-based filtering (requires drawing in canvas)
- Could not test right-click > Current value > Use as filter (requires physical mouse context menu)

### Suggestions for the platform
- Add JS API for programmatically setting structure filter values
- Consider adding a SMILES input field to the structure filter (like Scaffold Tree has)

### Suggestions for the scenario
- Scenario steps are numbered confusingly (multiple "1." restarts)
- Should use consistent step numbering
- Add expected filtered row counts for each step
