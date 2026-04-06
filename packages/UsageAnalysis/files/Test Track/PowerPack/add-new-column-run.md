# Add New Columns — Run Results

**Date**: 2026-03-09
**URL**: https://release-ec2.datagrok.ai/
**Status**: PASS

## Steps

| # | Step | Result | Playwright | Notes |
|---|------|--------|------------|-------|
| 1 | Open the Demog Dataset | PASS | N/A | Opened via grok.data.getDemoTable('demog.csv'), 5850 rows |
| 2 | Press "Add new column" icon, dialog opens | PASS | N/A | Edit > Add New Column... opened dialog correctly |
| 3 | UI Check: no overlapping, proper resize | PASS | N/A | Dialog clean: columns list, functions list, preview grid, no visual issues |
| 4 | Add column "New" with Round(${HEIGHT}+${WEIGHT}) | PASS | N/A | Column created with correct rounded integer values (234, 257, 221, 223) |
| 5 | Recent Activities: select last formula, autofill | PASS | N/A | History icon showed 3 entries; clicking one autofilled name, type, and formula |

## Summary

All 5 steps passed. The Add New Column dialog works correctly with formula entry, autocomplete highlighting, interactive preview, and formula history/recall.

## Retrospective

### What worked well
- Dialog layout is clean with no visual issues
- Formula preview updates in real time
- Column type auto-detection works (detected int for Round result)
- History feature stores and recalls previous formulas with full autofill

### What did not work
- Nothing — all steps passed

### Suggestions for the platform
- N/A

### Suggestions for the scenario
- Step 4 says "using autocomplete hints and columns drag-n-drop" but drag-n-drop is hard to test via automation; consider separating UI-specific interactions into their own steps
